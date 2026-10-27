// Copyright (c) 2022 Ultimaker B.V.
// CuraEngine is released under the terms of the AGPLv3 or higher.

#include <algorithm> //For std::partition_copy and std::min_element.
#include <cstddef>
#include <cstdint>
#include <memory>
#include <cassert>
#include <cmath>
#include <math.h>
#include <limits>
#include <unordered_set>

#include "WallToolPaths.hpp"

#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Arachne/utils/ExtrusionLine.hpp"
#include "libslic3r/Line.hpp"
#include "libslic3r/Arachne/utils/PolygonsPointIndex.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Arachne/BeadingStrategy/BeadingStrategyFactory.hpp"
#include "libslic3r/Arachne/utils/ExtrusionJunction.hpp"
#include "libslic3r/Arachne/utils/SparsePointGrid.hpp"
#include "libslic3r/Arachne/utils/SquareGrid.hpp"
#include "SkeletalTrapezoidation.hpp"
#include "../ClipperUtils.hpp"
#include "libslic3r/libslic3r.h"
#include "utils/linearAlg2D.hpp"
#include "EdgeGrid.hpp"
#include "utils/SparseLineGrid.hpp"
#include "Geometry.hpp"
#include "utils/PolylineStitcher.hpp"
#include "SVG.hpp"
#include "Utils.hpp"

#include <boost/log/trivial.hpp>
#include <utility>
#include <vector>

//#define ARACHNE_STITCH_PATCH_DEBUG

namespace Slic3r::Arachne
{

WallToolPathsParams make_paths_params(const int layer_id, const PrintObjectConfig &print_object_config, const PrintConfig &print_config)
{
    WallToolPathsParams input_params;
    {
        const double min_nozzle_diameter = *std::min_element(print_config.nozzle_diameter.values.begin(), print_config.nozzle_diameter.values.end());
        if (const auto &min_feature_size_opt = print_object_config.min_feature_size)
            input_params.min_feature_size = min_feature_size_opt.value * 0.01 * min_nozzle_diameter;

        if (const auto &min_wall_length_factor_opt = print_object_config.min_length_factor)
            input_params.min_length_factor = min_wall_length_factor_opt.value;
        else
            input_params.min_length_factor = 0.5f;

        // ORCA: in mm, or as a percentage of the nozzle diameter.
        if (layer_id == 0) {
            if (print_object_config.initial_layer_min_bead_width.value > 0.)
                input_params.min_bead_width = float(print_object_config.initial_layer_min_bead_width.get_abs_value(min_nozzle_diameter));
        } else {
            if (print_object_config.min_bead_width.value > 0.)
                input_params.min_bead_width = float(print_object_config.min_bead_width.get_abs_value(min_nozzle_diameter));
        }

        if (const auto &wall_transition_filter_deviation_opt = print_object_config.wall_transition_filter_deviation)
            input_params.wall_transition_filter_deviation = wall_transition_filter_deviation_opt.value * 0.01 * min_nozzle_diameter;

        if (const auto &wall_transition_length_opt = print_object_config.wall_transition_length)
            input_params.wall_transition_length = wall_transition_length_opt.value * 0.01 * min_nozzle_diameter;

        input_params.wall_transition_angle   = print_object_config.wall_transition_angle.value;
        input_params.wall_distribution_count = print_object_config.wall_distribution_count.value;

        input_params.is_top_or_bottom_layer = false; // Set to default value

        if (const auto& wall_maximum_resolution_opt = print_object_config.wall_maximum_resolution)
            input_params.wall_maximum_resolution = scaled<coord_t>(wall_maximum_resolution_opt.value);

        if (const auto& wall_maximum_deviation_opt = print_object_config.wall_maximum_deviation)
            input_params.wall_maximum_deviation = scaled<coord_t>(wall_maximum_deviation_opt.value);
    }

    return input_params;
}

WallToolPaths::WallToolPaths(const Polygons& outline, const coord_t bead_width_0, const coord_t bead_width_x,
                             const size_t inset_count, const coord_t wall_0_inset, const coordf_t layer_height, const WallToolPathsParams &params)
    : outline(outline)
    , bead_width_0(bead_width_0)
    , bead_width_x(bead_width_x)
    , inset_count(inset_count)
    , wall_0_inset(wall_0_inset)
    , layer_height(layer_height)
    , print_thin_walls(Slic3r::Arachne::fill_outline_gaps)
    , min_feature_size(scaled<coord_t>(params.min_feature_size))
    , min_bead_width(scaled<coord_t>(params.min_bead_width))
    , small_area_length(static_cast<double>(bead_width_0) / 2.)
    , wall_transition_filter_deviation(scaled<coord_t>(params.wall_transition_filter_deviation))
    , toolpaths_generated(false)
    , m_params(params)
{
}

void simplify(Polygon &thiss, const int64_t smallest_line_segment_squared, const int64_t allowed_error_distance_squared)
{
    if (thiss.size() < 3) {
        thiss.points.clear();
        return;
    }
    if (thiss.size() == 3)
        return;

    Polygon new_path;
    Point previous = thiss.points.back();
    Point previous_previous = thiss.points.at(thiss.points.size() - 2);
    Point current = thiss.points.at(0);

    /* When removing a vertex, we check the height of the triangle of the area
     being removed from the original polygon by the simplification. However,
     when consecutively removing multiple vertices the height of the previously
     removed vertices w.r.t. the shortcut path changes.
     In order to not recompute the new height value of previously removed
     vertices we compute the height of a representative triangle, which covers
     the same amount of area as the area being cut off. We use the Shoelace
     formula to accumulate the area under the removed segments. This works by
     computing the area in a 'fan' where each of the blades of the fan go from
     the origin to one of the segments. While removing vertices the area in
     this fan accumulates. By subtracting the area of the blade connected to
     the short-cutting segment we obtain the total area of the cutoff region.
     From this area we compute the height of the representative triangle using
     the standard formula for a triangle area: A = .5*b*h
     */
    int64_t accumulated_area_removed = int64_t(previous.x()) * int64_t(current.y()) - int64_t(previous.y()) * int64_t(current.x()); // Twice the Shoelace formula for area of polygon per line segment.

    for (size_t point_idx = 0; point_idx < thiss.points.size(); point_idx++) {
        current = thiss.points.at(point_idx % thiss.points.size());

        //Check if the accumulated area doesn't exceed the maximum.
        Point next;
        if (point_idx + 1 < thiss.points.size()) {
            next = thiss.points.at(point_idx + 1);
        } else if (point_idx + 1 == thiss.points.size() && new_path.size() > 1) { // don't spill over if the [next] vertex will then be equal to [previous]
            next = new_path[0]; //Spill over to new polygon for checking removed area.
        } else {
            next = thiss.points.at((point_idx + 1) % thiss.points.size());
        }
        const int64_t removed_area_next = int64_t(current.x()) * int64_t(next.y()) - int64_t(current.y()) * int64_t(next.x()); // Twice the Shoelace formula for area of polygon per line segment.
        const int64_t negative_area_closing = int64_t(next.x()) * int64_t(previous.y()) - int64_t(next.y()) * int64_t(previous.x()); // area between the origin and the short-cutting segment
        accumulated_area_removed += removed_area_next;

        const int64_t length2 = (current - previous).cast<int64_t>().squaredNorm();

        // Orca:
        // Checking if the segment's length is smaller than 5 microns (0.005mm).
        // The value of `length2` is scaled and squared, so we need to compare it with the squared value of 5 microns
        if (length2 < Slic3r::sqr(scaled<coord_t>(0.005))) {
            // We're allowed to always delete segments of less than 5 micron.
            continue;
        }

        const int64_t area_removed_so_far = accumulated_area_removed + negative_area_closing; // close the shortcut area polygon
        const int64_t base_length_2 = (next - previous).cast<int64_t>().squaredNorm();

        if (base_length_2 == 0) //Two line segments form a line back and forth with no area.
            continue; //Remove the vertex.
        //We want to check if the height of the triangle formed by previous, current and next vertices is less than allowed_error_distance_squared.
        //1/2 L = A           [actual area is half of the computed shoelace value] // Shoelace formula is .5*(...) , but we simplify the computation and take out the .5
        //A = 1/2 * b * h     [triangle area formula]
        //L = b * h           [apply above two and take out the 1/2]
        //h = L / b           [divide by b]
        //h^2 = (L / b)^2     [square it]
        //h^2 = L^2 / b^2     [factor the divisor]
        const int64_t height_2 = double(area_removed_so_far) * double(area_removed_so_far) / double(base_length_2);
        // Orca: The value of `height_2` is squared, so we need to compare it with the squared value
        if ((height_2 <= Slic3r::sqr(colinear_vertex_tolerance()) //Almost exactly colinear (barring rounding errors).
             && Line::distance_to_infinite(current, previous, next) <= double(colinear_vertex_tolerance()))) // make sure that height_2 is not small because of cancellation of positive and negative areas
            continue;

        if (length2 < smallest_line_segment_squared
            && height_2 <= allowed_error_distance_squared) // removing the vertex doesn't introduce too much error.)
        {
            const int64_t next_length2 = (current - next).cast<int64_t>().squaredNorm();
            if (next_length2 > 4 * smallest_line_segment_squared) {
                // Special case; The next line is long. If we were to remove this, it could happen that we get quite noticeable artifacts.
                // We should instead move this point to a location where both edges are kept and then remove the previous point that we wanted to keep.
                // By taking the intersection of these two lines, we get a point that preserves the direction (so it makes the corner a bit more pointy).
                // We just need to be sure that the intersection point does not introduce an artifact itself.
                Point intersection_point;
                bool has_intersection = Line(previous_previous, previous).intersection_infinite(Line(current, next), &intersection_point);
                if (!has_intersection
                    || Line::distance_to_infinite_squared(intersection_point, previous, current) > double(allowed_error_distance_squared)
                    || (intersection_point - previous).cast<int64_t>().squaredNorm() > smallest_line_segment_squared  // The intersection point is way too far from the 'previous'
                    || (intersection_point - next).cast<int64_t>().squaredNorm() > smallest_line_segment_squared)     // and 'next' points, so it shouldn't replace 'current'
                {
                    // We can't find a better spot for it, but the size of the line is more than 5 micron.
                    // So the only thing we can do here is leave it in...
                }
                else {
                    // New point seems like a valid one.
                    current = intersection_point;
                    // If there was a previous point added, remove it.
                    if(!new_path.empty()) {
                        new_path.points.pop_back();
                        previous = previous_previous;
                    }
                }
            } else {
                continue; //Remove the vertex.
            }
        }
        //Don't remove the vertex.
        accumulated_area_removed = removed_area_next; // so that in the next iteration it's the area between the origin, [previous] and [current]
        previous_previous = previous;
        previous = current; //Note that "previous" is only updated if we don't remove the vertex.
        new_path.points.push_back(current);
    }

    thiss = new_path;
}

/*!
     * Removes vertices of the polygons to make sure that they are not too high
     * resolution.
     *
     * This removes points which are connected to line segments that are shorter
     * than the `smallest_line_segment`, unless that would introduce a deviation
     * in the contour of more than `allowed_error_distance`.
     *
     * Criteria:
     * 1. Never remove a vertex if either of the connceted segments is larger than \p smallest_line_segment
     * 2. Never remove a vertex if the distance between that vertex and the final resulting polygon would be higher than \p allowed_error_distance
     * 3. The direction of segments longer than \p smallest_line_segment always
     * remains unaltered (but their end points may change if it is connected to
     * a small segment)
     *
     * Simplify uses a heuristic and doesn't neccesarily remove all removable
     * vertices under the above criteria, but simplify may never violate these
     * criteria. Unless the segments or the distance is smaller than the
     * rounding error of 5 micron.
     *
     * Vertices which introduce an error of less than 5 microns are removed
     * anyway, even if the segments are longer than the smallest line segment.
     * This makes sure that (practically) colinear line segments are joined into
     * a single line segment.
     * \param smallest_line_segment Maximal length of removed line segments.
     * \param allowed_error_distance If removing a vertex introduces a deviation
     * from the original path that is more than this distance, the vertex may
     * not be removed.
 */
void simplify(Polygons &thiss, const int64_t smallest_line_segment = scaled<coord_t>(0.01), const int64_t allowed_error_distance = scaled<coord_t>(0.005))
{
    const int64_t allowed_error_distance_squared = int64_t(allowed_error_distance) * int64_t(allowed_error_distance);
    const int64_t smallest_line_segment_squared = int64_t(smallest_line_segment) * int64_t(smallest_line_segment);
    for (size_t p = 0; p < thiss.size(); p++)
    {
        simplify(thiss[p], smallest_line_segment_squared, allowed_error_distance_squared);
        if (thiss[p].size() < 3)
        {
            thiss.erase(thiss.begin() + p);
            p--;
        }
    }
}

typedef SparseLineGrid<PolygonsPointIndex, PolygonsPointIndexSegmentLocator> LocToLineGrid;
std::unique_ptr<LocToLineGrid>                                               createLocToLineGrid(const Polygons &polygons, int square_size)
{
    unsigned int n_points = 0;
    for (const auto &poly : polygons)
        n_points += poly.size();

    auto ret = std::make_unique<LocToLineGrid>(square_size, n_points);

    for (unsigned int poly_idx = 0; poly_idx < polygons.size(); poly_idx++)
        for (unsigned int point_idx = 0; point_idx < polygons[poly_idx].size(); point_idx++)
            ret->insert(PolygonsPointIndex(&polygons, poly_idx, point_idx));
    return ret;
}

/* Note: Also tries to solve for near-self intersections, when epsilon >= 1
 */
void fixSelfIntersections(const coord_t epsilon, Polygons &thiss)
{
    if (epsilon < 1)
        return;

    const int64_t half_epsilon = (epsilon + 1) / 2;

    // Points too close to line segments should be moved a little away from those line segments, but less than epsilon,
    //   so at least half-epsilon distance between points can still be guaranteed.
    const coord_t grid_size  = scaled<coord_t>(2.);
    auto              query_grid = createLocToLineGrid(thiss, grid_size);

    const auto    move_dist         = std::max<int64_t>(2L, half_epsilon - 2);
    const int64_t half_epsilon_sqrd = half_epsilon * half_epsilon;

    const size_t n = thiss.size();
    for (size_t poly_idx = 0; poly_idx < n; poly_idx++) {
        const size_t pathlen = thiss[poly_idx].size();
        for (size_t point_idx = 0; point_idx < pathlen; ++point_idx) {
            Point &pt = thiss[poly_idx][point_idx];
            for (const auto &line : query_grid->getNearby(pt, epsilon)) {
                const size_t line_next_idx = (line.point_idx + 1) % thiss[line.poly_idx].size();
                if (poly_idx == line.poly_idx && (point_idx == line.point_idx || point_idx == line_next_idx))
                    continue;

                const Line segment(thiss[line.poly_idx][line.point_idx], thiss[line.poly_idx][line_next_idx]);
                Point      segment_closest_point;
                segment.distance_to_squared(pt, &segment_closest_point);

                if (half_epsilon_sqrd >= (pt - segment_closest_point).cast<int64_t>().squaredNorm()) {
                    const Point  &other = thiss[poly_idx][(point_idx + 1) % pathlen];
                    const Vec2i64 vec   = (LinearAlg2D::pointIsLeftOfLine(other, segment.a, segment.b) > 0 ? segment.b - segment.a : segment.a - segment.b).cast<int64_t>();
                    assert(Slic3r::sqr(double(vec.x())) < double(std::numeric_limits<int64_t>::max()));
                    assert(Slic3r::sqr(double(vec.y())) < double(std::numeric_limits<int64_t>::max()));
                    const int64_t len   = vec.norm();
                    pt.x() += (-vec.y() * move_dist) / len;
                    pt.y() += (vec.x() * move_dist) / len;
                }
            }
        }
    }
}

/*!
     * Removes overlapping consecutive line segments which don't delimit a positive area.
 */
void removeDegenerateVerts(Polygons &thiss)
{
    for (size_t poly_idx = 0; poly_idx < thiss.size(); poly_idx++) {
        Polygon &poly = thiss[poly_idx];
        Polygon  result;

        auto isDegenerate = [](const Point &last, const Point &now, const Point &next) {
            Vec2i64 last_line = (now - last).cast<int64_t>();
            Vec2i64 next_line = (next - now).cast<int64_t>();
            return last_line.dot(next_line) == -1 * last_line.norm() * next_line.norm();
        };
        bool isChanged = false;
        for (size_t idx = 0; idx < poly.size(); idx++) {
            const Point &last = (result.size() == 0) ? poly.back() : result.back();
            if (idx + 1 == poly.size() && result.size() == 0)
                break;

            const Point &next = (idx + 1 == poly.size()) ? result[0] : poly[idx + 1];
            if (isDegenerate(last, poly[idx], next)) { // lines are in the opposite direction
                // don't add vert to the result
                isChanged = true;
                while (result.size() > 1 && isDegenerate(result[result.size() - 2], result.back(), next))
                    result.points.pop_back();
            } else {
                result.points.emplace_back(poly[idx]);
            }
        }

        if (isChanged) {
            if (result.size() > 2) {
                poly = result;
            } else {
                thiss.erase(thiss.begin() + poly_idx);
                poly_idx--; // effectively the next iteration has the same poly_idx (referring to a new poly which is not yet processed)
            }
        }
    }
}

void removeSmallAreas(Polygons &thiss, const double min_area_size, const bool remove_holes)
{
    auto new_end = thiss.end();
    if (remove_holes) {
        for (auto it = thiss.begin(); it < new_end;) {
            // All polygons smaller than target are removed by replacing them with a polygon from the back of the vector.
            if (fabs(it->area()) < min_area_size) {
                --new_end;
                *it = std::move(*new_end);
                continue; // Don't increment the iterator such that the polygon just swapped in is checked next.
            }
            ++it;
        }
    } else {
        // For each polygon, computes the signed area, move small outlines at the end of the vector and keep pointer on small holes
        Polygons small_holes;
        for (auto it = thiss.begin(); it < new_end;) {
            if (double area = it->area(); fabs(area) < min_area_size) {
                if (area >= 0) {
                    --new_end;
                    if (it < new_end) {
                        std::swap(*new_end, *it);
                        continue;
                    } else { // Don't self-swap the last Path
                        break;
                    }
                } else {
                    small_holes.push_back(*it);
                }
            }
            ++it;
        }

        // Removes small holes that have their first point inside one of the removed outlines
        // Iterating in reverse ensures that unprocessed small holes won't be moved
        const auto removed_outlines_start = new_end;
        for (auto hole_it = small_holes.rbegin(); hole_it < small_holes.rend(); hole_it++)
            for (auto outline_it = removed_outlines_start; outline_it < thiss.end(); outline_it++)
                if (Polygon(*outline_it).contains(*hole_it->begin())) {
                    new_end--;
                    *hole_it = std::move(*new_end);
                    break;
                }
    }
    thiss.resize(new_end-thiss.begin());
}

void removeColinearEdges(Polygon &poly, const double max_deviation_angle)
{
    // TODO: Can be made more efficient (for example, use pointer-types for process-/skip-indices, so we can swap them without copy).
    size_t num_removed_in_iteration = 0;
    do {
        num_removed_in_iteration = 0;
        std::vector<bool> process_indices(poly.points.size(), true);

        bool go = true;
        while (go) {
            go = false;

            const auto  &rpath   = poly;
            const size_t pathlen = rpath.size();
            if (pathlen <= 3)
                return;

            std::vector<bool> skip_indices(poly.points.size(), false);

            Polygon new_path;
            for (size_t point_idx = 0; point_idx < pathlen; ++point_idx) {
                // Don't iterate directly over process-indices, but do it this way, because there are points _in_ process-indices that should nonetheless
                // be skipped:
                if (!process_indices[point_idx]) {
                    new_path.points.push_back(rpath[point_idx]);
                    continue;
                }

                // Should skip the last point for this iteration if the old first was removed (which can be seen from the fact that the new first was skipped):
                if (point_idx == (pathlen - 1) && skip_indices[0]) {
                    skip_indices[new_path.size()] = true;
                    go                            = true;
                    new_path.points.push_back(rpath[point_idx]);
                    break;
                }

                const Point &prev = rpath[(point_idx - 1 + pathlen) % pathlen];
                const Point &pt   = rpath[point_idx];
                const Point &next = rpath[(point_idx + 1) % pathlen];

                float angle = LinearAlg2D::getAngleLeft(prev, pt, next); // [0 : 2 * pi]
                if (angle >= float(M_PI)) { angle -= float(M_PI); }                    // map [pi : 2 * pi] to [0 : pi]

                // Check if the angle is within limits for the point to 'make sense', given the maximum deviation.
                // If the angle indicates near-parallel segments ignore the point 'pt'
                if (angle > max_deviation_angle && angle < M_PI - max_deviation_angle) {
                    new_path.points.push_back(pt);
                } else if (point_idx != (pathlen - 1)) {
                    // Skip the next point, since the current one was removed:
                    skip_indices[new_path.size()] = true;
                    go                            = true;
                    new_path.points.push_back(next);
                    ++point_idx;
                }
            }
            poly = new_path;
            num_removed_in_iteration += pathlen - poly.points.size();

            process_indices.clear();
            process_indices.insert(process_indices.end(), skip_indices.begin(), skip_indices.end());
        }
    } while (num_removed_in_iteration > 0);
}

void removeColinearEdges(Polygons &thiss, const double max_deviation_angle = 0.0005)
{
    for (int p = 0; p < int(thiss.size()); p++) {
        removeColinearEdges(thiss[p], max_deviation_angle);
        if (thiss[p].size() < 3) {
            thiss.erase(thiss.begin() + p);
            p--;
        }
    }
}

const std::vector<VariableWidthLines> &WallToolPaths::generate()
{
    if (this->inset_count < 1)
        return toolpaths;

    const coord_t smallest_segment = m_params.wall_maximum_resolution;
    const coord_t allowed_distance = m_params.wall_maximum_deviation;
    const coord_t epsilon_offset = (allowed_distance / 2) - 1;
    const double  transitioning_angle = Geometry::deg2rad(m_params.wall_transition_angle);
    const coord_t discretization_step_size = scaled<coord_t>(0.8);

    // Simplify outline for boost::voronoi consumption. Absolutely no self intersections or near-self intersections allowed:
    // TODO: Open question: Does this indeed fix all (or all-but-one-in-a-million) cases for manifold but otherwise possibly complex polygons?
    Polygons prepared_outline = offset(offset(offset(outline, -epsilon_offset), epsilon_offset * 2), -epsilon_offset);
    simplify(prepared_outline, smallest_segment, allowed_distance);
    fixSelfIntersections(epsilon_offset, prepared_outline);
    removeDegenerateVerts(prepared_outline);
    removeColinearEdges(prepared_outline, 0.005);
    // Removing collinear edges may introduce self intersections, so we need to fix them again
    fixSelfIntersections(epsilon_offset, prepared_outline);
    removeDegenerateVerts(prepared_outline);
    removeSmallAreas(prepared_outline, small_area_length * small_area_length, false);

    // The functions above could produce intersecting polygons that could cause a crash inside Arachne.
    // Applying Clipper union should be enough to get rid of this issue.
    // Clipper union also fixed an issue in Arachne that in post-processing Voronoi diagram, some edges
    // didn't have twin edges. (a non-planar Voronoi diagram probably caused this).
    prepared_outline = union_(prepared_outline);

    if (area(prepared_outline) <= 0) {
        assert(toolpaths.empty());
        return toolpaths;
    }

    const float external_perimeter_extrusion_width = Flow::rounded_rectangle_extrusion_width_from_spacing(unscale<float>(bead_width_0), float(this->layer_height));
    const float perimeter_extrusion_width          = Flow::rounded_rectangle_extrusion_width_from_spacing(unscale<float>(bead_width_x), float(this->layer_height));

    const coord_t wall_transition_length = scaled<coord_t>(this->m_params.wall_transition_length);
	
	const double wall_split_middle_threshold = m_params.wall_split_middle_threshold > 0. ? m_params.wall_split_middle_threshold :
        std::clamp(2. * unscaled<double>(this->min_bead_width) / external_perimeter_extrusion_width - 1., 0.01, 0.99); // For an uneven nr. of lines: When to split the middle wall into two.
    const double wall_add_middle_threshold   = m_params.wall_add_middle_threshold > 0. ? m_params.wall_add_middle_threshold :
        std::clamp(unscaled<double>(this->min_bead_width) / perimeter_extrusion_width, 0.01, 0.99); // For an even nr. of lines: When to add a new middle in between the innermost two walls.
    
    const int wall_distribution_count = this->m_params.wall_distribution_count;
    const size_t max_bead_count = (inset_count < std::numeric_limits<coord_t>::max() / 2) ? 2 * inset_count : std::numeric_limits<coord_t>::max();
    const auto beading_strat = BeadingStrategyFactory::makeStrategy
        (
            bead_width_0,
            bead_width_x,
            wall_transition_length,
            transitioning_angle,
            print_thin_walls,
            min_bead_width,
            min_feature_size,
            wall_split_middle_threshold,
            wall_add_middle_threshold,
            max_bead_count,
            wall_0_inset,
            wall_distribution_count,
            0.5,
            m_params.overflow_early,
            m_params.transition_thresholds
        );
    const coord_t transition_filter_dist   = scaled<coord_t>(100.f);
    const coord_t allowed_filter_deviation = wall_transition_filter_deviation;
    SkeletalTrapezoidation wall_maker
    (
        prepared_outline,
        *beading_strat,
        beading_strat->getTransitioningAngle(),
        discretization_step_size,
        transition_filter_dist,
        allowed_filter_deviation,
        wall_transition_length
    );
    wall_maker.generateToolpaths(toolpaths);

    stitchToolPaths(toolpaths, this->bead_width_x);

    removeSmallLines(toolpaths);

    separateOutInnerContour();

    simplifyToolPaths(toolpaths, m_params);

    removeEmptyToolPaths(toolpaths);
    assert(std::is_sorted(toolpaths.cbegin(), toolpaths.cend(),
                          [](const VariableWidthLines& l, const VariableWidthLines& r)
                          {
                              return l.front().inset_idx < r.front().inset_idx;
                          }) && "WallToolPaths should be sorted from the outer 0th to inner_walls");
    toolpaths_generated = true;
    return toolpaths;
}

void WallToolPaths::stitchToolPaths(std::vector<VariableWidthLines> &toolpaths, const coord_t bead_width_x)
{
    const coord_t stitch_distance = bead_width_x - 1; //In 0-width contours, junctions can cause up to 1-line-width gaps. Don't stitch more than 1 line width.

    for (unsigned int wall_idx = 0; wall_idx < toolpaths.size(); wall_idx++) {
        VariableWidthLines& wall_lines = toolpaths[wall_idx];

        VariableWidthLines stitched_polylines;
        VariableWidthLines closed_polygons;
        PolylineStitcher<VariableWidthLines, ExtrusionLine, ExtrusionJunction>::stitch(wall_lines, stitched_polylines, closed_polygons, stitch_distance);
#ifdef ARACHNE_STITCH_PATCH_DEBUG
        for (const ExtrusionLine& line : stitched_polylines) {
            if ( ! line.is_odd && line.polylineLength() > 3 * stitch_distance && line.size() > 3) {
                BOOST_LOG_TRIVIAL(error) << "Some even contour lines could not be closed into polygons!";
                assert(false && "Some even contour lines could not be closed into polygons!");
                BoundingBox aabb;
                for (auto line2 : wall_lines)
                    for (auto j : line2)
                        aabb.merge(j.p);
                {
                    static int iRun = 0;
                    SVG svg(debug_out_path("contours_before.svg-%d.png", iRun), aabb);
                    std::array<const char *, 8> colors    = {"gray", "black", "blue", "green", "lime", "purple", "red", "yellow"};
                    size_t                      color_idx = 0;
                    for (auto& inset : toolpaths)
                        for (auto& line2 : inset) {
                            // svg.writePolyline(line2.toPolygon(), col);

                            Polygon poly = line2.toPolygon();
                            Point last = poly.front();
                            for (size_t idx = 1 ; idx < poly.size(); idx++) {
                                Point here = poly[idx];
                                svg.draw(Line(last, here), colors[color_idx]);
//                                svg.draw_text((last + here) / 2, std::to_string(line2.junctions[idx].region_id).c_str(), "black");
                                last = here;
                            }
                            svg.draw(poly[0], colors[color_idx]);
                            // svg.nextLayer();
                            // svg.writePoints(poly, true, 0.1);
                            // svg.nextLayer();
                            color_idx = (color_idx + 1) % colors.size();
                        }
                }
                {
                    static int iRun = 0;
                    SVG svg(debug_out_path("contours-%d.svg", iRun), aabb);
                    for (auto& inset : toolpaths)
                        for (auto& line2 : inset)
                            svg.draw_outline(line2.toPolygon(), "gray");
                    for (auto& line2 : stitched_polylines) {
                        const char *col = line2.is_odd ? "gray" : "red";
                        if ( ! line2.is_odd)
                            std::cerr << "Non-closed even wall of size: " << line2.size()  << " at " << line2.front().p << "\n";
                        if ( ! line2.is_odd)
                            svg.draw(line2.front().p);
                        Polygon poly = line2.toPolygon();
                        Point last = poly.front();
                        for (size_t idx = 1 ; idx < poly.size(); idx++)
                        {
                            Point here = poly[idx];
                            svg.draw(Line(last, here), col);
                            last = here;
                        }
                    }
                    for (auto line2 : closed_polygons)
                        svg.draw(line2.toPolygon());
                }
            }
        }
#endif // ARACHNE_STITCH_PATCH_DEBUG
        wall_lines = stitched_polylines; // replace input toolpaths with stitched polylines

        for (ExtrusionLine& wall_polygon : closed_polygons)
        {
            if (wall_polygon.junctions.empty())
            {
                continue;
            }

            // PolylineStitcher, in some cases, produced closed extrusion (polygons),
            // but the endpoints differ by a small distance. So we reconnect them.
            // FIXME Lukas H.: Investigate more deeply why it is happening.
            if (wall_polygon.junctions.front().p != wall_polygon.junctions.back().p &&
                (wall_polygon.junctions.back().p - wall_polygon.junctions.front().p).cast<double>().norm() < stitch_distance) {
                wall_polygon.junctions.emplace_back(wall_polygon.junctions.front());
            }
            wall_polygon.is_closed = true;
            wall_lines.emplace_back(std::move(wall_polygon)); // add stitched polygons to result
        }
#ifdef DEBUG
        for (ExtrusionLine& line : wall_lines)
        {
            assert(line.inset_idx == wall_idx);
        }
#endif // DEBUG
    }
}

template<typename T> bool shorterThan(const T &shape, const coord_t check_length)
{
    const auto *p0     = &shape.back();
    int64_t     length = 0;
    for (const auto &p1 : shape) {
        length += (*p0 - p1).template cast<int64_t>().norm();
        if (length >= check_length)
            return false;
        p0 = &p1;
    }
    return true;
}

void WallToolPaths::removeSmallLines(std::vector<VariableWidthLines> &toolpaths)
{
    for (VariableWidthLines &inset : toolpaths) {
        for (size_t line_idx = 0; line_idx < inset.size(); line_idx++) {
            ExtrusionLine &line      = inset[line_idx];
            coord_t        min_width = std::numeric_limits<coord_t>::max();
            for (const ExtrusionJunction &j : line)
                min_width = std::min(min_width, j.w);
            // Only use min_length_factor for non-topmost, to prevent top gaps. Otherwise use default value.
            if (line.is_odd && !line.is_closed && shorterThan(line, m_params.is_top_or_bottom_layer ? (min_width / 2) : (min_width * m_params.min_length_factor))) { // remove line
                line = std::move(inset.back());
                inset.erase(--inset.end());
                line_idx--; // reconsider the current position
            }
        }
    }
}

void WallToolPaths::simplifyToolPaths(std::vector<VariableWidthLines>& toolpaths, const WallToolPathsParams& params)
{
    const int64_t maximum_resolution = params.wall_maximum_resolution;
    const int64_t maximum_deviation  = params.wall_maximum_deviation;

    const int64_t smallest_line_segment_squared  = maximum_resolution * maximum_resolution;
    const int64_t allowed_error_distance_squared = maximum_deviation * maximum_deviation;

    const int64_t maximum_extrusion_area_deviation = Slic3r::Arachne::meshfix_maximum_extrusion_area_deviation(); // unit: μm²

    for (VariableWidthLines& lines : toolpaths)
    {
        for (ExtrusionLine& line : lines)
        {
            line.simplify(smallest_line_segment_squared, allowed_error_distance_squared, maximum_extrusion_area_deviation);
        }
    }
}

const std::vector<VariableWidthLines> &WallToolPaths::getToolPaths()
{
    if (!toolpaths_generated)
        return generate();
    return toolpaths;
}

void WallToolPaths::separateOutInnerContour()
{
    //We'll remove all 0-width paths from the original toolpaths and store them separately as polygons.
    std::vector<VariableWidthLines> actual_toolpaths;
    actual_toolpaths.reserve(toolpaths.size()); //A bit too much, but the correct order of magnitude.
    std::vector<VariableWidthLines> contour_paths;
    contour_paths.reserve(toolpaths.size() / inset_count);
    inner_contour.clear();
    for (const VariableWidthLines &inset : toolpaths) {
        if (inset.empty())
            continue;
        bool is_contour = false;
        for (const ExtrusionLine &line : inset) {
            for (const ExtrusionJunction &j : line) {
                if (j.w == 0)
                    is_contour = true;
                else
                    is_contour = false;
                break;
            }
        }

        if (is_contour) {
#ifdef DEBUG
            for (const ExtrusionLine &line : inset)
                for (const ExtrusionJunction &j : line)
                    assert(j.w == 0);
#endif // DEBUG
            for (const ExtrusionLine &line : inset) {
                if (line.is_odd)
                    continue;            // odd lines don't contribute to the contour
                else if (line.is_closed) // sometimes an very small even polygonal wall is not stitched into a polygon
                    inner_contour.emplace_back(line.toPolygon());
            }
        } else {
            actual_toolpaths.emplace_back(inset);
        }
    }
    if (!actual_toolpaths.empty())
        toolpaths = std::move(actual_toolpaths); // Filtered out the 0-width paths.
    else
        toolpaths.clear();

    //The output walls from the skeletal trapezoidation have no known winding order, especially if they are joined together from polylines.
    //They can be in any direction, clockwise or counter-clockwise, regardless of whether the shapes are positive or negative.
    //To get a correct shape, we need to make the outside contour positive and any holes inside negative.
    //This can be done by applying the even-odd rule to the shape. This rule is not sensitive to the winding order of the polygon.
    //The even-odd rule would be incorrect if the polygon self-intersects, but that should never be generated by the skeletal trapezoidation.
    inner_contour = union_(inner_contour, pftEvenOdd);
}

const Polygons& WallToolPaths::getInnerContour()
{
    if (!toolpaths_generated && inset_count > 0)
    {
        generate();
    }
    else if(inset_count == 0)
    {
        return outline;
    }
    return inner_contour;
}

bool WallToolPaths::removeEmptyToolPaths(std::vector<VariableWidthLines> &toolpaths)
{
    toolpaths.erase(std::remove_if(toolpaths.begin(), toolpaths.end(), [](const VariableWidthLines& lines)
                                   {
                                       return lines.empty();
                                   }), toolpaths.end());
    return toolpaths.empty();
}

/*!
     * Get the order constraints of the insets when printing walls per region / hole.
     * Each returned pair consists of adjacent wall lines where the left has an inset_idx one lower than the right.
     *
     * Odd walls follow the same order as even ones, so that the centre line of a narrow wall is an inner wall
     * like any other when printing inner to outer.
     *
     * \param outer_to_inner Whether the wall polygons with a lower inset_idx should go before those with a higher one.
 */
WallToolPaths::ExtrusionLineSet WallToolPaths::getRegionOrder(const std::vector<ExtrusionLine *> &input, const bool outer_to_inner)
{
    ExtrusionLineSet order_requirements;
    // We build a grid where we map toolpath vertex locations to toolpaths,
    // so that we can easily find which two toolpaths are next to each other,
    // which is the requirement for there to be an order constraint.
    //
    // We use a PointGrid rather than a LineGrid to save on computation time.
    // In very rare cases two insets might lie next to each other without having neighboring vertices, e.g.
    //  \            .
    //   |  /        .
    //   | /         .
    //   ||          .
    //   | \         .
    //   |  \        .
    //  /            .
    // However, because of how Arachne works this will likely never be the case for two consecutive insets.
    // On the other hand one could imagine that two consecutive insets of a very large circle
    // could be simplify()ed such that the remaining vertices of the two insets don't align.
    // In those cases the order requirement is not captured,
    // which means that the PathOrderOptimizer *might* result in a violation of the user set path order.
    // This problem is expected to be not so severe and happen very sparsely.

    coord_t max_line_w = 0u;
    for (const ExtrusionLine *line : input) // compute max_line_w
        for (const ExtrusionJunction &junction : *line)
            max_line_w = std::max(max_line_w, junction.w);
    if (max_line_w == 0u)
        return order_requirements;

    struct LineLoc
    {
        ExtrusionJunction    j;
        const ExtrusionLine *line;
    };
    struct Locator
    {
        Point operator()(const LineLoc &elem) { return elem.j.p; }
    };

    // How much farther two verts may be apart due to corners.
    // This distance must be smaller than 2, because otherwise
    // we could create an order requirement between e.g.
    // wall 2 of one region and wall 3 of another region,
    // while another wall 3 of the first region would lie in between those two walls.
    // However, higher values are better against the limitations of using a PointGrid rather than a LineGrid.
    constexpr float diagonal_extension = 1.9f;
    const auto      searching_radius   = coord_t(max_line_w * diagonal_extension);
    using GridT                        = SparsePointGrid<LineLoc, Locator>;
    GridT grid(searching_radius);

    for (const ExtrusionLine *line : input)
        for (const ExtrusionJunction &junction : *line) grid.insert(LineLoc{junction, line});
    for (const std::pair<const SquareGrid::GridPoint, LineLoc> &pair : grid) {
        const LineLoc       &lineloc_here = pair.second;
        const ExtrusionLine *here         = lineloc_here.line;
        Point                loc_here     = pair.second.j.p;
        std::vector<LineLoc> nearby_verts = grid.getNearby(loc_here, searching_radius);
        for (const LineLoc &lineloc_nearby : nearby_verts) {
            const ExtrusionLine *nearby = lineloc_nearby.line;
            if (nearby == here)
                continue;
            if (nearby->inset_idx == here->inset_idx)
                continue;
            if (nearby->inset_idx > here->inset_idx + 1)
                continue; // not directly adjacent
            if (here->inset_idx > nearby->inset_idx + 1)
                continue; // not directly adjacent
            if (!shorter_then(loc_here - lineloc_nearby.j.p, (lineloc_here.j.w + lineloc_nearby.j.w) / 2 * diagonal_extension))
                continue; // points are too far away from each other
            if ((nearby->inset_idx < here->inset_idx) == outer_to_inner) {
                order_requirements.emplace(std::make_pair(nearby, here));
            } else {
                assert((nearby->inset_idx > here->inset_idx) == outer_to_inner);
                order_requirements.emplace(std::make_pair(here, nearby));
            }
        }
    }
    return order_requirements;
}


// ORCA: nominal walls, see WallToolPaths.hpp.

// Width at each point of a thick polyline: the two segment ends meeting there, averaged.
static std::vector<double> point_widths(const ThickPolyline &tp)
{
    const size_t        n = tp.points.size();
    std::vector<double> w(n);
    for (size_t i = 0; i < n; ++ i)
        w[i] = i == 0 ? tp.width[0] : i == n - 1 ? tp.width[2 * (n - 2) + 1] : 0.5 * (tp.width[2 * (i - 1) + 1] + tp.width[2 * i]);
    return w;
}

// Where the width w of a center line passes cap: a wall of width opt on its left, at (w - wall width) / 2 from the old
// centerline, and a filler on its right, at -wall width / 2, which together span exactly the old line's [-w/2, w/2].
// Over ramp from where the width crosses cap the wall narrows from w to opt as it slides left and the filler grows
// from nothing, so the wall never jumps sideways. Where the filler would be narrower than min_filler it is left out.
static void split_center_line(const ThickPolyline &tp, double opt, double cap, double ramp, double min_filler,
                              ThickPolylines &walls, ThickPolylines &fillers)
{
    struct V { Vec2d p; double w; };
    std::vector<V> v;
    {
        const std::vector<double> w = point_widths(tp);
        for (size_t i = 0; i < tp.points.size(); ++ i) {
            const Vec2d p = tp.points[i].cast<double>();
            if (! v.empty() && (v.back().p - p).norm() < 1.)
                continue;
            v.push_back({ p, w[i] });
        }
    }
    const bool closed = tp.points.size() > 2 && tp.points.front() == tp.points.back();
    if (v.size() < 2) {
        walls.push_back(tp);
        return;
    }
    if (closed && (v.front().p - v.back().p).norm() > 1.)
        v.push_back(v.front());
    bool any = false, all = true;
    for (const V &x : v) {
        any |= x.w > cap;
        all &= x.w > cap;
    }
    if (! any) {
        walls.push_back(tp);
        return;
    }
    if (closed && ! all) {
        // Start a closed line where it is not split, so that no split runs across its start.
        v.pop_back();
        size_t k = 0;
        while (v[k].w > cap)
            ++ k;
        std::rotate(v.begin(), v.begin() + k, v.end());
        v.push_back(v.front());
    }
    // The points where the width crosses cap.
    std::vector<V> r{ v.front() };
    for (size_t i = 1; i < v.size(); ++ i) {
        const V &a = v[i - 1], &b = v[i];
        if ((a.w > cap) != (b.w > cap)) {
            const double t = (cap - a.w) / (b.w - a.w);
            if (t > 1e-6 && t < 1. - 1e-6)
                r.push_back({ a.p + (b.p - a.p) * t, cap });
        }
        r.push_back(b);
    }
    // Runs over cap, as first and last point, the crossings included.
    auto runs_of = [cap](const std::vector<V> &pts) {
        std::vector<std::pair<size_t, size_t>> runs;
        for (size_t i = 0; i < pts.size(); ++ i) {
            if (pts[i].w <= cap + 0.5)
                continue;
            size_t a = i, b = i;
            while (a > 0 && pts[a - 1].w > cap - 0.5)
                -- a;
            while (b + 1 < pts.size() && pts[b + 1].w > cap - 0.5)
                ++ b;
            runs.emplace_back(a, b);
            i = b;
        }
        return runs;
    };
    auto arc_lengths = [](const std::vector<V> &pts) {
        std::vector<double> s(pts.size(), 0.);
        for (size_t i = 1; i < pts.size(); ++ i)
            s[i] = s[i - 1] + (pts[i].p - pts[i - 1].p).norm();
        return s;
    };
    // A ramp only starts or ends at a crossing; a run that reaches the end of an open line stays split to its end.
    auto is_crossing = [cap](const std::vector<V> &pts, size_t i) { return pts[i].w <= cap + 0.5; };
    {
        // Points where the ramps end, so that they end exactly a ramp length in.
        const std::vector<double> s = arc_lengths(r);
        std::vector<double>       cuts;
        if (! (all && closed))
            for (auto [a, b] : runs_of(r)) {
                if (is_crossing(r, a) && s[a] + ramp < s[b])
                    cuts.push_back(s[a] + ramp);
                if (is_crossing(r, b) && s[b] - ramp > s[a])
                    cuts.push_back(s[b] - ramp);
            }
        std::sort(cuts.begin(), cuts.end());
        std::vector<V> rr{ r.front() };
        size_t         c = 0;
        for (size_t i = 1; i < r.size(); ++ i) {
            for (; c < cuts.size() && cuts[c] < s[i] - 1.; ++ c)
                if (cuts[c] > s[i - 1] + 1.) {
                    const double t = (cuts[c] - s[i - 1]) / (s[i] - s[i - 1]);
                    rr.push_back({ r[i - 1].p + (r[i].p - r[i - 1].p) * t, r[i - 1].w + (r[i].w - r[i - 1].w) * t });
                }
            rr.push_back(r[i]);
        }
        r.swap(rr);
    }
    const size_t              m = r.size();
    const std::vector<double> s = arc_lengths(r);
    // Offset direction at each point, to the left, mitred.
    std::vector<Vec2d> sn(m - 1), vn(m);
    for (size_t i = 0; i + 1 < m; ++ i) {
        const Vec2d d = r[i + 1].p - r[i].p;
        sn[i] = Vec2d(- d.y(), d.x()) / d.norm();
    }
    for (size_t i = 0; i < m; ++ i) {
        const Vec2d  n0  = i > 0 ? sn[i - 1] : (closed ? sn[m - 2] : sn[0]);
        const Vec2d  n1  = i + 1 < m ? sn[i] : (closed ? sn[0] : sn[m - 2]);
        const Vec2d  sum = n0 + n1;
        const double len = sum.norm();
        vn[i] = len < 1e-9 ? n1 : Vec2d(sum / len / std::max((sum / len).dot(n1), 0.5));
    }
    // How far into its split each point is: 0 outside or at a crossing, 1 a ramp length or more inside.
    const auto          runs = runs_of(r);
    std::vector<double> beta(m, 0.);
    for (auto [a, b] : runs)
        for (size_t i = a; i <= b; ++ i) {
            const double da = is_crossing(r, a) ? s[i] - s[a] : std::numeric_limits<double>::max();
            const double db = is_crossing(r, b) ? s[b] - s[i] : std::numeric_limits<double>::max();
            beta[i] = (all && closed) ? 1. : std::clamp(std::min(da, db) / ramp, 0., 1.);
        }
    // The wall narrows to opt over the ramp, and is never wider than cap: where the center line widens faster than the
    // ramp narrows the wall, the filler takes the difference.
    auto wall_width = [&](size_t i) {
        double ww = std::min(cap, r[i].w - beta[i] * std::max(r[i].w - opt, 0.));
        // Past the cap the filler is never narrower than min_filler: the wall gives up the difference.
        if (r[i].w > cap - 0.5 && r[i].w - ww < min_filler)
            ww = r[i].w - min_filler;
        return ww;
    };
    ThickPolyline wall;
    double        last_w = 0.;
    for (size_t i = 0; i < m; ++ i) {
        const double ww = wall_width(i);
        const Vec2d  p  = r[i].p + vn[i] * (0.5 * (r[i].w - ww));
        if (! wall.points.empty() && (p - wall.points.back().cast<double>()).norm() < 1.)
            continue;
        if (! wall.points.empty()) {
            wall.width.push_back(last_w);
            wall.width.push_back(ww);
        }
        wall.points.push_back(p.cast<coord_t>());
        last_w = ww;
    }
    if (wall.points.size() >= 2)
        walls.push_back(std::move(wall));
    for (auto [a, b] : runs) {
        ThickPolyline f;
        double        fw_prev = 0.;
        for (size_t i = a; i <= b; ++ i) {
            const double ww = wall_width(i);
            const double fw = r[i].w - ww;
            if (fw < min_filler) {
                // Narrower than a filler may be: left out, ending this filler or not yet starting one.
                if (f.points.size() >= 2)
                    fillers.push_back(std::move(f));
                f = ThickPolyline();
                continue;
            }
            const Vec2d  p  = r[i].p - vn[i] * (0.5 * ww);
            if (! f.points.empty()) {
                if ((p - f.points.back().cast<double>()).norm() < 1.)
                    continue;
                f.width.push_back(fw_prev);
                f.width.push_back(fw);
            }
            f.points.push_back(p.cast<coord_t>());
            fw_prev = fw;
        }
        if (f.points.size() >= 2)
            fillers.push_back(std::move(f));
    }
}

static ExtrusionLine to_extrusion_line(const ThickPolyline &tp, size_t inset_idx, bool is_odd, size_t perimeter_index)
{
    ExtrusionLine             line(inset_idx, is_odd);
    const std::vector<double> w = point_widths(tp);
    for (size_t i = 0; i < tp.points.size(); ++ i)
        line.junctions.emplace_back(tp.points[i], coord_t(w[i]), perimeter_index);
    line.is_closed = tp.points.size() > 2 && tp.points.front() == tp.points.back();
    return line;
}

// The runs of a thick polyline at least min_width wide, each a polyline of its own.
static ThickPolylines thick_runs_at_least(const ThickPolyline &tp, double min_width)
{
    ThickPolylines            out;
    const std::vector<double> w = point_widths(tp);
    ThickPolyline             run;
    for (size_t i = 0; i < tp.points.size(); ++ i) {
        if (w[i] < min_width) {
            if (run.points.size() >= 2)
                out.push_back(std::move(run));
            run = ThickPolyline();
            continue;
        }
        if (! run.points.empty()) {
            run.width.push_back(w[i - 1]);
            run.width.push_back(w[i]);
        }
        run.points.push_back(tp.points[i]);
    }
    if (run.points.size() >= 2)
        out.push_back(std::move(run));
    return out;
}

// The runs of a line at least min_width wide, each a line of its own.
static VariableWidthLines runs_at_least(const ExtrusionLine &line, coord_t min_width)
{
    VariableWidthLines out;
    ExtrusionLine      run(line.inset_idx, line.is_odd);
    for (const ExtrusionJunction &j : line.junctions) {
        if (j.w >= min_width) {
            run.junctions.push_back(j);
            continue;
        }
        if (run.junctions.size() >= 2)
            out.push_back(run);
        run.junctions.clear();
    }
    if (run.junctions.size() >= 2)
        out.push_back(run);
    for (ExtrusionLine &l : out)
        l.is_closed = false;
    return out;
}

ThickPolylines split_center_walls(std::vector<VariableWidthLines> &toolpaths, coord_t nominal_width, coord_t max_width, coord_t ramp_length,
                                  coord_t min_filler_width, coord_t split_wall_width, bool narrow_center_lines_are_fillers)
{
    ThickPolylines fillers;
    const double   opt        = double(nominal_width);
    const double   cap        = std::max(double(max_width), opt);
    const double   min_filler = std::max(double(min_filler_width), scaled<double>(0.02));
    auto is_center = [](const ExtrusionLine &l) { return l.is_odd && l.inset_idx > 0 && l.junctions.size() >= 2; };
    auto key       = [](const Point &p) { return std::make_pair(p.x(), p.y()); };

    // Where center lines end, taken before their thin ends are dropped below: lines meet there.
    Points center_end_points;
    for (const VariableWidthLines &lines : toolpaths)
        for (const ExtrusionLine &l : lines)
            if (is_center(l)) {
                center_end_points.push_back(l.junctions.front().p);
                center_end_points.push_back(l.junctions.back().p);
            }
    const double match_dist = 0.25 * opt;
    auto at_center_end = [&](const Point &p) {
        for (const Point &q : center_end_points)
            if ((q - p).cast<double>().norm() < match_dist)
                return true;
        return false;
    };

    // With the slack all in the center line, Arachne's transitions thin a center line towards nothing, even below:
    // the parts narrower than a filler may be are dropped.
    for (VariableWidthLines &lines : toolpaths) {
        bool thin = false;
        for (const ExtrusionLine &l : lines)
            if (is_center(l))
                for (const ExtrusionJunction &j : l.junctions)
                    thin |= j.w < coord_t(min_filler);
        if (! thin)
            continue;
        VariableWidthLines out;
        for (ExtrusionLine &l : lines) {
            bool l_thin = false;
            if (is_center(l))
                for (const ExtrusionJunction &j : l.junctions)
                    l_thin |= j.w < coord_t(min_filler);
            if (l_thin)
                append(out, runs_at_least(l, coord_t(min_filler)));
            else
                out.push_back(std::move(l));
        }
        lines = std::move(out);
    }

    // The lines that take over where a center line ends start there too, flaring from its full width: they start a
    // junction later, in their own lanes. A closed line that started and ended there becomes open.
    std::vector<std::pair<coord_t, coord_t>> arm_ends;
    for (VariableWidthLines &lines : toolpaths)
        for (ExtrusionLine &l : lines) {
            if (l.is_odd || l.inset_idx == 0 || l.junctions.size() < 3)
                continue;
            const bool trim_front = at_center_end(l.junctions.front().p);
            const bool trim_back  = at_center_end(l.junctions.back().p);
            if (trim_back && l.junctions.size() >= 3) {
                l.junctions.pop_back();
                l.is_closed = false;
                arm_ends.push_back(key(l.junctions.back().p));
            }
            if (trim_front && l.junctions.size() >= 3) {
                l.junctions.erase(l.junctions.begin());
                l.is_closed = false;
                arm_ends.push_back(key(l.junctions.front().p));
            }
            // The flare runs on over the transition: near a trimmed end, a line is no wider than it settles to.
            for (int end = 0; end < 2; ++ end) {
                if (! (end ? trim_back : trim_front))
                    continue;
                const size_t n = l.junctions.size();
                auto at = [&](size_t i) -> ExtrusionJunction & { return l.junctions[end ? n - 1 - i : i]; };
                double dist = 0.;
                size_t i    = 1;
                for (; i < n && dist < 2. * opt; ++ i)
                    dist += (at(i).p - at(i - 1).p).cast<double>().norm();
                const coord_t settled = at(std::min(i, n - 1)).w;
                for (size_t k = 0; k < std::min(i, n); ++ k)
                    at(k).w = std::min(at(k).w, settled);
            }
        }
    // A pair can also run into a center line's end and straight back out, as one line folded there: the fold takes the
    // center line's full width and flares out on both sides. It is capped the same way, on each side.
    for (VariableWidthLines &lines : toolpaths)
        for (ExtrusionLine &l : lines) {
            const size_t n = l.junctions.size();
            if (l.is_odd || l.inset_idx == 0 || n < 3)
                continue;
            for (size_t f = 1; f + 1 < n; ++ f) {
                if (! at_center_end(l.junctions[f].p))
                    continue;
                coord_t settled[2];
                size_t  reach[2];
                for (int dir = 0; dir < 2; ++ dir) {
                    double dist = 0.;
                    size_t i    = f;
                    while (dist < 2. * opt && (dir ? i + 1 < n : i > 0)) {
                        const size_t next = dir ? i + 1 : i - 1;
                        dist += (l.junctions[next].p - l.junctions[i].p).cast<double>().norm();
                        i = next;
                    }
                    settled[dir] = l.junctions[i].w;
                    reach[dir]   = i;
                }
                for (size_t i = reach[0]; i < f; ++ i)
                    l.junctions[i].w = std::min(l.junctions[i].w, settled[0]);
                for (size_t i = f + 1; i <= reach[1]; ++ i)
                    l.junctions[i].w = std::min(l.junctions[i].w, settled[1]);
                l.junctions[f].w = std::min(l.junctions[f].w, std::max(settled[0], settled[1]));
            }
        }
    std::sort(arm_ends.begin(), arm_ends.end());

    // Split the center lines. Walls out of a split are marked, to join them below.
    std::vector<std::vector<bool>> split_wall(toolpaths.size());
    for (size_t k = 0; k < toolpaths.size(); ++ k) {
        VariableWidthLines out;
        for (ExtrusionLine &l : toolpaths[k]) {
            if (! is_center(l)) {
                out.push_back(std::move(l));
                split_wall[k].push_back(false);
                continue;
            }
            ThickPolylines walls, f;
            split_center_line(to_thick_polyline(l), std::min(double(split_wall_width), cap), cap, double(ramp_length), min_filler, walls, f);
            append(fillers, std::move(f));
            for (ThickPolyline &w : walls) {
                bool narrow = narrow_center_lines_are_fillers;
                for (const coordf_t width : w.width)
                    narrow &= width < opt;
                if (narrow) {
                    fillers.push_back(std::move(w));
                    continue;
                }
                out.push_back(to_extrusion_line(w, l.inset_idx, true, l.junctions.front().perimeter_index));
                split_wall[k].push_back(true);
            }
        }
        toolpaths[k] = std::move(out);
    }

    // A wall out of a split continues, where its center line ended at a transition, into the line that takes over in
    // its lane: the nearest trimmed end, within a line width.
    const double join_dist = 1.5 * opt;
    for (size_t k = 0; k < toolpaths.size(); ++ k) {
        VariableWidthLines &lines = toolpaths[k];
        for (size_t i = 0; i < lines.size(); ++ i) {
            if (! split_wall[k][i] || lines[i].junctions.size() < 2 || lines[i].is_closed)
                continue;
            for (int end = 0; end < 2; ++ end) {
                ExtrusionLine &w  = lines[i];
                const Point    we = end ? w.junctions.back().p : w.junctions.front().p;
                double         best = join_dist;
                size_t         bj   = 0;
                bool           bfront = true;
                for (size_t j = 0; j < lines.size(); ++ j) {
                    if (j == i || lines[j].is_odd || lines[j].junctions.size() < 2)
                        continue;
                    for (int ae = 0; ae < 2; ++ ae) {
                        const Point &p = ae ? lines[j].junctions.back().p : lines[j].junctions.front().p;
                        if (! std::binary_search(arm_ends.begin(), arm_ends.end(), key(p)))
                            continue;
                        const double d = (p - we).cast<double>().norm();
                        if (d < best) {
                            best   = d;
                            bj     = j;
                            bfront = ae == 0;
                        }
                    }
                }
                if (best >= join_dist)
                    continue;
                ExtrusionLine &a = lines[bj];
                if (! bfront)
                    std::reverse(a.junctions.begin(), a.junctions.end());
                if (! end)
                    std::reverse(w.junctions.begin(), w.junctions.end());
                w.junctions.insert(w.junctions.end(), a.junctions.begin(), a.junctions.end());
                w.is_odd = false;
                a.junctions.clear();
                if (! end)
                    std::reverse(w.junctions.begin(), w.junctions.end());
            }
        }
        lines.erase(std::remove_if(lines.begin(), lines.end(), [](const ExtrusionLine &l) { return l.junctions.size() < 2; }), lines.end());
    }
    // No wall inside the outer walls wider than the maximum: where Arachne still flares one past it, the flare lies
    // over lines that cover that place already.
    for (VariableWidthLines &lines : toolpaths)
        for (ExtrusionLine &l : lines)
            if (l.inset_idx > 0)
                for (ExtrusionJunction &j : l.junctions)
                    j.w = std::min(j.w, coord_t(cap));
    // Nor any filler narrower than the minimum.
    ThickPolylines out;
    for (const ThickPolyline &f : fillers)
        append(out, thick_runs_at_least(f, min_filler));
    return out;
}

ThickPolylines narrow_infill_strips(ExPolygons &infill_contour, coord_t min_infill_width, double min_width, double max_width)
{
    // Only regions narrow all over: the narrow corners of a region with room for infill are left to the infill.
    ThickPolylines out;
    ExPolygons     keep;
    for (ExPolygon &region : infill_contour) {
        if (offset_ex(region, - float(min_infill_width / 2.)).empty()) {
            // Simplified first, as for gap fill: the inner contour's fine detail breaks the medial axis into pieces.
            ExPolygon simplified = region;
            simplified.douglas_peucker(scaled<double>(0.0125));
            simplified.medial_axis(min_width, max_width, &out);
        }
        else
            keep.push_back(std::move(region));
    }
    infill_contour = std::move(keep);
    ThickPolylines trimmed;
    for (const ThickPolyline &t : out)
        append(trimmed, thick_runs_at_least(t, min_width));
    return trimmed;
}


// The area a thick polyline's lines cover: a quad per segment and a disc at each point, so that lines meeting at an
// angle or ending leave no hairline open.
static Polygons thick_polyline_cover(const ThickPolyline &tp)
{
    Polygons                  out;
    const std::vector<double> w = point_widths(tp);
    for (size_t i = 0; i + 1 < tp.points.size(); ++ i) {
        const Vec2d  a = tp.points[i].cast<double>(), b = tp.points[i + 1].cast<double>(), d = b - a;
        const double len = d.norm();
        if (len < 1.)
            continue;
        const Vec2d n(- d.y() / len, d.x() / len);
        Polygon quad({ (a + n * (0.5 * w[i])).cast<coord_t>(), (b + n * (0.5 * w[i + 1])).cast<coord_t>(),
                       (b - n * (0.5 * w[i + 1])).cast<coord_t>(), (a - n * (0.5 * w[i])).cast<coord_t>() });
        quad.make_counter_clockwise();
        out.emplace_back(std::move(quad));
    }
    for (size_t i = 0; i < tp.points.size(); ++ i)
        if (w[i] > 0.) {
            Polygon disc;
            for (int k = 0; k < 12; ++ k) {
                const double angle = 2. * PI * k / 12.;
                disc.points.emplace_back((tp.points[i].cast<double>() + Vec2d(cos(angle), sin(angle)) * (0.5 * w[i])).cast<coord_t>());
            }
            out.emplace_back(std::move(disc));
        }
    return out;
}

ThickPolylines fill_uncovered(const Polygons &outline, const std::vector<VariableWidthLines> &toolpaths, const ThickPolylines &fillers,
                              const ExPolygons &infill_contour, double min_width, double max_width)
{
    Polygons covered = to_polygons(infill_contour);
    for (const VariableWidthLines &lines : toolpaths)
        for (const ExtrusionLine &l : lines)
            if (l.junctions.size() >= 2)
                append(covered, thick_polyline_cover(to_thick_polyline(l)));
    for (const ThickPolyline &f : fillers)
        append(covered, thick_polyline_cover(f));
    // Grown a little, so that the seams between lines are not taken for unfilled; the opening drops the hairlines left.
    ExPolygons     unfilled = opening_ex(diff_ex(outline, offset(covered, float(scaled<double>(0.01)))), float(min_width / 2.));
    ThickPolylines lines;
    for (ExPolygon &area : unfilled) {
        area.douglas_peucker(scaled<double>(0.0125));
        // Not only where the sides are near parallel: what is left unfilled is often a wedge.
        area.medial_axis(min_width, max_width, &lines, false);
    }
    ThickPolylines out;
    for (const ThickPolyline &t : lines)
        append(out, thick_runs_at_least(t, min_width));
    return out;
}

} // namespace Slic3r::Arachne
