// Copyright (c) 2020 Ultimaker B.V.
// CuraEngine is released under the terms of the AGPLv3 or higher.

#ifndef CURAENGINE_WALLTOOLPATHS_H
#define CURAENGINE_WALLTOOLPATHS_H

#include <cstddef>
#include <boost/container_hash/hash.hpp>
#include <memory>
#include <ankerl/unordered_dense.h>
#include <vector>
#include <utility>

#include "libslic3r/libslic3r.h"
#include "libslic3r/Point.hpp"
#include "utils/ExtrusionLine.hpp"
#include "../Polygon.hpp"
#include "../PrintConfig.hpp"

namespace Slic3r::Arachne
{

constexpr bool    fill_outline_gaps = true;
inline coord_t    meshfix_maximum_resolution() { return scaled<coord_t>(0.5); }
inline coord_t    meshfix_maximum_deviation() { return scaled<coord_t>(0.025); }
inline coord_t    meshfix_maximum_extrusion_area_deviation() { return scaled<coord_t>(2.); }

class WallToolPathsParams
{
public:
    float   min_bead_width                   = 0.f;
    float   min_feature_size                 = 0.f;
    float   min_length_factor                = 0.5f;
    float   wall_transition_length           = 0.f;
    float   wall_transition_angle            = 10.f;
    float   wall_transition_filter_deviation = 0.f;
    int     wall_distribution_count          = 1;
    bool    is_top_or_bottom_layer           = false;
    // ORCA: when the bead count steps up, as fractions of the bead width: from one bead to two (split) and from an
    // even count to the next odd one (add). Left at 0 they follow min_bead_width, which is also the width thin
    // features are widened to; a caller that wants the steps elsewhere without widening sets them here.
    double  wall_split_middle_threshold      = 0.;
    double  wall_add_middle_threshold        = 0.;
    // ORCA: at the bead count limit, leave extra thickness to the inner contour (infill) instead of widening the walls.
    bool    overflow_early                   = false;
    // ORCA: thresholds by inner bead count n, for the transition from n inner beads to n + 1 (see
    // BeadingStrategy::setCountThresholds()). Empty: the two thresholds above.
    std::vector<double> transition_thresholds;

    coord_t wall_maximum_resolution = meshfix_maximum_resolution();
    coord_t wall_maximum_deviation  = meshfix_maximum_deviation();
};

WallToolPathsParams make_paths_params(const int layer_id, const PrintObjectConfig &print_object_config, const PrintConfig &print_config);

// ORCA: nominal walls. Arachne is set to keep its walls at nominal width and put the slack into one center line; these
// shape that center line. split_center_walls():
//  - Where a center line ends at a bead count transition, the two lines that take over start there too, each flaring
//    from the center line's full width: they start a junction later instead, in their own lanes.
//  - A center line wider than max_width becomes a wall of split_wall_width on its left and a filler on its right, which
//    together cover exactly what it covered; the split ramps in and out over ramp_length.
//  - Such a wall continues, at a transition, into the line that takes over in its lane.
//  - With narrow_center_lines_are_fillers, center lines narrower than nominal_width all along are fillers.
//  - No wall inside the outer walls is left wider than max_width, nor any filler narrower than min_filler_width:
//    what Arachne still flares past the maximum covers what is covered already, and thinner fillers are left out.
// Fillers are taken out of the toolpaths and returned, to be printed after the walls as gap fill.
// All widths are spacings, as Arachne's are.
ThickPolylines split_center_walls(std::vector<VariableWidthLines> &toolpaths, coord_t nominal_width, coord_t max_width, coord_t ramp_length,
                                  coord_t min_filler_width, coord_t split_wall_width, bool narrow_center_lines_are_fillers);

// ORCA: nominal walls leave extra thickness at the wall count limit to the inner contour, as strips too narrow for
// infill. Takes the regions of infill_contour narrower than min_infill_width all over out of it and returns lines
// filling them, as gap fill.
ThickPolylines narrow_infill_strips(ExPolygons &infill_contour, coord_t min_infill_width, double min_width, double max_width);

// ORCA: nominal walls put a single center line where Arachne would put a pair, and Arachne stops a single line short of
// the narrow end of a tapering region, so the tips of such regions are left empty. Returns gap fill for what lies inside
// outline but under none of the toolpaths, the fillers and the infill: no line narrower than min_width or wider than
// max_width.
ThickPolylines fill_uncovered(const Polygons &outline, const std::vector<VariableWidthLines> &toolpaths, const ThickPolylines &fillers,
                              const ExPolygons &infill_contour, double min_width, double max_width);

class WallToolPaths
{
public:
    /*!
     * A class that creates the toolpaths given an outline, nominal bead width and maximum amount of walls
     * \param outline An outline of the area in which the ToolPaths are to be generated
     * \param bead_width_0 The bead width of the first wall used in the generation of the toolpaths
     * \param bead_width_x The bead width of the inner walls used in the generation of the toolpaths
     * \param inset_count The maximum number of parallel extrusion lines that make up the wall
     * \param wall_0_inset How far to inset the outer wall, to make it adhere better to other walls.
     */
    WallToolPaths(const Polygons& outline, coord_t bead_width_0, coord_t bead_width_x, size_t inset_count, coord_t wall_0_inset, coordf_t layer_height, const WallToolPathsParams &params);

    /*!
     * Generates the Toolpaths
     * \return A reference to the newly create  ToolPaths
     */
    const std::vector<VariableWidthLines> &generate();

    /*!
     * Gets the toolpaths, if this called before \p generate() it will first generate the Toolpaths
     * \return a reference to the toolpaths
     */
    const std::vector<VariableWidthLines> &getToolPaths();

    /*!
     * Compute the inner contour of the walls. This contour indicates where the walled area ends and its infill begins.
     * The inside can then be filled, e.g. with skin/infill for the walls of a part, or with a pattern in the case of
     * infill with extra infill walls.
     */
    void separateOutInnerContour();

    /*!
     * Gets the inner contour of the area which is inside of the generated tool
     * paths.
     *
     * If the walls haven't been generated yet, this will lazily call the
     * \p generate() function to generate the walls with variable width.
     * The resulting polygon will snugly match the inside of the variable-width
     * walls where the walls get limited by the LimitedBeadingStrategy to a
     * maximum wall count.
     * If there are no walls, the outline will be returned.
     * \return The inner contour of the generated walls.
     */
    const Polygons& getInnerContour();

    /*!
     * Removes empty paths from the toolpaths
     * \param toolpaths the VariableWidthPaths generated with \p generate()
     * \return true if there are still paths left. If all toolpaths were removed it returns false
     */
    static bool removeEmptyToolPaths(std::vector<VariableWidthLines> &toolpaths);

    using ExtrusionLineSet = ankerl::unordered_dense::set<std::pair<const ExtrusionLine *, const ExtrusionLine *>, boost::hash<std::pair<const ExtrusionLine *, const ExtrusionLine *>>>;

    /*!
     * Get the order constraints of the insets when printing walls per region / hole.
     * Each returned pair consists of adjacent wall lines where the left has an inset_idx one lower than the right.
     *
     * Odd walls should always go after their enclosing wall polygons.
     *
     * \param outer_to_inner Whether the wall polygons with a lower inset_idx should go before those with a higher one.
     */
    static ExtrusionLineSet getRegionOrder(const std::vector<ExtrusionLine *> &input, bool outer_to_inner);

protected:
    /*!
     * Stitch the polylines together and form closed polygons.
     *
     * Works on both toolpaths and inner contours simultaneously.
     */
    static void stitchToolPaths(std::vector<VariableWidthLines> &toolpaths, coord_t bead_width_x);

    /*!
     * Remove polylines shorter than half the smallest line width along that polyline.
     */
    void removeSmallLines(std::vector<VariableWidthLines> &toolpaths);

    /*!
     * Simplifies the variable-width toolpaths by calling the simplify on every line in the toolpath using the provided
     * settings.
     * \param toolpaths The toolpaths vector to simplify
     * \param params The settings as provided by the user
     * \return
     */
    static void simplifyToolPaths(std::vector<VariableWidthLines>& toolpaths, const WallToolPathsParams& params);

private:
    const Polygons& outline; //<! A reference to the outline polygon that is the designated area
    coord_t bead_width_0; //<! The nominal or first extrusion line width with which libArachne generates its walls
    coord_t bead_width_x; //<! The subsequently extrusion line width with which libArachne generates its walls if WallToolPaths was called with the nominal_bead_width Constructor this is the same as bead_width_0
    size_t inset_count; //<! The maximum number of walls to generate
    coord_t wall_0_inset; //<! How far to inset the outer wall. Should only be applied when printing the actual walls, not extra infill/skin/support walls.
    coordf_t layer_height;
    bool print_thin_walls; //<! Whether to enable the widening beading meta-strategy for thin features
    coord_t min_feature_size; //<! The minimum size of the features that can be widened by the widening beading meta-strategy. Features thinner than that will not be printed
    coord_t min_bead_width;  //<! The minimum bead size to use when widening thin model features with the widening beading meta-strategy
    double small_area_length; //<! The length of the small features which are to be filtered out, this is squared into a surface
    coord_t wall_transition_filter_deviation; //!< The allowed line width deviation induced by filtering
    bool toolpaths_generated; //<! Are the toolpaths generated
    std::vector<VariableWidthLines> toolpaths; //<! The generated toolpaths
    Polygons inner_contour;  //<! The inner contour of the generated toolpaths
    const WallToolPathsParams m_params;
};

} // namespace Slic3r::Arachne

#endif // CURAENGINE_WALLTOOLPATHS_H
