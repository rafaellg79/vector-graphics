// Based on "Massively-Parallel Vector Graphics", Ganacim et al., Nov 2014
#ifndef SHORTCUT_TREE_H
#define SHORTCUT_TREE_H

#include <stdlib.h>
#include <math.h>

#include "vec.h"
#include "bounding_box.h"
#include "curve.h"
#include "path.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ShortcutNode{
    struct ShortcutNode *cells;
    Element *elements;
    int total_segments, size;
    BoundingBox BB;
} ShortcutNode;

typedef struct ShortcutTree {
    Shape* shapes;
    Curve* shortcuts;
    ShortcutNode * root;
    int size;
    int num_shortcuts, max_shortcuts;
    int* max_segs;
    int max_depth;
} ShortcutTree;

inline ShortcutNode* get(ShortcutNode* tree, int node_id, int max_depth) {
    if (node_id < 0) return NULL;
    else if (node_id == 0) return tree;
    int* path = (int*)malloc(max_depth * sizeof(int));
    int descendants = 0;
    int current_id = node_id;
    while(descendants < max_depth && current_id > 0){
        current_id -= 1;
        path[descendants++] = current_id & 3;
        current_id >>= 2;
    }
    ShortcutNode* node = tree;
    while (descendants && node != NULL) {
        node = node->cells + (path[--descendants]);
    }
    free(path);
    return node;
}

void build_tree(ShortcutTree* tree, real cellsize);
void build_tree_rec(ShortcutTree* tree, ShortcutNode* nodes, const int depth);
void delete_nodes(ShortcutNode* node);
void delete_tree(ShortcutTree* tree);

#ifdef __cplusplus
}
#endif

// Horizontal ray test returns
// a positive value if the point (x, y) is to the left of the segment,
// a negative value if the point (x, y) is to the right of the segment and
// 0 if the point (x, y) is over the segment.
// For this test, the segment is considered closed on both sides.
real hhit(Curve curve, real x, real y) {
    if (x < curve.BB.xmin || y <= curve.BB.ymin || y > curve.BB.ymax)
        return -1.;
    if (x > curve.BB.xmax)
        return 1.;
    real lv = curve.l.x * x + curve.l.y * y + curve.l.z;
    switch (curve.type) {
    case LINE:
        return lv;

    case QUADRATIC: {
        if (lv > 0.)
            return curve.convex ? 1. : -1.;
        real dx = x - curve.x0, dy = y - curve.y0;
        real qv = dx * (dx * curve.F[1] + dy * curve.F[2] + curve.F[0])
            + dy * (dy * curve.F[4] + curve.F[3]);
        return qv;
    }

    case CUBIC: {
        if (lv > 0.)
            return curve.convex ? 1. : -1.;
        if (curve.m.x * x + curve.m.y * y + curve.m.z > 0. ||
            curve.n.x * x + curve.n.y * y + curve.n.z > 0.)
            return (!curve.convex) ? 1. : -1.;
        real dx = x - curve.x0, dy = y - curve.y0;
        real cv = dx * (dx * (dx * curve.F[8] + dy * curve.F[7] + curve.F[2])
            + dy * (dy * curve.F[6] + curve.F[3]) + curve.F[0])
            + dy * (dy * (dy * curve.F[5] + curve.F[4]) + curve.F[1]);
        return cv;
    }

    default:
        // Unknown curve type
        return NAN;
    }
}

// Vertical ray test returns
// a positive value if the point (x, y) is below the segment,
// a negative value if the point (x, y) is above the segment and
// 0 if the point (x, y) is over the segment.
// For this test, the segment is considered closed on both sides.
real vhit(Curve curve, real x, real y) {
    if (x < curve.BB.xmin || x > curve.BB.xmax || y < curve.BB.ymin)
        return -1.;
    if (y > curve.BB.ymax)
        return 1.;
    real lv = curve.l.x * x + curve.l.y * y + curve.l.z;
    real q = curve.asc ? -1. : 1.;
    switch (curve.type) {
    case LINE:
        return lv * q;

    case QUADRATIC: {
        if (lv > 0.)
            return (curve.convex == (q > 0)) ? 1. : -1.;
        real dx = x - curve.x0, dy = y - curve.y0;
        real qv = dx * (dx * curve.F[1] + dy * curve.F[2] + curve.F[0])
            + dy * (dy * curve.F[4] + curve.F[3]);
        return qv * q;
    }

    case CUBIC: {
        if (lv > 0.)
            return (curve.convex == (q > 0)) ? 1. : -1.;
        if (curve.m.x * x + curve.m.y * y + curve.m.z > 0. ||
            curve.n.x * x + curve.n.y * y + curve.n.z > 0.)
            return (!curve.convex == (q > 0)) ? 1. : -1.;
        real dx = x - curve.x0, dy = y - curve.y0;
        real cv = dx * (dx * (dx * curve.F[8] + dy * curve.F[7] + curve.F[2])
            + dy * (dy * curve.F[6] + curve.F[3]) + curve.F[0])
            + dy * (dy * (dy * curve.F[5] + curve.F[4]) + curve.F[1]);
        return cv * q;
    }

    default:
        fprintf(stderr, "Unknown curve type: %d\n", curve.type);
        return NAN;
    }
}

// Subcell combination identifiers
// BL   Bottom-Left
// BR   Bottom-Right
// TL   Top-Left
// TR   Top-Right
// B    BL and BR
// T    TL and TR
// L    BL and TL
// R    BR and TR
// A    BL and TR
// NBL  BR, TL and TR
// NBR  BL, TL and TR
// NTL  BL, BR and TR
// NTR  BL, BR and TL
enum subcells {
    BL = 1, BR = 2, TL = 4, TR = 8,
    B = 3, T = 12, L = 5, R = 10, A = 9,
    NBL = 14, NBR = 13, NTL = 11, NTR = 7
};

// Determines on which subcells a segment is contained
// Assumes the segment is either a monotonic segment
// contained on the cell or a shortcut
char classify_seg(Curve seg, BoundingBox BB) {
    real ymid = (BB.ymax + BB.ymin) * .5;
    real xmid = (BB.xmax + BB.xmin) * .5;

    if (seg.BB.xmax <= BB.xmin) // Outside the cell
        // Shortcut
        return 0;

    // Bounding box intersection with halves
    bool l = seg.BB.xmin < xmid || ((seg.BB.xmin == xmid) && !seg.asc);
    bool r = seg.BB.xmax > xmid;
    bool b = seg.BB.ymin < ymid;
    bool t = seg.BB.ymax > ymid;

    // Horizontal segment on middle boundary
    if (!(b || t))
        return 0;

    // If bounding box overlaps a single subcell,
    // return that subcell
    if (!r && !t)
        return BL;
    if (!l && !t)
        return BR;
    if (!r && !b)
        return TL;
    if (!l && !b)
        return TR;

    // If bounding box overlaps 2 subcells,
    // differentiate cases using single ray test
    real v;
    if (!r) {
        // left
        v = hhit(seg, BB.xmin, ymid);
        if (v < 0.) return L;
        else return seg.asc ? TL : BL;
    }
    if (!l) {
        // right
        v = hhit(seg, BB.xmax, ymid);
        if (v > 0.) return R;
        else if (v == 0.)
            return seg.asc ? BR : R;
        else
            return seg.asc ? BR : TR;
    }
    if (!t) {
        // bottom
        v = vhit(seg, xmid, BB.ymin);
        if (v < 0.) return B;
        else return seg.asc ? BR : BL;
    }
    if (!b) {
        // top
        v = vhit(seg, xmid, BB.ymax);
        if (v > 0.) return T;
        else if (v == 0.)
            return seg.asc ? TL : T;
        else
            return seg.asc ? TL : TR;
    }

    // If bounding box overlaps all subcells,
    // differentiate cases using three ray tests
    real hc = hhit(seg, xmid, ymid);
    if (seg.asc) {
        if (hc < 0.) {
            real hr = hhit(seg, BB.xmax, ymid);
            real vb = vhit(seg, xmid, BB.ymin);
            if (hr <= 0.) {
                if (vb >= 0.)
                    return BR;
                else
                    return B;
            }
            else {
                if (vb >= 0.)
                    return R;
                else
                    return NTL;
            }
        }
        else if (hc > 0.) {
            real vt = vhit(seg, xmid, BB.ymax);
            real hl = hhit(seg, BB.xmin, ymid);
            if (vt <= 0.) {
                if (hl >= 0.)
                    return TL;
                else
                    return L;
            }
            else {
                if (hl >= 0)
                    return T;
                else
                    return NBR;
            }
        }
        else
            return A;
    }
    else {
        if (hc < 0.) {
            real hr = hhit(seg, BB.xmax, ymid);
            real vt = vhit(seg, xmid, BB.ymax);
            if (hr < 0.) {
                if (vt >= 0.)
                    return T;
                else
                    return TR;
            }
            else {
                if (vt >= 0.)
                    return NBL;
                else
                    return R;
            }
        }
        else {
            real vb = vhit(seg, xmid, BB.ymin);
            real hl = hhit(seg, BB.xmin, ymid);
            if (vb >= 0.) {
                if (hl >= 0.)
                    return BL;
                else
                    return L;
            }
            else {
                if (hl >= 0)
                    return B;
                else
                    return NTR;
            }
        }
    }
}

// Convert path to shortcuted path
void cut_path(ShortcutTree* tree, Shape shape, Element* element, BoundingBox BB) {
    Curve* curves = shape.path;
    int i;
    for (i = 0; i < shape.num_segments; i++) {
        Curve seg = curves[i];

        if (seg.BB.xmax > BB.xmin && seg.BB.xmin < BB.xmax &&
            seg.BB.ymax > BB.ymin && seg.BB.ymin < BB.ymax) {
            // Bounding box intersection

            Curve* curve = NULL;
            bool computed_lcross = false;
            bool lcross = false;

            if ((BB.xmin <= seg.x0 && BB.ymin <= seg.y0 && seg.x0 <= BB.xmax && seg.y0 <= BB.ymax) ||
                (BB.xmin <= seg.xf && BB.ymin <= seg.yf && seg.xf <= BB.xmax && seg.yf <= BB.ymax)) {
                // Endpoint inside
                curve = curves + i;
            }
            else {
                // Endpoint outside
                computed_lcross = true;
                if (seg.asc) {
                    lcross = (vhit(seg, BB.xmin, BB.ymax) > 0) && (vhit(seg, BB.xmin, BB.ymin) <= 0);
                    if (lcross
                        // Crosses left, top or right borders
                        || (vhit(seg, BB.xmax, BB.ymax) >= 0 && vhit(seg, BB.xmax, BB.ymin) < 0)
                        || (hhit(seg, BB.xmax, BB.ymax) >= 0 && hhit(seg, BB.xmin, BB.ymax) < 0)) {
                        curve = curves + i;
                    }
                }
                else {
                    lcross = (vhit(seg, BB.xmin, BB.ymax) >= 0) && (vhit(seg, BB.xmin, BB.ymin) < 0);
                    if (lcross
                        // Crosses left, top or right borders
                        || (vhit(seg, BB.xmax, BB.ymax) >= 0 && vhit(seg, BB.xmax, BB.ymin) <= 0)
                        || (hhit(seg, BB.xmax, BB.ymax) >= 0 && hhit(seg, BB.xmin, BB.ymax) <= 0)) {
                        curve = curves + i;
                    }
                }
            }

            if (curve) {
                push_element(element, curve);
                // Calculate left crossing if not already done
                if (!computed_lcross) {
                    if (seg.asc)
                        lcross = (vhit(seg, BB.xmin, BB.ymax) > 0) && (vhit(seg, BB.xmin, BB.ymin) <= 0);
                    else
                        lcross = (vhit(seg, BB.xmin, BB.ymax) >= 0) && (vhit(seg, BB.xmin, BB.ymin) < 0);
                }

                // Create shortcut if necessary
                if (lcross) {
                    if (seg.x0 <= BB.xmin) {
                        tree->shortcuts[tree->num_shortcuts] = linear_segment(seg.x0, BB.ymin, seg.x0, seg.y0);
                        curve->sh = tree->shortcuts + tree->num_shortcuts++;
                        push_element(element, curve->sh);
                    }
                    else if (seg.xf <= BB.xmin) {
                        tree->shortcuts[tree->num_shortcuts] = linear_segment(seg.xf, seg.yf, seg.xf, BB.ymin);
                        curve->sh = tree->shortcuts + tree->num_shortcuts++;
                        push_element(element, curve->sh);
                    }
                }
            }
        }
        // Add to initial winding number, independent of intersection
        element->init_w += winding(&seg, BB.xmin, BB.ymax);
    }
}

// Subdivide path into cells
void subdivide_path(ShortcutTree* tree, Element path, ShortcutNode* cell) {
    real xmin = cell->BB.xmin, ymin = cell->BB.ymin;
    real xmax = cell->BB.xmax, ymax = cell->BB.ymax;
    real xmid = (xmin + xmax) * .5, ymid = (ymin + ymax) * .5;

    Element elements[4];
    int i;
    for (i = 0; i < 4; i++) {
        elements[i].shape = path.shape;
        elements[i].curves = cell->cells[i].elements[0].curves + cell->cells[i].total_segments;
        elements[i].size = 0;
        elements[i].init_w = path.init_w;
    }
    for (i = 0; i < path.size; i++) {
        Curve seg = *(path.curves[i]);
        char Class = classify_seg(seg, cell->BB);

        // Left crossing vertical ray tests
        // Ordered as:
        //    4....5....x
        //    | TL | TR |
        //    2....3....x
        //    | BL | BR |
        //    0....1....x
        // x = not computed
        real hits[6] = { 0.0 };
        if (Class & L) {
            hits[2] = vhit(seg, xmin, ymid);
            if (Class & BL)
                hits[0] = vhit(seg, xmin, ymin);
            if (Class & TL)
                hits[4] = vhit(seg, xmin, ymax);
        }
        if (Class & R) {
            hits[3] = vhit(seg, xmid, ymid);
            if (Class & BR)
                hits[1] = vhit(seg, xmid, ymin);
            if (Class & TR)
                hits[5] = vhit(seg, xmid, ymax);
        }

        // Shortcuts
        real xs, ys0, ys1;
        Curve* sh = seg.sh;
        if (!sh) {
            if (seg.x0 <= seg.xf)
                xs = seg.x0, ys0 = ymin, ys1 = seg.y0;
            else
                xs = seg.xf, ys0 = seg.yf, ys1 = ymin;
        }

        int j = 0, subcell = 1;
        for (; subcell < (BL | BR | TL | TR); subcell <<= 1, j += 1) {
            if (Class & subcell) {
                push_element(elements + j, path.curves[i]);
                if (hits[j] <= 0 && hits[j + 2] >= 0) {
                    if (!sh) {
                        tree->shortcuts[tree->num_shortcuts] = linear_segment(xs, ys0, xs, ys1);
                        sh = tree->shortcuts + tree->num_shortcuts++;
                    }
                    push_element(elements + j, sh);
                }
            }
        }
        path.curves[i]->sh = sh;

        if ((Class == 0) || (Class & L)) {
            elements[0].init_w += winding(&seg, xmin, ymid);
            elements[1].init_w += winding(&seg, xmid, ymid);
            elements[3].init_w += winding(&seg, xmid, ymax);
        }
    }

    for (i = 0; i < 4; i++) {
        if (elements[i].size || elements[i].init_w) {
            cell->cells[i].elements[cell->cells[i].size++] = elements[i];
            cell->cells[i].total_segments += elements[i].size;
        }
    }
}

// Divide node into 4 cells. Each cell represents a quadrant of the node
// bounding box cut by the x and y axes centered at the box midpoint, computed
// as the average of the box corner vertices.
void build_tree_rec(ShortcutTree* tree, ShortcutNode* node, int depth) {
    if ((depth >= tree->max_depth) || (node->total_segments <= tree->max_segs[depth]) || node->size == 0) {
        return;
    }

    real xmid = (node->BB.xmin + node->BB.xmax) * .5;
    real ymid = (node->BB.ymin + node->BB.ymax) * .5;

    node->cells = (ShortcutNode*)malloc(sizeof(ShortcutNode) * 4);
    ShortcutNode* cells = node->cells;
    cells[0].BB = create_BoundingBox(node->BB.xmin, node->BB.ymin, xmid, ymid);
    cells[1].BB = create_BoundingBox(xmid, node->BB.ymin, node->BB.xmax, ymid);
    cells[2].BB = create_BoundingBox(node->BB.xmin, ymid, xmid, node->BB.ymax);
    cells[3].BB = create_BoundingBox(xmid, ymid, node->BB.xmax, node->BB.ymax);

    int i;

    for (i = 0; i < 4; i++) {
        cells[i].cells = NULL;
        cells[i].elements = (Element*)malloc(sizeof(Element) * node->size);
        cells[i].elements[0].curves = (Curve**)malloc(sizeof(Curve*) * 2 * node->total_segments);
        cells[i].total_segments = 0;
        cells[i].size = 0;
    }

    for (i = 0; i < node->size; i++) {
        subdivide_path(tree, node->elements[i], node);
    }

    for (i = 0; i < 4; i++) {
        build_tree_rec(tree, node->cells + i, depth + 1);
    }
}

// Build the shortcut tree inner nodes recursively.
// The node elements are formed based on the list of shapes,
// the bounding box of the nodes is given by cell_size and
// the recursion continues until the tree->max_depth is
// reached or node->total_segments <= tree->max_segs[depth].
void build_tree(ShortcutTree* tree, real cell_size) {
    tree->root = (ShortcutNode*)malloc(sizeof(ShortcutNode));
    memset(tree->root, 0, sizeof(ShortcutNode));
    tree->root->BB.xmax = cell_size;
    tree->root->BB.ymax = cell_size;

    if (tree->size <= 0) {
        tree->shortcuts = NULL;
        tree->num_shortcuts = 0;
        tree->max_shortcuts = 0;
        return;
    }

    Shape last_shape = tree->shapes[tree->size - 1];
    int total_curves = (int)(last_shape.path + last_shape.num_segments - tree->shapes[0].path);
    tree->root->elements = (Element*)malloc(sizeof(Element) * tree->size);
    tree->root->elements[0].curves = (Curve**)malloc(sizeof(Curve*) * 3 * total_curves);
    tree->shortcuts = (Curve*)malloc(sizeof(Curve) * total_curves);
    tree->num_shortcuts = 0;
    tree->max_shortcuts = total_curves;
    int i, n;
    for (i = 0, n = 0; i < tree->size; i++) {
        Element e;
        e.shape = tree->shapes + i;
        e.curves = tree->root->elements[0].curves + tree->root->total_segments;
        e.init_w = 0;
        e.size = 0;
        cut_path(tree, tree->shapes[i], &e, tree->root->BB);
        if (e.size || e.init_w) {
            tree->root->elements[n++] = e;
            tree->root->total_segments += e.size;
            tree->root->size++;
        }
    }

    build_tree_rec(tree, tree->root, 0);
}

void delete_nodes(ShortcutNode* node) {
    int i;
    if (node == NULL) {
        return;
    }
    if (node->cells != NULL) {
        for (i = 0; i < 4; i++)
            delete_nodes(node->cells + i);
        free(node->cells);
        node->cells = NULL;
    }
    free(node->elements->curves);
    free(node->elements);
    node->elements = NULL;
}

void delete_tree(ShortcutTree* tree) {
    int i;
    delete_nodes(tree->root);
    free(tree->shortcuts);

    tree->root = NULL;
    tree->shortcuts = NULL;

    free(tree->shapes[0].path);
    for (i = 0; i < tree->size; i++) {
        Paint paint = tree->shapes[i].paint;
        if (paint.type == PAINT_LINEAR)
            free(paint.linear.ramp.stops);
        else if (paint.type == PAINT_RADIAL)
            free(paint.radial.ramp.stops);
    }
    free(tree->shapes);
    free(tree->max_segs);

    tree->shapes = NULL;
    tree->max_segs = NULL;
}

#endif //SHORTCUT_TREE_H