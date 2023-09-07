//Standard libraries
#include <stdio.h>
#include <math.h>

//External libraries
#include "nanosvg.h"

//Headers
#include "shortcut_tree.h"
#include "util.h"
#include "vec.h"
#include "filter.h"
#include "raster.h"

#define NUM_THREADS 16

unsigned int pack(unsigned char r, unsigned char g, unsigned char b, unsigned char a) {
    return (a << 24) | (b << 16) | (g << 8) | r;
}

unsigned int pack_vec(vec4 c) {
    unsigned char r = (char)(255 * c.x), g = (char)(255 * c.y), b = (char)(255 * c.z), a = (char)(255 * c.w);
    return pack(r, g, b, a);
}

vec4 unpack(unsigned int c) {
    vec4 color;
    color.x = (c & 0x000000FF) / 255.;
    color.y = ((c & 0x0000FF00) >> 8) / 255.;
    color.z = ((c & 0x00FF0000) >> 16) / 255.;
    color.w = ((c & 0xFF000000) >> 24) / 255.;
    return color;
}

real* read_pattern(int* n, const char* filename) {
    FILE* file = fopen(filename, "r");
    if (file == NULL || fscanf(file, "%i\n", n) != 1)
        return NULL;
    real* pattern = (real*)malloc(sizeof(real) * 2 * (*n));
    int i;
    for (i = 0; i < 2 * (*n); i += 2) {
        if (fscanf(file, "%lf %lf\n", pattern + i, pattern + i + 1) != 2) {
            fprintf(stderr, "Error reading pattern");
            free(pattern);
            return NULL;
        }
    }
    return pattern;
}

// Assumes shape->path has enough space to store the curves after implicitization
// It should be enough if shape->path[shape->num_segments + npts] is a valid memory address
void fill_shape(Shape * shape, real * pts, int npts){
    int i, j;
    for (i = 0; i + 3 < npts; i += 3) {
        real curve_pts[8];
        for (j = 0; j < 8; j += 2) {
            curve_pts[j] = pts[i * 2 + j];
            curve_pts[j + 1] = pts[i * 2 + j + 1];
        }
        enum Cubic type = regenerate(curve_pts);
        monotonize(shape, curve_pts, type);
    }
}

Shape loop(){
    real pts[] = { 195.500000,  15.500000, 195.500000,  15.500000,  195.500000, 195.500000, // right line
                   195.500000, 195.500000, 195.500000, 195.500000,   15.500000, 195.500000, // bottom line
                    15.500000, 195.500000,  15.500000, 195.500000,   15.500000,  15.500000, // left line
                    15.500000,  15.500000, 345.500000, 255.500000, -134.500000, 255.500000, // loop
                   195.500000,  15.500000 };
    int npts = 13;

    Paint paint = { PAINT_COLOR, 1.0, {1.0, 0.0, 0.0, 1.0} };
    Shape s;
    s.path = (Curve*)malloc(sizeof(Curve) * 16);
    s.paint = paint;
    s.winding_test = evenodd;
    s.BB = create_BoundingBox(15.5, 15.5, 195.5, 195.5);
    s.num_segments = 0;
    fill_shape(&s, pts, npts);
    return s;
}

winding_function get_winding_test(char rule) {
    switch (rule) {
    case NSVG_FILLRULE_NONZERO:
        return nonzero;
    case NSVG_FILLRULE_EVENODD:
        return evenodd;
    default:
        fprintf(stderr, "Unknown fill rule: %c\n", rule);
        return NULL;
    }
}

spread_function get_spread_function(char spread) {
    switch (spread) {
    case NSVG_SPREAD_PAD:
        return pad;
    case NSVG_SPREAD_REFLECT:
        return reflect;
    case NSVG_SPREAD_REPEAT:
        return repeat;
    default:
        fprintf(stderr, "Unknown spread: %c\n", spread);
        return NULL;
    }
}

Paint NSVGpaint_to_paint(NSVGpaint nsvg_paint) {
    Paint paint;
    switch (nsvg_paint.type) {
    case NSVG_PAINT_NONE:
        paint.type = PAINT_NONE;
        return paint;
    case NSVG_PAINT_COLOR:
        paint.type = PAINT_COLOR;
        paint.color = unpack(nsvg_paint.color);
        return paint;
    case NSVG_PAINT_LINEAR_GRADIENT:
    {
        paint.type = PAINT_LINEAR;
        NSVGgradient* gradient = nsvg_paint.gradient;
        Ramp ramp;
        ramp.nstops = gradient->nstops;
        ramp.stops = (RampStop*)malloc(sizeof(RampStop) * ramp.nstops);
        int i;
        for (i = 0; i < ramp.nstops; i++) {
            ramp.stops[i].offset = gradient->stops[i].offset;
            ramp.stops[i].color = unpack(gradient->stops[i].color);
        }
        spread_function spread = get_spread_function(gradient->spread);

        LinearGradient* linear = &paint.linear;
        linear->a = (real)gradient->xform[1];
        linear->b = (real)gradient->xform[3];
        linear->c = (real)gradient->xform[5];
        linear->spread = spread;
        linear->ramp = ramp;
        return paint;
    }
    case NSVG_PAINT_RADIAL_GRADIENT:
    {
        paint.type = PAINT_RADIAL;
        NSVGgradient* gradient = nsvg_paint.gradient;
        Ramp ramp;
        ramp.nstops = gradient->nstops;
        ramp.stops = (RampStop*)malloc(sizeof(RampStop) * ramp.nstops);
        int i;
        for (i = 0; i < ramp.nstops; i++) {
            ramp.stops[i].offset = gradient->stops[i].offset;
            ramp.stops[i].color = unpack(gradient->stops[i].color);
        }
        spread_function spread = get_spread_function(gradient->spread);

        real xform[6];
        for (i = 0; i < 6; i++)
            xform[i] = (real)gradient->xform[i];
        paint.radial = radial_gradient(gradient->fx, gradient->fy, xform, spread, ramp);
        return paint;
    }
    default:
        fprintf(stderr, "No valid paint type");
        memset(&paint, 0, sizeof(Paint));
        return paint;
    }
}

// Generates a shortcut tree from the scene
ShortcutTree preprocess(NSVGimage* scene, int max_depth, int* max_segs) {
    int n = 0, m = 0;

    //Calculate the maximum amount of points that may be allocate after implicitization
    for (NSVGshape* shape = scene->shapes; shape != NULL; shape = shape->next) {
        for (NSVGpath* path = shape->paths; path != NULL; path = path->next) {
            n += (path->npts - 1) * 2;
        }
        m++;
    }

    ShortcutTree tree;
    if (n == 0) {
        tree.shapes = NULL;
        tree.size = 0;
        return tree;
    }

    tree.shapes = (Shape*)malloc(sizeof(Shape) * m);
    tree.shapes[0].path = (Curve*)malloc(sizeof(Curve) * n);
    tree.size = 0;
    tree.max_depth = max_depth;
    tree.max_segs = max_segs;

    for (NSVGshape* shape = scene->shapes; shape != NULL; shape = shape->next) {
        Shape* s = tree.shapes + tree.size;
        if (tree.size++) {
            s->path = (s - 1)->path + (s - 1)->num_segments;
        }
        s->paint = NSVGpaint_to_paint(shape->fill);
        s->paint.opacity = shape->opacity;
        s->winding_test = get_winding_test(shape->fillRule);
        s->BB = create_BoundingBox(shape->bounds[0], shape->bounds[1], shape->bounds[2], shape->bounds[3]);
        s->num_segments = 0;
        int i, j;
        for (NSVGpath* path = shape->paths; path != NULL; path = path->next) {
            for (i = 0; i + 3 < path->npts; i += 3) {
                real pts[8];
                for (j = 0; j < 8; j += 2) {
                    pts[j] = path->pts[i * 2 + j];
                    pts[j + 1] = path->pts[i * 2 + j + 1];
                }
                enum Cubic type = regenerate(pts);
                monotonize(s, pts, type);
            }
        }
    }

    real cell_size = pow(2., ceil(log2(fmax(scene->height, scene->width))));
    build_tree(&tree, cell_size);

    return tree;
}