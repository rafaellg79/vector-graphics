//Standard libraries
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

//External libraries
//nanosvg
#define NANOSVG_IMPLEMENTATION
#include "nanosvg.h"

//stb_image_write
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

//Headers
#include "shortcut_tree.h"
#include "util.h"
#include "vec.h"
#include "filter.h"
#include "exporter.h"
#include "raster.h"

unsigned int pack(unsigned char r, unsigned char g, unsigned char b, unsigned char a) {
    return (a << 24) | (b << 16) | (g << 8) | r;
}

unsigned int pack_vec(vec4 c){
    unsigned char r = 255 * c.x, g = 255 * c.y, b = 255 * c.z, a = 255 * c.w;
    return pack(r, g, b, a);
}

vec4 unpack(unsigned int c){
    vec4 color;
    color.x =  (c & 0x000000FF)        / 255.;
    color.y = ((c & 0x0000FF00) >>  8) / 255.;
    color.z = ((c & 0x00FF0000) >> 16) / 255.;
    color.w = ((c & 0xFF000000) >> 24) / 255.;
    return color;
}

real* read_pattern(int *n, const char * filename){
    FILE *file = fopen(filename, "r");
    if(file == NULL || fscanf(file, "%i\n", n) != 1)
        return NULL;
    real* pattern = (real*)malloc(sizeof(real) * 2 * (*n));
    int i;
    for(i = 0; i < 2*(*n); i+=2){
        if (fscanf(file, "%lf %lf\n", pattern + i, pattern + i + 1) != 2) {
            fprintf(stderr, "Error reading pattern");
            free(pattern);
            return NULL;
        }
    }
    return pattern;
}

winding_function get_winding_test(char rule){
    switch(rule){
        case NSVG_FILLRULE_NONZERO:
            return nonzero;
        case NSVG_FILLRULE_EVENODD:
            return evenodd;
        default:
            fprintf(stderr, "Unknown fill rule: %c\n", rule);
            return NULL;
    }
}

spread_function get_spread_function(char spread){
    switch(spread){
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

Paint NSVGpaint_to_paint(NSVGpaint nsvg_paint){
    Paint paint;
    switch(nsvg_paint.type){
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
            NSVGgradient *gradient = nsvg_paint.gradient;
            Ramp ramp;
            ramp.nstops = gradient->nstops;
            ramp.stops = (RampStop *)malloc(sizeof(RampStop) * ramp.nstops);
            int i;
            for(i=0; i<ramp.nstops; i++){
                ramp.stops[i].offset = gradient->stops[i].offset;
                ramp.stops[i].color = unpack(gradient->stops[i].color);
            }
            spread_function spread = get_spread_function(gradient->spread);
            
            LinearGradient *linear = &paint.linear;
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
            NSVGgradient *gradient = nsvg_paint.gradient;
            Ramp ramp;
            ramp.nstops = gradient->nstops;
            ramp.stops = (RampStop *)malloc(sizeof(RampStop) * ramp.nstops);
            int i;
            for(i=0; i<ramp.nstops; i++){
                ramp.stops[i].offset = gradient->stops[i].offset;
                ramp.stops[i].color = unpack(gradient->stops[i].color);
            }
            spread_function spread = get_spread_function(gradient->spread);
            
            real xform[6];
            for(i = 0; i < 6; i++)
                xform[i] = (real)gradient->xform[i];
            paint.radial = radial_gradient(gradient->fx, gradient->fy, xform, spread, ramp);
            return paint;
        }
        default:
            fprintf(stderr, "No valid paint type");
            return paint;
    }
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

// Generates a shortcut tree from the scene
ShortcutTree preprocess(NSVGimage * scene, int max_depth, int * max_segs){
    int n = 0, m = 0;
    
    //Calculate the maximum amount of points that may be allocate after implicitization
    for (NSVGshape * shape = scene->shapes; shape != NULL; shape = shape->next) {
        for (NSVGpath * path = shape->paths; path != NULL; path = path->next) {
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
    
    tree.shapes = (Shape *)malloc(sizeof(Shape) * m);
    tree.shapes[0].path = (Curve *)malloc(sizeof(Curve) * n);
    tree.size = 0;
    tree.max_depth = max_depth;
    tree.max_segs = max_segs;
    
    for (NSVGshape * shape = scene->shapes; shape != NULL; shape = shape->next) {
        Shape * s = tree.shapes + tree.size;
        if(tree.size++){
            s->path = (s-1)->path+(s-1)->num_segments;
        }
        s->paint = NSVGpaint_to_paint(shape->fill);
        s->paint.opacity = shape->opacity;
        s->winding_test = get_winding_test(shape->fillRule);
        s->BB = create_BoundingBox(shape->bounds[0], shape->bounds[1], shape->bounds[2], shape->bounds[3]);
        s->num_segments = 0;
        for (NSVGpath * path = shape->paths; path != NULL; path = path->next) {
            real * pts = (real *)malloc(sizeof(real) * 2 * path->npts);
            for(int i = 0; i < 2*path->npts; i++) pts[i] = path->pts[i];
            fill_shape(s, pts, path->npts);
        }
    }
    
    real cell_size = pow(2., ceil(log2(fmax(scene->height, scene->width))));
    build_tree(&tree, cell_size);
    
    return tree;
}

int main(int argc, char * argv[]){
    
    if (argc < 3) {
        printf("Usage: %s <input_file> <output_file> [options]\n", argv[0]);
        return 0;
    }
    
    int max_depth = 0, num_samples = 0;
    bool dump = false;
    int dump_cell_id = 0;
    real* pattern = NULL;
    char * file_name;
    int i;

    // Simple parsing for example applications do NOT use in real applications
    // Use cargs, getopt (POSIX), Argp (GNU), Boost.Program_options (Boost, if using C++), etc
    for(i = 3; i < argc; i++){
        const char * arg = argv[i];
        if(strcmp(arg, "-dump") == 0){
            dump = true;
        }else if(strcmp(arg, "-dump-cell") == 0){
            dump = true;
            if(i+1 < argc){
                dump_cell_id = atoi(argv[++i]);
            }
        }else if(strcmp(arg, "-depth") == 0){
            if(i+1 < argc){
                max_depth = atoi(argv[++i]);
            }
        }else if(strcmp(arg, "-pattern") == 0){
            if(i+1 < argc){
                pattern = read_pattern(&num_samples, argv[++i]);
            }
        }else{
            fprintf(stderr, "ignored parameter: %s\n", arg);
        }
    }
    
	NSVGimage* nsvg_image;
	nsvg_image = nsvgParseFromFile(argv[1], "px", 96);
    if(nsvg_image == NULL){
        fprintf(stderr, "Error reading svg file!");
        return 1;
    }
    int width = ceil(nsvg_image->width), height = ceil(nsvg_image->height);
    unsigned int * image = (unsigned int *)malloc(sizeof(unsigned int) * width * height);

    if (pattern == NULL) {
        pattern = (real*)malloc(sizeof(real) * 2);
        pattern[0] = 0.0; pattern[1] = 0.0;
        num_samples = 1;
    }
    
    // the depth of a cell with a 1 px^2 area
    int pixel_depth = fmax(0., ceil(log2(fmax(width, height))));
    if(!max_depth){
        max_depth = pixel_depth;
    }
    
    int * max_segs = (int *)malloc(sizeof(int) * max_depth);
    int depth;
    real inv_cell_samples = 1.0 / (((1 << pixel_depth) << pixel_depth) * num_samples);
    for(depth = 0; depth < max_depth; depth++){
        max_segs[depth] = (int)(pow(1.5, depth) * inv_cell_samples);
        inv_cell_samples *= 4;
    }
    
    //preprocess...
    printf("Preprocessing...\n");
    clock_t start = clock();
    ShortcutTree tree = preprocess(nsvg_image, max_depth, max_segs);
    clock_t finish = clock();
    real elapsed = (real)(finish - start) / CLOCKS_PER_SEC;
    fprintf(stderr, "Preprocessing elapsed time: %lfs\n", elapsed);
    printf("Finished preprocessing.\n");
    
    //rendering nsvgimage...
    
    printf("Rendering...\n");
    
    start = clock();
	vec4 white = {1.0, 1.0, 1.0, 1.0};

    int processed_lines = 0;
    for (i = 0; i < height; i++) {
        fprintf(stderr, "\r%.2f%%", 100. * i / height);
        real y = i + 0.5;
        int j;
        for (j = 0; j < width; j++) {
            real x = j + 0.5;
            image[i * width + j] = pack_vec(supersample(tree.root, pattern, num_samples, x, y, white));
        }
    }
    fprintf(stderr, "\r100.00%%");

    finish = clock();
    
    elapsed = (real)(finish - start) / CLOCKS_PER_SEC;
    fprintf(stderr, "\nRendering elapsed time: %lfs\n", elapsed);
    printf("Finished rendering.\n");
    
    printf("Saving to file...\n");
    
    int error_code = stbi_write_png(argv[2], width, height, 4, (const unsigned char *)image, 0);
    if(error_code == 0){
        fprintf(stderr, "error: %i\n", error_code);
        return 1;
    }
    printf("File saved successfully.\n");

    free(image);
    if(dump){
        printf("dumping tree\n");

        int length = strlen(argv[2])+1;
        char * filename = (char *)malloc(sizeof(char) * (length+4));
        char * extension;
        strcpy(filename, argv[2]);
        extension = filename + length;
        while (--extension != filename && *extension != '.');
        strcpy(extension, ".svg");

        // TODO: escape percentage from svg file when writing
        char * svg = cell_to_svg(tree.root, dump_cell_id, 10);

        if (write(filename, svg)) {
            printf("finished dumping tree\n");
        }
        else {
            printf("unable to write dump file\n");
        }

        free(filename);
        free(svg);
    }
    
	nsvgDelete(nsvg_image);
    delete_tree(&tree);
    free(pattern);

    return 0;
}