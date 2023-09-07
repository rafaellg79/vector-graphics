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
#include "common.h"

int main(int argc, char * argv[]){
    
    if (argc < 3) {
        printf("Usage: %s <input_file> <output_file> [options]\n", argv[0]);
        return 0;
    }
    
    int max_depth = 0, num_samples = 0;
    bool dump = false;
    int dump_cell_id = 0;
    real* pattern = NULL;
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
    int width = (int)ceil(nsvg_image->width), height = (int)ceil(nsvg_image->height);
    unsigned int * image = (unsigned int *)malloc(sizeof(unsigned int) * width * height);

    if (pattern == NULL) {
        pattern = (real*)malloc(sizeof(real) * 2);
        pattern[0] = 0.0; pattern[1] = 0.0;
        num_samples = 1;
    }
    
    // the depth of a cell with a 1 px^2 area
    int pixel_depth = (int)fmax(0., ceil(log2(fmax(width, height))));
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

        int length = (int)strlen(argv[2])+1;
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