//Standard libraries
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

//External libraries
//stb_image_write
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

//Headers
#include "shortcut_tree.h"
#include "util.h"
#include "vec.h"
#include "filter.h"
#include "raster.h"
#include "common.h"

int main(int argc, char* argv[]) {

    if (argc < 2) {
        printf("Usage: %s <output_file>\n", argv[0]);
        return 0;
    }

    int i, j;
    int max_depth = 10, num_samples = 1;
    real pattern[] = {0.0, 0.0};

    int width = 210, height = 210;
    unsigned int* image = (unsigned int*)malloc(sizeof(unsigned int) * width * height);

    if (!max_depth) {
        max_depth = fmax(0., ceil(log2(fmax(width, height))));
    }

    int* max_segs = (int*)malloc(sizeof(int) * max_depth);
    int depth;
    for (depth = 0; depth < max_depth; depth++) {
        max_segs[depth] = (int)fmax(1, depth * depth);
    }

    // Tree creation
    ShortcutTree tree;
    tree.shapes = (Shape *)malloc(sizeof(Shape));
    tree.shapes[0] = loop();
    tree.size = 1;
    tree.max_depth = max_depth;
    tree.max_segs = max_segs;

    real cell_size = pow(2., ceil(log2(fmax(height, width))));
    build_tree(&tree, cell_size);

    vec4 white = { 1.0, 1.0, 1.0, 1.0 };

    printf("Rendering...\n");
    for (i = 0; i < height; i++) {
        real y = i + 0.5;
        int j;
        for (j = 0; j < width; j++) {
            real x = j + 0.5;
            image[i * width + j] = pack_vec(supersample(tree.root, pattern, num_samples, x, y, white));
        }
        fprintf(stderr, "\r%.2f%%", 100. * i / height);
    }
    fprintf(stderr, "\r100.00%%");

    printf("Finished rendering.\n");

    printf("Saving to file...\n");

    int error_code = stbi_write_png(argv[1], width, height, 4, (const unsigned char*)image, 0);
    if (error_code == 0) {
        fprintf(stderr, "error: %i\n", error_code);
        return 1;
    }
    printf("File saved successfully.\n");

    free(image);
    delete_tree(&tree);

    return 0;
}