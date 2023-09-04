cmake_minimum_required(VERSION 3.9)
project(vector_graphics)

include_directories(examples/external examples/external/nanosvg/src src)

set(VECTOR_GRAPHICS_SOURCE
    src/vec.h
    src/util.h
    src/bounding_box.h
    src/quadratic.h
    src/bezier.h
    src/curve.h
    src/path.h
    src/filter.h
    src/shortcut_tree.h
    src/raster.h
)

set(EXTERNAL
    examples/external/stb_image_write.h
    examples/external/nanosvg/src/nanosvg.h
)

add_executable(example1
    ${VECTOR_GRAPHICS_SOURCE}
    ${EXTERNAL}
    examples/example1.c
)
    
add_executable(example2
    ${VECTOR_GRAPHICS_SOURCE}
    ${EXTERNAL}
    src/exporter.h
    examples/example2.c
)