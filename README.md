# Vector Graphics
A 2D vector graphics renderer written in C.

## Building
To download all necessary files for this project `clone` this repository with the `--recurse-submodules` parameter.
After downloading all necessary files you may build the project with `CMake`.

The below commands will download this repository and build the project in a `build` directory:
```
$ git clone --recurse-submodules git@github.com:rafaellg79/vector-graphics.git
$ cmake -B build -S vector-graphics
$ cmake --build build
```

The resulting executables may then be found in the `build` directory.

### Running
After building, two example executables will be generated.

The `example1` asks for an output filename parameter and render a predefined loop with 210x210 resolution.
For instance:
```
$ example1 cubic.png
```

Will render the loop and save into a `cubic.png` file.
Note that the path to the file should exist before running the executable.

The `example2` requires two files path, the first an input file in *.svg* format and the second an output file in a supported raster format.
It will read the input file and save a raster version in the output file path.
Most features of the SVG specification are not supported.
