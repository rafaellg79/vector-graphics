#ifndef VEC_H
#define VEC_H

#include <stdio.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef double real;

typedef struct vec1{
    real x;
} vec1;

typedef struct vec2{
    real x, y;
} vec2;

typedef struct vec3{
    real x, y, z;
} vec3;

typedef struct vec4{
    real x, y, z, w;
} vec4;

#ifdef __cplusplus
}
#endif

#endif //VEC_H