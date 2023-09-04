#ifndef PATH_H
#define PATH_H

#include <stdlib.h>

#include "vec.h"
#include "curve.h"
#include "paint.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef bool (*winding_function)(int);

typedef struct Shape{
    BoundingBox BB;
    Curve * path;
    int num_segments;
    Paint paint;
    winding_function winding_test;
} Shape;

typedef struct Element{
    Shape * shape;
    Curve ** curves;
    int size;
    int init_w;
} Element;

inline void push_element(Element* e, Curve* c);
inline void push_shape(Shape* s, Curve c);
int inside_shape(Shape shape, real x, real y);
int inside_element(Element* element, real x, real y);

#ifdef __cplusplus
}
#endif

inline void push_element(Element * e, Curve * c){
    e->curves[e->size++] = c;
}

inline void push_shape(Shape * s, Curve c){
    s->path[s->num_segments++] = c;
}

int inside_shape(Shape shape, real x, real y){
    int w = 0;
    Curve * curves = shape.path;
    int i;
    for(i = 0; i < shape.num_segments; i++){
        w = w + winding(curves+i, x, y);
    }
    return w;
}

int inside_element(Element * element, real x, real y){
    int i = 0, w = 0;
    for(i = 0; i < element->size; i++){
        w = w + winding(element->curves[i], x, y);
    }
    return w + element->init_w;
}

#endif //PATH_H