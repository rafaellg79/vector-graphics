#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include <stdlib.h>

#include "vec.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct BoundingBox{
    real xmin, ymin, xmax, ymax;
} BoundingBox;

inline BoundingBox create_BoundingBox(const real xmin, const real ymin, const real xmax, const real ymax);

#ifdef __cplusplus
}
#endif

// Creates a BoundingBox object representing the bounding box given by 2 corner points
// p0 = (min(x0, x1), min(y0, y1)) and p1 = (max(x0, x1), max(y0, y1))
inline BoundingBox create_BoundingBox(const real x0, const real y0, const real x1, const real y1){
    BoundingBox b;
    if(x0 < x1){
        b.xmin = x0;
        b.xmax = x1;
    }
    else{
        b.xmin = x1;
        b.xmax = x0;
    }
    if(y0 < y1){
        b.ymin = y0;
        b.ymax = y1;
    }
    else{
        b.ymin = y1;
        b.ymax = y0;
    }
    return b;
}

#endif //BOUNDING_BOX_H