#ifndef FILTER_H
#define FILTER_H

#include <stdlib.h>

#include "vec.h"
#include "util.h"
#include "quadratic.h"
#include "bezier.h"
#include "path.h"

#ifdef __cplusplus
extern "C" {
#endif

void monotonize(Shape* shape, real* pts, enum Cubic type);

#ifdef __cplusplus
}
#endif

int st_float_cmp(const void * a, const void * b){
    const real* f_a = (const real*)a;
    const real* f_b = (const real*)b;
    return (*f_a > *f_b) - (*f_a < *f_b);
}

inline void uniq(real * t, int * n){
    qsort(t, *n, sizeof(real), st_float_cmp);
    int i = 0, j = 1;
    while (j < *n) {
        if (t[i] == t[j]) {
            j = j + 1;
        }
        else{
            i = i + 1;
            t[i] = t[j];
            j = j + 1;
        }
    }
    *n = i + 1;
}

inline void filterroot(real * ts, int * m, real t, real s){
    real ss = sign(s);
    if(0.f < ss * t && ss * t < ss * s){
        ts[(*m)++] = t / s;
    }
}

void monotonize_quadratic(Shape * shape, real x0, real y0, real x1, real y1, real x2, real y2){
    real ts[4];
    ts[0] = 0.f;
    int m = 1;
    
    // monotonize in x
    real t = x0 - x1;
    real s = x0 - 2.f*x1 + x2;
    filterroot(ts, &m, t, s);
    
    // monotonize in y
    t = y0 - y1;
    s = y0 - 2.f*y1 + y2;
    filterroot(ts, &m, t, s);
    
    ts[m++] = 1.;
    uniq(ts, &m);
    real xp = x0, yp = y0;
    for (int i = 1; i < m; i++){
        vec2 cp[3];
        cut2(cp, ts[i-1], ts[i], x0, y0, x1, y1, x2, y2);
        Curve qs = quadratic_segment(xp, yp, cp[1].x, cp[1].y, cp[2].x, cp[2].y);
        push_shape(shape, qs);
        xp = cp[2].x, yp = cp[2].y;
    }
}

//Break a cubic bezier curve segment into several monotonic cubic segments
//also removes double and inflection points
void monotonize_cubic(Shape * shape, real x0, real y0, real x1, real y1, real x2, real y2, real x3, real y3){
    vec4 d = crosspmatrix(x0, y0, x1, y1, x2, y2, x3, y3);
    d.x = 3.f*d.z*d.z - 4.f*d.y*d.w;
    real ts[8];
    ts[0] = 0.f;
    int n = 0, m = 1;
    
    // monotonize in x
    real a = -x0 + 3.f*(x1 - x2) + x3;
    real b = 2.f*(x0 - 2.f*x1 + x2);
    real c = x1 - x0;
    vec4 t = solve_quadratic(&n, a, b, c);
    if (n > 0) filterroot(ts, &m, t.x, t.y);
    if (n > 1) filterroot(ts, &m, t.z, t.w);
    
    // monotonize in y
    a = -y0 + 3.f*(y1 - y2) + y3;
    b = 2.f*(y0 - 2.f*y1 + y2);
    c = y1 - y0;
    t = solve_quadratic(&n, a, b, c);
    if (n > 0) filterroot(ts, &m, t.x, t.y);
    if (n > 1) filterroot(ts, &m, t.z, t.w);
    
    // break at double point
    if (!is_almost_zero(d.y) && d.x < 0.f) {
        t = solve_quadratic_precomputed_delta(&n, d.y*d.y, -d.y*d.z, d.z*d.z-d.y*d.w, -.25*d.y*d.y*d.x);
    }
    // break at inflections
    else {
        t = solve_quadratic_precomputed_delta(&n, -3.f*d.y, 3.f*d.z, -d.w, .25*3.*d.x);
    }
    if (n > 0) filterroot(ts, &m, t.x, t.y);
    if (n > 1) filterroot(ts, &m, t.z, t.w);
    
    ts[m++] = 1.;
    uniq(ts, &m);
    real xp = x0, yp = y0;
    int i;
    for(i = 1; i < m; i++){
        vec2 cp[4];
        cut3(cp, ts[i-1], ts[i], x0, y0, x1, y1, x2, y2, x3, y3);
        Curve cs = cubic_segment(xp, yp, cp[1].x, cp[1].y, cp[2].x, cp[2].y, cp[3].x, cp[3].y);
        push_shape(shape, cs);
        xp = cp[3].x, yp = cp[3].y;
    }
}

void monotonize(Shape * shape, real* pts, enum Cubic type){
    if (type == DEGENERATED_QUADRATIC) {
        monotonize_quadratic(shape, pts[0], pts[1], pts[2], pts[3], pts[4], pts[5]);
    } else if (type == LINE_OR_POINT) {
        if(!is_almost_equal(pts[0], pts[2]) || !is_almost_equal(pts[1], pts[3])){
            Curve ls = linear_segment(pts[0], pts[1], pts[2], pts[3]);
            push_shape(shape, ls);
        }
    } else {
        monotonize_cubic(shape, pts[0], pts[1], pts[2], pts[3], pts[4], pts[5], pts[6], pts[7]);
    }
}

#endif // FILTER_H