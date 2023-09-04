#ifndef BEZIER_H
#define BEZIER_H

#include "vec.h"
#include "util.h"

#ifdef __cplusplus
extern "C" {
#endif

enum Cubic{
    //Cubic bezier curve classes
    SERPENTINE,
    LOOP,
    CUSP_WITH_INFLECTION_AT_INFINITY,
    CUSP_WITH_CUSP_AT_INFINITY,
    //Degenerate cases
    DEGENERATED_QUADRATIC,
    LINE_OR_POINT
};

void cut2(vec2 *p, const real  a, const real  b,
                   const real x0, const real y0, 
                   const real x1, const real y1, 
                   const real x2, const real y2);
void cut3(vec2 *p, const real  a, const real  b, 
                   const real x0, const real y0, 
                   const real x1, const real y1, 
                   const real x2, const real y2, 
                   const real x3, const real y3);
inline vec2 at1(real t, real x0, real y0, real x1, real y1);
vec2 at2(real t, real x0, real y0, real x1, real y1, real x2, real y2);
vec2 at3(real t, real x0, real y0, real x1, real y1, real x2, real y2, real x3, real y3);
vec4 crosspmatrix(real x0, real y0, real x1, real y1, real x2, real y2, real x3, real y3);
enum Cubic classify3(real x0, real y0, real x1, real y1, real x2, real y2, real x3, real y3);
enum Cubic regenerate(real * pts);

#ifdef __cplusplus
}
#endif

inline real lerp1(real const a, real const x0, real const x1){
    if (a == 1.0) return x1;
    return x0 + a * (x1 - x0);
}

inline real lerp2(real const a, real const b, real const x0, real const x1, real const x2){
    real x00 = lerp1(a, x0, x1);
    real x01 = lerp1(a, x1, x2);
    
    return lerp1(b, x00, x01);
}

inline real lerp3(real const a, real const b, real const c, real const x0, real const x1, real const x2, real const x3){
    real x00 = lerp1(a, x0, x1);
    real x01 = lerp1(a, x1, x2);
    real x02 = lerp1(a, x2, x3);
    
    return lerp2(b, c, x00, x01, x02);
}
/*
// generalized lerp using the recursive form of the De Casteljau's algorithm
real lerp_rec(int n, real * t, real * x){
    if(n == 1) return lerp1(t[0], x[0], x[1]);
    real a = lerp_rec(n-1, t, x);
    real b = lerp_rec(n-1, t, x+1);
    return lerp1(t[n-1], a, b);
}

// generalized lerp using the iterative form of the De Casteljau's algorithm
real lerp_iter(int n, real * t, real * x){
    real * xn = (real *)malloc(sizeof(real) * (n+1));

    int i;
    for(i = 0; i < n+1; i++){
        xn[i] = x[i];
    }
    
    for(i = n; i >= 0; i--){
        int j;
        for(j = 0; j < i; j++){
            xn[j] = lerp1(t[j], xn[j], xn[j+1]);
        }
    }
    
    real result = xn[0];
    free(xn);
    return result;
}*/

static inline vec3 cut2_base(real a, real b, real x0, real x1, real x2){
    real x00 = lerp1(a, x0, x1);
    real x01 = lerp1(a, x1, x2);
    
    vec3 u;
    u.x = lerp1(a, x00, x01);
    u.y = lerp1(b, x00, x01);
    u.z = lerp2(b, b, x0, x1, x2);
    
    return u;
}

static inline vec4 cut3_base(real a, real b, real x0, real x1, real x2, real x3){
    real x00 = lerp1(a, x0, x1);
    real x01 = lerp1(a, x1, x2);
    real x02 = lerp1(a, x2, x3);
    
    real x10 = lerp1(a, x00, x01);
    real x11 = lerp1(a, x01, x02);
    
    vec4 u;
    
    u.x = lerp1(a, x10, x11);
    u.y = lerp1(b, x10, x11);
    
    x10 = lerp1(b, x00, x01);
    x11 = lerp1(b, x01, x02);
    u.z = lerp1(b, x10, x11);
    
    u.w = lerp3(b, b, b, x0, x1, x2, x3);
    return u;
}

inline vec3 power2(real x0, real x1, real x2){
    vec3 v = {x0, 2.f*(x1-x0), (x2-x1)-(x1-x0)};
    return v;
}

inline vec4 power3(real x0, real x1, real x2, real x3){
    vec4 v = {x0, 3.f*(x1-x0), 3.f*(x2-2.f*x1+x0), x3-x0-3.f*(x2-x1)};
    return v;
}

void cut2(vec2 * p, const real  a, const real  b,
                    const real x0, const real y0,
                    const real x1, const real y1,
                    const real x2, const real y2){
    vec3 u = cut2_base(a, b, x0, x1, x2);
    vec3 v = cut2_base(a, b, y0, y1, y2);
    p[0].x = u.x; p[0].y = v.x;
    p[1].x = u.y; p[1].y = v.y;
    p[2].x = u.z; p[2].y = v.z;
}

void cut3(vec2 * p, const real a, const real b,
                       const real x0, const real y0,
                       const real x1, const real y1,
                       const real x2, const real y2,
                       const real x3, const real y3){
    vec4 u = cut3_base(a, b, x0, x1, x2, x3);
    vec4 v = cut3_base(a, b, y0, y1, y2, y3);
    p[0].x = u.x; p[0].y = v.x;
    p[1].x = u.y; p[1].y = v.y;
    p[2].x = u.z; p[2].y = v.z;
    p[3].x = u.w; p[3].y = v.w;
}

vec2 at1(real t, real x0, real y0, real x1, real y1){
    vec2 p;
    p.x = lerp1(t, x0, x1);
    p.y = lerp1(t, y0, y1);
    return p;
}

vec2 at2(real t, real x0, real y0, real x1, real y1, real x2, real y2){
    vec2 p;
    p.x = lerp2(t, t, x0, x1, x2);
    p.y = lerp2(t, t, y0, y1, y2);
    return p;
}

vec2 at3(real t, real x0, real y0, real x1, real y1, real x2, real y2, real x3, real y3){
    vec2 p;
    p.x = lerp3(t, t, t, x0, x1, x2, x3);
    p.y = lerp3(t, t, t, y0, y1, y2, y3);
    return p;
}

// Returns the tangents lines to a cubic bezier at t=0 and t=1
// Assumes the bezier curve is not degenerate to a point
void cubic_tan(vec3 * d0, vec3 * d1, const real x0, const real y0, const real x1, const real y1, const real x2, const real y2, const real x3, const real y3){
    real xt = x1, yt = y1;
    if (is_almost_zero(norm2(xt-x0, yt-y0))) {
        xt = x2; yt = y2;
        if (is_almost_zero(norm2(xt-x0, yt-y0))) {
            xt = x3; yt = y3;
        }
    }
    *d0 = cross3(xt, yt, 1., x0, y0, 1.);

    xt = x2; yt = y2;
    if (is_almost_zero(norm2(xt-x3, yt-y3))) {
        xt = x1; yt = y1;
        if (is_almost_zero(norm2(xt-x3, yt-y3))) {
            xt = x0; yt = y0;
        }
    }
    *d1 = cross3(x3, y3, 1., xt, yt, 1.);
}

vec4 crosspmatrix(real x0, real y0, real x1, real y1, real x2, real y2, real x3, real y3){
    vec4 u = power3(x0, x1, x2, x3);
    vec4 v = power3(y0, y1, y2, y3);
    vec4 d;
    d.x = 0.f;
    d.y = det2(u.z, u.w, v.z, v.w);
    d.z = -det2(u.y, u.w, v.y, v.w);
    d.w = det2(u.y, u.z, v.y, v.z);
    return d;
}

enum Cubic classify3(real x0, real y0, real x1, real y1, real x2, real y2, real x3, real y3){
    vec4 d = crosspmatrix(x0, y0, x1, y1, x2, y2, x3, y3);
    d.x = 3.f*d.z*d.z-4.f*d.y*d.w;
    if (!is_almost_zero(d.y)) {
        if(d.x > 0.f)
            return SERPENTINE;
        else if(d.x < 0.f)
            return LOOP;
        else
            return CUSP_WITH_INFLECTION_AT_INFINITY;
    }
    else if(!is_almost_zero(d.z))
        return CUSP_WITH_CUSP_AT_INFINITY;
    else if(!is_almost_zero(d.w))
        return DEGENERATED_QUADRATIC;
    else
        return LINE_OR_POINT;
}

//Receives the control points from a cubic bezier curve,
//then if the bezier is degenerated to a lower order curve
//modify the points to store only the necessary points
//and return the type of the curve
enum Cubic regenerate(real* pts){
    enum Cubic type = classify3(pts[0], pts[1], pts[2], pts[3], pts[4], pts[5], pts[6], pts[7]);
    if (type == LINE_OR_POINT) {
        pts[2] = pts[6]; pts[3] = pts[7];
    }
    else if (type == DEGENERATED_QUADRATIC) {
        vec3 d0, d1;
        cubic_tan(&d0, &d1, pts[0], pts[1], pts[2], pts[3], pts[4], pts[5], pts[6], pts[7]);
        vec3 p = cross3(d0.x, d0.y, d0.z, d1.x, d1.y, d1.z);
        
        pts[2] = p.x / p.z; pts[3] = p.y / p.z;
        pts[4] = pts[6]; pts[5] = pts[7];
    }
    return type;
}

#endif //BEZIER_H