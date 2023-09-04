#ifndef UTIL_H
#define UTIL_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#define MAX_REL (8*FLT_EPSILON)
#define MAX_ABS (8*FLT_EPSILON)

#ifdef __cplusplus
extern "C" {
#endif

inline bool is_almost_equal(real a, real b);
inline bool is_almost_zero(real a);
inline bool is_almost_one(real f);
inline real chop(real v);
inline real det2(real a, real b, real c, real d);
bool xform_inverse(real* inv, real* t);
inline real norm2(real x, real y);
inline vec3 cross3(real x0, real y0, real w0, real x1, real y1, real w1);
inline real sign(real v);
inline bool write(const char* file, const char* text);

#ifdef __cplusplus
}
#endif

inline bool is_almost_equal(real a, real b){
    real diff = fabs(a-b);
    if (diff < MAX_ABS) return true;
    a = fabs(a); b = fabs(b);
    return diff <= fmax(a, b) * MAX_REL;
}

inline bool is_almost_zero(real a){
    return fabs(a) < MAX_ABS;
}

inline bool is_almost_one(real f){
    real diff = fabs(f-1.);
    if (diff < MAX_ABS) return true;
    f = fabs(f);
    return diff <= fmax(f,1.) * MAX_REL;
}

inline real chop(real v){
    return is_almost_zero(v) ? 0 : v;
}

inline real det2(real a, real b, real c, real d){
    return a * d - b * c;
}

// If t is singular then return false otherwise
// store the inverse of t in inv and return true.
// The inversion algorithm consists in finding the inverse
// determinant and multiplying by the adjugate of t.
//
// Note that t must be of the form:
// t0 t2 t4
// t1 t3 t5
//  0  0  1
bool xform_inverse(real* inv, real* t){
    real det = det2(t[0], t[1], t[2], t[3]);
	if (is_almost_zero(det)) {
		return false;
	}
    real invdet = 1.0 / det;
	inv[0] = t[3] * invdet;
	inv[2] = -t[2] * invdet;
	inv[4] = inv[2] * t[5] - inv[0] * t[4];
	inv[1] = -t[1] * invdet;
	inv[3] = t[0] * invdet;
	inv[5] = inv[1] * t[4] - inv[3] * t[5];
    return true;
}
inline real norm2(real x, real y){
    return x*x + y*y;
}

inline vec3 cross3(real x0, real y0, real w0, real x1, real y1, real w1){
    vec3 v = {y0*w1 - y1*w0, w0*x1 - w1*x0, x0*y1 - x1*y0};
    return v;
}

inline real sign(real v){
    if(v < 0.) return -1.;
    else if(v > 0.) return 1.;
    else return 0.;
}

inline bool write(const char *file, const char *text){
    FILE *f = fopen(file, "w+");
    if (f == NULL)
        return false;
    fprintf(f, text);
    fclose(f);
    return true;
}

#endif //UTIL_H