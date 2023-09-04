// From "How to Solve a Quadratic Equation"
// Jim Blinn's Corner, Nov/Dec 2005
#ifndef QUADRATIC_H
#define QUADRATIC_H

#include <math.h>
#include "vec.h"

#ifdef __cplusplus
extern "C" {
#endif

extern vec4 solve_quadratic_precomputed_delta(int* n, real a, real b, real c, real delta);
extern vec4 solve_quadratic(int * n, real a, real b, real c);

#ifdef __cplusplus
}
#endif

vec4 solve_quadratic_precomputed_delta(int* n, real a, real b, real c, real delta) {
    b *= .5;
    if (delta >= 0) {
        real d = sqrt(delta);
        *n = 2;
        if (b > 0) {
            real e = b + d;
            vec4 v = { -c, e, e, -a };
            return v;
        }
        else if (b < 0.f) {
            real e = -b + d;
            vec4 v = { e, a, c, e };
            return v;
        }
        else if (fabs(a) > fabs(c)) {
            vec4 v = { d, a, -d, a };
            return v;
        }
        else {
            vec4 v = { -c, d, c, d };
            return v;
        }
    }
    else {
        *n = 0;
        vec4 v = { 0, 0, 0, 0 };
        return v;
    }
}

// returns a vec4 with x=t0, y=s0, z=t1 and w=s1
// and stores in n the number of real roots of a*x^2 + b*x + c == 0
// such that each root i in 0...n-1 is given by ti/si
inline vec4 solve_quadratic(int* n, real a, real b, real c) {
    return solve_quadratic_precomputed_delta(n, a, b, c, 0.25f * b * b - a * c);
}

#endif // QUADRATIC_H