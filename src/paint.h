#ifndef PAINT_H
#define PAINT_H

#include <stdlib.h>

#include "vec.h"
#include "util.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef real (*spread_function)(real);

enum ColorType { PAINT_NONE, PAINT_COLOR, PAINT_LINEAR, PAINT_RADIAL };

typedef struct RampStop{
	real offset;
	vec4 color;
} RampStop;

typedef struct Ramp{
	RampStop * stops;
    int nstops;
} Ramp;

typedef struct LinearGradient{
    real a, b, c;
    spread_function spread;
    Ramp ramp;
} LinearGradient;

typedef struct RadialGradient{
    real xx, xy, x, yy, y, w;
    real fx, fy;
    spread_function spread;
    Ramp ramp;
} RadialGradient;

typedef struct Paint{
    char type;
    real opacity;
    union {
        vec4 color;
        LinearGradient linear;
        RadialGradient radial;
    };
} Paint;

real pad(real t);
real reflect(real t);
real repeat(real t);
LinearGradient linear_gradient(const real x1, const real y1,
                               const real x2, const real y2,
                               spread_function spread,
                               Ramp ramp,
                               real * xform);
RadialGradient radial_gradient(const real fx, const real fy,
                               real * T,
                               spread_function spread,
                               Ramp ramp);

#ifdef __cplusplus
}
#endif

real pad(real t){
    if(t<0.0)
        return 0.0;
    else if(t>1.0)
        return 1.0;
    else
        return t;
}

real reflect(real t){
    return 2 * fabs(0.5*t - floor(0.5*t+0.5));
}

real repeat(real t){
    return t - floor(t);
}

LinearGradient linear_gradient(real x1, real y1,
                               real x2, real y2,
                               spread_function spread, Ramp ramp,
                               real * xform){
    real xformed_x1 = x1*xform[0] + y1*xform[2] + xform[4];
    real xformed_y1 = x1*xform[1] + y1*xform[3] + xform[5];
    real xformed_x2 = x2*xform[0] + y2*xform[2] + xform[4];
    real xformed_y2 = x2*xform[1] + y2*xform[3] + xform[5];
    
    real nx = (y1 - y2)*xform[0] + (x2 - x1)*xform[2];
    real ny = (y1 - y2)*xform[1] + (x2 - x1)*xform[3];
    real dx = xformed_x2 - xformed_x1;
    real dy = xformed_y2 - xformed_y1;
    
    real T[] = {dx, nx, dy, ny, 0, 0};
    real inv[6];
    if(!xform_inverse(inv, T)){
		// ERROR: No inverse
        fprintf(stderr, "LinearGradient without inverse!");
        LinearGradient gradient;
        memset(&gradient, 0, sizeof(LinearGradient));
        return gradient;
    }
	
    LinearGradient gradient = {
        inv[0],
        inv[2],
        -inv[0]*x1-inv[2]*y1+inv[4],
        spread,
        ramp
    };
    return gradient;
}

// Creates a radial gradient with focus (fx, fy) in unit circle,
// given a transform T from unit circle to viewport space,
// a color ramp and a spread function
RadialGradient radial_gradient(real fx, real fy,
                               real * T,
                               spread_function spread,
                               Ramp ramp){
    real Tfx = T[0]*fx+T[2]*fy+T[4];
    real Tfy = T[1]*fx+T[3]*fy+T[5];
    
    // Q is a symmetric 3x3 matrix stored in the form
    // Q[0] Q[1] Q[3]
    // Q[1] Q[2] Q[4]
    // Q[3] Q[4] Q[5]
    real Q[6];
    Q[0] = T[0]*T[0] + T[1]*T[1];
    Q[1] = T[0]*T[2] + T[1]*T[3];
    Q[2] = T[2]*T[2] + T[3]*T[3];
    Q[3] = T[0]*Tfx  + T[1]*Tfy;
    Q[4] = T[2]*Tfx  + T[3]*Tfy;
    Q[5] = Tfx *Tfx  + Tfy *Tfy - 1;
    
    RadialGradient gradient;
    // Quadric equation (x'Qx)
    gradient.xx =   Q[0];
    gradient.xy = 2*Q[1];
    gradient.x  = 2*Q[3];
    gradient.yy =   Q[2];
    gradient.y  = 2*Q[4];
    gradient.w  =   Q[5];
    // Focus
    gradient.fx = fx;
    gradient.fy = fy;
    // Spread function
    gradient.spread = spread;
    // Ramp
    gradient.ramp = ramp;
    return gradient;
}

#endif //PAINT_H