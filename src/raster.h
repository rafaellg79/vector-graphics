#ifndef RASTER_H
#define RASTER_H

#include <math.h>

#include "util.h"
#include "quadratic.h"
#include "shortcut_tree.h"

#define GAMMA 2.2
#define INVGAMMA (1. / 2.2)

#ifdef __cplusplus
extern "C" {
#endif

bool nonzero(int wn);
bool evenodd(int wn);
vec4 sample(ShortcutNode *cell, real x, real y, vec4 background_color);
vec4 supersample(ShortcutNode *cell, real* pattern, int n, real x, real y, vec4 background_color);

#ifdef __cplusplus
}
#endif

bool nonzero(int wn){
    return wn != 0;
}

bool evenodd(int wn){
    return wn % 2;
}

// blend two colors based on alpha values
vec4 blend(vec4 c1, vec4 c2){
    c2.w = (1 - c1.w) * c2.w;
    
    c1.x = c1.x + c2.x * c2.w;
    c1.y = c1.y + c2.y * c2.w;
    c1.z = c1.z + c2.z * c2.w;
    c1.w = c1.w + c2.w;
    return c1;
}

vec4 ramp_color(Ramp ramp, real t){
    RampStop *stops = ramp.stops;
    if(t <= stops[0].offset){
        return stops[0].color;
    } else if(t >= stops[ramp.nstops-1].offset){
        return stops[ramp.nstops-1].color;
    }
    
    int i;
    for(i = 1; i < ramp.nstops; i++){
        if(stops[i].offset >= t){
            real t1 = stops[i-1].offset;
            real t2 = stops[i].offset;
            vec4 c1 = stops[i-1].color;
            vec4 c2 = stops[i].color;
            
            if(is_almost_equal(t1, t2)){
                return c1;
            }
            
            real m = (t2 - t) / (t2 - t1);
            c1.x = lerp1(m, c2.x, c1.x);
            c1.y = lerp1(m, c2.y, c1.y);
            c1.z = lerp1(m, c2.z, c1.z);
            c1.w = lerp1(m, c2.w, c1.w);
            return c1;
        }
    }
    return stops[ramp.nstops-1].color;
}

vec4 color_at(Paint paint, real x, real y){
    vec4 color = {0.0, 0.0, 0.0, 0.0};
    switch(paint.type){
        case PAINT_COLOR:
            color = paint.color;
            color.w *= paint.opacity;
            return color;
        case PAINT_LINEAR:{
            LinearGradient linear = paint.linear;
            real t = linear.spread(linear.a * x + linear.b * y + linear.c);
            color = ramp_color(linear.ramp, t);
            color.w *= paint.opacity;
            return color;
        }
        case PAINT_RADIAL:{
            RadialGradient radial = paint.radial;
            real dx = x - radial.fx;
            real dy = y - radial.fy;
            int n;
            
            real a = radial.w;
            real b = dx*radial.x + dy*radial.y;
            real c = dx*(dx*radial.xx + dy*radial.xy) + dy*dy*radial.yy;
            if( is_almost_zero(a) && b > 0 && !is_almost_zero(b))
                return color;
            
            vec4 sol = solve_quadratic(&n, a, b, c);
            
            real t;
            if(n > 0 && sol.x * sol.y > 0){
                t = sol.x / sol.y;
            } else if(n > 0 && sol.z * sol.w > 0){
                t = sol.z / sol.w;
            } else if(is_almost_zero(dx) && is_almost_zero(dy)){
                t = 0;
            } else{
                return color;
            }
            
            t = paint.radial.spread(t);
            color = ramp_color(radial.ramp, t);
            color.w *= paint.opacity;
            return color;
        }
        default:
            return color;
    }
}

vec4 sample(ShortcutNode *cell, real x, real y, vec4 background_color){
    // Find subcell containing (x, y)
    ShortcutNode* sub = cell->cells;
    
    while(sub){
        bool top = (y > sub->BB.ymax);
        bool right = (x > sub->BB.xmax);
        cell = sub + (2 * top + right);
        sub = cell->cells;
    }
    
    vec4 color = {0.0, 0.0, 0.0, 0.0};
    for (int i = cell->size-1; i >= 0; i--) {
        //fill test
        Element * element = cell->elements + i;
        Paint paint = element->shape->paint;
        if ((paint.type != PAINT_NONE) && element->shape->winding_test(inside_element(element, x, y))) {
            color = blend(color, color_at(paint, x, y));
            if (is_almost_one(color.w)) {
                color.w = 1.;
                return color;
            }
        }
        //TODO: stroke test
    }
    return blend(color, background_color);
}

// Supersample pixel
vec4 supersample(ShortcutNode *cell, real* pattern, int n, real x, real y, vec4 background_color){
    // Find subcell containing (x, y)
    ShortcutNode* sub = cell->cells;
    
    while(sub && ((sub->BB.xmax - sub->BB.xmin) >= 1.0)){
        bool top = (y > sub->BB.ymax);
        bool right = (x > sub->BB.xmax);
        cell = sub + (2 * top + right);
        sub = cell->cells;
    }
    
    // Blend sample colors
    real inv_n = 1. / n;
    vec4 final_color = {0.0, 0.0, 0.0, 0.0};
    int i;
    for(i = 0; i < 2*n; i+=2){
        real dx = pattern[i], dy = pattern[i+1];
        vec4 sample_color = sample(cell, x+dx, y+dy, background_color);
        final_color.x += pow(sample_color.x, GAMMA);
        final_color.y += pow(sample_color.y, GAMMA);
        final_color.z += pow(sample_color.z, GAMMA);
        final_color.w += sample_color.w;
    }
    
    final_color.x = pow(final_color.x*inv_n, INVGAMMA);
    final_color.y = pow(final_color.y*inv_n, INVGAMMA);
    final_color.z = pow(final_color.z*inv_n, INVGAMMA);
    final_color.w *= inv_n;
    
    return final_color;
}

#endif //RASTER_H