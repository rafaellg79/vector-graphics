#ifndef CURVE_H
#define CURVE_H

#include <stdlib.h>

#include "vec.h"
#include "util.h"
#include "bounding_box.h"
#include "bezier.h"

#ifdef __cplusplus
extern "C" {
#endif

enum CurveOrder{
    POINT = 0,
    LINE = 1,
    QUADRATIC = 2,
    CUBIC = 3,
    RQUADRATIC = 4
};

// An implicit curve
typedef struct Curve{
    BoundingBox BB;
    enum CurveOrder type;
    bool asc;
    bool convex;
    vec3 l, m, n;
    real F[9];
    int s;
    real x1, y1, x2, y2;
    real x0, y0, xf, yf;
    struct Curve * sh;
} Curve;

int winding(const Curve* curve, real x, real y);
Curve linear_segment(const real x0, const real y0,
                   const real x1, const real y1);
Curve quadratic_segment(const real x0, const real y0,
                        const real x1, const real y1,
                        const real x2, const real y2);
Curve cubic_segment(const real x0, const real y0,
                    const real x1, const real y1,
                    const real x2, const real y2,
                    const real x3, const real y3);

#ifdef __cplusplus
}
#endif

int winding(const Curve* curve, real x, real y) {
    // outside bounding box and non-intersecting
    if (y <= curve->BB.ymin || y > curve->BB.ymax || x < curve->BB.xmin)
        return 0;
    // outside bounding box and intersecting
    if (x >= curve->BB.xmax)
        return curve->s;
    real lv = curve->l.x * x + curve->l.y * y + curve->l.z;
    switch (curve->type) {
    case LINE:
        if (lv >= 0.)
            return curve->s;
        return 0;

    case QUADRATIC: {
        // triangle test
        if (lv > 0.)
            return curve->convex ? curve->s : 0;
        // implict test
        real dx = x - curve->x0, dy = y - curve->y0;
        real qv = dx * (dx * curve->F[1] + dy * curve->F[2] + curve->F[0])
            + dy * (dy * curve->F[4] + curve->F[3]);
        if (qv >= 0.)
            return curve->s;
        return 0;
    }

    case CUBIC: {
        // triangle test
        if (lv > 0.)
            return curve->convex ? curve->s : 0;
        // tangents test
        if (curve->m.x * x + curve->m.y * y + curve->m.z > 0. ||
            curve->n.x * x + curve->n.y * y + curve->n.z > 0.)
            return (!curve->convex) ? curve->s : 0;
        // implicit test
        real dx = x - curve->x0, dy = y - curve->y0;
        real cv = dx * (dx * (dx * curve->F[8] + dy * curve->F[7] + curve->F[2])
            + dy * (dy * curve->F[6] + curve->F[3]) + curve->F[0])
            + dy * (dy * (dy * curve->F[5] + curve->F[4]) + curve->F[1]);
        if (cv >= 0.)
            return curve->s;
        return 0;
    }

    default:
        // Unknown curve type
        return 0;
    }
}

// Creates a Curve object representing the line segment C that
// starts at point p0 = (x0, y0) and ends at point p1 = (x1, y1).
Curve linear_segment(real x0, real y0, real x1, real y1){
    Curve c = {
        create_BoundingBox(x0, y0, x1, y1), //Bounding box
        LINE,                               //type of bezier curve
        (x0 - x1) * (y0 - y1) > 0.,         //is the curve ascending
    };
    if (y0 <= y1) {
        c.s = 1;
        c.l = cross3(x1, y1, 1., x0, y0, 1.);
    } else {
        c.s = -1;
        c.l = cross3(x0, y0, 1., x1, y1, 1.);
    }
    
    c.x0 = x0; c.y0 = y0;
    c.xf = x1; c.yf = y1;
    c.sh = NULL;
    return c;
}

// Creates a Curve object representing the quadratic bezier curve C given by
// the control points p0 = (x0, y0), p1 = (x1, y1) and p2 = (x2, y2).
//
// Note: Assumes that C is a monotonic curve.
Curve quadratic_segment(real x0, real y0, real x1, real y1, real x2, real y2){
    Curve c = {
        create_BoundingBox(x0, y0, x2, y2), //Bounding box
        QUADRATIC,                          //type of bezier curve
        (x0 - x2) * (y0 - y2) > 0.,         //is the curve ascending
    };
    
    if (y0 <= y2) {
        c.s = 1;
        c.l = cross3(x2, y2, 1., x0, y0, 1.);
    } else {
        c.s = -1;
        c.l = cross3(x0, y0, 1., x2, y2, 1.);
    }
    
    vec2 mid = at2(.5, x0, y0, x1, y1, x2, y2);
    c.convex = c.l.x*mid.x + c.l.y*mid.y + c.l.z < 0.;
    if (!c.convex) {
        c.l.x = -c.l.x;
        c.l.y = -c.l.y;
        c.l.z = -c.l.z;
    }
    
    real dx1 = x1-x0, dy1 = y1-y0;
    real dx2 = x2-x0, dy2 = y2-y0;

    real sgn = 2*dy2*(dx1*dy2 - dx2*dy1);
    real inv = (sgn < 0.) ? -1. : 1.;

    //implicit function coefficients(order: x, xy, xx, y, yy)
    c.F[0] = (4*dy1*(dx1*dy2 - dx2*dy1))*inv;
    c.F[1] = (4*dy1*(dy1 - dy2) + dy2*dy2)*inv;
    c.F[2] = (4*dx1*(dy2 - 2*dy1) + 2*dx2*(2*dy1 - dy2))*inv;
    c.F[3] = (4*dx1*(dx2*dy1 - dx1*dy2))*inv;
    c.F[4] = (4*dx1*(dx1 - dx2) + dx2*dx2)*inv;
    
    c.x0 = x0; c.y0 = y0;
    c.xf = x2; c.yf = y2;
    c.x1 = x1; c.y1 = y1;
    c.sh = NULL;
    return c;
}

// Creates a Curve object representing the cubic bezier curve C given by
// the control points p0 = (x0, y0), p1 = (x1, y1), p2 = (x2, y2) and p3 = (x3, y3).
//
// Note: Assumes that C is a monotonic curve.
Curve cubic_segment(real x0, real y0,
                    real x1, real y1,
                    real x2, real y2,
                    real x3, real y3){
    Curve c = {
        create_BoundingBox(x0, y0, x3, y3), //Bounding box
        CUBIC,                              //type of bezier curve
        (x0 - x3) * (y0 - y3) > 0.,         //is the curve ascending
    };

    // midline
    vec3 l;
    if (y0 <= y3) {
        c.s = 1;
        l = cross3(x3, y3, 1., x0, y0, 1.);
    } else {
        c.s = -1;
        l = cross3(x0, y0, 1., x3, y3, 1.);
    }

    // convexity
    vec2 mid = at3(.5, x0, y0, x1, y1, x2, y2, x3, y3);
    bool convex = l.x*mid.x + l.y*mid.y + l.z < 0.;
    if (!convex) {
        l.x = -l.x;
        l.y = -l.y;
        l.z = -l.z;
    }

    // tangents
    vec3 m, n;
    cubic_tan(&m, &n, x0, y0, x1, y1, x2, y2, x3, y3);
    if ((c.s > 0) == convex) {
        m.x = -m.x; n.x = -n.x;
        m.y = -m.y; n.y = -n.y;
        m.z = -m.z; n.z = -n.z;
    }

    // implicit form coefficients
    real dx1 = x1 - x0, dy1 = y1 - y0;
    real dx2 = x2 - x0, dy2 = y2 - y0;
    real dx3 = x3 - x0, dy3 = y3 - y0;

    real sig = (dy1 - dy2 - dy3)*(-dx3*dx3*(4*dy1*dy1 - 2*dy1*dy2 + dy2*dy2)
               + dx1*dx1*(9*dy2*dy2 - 6*dy2*dy3 - 4*dy3*dy3)
               + dx2*dx2*(9*dy1*dy1 - 12*dy1*dy3 - dy3*dy3)
               + 2*dx1*dx3*(-dy2*(6*dy2 + dy3)
               + dy1*(3*dy2 + 4*dy3))
               - 2*dx2*(dx3*(3*dy1*dy1 - dy2*dy3 + dy1*(-6*dy2 + dy3))
               + dx1*(dy1*(9*dy2 - 3*dy3) - dy3*(6*dy2 + dy3))));

    real inv = (sig < 0.f) ? -1.f : 1.f;

    //implicit function coefficients(order: x, y, xx, xy, yy, yyy, xyy, xxy, xxx)
    c.F[0] = (27*dy1*(dx1*(dy3*(dx1*dy3-2*dx3*dy1-3*dx2*dy2)+3*dx3*dy2*dy2)
           + dy1*(3*dx2*(dx2*dy3-dx3*dy2)+dx3*dx3*dy1)))*inv;

    c.F[1] = (27*dx1*(dx1*(dx3*(2*dy1*dy3-3*dy2*dy2)+dy3*(3*dx2*dy2-dx1*dy3))
           + dy1*(3*dx2*(dx3*dy2 - dx2*dy3) - dx3*dx3*dy1)))*inv;

    c.F[2] = (9*(dy1*(dy3*(dx1*(6*(dy1-dy3)+9*dy2)+dx2*(9*dy2+2*dy3-18*dy1)
           + 3*dx3*(2*dy1-dy2))+9*dy2*(dx3*(dy1-dy2)+dx2*dy1-dx1*dy2)
           - 6*dx3*dy1*dy1)+dy2*(3*dy2*(dx3*dy2-dx2*dy3)+dx1*dy3*dy3)))*inv;

    c.F[3] = (dx1*(dy3*(dx1*(54*dy3-81*dy2-108*dy1)+27*dx2*(9*dy1-3*dy2
           - dy3)+(9*dx3*dy2))+dx3*(dy1*(108*dy1-243*dy2)+(81*dy2*dy2))
           + (81*dx1*dy2*dy2))+dx2*(9*dy1*(dx3*(9*(dy1+dy2)-dy3)-9*dx2*(dy1+dy3))
           + 54*dy2*(dx2*dy3-dx3*dy2))+27*dx3*dx3*dy1*(dy2-2*dy1))*inv;

    c.F[4] = (dx1*(dx3*(54*dx1*(3*dy2-dy1-dy3)+27*dx2*(dy3-3*(dy1+dy2))
           + 18*dx3*(3*dy1-dy2))+81*dx2*(dx2*(dy1+dy3)-dx1*(dy2+dy3))
           + 54*dx1*dx1*dy3) + 9*dx2*(3*dx2*(dx3*dy2-dx2*dy3)-dx3*dx3*dy1))*inv;

    c.F[5] = (dx1*(dx1*(81*dx2 - 27*dx3 - 27*dx1) + dx2*(54*dx3 - 81*dx2)
           - 9*dx3*dx3)+ dx2*(27*dx2*(dx2 - dx3) + 9*dx3*dx3)- dx3*dx3*dx3)*inv;

    c.F[6] = ((3*(dy2-dy1)-dy3)*(dx1*(54*dx2-27*dx1-18*dx3)+dx2*(18*dx3
           - 27*dx2)-3*dx3*dx3))*inv;
    
    c.F[7] = ((3*(dx2-dx1)-dx3)*(dy1*(27*dy1-54*dy2+18*dy3)+dy2*(27*dy2
           - 18*dy3)+3*dy3*dy3))*inv;
    
    c.F[8] = (dy1*(27*dy1*(dy1-3*dy2+dy3)+dy2*(81*dy2-54*dy3)+9*dy3*dy3)
           + dy2*(27*dy2*(dy3-dy2)-9*dy3*dy3)+dy3*dy3*dy3)*inv;

    c.l = l;
    c.m = m;
    c.n = n;
    c.convex = convex;
    c.x0 = x0; c.y0 = y0;
    c.xf = x3; c.yf = y3;
    c.x1 = x1; c.y1 = y1;
    c.x2 = x2; c.y2 = y2;
    c.sh = NULL;
    return c;
}

#endif //CURVE_H