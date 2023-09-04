#ifndef EXPORTER_H
#define EXPORTER_H

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "util.h"
#include "path.h"
#include "shortcut_tree.h"

typedef struct string {
    char* data;
    unsigned int size, capacity;
}string;

int realloc_to_fit(string* str, unsigned int new_size){
    while (str->capacity <= new_size) {
        str->capacity <<= 1;
    }
    str->data = (char*)realloc(str->data, str->capacity * sizeof(char));
    if (str->data == NULL){
        str->capacity = 0;
        return 1;
    }
    return 0;
}

int append(string* dst, char* src){
    int src_size = strlen(src);
    if (realloc_to_fit(dst, dst->size + src_size))
        return 1;
    strcpy(dst->data + dst->size, src);
    dst->size += src_size;
    return 0;
}

int curve_to_svg(Curve c, char *str){
    int n = 0;
    switch (c.type) {
        case LINE:
            n = sprintf(str, "L%f %f ", chop(c.xf), chop(c.yf));
            break;
        case QUADRATIC:
            n = sprintf(str, "Q%f %f %f %f ", chop(c.x1), chop(c.y1), chop(c.xf), chop(c.yf));
            break;
        case CUBIC:
            n = sprintf(str, "C%f %f %f %f %f %f ", chop(c.x1), chop(c.y1), chop(c.x2), chop(c.y2), chop(c.xf), chop(c.yf));
            break;
    }
    return n;
}

int element_to_svg(Element e, string * path){
    if(e.size == 0)
        return 0;
    Curve ** curves = e.curves;
    if (realloc_to_fit(path, path->size + (e.size + 1) * 150)) {
        return 1;
    }
    append(path, "<path class=\"curve\" marker-start=\"url(#curve_arrow)\" marker-mid=\"url(#curve_arrow)\" marker-end=\"url(#curve_arrow)\" d=\"");

    int i;
    for(i = 0; i < e.size; i++){
        Curve c = (*curves[i]);

        path->size += sprintf(path->data + path->size, "M%f %f ", c.x0, c.y0);
        path->size += curve_to_svg(c, path->data + path->size);
    }
    return append(path, "\"/>\n");
}

int shortcuts_to_svg(Element e, string * path){
    if(e.size == 0)
        return 0;
    Curve ** curves = e.curves;
    if (realloc_to_fit(path, path->size + (e.size + 1) * 150) ||
        append(path, "<path class=\"shortcut\" marker-end=\"url(#shortcut_arrow)\" d=\"")){
        return 1;
    }

    int i;
    for(i = 0; i < e.size; i++){
        if(curves[i]->sh){
            Curve c = *(curves[i]->sh);

            path->size += sprintf(path->data + path->size, "M%f %f ", c.x0, c.y0);
            path->size += curve_to_svg(c, path->data + path->size);
        }
    }
    return append(path, "\"/>\n");
}

int cell_to_svg_rec(ShortcutNode *cell, string *buffer, int id, int recursion_level){
    BoundingBox BB = cell->BB;
    int i;
    if (buffer->capacity - buffer->size < 200) {
        if (buffer->capacity >= 256) buffer->capacity <<= 1;
        else buffer->capacity = 512;
        buffer->data = (char*)realloc(buffer->data, buffer->capacity * sizeof(char));
    }
    buffer->size += sprintf(buffer->data+buffer->size, "<rect class=\"cell\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\"/>\n", BB.xmin, BB.ymin, fabs(BB.xmax - BB.xmin), fabs(BB.ymax - BB.ymin));
    if(cell->cells == NULL){
        double cell_size = (BB.xmax-BB.xmin);
        buffer->size += sprintf(buffer->data + buffer->size, "<text x=\"%f\" y=\"%f\" font-size=\"%fpx\">%d, %d</text>\n", BB.xmin + cell_size * 0.1, BB.ymin + cell_size * 0.3, cell_size * 0.2, id, cell->total_segments);
    }
    if(recursion_level && cell->cells)
        for(i = 0; i < 4; i++){
            if(cell_to_svg_rec(cell->cells+i, buffer, (id << 2) + i + 1, recursion_level-1)){
                return 1;
            }
        }
    else if (cell->elements) {
        for (i = 0; i < cell->size; i++) {
            if(element_to_svg(cell->elements[i], buffer)){
                return 1;
            }
        }
        for (i = 0; i < cell->size; i++) {
            if(shortcuts_to_svg(cell->elements[i], buffer)){
                return 1;
            }
        }
    }
    return 0;
}

char * cell_to_svg(ShortcutNode *cell, int id, int recursion_level){
    if (cell == NULL)
        return "";
    char buffer[256];
    
    BoundingBox BB = cell->BB;
    string svg;
    svg.data = (char *)malloc(sizeof(char) * 2048);
    svg.capacity = 2048;
    svg.size = sprintf(svg.data, "\
<?xml version=\"1.0\" standalone=\"no\"?>\n\
<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20010904//EN\"\n\
    \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n\
<svg viewBox=\"%f %f %f %f\" xmlns=\"http://www.w3.org/2000/svg\">\n\
  <style>\n\
    *{stroke-width:%f}\n\
    rect.cell { stroke: #000000; fill: transparent;}\n\
    path.linear { stroke: #FFFFB3; fill: transparent;}\n\
    path.quadratic { stroke: #00B3FF; fill: transparent;}\n\
    path.cubic { stroke: #FFB300; fill: transparent;}\n\
    path.curve { stroke: #FFB300; fill: transparent;}\n\
    path.shortcut { stroke: #336600; fill: transparent;}\n\
    text { fill: #FF0000; }\n\
  </style>\n\
  <defs>\n\
    <marker style=\"stroke: #000000; fill: #FFB300\" id=\"curve_arrow\" markerWidth = \"10\" markerHeight = \"10\" refY = \"5\" orient = \"auto\">\n\
      <polygon points = \"0 0 10 5 0 10\"/>\n\
    </marker>\n\
    <marker style=\"stroke: #000000; fill: #336600\" id=\"shortcut_arrow\" markerWidth = \"10\" markerHeight = \"10\" refX=\"50%%\" refY = \"5\" orient = \"auto\">\n\
      <polygon points = \"0 0 10 5 0 10\"/>\n\
    </marker>\n\
  </defs>\n", BB.xmin, BB.ymin, BB.xmax - BB.xmin, BB.ymax - BB.ymin, fmin(10.0, BB.xmax - BB.xmin) / 20.0);
    
    if(cell_to_svg_rec(get(cell, id, recursion_level), &svg, id, recursion_level)){
        free(svg.data);
        return NULL;
    }
    append(&svg, "</svg>\n");
    return svg.data;
}

#endif //EXPORTER_H