#ifndef KM_COLOR
#define KM_COLOR

#include <stdint.h>
#include "kvec.h"
#include "sdict.h"
#include "asg.h"

typedef struct {
	uint64_t c1, c2;
} km_color_t;

// for representing multi-coloring for folding
typedef kvec_t(uint64_t) km_multicolor_t;

km_color_t *km_colors_read(const char *fn, sdict_t *d);
int km_cut_cross(asg_t *g, km_color_t *c);

km_color_t *km_intervalize(km_multicolor_t *colors, size_t n_reads);

#endif
