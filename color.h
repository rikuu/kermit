#ifndef KM_COLOR
#define KM_COLOR

#include <stdint.h>
#include "kvec.h"
#include "sdict.h"
#include "asg.h"

typedef struct {
	uint16_t c1, c2;
} km_color_t;

// for representing multi-coloring for folding
typedef kvec_t(uint16_t) km_multicolor_t;

km_color_t *km_colors_read(const char *fn, sdict_t *d);
int km_cut_cross(asg_t *g, km_color_t *c);

void km_fold(km_multicolor_t *colors, size_t n_reads, uint16_t n_bins, uint16_t min_coverage);
km_color_t *km_filter_multi(km_multicolor_t *colors, size_t n_reads);

#endif
