#ifndef KM_COLOR
#define KM_COLOR

#include <stdint.h>

#include "kvec.h"
#include "sdict.h"
#include "asg.h"

typedef struct {
	uint64_t c1, c2;
} km_color_t;

// Queue structure for propagation
typedef struct {
	uint32_t node;
	int depth;
} queue_item_t;

typedef kvec_t(km_color_t) km_color_v;
typedef kvec_t(uint64_t) km_multicolor_t;

#define COLORED(n) (n.c1 != 0 || n.c2 != 0)

void km_cf_print(sdict_t *d, km_color_t *c);
km_color_t *km_colors_read(const char *fn, sdict_t *d);
sdict_t *km_exclude(const char *fn, uint64_t color);

int km_cut_cross(asg_t *g, km_color_t *c);

km_color_t *km_intervalize(km_multicolor_t *colors, size_t n_reads);
void km_propagate(asg_t *g, km_color_t *colors, int max_depth);

#endif
