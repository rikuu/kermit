#ifndef KM_CF
#define KM_CF

#include <stdint.h>

#include "paf.h"

typedef paf_file_t color_file_t;

typedef struct {
	const char *qn;
	uint16_t c1, c2;
} color_rec_t;

color_file_t *cf_open(const char *fn);
int cf_close(color_file_t *pf);
int color_read(color_file_t *pf, color_rec_t *r);

#endif
