#ifndef KM_MARKER
#define KM_MARKER

#include <stdint.h>

#include "paf.h"

typedef paf_file_t marker_file_t;

typedef struct {
	const char *n;
	uint32_t p;
	uint16_t bin;
} marker_rec_t;

marker_file_t *marker_open(const char *fn);
int marker_close(marker_file_t *pf);
int marker_read(marker_file_t *pf, marker_rec_t *r);

#endif
