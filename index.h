#ifndef KM_INDEX
#define KM_INDEX

#include "sdict.h"

typedef struct {
	int32_t qn;
	uint32_t s, e;
} km_hit_t;

typedef kvec_t(km_hit_t) km_hit_v;

typedef struct {
	uint32_t n_bins;
	km_hit_v *bins;
} km_target_t;

typedef struct  {
	uint32_t length;
	sdict_t *d;
	km_target_t *targets;
} km_idx_t;

void km_idx_destroy(km_idx_t *idx);
km_hit_v km_pileup(km_idx_t *idx, const char *contig, const uint32_t pos);
km_idx_t *km_build_idx(const char *fn, sdict_t *d, const uint32_t length, const uint32_t max_overhang);

#endif
