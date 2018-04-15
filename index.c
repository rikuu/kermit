#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "kvec.h"
#include "sdict.h"
#include "paf.h"
#include "index.h"

void km_idx_destroy(km_idx_t *idx)
{
	if (idx == 0) return;
	for (uint32_t i = 0; i < idx->d->n_seq; i++) {
		for (uint32_t j = 0; j < idx->targets[i].n_bins; j++)
			free(idx->targets[i].bins[j].a);
		free(idx->targets[i].bins);
	}
	free(idx->targets);
	sd_destroy(idx->d);
	free(idx);
}

km_hit_v km_pileup(km_idx_t *idx, const char *contig, const uint32_t pos) {
	km_hit_v v = {0,0,0};

	int32_t tid = sd_get(idx->d, contig);
	size_t id = pos / idx->length;
	if (tid < 0 || id >= idx->targets[tid].n_bins) return v;
	km_hit_v *bin = &idx->targets[tid].bins[id];
	for (size_t i = 0; i < bin->n; i++) {
		const km_hit_t h = bin->a[i];
		if (h.s > pos || h.e < pos) continue;
		kv_push(km_hit_t, v, h);
	}
	return v;
}

static inline uint32_t min(const uint32_t a, const uint32_t b) {
	return a < b ? a : b;
}

typedef struct {
	int32_t qn, tn;
	uint32_t s, e, tl;
} km_hit2_t;

km_idx_t *km_build_idx(const char *fn, sdict_t *d, const uint32_t length, const uint32_t max_overhang) {
	km_idx_t *idx = (km_idx_t *) calloc(1, sizeof(km_idx_t));
	idx->length = length;
	idx->d = sd_init();

	paf_file_t *fp = paf_open(fn);
	if (fp == 0) {
		fprintf(stderr, "ERROR: Failed to open \"%s\"\n", fn);
		exit(1);
	}

	// Find longest matches for reads
	kvec_t(km_hit2_t) v = {0,0,0};
	uint32_t tot = 0;
	paf_rec_t r;
	while (paf_read(fp, &r) >= 0) {
		tot++;
		uint32_t tl = r.te - r.ts;
		int32_t id = sd_get(d, r.qn);
		if (id > 0 && d->seq[id].len < tl) continue;
		km_hit2_t *h;
		kv_pushp(km_hit2_t, v, &h);
		h->qn = sd_put(d, r.qn, tl);
		h->tn = sd_put(idx->d, r.tn, r.tl);
		h->tl = d->seq[h->qn].len = tl;
		h->s = r.ts > min(r.qs, max_overhang) ? (r.ts - min(r.qs, max_overhang)) : 0;
    h->e = min(r.tl, r.te + min(r.ql - r.qe, max_overhang));
	}
	paf_close(fp);

	// Allocate bins
	idx->targets = (km_target_t*) calloc(idx->d->n_seq, sizeof(km_target_t));
	for (uint32_t i = 0; i < idx->d->n_seq; i++) {
		if (idx->d->seq[i].len == 0) continue;
		uint32_t n_bins = (idx->d->seq[i].len / length) + 2;
		idx->targets[i].n_bins = n_bins;
		idx->targets[i].bins = (km_hit_v*) calloc(n_bins, sizeof(km_hit_v));
	}

	// Add hits to bins
	uint32_t stored = 0;
	for (size_t i = 0; i < v.n; i++) {
		km_hit2_t h2 = v.a[i];
		if (d->seq[h2.qn].len != h2.tl) {
			// TODO: Squeeze dict; affects qns?
			// d->seq[h2->qn].del = 1;
			continue;
		}

		stored++;
		km_hit_v *bins = idx->targets[h2.tn].bins;
		size_t sid = (size_t) floorf(h2.s / (float) length), eid = (size_t) ceilf(h2.e / (float) length);
		assert(sid < idx->targets[h2.tn].n_bins);
		assert(eid < idx->targets[h2.tn].n_bins);
		assert(sid <= eid);

		km_hit_t *h;
		kv_pushp(km_hit_t, bins[sid], &h);
		h->qn = h2.qn;
		h->s = h2.s;
		h->e = h2.e;

		for (size_t i = sid+1; i < eid+1; i++) {
			kv_push(km_hit_t, bins[i], *h);
		}
	}

	// sd_squeeze(d);
	free(v.a);
	fprintf(stderr, "[M::%s] read %u hits; stored %d hits\n", __func__, tot, stored);
	return idx;
}
