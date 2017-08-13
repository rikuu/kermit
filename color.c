#include <string.h>
#include <assert.h>
#include <stdio.h>

#include "asg.h"
#include "kvec.h"

#include "cf.h"
#include "color.h"

km_color_t *km_colors_read(const char *fn, sdict_t *d)
{
	color_file_t *fp = cf_open(fn);
	if (fp == 0) {
		fprintf(stderr, "[E::%s] could not open color file %s\n", __func__, fn);
		exit(1);
	}

	km_color_t *colors = (km_color_t*) calloc(d->n_seq, sizeof(km_color_t));
	color_rec_t r;
	size_t ones = 0, twos = 0, tot = 0;
	while (color_read(fp, &r) >= 0) {
		++tot;
		int32_t id = sd_get(d, r.qn);
		if (id < 0) continue;
		colors[id].c1 = r.c1, colors[id].c2 = r.c2;

		// TODO: count unique colors?
		assert(r.c1 != 0);
		if (r.c2 == 0) ones++;
		else twos++;

		assert(strcmp(d->seq[id].name, r.qn) == 0);
	}
	cf_close(fp);
	fprintf(stderr, "[M::%s] read %ld hits with %ld ones, %ld twos\n", __func__, tot, ones, twos);
	return colors;
}

// remove color crossing arcs
int km_cut_cross(asg_t *g, km_color_t *c) {
	uint32_t e, n_cross = 0;
	for (e = 0; e < g->n_arc; ++e) {
		uint32_t v = g->arc[e].ul>>33, u = g->arc[e].v>>1;
		uint16_t cv1 = c[v].c1, cv2 = c[v].c2;
		uint16_t cu1 = c[u].c1, cu2 = c[u].c2;
		if (cv1 == 0 || cu1 == 0 || cv1 == cu1) continue;
		if (cv2 != 0 && cv2 == cu1) continue;
		if ((cv2 != 0 && cu2 != 0) && (cv1 == cu2 || cv2 == cu2)) continue;
		g->arc[e].del = 1, ++n_cross;
	}
	fprintf(stderr, "[M::%s] removed %d color crossing arcs\n", __func__, n_cross);
	if (n_cross) {
		asg_cleanup(g);
		asg_symm(g);
	}
	return n_cross;
}

// filter out reads with >2 colors
km_color_t *km_filter_multi(km_multicolor_t *colors, size_t n_reads) {
	size_t n_filt = 0, ones = 0, twos = 0;
	km_color_t *filtered = (km_color_t*) calloc(n_reads, sizeof(km_color_t));
	for (size_t i = 0; i < n_reads; i++) {
		if (colors[i].n == 1) filtered[i].c1 = colors[i].a[0], ones++;
		else if (colors[i].n == 2)
			filtered[i].c1 = colors[i].a[0], filtered[i].c2 = colors[i].a[1], twos++;
		else n_filt++;
	}
	fprintf(stderr, "[M::%s] filtered %zu reads; %zu ones, %zu twos\n", __func__, n_filt, ones, twos);
	return filtered;
}

// reduce colors by folding unnecessary colors
void km_fold(km_multicolor_t *colors, size_t n_reads, uint16_t n_bins, uint16_t min_coverage) {
	// Count number of reads crossing each bin
	size_t n_cross = 0;
	uint16_t *crossing = (uint16_t*) calloc(n_bins, sizeof(uint16_t));
	for (size_t i = 0; i < n_reads; i++) {
		if (colors[i].n < 3) continue;
		n_cross++;
		for (size_t j = 0; j < colors[i].n; j++) {
			uint16_t color = colors[i].a[j];
			assert(color < n_bins);
			if (crossing[color] == 255) continue;
			crossing[color]++;
		}
	}

	fprintf(stderr, "[M::%s] %zu crossing reads\n", __func__, n_cross);
	if (n_cross == 0) {
		free(crossing);
		return;
	}

	// Find all stretches of crossed bins with high coverage
	// We repurpose crossing here to quickly find the start and end of folds
	size_t n_fold = 0;
	for (uint16_t i = 0; i < n_bins; i++) {
		if (crossing[i] >= min_coverage) {
			uint16_t j = i;
			while (j < n_bins && crossing[j] >= min_coverage) {
				crossing[j] = 0;
				j++;
			}
			if (j-1 - i >= 3) {
				n_fold += j-1 - i;
				crossing[i] = i;
				crossing[i+1] = j-1;
			}
			i = j;
		} else {
			crossing[i] = 0;
		}
	}

	for (size_t i = 0; i < n_reads; i++) {
		for (size_t j = 0; j < colors[i].n; j++) {
			// This assumes all reads that use colors inside a fold are contained
			// TODO: Check for reads being left or right of the fold and use appropriate colors
			uint16_t c = colors[i].a[j];
			if (crossing[c] != 0) {
				colors[i].n = 2, colors[i].m = 2;
				// colors[i].a = realloc(colors[i].a, sizeof(uint16_t) * 2);
				colors[i].a[0] = crossing[c];
				colors[i].a[1] = crossing[c+1];
				break;
			}
		}
	}
	free(crossing);
	fprintf(stderr, "[M::%s] folded %zu colors\n", __func__, n_fold);
}
