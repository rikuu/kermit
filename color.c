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
	uint32_t n_cross = 0;
	for (uint32_t e = 0; e < g->n_arc; ++e) {
		uint32_t v = g->arc[e].ul>>33, u = g->arc[e].v>>1;
		if (c[v].c1 <= c[u].c2 && c[u].c1 <= c[v].c2) continue;
		g->arc[e].del = 1, ++n_cross;
	}
	fprintf(stderr, "[M::%s] removed %d color crossing arcs\n", __func__, n_cross);
	if (n_cross) {
		asg_cleanup(g);
		asg_symm(g);
	}
	return n_cross;
}

km_color_t *km_intervalize(km_multicolor_t *colors, size_t n_reads) {
	km_color_t *intervals = (km_color_t*) calloc(n_reads, sizeof(km_color_t));
	for (size_t i = 0; i < n_reads; i++) {
		if (colors[i].n == 0) continue;
		uint64_t min = colors[i].a[0], max = colors[i].a[0];
		for (size_t j = 1; j < colors[i].n; j++) {
			min = (colors[i].a[j] < min) ? colors[i].a[j] : min;
			max = (colors[i].a[j] > max) ? colors[i].a[j] : max;
		}
		intervals[i].c1 = min;
		intervals[i].c2 = max;
	}
	return intervals;
}
