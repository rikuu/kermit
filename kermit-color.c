#include <stdint.h>
#include <sys/types.h>
#include <zlib.h>
#include <stdio.h>
#include <string.h>

#include "sys.h"

#include "color.h"
#include "cf.h"
#include "index.h"
#include "marker.h"

static void print_colors(const sdict_t *d, const km_color_t *c)
{
	for (size_t i = 0; i < d->n_seq; ++i) {
		if (c[i].c1 == 0) continue;
		else if (c[i].c2 == 0) printf("%s\t%d\n", d->seq[i].name, c[i].c1);
		else printf("%s\t%d\t%d\n", d->seq[i].name, c[i].c1, c[i].c2);
	}
}

static void print_multicolors(const sdict_t *d, const km_multicolor_t *c) {
	for (size_t i = 0; i < d->n_seq; i++) {
		if (c[i].n == 0) continue;
		printf("%s", d->seq[i].name);
		for (size_t j = 0; j < c[i].n; j++) {
			printf("\t%u", c[i].a[j]);
		}
		printf("\n");
	}
}

km_multicolor_t *km_align_markers(char **map_fns, uint16_t n_maps, km_idx_t *idx, size_t n_reads, uint16_t *n_bins)
{
	uint32_t n_markers=0;
	km_multicolor_t *colors = calloc(n_reads, sizeof(km_multicolor_t));
	for (uint16_t i = 0; i < n_maps; i++) {
		marker_file_t *fp = marker_open(map_fns[i]);
		if (fp == 0) {
			fprintf(stderr, "ERROR: Failed to open \"%s\"\n", map_fns[i]);
			exit(1);
		}

		marker_rec_t r;
		while (marker_read(fp, &r) >= 0) {
			n_markers++;
			const uint16_t color = *n_bins = i << 8 | (r.bin + 1);
			const km_hit_v h = km_pileup(idx, r.n, r.p);
			for (size_t j = 0; j < h.n; j++) {
				size_t id = h.a[j].qn;
				if (id >= n_reads || (colors[id].n > 0 && colors[id].a[colors[id].n - 1] == color)) continue;
				kv_push(uint16_t, colors[id], color);
			}
		}
		marker_close(fp);
	}

	size_t ones=0, twos=0, multis=0;
	for (size_t i = 0; i < n_reads; i++) {
		if (colors[i].n == 1) ones++;
		else if (colors[i].n == 2) twos++;
		else if (colors[i].n > 2) multis++;
	}
	fprintf(stderr, "[M::%s] aligned %u markers; %zu ones, %zu twos, %zu multis\n", __func__, n_markers, ones, twos, multis);
	return colors;
}

int main(int argc, char *argv[])
{
	char *paf_fn = 0, **map_fns = 0;
	uint16_t n_maps, min_coverage = 3, bin_length = 16;
	uint32_t max_overhang = 250;
	int no_fold = 0, no_merge = 1, no_colors = 0, c;

	while ((c = getopt(argc, argv, "c:o:l:pfmV")) >= 0) {
		if (c == 'c') min_coverage = atoi(optarg);
		else if (c == 'o') max_overhang = atoi(optarg);
		else if (c == 'l') bin_length = atoi(optarg);
		else if (c == 'f') no_fold = 1;
		else if (c == 'm') no_merge = 0;
		else if (c == 'p') no_colors = 1;
		else if (c == 'V') {
			// printf("%s\n", KM_VERSION);
			return 0;
		}
	}

	if (argc - optind >= 2) {
		paf_fn = argv[optind];

		n_maps = argc - optind - 1;
		map_fns = malloc(sizeof(char*) * n_maps);
		memcpy(map_fns, argv + (optind + 1), sizeof(char*) * n_maps);
	} else {
		fprintf(stderr, "Usage: kermit-color [options] <in.paf> <markers> [<markers2>,..]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, " -c INT      min fold coverage [%d]\n", min_coverage);
		fprintf(stderr, " -f          disable color folding\n");
		// fprintf(stderr, " -p          output multicolored\n");
		// fprintf(stderr, " -m          merge similar colors\n");
		fprintf(stderr, " -V          print version number\n");
		return 1;
	}

	sys_init();
	sdict_t *d = sd_init();

	fprintf(stderr, "[M::%s] ===> Step 1: indexing read mappings <===\n", __func__);
	km_idx_t *idx = km_build_idx(paf_fn, d, bin_length, max_overhang);

	// TODO: Support reference-only coloring
	fprintf(stderr, "[M::%s] ===> Step 2: mapping markers to reads <===\n", __func__);
	size_t n_reads = d->n_seq; uint16_t n_bins = 0;
	km_multicolor_t *multicolors = km_align_markers(map_fns, n_maps, idx, n_reads, &n_bins);
	km_idx_destroy(idx);

	fprintf(stderr, "[M::%s] %zu reads mapped to %u bins\n", __func__, n_reads, n_bins);

	if (!no_colors) {
		if (!no_fold) {
			fprintf(stderr, "[M::%s] ===> Step 3: folding colors <===\n", __func__);
			km_fold(multicolors, n_reads, n_bins, min_coverage);
		}
		km_color_t *colors = km_filter_multi(multicolors, n_reads);
		print_colors(d, colors);
		free(colors);
	} else {
		print_multicolors(d, multicolors);
	}

	sd_destroy(d);
	free(multicolors);

	//fprintf(stderr, "[M::%s] Version: %s\n", __func__, KM_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (int i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, sys_realtime(), sys_cputime());
	return 0;
}
