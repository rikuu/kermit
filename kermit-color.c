#include <stdint.h>
#include <sys/types.h>
#include <zlib.h>
#include <stdio.h>
#include <string.h>

#include "sys.h"

#include "kermit.h"
#include "color.h"
#include "cf.h"
#include "index.h"
#include "marker.h"

static void color_stats(const km_multicolor_t *colors, const size_t n_reads)
{
	size_t uncolored=0, ones=0, twos=0, multis=0;
	for (size_t i = 0; i < n_reads; i++) {
		if (colors[i].n == 0) uncolored++;
		else if (colors[i].n == 1) ones++;
		else if (colors[i].n == 2) twos++;
		else multis++;
	}
	fprintf(stderr, "[M::%s] %zu ones, %zu twos, %zu multis, %zu uncolored\n", __func__, ones, twos, multis, uncolored);
}

km_multicolor_t *km_align_markers(char **map_fns, size_t n_maps, km_idx_t *idx, size_t n_reads, int unicolor)
{
	size_t n_markers = 0;
	km_multicolor_t *colors = (km_multicolor_t*) calloc(n_reads, sizeof(km_multicolor_t));
	for (size_t i = 0; i < n_maps; i++) {
		marker_file_t *fp = marker_open(map_fns[i]);
		if (fp == 0) {
			fprintf(stderr, "[E::%s] could not open map file %s\n", __func__, map_fns[i]);
			exit(1);
		}

		marker_rec_t r;
		while (marker_read(fp, &r) >= 0) {
			n_markers++;
			const uint64_t color = unicolor ? (i << 16) : (i << 16 | (r.bin + 1));
			km_hit_v h = km_pileup(idx, r.n, r.p);
			for (size_t j = 0; j < h.n; j++) {
				int32_t id = h.a[j].qn;
				if (id == -1 || id >= (int32_t) n_reads || (colors[id].n > 0 && colors[id].a[colors[id].n - 1] == color)) continue;
				kv_push(uint64_t, colors[id], color);
			}
			free(h.a);
		}
		marker_close(fp);
	}
	fprintf(stderr, "[M::%s] aligned %zu markers\n", __func__, n_markers);
	color_stats(colors, n_reads);
	return colors;
}

km_multicolor_t *km_align_reference(km_idx_t *idx, size_t n_reads, int unicolor)
{
	km_multicolor_t *colors = (km_multicolor_t*) calloc(n_reads, sizeof(km_multicolor_t));
	for (size_t i = 0; i < idx->d->n_seq; i++) {
		km_target_t *target = &idx->targets[i];
		for (size_t j = 0; j < target->n_bins; j++) {
			const uint64_t color = unicolor ? (i << 16) : (i << 16 | (j + 1));
			for (size_t k = 0; k < target->bins[j].n; k++) {
				int32_t id = target->bins[j].a[k].qn;
				if (id == -1 || id >= (int32_t) n_reads || (colors[id].n > 0 && colors[id].a[colors[id].n - 1] == color)) continue;
				kv_push(uint64_t, colors[id], color);
			}
		}
	}
	color_stats(colors, n_reads);
	return colors;
}

int main(int argc, char *argv[])
{
	char *paf_fn = 0, **map_fns = 0;
	uint16_t n_maps = 0;
	uint32_t max_overhang = 250, bin_length = 10000;
	int reference_only = 1, unicolor = 0, c;

	while ((c = getopt(argc, argv, "o:l:Vu")) >= 0) {
		if (c == 'o') max_overhang = atoi(optarg);
		else if (c == 'l') bin_length = atoi(optarg);
		else if (c == 'u') unicolor = 1;
		else if (c == 'V') {
			printf("%s\n", KM_VERSION);
			return 0;
		}
	}

	if (argc - optind >= 1) {
		paf_fn = argv[optind];

		if (argc - optind > 1) {
			reference_only = 0;
			n_maps = argc - optind - 1;
			map_fns = (char**) malloc(sizeof(char*) * n_maps);
			memcpy(map_fns, argv + (optind + 1), sizeof(char*) * n_maps);
		}
	} else {
		fprintf(stderr, "Usage: kermit-color [options] <in.paf> <markers> [<markers2>,..]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, " -V          print version number and exit\n");
		fprintf(stderr, " -o [INT]    max extension of mapping\n");
		fprintf(stderr, " -l [INT]    index bin length\n");
		fprintf(stderr, " -u          assign single color per contig\n");
		return 1;
	}

	sys_init();
	sdict_t *d = sd_init();

	fprintf(stderr, "[M::%s] ===> Step 1: indexing read mappings <===\n", __func__);
	km_idx_t *idx = km_build_idx(paf_fn, d, bin_length, max_overhang);

	km_multicolor_t *multicolors = 0;
	if (!reference_only) {
		fprintf(stderr, "[M::%s] ===> Step 2: mapping markers to reads <===\n", __func__);
		multicolors = km_align_markers(map_fns, n_maps, idx, d->n_seq, unicolor);
		free(map_fns);
	} else {
		fprintf(stderr, "[M::%s] ===> Step 2: reversing index <===\n", __func__);
		multicolors = km_align_reference(idx, d->n_seq, unicolor);
	}

	km_idx_destroy(idx);

	km_color_t *colors = km_intervalize(multicolors, d->n_seq);
	km_cf_print(d, colors);
	free(colors);

	for (size_t i = 0; i < d->n_seq; i++)
		free(multicolors[i].a);
	free(multicolors);

	sd_destroy(d);

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, KM_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (int i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, sys_realtime(), sys_cputime());
	return 0;
}
