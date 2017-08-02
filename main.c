#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kvec.h"
#include "sys.h"
#include "paf.h"
#include "sdict.h"
#include "miniasm.h"
#include "color.h"

#define MA_VERSION "0.2-r168-dirty"
#define KM_VERSION "0.1"

int main(int argc, char *argv[])
{
	ma_opt_t opt;
	int i, c, stage = 100, no_first = 0, no_second = 0, bi_dir = 1, o_set = 0, no_cont = 0, no_cross = 0;
	sdict_t *d, *excl = 0;
	ma_sub_t *sub = 0;
	ma_hit_t *hit;
	size_t n_hits;
	float cov = 40.0;
	char *fn_reads = 0, *fn_colors = 0, *outfmt = "ug";

	ma_opt_init(&opt);
	while ((c = getopt(argc, argv, "n:m:s:c:S:i:d:g:o:h:I:r:f:e:p:C:12VBRbF:")) >= 0) {
		if (c == 'm') opt.min_match = atoi(optarg);
		else if (c == 'i') opt.min_iden = atof(optarg);
		else if (c == 's') opt.min_span = atoi(optarg);
		else if (c == 'c') opt.min_dp = atoi(optarg);
		else if (c == 'o') opt.min_ovlp = atoi(optarg), o_set = 1;
		else if (c == 'S') stage = atoi(optarg);
		else if (c == 'd') opt.bub_dist = atoi(optarg);
		else if (c == 'g') opt.gap_fuzz = atoi(optarg);
		else if (c == 'h') opt.max_hang = atoi(optarg);
		else if (c == 'I') opt.int_frac = atof(optarg);
		else if (c == 'e') opt.max_ext = atoi(optarg);
		else if (c == 'f') fn_reads = optarg;
		else if (c == 'p') outfmt = optarg;
		else if (c == '1') no_first = 1;
		else if (c == '2') no_second = 1;
		else if (c == 'n') opt.n_rounds = atoi(optarg) - 1;
		else if (c == 'B') bi_dir = 1;
		else if (c == 'b') bi_dir = 0;
		else if (c == 'R') no_cont = 1;
		else if (c == 'F') opt.final_ovlp_drop_ratio = atof(optarg);
		else if (c == 'C') fn_colors = optarg;
		else if (c == 'V') {
			printf("%s\n", KM_VERSION);
			return 0;
		} else if (c == 'r') {
			char *s;
			opt.max_ovlp_drop_ratio = strtod(optarg, &s);
			if (*s == ',') opt.min_ovlp_drop_ratio = strtod(s + 1, &s);
		}
	}
	if (o_set == 0) opt.min_ovlp = opt.min_span;
	if (argc == optind) {
		fprintf(stderr, "Usage: kermit [options] <in.paf>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  Pre-selection:\n");
		fprintf(stderr, "    -R          prefilter clearly contained reads (2-pass required)\n");
		fprintf(stderr, "    -m INT      min match length [%d]\n", opt.min_match);
		fprintf(stderr, "    -i FLOAT    min identity [%.2g]\n", opt.min_iden);
		fprintf(stderr, "    -s INT      min span [%d]\n", opt.min_span);
		fprintf(stderr, "    -c INT      min coverage [%d]\n", opt.min_dp);
		fprintf(stderr, "  Overlap:\n");
		fprintf(stderr, "    -o INT      min overlap [same as -s]\n");
		fprintf(stderr, "    -h INT      max over hang length [%d]\n", opt.max_hang);
		fprintf(stderr, "    -I FLOAT    min end-to-end match ratio [%.2g]\n", opt.int_frac);
		fprintf(stderr, "  Layout:\n");
		fprintf(stderr, "    -g INT      max gap differences between reads for trans-reduction [%d]\n", opt.gap_fuzz);
		fprintf(stderr, "    -d INT      max distance for bubble popping [%d]\n", opt.bub_dist);
		fprintf(stderr, "    -e INT      small unitig threshold [%d]\n", opt.max_ext);
		fprintf(stderr, "    -f FILE     read sequences []\n");
		fprintf(stderr, "    -n INT      rounds of short overlap removal [%d]\n", opt.n_rounds + 1);
		fprintf(stderr, "    -r FLOAT[,FLOAT]\n");
		fprintf(stderr, "                max and min overlap drop ratio [%.2g,%.2g]\n", opt.max_ovlp_drop_ratio, opt.min_ovlp_drop_ratio);
		fprintf(stderr, "    -F FLOAT    aggressive overlap drop ratio in the end [%.2g]\n", opt.final_ovlp_drop_ratio);
		fprintf(stderr, "    -C FILE     vertex colors []\n");
		fprintf(stderr, "  Miscellaneous:\n");
		fprintf(stderr, "    -p STR      output information: sg or ug [%s]\n", outfmt);
		fprintf(stderr, "    -b          both directions of an arc are present in input\n");
		fprintf(stderr, "    -1          skip 1-pass read selection\n");
		fprintf(stderr, "    -2          skip 2-pass read selection\n");
		fprintf(stderr, "    -V          print version number\n");
		return 1;
	}

	sys_init();
	d = sd_init();

	if (no_cont || (fn_colors && no_cross)) {
		fprintf(stderr, "[M::%s] ===> Step 0: removing reads <===\n", __func__);
		if (fn_colors && no_cross) {
			fprintf(stderr, "[M::%s] ===> Step 0.1: removing colored reads <===\n", __func__);
			// TODO: Exclude colors not specified
			// km_color_t *colors = km_colors_read(fn_colors, d);
			// excl = km_hit_prune(sg, colors);
			// free(colors);
		}

		if (no_cont) {
			fprintf(stderr, "[M::%s] ===> Step 0.2: removing contained reads <===\n", __func__);
			excl = ma_hit_no_cont(argv[optind], opt.min_span, opt.min_match, opt.max_hang, opt.int_frac);
		}
	}

	fprintf(stderr, "[M::%s] ===> Step 1: reading read mappings <===\n", __func__);
	hit = ma_hit_read(argv[optind], opt.min_span, opt.min_match, d, &n_hits, bi_dir, excl);
	if (excl) sd_destroy(excl);

	if (!no_first) {
		fprintf(stderr, "[M::%s] ===> Step 2: 1-pass (crude) read selection <===\n", __func__);
		if (stage >= 2) {
			sub = ma_hit_sub(opt.min_dp, opt.min_iden, 0, n_hits, hit, d->n_seq);
			n_hits = ma_hit_cut(sub, opt.min_span, n_hits, hit);
		}
		if (stage >= 3) n_hits = ma_hit_flt(sub, opt.max_hang * 1.5, opt.min_ovlp * .5, n_hits, hit, &cov);
	}

	if (!no_second) {
		fprintf(stderr, "[M::%s] ===> Step 3: 2-pass (fine) read selection <===\n", __func__);
		if (stage >= 4) {
			ma_sub_t *sub2;
			sub2 = ma_hit_sub(opt.min_dp, opt.min_iden, opt.min_span/2, n_hits, hit, d->n_seq);
			n_hits = ma_hit_cut(sub2, opt.min_span, n_hits, hit);
			if (sub != 0) {
				ma_sub_merge(d->n_seq, sub, sub2);
				free(sub2);
			} else {
				sub = sub2;
			}
		}
		if (stage >= 5) n_hits = ma_hit_contained(&opt, d, sub, n_hits, hit);
	}

	hit = (ma_hit_t*)realloc(hit, n_hits * sizeof(ma_hit_t));

	fprintf(stderr, "[M::%s] ===> Step 4: 1-pass graph cleaning <===\n", __func__);
	asg_t *sg = ma_sg_gen(&opt, d, sub, n_hits, hit);
	if (fn_colors && !no_cross) {
		fprintf(stderr, "[M::%s] ===> Step 4.1: reading graph coloring <===\n", __func__);
		km_color_t *colors = km_colors_read(fn_colors, d);

		fprintf(stderr, "[M::%s] ===> Step 4.2: removing color crossing arcs <===\n", __func__);
		km_cut_cross(sg, colors);
		free(colors);
	}

	fprintf(stderr, "[M::%s] ===> Step 5: 2-pass graph cleaning <===\n", __func__);
	if (stage >= 6) {
		fprintf(stderr, "[M::%s] ===> Step 5.1: transitive reduction <===\n", __func__);
		asg_arc_del_trans(sg, opt.gap_fuzz);
	}
	if (stage >= 7) {
		fprintf(stderr, "[M::%s] ===> Step 5.2: initial tip cutting and bubble popping <===\n", __func__);
		asg_cut_tip(sg, opt.max_ext);
		asg_pop_bubble(sg, opt.bub_dist);
	}
	if (stage >= 9) {
		fprintf(stderr, "[M::%s] ===> Step 5.3: cutting short overlaps (%d rounds in total) <===\n", __func__, opt.n_rounds + 1);
		for (i = 0; i <= opt.n_rounds; ++i) {
			float r = opt.min_ovlp_drop_ratio + (opt.max_ovlp_drop_ratio - opt.min_ovlp_drop_ratio) / opt.n_rounds * i;
			if (asg_arc_del_short(sg, r) != 0) {
				asg_cut_tip(sg, opt.max_ext);
				asg_pop_bubble(sg, opt.bub_dist);
			}
		}
	}
	if (stage >= 10) {
		fprintf(stderr, "[M::%s] ===> Step 5.4: removing short internal sequences and bi-loops <===\n", __func__);
		asg_cut_internal(sg, 1);
		asg_cut_biloop(sg, opt.max_ext);
		asg_cut_tip(sg, opt.max_ext);
		asg_pop_bubble(sg, opt.bub_dist);
	}
	if (stage >= 11) {
		fprintf(stderr, "[M::%s] ===> Step 5.5: aggressively cutting short overlaps <===\n", __func__);
		if (asg_arc_del_short(sg, opt.final_ovlp_drop_ratio) != 0) {
			asg_cut_tip(sg, opt.max_ext);
			asg_pop_bubble(sg, opt.bub_dist);
		}
	}

	if (strcmp(outfmt, "ug") == 0) {
		fprintf(stderr, "[M::%s] ===> Step 6: generating unitigs <===\n", __func__);
		ma_ug_t *ug = ma_ug_gen(sg);
		if (fn_reads) ma_ug_seq(ug, d, sub, fn_reads);
		ma_ug_print(ug, d, sub, stdout);
		ma_ug_destroy(ug);
	} else ma_sg_print(sg, d, sub, stdout);

	asg_destroy(sg);
	free(sub); free(hit);
	sd_destroy(d);

	fprintf(stderr, "[M::%s] Version: %s (miniasm %s)\n", __func__, KM_VERSION, MA_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, sys_realtime(), sys_cputime());
	return 0;
}
