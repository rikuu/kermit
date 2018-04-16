#include <string.h>
#include <assert.h>
#include <stdio.h>

#include "asg.h"
#include "kvec.h"
#include "ksort.h"

#include "cf.h"
#include "color.h"

#define km_color_key(a) ((a).c1)
#define km_multicolor_key(a) ((a))
KRADIX_SORT_INIT(color, km_color_t, km_color_key, 8)
KRADIX_SORT_INIT(multicolor, uint64_t, km_multicolor_key, 8)

void km_cf_print(sdict_t *d, km_color_t *c) {
	for (size_t i = 0; i < d->n_seq; ++i) {
		if (c[i].c1 == 0) continue;
		else if (c[i].c2 == 0) printf("%s\t%llu\n", d->seq[i].name, c[i].c1);
		else printf("%s\t%llu\t%llu\n", d->seq[i].name, c[i].c1, c[i].c2);
	}
}

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
		int32_t id = sd_get(d, r.qn);
		if (id < 0 || r.c1 == 0) continue;
		colors[id].c1 = r.c1, colors[id].c2 = r.c2;

		++tot;
		if (r.c2 - r.c1 == 0) ++ones;
		else if (r.c2 - r.c1 == 1) ++twos;
	}
	cf_close(fp);
	fprintf(stderr, "[M::%s] stored %ld colors with %ld ones, %ld twos, %ld multis, %ld uncolored\n",
		__func__, tot, ones, twos, tot-ones-twos, d->n_seq-tot);
	return colors;
}

// remove color crossing arcs
int km_cut_cross(asg_t *g, km_color_t *c) {
	size_t n_cross = 0;
	for (uint32_t e = 0; e < g->n_arc; ++e) {
		uint32_t v = g->arc[e].ul>>33, u = g->arc[e].v>>1;
		if (c[v].c2 >= c[u].c1 && c[u].c2 >= c[v].c1) continue;
		g->arc[e].del = 1, ++n_cross;
	}
	fprintf(stderr, "[M::%s] removed %ld color crossing arcs\n", __func__, n_cross);
	if (n_cross) {
		asg_cleanup(g);
		asg_symm(g);
	}
	return n_cross;
}

// if two color intervals overlap or are consecutive
static inline int overlap(const km_color_t a, const km_color_t b) {
	uint32_t lower = (a.c2 < b.c2) ? a.c2 : b.c2;
	uint32_t upper = (a.c1 > b.c1) ? a.c1 : b.c1;
	return (lower+1 >= upper);
}

static km_color_t merge_colors(km_color_v colors) {
	// TODO: Some heuristic for chimericity
	radix_sort_color(colors.a, colors.a + colors.n);
	km_color_t result = colors.a[0];
	for (size_t i = 1; i < colors.n; i++) {
		if (overlap(result, colors.a[i])) {
			result.c1 = (colors.a[i].c1 < result.c1) ? colors.a[i].c1 : result.c1;
			result.c1 = (colors.a[i].c2 > result.c2) ? colors.a[i].c2 : result.c2;
		}
	}
	return result;
}

// Move multicoloured nodes to interval coloring, i.e. collapse the colors
km_color_t *km_intervalize(km_multicolor_t *colors, size_t n_reads) {
	km_color_t *intervals = (km_color_t*) calloc(n_reads, sizeof(km_color_t));
	for (size_t i = 0; i < n_reads; i++) {
		if (colors[i].n == 0) continue;
		// TODO: Combine with color merging
		radix_sort_multicolor(colors[i].a, colors[i].a + colors[i].n);
		uint64_t min = colors[i].a[0], max = colors[i].a[0];
		for (size_t j = 1; j < colors[i].n; j++) {
			min = (colors[i].a[j] == min-1) ? colors[i].a[j] : min;
			max = (colors[i].a[j] == max+1) ? colors[i].a[j] : max;
		}
		intervals[i].c1 = min;
		intervals[i].c2 = max;
	}
	return intervals;
}

static inline int contains(const km_color_v v, const km_color_t x) {
	for (size_t i = 0; i < v.n; i++) {
		if (v.a[i].c1 == x.c1 && v.a[i].c2 == x.c2)
			return 1;
	}
	return 0;
}

typedef struct {
	uint32_t node;
	int depth;
} queue_t;

// Propagate colors through the graph
void km_propagate(asg_t *g, km_color_t *colors, int max_depth) {
	asg_cleanup(g); // sort and index
	uint32_t n_vtx = g->n_seq * 2;
	size_t n_uncolored = 0, n_colored = 0;
	int *visited = (int*) calloc(n_vtx, sizeof(int));
	km_color_v reachable = {0, 0, 0};
	kvec_t(queue_t) queue = {0, 0, 0};
	for (uint32_t v = 0; v < n_vtx; ++v) {
		if (colored(colors[v])) continue;
		n_uncolored++;
		queue_t s = {v, 0};
		kv_push(queue_t, queue, s);
		while (queue.n > 0) {
			s = queue.a[--(queue.n)];
			visited[s.node] = 1;
			if (s.depth > max_depth) continue;
			if (colored(colors[s.node])) {
				if (!contains(reachable, colors[s.node]))
					kv_push(km_color_t, reachable, colors[s.node]);
			} else {
				// push all neighbors to stack
				uint32_t nv = asg_arc_n(g, s.node);
				asg_arc_t *av = asg_arc_a(g, s.node);
				for (uint32_t j = 0; j < nv; ++j) {
					queue_t neighbor = {av[j].v^1, s.depth+1};
					if (!av[j].del && !visited[neighbor.node])
						kv_push(queue_t, queue, neighbor);
				}
			}
		}
		if (reachable.n == 0) continue;
		km_color_t r = merge_colors(reachable);
		if (colored(r)) {
			colors[v] = r, ++n_colored;
		} else asg_seq_del(g, v);

		memset(visited, 0, g->n_seq);
		reachable.n = 0, queue.n = 0;
	}
	fprintf(stderr, "[M::%s] colored %ld, removed %ld uncolored nodes\n", __func__, n_colored, n_uncolored - n_colored);
	free(visited);
	kv_destroy(reachable);
	kv_destroy(queue);
	if (n_uncolored != n_colored) asg_cleanup(g);
}
