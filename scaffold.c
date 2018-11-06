#include "assert.h"

#include "kvec.h"
#include "ksort.h"

#include "color.h"
#include "scaffold.h"

#define km_utgc_key(a) ((a).c1)
KRADIX_SORT_INIT(utgc, km_utgc_t, km_utgc_key, 8)

double km_estimate_length(km_ugc_t *ugc) {
    const ma_ug_t *ug = ugc->ug;
    const km_utgc_t *utgc = ugc->utgc;
    
    double avg = 0.;
    for (uint32_t i = 0; i < ug->u.n; ++i) {
        const uint32_t l = ug->u.a[utgc[i].utg].len;
        const uint64_t c = (utgc[i].c2 > utgc[i].c1) ? utgc[i].c2 - utgc[i].c1 : utgc[i].c1 - utgc[i].c2;
        avg += (double)l / (double)c;
    }
    return avg / (double)ug->u.n;
}

km_ugc_t* km_ug_color(const ma_ug_t *ug, const km_color_t *colors)
{
    km_utgc_t *utgc = (km_utgc_t *) calloc(ug->u.n, sizeof(km_utgc_t));
	for (uint32_t i = 0; i < ug->u.n; ++i) {
		const ma_utg_t *p = &ug->u.a[i];
        utgc[i].utg = i;

        utgc[i].c1 = colors[p->a[0] >> 33].c1;
        utgc[i].c2 = colors[p->a[p->n - 1] >> 33].c2;
        
        /*uint64_t min, max;
        min = colors[p->a[0] >> 33].c1;
        max = colors[p->a[0] >> 33].c2;
        for (uint32_t j = 1; j < p->n; j++) {
            const km_color_t c = colors[p->a[j] >> 33];
            min = min < c.c1 ? min : c.c1;
            max = max > c.c2 ? max : c.c2;
        }
        utgc[i].c1 = min, utgc[i].c2 = max;*/
	}

    radix_sort_utgc(utgc, utgc + ug->u.n);
    
    km_ugc_t *ugc = (km_ugc_t *) calloc(1, sizeof(km_ugc_t));
    ugc->c = colors;
    ugc->ug = ug;
    ugc->utgc = utgc;

    return ugc;
}

void km_print_ugc(const km_ugc_t *ugc, FILE *fp) {
    const size_t n_utg = ugc->ug->u.n;
    for (size_t i = 0; i < n_utg; i++) {
        const km_utgc_t *p = &ugc->utgc[i];
        fprintf(fp, "#\tutg%.6dl\t%d\t%d\n", p->utg+1, p->c1, p->c2);
    }
}

scaffold_t km_scaffold(const km_ugc_t *ugc, uint64_t max_dist)
{
    const size_t n_utg = ugc->ug->u.n;

    scaffold_t scaffolds = {0, 0, 0};
    size_t i = 0;
    while (i < n_utg) {
        size_t j = i+1;
        while (j < n_utg) {
            const uint64_t d = ugc->utgc[j].c1 - ugc->utgc[j-1].c2;
            if (d > max_dist) break;
            // TODO: fix case with unitigs overlapping
            j++;
        }
        kv_push(size_t, scaffolds, j-1);
        i = j;
    }
    fprintf(stderr, "[M::%s] sorted %zu unitigs to %zu scaffolds\n", __func__, n_utg, scaffolds.n);
    return scaffolds;
}

typedef struct {
	uint32_t node;
    size_t length;
} queue_item_t;

void fill(asg_t *g, km_color_t *color, uint32_t start, uint32_t end, km_color_t path_color)
{
    const uint32_t n_vtx = g->n_seq*2;

    int8_t *visited = (int8_t*) calloc(n_vtx, sizeof(int8_t));
	kvec_t(queue_item_t) queue = {0, 0, 0};
    queue_item_t s = {start, 0};
    kv_push(queue_item_t, queue, s);
    while (queue.n > 0) {
        s = queue.a[--(queue.n)];
        visited[s.node] = 1;

        km_color_t nc = color[s.node];
        fprintf(stderr, "visiting node (%llu, %llu) (%llu, %llu, %zu)\n", nc.c1, nc.c2, path_color.c1, path_color.c2, s.length);

        if (s.node == end) {
            fprintf(stderr, "path found! %llu, %llu, %zu\n", path_color.c1, path_color.c2, s.length);
            return;
        }
        
        // push all neighbors to stack
        uint32_t nv = asg_arc_n(g, s.node);
        asg_arc_t *av = asg_arc_a(g, s.node);
        for (uint32_t i = 0; i < nv; ++i) {
            queue_item_t neighbor = {av[i].v^1, s.length + asg_arc_len(av[i])};
            if (!color_overlap(path_color, color[neighbor.node], 1)) continue;
            if (av[i].del || visited[neighbor.node]) continue;
            kv_push(queue_item_t, queue, neighbor);
        }
    }
	kv_destroy(queue);
    free(visited);
}

void km_scaffold_fill(km_ugc_t *ugc, asg_t *g, scaffold_t scaf, sdict_t *d)
{
    const km_utgc_t *utgc = ugc->utgc;
    const double estimate = km_estimate_length(ugc);

    size_t s = 0;
    for (size_t i = 0; i < scaf.n; i++) {
        fprintf(stderr, "utg%.6lul\t%llu\t%llu\t%u\n",
            utgc[s].utg+1, utgc[s].c1, utgc[s].c2, ugc->ug->u.a[utgc[s].utg].len);

        size_t e = scaf.a[i];
        while (s < e) {
            const ma_utg_t *utg_s = &ugc->ug->u.a[utgc[s].utg];
            const ma_utg_t *utg_e = &ugc->ug->u.a[utgc[s+1].utg];
            uint32_t read_s = utg_s->a[utg_s->n - 1] >> 33;
            uint32_t read_e = utg_e->a[0] >> 33;
            km_color_t color = {ugc->c[read_s].c1, ugc->c[read_e].c2};
            
            fprintf(stderr, "-> utg%.6lul\t%llu\t%llu\t%u\n",
                utgc[s+1].utg+1, utgc[s+1].c1, utgc[s+1].c2, ugc->ug->u.a[utgc[s+1].utg].len);
            fprintf(stderr, "%s\n%s\n", d->seq[read_s].name, d->seq[read_e].name);
            fprintf(stderr, "%llu\n", (color.c2 - color.c1) * estimate);

            fill(g, ugc->c, read_s, read_e, color);
            s++;
        }
        s = e+1;
    }
}
