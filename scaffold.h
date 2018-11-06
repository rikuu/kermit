#include "miniasm.h"

#include "color.h"

typedef struct {
    uint64_t c1, c2;
    size_t utg;
} km_utgc_t;

typedef struct {
    ma_ug_t *ug;
    km_utgc_t *utgc;
    km_color_t *c;
} km_ugc_t;

inline void km_ugc_destroy(km_ugc_t *ugc) {
    free(ugc->utgc);
}

// TODO: Merge to ugc
typedef struct { size_t n, m; size_t *a; } scaffold_t;

km_ugc_t* km_ug_color(const ma_ug_t *ug, const km_color_t *colors);
scaffold_t km_scaffold(const km_ugc_t *ugc, uint64_t max_dist);
void km_scaffold_fill(km_ugc_t *ugc, asg_t *g, scaffold_t scaf, sdict_t *d);

void km_print_ugc(const km_ugc_t *ugc, FILE* fp);
