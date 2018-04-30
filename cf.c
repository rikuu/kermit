#include <stdint.h>
#include <stdio.h>
#include <zlib.h>

#include "kseq.h"
#include "paf.h"

#include "cf.h"

KSTREAM_INIT(gzFile, gzread, 0x10000)

color_file_t *cf_open(const char *fn)
{
	gzFile fp = gzopen(fn, "r");
	if (fp == 0) return 0;
	color_file_t *pf = (color_file_t*) calloc(1, sizeof(color_file_t));
	pf->fp = ks_init(fp);
	return pf;
}

int cf_close(color_file_t *pf)
{
	if (pf == 0) return 0;
	free(pf->buf.s);
	kstream_t *ks = (kstream_t*) pf->fp;
	gzclose(ks->f);
	ks_destroy(ks);
	free(pf);
	return 0;
}

static int color_parse(int l, char *s, color_rec_t *pr)
{
	pr->qn=0, pr->c1=0, pr->c2=0;
	char *q, *r;
	int i, t;
	for (i = t = 0, q = s; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) pr->qn = q;
		else if (t == 1) pr->c1 = strtol(q, &r, 10);
		else if (t == 2) pr->c2 = strtol(q, &r, 10);
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 1) return -1;
	return 0;
}

int color_read(color_file_t *pf, color_rec_t *r)
{
	int ret, dret;
file_read_more:
	ret = ks_getuntil((kstream_t*)pf->fp, KS_SEP_LINE, &pf->buf, &dret);
	if (ret < 0) return ret;
	ret = color_parse(pf->buf.l, pf->buf.s, r);
	if (ret < 0) goto file_read_more;
	return ret;
}
