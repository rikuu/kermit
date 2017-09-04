#include <stdio.h>
#include <zlib.h>

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)
#include "marker.h"

marker_file_t *marker_open(const char *fn)
{
	gzFile fp = gzopen(fn, "r");
	if (fp == 0) return 0;
	marker_file_t *pf = (marker_file_t *) calloc(1, sizeof(marker_file_t));
	pf->fp = ks_init(fp);
	return pf;
}

int marker_close(marker_file_t *pf)
{
	if (pf == 0) return 0;
	free(pf->buf.s);
	kstream_t *ks = (kstream_t *) pf->fp;
	gzclose(ks->f);
	ks_destroy(ks);
	free(pf);
	return 0;
}

static int marker_parse(int l, char *s, marker_rec_t *pr)
{
	char *q, *r;
	int i, t;
	for (i = t = 0, q = s; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) pr->n = q;
		else if (t == 1) pr->p = strtol(q, &r, 10);
		else if (t == 2) pr->bin = strtol(q, &r, 10);
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 2) return -1;
	return 0;
}

int marker_read(marker_file_t *pf, marker_rec_t *r)
{
	int ret, dret;
file_read_more:
	ret = ks_getuntil((kstream_t*)pf->fp, KS_SEP_LINE, &pf->buf, &dret);
	if (ret < 0) return ret;
	ret = marker_parse(pf->buf.l, pf->buf.s, r);
	if (ret < 0) goto file_read_more;
	return ret;
}
