// cram/cram.h — stub for htslib internal CRAM API
// We don't support CRAM; these stubs satisfy hts.c linkage.
#ifndef CRAM_CRAM_H
#define CRAM_CRAM_H

#include <stdarg.h>
#include <stdint.h>
#include "../htslib/cram.h"

#ifdef __cplusplus
extern "C" {
#endif

struct hFILE;

// ---- Types used by hts.c internally ----

// cram_index is a struct used as a linked-list node for CRAM index entries.
// hts.c accesses members: offset, slice, len, e_next.
typedef struct cram_index {
    int64_t  offset;
    int      slice;
    int      len;
    struct cram_index *e_next;
} cram_index;

// cram_range is a struct representing a query region.
// hts.c accesses members: refid, start (via initializer), end.
typedef struct cram_range {
    int       refid;
    int64_t   start;
    int64_t   end;
} cram_range;

// ---- Functions referenced by hts.c — implemented as stubs in cram_stubs.c ----

struct cram_fd *cram_dopen(struct hFILE *fp, const char *fn, const char *mode);
int cram_close(struct cram_fd *fd);
int cram_flush(struct cram_fd *fd);
int cram_eof(struct cram_fd *fd);
int cram_check_EOF(struct cram_fd *fd);
int cram_set_option(struct cram_fd *fd, enum hts_fmt_option opt, ...);
int cram_set_voption(struct cram_fd *fd, enum hts_fmt_option opt, va_list args);
struct hFILE *cram_hfile(struct cram_fd *fd);
int64_t cram_ptell(struct cram_fd *fd);

// Index functions
void cram_index_free(struct cram_fd *fd);
cram_index *cram_index_last(struct cram_fd *fd, int refid, cram_index *last);
cram_index *cram_index_query(struct cram_fd *fd, int refid, int64_t pos, cram_index *last);
cram_index *cram_index_query_last(struct cram_fd *fd, int refid, int64_t end);

#ifdef __cplusplus
}
#endif

#endif // CRAM_CRAM_H
