// cram_stubs.c — stub implementations for CRAM functions referenced by hts.c
// We don't support CRAM; these return error/NULL for any call path that
// attempts to use CRAM.
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include "htslib/hts.h"
#include "cram/cram.h"

// Forward declarations
struct cram_fd;
struct hFILE;

struct cram_fd *cram_dopen(struct hFILE *fp, const char *fn, const char *mode) {
    (void)fp; (void)fn; (void)mode;
    return NULL;  // CRAM not supported
}

int cram_close(struct cram_fd *fd) {
    (void)fd;
    return -1;
}

int cram_flush(struct cram_fd *fd) {
    (void)fd;
    return -1;
}

int cram_eof(struct cram_fd *fd) {
    (void)fd;
    return -1;
}

int cram_check_EOF(struct cram_fd *fd) {
    (void)fd;
    return 0;
}

int cram_set_option(struct cram_fd *fd, enum hts_fmt_option opt, ...) {
    (void)fd; (void)opt;
    return -1;
}

int cram_set_voption(struct cram_fd *fd, enum hts_fmt_option opt, va_list args) {
    (void)fd; (void)opt; (void)args;
    return -1;
}

struct hFILE *cram_hfile(struct cram_fd *fd) {
    (void)fd;
    return NULL;
}

int64_t cram_ptell(struct cram_fd *fd) {
    (void)fd;
    return -1;
}

void cram_index_free(struct cram_fd *fd) {
    (void)fd;
}

cram_index *cram_index_last(struct cram_fd *fd, int refid, cram_index *last) {
    (void)fd; (void)refid; (void)last;
    return NULL;
}

cram_index *cram_index_query(struct cram_fd *fd, int refid, int64_t pos, cram_index *last) {
    (void)fd; (void)refid; (void)pos; (void)last;
    return NULL;
}

cram_index *cram_index_query_last(struct cram_fd *fd, int refid, int64_t end) {
    (void)fd; (void)refid; (void)end;
    return NULL;
}
