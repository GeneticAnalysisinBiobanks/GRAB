// Redirect to the htscodecs copy
#include "../htscodecs/htscodecs/pooled_alloc.h"

// pool_free is commented out in the htscodecs header but needed by htslib
#ifndef CRAM_POOL_FREE_DEFINED
#define CRAM_POOL_FREE_DEFINED
static inline void pool_free(pool_alloc_t *p, void *ptr) {
    *(void **)ptr = p->free;
    p->free = ptr;
}
#endif
