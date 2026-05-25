// cram/string_alloc.h — minimal string pool allocator for htslib internal use
// Based on the original htslib cram/string_alloc API.
#ifndef CRAM_STRING_ALLOC_H
#define CRAM_STRING_ALLOC_H

#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct string_alloc_s {
    size_t  block_size;
    char   *cur_block;
    size_t  cur_used;
    size_t  cur_alloc;
    // Simple linked list of allocated blocks for cleanup
    struct sa_block {
        struct sa_block *next;
        char data[];
    } *blocks;
} string_alloc_t;

static inline string_alloc_t *string_pool_create(size_t block_size) {
    string_alloc_t *a = (string_alloc_t *)calloc(1, sizeof(*a));
    if (!a) return NULL;
    a->block_size = block_size ? block_size : 65536;
    return a;
}

static inline void string_pool_destroy(string_alloc_t *a) {
    if (!a) return;
    struct sa_block *b = a->blocks;
    while (b) {
        struct sa_block *n = b->next;
        free(b);
        b = n;
    }
    free(a);
}

static inline char *string_alloc(string_alloc_t *a, size_t len) {
    if (!a) return NULL;
    if (a->cur_used + len > a->cur_alloc) {
        size_t bs = a->block_size;
        if (len > bs) bs = len;
        struct sa_block *nb = (struct sa_block *)malloc(sizeof(struct sa_block) + bs);
        if (!nb) return NULL;
        nb->next = a->blocks;
        a->blocks = nb;
        a->cur_block = nb->data;
        a->cur_alloc = bs;
        a->cur_used  = 0;
    }
    char *p = a->cur_block + a->cur_used;
    a->cur_used += len;
    return p;
}

static inline char *string_dup(string_alloc_t *a, const char *s) {
    size_t len = strlen(s) + 1;
    char *p = string_alloc(a, len);
    if (p) memcpy(p, s, len);
    return p;
}

static inline char *string_ndup(string_alloc_t *a, const char *s, size_t len) {
    char *p = string_alloc(a, len + 1);
    if (p) {
        memcpy(p, s, len);
        p[len] = '\0';
    }
    return p;
}

#ifdef __cplusplus
}
#endif

#endif // CRAM_STRING_ALLOC_H
