// sam_stubs.c — stub implementations for SAM/FASTQ functions referenced by
// hts.c and vcf.c.  We only use htslib for VCF/BCF; SAM/CRAM/FASTQ code
// paths will never execute but the linker needs these symbols.
#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "sam_internal.h"

// --- sam_hdr_t lifecycle ---

void sam_hdr_destroy(sam_hdr_t *h) {
    if (!h) return;
    if (h->target_name) {
        int32_t i;
        for (i = 0; i < h->n_targets; i++)
            free(h->target_name[i]);
        free(h->target_name);
    }
    free(h->target_len);
    free(h->text);
    free(h);
}

// --- threading (no-ops) ---

int sam_set_thread_pool(htsFile *fp, htsThreadPool *p) {
    (void)fp; (void)p;
    return 0;
}

int sam_set_threads(htsFile *fp, int nthreads) {
    (void)fp; (void)nthreads;
    return 0;
}

// --- SAM state (no-ops) ---

int sam_state_destroy(samFile *fp) {
    (void)fp;
    return 0;
}

// --- FASTQ state (no-ops) ---

int fastq_state_set(samFile *fp, enum hts_fmt_option opt, ...) {
    (void)fp; (void)opt;
    return -1;
}

void fastq_state_destroy(samFile *fp) {
    (void)fp;
}

// --- index ---

int sam_idx_save(htsFile *fp) {
    (void)fp;
    return -1;
}

// --- bam1_t data reallocation ---

int sam_realloc_bam_data(bam1_t *b, size_t desired) {
    (void)b; (void)desired;
    return -1;
}
