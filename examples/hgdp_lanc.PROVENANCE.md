# `examples/hgdp_lanc.*` — provenance

Tiny RFMix → `.lanc` fixture built per
`dev-notes/methods/recode_rfmix/02_tiny_example.md` (Plan B) and regenerated for
the merged single-file `.lanc` format per
`dev-notes/methods/recode_rfmix/03_merged_and_fixes.md` (items #1 and #4). It
exercises the full RFMix → `--make-lanc` pipeline and is the verification ground
truth for `dev-notes/methods/recode_rfmix/01_lanc_format_and_recoding.md` §8.

## Files shipped under `examples/`

- **`hgdp_lanc.lanc`** — ONE merged plane-separated `.lanc` (all chromosomes as
  contiguous segments, in chromosome order; 32-byte header, per-segment footer
  directory). No per-chromosome `.lanc`/`.bim` files.
- **`hgdp_lanc.bim`** — one merged PLINK BIM, all 3000 markers in chromosome
  order (chr20, chr21, chr22).
- **`hgdp_lanc.fam`** — shared query-sample list (N = 925).
- **`rfmix.chr{20,21,22}.msp.tsv`** — the RFMix output, kept as verification
  ground truth.

## Source data

- Merged HGDP+1000 Genomes phased BCFs (SHAPEIT5), b38, per chromosome:
  `/home/share/master_data1/1000Genomes/hgdp_1kg_phased_haplotypes_v2/hgdp1kgp_chr{20,21,22}.filtered.SNV_INDEL.phased.shapeit5.bcf`
  (4091 samples).
- Sample metadata: `gnomad_meta_updated.tsv` in the same directory
  (col 1 `s`, col 174 `hgdp_tgp_meta.Project`, col 177 `hgdp_tgp_meta.Genetic.region`).
- Genetic maps (b38): `/home/share/master_data1/GeneticMaps/b38/chr{20,21,22}.b38.gmap.gz`.

## Reference panel (RFMix-internal; NOT stored in `.lanc`)

All 1000-Genomes-project samples with `Genetic.region ∈ {AFR, EUR, EAS, CSA}`
(AMR excluded ⇒ K = 4), intersected with the BCF sample list:

| Ancestry | Code | Samples |
|---|---|---|
| AFR | 0 | 879 |
| CSA | 1 | 599 |
| EAS | 2 | 583 |
| EUR | 3 | 618 |
| **Total** | | **2679** |

(1000 Genomes "SAS" is labeled `CSA` in this metadata.)

## Query cohort (stored in `.lanc`)

All HGDP samples intersected with the BCF sample list: **N = 925**, all regions
(AFR 107, AMR 62, CSA 183, EAS 233, EUR 153, MID 157, OCE 30). The MID/OCE/AMR
HGDP individuals are force-assigned by RFMix to the nearest of the 4 reference
ancestries.

## Markers

- Chromosomes: chr20, chr21, chr22 → **3 segments** in the merged `.lanc`.
- Per chromosome: **exactly 1000 markers**; **total M = 3000**.
- Marker class: biallelic SNVs, MAF ≥ 0.05 (MAF computed over all 4091 BCF
  samples via `bcftools +fill-tags -t MAF`), restricted to positions hosting
  exactly one biallelic SNV in the BCF (so every extracted site is unambiguous),
  then stride-thinned to be evenly spread across the selection region.

### Sub-region restriction (Plan B §3 sparse-input contingency)

1000 SNVs spread across a whole chromosome are too sparse for RFMix. Sites were
therefore drawn from a dense **20–40 Mb sub-region** on each chromosome (well
inside the genetic-map-covered range: chr20 81 kb–64.3 Mb, chr21 10.3–46.7 Mb,
chr22 15.3–50.8 Mb). The 1000 selected sites span, per chromosome:

| Chromosome | First POS | Last POS | Windows |
|---|---|---|---|
| chr20 | 20,001,903 | 39,620,640 | 33 |
| chr21 | 20,001,213 | 39,863,737 | 62 |
| chr22 | 20,000,241 | 39,524,652 | 77 |

### Window coordinate convention (item #4 fixed — NO sacrificial first site)

RFMix writes each window's `spos`/`epos` as **1-based** VCF positions. Internal
windows tile half-open (`epos_i == spos_{i+1}`); the **last window of each
chromosome reports `epos` = the last SNP's position (INCLUSIVE)**. The converter
now stores windows 0-based (subtracting 1) and treats each chromosome's last
window as inclusive, so the membership test operates in htslib's 0-based
`rec->pos` space for every site. Consequently **all 1000 selected sites per
chromosome map into a window** — the first-window first SNP (`POS == spos0`) is
kept, window-boundary SNPs land in the correct later window, and each
chromosome's final SNP (`POS == last-window epos`) is kept. There is **no
sacrificial first site**; exactly 1000 sites are supplied and 1000 are stored.

## RFMix

- Version: **RFMix v2.03-r0**.
- Command (per chromosome), all long options in `--opt=value` form (RFMix's
  long-option parser calls `atoi(optarg)` and mishandles the space-separated
  form):

  ```
  rfmix -f query.chrN.vcf.gz -r ref.chrN.vcf.gz -m ref.map -g genmap.chrN.tsv \
        -o rfmix.chrN --chromosome=chrN --n-threads=32 --random-seed=20240722
  ```

- No `-s` / `-c` / `-G` / `-e` overrides — RFMix defaults were used.
- Genetic map supplied to RFMix in `chrom  physical_position  cM` order with
  `chr`-prefixed contig names.
- Reference and query VCFs extracted at the **same** 1000 positions per
  chromosome (so RFMix's site intersection is exact) and fed unchanged to both
  RFMix and the converter (sample order preserved).
- K = 4 (`#Subpopulation order/codes: AFR=0  CSA=1  EAS=2  EUR=3`).

## Converter

```
grab2 --make-lanc --rfmix-msp <WORK>/rfmix --vcf <WORK>/query \
      --out examples/hgdp_lanc --compression-level 3 --threads 100
```

`--make-lanc` uses `--compression-level` for its internal zstd frames (default 3)
and ignores `--compression`, so no `--compression` flag is required (CLI item #5).
The `.lanc` output is byte-reproducible for a fixed GRAB binary.

## Verification (Plan A §8 / item #4 anchors — all pass)

- **(a)** `LancData("examples/hgdp_lanc")` opens: N = 925, K = 4, nSeg = 3, one
  `.lanc` file, total 3000 markers.
- **(b)** Per-window assigned SNV count **== RFMix `n snps`** for every window,
  per-window and summed: chr20 1000/1000, chr21 1000/1000, chr22 1000/1000
  (total 3000/3000). This is the independent anchor of the item-#4 fix (does not
  depend on the converter's own convention).
- **(c)** Decode invariant Σ_k Σ_i hapcount = 2N (= 1850) for every marker; max
  deviation 0; 0 unassigned (0xF) cells.
- **(d)** Determinism: `--make-lanc` run twice yields byte-identical
  `.lanc` / `.bim` / `.fam`.
- **(e)** Ground-truth reconstruction from `query.chr*.vcf.gz` + `rfmix.chr*.msp.tsv`
  (0-based `[spos, epos)` with last-window-inclusive semantics) vs the `.lanc`
  decode via `getAllAncestries`: **max absolute difference 0** over all
  (marker, sample, ancestry) cells.
