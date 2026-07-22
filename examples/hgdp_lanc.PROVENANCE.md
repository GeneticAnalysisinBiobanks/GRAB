# `examples/hgdp_lanc.*` — provenance

Tiny RFMix → `.lanc` fixture built per
`dev-notes/methods/recode_rfmix/02_tiny_example.md` (Plan B). It exercises the
full RFMix → `--make-lanc` pipeline and is the verification ground truth for
`dev-notes/methods/recode_rfmix/01_lanc_format_and_recoding.md` §8.

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

- Chromosomes: chr20, chr21, chr22.
- Per chromosome: **1000 markers** in the `.lanc`; **total M = 3000**.
- Marker class: biallelic SNVs, MAF ≥ 0.05 (MAF computed over all 4091 BCF
  samples via `bcftools +fill-tags -t MAF`), restricted to positions hosting
  exactly one biallelic SNV in the BCF (so every extracted site is unambiguous),
  then stride-thinned to be evenly spread across the selection region.

### Sub-region restriction (Plan B §3 sparse-input contingency)

1000 SNVs spread across a whole chromosome are too sparse for RFMix. Sites were
therefore drawn from a dense **20–40 Mb sub-region** on each chromosome (well
inside the genetic-map-covered range: chr20 81 kb–64.3 Mb, chr21 10.3–46.7 Mb,
chr22 15.3–50.8 Mb). Genetic span of the 20–40 Mb region: chr20 ≈ 10.9 cM,
chr21 ≈ 31.9 cM, chr22 ≈ 39.0 cM — sufficient for multiple RFMix windows per
chromosome.

### Sacrificial first site (exact-1000 accounting)

The converter's window cursor requires `window.spos <= pos0` where `pos0` is the
htslib 0-based POS, while RFMix writes each window's `spos` as the 1-based POS of
its first SNP. Consequently the very first SNP of each chromosome (the one at the
first window's `spos`) fails the test and is dropped — exactly one marker per
chromosome. To land exactly 1000 markers per chromosome, **1001 sites** were
supplied to RFMix per chromosome; the converter drops the first, leaving 1000.
The dropped (sacrificial) positions were chr20:20001903, chr21:20001213,
chr22:20000241.

## RFMix

- Version: **RFMix v2.03-r0**.
- Command (per chromosome), all long options in `--opt=value` form (RFMix's
  long-option parser calls `atoi(optarg)` and segfaults on the space-separated
  form):

  ```
  rfmix -f query.chrN.vcf.gz -r ref.chrN.vcf.gz -m ref.map -g genmap.chrN.tsv \
        -o rfmix.chrN --chromosome=chrN --n-threads=100 --random-seed=20240722
  ```

- No `-s` / `-c` / `-G` / `-n` overrides — RFMix defaults were used.
- Genetic map supplied to RFMix in `chrom  physical_position  cM` order with
  `chr`-prefixed contig names.
- Reference and query VCFs extracted at the **same** 1001 positions per
  chromosome (so RFMix's site intersection is exact) and fed unchanged to both
  RFMix and the converter (sample order preserved).

- K = 4 (`#Subpopulation order/codes: AFR=0  CSA=1  EAS=2  EUR=3`).
- RFMix windows (msp data rows) per chromosome:

  | Chromosome | Windows |
  |---|---|
  | chr20 | 29 |
  | chr21 | 66 |
  | chr22 | 74 |

## Converter

```
grab2 --make-lanc --rfmix-msp <WORK>/rfmix --vcf <WORK>/query \
      --out examples/hgdp_lanc --compression zst --compression-level 3 --threads 100
```

Deviation from the Plan B recipe: the recipe lists `--compression-level 3`
alone, but this binary's generic validation requires `--compression gz|zst`
whenever `--compression-level` is explicit. `--compression zst` is inert for
`--make-lanc` (the branch ignores it; the `.lanc` is always internally
zstd-framed and the level is the frame level), so `--compression zst
--compression-level 3` produces exactly the intended level-3 `.lanc`. Omitting
`--compression-level` would yield a byte-identical file (default frame level 3).

## Result and self-check

- `.lanc` dimensions (confirmed by the reader via `--cal-phi`):
  **4 ancestries × 3000 markers × 925 samples**, 3 chromosome files.
- Round-trip determinism: running `--make-lanc` twice yields byte-identical
  `.lanc` / `.bim` / `.fam`.
- Files shipped under `examples/`:
  `hgdp_lanc.chr{20,21,22}.lanc`, `hgdp_lanc.chr{20,21,22}.bim`,
  `hgdp_lanc.fam`, and the ground-truth `rfmix.chr{20,21,22}.msp.tsv`.
