# LEAF

LEAF (Local Ethnicity-Aware Frequencies) stratifies the analysis cohort
into $K$ ancestry-homogeneous clusters and performs a cluster-specific
WtCoxG marker test in each cluster using an ancestry-matched synthetic
external allele frequency. The cluster-level results are pooled by a
fixed-effects score-based meta-analysis. Compared to a global WtCoxG run
with a single external reference, LEAF aims to (i) avoid spurious
batch-effect rejections caused by ancestry mismatch between the cohort and
the reference panel, and (ii) recover the precision gain of the external
information separately in each ancestral stratum.

Notation:

- $K$: number of ancestry clusters, either supplied via
  `--leaf-nclusters`, inferred from a precomputed
  `--leaf-cluster-file`, or defaulted to the number of reference
  populations.
- $P$: number of phenotypes processed in a single LEAF invocation.
- $J$: number of reference populations whose allele frequencies are
  loaded from plink2 `.afreq` files (one per population, supplied as a
  comma-separated list to `--ref-af`; $J \le 8$).
- $n$: total number of subjects in the union of all cluster slices
  (i.e. retained after `dropNaInColumns` and the `--keep`/`--remove`
  filters); $n_c$ denotes the size of cluster $c$.
- $\mathbf{R}^{(p)} \in \mathbb{R}^n$: residual vector from the global
  null-model fit for phenotype $p$; the cluster slice is
  $\mathbf{R}^{(p,c)} \in \mathbb{R}^{n_c}$.
- $\mathbf{w}^{(p)} \in \mathbb{R}^n$: per-subject weight vector for
  phenotype $p$; the cluster slice is $\mathbf{w}^{(p,c)}$.
- $\mathbf{f}_j \in \mathbb{R}^M$: per-marker reference allele
  frequency for population $j$, with allele orientation matched against
  the `.bim` REF/ALT pair; $M$ is the number of markers present in
  every reference file.
- $\mathrm{obs\_ct}_{m,j}$: total reference allele count for marker $m$
  in population $j$ (`OBS_CT` column of the plink2 `.afreq` file).
- $\hat{\boldsymbol{\alpha}}^{(p,c)} \in \mathbb{R}^J$: per-cluster,
  per-phenotype ancestry proportions estimated by Summix; they satisfy
  $\hat\alpha^{(p,c)}_j \ge 0$ and $\sum_j \hat\alpha^{(p,c)}_j = 1$.
- $\pi$: reference disease prevalence (`--prevalence`).

The null hypothesis is $H_0: \beta = 0$ in each cluster, with the
meta-analysis combining the per-cluster scores under a common $\beta$.
All LEAF-specific code lives in
[src/wtcoxg/leaf.cpp](../../src/wtcoxg/leaf.cpp) (declarations in
[src/wtcoxg/leaf.hpp](../../src/wtcoxg/leaf.hpp)); the per-cluster marker
test delegates to `WtCoxGMethod` from
[src/wtcoxg/wtcoxg.cpp](../../src/wtcoxg/wtcoxg.cpp), whose marker-level
formulas are documented in [wtcoxg.md](wtcoxg.md).

## 1. Cluster assignment

LEAF determines the cluster label of each retained subject in one of two
ways.

**K-means on principal components.** When `--pc-cols` supplies the
column names of $d$ principal components $\mathbf{X} \in
\mathbb{R}^{n \times d}$, LEAF runs $K$-means with $K$-means++
initialisation and `nstart = 25` random restarts
(`kmeansCluster` in
[src/wtcoxg/leaf.cpp](../../src/wtcoxg/leaf.cpp)). The objective is

$$
\min_{\mathbf{c}_1, \ldots, \mathbf{c}_K,\;
      \ell : \{1,\ldots,n\} \to \{1,\ldots,K\}}
\sum_{i=1}^n \| \mathbf{X}_i - \mathbf{c}_{\ell(i)} \|^2, \tag{1}
$$

where $\mathbf{c}_k$ is the cluster centroid. Lloyd iterations are
executed up to a maximum of 300 sweeps and terminate as soon as no label
changes during an assignment step; the assignment distances use the
factorisation
$\|\mathbf{x}_i - \mathbf{c}_j\|^2 = \|\mathbf{x}_i\|^2 - 2\mathbf{x}_i^\top
\mathbf{c}_j + \|\mathbf{c}_j\|^2$ so that the $n \times K$ distance
matrix is computed by a single Eigen / BLAS `dgemm`. The restart with the
lowest inertia is kept and the labels are renumbered to
$\{1, \ldots, K\}$ to match the R convention.

**Pre-computed cluster file.** When `--leaf-cluster-file` is supplied,
`parseLeafClusterFile` in
[src/wtcoxg/leaf.cpp](../../src/wtcoxg/leaf.cpp) reads a tab- or
whitespace-separated file with a mandatory header line, picks out the
`IID` (or `#IID`) and `cluster` columns by name, and returns the label
vector aligned to `usedIIDs`. The inferred $K = \max(\text{cluster})$ is
cross-checked against `--leaf-nclusters` when both are provided; missing
IIDs or non-integer cluster values are rejected. This path is the
recommended interface when an ancestry assignment is already available
from upstream tooling.

Either path yields the cluster index $c(i) \in \{1, \ldots, K\}$ for each
subject $i$, the per-cluster bitmask
$\mathbf{1}_{\mathcal C_c} \in \{0,1\}^n$, and the cluster-to-union index
map $\mathrm{idx}_c$ used by the marker engine. The genotype reader is
loaded once on the union of all $K$ cluster bitmasks, so each marker is
decoded a single time per chunk regardless of $K$.

## 2. Global null model and cluster slicing

For each phenotype $p$, LEAF fits a single global null model on all $n$
subjects (Section 2.1) and then slices the resulting residual and weight
vectors to the $K$ cluster memberships (Section 2.2). This contrasts with
a per-cluster refit, in which the per-cluster covariate effects would
absorb any residual ancestry signal that the global covariate adjustment
left behind. Following the LEAF reference R implementation, the global
fit is preferred because it produces residuals whose mean-zero structure
holds across all clusters jointly.

### 2.1 Weight construction

LEAF uses a fixed control-to-case sampling-weight ratio that depends only
on the reference prevalence $\pi$, regardless of the empirical
case/control split in the cohort:

$$
w_i^{(p)}
= \begin{cases}
1, & \text{if } \text{event}_i^{(p)} = 1, \\[2pt]
\dfrac{1 - \pi}{\pi}, & \text{if } \text{event}_i^{(p)} = 0.
\end{cases} \tag{2}
$$

This weight is computed by `leafRegrWeight` in
[src/wtcoxg/leaf.cpp](../../src/wtcoxg/leaf.cpp) and differs from
WtCoxG's `regression::calRegrWeight`, which rescales the control weight
by the empirical case/control ratio. The LEAF form matches the reference
R pipeline, in which controls always carry the fixed weight
$(1-\pi)/\pi$.

### 2.2 Residual fitting

The residual vector $\mathbf{R}^{(p)}$ depends on the phenotype type
inferred from `--pheno-name`:

- Binary phenotype (column $y \in \{0,1\}$): weighted logistic
  regression with design matrix $[\mathbf{1}_n \mid \mathbf{X}_{\text{cov}}]$
  and weights from (2), returning the response residuals $y_i - \hat\mu_i$
  via `regression::logisticResiduals`.
- Survival phenotype (`TIME:EVENT`): weighted Cox proportional-hazards
  regression with design matrix $\mathbf{X}_{\text{cov}}$ (no intercept
  column) and weights from (2), returning the martingale residuals
  $\delta_i - \hat\Lambda_0(t_i) e^{\hat\beta^\top \mathbf{X}_i}$
  via `regression::coxResiduals`.

Both fitters live in
[src/util/regression.cpp](../../src/util/regression.cpp); the formal
definitions are described in [residuals.md](residuals.md). The cluster
slice is

$$
\mathbf{R}^{(p,c)} = (R_i^{(p)})_{i \in \mathcal C_c},
\qquad
\mathbf{w}^{(p,c)} = (w_i^{(p)})_{i \in \mathcal C_c}, \tag{3}
$$

and the global weight sum $S_w^{(p)} = \sum_{i=1}^n w_i^{(p)}$ is
retained for use in (8).

## 3. Per-cluster ancestry-matched reference

LEAF replaces the user-supplied per-marker external allele frequency
$\hat f_{\mathrm{ext}}$ with a per-(phenotype, cluster) synthesis from $J$
reference-population frequencies, weighted by Summix-estimated ancestry
proportions.

### 3.1 Per-cluster, per-phenotype observed allele frequency

For each (phenotype $p$, cluster $c$, marker $m$), LEAF computes case and
control means $\hat f_1^{(p,c,m)}, \hat f_0^{(p,c,m)}$ over the
cluster's subject slice in a single genotype scan, then forms the
equal-blend observed allele frequency

$$
\hat f_{\mathrm{obs}}^{(p,c,m)}
= \tfrac{1}{2}\!\left(\hat f_0^{(p,c,m)} + \hat f_1^{(p,c,m)}\right) \tag{4}
$$

whenever both $n_0^{(p,c,m)} > 0$ and $n_1^{(p,c,m)} > 0$. Markers
without any case or any control in the cluster contribute NaN and are
filtered out by Summix. The folded value
$\hat f_{\mathrm{int}}^{(p,c,m)} = \min(\hat f_{\mathrm{obs}}^{(p,c,m)},
1 - \hat f_{\mathrm{obs}}^{(p,c,m)})$ is stored separately for use as the
MAF-group key in the WtCoxG Phase 2 estimator.

### 3.2 Summix ancestry estimation

For each (phenotype, cluster) pair, Summix solves

$$
\hat{\boldsymbol{\alpha}}^{(p,c)}
= \arg\min_{\boldsymbol{\alpha} \in \Delta^{J-1}}
  \sum_{m=1}^M
  \left(\hat f_{\mathrm{obs}}^{(p,c,m)}
        - \sum_{j=1}^J \alpha_j f_{m,j}\right)^2, \tag{5}
$$

where $\Delta^{J-1} = \{\boldsymbol{\alpha} : \alpha_j \ge 0,
\sum_j \alpha_j = 1\}$. The implementation enumerates all
$2^J - 1$ non-empty active sets of references, solves the equality-
constrained least-squares problem on each by the $|\,\mathcal S\,|+1$
dimensional KKT system

$$
\begin{pmatrix}
  2\,\mathbf{D}_{\mathcal S}^\top \mathbf{D}_{\mathcal S} & \mathbf{1} \\[2pt]
  \mathbf{1}^\top & 0
\end{pmatrix}
\begin{pmatrix}
  \boldsymbol{\alpha}_{\mathcal S} \\[2pt] \lambda
\end{pmatrix}
=
\begin{pmatrix}
  2\,\mathbf{D}_{\mathcal S}^\top \hat{\mathbf{f}}_{\mathrm{obs}} \\[2pt] 1
\end{pmatrix}, \tag{6}
$$

retains only solutions with $\boldsymbol{\alpha}_{\mathcal S} \ge 0$, and
keeps the lowest-objective solution. Following the Summix reference (R
package, [github.com/hendriau/Summix](https://github.com/hendriau/Summix)),
LEAF then drops any reference whose pass-1 proportion is below $0.01$ and
refits on the surviving subset; this is the `<1%` removal step.
Implemented in `summixEstimate` /
`summixCalcConstrainedLS` in
[src/wtcoxg/leaf.cpp](../../src/wtcoxg/leaf.cpp). The cap $J \le 8$ keeps
the enumeration at $\le 255$ subsets per fit.

### 3.3 Ancestry-matched allele frequency and effective sample size

The synthesised per-cluster, per-phenotype external allele frequency at
marker $m$ is the weighted average

$$
\hat f_{\mathrm{ext}}^{(p,c,m)}
= \sum_{j=1}^J \hat\alpha_j^{(p,c)}\, f_{m,j}, \tag{7}
$$

and the effective external sample size is the harmonic-style aggregator

$$
\mathrm{obs\_ct}^{(p,c,m)}
= \left(
    \sum_{j : \hat\alpha_j^{(p,c)} > 0,\; \mathrm{obs\_ct}_{m,j} > 0}
    \frac{(\hat\alpha_j^{(p,c)})^2}{\mathrm{obs\_ct}_{m,j}}
  \right)^{-1}, \tag{8}
$$

with the convention $\mathrm{obs\_ct}^{(p,c,m)} = 0$ when the
denominator is zero. Equation (8) is the variance-minimising effective
$n$ for the weighted sum
$\sum_j \hat\alpha_j^{(p,c)} \hat f_{m,j}$ under the independent-panel
binomial variance approximation
$\mathbb V(\hat f_{m,j}) \propto 1/\mathrm{obs\_ct}_{m,j}$. The
resulting pair
$(\hat f_{\mathrm{ext}}^{(p,c,m)}, \mathrm{obs\_ct}^{(p,c,m)})$ is the
single per-marker external block that the per-cluster `WtCoxGMethod`
consumes.

## 4. Per-cluster WtCoxG marker test

For each (phenotype $p$, cluster $c$), LEAF instantiates an independent
`WtCoxGMethod` from the cluster-sliced residuals
$\mathbf{R}^{(p,c)}$, weights $\mathbf{w}^{(p,c)}$, and the synthesised
external reference of Section 3. The instantiation runs the full WtCoxG
Phase 2 estimator on the cluster slice (per-MAF-group $\tau$, $\sigma^2$,
$b$, and variance ratios; see [wtcoxg.md](wtcoxg.md) Section 3) with one
calibration that is specific to LEAF.

### 4.1 Global weight normalisation for $\widetilde{\mathbf{w}}$

`testBatchEffects` accepts an optional `globalSumWeight` argument that
overrides the cluster-local denominator in the normalisation
$\widetilde w_i = w_i / (2 \sum_j w_j)$ of equation (7) in
[wtcoxg.md](wtcoxg.md). LEAF passes the global weight sum
$S_w^{(p)} = \sum_{i=1}^n w_i^{(p)}$ for the **outer** quantities
$\widetilde{\mathbf{w}}$ that feed the variance-ratio correction
$\rho_{w0}$ and the per-MAF-group $V_{\mathrm{Sbat}}$ used by the
$(\tau, \sigma^2)$ estimator, so that the GRM half-storage bilinear form
$H_{\boldsymbol{\Phi}}(\widetilde{\mathbf{w}}, \widetilde{\mathbf{w}})$
(defined in [wtcoxg.md](wtcoxg.md); off-diagonal entries of
$\boldsymbol{\Phi}$ counted once) remains commensurate with the
cluster-shared external-reference variance term $1 / \mathrm{obs\_ct}$.
The **inner** quantities used inside the $b$-estimator's
`fun.optimalWeight` retain the cluster-local normalisation
$\widetilde w_i^{\mathrm{loc}} = w_i / (2 \sum_{j \in \mathcal C_c} w_j)$,
mirroring the LEAF.R reference. The two scales coincide when there is
exactly one cluster, in which case LEAF reduces exactly to WtCoxG.

The departure from the symmetric quadratic-form derivation of
Li et al. (2025) described in the notation section of
[wtcoxg.md](wtcoxg.md) propagates to LEAF directly: the per-cluster
$\rho_{w0}$, $\rho_{\mathrm{int}}$, and $\rho_{\mathrm{ext}}$ are
evaluated with the half-storage bilinear form $H_{\boldsymbol{\Phi}}$
on the cluster slice, not with the symmetric quadratic form
$\mathbf{x}^\top \boldsymbol{\Phi}_c \mathbf{x}$. The magnitude of
the resulting numerical discrepancy is larger in LEAF than in
single-cluster WtCoxG because intra-cluster relatedness is, by
construction, denser than relatedness in the union cohort.

### 4.2 Per-cluster outputs

For each marker, the per-cluster `WtCoxGMethod::computeDual` returns the
quadruple $(p_{\mathrm{ext}}^{(c)}, p_{\mathrm{noext}}^{(c)},
S_{\mathrm{ext}}^{(c)}, S_{\mathrm{noext}}^{(c)})$ via the same
mechanism as the standalone WtCoxG run; the per-cluster batch-effect
p-value $p_{\mathrm{bat}}^{(c)}$ is read from the cluster's
`WtCoxGRefInfo` table. The cluster results are not written to disk
individually; only the meta-analysed p-values and the per-cluster
diagnostics enter the final output (Section 6).

## 5. Fixed-effects meta-analysis

For each marker, LEAF combines the per-cluster scores under a common
$\beta$. Let $S_c$ denote the per-cluster score and $p_c$ its WtCoxG
p-value (either the external-conditional p-value for the meta-extended
column or the internal-only p-value for the meta-internal column). From
the chi-square calibration

$$
\widehat{\mathbb V}(S_c)
= \frac{S_c^2}{\chi^2_{1,\, 1 - p_c}}, \tag{9}
$$

where $\chi^2_{1, 1-p}$ is the upper-$p$ quantile of a one-degree-of-
freedom chi-square distribution (`math::qchisq(p, 1, false, false)` in
the code), the fixed-effects meta z-score is

$$
Z_{\mathrm{meta}}
= \frac{\sum_{c : \mathcal V_c} S_c}{\sqrt{\sum_{c : \mathcal V_c} \widehat{\mathbb V}(S_c)}}, \tag{10}
$$

where the sum runs over clusters $c$ satisfying $S_c \ne \mathrm{NaN}$,
$p_c \in (0, 1)$, and $\widehat{\mathbb V}(S_c) \ne \mathrm{NaN}$ (the
condition is denoted $\mathcal V_c$ in (10)). Clusters in which the
external branch fails the QC gates of [wtcoxg.md](wtcoxg.md) Section 6
are silently excluded from the corresponding meta column. The two-sided
meta p-value is

$$
p_{\mathrm{meta}}
= \mathrm{erfc}\!\left(\frac{|Z_{\mathrm{meta}}|}{\sqrt{2}}\right). \tag{11}
$$

The chi-square calibration in (9) is the same device used internally by
WtCoxG to back out an SPA-calibrated score variance (equation (27) of
[wtcoxg.md](wtcoxg.md)). At the meta stage it serves a different
purpose: rather than calibrating a single test, it propagates the
saddlepoint-corrected scale of each cluster's score into the variance
denominator of the pooled z-score, so that the meta p-value inherits the
SPA correction component-wise. A pure normal-approximation variance would
be too small in the tails and would produce inflated $|Z_{\mathrm{meta}}|$
under heavy ascertainment.

If $\sum_c \widehat{\mathbb V}(S_c) \le 0$ (e.g. all clusters dropped due
to QC), the meta p-value is reported as `NaN`.

## 6. Output columns

For each marker, `LEAFMethod::getResultVec` appends $2 + 3K$ columns
after the standard meta block
(`CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC HWE_P`):

| Column           | Definition                                                                 |
|------------------|----------------------------------------------------------------------------|
| `meta.p_ext`     | Fixed-effects meta p-value (11) over the per-cluster external-conditional scores. |
| `meta.p_noext`   | Fixed-effects meta p-value (11) over the per-cluster internal-only scores.        |
| `cl<c>.p_ext`    | Cluster $c$ external-conditional p-value $p_{\mathrm{ext}}^{(c)}$ (cf. `p_ext` in [wtcoxg.md](wtcoxg.md)). |
| `cl<c>.p_noext`  | Cluster $c$ internal-only p-value $p_{\mathrm{noext}}^{(c)}$ (cf. `p_noext` in [wtcoxg.md](wtcoxg.md)).   |
| `cl<c>.p_batch`  | Cluster $c$ batch-effect p-value $p_{\mathrm{bat}}^{(c)}$ (cf. `p_batch` in [wtcoxg.md](wtcoxg.md)).      |

The cluster indices $c$ run from $1$ to $K$ in the order returned by the
K-means algorithm (or the order encoded in `--leaf-cluster-file`).
Per-cluster z-scores are not written; the meta-analysis uses the raw
per-cluster scores $S_c$ for combination, and the per-cluster `p_ext` and
`p_noext` columns provide the diagnostic detail.

## 7. Engine integration and parallelism

The dispatcher entry points are `runLEAFPheno` (single phenotype) and
`runLEAF` (multi-phenotype) in
[src/wtcoxg/leaf.cpp](../../src/wtcoxg/leaf.cpp). Five pipeline phases
run in order:

- **A. Shared loading.** Subject data, K-means or cluster file,
  per-cluster bitmasks, union bitmask, genotype reader, reference
  `.afreq` files, and per-cluster sparse GRMs are all loaded once on the
  main thread.
- **B. Parallel global null-model fits.** $P$ global fits are dispatched
  to $\min(\text{nthreads}, P)$ workers via an atomic task counter.
- **C. Per-(phenotype, cluster) genotype scan and Summix.** A single
  marker-by-marker pass populates the $P \times K \times M$ table of
  per-stratum allele frequencies; the subsequent $P \times K$ Summix fits
  are parallelised with $\min(\text{nthreads}, P \cdot K)$ workers.
- **D. Per-(phenotype, cluster) batch-effect testing.** The
  `testBatchEffects` of [wtcoxg.md](wtcoxg.md) Section 3 is invoked
  independently for each of the $P \cdot K$ slices, parallelised the
  same way.
- **E. Single multi-phenotype marker engine call.** The $P$ LEAF tasks
  are passed to the shared `multiPhenoEngine`
  ([src/engine/marker.cpp](../../src/engine/marker.cpp)), which decodes
  each genotype chunk once and dispatches it to the per-phenotype
  `LEAFMethod` clones.

The marker engine evaluates the per-cluster scores via the fused-GEMM
path: `LEAFMethod::fillUnionResiduals` populates $2K$ residual columns
per phenotype (one zero-padded residual column and one zero-padded
indicator column per cluster), so that the engine's union-level matmul
$\mathbf{A}^\top \mathbf{G}_{\mathrm{batch}}$ delivers, per marker and
per cluster, the pair $(\mathbf{R}^{(c)\top} \mathbf{g}^{(c)},
\sum_{i \in \mathcal C_c} g_i)$ needed by
`WtCoxGMethod::computeDualFromScalars`. The per-cluster gather of the
genotype vector is therefore avoided in the per-marker inner loop, and
the meta-analysis (11) is computed once per marker per phenotype on the
$K$ scalar pairs returned by the GEMM.

## 8. Required and optional command-line flags

The minimal invocation supplies a phenotype file with the response and
covariate columns named via `--pheno-name`, one `.afreq` file per
reference population via a comma-separated `--ref-af`, the disease
prevalence via `--prevalence`, and either the principal-component column
names via `--pc-cols` (for K-means clustering) or a precomputed cluster
file via `--leaf-cluster-file`. The optional flag `--leaf-nclusters`
overrides the default cluster count
$K = $ number of `.afreq` files when K-means is run; it is
cross-checked against the maximum cluster value in
`--leaf-cluster-file` when both are provided. The remaining flags
(`--batch-effect-p-threshold`, `--spa-z-threshold`,
`--outlier-iqr-threshold`, `--threads`, `--chunk-size`, `--geno`,
`--maf`, `--mac`, `--hwe`) behave identically to the WtCoxG
documentation; see [wtcoxg.md](wtcoxg.md) Section 6 and
[src/cli/flags.hpp](../../src/cli/flags.hpp).
