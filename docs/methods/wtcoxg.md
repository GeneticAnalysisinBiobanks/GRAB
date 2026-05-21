# WtCoxG

WtCoxG (Weighted Cox-type G-test, Li et al. 2025) performs marker-level
score tests in case-control or survival GWAS samples whose sampling design
departs from the source population prevalence. For each marker it
optionally borrows an external allele-frequency estimate from a reference
panel to recover the precision lost to the unbalanced ascertainment, while
guarding against allele-frequency drift between the internal cohort and the
reference panel via a separate batch-effect screen.

Notation:

- $n$: number of retained subjects in the internal cohort.
- $\mathbf{G} = (G_1, \ldots, G_n)^\top$: per-subject alternate-allele
  count for the marker under test; missing entries are mean-imputed to
  $2\hat f_{\mathrm{int}}$ (see below). WtCoxG does not flip alleles to
  enforce $\text{MAF} \le 0.5$.
- $\mathbf{R} = (R_1, \ldots, R_n)^\top$: residuals from the upstream null
  fit. The fit is a weighted logistic regression (binary phenotype) or a
  weighted Cox proportional-hazards regression (survival phenotype), with
  per-subject case-control weights
  $w_i = 1$ for cases and $w_i = (n_{\mathrm{case}}/n_{\mathrm{ctrl}}) /
  [\pi/(1-\pi)]$ for controls, where $\pi$ is the user-supplied reference
  prevalence (see [residuals.md](residuals.md) for the residual
  definitions; the weight is implemented in `regression::calRegrWeight` in
  [src/wtcoxg/regression.cpp](../../src/wtcoxg/regression.cpp)).
- $\mathbf{w} = (w_1, \ldots, w_n)^\top$: the same sampling weights, also
  reused by the batch-effect screen.
- $\boldsymbol{\Phi}$: optional sparse GRM, $n \times n$, symmetric.
  When supplied, WtCoxG computes three variance-ratio corrections from
  its file-stored entries. Let $\mathcal{S}(\boldsymbol{\Phi})$ denote
  the set of index pairs $(i, j)$ at which $\boldsymbol{\Phi}$ is
  stored: by convention each diagonal entry $(i, i)$ is stored once and
  each off-diagonal pair $\{(i, j), (j, i)\}$ is stored only on one
  side (the lower-triangular file representation). Define the
  half-storage bilinear form
  $$
  H_{\boldsymbol{\Phi}}(\mathbf{x}, \mathbf{y})
  := \sum_{(i, j) \in \mathcal{S}(\boldsymbol{\Phi})} \Phi_{ij}\, x_i\, y_j
   = \sum_{i} \Phi_{ii} x_i y_i
   + \sum_{\substack{i > j \\ (i, j) \in \mathcal{S}(\boldsymbol{\Phi})}}
     \Phi_{ij} x_i y_j,
  $$
  which is what `SparseGRM::halfStorageSum`
  ([src/io/sparse_grm.hpp](../../src/io/sparse_grm.hpp)) returns. The
  off-diagonal entries are counted **once**, not twice, so
  $H_{\boldsymbol{\Phi}}(\mathbf{x}, \mathbf{x})$ is **not** the
  symmetric quadratic form $\mathbf{x}^\top \boldsymbol{\Phi}
  \mathbf{x}$; the two differ by a factor of two on every off-diagonal
  contribution. WtCoxG and LEAF use $H_{\boldsymbol{\Phi}}$ for the
  variance-ratio numerators in order to reproduce the LEAF.R reference
  exactly. *Departure from the published derivation.* The methods
  section of Li et al. (2025, Nat. Comput. Sci. 5: 1064-1079) writes
  every variance-ratio numerator as the symmetric quadratic form
  $\mathbf{x}^\top \boldsymbol{\Phi} \mathbf{x}$ in its equations
  (9), (10), (21), and (31), in which each off-diagonal kinship
  entry contributes twice by symmetry. The implementation in this
  repository instead computes $H_{\boldsymbol{\Phi}}(\mathbf{x},
  \mathbf{x})$, in which each off-diagonal entry contributes once.
  The two conventions coincide on diagonal contributions and on
  isolated kinship pairs but differ by a factor approaching two on
  every off-diagonal whenever intra-cluster relatedness is dense; a
  24 % relative discrepancy in $\rho_{\mathrm{int}}$ between the two
  conventions has been observed in one LEAF cluster
  (see [src/wtcoxg/wtcoxg.cpp](../../src/wtcoxg/wtcoxg.cpp), inline
  comment in `testBatchEffects`). The half-storage convention is
  retained here so that the WtCoxG and LEAF outputs reproduce the
  LEAF.R reference bit-for-bit; reproducing the closed-form
  expressions of Li et al. (2025) would require replacing every
  occurrence of $H_{\boldsymbol{\Phi}}$ in Sections 2 and 3.3 below
  with the symmetric quadratic form
  $\mathbf{x}^\top \boldsymbol{\Phi} \mathbf{x}$.
- $\hat f_0, \hat f_1$: control and case mean genotypes divided by two,
  computed from the internal cohort on non-missing entries only.
- $n_0, n_1$: non-missing control and case counts at the marker.
- $\hat f_{\mathrm{int}} = (\sum_{i: G_i \ne \mathrm{NA}} G_i) / (2(n_0 +
  n_1))$: internal allele-frequency estimate after the per-stratum
  averaging.
- $\hat f_{\mathrm{ext}}$ and $\mathrm{obs\_ct}$: external allele
  frequency and external total allele count parsed from the plink2
  `.afreq` file (`ALT_FREQS` and `OBS_CT` columns; allele orientation is
  resolved against the `.bim` so that $\hat f_{\mathrm{ext}}$ refers to
  the same allele as $\hat f_{\mathrm{int}}$). The external sample size
  proxy is $n_{\mathrm{ext}} = \mathrm{obs\_ct}/2$.
- $b \in [0, 1]$: the external-information weight, estimated per
  MAF-group (Section 3.3 below).
- $\hat f = (1-b)\hat f_{\mathrm{int}} + b \hat f_{\mathrm{ext}}$:
  combined allele-frequency estimate used to centre the genotype.
- $p_{\mathrm{cut}}$: batch-effect p-value threshold (default 0.05).

The null hypothesis is $H_0: \beta = 0$, where $\beta$ is the per-allele
log-odds-ratio (logistic) or log-hazard-ratio (Cox) of the marker on the
phenotype. All marker-level computations are implemented in
[src/wtcoxg/wtcoxg.cpp](../../src/wtcoxg/wtcoxg.cpp); the corresponding
header declarations are in
[src/wtcoxg/wtcoxg.hpp](../../src/wtcoxg/wtcoxg.hpp).

## 1. Marker matching and per-marker stratum statistics

For each marker present in the genotype file, WtCoxG looks up an entry in
the reference allele-frequency file by `(CHROM, ID)` and accepts the match
only when the `(REF, ALT)` alleles either coincide with `(bim.REF,
bim.ALT)` or are the reversed pair. In the latter case
$\hat f_{\mathrm{ext}}$ is set to $1 - \mathtt{ALT\_FREQS}$ so that the
external estimate refers to the `bim`-side `ALT` allele. A separate
two-column numeric fallback path is provided for the case where the
reference file has no header and rows are assumed to be in `.bim` order;
that path stores `ALT_FREQS` and `OBS_CT` directly.

A single sequential pass over the matched markers populates the stratum
statistics $(\hat f_0, \hat f_1, n_0, n_1, \hat f_{\mathrm{int}})$ from
the genotype matrix and the case-control indicator vector (see
`computeMarkerStats` in [src/wtcoxg/wtcoxg.cpp](../../src/wtcoxg/wtcoxg.cpp)).

## 2. Batch-effect screen

The batch-effect screen tests whether the internally weighted allele
frequency departs significantly from the external estimate. Define the
case weight $w_1 = 1$ and the control weight

$$
w_0 = \frac{1 - \pi}{\pi} \cdot \frac{e_r}{1 - e_r},
\qquad
e_r = \frac{n_1}{n_0 + n_1}. \tag{1}
$$

The internally weighted allele frequency is

$$
\hat f_{\mathrm{w}}
= \frac{\hat f_0\, w_0\, n_0 + \hat f_1\, w_1\, n_1}
       {w_0\, n_0 + w_1\, n_1}, \tag{2}
$$

and the pooled estimate that enters the batch-effect variance is

$$
\hat f_{\mathrm{pool}}
= \frac{\hat f_0\, w_0\, n_0 + \hat f_1\, w_1\, n_1 + \hat f_{\mathrm{ext}}\,
        n_{\mathrm{ext}}\, w_0}
       {w_0\, n_0 + w_1\, n_1 + w_0\, n_{\mathrm{ext}}}. \tag{3}
$$

The unadjusted batch-effect z-statistic is

$$
Z_{\mathrm{bat}}
= \frac{\hat f_{\mathrm{w}} - \hat f_{\mathrm{ext}}}{\sqrt{V_{\mathrm{bat}}}}
\tag{4}
$$

with

$$
V_{\mathrm{bat}}
= \left\{
    \frac{n_1 w_1^2 + n_0 w_0^2}{2(n_1 w_1 + n_0 w_0)^2}
    + \frac{1}{\mathrm{obs\_ct}}
  \right\}
  \hat f_{\mathrm{pool}}(1 - \hat f_{\mathrm{pool}}). \tag{5}
$$

When a sparse GRM is supplied, WtCoxG applies a per-marker variance-ratio
correction $\rho_{w0}$ to $Z_{\mathrm{bat}}$ before computing the p-value:

$$
Z_{\mathrm{bat}}^{\mathrm{adj}}
= \frac{Z_{\mathrm{bat}}}{\sqrt{\rho_{w0}}},
\qquad
p_{\mathrm{bat}} = 2\,\Phi\!\left(-\left|Z_{\mathrm{bat}}^{\mathrm{adj}}\right|\right). \tag{6}
$$

The correction is

$$
\rho_{w0}
= \frac{H_{\boldsymbol{\Phi}}(\widetilde{\mathbf{w}}, \widetilde{\mathbf{w}})
        + 1/(2\,\mathrm{obs\_ct})}
       {\sum_{i=1}^n \widetilde w_i^2
        + 1/(2\,\mathrm{obs\_ct})},
\qquad
\widetilde w_i = \frac{w_i}{2\sum_{j=1}^n w_j}. \tag{7}
$$

The numerator uses the half-storage bilinear form
$H_{\boldsymbol{\Phi}}$ defined in the notation above (each off-diagonal
GRM entry contributes once, not twice). The denominator is the plain
element-wise sum of squares $\sum_i \widetilde w_i^2$, which does not
depend on $\boldsymbol{\Phi}$. Both quantities are evaluated by
`SparseGRM::halfStorageSum` and `w1.array().square().sum()` respectively
in [src/wtcoxg/wtcoxg.cpp](../../src/wtcoxg/wtcoxg.cpp), matching the
LEAF.R reference. If no GRM is supplied, $\rho_{w0} = 1$.

A marker passes the batch-effect screen when
$p_{\mathrm{bat}} \ge p_{\mathrm{cut}}$, and the external reference is
incorporated only for those markers. Markers with
$p_{\mathrm{bat}} < p_{\mathrm{cut}}$ are routed through the
batch-detected branch described in Section 4.1, which integrates the
conditional p-value over the complementary tail region of the screen.

## 3. Per-MAF-group parameter estimation

Three parameters of the conditional p-value are not fixed a priori but
estimated from the empirical distribution of batch-effect p-values within
each internal-MAF bin: the mixture weight $\tau$ of the batch-affected
component, the per-MAF residual variance $\sigma^2$ of the external
estimate, and the external-information weight $b$. WtCoxG groups markers
by $\hat f_{\mathrm{int}}$ into the bins
$(-10^{-5}, 0.05], (0.05, 0.10], \ldots, (0.35, 0.40], (0.40,
\max \hat f_{\mathrm{int}}]$ and runs the three estimators independently
in each bin.

### 3.1 Estimation of $\tau$ and $\sigma^2$

Let $V_{\mathrm{Sbat}}$ denote the model-based variance of the batch
statistic at the bin centre $\mu = (\mathrm{lo}+\mathrm{hi})/2$:

$$
V_{\mathrm{Sbat}}
= \rho_{w0} \left\{
    \left(\sum_i \widetilde w_i^2\right) \cdot 2\mu(1-\mu)
    + \frac{\mu(1-\mu)}{\mathrm{obs\_ct}}
  \right\}, \tag{8}
$$

with $\rho_{w0} = 1$ when no GRM is supplied. For each cutoff
$c_j \in \{0.01, 0.11, 0.21, 0.31\}$, the empirical pass rate over a
widened MAF window $(\mathrm{lo} - 0.1, \mathrm{hi} + 0.1)$ is

$$
\hat p_{\mathrm{den}, j}
= \frac{1}{m} \sum_{k=1}^m \mathbb{I}\!\left(p_{\mathrm{bat}, k} > c_j\right). \tag{9}
$$

The corresponding model probability is

$$
p_j(\tau, \sigma^2)
= \tau \left[\Phi\!\left(\frac{u_j}{\sqrt{V_{\mathrm{Sbat}} + \sigma^2}}\right)
            - \Phi\!\left(\frac{\ell_j}{\sqrt{V_{\mathrm{Sbat}} + \sigma^2}}\right)
       \right]
+ (1-\tau)(1-c_j), \tag{10}
$$

where the integration limits are

$$
\ell_j = -z_{1-c_j/2}\sqrt{V_{\mathrm{Sbat}}},
\qquad
u_j = z_{1-c_j/2}\sqrt{V_{\mathrm{Sbat}}}. \tag{11}
$$

The implementation estimates $(\tau, \sigma^2)$ by Nelder-Mead
minimisation of

$$
Q(\tau, \sigma^2)
= \sum_j
  \left(
    \frac{\hat p_{\mathrm{den}, j} - p_j(\tau, \sigma^2)}{\hat p_{\mathrm{den}, j}}
  \right)^2,
\qquad (\tau, \sigma^2) \in [0,1]\times[0,1], \tag{12}
$$

with soft penalty walls at the box boundaries and a post-hoc clamp to
$[0,1]^2$. Markers in the bin all share the resulting $(\tau, \sigma^2)$.

### 3.2 Estimation of $b$

The external-information weight $b$ is selected to maximise the local
detection sensitivity at the genome-wide threshold $5 \times 10^{-8}$. For
each candidate $b \in [0, 1]$ and each trial case allele frequency $\mu_1$,
WtCoxG evaluates the conditional p-value

$$
p_{\mathrm{con}}(b, \mu_1)
= \frac{2\left\{\tau\, p_1(b, \mu_1) + (1-\tau)\, p_0(b, \mu_1)\right\}}
       {p_{\mathrm{deno}}(b, \mu_1)}, \tag{13}
$$

where $p_0(b, \mu_1)$ and $p_1(b, \mu_1)$ are the bivariate normal
rectangle probabilities of the pair $(S, S_{\mathrm{bat}})$ defined in
Section 4 evaluated at the trial $\mu_1$ (with $\sigma^2 = 0$ and
$\sigma^2 > 0$ respectively), the control allele frequency is fixed at the
bin centre, and $p_{\mathrm{deno}}(b, \mu_1)$ is the selection probability
of the screen. For each fixed $b$, the implementation finds
$\mu_1^\star(b)$ as the smallest $\mu_1$ that solves
$p_{\mathrm{con}}(b, \mu_1) = 5 \times 10^{-8}$ by Brent root-finding, and
then chooses

$$
\hat b = \arg\min_{0 \le b \le 1} \mu_1^\star(b) \tag{14}
$$

by one-dimensional Brent minimisation. The interpretation is that $\hat b$
minimises the case allele frequency required to declare genome-wide
significance in the bin.

### 3.3 Variance ratios

When a sparse GRM is supplied, three variance ratios are formed once per
MAF bin (or once per marker for $\rho_{w0}$). With
$\widetilde R_i^{(\mathrm{int})} = R_i - \bar R$ and
$\widetilde R_i^{(\mathrm{ext})}(b) = R_i - b\,\bar R$, the half-storage
bilinear forms

$$
Q_w = H_{\boldsymbol{\Phi}}(\widetilde{\mathbf{w}}, \widetilde{\mathbf{w}}),
\quad
Q_R = H_{\boldsymbol{\Phi}}(\widetilde{\mathbf{R}}^{(\mathrm{int})},
                            \widetilde{\mathbf{R}}^{(\mathrm{int})}),
\quad
Q_R^{\mathrm{ext}}(b)
= H_{\boldsymbol{\Phi}}(\widetilde{\mathbf{R}}^{(\mathrm{ext})}(b),
                         \widetilde{\mathbf{R}}^{(\mathrm{ext})}(b)) \tag{15}
$$

yield

$$
\rho_{\mathrm{int}}
= \frac{Q_R}{\sum_i (\widetilde R_i^{(\mathrm{int})})^2},
\qquad
\rho_{\mathrm{ext}}
= \frac{Q_R^{\mathrm{ext}}(\hat b) + 2\hat b^2 (\sum_i R_i)^2 / \mathrm{obs\_ct}}
       {\sum_i (\widetilde R_i^{(\mathrm{ext})}(\hat b))^2
        + 2\hat b^2 (\sum_i R_i)^2 / \mathrm{obs\_ct}}, \tag{16}
$$

with $\rho_{w0}$ as in (7). Each numerator uses the half-storage
bilinear form $H_{\boldsymbol{\Phi}}$ (off-diagonal entries counted
once); the denominators are plain element-wise sums of squares of the
centred residuals and do not involve $\boldsymbol{\Phi}$. When no GRM
is supplied, all three ratios are set to $1$.

## 4. Score statistic and conditional p-value

For markers that pass the QC thresholds (Section 6), WtCoxG forms the
combined allele frequency $\hat f = (1-b)\hat f_{\mathrm{int}} + b\hat
f_{\mathrm{ext}}$ and the score statistic

$$
S = \sum_{i=1}^n R_i (G_i - 2\hat f)
  = \mathbf{R}^\top \mathbf{G} - 2\hat f \sum_i R_i. \tag{17}
$$

The internal-only branch is the special case $b = 0$, in which case
$\hat f = \hat f_{\mathrm{int}}$ and $S$ reduces to the standard
centred-genotype score.

The conditional p-value is derived from the joint distribution of $S$ and
the batch-effect score statistic

$$
S_{\mathrm{bat}} = \sum_{i=1}^n \widetilde w_i G_i, \tag{18}
$$

which, under the diploid binomial model
$G_i \sim \mathrm{Binomial}(2, \hat f)$ independent across $i$ (extended by
the variance-ratio correction when relatedness is modelled), has variance

$$
\widehat{\mathbb{V}}(S_{\mathrm{bat}})
= \left(\sum_i \widetilde w_i^2\right) \cdot 2\hat f(1-\hat f)
  + \widehat{\mathbb{V}}(\hat f_{\mathrm{ext}}),
\qquad
\widehat{\mathbb{V}}(\hat f_{\mathrm{ext}})
= \frac{\hat f(1-\hat f)}{\mathrm{obs\_ct}} + \sigma^2. \tag{19}
$$

The variance contributed by the external estimate is set to zero in the
internal-only branch and in the $\sigma^2 = 0$ component of the conditional
mixture.

### 4.1 Branch A: batch effect detected ($p_{\mathrm{bat}} < p_{\mathrm{cut}}$)

When the batch-effect screen rejects, WtCoxG drops the external
information ($b$ is treated as $0$, so $\hat f = \hat f_{\mathrm{int}}$ and
the variance contribution $\widehat{\mathbb{V}}(\hat f_{\mathrm{ext}})$ is
absorbed only as a finite-reference panel term in the SPA CGF; see Section
5). It then evaluates the conditional p-value by integrating over the
complementary region $\{|S_{\mathrm{bat}}| > u\}$:

$$
p_{\mathrm{con}}^{\mathrm{A}}
= \frac{2\left\{\tau\, p_1^{\mathrm{A}} + (1-\tau)\, p_0^{\mathrm{A}}\right\}}
       {p_{\mathrm{deno}}^{\mathrm{A}}}, \tag{20}
$$

with the bivariate normal rectangle probabilities

$$
p_\bullet^{\mathrm{A}}
= \Pr\!\left(S/\sqrt{\rho_\bullet} \le -|S|/\sqrt{\rho_\bullet},\;
             S_{\mathrm{bat}}/\sqrt{\rho_{w\bullet}} \le \ell/\sqrt{\rho_{w0}}\right)
+ \Pr\!\left(S/\sqrt{\rho_\bullet} \le -|S|/\sqrt{\rho_\bullet},\;
             S_{\mathrm{bat}}/\sqrt{\rho_{w\bullet}} \ge u/\sqrt{\rho_{w0}}\right), \tag{21}
$$

evaluated under
$\boldsymbol{\Sigma}_\bullet =
\begin{pmatrix}\widehat{\mathbb{V}}_\bullet(S) & \widehat{\mathrm{Cov}}(S_{\mathrm{bat}}, S) \\
\widehat{\mathrm{Cov}}(S_{\mathrm{bat}}, S) & \widehat{\mathbb{V}}(S_{\mathrm{bat}}) + \sigma^2 \mathbb{I}\{\bullet = 1\}\end{pmatrix}$
for $\bullet \in \{0, 1\}$, with $\rho_\bullet = \rho_{\mathrm{int}}$ for
the score axis and $\rho_{w\bullet} = \rho_{w0}$ for the batch axis in
both components. The selection probability is

$$
p_{\mathrm{deno}}^{\mathrm{A}}
= 2\tau\, \Phi\!\left(\frac{\ell}{\sqrt{\rho_{w0}}\,\sqrt{\widehat{\mathbb{V}}(S_{\mathrm{bat}}) + \sigma^2}}\right)
+ (1-\tau)\, p_{\mathrm{cut}}. \tag{22}
$$

### 4.2 Branch B: batch effect not detected ($p_{\mathrm{bat}} \ge p_{\mathrm{cut}}$)

When the screen does not reject, WtCoxG retains the external information
and integrates over the interior region $\{\ell \le S_{\mathrm{bat}} \le u\}$:

$$
p_{\mathrm{con}}^{\mathrm{B}}
= \frac{2\left\{\tau\, p_1^{\mathrm{B}} + (1-\tau)\, p_0^{\mathrm{B}}\right\}}
       {p_{\mathrm{deno}}^{\mathrm{B}}}, \tag{23}
$$

with

$$
p_\bullet^{\mathrm{B}}
= \Pr\!\left(S/\sqrt{\rho_\bullet} \le -|S|/\sqrt{\rho_\bullet},\;
             \frac{\ell}{\sqrt{\rho_{w\bullet}}} \le S_{\mathrm{bat}} \le \frac{u}{\sqrt{\rho_{w\bullet}}}\right), \tag{24}
$$

evaluated under the same covariance matrices
$\boldsymbol{\Sigma}_\bullet$. The score-axis variance ratio is
$\rho_0 = \rho_1 = \rho_{\mathrm{ext}}$ when an external estimate is used,
and both batch-axis ratios are $\rho_{w0}$. The selection probability is

$$
p_{\mathrm{deno}}^{\mathrm{B}}
= \tau \left[
    \Phi\!\left(\frac{u/\sqrt{\rho_{w0}}}{\sqrt{\widehat{\mathbb{V}}(S_{\mathrm{bat}}) + \sigma^2}}\right)
    - \Phi\!\left(\frac{\ell/\sqrt{\rho_{w0}}}{\sqrt{\widehat{\mathbb{V}}(S_{\mathrm{bat}}) + \sigma^2}}\right)
  \right]
+ (1 - \tau)(1 - p_{\mathrm{cut}}). \tag{25}
$$

The integration limits in both branches are

$$
\ell = -z_{1-p_{\mathrm{cut}}/2}\,
       \sqrt{\widehat{\mathbb{V}}(S_{\mathrm{bat}})}\,
       \sqrt{\rho_{w0}},
\qquad
u = -\ell. \tag{26}
$$

### 4.3 SPA-calibrated score variance for the rectangle integration

The implementation does not insert the closed-form score variance into the
bivariate normal calculation. Instead, it first computes the marginal SPA
p-value $p_{\mathrm{SPA}}^{\mathrm{int}}$ of the internal score $S$ (see
Section 5) and inverts the chi-square calibration

$$
\widehat{\mathbb{V}}_\bullet(S)
= \frac{S^2}{\rho_\bullet\, \chi^2_{1, 1-p_{\mathrm{SPA}, \bullet}}}, \tag{27}
$$

where $\chi^2_{1, 1-p}$ denotes the upper-$p$ quantile of a one-degree-of-
freedom chi-square distribution and
$p_{\mathrm{SPA}, 0}, p_{\mathrm{SPA}, 1}$ are the marginal SPA p-values
under $\sigma^2 = 0$ and $\sigma^2 > 0$ respectively. The covariance
between $S$ and $S_{\mathrm{bat}}$ that enters
$\boldsymbol{\Sigma}_\bullet$ is rescaled by
$\sqrt{\widehat{\mathbb{V}}_\bullet(S) / V_{\mathrm{int},\bullet}}$ from its
closed-form value

$$
\widehat{\mathrm{Cov}}_\bullet(S_{\mathrm{bat}}, S)
= \left(\sum_i \widetilde w_i \widetilde R_i^{(b)}\right) \cdot 2\hat f(1-\hat f)
  + 2b\left(\sum_i R_i\right)
    \left\{\widehat{\mathbb{V}}(\hat f_{\mathrm{ext}}) + \sigma^2 \mathbb{I}\{\bullet = 1\}\right\}, \tag{28}
$$

where $\widetilde R_i^{(b)} = R_i - (1-b)\bar R$ and

$$
V_{\mathrm{int},\bullet}
= \sum_i (\widetilde R_i^{(b)})^2 \cdot 2\hat f(1-\hat f)
  + 4b^2\left(\sum_i R_i\right)^2
    \left\{\widehat{\mathbb{V}}(\hat f_{\mathrm{ext}}) + \sigma^2 \mathbb{I}\{\bullet = 1\}\right\} \tag{29}
$$

is the closed-form score variance under the same component. The rescaling
preserves the correlation between $S$ and $S_{\mathrm{bat}}$ while passing
the SPA-corrected scale of $S$ into the rectangle probability. Without
this step the bivariate normal tail would mix a Gaussian variance for $S$
with the saddlepoint-corrected mass that the user actually wants to
condition on.

## 5. Saddlepoint approximation with outlier split

The marginal p-value of $S$ that feeds (27) is obtained by saddlepoint
approximation under a diploid binomial model with allele frequency
$\hat f$, augmented by Gaussian terms that capture the additional sources
of variance entering $S$ from the external panel.

### 5.1 Residual centring and outlier split

WtCoxG centres the residuals once per evaluation by $bm = (1-b)\bar R$
and splits them by the $1.5 \times \mathrm{IQR}$ rule
(`detectOutliers` in
[src/util/outlier.hpp](../../src/util/outlier.hpp)) into the index sets
$\mathcal O$ (outliers) and $\mathcal N$ (non-outliers). The
non-outlier residuals are absorbed into a closed-form Gaussian
approximation; only the outlier residuals enter the per-iteration
empirical CGF sum.

### 5.2 Cumulant generating function

For a single non-missing genotype with $G \sim \mathrm{Binomial}(2,
\hat f)$, the moment generating function is
$M_G(t) = [(1-\hat f) + \hat f\, e^t]^2$ and the cumulant generating
function $K_G(t) = \log M_G(t)$ has derivatives that the implementation
evaluates in a single fused loop (see `outlierCgf_*` in
[src/wtcoxg/wtcoxg.cpp](../../src/wtcoxg/wtcoxg.cpp); scalar, AVX2, and
AVX-512 variants dispatch at runtime via `simdLevel()`). Writing
$r_i = R_i - bm$, the outlier contribution returns the triple

$$
K_0^{\mathcal O}(t) = \sum_{i \in \mathcal O} K_G(t r_i), \quad
K_1^{\mathcal O}(t) = \sum_{i \in \mathcal O} r_i K_G'(t r_i), \quad
K_2^{\mathcal O}(t) = \sum_{i \in \mathcal O} r_i^2 K_G''(t r_i). \tag{30}
$$

The non-outlier residuals contribute via their Gaussian second-order
expansion of $K_G$, which in closed form is

$$
K_0^{\mathcal N}(t) = m_{\mathcal N} t + \tfrac{1}{2} v_{\mathcal N} t^2,
\quad
K_1^{\mathcal N}(t) = m_{\mathcal N} + v_{\mathcal N} t,
\quad
K_2^{\mathcal N}(t) = v_{\mathcal N}, \tag{31}
$$

with the moments

$$
m_{\mathcal N} = 2\hat f \sum_{i \in \mathcal N}(R_i - bm)
              - 2b\hat f \sum_i R_i,
\quad
v_{\mathcal N}
= 2\hat f(1-\hat f) \sum_{i \in \mathcal N}(R_i - bm)^2. \tag{32}
$$

WtCoxG then augments the CGF with two scalar Gaussian corrections that
arise from the borrowed external estimate, written separately for the
first and second derivatives because they have distinct physical sources:

$$
v_{\mathrm{batch}}
= 4b^2 \left(\sum_i R_i\right)^2 \widehat{\mathbb{V}}(\hat f_{\mathrm{ext}}),
\quad
v_{\mathrm{finite}}
= \mathrm{obs\_ct} \left(\frac{\sum_i R_i}{n + n_{\mathrm{ext}}}\right)^2
  \hat f(1-\hat f). \tag{33}
$$

The total CGF and its first two derivatives are

$$
H(t) = K_0^{\mathcal O}(t) + m_{\mathcal N} t
     + \tfrac{1}{2}(v_{\mathcal N} + v_{\mathrm{batch}}) t^2, \tag{34}
$$

$$
H'(t) = K_1^{\mathcal O}(t) + m_{\mathcal N}
       + (v_{\mathcal N} + v_{\mathrm{batch}}) t, \tag{35}
$$

$$
H''(t) = K_2^{\mathcal O}(t) + v_{\mathcal N} + v_{\mathrm{finite}}. \tag{36}
$$

The asymmetry between $H'$ and $H''$ (one carries $v_{\mathrm{batch}}$,
the other $v_{\mathrm{finite}}$) is intentional: the batch-effect
variance enters the score expectation through the centred allele
frequency, whereas the finite-reference panel correction enters the
score variance through the per-marker allele-count uncertainty.

### 5.3 Saddlepoint and Lugannani-Rice tail

For a target $s = S/\rho_\bullet$, WtCoxG solves $H'(\zeta) = s$ by a
damped Newton-Raphson iteration that reuses the cached non-outlier
quantities (initial $\zeta_0 = 0$, tolerance $10^{-3}$, maximum
$100$ iterations). The Lugannani-Rice tail probability is

$$
\Pr(S/\rho_\bullet \le s)
\simeq \Phi\!\left(w + \frac{1}{w}\log\frac{v}{w}\right),
\quad
w = \mathrm{sign}(\zeta) \sqrt{2(\zeta s - H(\zeta))},
\quad
v = \zeta \sqrt{H''(\zeta)}. \tag{37}
$$

The two-sided marginal SPA p-value is

$$
p_{\mathrm{SPA}}
= \Pr\!\left(S/\rho_\bullet \ge |S|/\rho_\bullet\right)
+ \Pr\!\left(S/\rho_\bullet \le -|S|/\rho_\bullet\right), \tag{38}
$$

returned by `spaGOneSnpHomoFromScalars` in
[src/wtcoxg/wtcoxg.cpp](../../src/wtcoxg/wtcoxg.cpp). When
$|Z| = |S/\rho_\bullet| / \sqrt{V_{\mathrm{int},\bullet}}$ falls below
the SPA cutoff (`--spa-z-threshold`, default 2.0), the normal
approximation $p_{\mathrm{norm}} = 2\Phi(-|Z|)$ is returned in place of
the saddlepoint approximation.

## 6. Quality control gates

The external-reference branch is evaluated only when both alleles have
non-missing count $\ge 10$ at the marker (i.e. $\sum_i G_i \ge 10$ and
$2n - \sum_i G_i \ge 10$) and a batch-effect p-value is available. The
non-extended chunk-level QC filters
`--geno` (per-marker missing-rate cap),
`--maf`, `--mac`, and `--hwe` apply uniformly via the shared marker
engine; see [src/engine/marker.hpp](../../src/engine/marker.hpp).

## 7. Output columns

For each marker, `WtCoxGMethod::getResultVec` appends five
method-specific columns after the standard meta block
(`CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC HWE_P`):

| Column     | Definition                                                                                            |
|------------|-------------------------------------------------------------------------------------------------------|
| `p_ext`    | Conditional p-value using the external reference: $p_{\mathrm{con}}^{\mathrm{B}}$ from (23) when the marker passes the screen, $p_{\mathrm{con}}^{\mathrm{A}}$ from (20) when it does not. Reported as `NA` if the external branch is skipped (missing $p_{\mathrm{bat}}$, MAC below the threshold of Section 6, or $\hat V(S) \le 0$). |
| `p_noext`  | Internal-only p-value: marginal SPA p-value $p_{\mathrm{SPA}}$ from (38) computed at $b = 0$ with $\widehat{\mathbb{V}}(\hat f_{\mathrm{ext}}) = 0$ and $\rho_\bullet = \rho_{\mathrm{int}}$. |
| `z_ext`    | Score-axis z-score $Z = S / \sqrt{\widehat{\mathbb{V}}_1(S)}$ associated with `p_ext` (computed at the combined $\hat f$). |
| `z_noext`  | Score-axis z-score $Z = S / \sqrt{V_{\mathrm{int},0}}$ associated with `p_noext`. |
| `p_batch`  | Per-marker batch-effect p-value $p_{\mathrm{bat}}$ from (6).                                          |

The dispatcher entry points are `runWtCoxGPheno` (single phenotype) and
`runWtCoxG` (multi-phenotype) in
[src/wtcoxg/wtcoxg.cpp](../../src/wtcoxg/wtcoxg.cpp); both invoke the
shared multi-phenotype marker engine
(`multiPhenoEngine` in [src/engine/marker.cpp](../../src/engine/marker.cpp))
after the per-phenotype Phase 2 estimation has populated the per-marker
`WtCoxGRefInfo` table.
