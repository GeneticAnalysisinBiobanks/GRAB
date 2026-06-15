# SPAsqr

SPAsqr (Saddlepoint Approximation for smoothed Quantile Regression,
"SPA-squared") performs score tests for genetic variants on a quantitative
phenotype across a set of quantile levels $\tau$, then combines the
per-$\tau$ p-values into a single omnibus test. Relative to a mean-model
score test (SPACox / SPAmix / SPAGRM), it probes association at multiple
points of the conditional outcome distribution, which improves power when
the genetic effect is heteroscedastic or confined to a tail.

Let:

- $N \times 1$ vector $\mathbf{G}$ represent genotypes (allele counts) for a variant to be tested;
- $N \times p$ matrix $\mathbf{X}$ represent $p$ covariates (excluding intercept), and $\mathbf{Z} = [\mathbf{1}_N, \mathbf{X}]$ the design with intercept;
- $N \times 1$ vector $\mathbf{Y}$ represent the quantitative phenotype;
- $\boldsymbol{\tau} = (\tau_1, \ldots, \tau_m)$ represent the $m$ quantile levels (CLI `--spasqr-taus`);
- $N \times N$ sparse matrix $\mathbf{\Phi}$ represent the GRM, encoding cryptic relatedness;
- $\beta$ denote the genetic effect to be tested.

The null hypothesis at level $\tau$ is $H_0^{(\tau)}: \beta(\tau) = 0$, and the
omnibus null is the intersection $\bigcap_\tau H_0^{(\tau)}$.

This implementation is a pure C++17 / Eigen / Boost port of
`ref_code/src/mtSPAsqr.h` + `SPAsqr.R`. The score-mode path reuses the
SPAGRM saddlepoint machinery (`src/spasqr/spasqr.cpp`), one residual column
per $\tau$; the null model is fit by the QMME solver (`src/spasqr/qmme.cpp`).

## 1. Smoothed Quantile Regression Null Model (QMME)

### 1.1 Smoothed check loss

The classical quantile regression objective minimizes the check (pinball)
loss $\rho_\tau(u) = u(\tau - \mathbb{1}\{u < 0\})$, which is non-smooth at
$u = 0$. SPAsqr replaces it with a Gaussian-convolution-smoothed surrogate
of bandwidth $h$ (Heng and Wang, 2025):

$$
\ell_{h,\tau}(u) = \frac{h}{\sqrt{2\pi}} e^{-u^2/(2h^2)} + \frac{u}{2}\left[1 - 2\Phi\!\left(-\frac{u}{h}\right)\right] + \left(\tau - \tfrac{1}{2}\right) u \tag{1}
$$

whose first derivative is bounded, continuous, and recovers the check-loss
subgradient as $h \to 0$:

$$
\ell_{h,\tau}'(u) = \tau - \Phi\!\left(-\frac{u}{h}\right) \tag{2}
$$

The null objective for quantile level $\tau$, with residual
$r_i = Y_i - \mathbf{Z}_i^\top \boldsymbol{\beta}$, is the average loss

$$
f(\boldsymbol{\beta}) = \frac{1}{N} \sum_{i=1}^N \ell_{h,\tau}(Y_i - \mathbf{Z}_i^\top \boldsymbol{\beta}), \qquad
\nabla f(\boldsymbol{\beta}) = \frac{1}{N} \mathbf{Z}^\top\!\left[\Phi\!\left(-\frac{\mathbf{r}}{h}\right) - \tau\right] \tag{3}
$$

### 1.2 Quadratic Majorization–Minimization with Extrapolation

Because $\ell_{h,\tau}''(u) = \tfrac{1}{h}\phi(u/h) \le \tfrac{1}{h\sqrt{2\pi}}$,
the Hessian of $f$ is uniformly upper-bounded by a constant matrix, giving a
majorizer that is rebuilt only when the bandwidth changes:

$$
\mathbf{H} = \frac{1}{N\sqrt{2\pi}\,h}\,\mathbf{Z}^\top \mathbf{Z} + \delta \mathbf{I}, \qquad \delta = 10^{-6} \tag{4}
$$

The minimizer is found by MM with Nesterov-type extrapolation (paper
Algorithm 1; `qmme::SqrSolver::solve`):

$$
\mathbf{y}^{k} = \boldsymbol{\beta}^{k} + \frac{l}{l+2}\left(\boldsymbol{\beta}^{k} - \boldsymbol{\beta}^{k-1}\right), \qquad
\boldsymbol{\beta}^{k+1} = \mathbf{y}^{k} - \mathbf{H}^{-1}\nabla f(\mathbf{y}^{k}) \tag{5}
$$

The extrapolation counter $l$ restarts ($l \leftarrow 1$) whenever a step
ascends ($f(\boldsymbol{\beta}^{k+1}) > f(\boldsymbol{\beta}^{k})$) or after
`restartPeriod` $= 50$ steps; on ascent the iterate is additionally replaced
by a plain monotone MM step $\boldsymbol{\beta}^{k+1} = \boldsymbol{\beta}^{k} - \mathbf{H}^{-1}\nabla f(\boldsymbol{\beta}^{k})$,
preserving guaranteed descent. Internally $\mathbf{X}$ is standardized to
z-scores and $\mathbf{Y}$ is centered; the single Cholesky factor of
$\mathbf{H}$ is cached per phenotype and shared read-only across $\tau$
workers. Convergence is declared at $\lVert\nabla f\rVert_\infty \le$ `--spasqr-tol`.

### 1.3 Bandwidth selection

The bandwidth $h$ is set once per phenotype (`--spasqr-h` to fix it, else
auto):

$$
h = \frac{\operatorname{IQR}(\mathbf{Y})}{c}, \qquad c = \texttt{--spasqr-h-scale} \tag{6}
$$

with default $c = 3$ in score mode (and $c = 10$ in wald mode, §7, where the
kernel doubles as a conditional-density estimator). If the IQR-based value is
non-positive, a sample-size-dependent fallback $h = \max\!\big(((\ln N + p)/N)^{0.4},\, 0.05\big)$
is used.

## 2. Smoothed Residuals

For each $\tau$, the QMME fit yields $\hat{\boldsymbol{\beta}}(\tau)$ and raw
residuals $r_i = Y_i - \mathbf{Z}_i^\top \hat{\boldsymbol{\beta}}(\tau)$. The
quantity fed to the score test is the **smoothed residual**, equal to the
per-observation loss derivative (2):

$$
R(i, \tau) = \tau - \Phi\!\left(-\frac{r_i}{h}\right) = \ell_{h,\tau}'(r_i) \tag{7}
$$

Stacking over $\tau$ gives the residual matrix
$\mathbf{R} \in \mathbb{R}^{N \times m}$, $\mathbf{R}_{\cdot t} = R(\cdot, \tau_t)$.
At the QMME stationary point the gradient (3) vanishes, so each column is
approximately orthogonal to the design,

$$
\mathbf{Z}^\top \mathbf{R}_{\cdot t} \approx \mathbf{0}
\quad\Longrightarrow\quad
\sum_{i=1}^N R(i,\tau_t) \approx 0, \quad \mathbf{X}^\top \mathbf{R}_{\cdot t} \approx \mathbf{0} \tag{8}
$$

This is the SPAsqr analogue of the martingale-residual restriction
$\tilde{\mathbf{X}}^\top \mathbf{R} = \mathbf{0}$ in SPACox (cf. `spacox.md`
§1); it is what makes the score statistic (9) mean-zero under $H_0$ without
explicit per-marker centering.

## 3. Score Statistic and Variance (per $\tau$)

For each $\tau$ and each marker, the score statistic is the inner product of
the genotype and the smoothed residual column,

$$
S_\tau = \sum_{i=1}^N R(i,\tau)\, G_i = \mathbf{R}_{\cdot\tau}^\top \mathbf{G} \tag{9}
$$

computed for a batch of markers as the fused product
$\mathbf{R}^\top \mathbf{G}_{\text{batch}}$. Treating the genotype as a random
vector with covariance proportional to the GRM,
$\operatorname{Cov}(G_i, G_j) = \widehat{\operatorname{Var}}(G)\,\Phi_{ij}$,
the score variance is the residual quadratic form

$$
\widehat{\mathbb{V}}(S_\tau) = \widehat{\operatorname{Var}}(G)\,\big(\mathbf{R}_{\cdot\tau}^\top \mathbf{\Phi}\, \mathbf{R}_{\cdot\tau}\big)
\equiv G_{\text{var}} \cdot R\!\Phi R \tag{10}
$$

where $\widehat{\operatorname{Var}}(G)$ (`G_var`) is the empirical variance of
the post-imputation genotypes, and

$$
R\Phi R = \sum_{(i,j)\,\in\,\mathbf{\Phi}} c_{ij}\,\Phi_{ij}\,R(i,\tau)R(j,\tau), \qquad c_{ij} = \begin{cases}1 & i = j \\ 2 & i \neq j\end{cases} \tag{11}
$$

is accumulated once per $\tau$ over the stored lower-triangular GRM entries
(`SPAsqrPerTau::R_GRM_R`). The standardized score (z-statistic) is

$$
Z_\tau = \frac{S_\tau}{\sqrt{\widehat{\mathbb{V}}(S_\tau)}} \tag{12}
$$

Unlike SPAmix, SPAsqr uses a single overall minor-allele frequency
$\text{MAF} = \min(\widehat{f}, 1 - \widehat{f})$ per marker (homogeneous
population, as in SPAGRM), not individual-specific frequencies. Markers with
$\widehat{\mathbb{V}}(S_\tau) \le 0$ or $\text{MAF} \le 0$ return `NaN`.

## 4. Outlier Partitioning (per $\tau$ column)

For each residual column the sample is split into outliers $\mathcal{O}_\tau$
and non-outliers $\mathcal{N}_\tau$ by the Tukey IQR rule with an additional
absolute clamp specific to SPAsqr:

$$
\text{cutLo} = \max\!\big(Q_1 - k\cdot\operatorname{IQR},\, -B\big), \qquad
\text{cutHi} = \min\!\big(Q_3 + k\cdot\operatorname{IQR},\, +B\big) \tag{13}
$$

with $k = $ `--outlier-iqr-multiplier` (default $1.5$) and $B = $
`--spasqr-outlier-abs-bound` (default $0.55$); individual $i$ is an outlier in
column $\tau$ iff $R(i,\tau) < \text{cutLo}$ or $R(i,\tau) > \text{cutHi}$.
The absolute bound $B$ reflects that the smoothed residual (7) lives on a
bounded scale: $R(i,\tau) \in (\tau - 1, \tau)$. The partition is fixed across
all markers; per $\tau$ the code precomputes the outlier residual vector and
the non-outlier reductions
$\sum_{\mathcal{N}} R$ (`sum_R_nonOutlier`),
$R\Phi R\big|_{\mathcal{N}}$ (`R_GRM_R_nonOutlier`, restricting (11) to entries
with both endpoints non-outlier), and
$\lVert R_{\mathcal{O}}\rVert^2$ (`sum_unrelated_outliers2`).

## 5. Saddlepoint Approximation (per $\tau$)

### 5.1 Normal approximation in the bulk

If $|Z_\tau| \le$ `--spa-z-threshold` (default SPA cutoff $= 2$), the
two-sided p-value is the normal-approximation tail

$$
p_\tau = 2\,\Phi(-|Z_\tau|) \tag{14}
$$

(SPAsqr uses the non-strict comparison $|Z_\tau| \le$ cutoff.)

### 5.2 Variance-ratio rescaling of the score

The CGF in §5.3 is built from a model-based binomial genotype variance
$2\,\text{MAF}(1-\text{MAF})$, whereas the observed score variance (10) uses
the empirical $G_{\text{var}}$. To place the observed score on the CGF's
natural scale, SPAsqr applies a SAIGE-style variance-ratio correction. With
the CGF curvature at the origin

$$
K''(0) = 2\,\text{MAF}(1-\text{MAF})\big(R\Phi R\big|_{\mathcal{N}} + \lVert R_{\mathcal{O}}\rVert^2\big) \equiv \text{EmpVar} \tag{15}
$$

the rescaled score is

$$
S_\tau^{\text{adj}} = S_\tau \sqrt{\frac{\text{EmpVar}}{\widehat{\mathbb{V}}(S_\tau)}} \tag{16}
$$

This drops the implicit Hardy–Weinberg assumption
$G_{\text{var}} \approx 2\,\text{MAF}(1-\text{MAF})$ that an uncorrected score
would require. SPACox / SPAGRM in score-only form do not apply this rescaling;
it is the SPAsqr counterpart of the SPAmixPlus variance ratio $\rho$
(`spamix.md` §5.1).

### 5.3 Outlier-exact / non-outlier-Gaussian CGF

The total CGF of $S_\tau^{\text{adj}}$ partitions as in SPAmix (`spamix.md`
§4.5): outliers contribute an exact empirical CGF; non-outliers are folded
into a single Gaussian piece.

For an outlier $i$ with residual $R_i$, the diploid genotype
$G_i \sim \text{Binomial}(2,\text{MAF})$ gives the per-individual MGF
$[(1-\text{MAF}) + \text{MAF}\,e^{t R_i}]^2$, hence

$$
K_{\mathcal{O}}(t) = \sum_{i \in \mathcal{O}} 2\ln\!\big[(1-\text{MAF}) + \text{MAF}\,e^{t R_i}\big] \tag{17}
$$

evaluated together with its first two derivatives in a single fused loop
(scalar / AVX2 / AVX-512 dispatch via `outlierCgf`). The non-outlier piece is
the Gaussian CGF with the GRM-weighted cumulants

$$
\mu_{\mathcal{N}} = 2\,\text{MAF}\sum_{i \in \mathcal{N}} R_i, \qquad
\sigma_{\mathcal{N}}^2 = 2\,\text{MAF}(1-\text{MAF})\,R\Phi R\big|_{\mathcal{N}} \tag{18}
$$

so that

$$
K(t) = K_{\mathcal{O}}(t) + \mu_{\mathcal{N}}\,t + \tfrac{1}{2}\sigma_{\mathcal{N}}^2 t^2, \quad
K''(t) = K_{\mathcal{O}}''(t) + \sigma_{\mathcal{N}}^2 \tag{19}
$$

The GRM enters the non-outlier variance through $R\Phi R\big|_{\mathcal{N}}$,
capturing familial correlation among the central majority of the sample.

### 5.4 Saddlepoint root and Lugannani–Rice tail

For each tail at $s = \pm |S_\tau^{\text{adj}}|$, the saddlepoint
$\zeta$ solves $K'(\zeta) = s$ by safeguarded Newton iteration
(`scalarFastGetRoot`, initial guess $\zeta_1 = \min(|S_\tau^{\text{adj}}|/\widehat{\mathbb{V}}(S_\tau),\, 1.2)$,
with bisection fallback on sign changes). The Lugannani–Rice tail probability
is

$$
w = \operatorname{sgn}(\zeta)\sqrt{2(\zeta s - K(\zeta))}, \quad
v = \zeta\sqrt{K''(\zeta)}, \quad
\Pr(S < s) \simeq \Phi\!\left(w + \frac{1}{w}\ln\frac{v}{w}\right) \tag{20}
$$

and the two-sided p-value sums the upper tail at $+|S_\tau^{\text{adj}}|$ and
the lower tail at $-|S_\tau^{\text{adj}}|$:

$$
p_\tau = \Pr\!\big(S > +|S_\tau^{\text{adj}}|\big) + \Pr\!\big(S < -|S_\tau^{\text{adj}}|\big) \tag{21}
$$

## 6. Cauchy Combination Across $\tau$

The per-$\tau$ p-values $p_{\tau_1}, \ldots, p_{\tau_m}$ are aggregated into a
single omnibus statistic by the Cauchy combination test (Liu and Xie, 2020),
which is robust to the dependence among the $\tau$-specific statistics:

$$
T = \frac{1}{m}\sum_{t=1}^{m} \tan\!\Big[\big(\tfrac{1}{2} - p_{\tau_t}\big)\pi\Big], \qquad
P_{\text{CCT}} = \frac{1}{2} - \frac{\arctan(T)}{\pi} \tag{22}
$$

`NaN` entries are skipped; any $p_{\tau_t} \le 0$ short-circuits
$P_{\text{CCT}} = 0$; entries $\ge 1$ are clamped to $0.999$; overflow-safe
branches replace $\tan$ near its poles. $P_{\text{CCT}}$ is the headline
p-value of the method.

## 7. Wald Mode

Score mode (§1–§6) is the default. The wald mode (`--spasqr-mode wald`,
`src/spasqr/spasqr_wald.cpp`) instead reports effect sizes. For each marker
and each $\tau$ it refits the smoothed QR on the augmented design
$\mathbf{Z}^+ = [\mathbf{1}\,|\,\mathbf{X}\,|\,\mathbf{G}]$ by QMME, taking
$\hat{\beta}_G(\tau)$ as the last coefficient. The standard error is the
M-estimation sandwich variance evaluated at the fitted residual
$e_i = Y_i - \mathbf{Z}_i^{+\top}\hat{\boldsymbol{\theta}}$:

$$
\mathbf{B} = \frac{1}{N}\sum_i K_i\,\mathbf{Z}_i^+ \mathbf{Z}_i^{+\top}, \;\;
K_i = \frac{1}{h\sqrt{2\pi}} e^{-e_i^2/(2h^2)}; \qquad
\mathbf{M} = \frac{1}{N}\sum_i R_i^2\,\mathbf{Z}_i^+ \mathbf{Z}_i^{+\top}, \;\;
R_i = \tau - \Phi(-e_i/h) \tag{23}
$$

$$
\mathbf{V} = \frac{1}{N}\,\mathbf{B}^{-1}\mathbf{M}\,\mathbf{B}^{-1}, \qquad
\operatorname{SE}(\hat\beta_G) = \sqrt{V_{GG}}, \qquad
Z_\tau = \frac{\hat\beta_G(\tau)}{\operatorname{SE}(\hat\beta_G)} \tag{24}
$$

with $p_\tau = 2\big(1 - \Phi(|Z_\tau|)\big)$ and the same Cauchy combination
(22). The bread $\mathbf{B}$ is the kernel estimate of the conditional density
at the $\tau$-quantile (hence the larger default bandwidth scale $c = 10$);
the meat $\mathbf{M}$ is the empirical outer product of the score residuals.
No GRM is used. Because the per-marker QR refit is appreciably slower than the
score-mode dot product, the regression baseline restricts wald mode to 100
variants via `--extract` (see project `CLAUDE.md`).

## 8. Output Columns

For each marker, SPAsqr writes a wide, one-line-per-marker record appended to
the standard meta block
(`CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC HWE_P`). Column labels embed the
$\tau$ value (e.g. `P_tau0.5`).

| Column            | Definition |
|-------------------|------------|
| `P_CCT`           | Omnibus Cauchy-combined p-value (22) over all $\tau$. The headline result. |
| `P_tau{τ}`        | Per-$\tau$ two-sided p-value: normal approximation (14) when $\lvert Z_{\mathrm{Norm},\tau}\rvert \le$ cutoff, otherwise the saddlepoint p-value (21). |
| `Z_tau{τ}`        | Per-$\tau$ z-score consistent with `P_tau{τ}`: $\operatorname{sign}(Z_{\mathrm{Norm},\tau})\,\Phi^{-1}(1 - P_\tau/2)$ in score mode (equals `Z_Norm_tau{τ}` outside the SPA tail), or the Wald z-statistic (24) in wald mode. |
| `Z_Norm_tau{τ}`   | (score mode only) Per-$\tau$ raw standardized score (12) $S_\tau / \sqrt{\widehat{\mathbb{V}}(S_\tau)}$; not altered by SPA. |
| `BETA_tau{τ}`     | (wald mode only) Quantile-specific effect estimate $\hat\beta_G(\tau)$. |
| `SE_tau{τ}`       | (wald mode only) Sandwich standard error $\operatorname{SE}(\hat\beta_G)$ from (24). |

In score mode the emitted columns are `P_CCT`, then all `P_tau{τ}`, then all
`Z_tau{τ}`, then all `Z_Norm_tau{τ}` (grouped by statistic). SPA re-calibrates
only the per-$\tau$ `P` and the derived `Z` in the tails and does not alter
`Z_Norm`. Wald mode replaces the `Z_Norm_tau{τ}` group with `BETA_tau{τ}` and
`SE_tau{τ}` blocks; there `Z_tau{τ}` is already the p-consistent Wald z, so no
separate `Z_Norm` is emitted. A degenerate marker
($\widehat{\mathbb{V}}(S_\tau) \le 0$, monomorphic, or rejected by QC) reports
`NA` for the affected $\tau$ entries, which the Cauchy combination (22) skips.

## References

- Heng, Q. and Wang, Y. (2025). Smoothed quantile regression via Quadratic
  Majorization–Minimization with Extrapolation. *(QMME null-model solver.)*
- Liu, Y. and Xie, J. (2020). Cauchy combination test: a powerful test with
  analytic p-value calculation under arbitrary dependency structures.
  *J. Amer. Statist. Assoc.* 115(529), 393–402.
- Barndorff-Nielsen, O. E. Saddlepoint (Lugannani–Rice) tail approximation;
  see `spacox.md` §4 and `spamix.md` §5 for the diploid-genotype CGF reused
  here.
