# Null-model residuals in GRAB

All score-test methods in GRAB (SPACox, SPAGRM, SPAmix, SPAmixPlus,
SPAmixLocalPlus, WtCoxG, LEAF) consume a per-subject residual vector
$\mathbf{R} = (R_1, \ldots, R_n)^\top$ produced by an upstream null-model fit.
The score statistic for marker $j$ is

$$
S_j = \sum_{i=1}^n R_i\, G_{ij}, \tag{0}
$$

and the asymptotic null distribution of $S_j$ depends on the choice of
residual.  This document gives precise definitions of every residual that
GRAB constructs, together with the optimisation problem each one solves and
the analytical properties (mean-zero conditions, asymptotic mean and
variance under $H_0$) that justify its use in (0).

The implementations live in
[src/util/regression.{hpp,cpp}](../../src/util/regression.cpp) for the four
generic GLM-family residuals and in
[src/wtcoxg/regression.{hpp,cpp}](../../src/wtcoxg/regression.cpp) for the
case–control reweighting helper used by WtCoxG and LEAF.  Throughout this
document we let

- $n$ — number of subjects retained after NaN filtering;
- $p$ — number of columns in the design matrix $\mathbf{X}$ (intercept
  included when the model requires one);
- $w_i > 0$ — per-subject case weight passed to the fit;
- $\mathbf{X} \in \mathbb{R}^{n \times p}$ — covariate matrix
  (column 0 is the intercept whenever an intercept is part of the model);
- $\beta \in \mathbb{R}^{p}$ — regression coefficients estimated under the
  null;
- $y_i$ — phenotype value, with type depending on the residual family.

All four generic fits drop rows for which $y_i$, $w_i$, or any covariate is
$\mathrm{NaN}$ before fitting.  Constant covariate columns (other than the
intercept) and non-positive weights are rejected with a runtime error.

## 1. Linear residuals (quantitative trait)

Function: `regression::linearResiduals`
([src/util/regression.cpp:148](../../src/util/regression.cpp#L148)).

### 1.1 Model

For a quantitative response $y_i \in \mathbb{R}$ the null model is the
weighted linear regression

$$
y_i = \mathbf{X}_i^\top \beta + \varepsilon_i,
\qquad
\varepsilon_i \sim \big(0,\, \sigma^2 / w_i\big),
\qquad i = 1, \ldots, n. \tag{1.1}
$$

### 1.2 Estimator

The coefficient vector is the weighted ordinary-least-squares solution

$$
\hat\beta
\;=\;
\arg\min_{\beta \in \mathbb{R}^{p}}
\;\sum_{i=1}^n w_i \big(y_i - \mathbf{X}_i^\top \beta\big)^2
\;=\;
\big(\mathbf{X}^\top \mathbf{W} \mathbf{X}\big)^{-1}
\mathbf{X}^\top \mathbf{W}\, \mathbf{y}, \tag{1.2}
$$

with $\mathbf{W} = \operatorname{diag}(w_1, \ldots, w_n)$, computed in
closed form by an $\mathrm{LDL}^\top$ factorisation of
$\mathbf{X}^\top \mathbf{W} \mathbf{X}$.

### 1.3 Residual

$$
R_i \;=\; y_i - \mathbf{X}_i^\top \hat\beta. \tag{1.3}
$$

### 1.4 Mean-zero property

The first-order optimality condition is
$\mathbf{X}^\top \mathbf{W} \mathbf{R} = \mathbf{0}$.  If column 0 of
$\mathbf{X}$ is the intercept (the all-ones vector $\mathbf{1}_n$), then

$$
\sum_{i=1}^n w_i R_i = 0. \tag{1.4}
$$

When $w_i \equiv 1$ this reduces to $\sum_i R_i = 0$.

### 1.5 Use in (0)

Under $H_0: \beta_G = 0$ and the assumption that $G_{ij}$ is uncorrelated
with the residual conditional on $\mathbf{X}$, $S_j$ in (0) is asymptotically
normal with mean zero and variance
$\operatorname{Var}(S_j) = \sigma^2 \sum_i w_i^{-1} G_{ij}^2$ in the equally
weighted case.

## 2. Logistic residuals (binary trait)

Function: `regression::logisticResiduals`
([src/util/regression.cpp:400](../../src/util/regression.cpp#L400)).

### 2.1 Model

For $y_i \in \{0, 1\}$ the null model is the weighted binomial GLM with
canonical (logit) link:

$$
\operatorname{logit} \mu_i \;=\; \log \frac{\mu_i}{1 - \mu_i}
\;=\; \mathbf{X}_i^\top \beta,
\qquad
\mu_i \equiv \mathbb{P}(y_i = 1 \mid \mathbf{X}_i). \tag{2.1}
$$

The design matrix must contain an intercept column.

### 2.2 Estimator

$\hat\beta$ maximises the weighted binomial log-likelihood

$$
\ell(\beta)
\;=\;
\sum_{i=1}^n w_i \Big[\, y_i \log \mu_i + (1 - y_i) \log(1 - \mu_i) \,\Big],
\qquad
\mu_i = \operatorname{logistic}(\mathbf{X}_i^\top \beta). \tag{2.2}
$$

Maximisation is by Fisher scoring (iteratively reweighted least squares,
IRLS): at iteration $k$ with $\eta_i^{(k)} = \mathbf{X}_i^\top \beta^{(k)}$
and $\mu_i^{(k)} = \operatorname{logistic}(\eta_i^{(k)})$, define
working weights and pseudo-responses

$$
\tilde w_i^{(k)} \;=\; w_i\, \mu_i^{(k)} \big(1 - \mu_i^{(k)}\big),
\qquad
z_i^{(k)} \;=\; \eta_i^{(k)} + \frac{y_i - \mu_i^{(k)}}
                                    {\mu_i^{(k)}\big(1 - \mu_i^{(k)}\big)}, \tag{2.3}
$$

and solve the weighted normal equations
$\big(\mathbf{X}^\top \tilde{\mathbf{W}}^{(k)} \mathbf{X}\big) \beta^{(k+1)}
 = \mathbf{X}^\top \tilde{\mathbf{W}}^{(k)} \mathbf{z}^{(k)}$.  Iteration
terminates on relative deviance change
$|D^{(k)} - D^{(k-1)}| / (|D^{(k)}| + 0.1) < \mathrm{tol}$.

### 2.3 Residual

GRAB returns the **response residual**

$$
R_i \;=\; y_i - \hat\mu_i,
\qquad
\hat\mu_i = \operatorname{logistic}(\mathbf{X}_i^\top \hat\beta). \tag{2.4}
$$

This matches `residuals(glm(..., family = binomial()), type = "response")`
in R.

### 2.4 Mean-zero property

The score equation for $\beta$ at convergence is
$\mathbf{X}^\top (\mathbf{y} - \boldsymbol{\hat\mu}) = \mathbf{0}$ when
$w_i \equiv 1$.  With the intercept in $\mathbf{X}$ this yields

$$
\sum_{i=1}^n R_i \;=\; 0 \quad (\text{unweighted}); \tag{2.5}
$$

for general $w_i$ the analogous identity is
$\sum_i w_i' R_i = 0$ where $w_i' = w_i / \big[\mu_i(1 - \mu_i)\big]$ is the
iteration-converged IRLS weight, so $\sum_i R_i$ is generally not exactly
zero when weights are non-constant.

### 2.5 Use in (0)

Under $H_0: \beta_G = 0$, $S_j = \sum_i R_i G_{ij}$ is asymptotically
normal with mean zero; its variance is well approximated by
$\sum_i \hat\mu_i(1 - \hat\mu_i) (G_{ij} - \bar G_j)^2$ in the equally
weighted case.  The saddlepoint approximation in SPACox / SPAGRM /
SPAmix retains the exact CGF of the score under $H_0$.

## 3. Cox martingale residuals (time-to-event trait)

Function: `regression::coxResiduals`
([src/util/regression.cpp:198](../../src/util/regression.cpp#L198)).

### 3.1 Model

For each subject $i$ observe $(t_i, \delta_i)$ with $t_i \ge 0$ the
follow-up time and $\delta_i \in \{0, 1\}$ the event indicator.  The
Cox proportional-hazards model assumes

$$
\lambda(t \mid \mathbf{X}_i) \;=\; \lambda_0(t)\, \exp\!\big(\mathbf{X}_i^\top \beta\big), \tag{3.1}
$$

with $\lambda_0(t)$ left unspecified.  $\mathbf{X}$ carries no intercept
because any time-invariant intercept is absorbed into $\lambda_0$.

### 3.2 Estimator

$\hat\beta$ maximises the weighted Breslow partial log-likelihood

$$
\ell_{\mathrm{P}}(\beta)
\;=\;
\sum_{i: \delta_i = 1} w_i
\Big[\, \mathbf{X}_i^\top \beta
       - \log \sum_{\ell \in \mathcal{R}(t_i)} w_\ell\,
                \exp\!\big(\mathbf{X}_\ell^\top \beta\big)
\Big], \tag{3.2}
$$

where $\mathcal{R}(t_i) = \{\ell : t_\ell \ge t_i\}$ is the risk set.
Newton–Raphson with cumulative-sum bookkeeping (one forward pass per
iteration in descending-time order) is used, with the same
relative-loglik convergence rule as `survival::coxph(...)$control` in R.

After convergence the **Breslow estimator** of the baseline cumulative
hazard is

$$
\hat\Lambda_0(t)
\;=\;
\sum_{j:\, t_j \le t,\; \delta_j = 1}
\;\frac{\sum_{m \in \mathcal{T}_j} w_m \delta_m}
       {S_0(t_j)},
\qquad
S_0(t)
\;=\;
\sum_{\ell \in \mathcal{R}(t)} w_\ell\,
   \exp\!\big(\mathbf{X}_\ell^\top \hat\beta\big), \tag{3.3}
$$

where $\mathcal{T}_j$ is the set of subjects sharing event time $t_j$
(tie group).

### 3.3 Residual

The martingale residual is

$$
R_i \;=\; \delta_i - \hat\Lambda_0(t_i)\,
          \exp\!\big(\mathbf{X}_i^\top \hat\beta\big), \tag{3.4}
$$

equivalent to `residuals(coxph(...), type = "martingale")` in R.

### 3.4 Mean-zero property

By construction of the Breslow baseline,

$$
\sum_{i=1}^n w_i R_i \;=\; 0, \tag{3.5}
$$

independently of $\mathbf{X}$.  When $w_i \equiv 1$ this gives
$\sum_i R_i = 0$.  Furthermore, because the score equation
$\mathbf{U}(\hat\beta) = \mathbf{0}$ implies
$\sum_i (\mathbf{X}_i - \bar{\mathbf{X}}(t_i))\, R_i = \mathbf{0}$,
the residuals are uncorrelated (in the appropriate weighted inner product)
with each column of $\mathbf{X}$.

### 3.5 Use in (0)

Under $H_0: \beta_G = 0$, $S_j = \sum_i R_i G_{ij}$ is asymptotically
normal with mean zero; the variance is well approximated by
$\sum_i (G_{ij} - \bar G_j(t_i))^2$ in the unweighted case.  See
[spacox.md](spacox.md) for the projection that underlies the SPACox SPA.

## 4. Surrogate residual (ordinal trait)

Function: `regression::cumulativeLogitFit`
([src/util/regression.cpp:491](../../src/util/regression.cpp#L491)).

### 4.1 Model

Let $y_i \in \{0, 1, \ldots, J - 1\}$ be an ordinal response with
$J \ge 2$ ordered categories, observed jointly with covariates
$\mathbf{X}_i \in \mathbb{R}^{p}$.  The fixed-effects proportional-odds
cumulative-logit model assumes the existence of a latent continuous
response

$$
Y_i^\star \;=\; \mathbf{X}_i^\top \beta + \xi_i,
\qquad
\xi_i \overset{\mathrm{iid}}{\sim}\;\operatorname{Logistic}(0, 1),
\tag{4.1}
$$

discretised through fixed cutpoints
$-\infty = \varepsilon_{-1} < \varepsilon_0 < \varepsilon_1 < \cdots
       < \varepsilon_{J-2} < \varepsilon_{J-1} = +\infty$:

$$
Y_i = j
\iff
\varepsilon_{j-1} < Y_i^\star \le \varepsilon_j,
\qquad j = 0, 1, \ldots, J - 1. \tag{4.2}
$$

The induced category probabilities are
$\mathbb{P}(Y_i \le j \mid \mathbf{X}_i) = \nu_{ij}$, with

$$
\nu_{ij} \;=\; \operatorname{logistic}(\varepsilon_j - \mathbf{X}_i^\top \beta),
\qquad
\nu_{i,-1} = 0,
\quad \nu_{i,J-1} = 1, \tag{4.3}
$$

$$
\mu_{ij} \;=\; \nu_{ij} - \nu_{i,j-1}
\;=\; \mathbb{P}(Y_i = j \mid \mathbf{X}_i). \tag{4.4}
$$

### 4.2 Estimator

Let $\theta = (\varepsilon_0, \ldots, \varepsilon_{J-2}, \beta^\top)^\top
\in \mathbb{R}^{J - 1 + p}$.  The conditional log-likelihood is

$$
\ell(\theta)
\;=\;
\sum_{i = 1}^{n}
\sum_{j = 0}^{J - 1} \mathbf{1}\{Y_i = j\}\, \log \mu_{ij}. \tag{4.5}
$$

GRAB maximises $\ell$ by Fisher scoring with the outer-product (OPG)
approximation to the Hessian
$\mathbf{H} \approx \sum_i \mathbf{s}_i \mathbf{s}_i^\top$, where
$\mathbf{s}_i = \nabla_\theta \ell_i$.  A diagonal ridge
$10^{-6}\, \mathbf{I}$ is added to $\mathbf{H}$ for numerical stability;
the threshold ordering constraint
$\varepsilon_0 < \cdots < \varepsilon_{J-2}$ is restored after each
Newton step by $\varepsilon_j \gets \max(\varepsilon_j,
\varepsilon_{j-1} + 10^{-2})$.  This optimisation target is identical to
`MASS::polr(method = "logistic")` and `ordinal::clm(link = "logit")` in
R; the optimiser differs (R uses BFGS) and yields the same identified
$\nu_{ij}$ surface up to ~$10^{-5}$ in the present implementation.

### 4.3 Surrogate residual

GRAB constructs the **conditional surrogate residual** of Liu & Zheng
(2018, *JASA* 113(522), 845–854).  Given the fitted parameters
$(\hat\beta, \hat\varepsilon)$ and the observed category
$Y_i = y_i$, draw a latent value $\tilde Y_i^\star$ from the conditional
posterior of $Y_i^\star$ under (4.1)–(4.2):

$$
\tilde Y_i^\star
\mid Y_i = y_i, \mathbf{X}_i
\;\sim\;
\operatorname{Logistic}\!\big(\mathbf{X}_i^\top \hat\beta,\, 1\big)
\;\text{truncated to}\;
(\hat\varepsilon_{y_i - 1},\, \hat\varepsilon_{y_i}\,]. \tag{4.6}
$$

Apply the probability integral transform with the logistic CDF
$F(t \mid \mathbf{X}_i) =
 \operatorname{logistic}(t - \mathbf{X}_i^\top \hat\beta)$ to define the
**raw surrogate residual**

$$
\tilde R_i
\;=\; F(\tilde Y_i^\star \mid \mathbf{X}_i) - \tfrac{1}{2}. \tag{4.7}
$$

Under (4.6), $F(\tilde Y_i^\star \mid \mathbf{X}_i)$ is Uniform on
$\big(\hat\nu_{i, y_i - 1},\, \hat\nu_{i, y_i}\big]$, so the equivalent
direct sampler is

$$
\tilde R_i
\;=\;
\hat\nu_{i, y_i - 1}
\;+\; U_i \big(\hat\nu_{i, y_i} - \hat\nu_{i, y_i - 1}\big)
\;-\; \tfrac{1}{2},
\qquad
U_i \overset{\mathrm{iid}}{\sim}\;\operatorname{Uniform}(0, 1). \tag{4.8}
$$

Finally GRAB subtracts the algebraic sample mean to enforce
$\sum_i R_i = 0$ exactly:

$$
R_i \;=\; \tilde R_i - \frac{1}{n}\sum_{i'=1}^{n} \tilde R_{i'}. \tag{4.9}
$$

### 4.4 Algorithm (reproducible specification)

The following pseudocode, given the fitted $(\hat\beta, \hat\varepsilon)$
and a 32-bit unsigned seed $s$, produces the exact residual vector
$R$ emitted by GRAB.  No quantity depends on the order in which
subjects are presented other than through the deterministic mapping
$i \mapsto U_i$.

```text
Input:  X (n × p), y (n × 1, integer in {0..J-1}),
        β̂ (p × 1), ε̂ (J-1 × 1, strictly increasing), seed s.
Output: residual vector R (n × 1).

1.  Initialise an MT19937 engine with seed s and let
        U_1, U_2, ..., U_n
    be the first n outputs of generate_canonical<double, 53, mt19937>
    (see §4.5 below for the bit-exact specification).

2.  For each i = 1, ..., n:
       η_i ← X_i^⊤ β̂.
       if y_i = 0:        F_low_i ← 0
       else:              F_low_i ← logistic(ε̂_{y_i - 1} − η_i).
       if y_i = J - 1:    F_hi_i  ← 1
       else:              F_hi_i  ← logistic(ε̂_{y_i}     − η_i).
       R_i ← F_low_i + U_i · (F_hi_i − F_low_i) − 1/2.

3.  R ← R − mean(R).

4.  Return R.
```

The implementation in [src/util/regression.cpp:587-606](../../src/util/regression.cpp#L587-L606) follows this pseudocode line-for-line.  The
reference R re-implementation in
[examples_1kg/vs_develop/grab_develop.R](../../examples_1kg/vs_develop/grab_develop.R)
(function `fit_ordinal_resid`) does the same and agrees with feat/cpp at
the level of optimiser tolerance
($\mathrm{med} |\Delta R| \approx 2 \times 10^{-6}$,
 $\max |\Delta R| \approx 2 \times 10^{-5}$ on the 1KG sample).

### 4.5 Bit-exact Uniform(0,1) sequence

The uniform stream $U_1, U_2, \ldots$ is the canonical reduction of
`std::mt19937(s)` (Knuth TGFSR initialisation: $w_0 = s$,
$w_k = (1812433253 \cdot (w_{k-1} \oplus (w_{k-1} \gg 30)) + k) \bmod 2^{32}$
for $k = 1, \ldots, 623$; then the standard MT19937 twist and tempering)
combined via `std::generate_canonical<double, 53, mt19937>`:

$$
U_i
\;=\;
\frac{d^{(2i-1)} + d^{(2i)} \cdot 2^{32}}{2^{64}},
\qquad i = 1, 2, \ldots, \tag{4.10}
$$

where $d^{(k)}$ is the $k$-th tempered 32-bit output of the MT19937
engine.  This matches libstdc++'s implementation of
`std::uniform_real_distribution<double>` on a 32-bit URBG with the
default 53-bit mantissa target ($\lceil 53 / 32 \rceil = 2$ draws per
sample, lower-weight draw first).  The R reference implementation
([examples_1kg/vs_develop/grab_develop.R](../../examples_1kg/vs_develop/grab_develop.R),
function `uniform_mt19937`) reproduces (4.10) bit-for-bit using only
base R arithmetic by storing the 32-bit state as `double` values in
$[0, 2^{32})$ and decomposing every XOR / AND into two 16-bit halves.

### 4.6 Properties under H₀

Under (4.1)–(4.2) with $(\beta, \varepsilon)$ at their true values, the
probability-integral-transform identity gives
$F(Y_i^\star \mid \mathbf{X}_i) \sim \operatorname{Uniform}(0, 1)$
unconditionally; the truncated draw $\tilde Y_i^\star$ defined in (4.6)
inherits this distribution by the tower property of conditional
expectation:

$$
\tilde R_i \mid \mathbf{X}_i
\;\sim\;
\operatorname{Uniform}\!\big(-\tfrac{1}{2},\, +\tfrac{1}{2}\big),
\qquad
\tilde R_i \perp\!\!\!\perp \mathbf{X}_i. \tag{4.11}
$$

In particular $\mathbb{E}[\tilde R_i \mid \mathbf{X}_i] = 0$ and
$\operatorname{Var}(\tilde R_i \mid \mathbf{X}_i) = 1/12$ exactly, and
$\tilde R_i$ is independent of any function of $\mathbf{X}_i$ —
including a candidate marker $G_i$ under the null of no genetic effect.
Consequently

$$
\mathbb{E}_{H_0}\!\big[S_j\big] = 0,
\qquad
\operatorname{Var}_{H_0}\!\big[S_j\big]
\;=\;
\tfrac{1}{12}\, \sum_{i = 1}^{n} G_{ij}^2
\;+\;
\big(\text{contribution from }\mathbf{X}\text{-projection of }G_j\big),
\tag{4.12}
$$

and $S_j$ is asymptotically normal under $H_0$ — the configuration the
SPACox / SPAGRM / SPAmix saddlepoint approximations downstream are
designed for.

### 4.7 Reproducibility caveats

Because (4.7)/(4.8) involves randomness, two GRAB runs with `seed = 0`
(the CLI default, mapping to `std::random_device`) produce ordinal
residuals — and therefore ordinal p-values — that differ by
$O(1 / \sqrt{n})$ across runs.  To obtain reproducible ordinal output,
the user must pass `--seed <s>` with $s \ne 0$; the same $s$ supplied to
the R reference implementation regenerates the same residual vector to
within the optimiser-convergence tolerance of $\hat\nu$.

The arithmetic centring (4.9) is required because $\sum_i \tilde R_i$ is
$O(1/\sqrt n)$ but not exactly zero in any finite sample.  Centring does
not bias the asymptotic null distribution since
$\sum_i \tilde R_i / n \xrightarrow{P} 0$ under (4.11).

## 5. Case–control sampling-weight correction (WtCoxG / LEAF)

Function: `regression::calRegrWeight`
([src/wtcoxg/regression.cpp:17](../../src/wtcoxg/regression.cpp#L17)).

WtCoxG and LEAF combine the residuals of §2 and §3 with a per-subject
sampling weight $w_i$ that compensates for case ascertainment.  Let
$\pi$ be the assumed population disease prevalence (CLI flag
`--prevalence`) and let $\delta_i \in \{0, 1\}$ be the indicator that
subject $i$ is a case.  Define

$$
n_{\text{case}} = \sum_i \delta_i,
\qquad
n_{\text{ctrl}} = n - n_{\text{case}},
\qquad
r_{\text{samp}} = \frac{n_{\text{case}}}{n_{\text{ctrl}}},
\qquad
r_{\text{pop}}  = \frac{\pi}{1 - \pi}. \tag{5.1}
$$

The per-subject weight is

$$
w_i
\;=\;
\begin{cases}
\;1 & \text{if } \delta_i = 1 \text{ (cases unweighted)}, \\[2pt]
\;\dfrac{r_{\text{samp}}}{r_{\text{pop}}} & \text{if } \delta_i = 0
   \text{ (controls down-weighted to match population odds)}.
\end{cases} \tag{5.2}
$$

The factor $r_{\text{samp}} / r_{\text{pop}}$ aligns the case-to-control
odds in the weighted sample with the assumed population odds, so that
weighted score quantities computed from the residuals of §2 or §3
approximate their population-level counterparts.

### 5.1 WtCoxG residual

WtCoxG ([src/wtcoxg/wtcoxg.cpp:1483](../../src/wtcoxg/wtcoxg.cpp#L1483))
fits a weighted Cox model on $(t_i, \delta_i, \mathbf{X}_i, w_i)$ using
$\texttt{regression::coxResiduals}$ of §3 with the weights (5.2), and
returns the resulting martingale residual $R_i$ of (3.4).  When the
phenotype is binary rather than time-to-event,
[src/wtcoxg/wtcoxg.cpp:1497](../../src/wtcoxg/wtcoxg.cpp#L1497) falls
back to weighted logistic regression and returns the response residual
of (2.4).

The mean-zero property of (3.5) holds in the WtCoxG-weighted form
$\sum_i w_i R_i = 0$, with $w_i$ given by (5.2).

### 5.2 LEAF residual

LEAF partitions the $n$ subjects into $K$ ancestry clusters by k-means on
principal components (cluster labels
$c_i \in \{1, \ldots, K\}$).  For each cluster $c \in \{1, \ldots, K\}$
let $\mathcal{I}_c = \{i : c_i = c\}$ be its members.  LEAF computes a
**per-cluster** weight vector
$(w_i)_{i \in \mathcal{I}_c}$ by applying (5.2) to the within-cluster
case / control indicator
$(\delta_i)_{i \in \mathcal{I}_c}$, then fits a within-cluster weighted
GLM (binary trait, §2) or Cox model (time-to-event trait, §3) on
$(\mathbf{X}_i)_{i \in \mathcal{I}_c}$ to obtain a within-cluster residual
$(R_i)_{i \in \mathcal{I}_c}$.  Implementation:
[src/wtcoxg/leaf.cpp:742](../../src/wtcoxg/leaf.cpp#L742) for the weight,
[src/wtcoxg/leaf.cpp:751](../../src/wtcoxg/leaf.cpp#L751) and
[src/wtcoxg/leaf.cpp:757](../../src/wtcoxg/leaf.cpp#L757) for the
within-cluster Cox / logistic fits.

Within each cluster the mean-zero property of §2 / §3 holds in the
weighted form
$\sum_{i \in \mathcal{I}_c} w_i R_i = 0$, but $\sum_i R_i$ across clusters
is not, in general, zero.  The meta-analysis combination in
[leaf.cpp](../../src/wtcoxg/leaf.cpp) handles this by combining
per-cluster score statistics and their per-cluster variances rather than
by pooling the residuals.

## 6. Comparison with the develop-branch R reference

The reference workflow in
[examples_1kg/vs_develop/grab_develop.R](../../examples_1kg/vs_develop/grab_develop.R)
fits each null model in R and writes the residual vector to disk; feat/cpp
emits the same per-phenotype vectors via `--save-resid` to
`<out>.null.resid`.  On the 1000 Genomes test fixture
(`examples_1kg/1kg_auto.bed`, $n = 2504$), with feat/cpp invoked under
`--method SPACox --pheno-name Normal,Uniform,Binary1,Binary2,Exp1:Binary1,Exp2:Binary2,Binomial,Poisson --seed 12345 --threads 1`, the two implementations agree as follows:

| Phenotype | Trait kind | R reference | $\mathrm{med}\,|\Delta R|$ | $\max\,|\Delta R|$ | $\operatorname{cor}(R_{\text{R}}, R_{\text{C++}})$ |
| --- | --- | --- | ---: | ---: | ---: |
| Normal       | Quantitative   | `residuals(lm(y ~ X))`                                   | $5.6 \times 10^{-17}$ | $2.3 \times 10^{-15}$ | $1.000000$ |
| Uniform      | Quantitative   | `residuals(lm(y ~ X))`                                   | $1.4 \times 10^{-17}$ | $4.7 \times 10^{-16}$ | $1.000000$ |
| Binary1      | Binary         | `residuals(glm(., family = binomial()), type = "response")` | $1.3 \times 10^{-13}$ | $2.3 \times 10^{-11}$ | $1.000000$ |
| Binary2      | Binary         | same as above                                            | $8.2 \times 10^{-15}$ | $1.9 \times 10^{-14}$ | $1.000000$ |
| Exp1\_Binary1 | Time-to-event | `residuals(coxph(Surv(t, δ) ~ X), type = "martingale")` | $2.4 \times 10^{-6}$ | $4.2 \times 10^{-3}$ | $1.000000$ |
| Exp2\_Binary2 | Time-to-event | same as above                                          | $4.2 \times 10^{-17}$ | $4.4 \times 10^{-16}$ | $1.000000$ |
| Binomial     | Ordinal        | `fit_ordinal_resid` per §4.4 + `MASS::polr`              | $1.7 \times 10^{-6}$ | $2.0 \times 10^{-5}$ | $1.000000$ |
| Poisson      | Ordinal        | same as above                                            | $1.7 \times 10^{-6}$ | $1.6 \times 10^{-5}$ | $1.000000$ |

The non-zero residuals on Binary1, Exp1\_Binary1, and the two ordinal
phenotypes are entirely attributable to optimiser convergence
differences (GLM IRLS, Cox Newton, and BFGS-vs-Newton-OPG for ordinal),
not to algorithmic disagreement.  Three implementation choices secure
the perfect rank correlation across all eight phenotypes:

1. The R reference's ordinal residual is the **surrogate residual** of
   §4.3, not the linear-scale $y_i - \mathbb{E}[Y_i \mid \mathbf{X}_i]$
   that earlier versions of `grab_develop.R` used.  The latter is a
   valid score-test input but is not the residual feat/cpp emits.

2. R's ordinal residual draws its uniform stream $U_i$ from a base-R
   re-implementation of `std::mt19937(s)` and
   `std::generate_canonical<double, 53, mt19937>` (see §4.5 and the
   `uniform_mt19937` function inside `grab_develop.R`), so for any
   non-zero seed both implementations consume the **same** 53-bit
   Uniform(0, 1) sequence.

3. `feat/cpp` exposes the seed through `--seed` (CLI) →
   `EngineOptions::seed` →
   `regression::cumulativeLogitFit(..., seed)`.  The R reference reads
   the same seed from the environment variable `CUMULATIVELOGITFIT_SEED`
   (default 12345), which is what `grab2.sh` is configured to pass.
