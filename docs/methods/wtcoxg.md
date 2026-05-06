# WtCoxG

WtCoxG performs marker-level score tests using residuals from an upstream weighted survival null model, with optional borrowing of external allele-frequency information.

Let:

- $n \times 1$ vector $\mathbf{G} = (G_1, \ldots, G_n)^\top$ represent genotypes (allele counts) for a variant;
- $n \times 1$ vector $\mathbf{R} = (R_1, \ldots, R_n)^\top$ represent residuals from the upstream null fit;
- $n \times 1$ vector $\mathbf{w} = (w_1, \ldots, w_n)^\top$ represent subject weights used in the batch-effect screen;
- $n \times n$ sparse matrix $\mathbf{\Phi}$ represent the sparse GRM when available;
- $\hat f_0$ and $\hat f_1$ denote the control and case allele-frequency estimates from the study sample;
- $n_0$ and $n_1$ denote the non-missing control and case sample counts for the marker;
- $\hat f_{\mathrm{int}} = \frac{1}{2n}\sum_{i=1}^n G_i$ denote the internal allele-frequency estimate after mean imputation of missing genotypes;
- $\hat f_{\mathrm{ext}}$ denote the external allele-frequency estimate from the reference panel;
- $b \in [0, 1]$ denote the external-information weight;
- $\hat f = (1-b)\hat f_{\mathrm{int}} + b\hat f_{\mathrm{ext}}$ denote the combined allele-frequency estimate;
- $n_{\mathrm{ext}}$ denote the external sample size proxy, so that $2n_{\mathrm{ext}}$ is the external allele count;
- $\rho$ denote the variance-ratio correction used in the score test.

The null hypothesis is $H_0: \beta = 0$.

Missing genotypes are imputed by the observed mean genotype, i.e. $G_i \leftarrow 2\hat f_{\mathrm{int}}$ for missing entries. WtCoxG does not flip alleles to force MAF $\le 0.5$ in the marker-level test.

## 1. Batch-effect Screen and Parameter Estimation

WtCoxG first performs a marker-level batch-effect screen that compares the internally weighted allele frequency against the external reference.

Let

$$
e_r = \frac{n_1}{n_0 + n_1},
\qquad
w_0 = \frac{1-\pi}{\pi} \cdot \frac{e_r}{1-e_r},
\qquad
w_1 = 1, \tag{1}
$$

where $\pi$ is the reference prevalence supplied by the user.

The internally weighted allele-frequency estimate is

$$
\hat f_{\mathrm{w}} = \frac{\hat f_0 w_0 n_0 + \hat f_1 w_1 n_1}{w_0 n_0 + w_1 n_1}. \tag{2}
$$

The pooled allele-frequency estimate used in the batch-effect variance is

$$
\hat f_{\mathrm{pool}} = \frac{\hat f_0 w_0 n_0 + \hat f_1 w_1 n_1 + \hat f_{\mathrm{ext}} w_0 n_{\mathrm{ext}}}{w_0 n_0 + w_1 n_1 + w_0 n_{\mathrm{ext}}}. \tag{3}
$$

The unadjusted batch-effect z-statistic is

$$
Z_{\mathrm{bat}} = \frac{\hat f_{\mathrm{w}} - \hat f_{\mathrm{ext}}}{\sqrt{V_{\mathrm{bat}}}}, \tag{4}
$$

with

$$
V_{\mathrm{bat}} = \left\{ \frac{n_1 w_1^2 + n_0 w_0^2}{2(n_1 w_1 + n_0 w_0)^2} + \frac{1}{2n_{\mathrm{ext}}} \right\}
\hat f_{\mathrm{pool}} (1-\hat f_{\mathrm{pool}}). \tag{5}
$$

If a sparse GRM is available, WtCoxG applies a variance-ratio correction $\rho_{w0}$ and uses

$$
Z_{\mathrm{bat}}^{\mathrm{adj}} = \frac{Z_{\mathrm{bat}}}{\sqrt{\rho_{w0}}},
\qquad
p_{\mathrm{bat}} = 2\Phi\left(-\left| Z_{\mathrm{bat}}^{\mathrm{adj}} \right|\right). \tag{6}
$$

The quantity $\rho_{w0}$ in (6) is computed as

$$
\rho_{w0} = \frac{\tilde{\mathbf{w}}^\top \mathbf{\Phi} \tilde{\mathbf{w}} + \frac{1}{2n_{\mathrm{ext}}}}{\tilde{\mathbf{w}}^\top \tilde{\mathbf{w}} + \frac{1}{2n_{\mathrm{ext}}}}, \qquad \tilde w_i = \frac{w_i}{2\sum_{j=1}^n w_j}. \tag{7}
$$

If no sparse GRM is used, then $\rho_{w0} = 1$ in (6).

Only markers with $p_{\mathrm{bat}} \ge p_{\mathrm{cut}}$ proceed to the external-reference test.

### 1.1 Estimation of $\tau$ and $\sigma^2$

Markers are grouped by internal allele frequency, and within each MAF group WtCoxG estimates $\tau$ and $\sigma^2$ by matching the empirical distribution of batch-effect p-values. Here $\tau$ is the mixture weight of the batch-affected component, corresponding to the code variable `TPR`.

For cutoffs $c_j \in \{0.01, 0.11, 0.21, 0.31\}$, define the empirical pass rates

$$
\hat p_{\mathrm{den}, j} = \frac{1}{m} \sum_{k=1}^m I\left(p_{\mathrm{bat}, k} > c_j\right). \tag{8}
$$

Let $V_{\mathrm{Sbat}}$ denote the model-based variance of the batch statistic in that MAF group. WtCoxG models the pass probability by

$$
p_j(\tau, \sigma^2) = \tau \left[ \Phi\left(\frac{u_j}{\sqrt{V_{\mathrm{Sbat}} + \sigma^2}}\right) - \Phi\left(\frac{\ell_j}{\sqrt{V_{\mathrm{Sbat}} + \sigma^2}}\right) \right] + (1-\tau)(1-c_j). \tag{9}
$$

where

$$
\ell_j = -z_{1-c_j/2}\sqrt{V_{\mathrm{Sbat}}}, \qquad u_j = z_{1-c_j/2}\sqrt{V_{\mathrm{Sbat}}}. \tag{10}
$$

The implementation estimates $(\tau, \sigma^2)$ by minimizing

$$
\sum_j \left( \frac{\hat p_{\mathrm{den}, j} - p_j(\tau, \sigma^2)}{\hat p_{\mathrm{den}, j}} \right)^2. \tag{11}
$$

using Nelder-Mead optimization, then truncates both estimates to $[0, 1]$.

Thus, $\sigma^2$ in (16) is not a closed-form estimator. It is a per-MAF-group variance component learned from the batch-effect p-value distribution through the optimization in (11).

### 1.2 Estimation of $b$

The external-information weight $b$ is also estimated separately in each MAF group. It is not fixed a priori.

For a candidate $b \in [0, 1]$, WtCoxG defines a power proxy based on the conditional external-reference p-value. More specifically, it searches for the smallest case allele frequency $\mu_1$ such that the conditional p-value reaches the genome-wide threshold $5 \times 10^{-8}$.

For each candidate pair $(b, \mu_1)$, let

- $p_0(b, \mu_1)$ denote the bivariate normal rectangle probability under the model with $\sigma^2 = 0$;
- $p_1(b, \mu_1)$ denote the corresponding bivariate normal rectangle probability under the model with the estimated $\sigma^2 > 0$;
- $p_{\mathrm{deno}}(b, \mu_1)$ denote the selection probability of passing the batch-effect screen under the same candidate pair.

These quantities are evaluated by the same formulas later written explicitly in Section 6, but here they are viewed as functions of the trial parameters $(b, \mu_1)$ rather than as the final values at the observed marker.

The conditional p-value is treated as a scalar function of $(b, \mu_1)$:

$$
p_{\mathrm{con}}(b, \mu_1) = \frac{2\{\tau p_1(b, \mu_1) + (1-\tau)p_0(b, \mu_1)\}}{p_{\mathrm{deno}}(b, \mu_1)}. \tag{12}
$$

Here the control allele frequency is fixed at the current MAF-group center, and only $(b, \mu_1)$ varies across the optimization.

For each fixed $b$, WtCoxG then solves the scalar equation in $\mu_1$

$$
p_{\mathrm{con}}(b, \mu_1) = 5 \times 10^{-8}. \tag{13}
$$

For each $b$, let $\mu_1^{\star}(b)$ denote the smallest root of (13), i.e. the smallest case allele frequency that makes the conditional external-reference p-value equal to $5 \times 10^{-8}$. If (13) has multiple roots in the admissible search interval, $\mu_1^{\star}(b)$ means the numerically smallest one. WtCoxG then chooses the $b$ that minimizes this required case allele frequency:

$$
\hat b = \arg\min_{0 \le b \le 1} \mu_1^{\star}(b). \tag{14}
$$

The optimization is again numerical. In the code, $b$ is obtained by one-dimensional Brent minimization.

## 2. Residual-free Quantities

Before introducing the residual-based score, WtCoxG uses two variance components that do not depend on $\mathbf{R}$.

Under the diploid binomial model, the genotype variance at allele frequency $\hat f$ is

$$
\widehat{\mathbb{V}}(G_i \mid \hat f) = 2\hat f(1-\hat f). \tag{15}
$$

WtCoxG models the uncertainty of the external frequency by

$$
\widehat{\mathbb{V}}(\hat f_{\mathrm{ext}}) = \frac{\hat f(1-\hat f)}{2n_{\mathrm{ext}}} + \sigma^2. \tag{16}
$$

Here $\sigma^2$ is the per-MAF-group estimate from (11). It is set to zero in the internal-only branch and in the first-stage external calibration.

## 3. Score Statistic

WtCoxG centers the genotype by the combined allele frequency and forms the score statistic

$$
S = \sum_{i=1}^n R_i (G_i - 2\hat f). \tag{17}
$$

Equivalently,

$$
S = \mathbf{R}^\top (\mathbf{G} - 2\hat f\mathbf{1}_n). \tag{18}
$$

The internal-only analysis is the special case $b = 0$, so that $\hat f = \hat f_{\mathrm{int}}$.

The implementation then applies a variance-ratio correction to the score numerator and uses

$$
\frac{S}{\rho}. \tag{19}
$$

where $\rho$ is:

- $\rho = \mathrm{var\_ratio}_{\mathrm{int}}$ for the internal-only test when sparse-GRM adjustment is available;
- $\rho = \mathrm{var\_ratio}_{\mathrm{ext}}$ for the external-reference test;
- $\rho = 1$ when no variance-ratio adjustment is applied.

The external-reference branch is evaluated only when the marker passes the batch-effect screen and both alleles have MAC at least 10. Otherwise the external result is reported as missing.

## 4. Score Variance

Define

$$
\bar R = \frac{1}{n}\sum_{i=1}^n R_i, \qquad \tilde R_i = R_i - (1-b)\bar R. \tag{20}
$$

The score variance used by the SPA routine is the estimated variance of the unadjusted score $S$ itself:

$$
\widehat{\mathbb{V}}(S) = \sum_{i=1}^n \tilde R_i^2 \cdot 2\hat f(1-\hat f) + 4b^2\left(\sum_{i=1}^n R_i\right)^2 \widehat{\mathbb{V}}(\hat f_{\mathrm{ext}}). \tag{21}
$$

The corresponding z-score is

$$
Z = \frac{S / \rho}{\sqrt{\widehat{\mathbb{V}}(S)}}. \tag{22}
$$

Therefore, $\widehat{\mathbb{V}}(S)$ is not the estimated variance of $S / \rho$. If one wanted the variance of the adjusted score under the same deterministic-ratio approximation, it would be $\widehat{\mathbb{V}}(S) / \rho^2$, but the implementation does not use that quantity explicitly.

For the external-reference conditional test, WtCoxG also constructs the batch-effect score variance

$$
\widehat{\mathbb{V}}(S_{\mathrm{bat}}) = \left(\sum_{i=1}^n \tilde w_i^2\right) \cdot 2\hat f(1-\hat f) + \widehat{\mathbb{V}}(\hat f_{\mathrm{ext}}). \tag{23}
$$

Here $\tilde w_i = \frac{w_i}{2\sum_{j=1}^n w_j}$ in (7) are the normalized subject weights used in the batch-effect screen.

The covariance between the score and the batch-effect statistic is

$$
\widehat{\mathrm{Cov}}(S_{\mathrm{bat}}, S)
= \left(\sum_{i=1}^n \tilde w_i \tilde R_i\right) 2\hat f(1-\hat f) + 2b\left(\sum_{i=1}^n R_i\right) \widehat{\mathbb{V}}(\hat f_{\mathrm{ext}}). \tag{24}
$$

In the implementation, this covariance is re-scaled by the SPA-implied score variance before entering the bivariate normal tail calculation. Its role is to specify the joint Gaussian approximation of $(S, S_{\mathrm{bat}})$ after selection on the batch-effect screen. Without this covariance term, the rectangle probabilities used in the conditional p-value would ignore dependence between the screening statistic and the association score.

## 5. Saddlepoint Approximation

If $|Z|$ is below the SPA cutoff, WtCoxG uses the normal approximation

$$
p_{\mathrm{norm}} = 2\Phi(-|Z|). \tag{25}
$$

Otherwise it applies SPA to the score statistic under a diploid binomial model with allele frequency $\hat f$.

### 5.1 Cumulant Generating Function

For a single genotype term with $G \sim \mathrm{Binomial}(2, \hat f)$, define the moment generating function

$$
M_G(t) = \left[(1-\hat f) + \hat f e^t\right]^2. \tag{26}
$$

and the cumulant generating function

$$
K_G(t) = \log M_G(t). \tag{27}
$$

WtCoxG evaluates the score CGF as

$$
H(t) = \sum_{i=1}^n K_G\left(t\tilde R_i\right) - 2b\left(\sum_{i=1}^n R_i\right)\hat f\, t + 2b^2\left(\sum_{i=1}^n R_i\right)^2 \widehat{\mathbb{V}}(\hat f_{\mathrm{ext}}) t^2. \tag{28}
$$

Its first derivative is

$$
H'(t) = \sum_{i=1}^n \tilde R_i K_G'\left(t\tilde R_i\right) - 2b\left(\sum_{i=1}^n R_i\right)\hat f + 4b^2\left(\sum_{i=1}^n R_i\right)^2 \widehat{\mathbb{V}}(\hat f_{\mathrm{ext}}) t. \tag{29}
$$

and the second derivative is

$$
H''(t) = \sum_{i=1}^n \tilde R_i^2 K_G''\left(t\tilde R_i\right) + 2n_{\mathrm{ext}} \left(\frac{\sum_{i=1}^n R_i}{n + n_{\mathrm{ext}}}\right)^2 \hat f(1-\hat f). \tag{30}
$$

The saddlepoint $\zeta$ solves

$$
H'(\zeta) = s. \tag{31}
$$

where $s = |S| / \rho$ or $s = -|S| / \rho$. The implementation brackets the root by interval expansion and then applies Brent's method.

### 5.2 Lugannani-Rice Tail Probability

Given the saddlepoint $\zeta$, WtCoxG computes

$$
w = \mathrm{sign}(\zeta)\sqrt{2\{\zeta s - H(\zeta)\}}, \qquad v = \zeta \sqrt{H''(\zeta)}. \tag{32}
$$

The one-sided tail probability is approximated by

$$
\Pr(S / \rho \le s) \simeq \Phi\left(w + \frac{1}{w}\log\frac{v}{w}\right). \tag{33}
$$

The two-sided SPA p-value is obtained by summing the upper and lower tails:

$$
p_{\mathrm{SPA}} = \Pr(S / \rho \ge |S| / \rho) + \Pr(S / \rho \le -|S| / \rho). \tag{34}
$$

The marker-level routine returns $p_{\mathrm{SPA}}$ when SPA is used, together with the unadjusted score $S$ and z-score $Z$.

## 6. Conditional External-Reference P-value

The external-reference branch does not report the raw SPA p-value directly. Instead it combines the score statistic with the batch-effect screen.

Let $p_{\mathrm{bat}}$ be the batch-effect p-value and let $p_{\mathrm{cut}}$ be the screening threshold. Markers with $p_{\mathrm{bat}} < p_{\mathrm{cut}}$ are excluded.

For markers that pass the screen, WtCoxG computes two bivariate normal tail probabilities:

- $p_0$ under $\sigma^2 = 0$;
- $p_1$ under the estimated $\sigma^2 > 0$.

Both use a rectangle probability of the form

$$
\Pr\left(S / \rho \le -|S| / \rho,\; \ell \le S_{\mathrm{bat}} \le u\right). \tag{35}
$$

with bounds

$$
\ell = -z_{1-p_{\mathrm{cut}}/2}\sqrt{\widehat{\mathbb{V}}(S_{\mathrm{bat}})}\sqrt{\rho_{w0}}, \qquad u = z_{1-p_{\mathrm{cut}}/2}\sqrt{\widehat{\mathbb{V}}(S_{\mathrm{bat}})}\sqrt{\rho_{w0}}. \tag{36}
$$

In the current implementation, $\rho_{w1} = \rho_{w0}$ for the external-reference path.

The corresponding covariance matrices are

$$
\Sigma_0 = \begin{pmatrix}
\widehat{\mathbb{V}}_0(S) & \widehat{\mathrm{Cov}}_0(S_{\mathrm{bat}}, S) \\
\widehat{\mathrm{Cov}}_0(S_{\mathrm{bat}}, S) & \widehat{\mathbb{V}}(S_{\mathrm{bat}})
\end{pmatrix},
\qquad
\Sigma_1 = \begin{pmatrix}
\widehat{\mathbb{V}}_1(S) & \widehat{\mathrm{Cov}}_1(S_{\mathrm{bat}}, S) \\
\widehat{\mathrm{Cov}}_1(S_{\mathrm{bat}}, S) & \widehat{\mathbb{V}}(S_{\mathrm{bat}}) + \sigma^2
\end{pmatrix}. \tag{37}
$$

The final conditional p-value function introduced in (12), now evaluated at the observed marker statistics, is

$$
p_{\mathrm{con}} = \frac{2\{\tau p_1 + (1-\tau)p_0\}}{p_{\mathrm{deno}}}. \tag{38}
$$

where

$$
p_{\mathrm{deno}} = \tau\left[\Phi_{\sigma^2}(u / \sqrt{\rho_{w1}}) - \Phi_{\sigma^2}(\ell / \sqrt{\rho_{w1}})\right] + (1-\tau)(1-p_{\mathrm{cut}}). \tag{39}
$$

and $\Phi_{\sigma^2}(\cdot)$ denotes the normal CDF with variance $\widehat{\mathbb{V}}(S_{\mathrm{bat}}) + \sigma^2$ after the corresponding variance-ratio scaling.

The internal-only branch skips this conditional adjustment and reports the direct SPA / normal p-value from Sections 3, 4, and 5.
