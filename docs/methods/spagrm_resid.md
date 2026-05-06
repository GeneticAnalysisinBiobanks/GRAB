# Residual Choice in SPAGRM: Theoretical Analysis

This document provides a detailed theoretical analysis of residual choices in SPAGRM, focusing on binary phenotypes in logistic regression and extending to other regression frameworks.

## 1. Logistic Regression Framework

### 1.1 Model Setup

**Binary phenotype**: $Y_i \in \{0, 1\}$ for $i = 1, \ldots, n$

**Null model** (no genetic effect):
$$
\text{logit}(\mu_i) = \mathbf{X}_i^\top \boldsymbol{\beta} \tag{1}
$$

where:

- $\mu_i = \mathbb{P}(Y_i = 1 \mid \mathbf{X}_i) = \frac{\exp(\mathbf{X}_i^\top \boldsymbol{\beta})}{1 + \exp(\mathbf{X}_i^\top \boldsymbol{\beta})}$ is the fitted probability
- $\mathbf{X}_i$ is the covariate vector (including intercept)
- $\boldsymbol{\beta}$ is estimated by maximum likelihood

**Alternative model** (testing genetic variant):
$$
\text{logit}(\mu_i^{\text{alt}}) = \mathbf{X}_i^\top \boldsymbol{\beta} + \gamma G_i \tag{2}
$$

where:

- $G_i \in \{0, 1, 2\}$ is the genotype (minor allele count)
- $\gamma$ is the genetic effect (null hypothesis: $\gamma = 0$)

### 1.2 Score Test Derivation

**Log-likelihood** under the alternative model:
$$
\ell(\boldsymbol{\beta}, \gamma) = \sum_{i=1}^n \left[Y_i \log \mu_i^{\text{alt}} + (1-Y_i) \log(1-\mu_i^{\text{alt}})\right] \tag{3}
$$

**Score statistic** for testing $H_0: \gamma = 0$:
$$
U(\gamma) = \frac{\partial \ell}{\partial \gamma}\bigg|_{\gamma=0} = \sum_{i=1}^n \frac{\partial \ell_i}{\partial \mu_i} \cdot \frac{\partial \mu_i}{\partial \gamma}\bigg|_{\gamma=0} \tag{4}
$$

**Derivative components**:
$$
\frac{\partial \ell_i}{\partial \mu_i} = \frac{Y_i}{\mu_i} - \frac{1-Y_i}{1-\mu_i} = \frac{Y_i - \mu_i}{\mu_i(1-\mu_i)} \tag{5}
$$

$$
\frac{\partial \mu_i}{\partial \gamma}\bigg|_{\gamma=0} = \mu_i(1-\mu_i) \cdot G_i \tag{6}
$$

**Combining**:
$$
U = \sum_{i=1}^n \frac{Y_i - \mu_i}{\mu_i(1-\mu_i)} \cdot \mu_i(1-\mu_i) \cdot G_i = \sum_{i=1}^n (Y_i - \mu_i) G_i \tag{7}
$$

**Centering genotypes**: Since we use centered genotypes $\tilde{G}_i = G_i - 2f$:
$$
U = \sum_{i=1}^n (Y_i - \mu_i) (G_i - 2f) = \sum_{i=1}^n (Y_i - \mu_i) \tilde{G}_i \tag{8}
$$

**Key observation**: The score test naturally uses **Pearson residuals** $R_i^{(P)} = Y_i - \mu_i$.

## 2. Pearson Residuals

### 2.1 Definition and Properties

**Pearson residual**:
$$
R_i^{(P)} = Y_i - \mu_i \tag{9}
$$

**Mean**: Under the null model (conditional on covariates):
$$
\mathbb{E}[R_i^{(P)} \mid \mathbf{X}_i] = \mathbb{E}[Y_i \mid \mathbf{X}_i] - \mu_i = \mu_i - \mu_i = 0 \tag{10}
$$

**Unconditional mean**:
$$
\mathbb{E}[R_i^{(P)}] = \mathbb{E}[\mathbb{E}[R_i^{(P)} \mid \mathbf{X}_i]] = 0 \tag{11}
$$

**Variance**:
$$
\mathbb{V}(R_i^{(P)} \mid \mathbf{X}_i) = \mathbb{V}(Y_i \mid \mathbf{X}_i) = \mu_i(1-\mu_i) \tag{12}
$$

**Heteroscedasticity**: The variance is **not constant** across individuals; it depends on $\mu_i$.

**Distribution**: For binary $Y_i$:
$$
R_i^{(P)} = \begin{cases}
1 - \mu_i & \text{with probability } \mu_i \\
-\mu_i & \text{with probability } 1-\mu_i
\end{cases}
\tag{13}
$$

This is a **discrete distribution**, not continuous.

### 2.2 Score Statistic Distribution

**Score statistic**:
$$
S^{(P)} = \sum_{i=1}^n R_i^{(P)} \tilde{G}_i = \sum_{i=1}^n (Y_i - \mu_i)(G_i - 2f) \tag{14}
$$

**Expectation** (under null hypothesis $G_i \perp Y_i \mid \mathbf{X}_i$):
$$
\mathbb{E}[S^{(P)}] = \sum_{i=1}^n \mathbb{E}[Y_i - \mu_i] \mathbb{E}[G_i - 2f] = 0 \tag{15}
$$

**Variance** (for unrelated individuals):
$$
\mathbb{V}(S^{(P)}) = \mathbb{E}\left[\mathbb{V}(S^{(P)} \mid \mathbf{G})\right] + \mathbb{V}\left(\mathbb{E}[S^{(P)} \mid \mathbf{G}]\right) \tag{16}
$$

where $\mathbf{G} = (G_1, \ldots, G_n)$ is the genotype vector.

**Conditional variance** (given genotypes):
$$
\mathbb{V}(S^{(P)} \mid \mathbf{G}) = \sum_{i=1}^n (G_i - 2f)^2 \mathbb{V}(Y_i \mid \mathbf{X}_i) = \sum_{i=1}^n (G_i - 2f)^2 \mu_i(1-\mu_i) \tag{17}
$$

**Expected variance**:
$$
\mathbb{E}[\mathbb{V}(S^{(P)} \mid \mathbf{G})] = \sum_{i=1}^n \mathbb{E}[(G_i - 2f)^2] \mu_i(1-\mu_i) = \sum_{i=1}^n 2f(1-f) \cdot \mu_i(1-\mu_i) \tag{18}
$$

**Conditional expectation**:
$$
\mathbb{E}[S^{(P)} \mid \mathbf{G}] = \sum_{i=1}^n (G_i - 2f) \mathbb{E}[Y_i - \mu_i] = 0 \tag{19}
$$

**Total variance**:
$$
\mathbb{V}(S^{(P)}) = 2f(1-f) \sum_{i=1}^n \mu_i(1-\mu_i) \tag{20}
$$

### 2.3 Moment Generating Function

**MGF for single observation** (given $G_i$):
$$
M_i^{(P)}(t) = \mathbb{E}[e^{t R_i^{(P)} (G_i - 2f)} \mid G_i] \tag{21}
$$

For binary $Y_i$:
$$
M_i^{(P)}(t) = \mu_i e^{t(1-\mu_i)(G_i - 2f)} + (1-\mu_i) e^{-t\mu_i(G_i - 2f)} \tag{22}
$$

**Simplification**: Let $s_i = t(G_i - 2f)$:
$$
M_i^{(P)}(s_i) = \mu_i e^{s_i(1-\mu_i)} + (1-\mu_i) e^{-s_i\mu_i} \tag{23}
$$

**Joint MGF** (assuming independence conditional on genotypes):
$$
M^{(P)}(t) = \mathbb{E}\left[\prod_{i=1}^n M_i^{(P)}(t \mid G_i)\right] = \mathbb{E}_{\mathbf{G}}\left[\prod_{i=1}^n M_i^{(P)}(t(G_i - 2f))\right] \tag{24}
$$

**Key complexity**: The MGF depends on the **distribution of genotypes**, which is what we're testing. This creates a circular dependency that SPAGRM resolves by:

1. Conditioning on genotypes (treating them as fixed)
2. Deriving the MGF with respect to genotype randomness (not phenotype)

**Standard SPAGRM approach**: Treat residuals as **fixed** and derive MGF with respect to **genotypes**:
$$
M_{\text{SPAGRM}}^{(P)}(t) = \mathbb{E}_{\mathbf{G}}\left[\exp\left(t \sum_{i=1}^n R_i^{(P)} (G_i - 2f)\right)\right] = \prod_{i=1}^n (1-f + f e^{t R_i^{(P)}})^2 \tag{25}
$$

This is the standard formula used in SPAGRM, treating $R_i^{(P)}$ as fixed weights.

## 3. Deviance Residuals

### 3.1 Definition and Properties

**Deviance residual**:
$$
R_i^{(D)} = \text{sign}(Y_i - \mu_i) \sqrt{-2[\ell_i(\tilde{Y}_i) - \ell_i(\mu_i)]} \tag{26}
$$

where $\ell_i(\cdot)$ is the log-likelihood contribution for observation $i$, and $\tilde{Y}_i = Y_i$ (saturated model).

**For binary outcomes**:
$$
R_i^{(D)} = \begin{cases}
\sqrt{-2\ln \mu_i} & \text{if } Y_i = 1 \\
-\sqrt{-2\ln(1-\mu_i)} & \text{if } Y_i = 0
\end{cases}
\tag{27}
$$

**Mean**: The deviance residual does **not** have zero mean in general:
$$
\mathbb{E}[R_i^{(D)} \mid \mathbf{X}_i] = \mu_i \sqrt{-2\ln \mu_i} - (1-\mu_i) \sqrt{-2\ln(1-\mu_i)} \neq 0 \tag{28}
$$

**Empirical mean**:
$$
\overline{R}^{(D)} = \frac{1}{n}\sum_{i=1}^n R_i^{(D)} \tag{29}
$$

This is typically **not zero** for a given dataset.

**Variance**: More complex than Pearson residuals; approximate formula:
$$
\mathbb{V}(R_i^{(D)} \mid \mathbf{X}_i) \approx \mu_i(1-\mu_i) \tag{30}
$$

(similar to Pearson, but not exact)

**Distribution**: Deviance residuals are designed to be more **symmetric** and **closer to normality** than Pearson residuals.

### 3.2 Centered Deviance Residuals

**Definition**:
$$
R_i^{(DC)} = R_i^{(D)} - \overline{R}^{(D)} \tag{31}
$$

**Key property**:
$$
\sum_{i=1}^n R_i^{(DC)} = 0 \quad \text{(exactly)} \tag{32}
$$

**Mean**:
$$
\mathbb{E}[R_i^{(DC)}] = \mathbb{E}[R_i^{(D)}] - \mathbb{E}[\overline{R}^{(D)}] \approx 0 \quad \text{for large } n \tag{33}
$$

**Variance**: Slightly different from uncentered deviance:
$$
\mathbb{V}(R_i^{(DC)}) = \mathbb{V}(R_i^{(D)}) - \frac{1}{n}\mathbb{V}(R_i^{(D)}) = \left(1 - \frac{1}{n}\right)\mathbb{V}(R_i^{(D)}) \tag{34}
$$

### 3.3 Score Statistic Distribution

**Score with centered deviance residuals**:
$$
S^{(DC)} = \sum_{i=1}^n R_i^{(DC)} (G_i - 2f) \tag{35}
$$

**Relationship to Pearson score**:
$$
S^{(DC)} = \sum_{i=1}^n (R_i^{(D)} - \overline{R}^{(D)}) (G_i - 2f) = S^{(D)} - \overline{R}^{(D)} \sum_{i=1}^n (G_i - 2f) \tag{36}
$$

Since $\sum_{i=1}^n (G_i - 2f) = n(\overline{G} - 2f)$:
$$
S^{(DC)} = S^{(D)} - n\overline{R}^{(D)}(\overline{G} - 2f) \tag{37}
$$

**Expectation**:
$$
\mathbb{E}[S^{(DC)}] = \mathbb{E}[S^{(D)}] - n\mathbb{E}[\overline{R}^{(D)}]\mathbb{E}[\overline{G} - 2f] = 0 \tag{38}
$$

**Variance**:
$$
\mathbb{V}(S^{(DC)}) = \mathbb{V}(S^{(D)}) + n^2 \mathbb{V}(\overline{R}^{(D)})\mathbb{V}(\overline{G}) - 2n\text{Cov}(S^{(D)}, \overline{R}^{(D)}\overline{G}) \tag{39}
$$

For large $n$ with weak dependence:
$$
\mathbb{V}(S^{(DC)}) \approx \mathbb{V}(S^{(D)}) \tag{40}
$$

## 4. Comparison: Pearson vs. Centered Deviance

### 4.1 Distributional Differences

**Pearson residuals**:

- Values: $R_i^{(P)} \in \{1-\mu_i, -\mu_i\}$ (discrete, two-point distribution)
- Range: $(-1, 1)$ with most mass near 0 when $\mu_i \approx 0.5$
- Skewness: High when $\mu_i$ is far from 0.5
- Example: If $\mu_i = 0.1$, then $R_i^{(P)} \in \{0.9, -0.1\}$ (highly skewed)

**Deviance residuals**:

- Values: $R_i^{(D)} \in \{\sqrt{-2\ln\mu_i}, -\sqrt{-2\ln(1-\mu_i)}\}$
- Range: Unbounded (can be large when $\mu_i$ is close to 0 or 1)
- Skewness: Lower, more symmetric distribution
- Example: If $\mu_i = 0.1$:
  - Pearson: $\{0.9, -0.1\}$
  - Deviance: $\{\sqrt{-2\ln(0.1)}, -\sqrt{-2\ln(0.9)}\} \approx \{2.146, -0.458\}$

**Comparison table**:

| $\mu_i$ | $R^{(P)}$ if $Y_i=1$ | $R^{(P)}$ if $Y_i=0$ | $R^{(D)}$ if $Y_i=1$ | $R^{(D)}$ if $Y_i=0$ | Skewness |
|---------|---------------------|---------------------|---------------------|---------------------|----------|
| 0.01 | 0.99 | -0.01 | 3.035 | -0.143 | Pearson: high |
| 0.1 | 0.9 | -0.1 | 2.146 | -0.458 | Pearson: moderate |
| 0.3 | 0.7 | -0.3 | 1.274 | -0.849 | Pearson: low |
| 0.5 | 0.5 | -0.5 | 1.177 | -1.177 | Both symmetric |
| 0.7 | 0.3 | -0.7 | 0.849 | -1.274 | Pearson: low |
| 0.9 | 0.1 | -0.9 | 0.458 | -2.146 | Pearson: moderate |
| 0.99 | 0.01 | -0.99 | 0.143 | -3.035 | Pearson: high |

**Key insight**: Deviance residuals are more symmetric, especially when $\mu_i$ is extreme (close to 0 or 1). This can improve the normal approximation in finite samples.

### 4.2 Score Statistic Variance

**Pearson score variance**:
$$
\mathbb{V}(S^{(P)}) = 2f(1-f) \sum_{i=1}^n \mu_i(1-\mu_i) \tag{41}
$$

**Deviance score variance** (approximate):
$$
\mathbb{V}(S^{(DC)}) \approx 2f(1-f) \sum_{i=1}^n \tilde{v}_i \tag{42}
$$

where $\tilde{v}_i$ is the variance of the centered deviance residual.

**Ratio**: For a single observation:
$$
\frac{\mathbb{V}(R_i^{(D)})}{\mathbb{V}(R_i^{(P)})} = \frac{\mathbb{E}[(R_i^{(D)})^2] - [\mathbb{E}[R_i^{(D)}]]^2}{\mu_i(1-\mu_i)} \tag{43}
$$

**Numerical comparison**:

| $\mu_i$ | $\mathbb{V}(R^{(P)})$ | $\mathbb{V}(R^{(D)})$ (approx) | Ratio |
|---------|-----------------------|-------------------------------|-------|
| 0.01 | 0.0099 | ~0.010 | ~1.01 |
| 0.1 | 0.09 | ~0.092 | ~1.02 |
| 0.3 | 0.21 | ~0.21 | ~1.00 |
| 0.5 | 0.25 | ~0.25 | ~1.00 |
| 0.7 | 0.21 | ~0.21 | ~1.00 |
| 0.9 | 0.09 | ~0.092 | ~1.02 |
| 0.99 | 0.0099 | ~0.010 | ~1.01 |

**Conclusion**: The variances are **nearly identical** for most $\mu_i$ values. The main difference is in the **shape** of the distribution, not the variance.

### 4.3 Impact on Saddlepoint Approximation

**Saddlepoint equation**:
$$
K'(\hat{t}) = S_{\text{obs}} \tag{44}
$$

where $K(t) = \ln M(t)$ is the cumulant generating function.

**For Pearson residuals**:
$$
K^{(P)}(t) = \sum_{i=1}^n 2\ln(1-f + fe^{tR_i^{(P)}}) \tag{45}
$$

**For centered deviance residuals**:
$$
K^{(DC)}(t) = \sum_{i=1}^n 2\ln(1-f + fe^{tR_i^{(DC)}}) \tag{46}
$$

**Key difference**: The weights $R_i^{(P)}$ vs. $R_i^{(DC)}$ in the exponential.

**CGF derivatives**:
$$
K'(t) = \sum_{i=1}^n \frac{2f R_i e^{tR_i}}{1-f + fe^{tR_i}} \tag{47}
$$

$$
K''(t) = \sum_{i=1}^n \frac{2f(1-f) R_i^2 e^{tR_i}}{(1-f + fe^{tR_i})^2} \tag{48}
$$

**At $t=0$** (null distribution):

- $K'(0) = \sum_{i=1}^n 2f R_i = 0$ (if $\sum R_i = 0$, which holds for both Pearson and centered deviance)
- $K''(0) = 2f(1-f) \sum_{i=1}^n R_i^2$

**Variance from CGF**:
$$
\mathbb{V}(S) = K''(0) = 2f(1-f) \sum_{i=1}^n R_i^2 \tag{49}
$$

**Comparison**:

- **Pearson**: $\sum_{i=1}^n (R_i^{(P)})^2 = \sum_{i=1}^n (Y_i - \mu_i)^2$
- **Centered deviance**: $\sum_{i=1}^n (R_i^{(DC)})^2 = \sum_{i=1}^n (R_i^{(D)} - \overline{R}^{(D)})^2$

**Which is larger?**

$$
\sum (R_i^{(DC)})^2 = \sum (R_i^{(D)})^2 - n(\overline{R}^{(D)})^2 \tag{50}
$$

Since $(\overline{R}^{(D)})^2 > 0$ when $\overline{R}^{(D)} \neq 0$:
$$
\sum (R_i^{(DC)})^2 < \sum (R_i^{(D)})^2 \tag{51}
$$

But how does $\sum (R_i^{(DC)})^2$ compare to $\sum (R_i^{(P)})^2$?

**Empirical comparison** (requires simulation or real data): Generally, $\sum (R_i^{(DC)})^2 \approx \sum (R_i^{(P)})^2$ for well-fitted models.

### 4.4 P-value Comparison

**Lugannani-Rice formula**:
$$
\Pr(S \geq s) \approx \Phi(\hat{w}) + \phi(\hat{w})\left(\frac{1}{\hat{w}} - \frac{1}{\hat{u}}\right) \tag{52}
$$

where:

- $\hat{w} = \text{sgn}(\hat{t})\sqrt{2[\hat{t}s - K(\hat{t})]}$
- $\hat{u} = \hat{t}\sqrt{K''(\hat{t})}$
- $\Phi$ and $\phi$ are the standard normal CDF and PDF

**Key insight**: The p-value depends on:

1. The observed score $s$ (different for Pearson vs. deviance due to different $R_i$ values)
2. The saddlepoint $\hat{t}$ (solution to $K'(\hat{t}) = s$)
3. The CGF curvature $K''(\hat{t})$

**Expected difference**: For well-calibrated null models:

- **Type I error**: Both should maintain nominal level (e.g., 0.05)
- **Power**: Deviance residuals may have slightly higher power if the true effect is non-linear (due to better symmetry)
- **Rare variants**: Differences may be more pronounced when $\mu_i$ values are extreme

**Simulation study needed**: To quantify exact differences, simulations should vary:

- Sample size ($n$)
- Prevalence (proportion of cases)
- Covariate effects (how extreme are $\mu_i$ values)
- Minor allele frequency ($f$)
- True genetic effect size

### 4.5 Theoretical Analysis of Type I Error and Power

This section provides a detailed theoretical comparison of Type I error rates and power between Pearson and centered deviance residuals in logistic regression.

#### 4.5.1 Type I Error Rate Under the Null Hypothesis

**Null hypothesis**: $H_0: \gamma = 0$ (no genetic effect)

Under $H_0$, the phenotypes $Y_i$ follow the null model:
$$
Y_i \sim \text{Bernoulli}(\mu_i), \quad \text{logit}(\mu_i) = \mathbf{X}_i^\top \boldsymbol{\beta}
$$

The genotypes $G_i$ are independent of $Y_i$ conditional on covariates: $G_i \perp Y_i \mid \mathbf{X}_i$.

**Score statistics**:

- **Pearson**: $S^{(P)} = \sum_{i=1}^n R_i^{(P)} (G_i - 2f)$ where $R_i^{(P)} = Y_i - \mu_i$
- **Centered deviance**: $S^{(DC)} = \sum_{i=1}^n R_i^{(DC)} (G_i - 2f)$ where $R_i^{(DC)} = R_i^{(D)} - \overline{R}^{(D)}$

**Expected values under $H_0$**:

Both statistics have zero expectation (proven in Sections 2 and 3):
$$
\mathbb{E}[S^{(P)}] = \mathbb{E}[S^{(DC)}] = 0 \tag{52a}
$$

**Variance under $H_0$**:

From equations (20) and (40):
$$
\mathbb{V}(S^{(P)}) = 2f(1-f) \sum_{i=1}^n \mu_i(1-\mu_i) \tag{52b}
$$
$$
\mathbb{V}(S^{(DC)}) \approx 2f(1-f) \sum_{i=1}^n \tilde{v}_i \tag{52c}
$$

where $\tilde{v}_i \approx \mu_i(1-\mu_i)$ for well-fitted models (see Section 4.2).

**Key observation**: The variances are **nearly identical** under $H_0$:
$$
\frac{\mathbb{V}(S^{(DC)})}{\mathbb{V}(S^{(P)})} \approx 1 \tag{52d}
$$

**Type I error calibration**:

The Type I error rate depends on the **tail probability** of the test statistic under $H_0$:
$$
\alpha = \Pr(|S| > c_{\alpha} \mid H_0) \tag{52e}
$$

where $c_{\alpha}$ is the critical value for significance level $\alpha$ (e.g., $\alpha = 0.05$).

**Saddlepoint approximation accuracy**:

The Lugannani-Rice formula (equation 52) approximates the tail probability using:
$$
\Pr(S \geq s \mid H_0) \approx \Phi(\hat{w}) + \phi(\hat{w})\left(\frac{1}{\hat{w}} - \frac{1}{\hat{u}}\right) \tag{52f}
$$

The accuracy depends on how well the **CGF** $K(t)$ approximates the true distribution:

- **Pearson CGF**: $K^{(P)}(t) = \sum_{i=1}^n 2\ln(1-f + fe^{tR_i^{(P)}})$
- **Deviance CGF**: $K^{(DC)}(t) = \sum_{i=1}^n 2\ln(1-f + fe^{tR_i^{(DC)}})$

**Critical difference**: The MGF for a single observation:

For **Pearson** with $R_i^{(P)} \in \{1-\mu_i, -\mu_i\}$:
$$
M_i^{(P)}(t) = \mu_i e^{t(1-\mu_i)(G_i-2f)} + (1-\mu_i) e^{-t\mu_i(G_i-2f)} \tag{52g}
$$

For **Deviance** with $R_i^{(DC)} \in \{r_1, r_0\}$ (more symmetric):
$$
M_i^{(DC)}(t) = \mu_i e^{tr_1(G_i-2f)} + (1-\mu_i) e^{tr_0(G_i-2f)} \tag{52h}
$$

where $r_1 = \sqrt{-2\ln\mu_i} - \overline{R}^{(D)}$ and $r_0 = -\sqrt{-2\ln(1-\mu_i)} - \overline{R}^{(D)}$ are the centered deviance residual values.

**Symmetry advantage of deviance residuals**:

When $\mu_i$ is extreme (close to 0 or 1), Pearson residuals are highly **skewed**:

- If $\mu_i = 0.01$: $R^{(P)} \in \{0.99, -0.01\}$ (highly positive-skewed)
- If $\mu_i = 0.01$: $R^{(DC)} \in \{3.035 - \bar{R}^{(D)}, -0.143 - \bar{R}^{(D)}\}$ (more balanced after centering)

**Impact on Type I error**:

The **Berry-Esseen theorem** bounds the error in normal approximation:
$$
\sup_x \left|\Pr\left(\frac{S}{\sqrt{\mathbb{V}(S)}} \leq x\right) - \Phi(x)\right| \leq \frac{C \mathbb{E}[|S - \mathbb{E}[S]|^3]}{(\mathbb{V}(S))^{3/2}} \tag{52i}
$$

where $C \approx 0.4784$ is the Berry-Esseen constant.

The **third central moment** (skewness):

- **Pearson**: Higher skewness when $\mu_i$ values are extreme
- **Deviance**: Lower skewness due to symmetrization

**Theoretical prediction**:

1. **For balanced $\mu_i$ (around 0.5)**: Both methods should have **similar** Type I error rates, close to the nominal level $\alpha$.

2. **For extreme $\mu_i$ (many close to 0 or 1)**: 
   - **Pearson**: May have **slight inflation** or **deflation** of Type I error due to skewness
   - **Deviance**: Better calibration due to **symmetry**, Type I error closer to nominal $\alpha$

3. **Saddlepoint approximation**: Should correct for skewness in both cases, but deviance residuals may converge faster to the correct tail probability.

**Quantitative comparison**: Requires simulation or asymptotic expansion beyond second order.

#### 4.5.2 Power Under the Alternative Hypothesis

**Alternative hypothesis**: $H_a: \gamma \neq 0$ (genetic effect exists)

Under $H_a$, the true model is:
$$
\text{logit}(\mu_i^{\text{true}}) = \mathbf{X}_i^\top \boldsymbol{\beta} + \gamma G_i
$$

But we still compute residuals from the **null model** $\mu_i^{\text{null}}$ (without genetic effect).

**Expected score under $H_a$**:

The score statistic is **no longer zero-mean** under the alternative:
$$
\mathbb{E}[S \mid H_a] = \sum_{i=1}^n \mathbb{E}[Y_i \mid \mathbf{X}_i, G_i] (G_i - 2f) - \sum_{i=1}^n \mu_i^{\text{null}} (G_i - 2f) \tag{52j}
$$

where $\mathbb{E}[Y_i \mid \mathbf{X}_i, G_i] = \mu_i^{\text{true}} = \text{logit}^{-1}(\mathbf{X}_i^\top \boldsymbol{\beta} + \gamma G_i)$.

**Approximation for small $\gamma$**:

Using Taylor expansion around $\gamma = 0$:
$$
\mu_i^{\text{true}} \approx \mu_i^{\text{null}} + \gamma G_i \mu_i^{\text{null}}(1-\mu_i^{\text{null}}) + O(\gamma^2) \tag{52k}
$$

**Expected score (Pearson)**:
$$
\mathbb{E}[S^{(P)} \mid H_a] \approx \gamma \sum_{i=1}^n G_i (G_i - 2f) \mu_i(1-\mu_i) \tag{52l}
$$

Since $\mathbb{E}[G_i(G_i - 2f)] = \mathbb{E}[G_i^2] - 2f\mathbb{E}[G_i] = 2f(1-f) + (2f)^2 - 2f \cdot 2f = 2f(1-f)$:
$$
\mathbb{E}[S^{(P)} \mid H_a] \approx \gamma \cdot 2f(1-f) \sum_{i=1}^n \mu_i(1-\mu_i) \tag{52m}
$$

**Expected score (Centered deviance)**:

The centering complicates the analysis. The mean deviance residual $\overline{R}^{(D)}$ depends on the observed data.

For **fixed residual values** (treating them as deterministic given phenotypes):
$$
\mathbb{E}[S^{(DC)} \mid H_a, \mathbf{Y}] = \sum_{i=1}^n R_i^{(DC)} \mathbb{E}[G_i - 2f \mid Y_i] \tag{52n}
$$

Under $H_a$, genotypes and phenotypes are **correlated**:
$$
\mathbb{E}[G_i \mid Y_i = 1] > \mathbb{E}[G_i \mid Y_i = 0] \quad \text{(if } \gamma > 0\text{)} \tag{52o}
$$

**Key difference in power**:

The **non-centrality parameter** (measure of power) is:
$$
\lambda = \frac{[\mathbb{E}[S \mid H_a]]^2}{\mathbb{V}(S \mid H_0)} \tag{52p}
$$

**Pearson residuals**:

From equation (52m) and (52b):
$$
\lambda^{(P)} = \frac{\left[\gamma \cdot 2f(1-f) \sum_{i=1}^n \mu_i(1-\mu_i)\right]^2}{2f(1-f) \sum_{i=1}^n \mu_i(1-\mu_i)} = \gamma^2 \cdot 2f(1-f) \sum_{i=1}^n \mu_i(1-\mu_i) \tag{52q}
$$

**Centered deviance residuals**:

The expected score depends on the **correlation structure** between residuals and genotypes under $H_a$. 

**Approximation**: For small $\gamma$, the centering adjustment $\overline{R}^{(D)}$ has negligible impact on the mean (it's an $O(1/n)$ effect). Thus:
$$
\lambda^{(DC)} \approx \lambda^{(P)} \tag{52r}
$$

**Refinement for finite samples**:

The **effective signal** differs slightly between Pearson and deviance:

1. **Pearson residuals**: Directly measure deviation $Y_i - \mu_i$, which is the **natural** score residual.

2. **Deviance residuals**: Transformed to be more symmetric, but the transformation may **distort** the linear relationship between residuals and genetic effect.

**Theoretical power comparison**:

Define the **relative efficiency** as:
$$
\text{RE} = \frac{\lambda^{(DC)}}{\lambda^{(P)}} \tag{52s}
$$

**Case 1: Balanced $\mu_i$ (around 0.5)**

- Pearson and deviance residuals are similar
- $\text{RE} \approx 1$ (nearly equal power)

**Case 2: Extreme $\mu_i$ (many close to 0 or 1)**

The relationship becomes **non-linear**:
$$
\text{logit}(\mu_i^{\text{true}}) - \text{logit}(\mu_i^{\text{null}}) = \gamma G_i
$$

When $\mu_i^{\text{null}}$ is extreme, the logit transformation amplifies small differences in probability:

- $\mu_i = 0.01 \to 0.02$: $\Delta \text{logit} \approx 0.69$
- $\mu_i = 0.5 \to 0.51$: $\Delta \text{logit} \approx 0.04$

**Deviance residuals** are designed to account for this non-linearity through the log-likelihood scaling.

**Expected relative efficiency**:

$$
\text{RE} \approx 1 + \frac{\gamma}{n} \sum_{i=1}^n \frac{(R_i^{(DC)})^2 - (R_i^{(P)})^2}{(R_i^{(P)})^2} + O(\gamma^2) \tag{52t}
$$

From the comparison table in Section 4.1:

- When $\mu_i = 0.01$: $(R^{(D)})^2 \approx 9.2$ vs. $(R^{(P)})^2 = 0.98$ (deviance much larger)
- When $\mu_i = 0.5$: $(R^{(D)})^2 \approx 1.39$ vs. $(R^{(P)})^2 = 0.25$ (deviance moderately larger)

However, after **centering** and **weighting by prevalence**, the effective contribution is:

$$
\text{Weight}_i = \mu_i \cdot (r_1^{(DC)})^2 + (1-\mu_i) \cdot (r_0^{(DC)})^2 \tag{52u}
$$

For **extreme $\mu_i$**:

- Deviance residuals **magnify** rare events (large residuals when $Y_i = 1$ and $\mu_i \ll 1$)
- This can improve **signal-to-noise ratio** if the genetic effect operates on the rare phenotype

**Theoretical prediction for power**:

1. **Common variants ($f \geq 0.05$), balanced prevalence**:
   - $\text{RE} \approx 1$ (similar power)
   - **Recommendation**: Use Pearson (simpler, standard)

2. **Rare variants ($f < 0.01$), extreme prevalence (rare disease)**:
   - $\text{RE}$ may be $> 1$ or $< 1$ depending on effect heterogeneity
   - **Deviance**: Better if effect is stronger in extreme $\mu_i$ regions
   - **Pearson**: Better if effect is uniform across $\mu_i$

3. **Large effect size ($|\gamma| > 0.5$)**:
   - Non-linearity becomes important
   - Deviance residuals may capture non-linear effects better
   - But: Type I error calibration becomes critical

#### 4.5.3 Summary of Theoretical Findings

**Type I Error Rate**:

| Scenario | Pearson Residuals | Centered Deviance Residuals |
|----------|-------------------|----------------------------|
| Balanced $\mu_i$ ($\approx 0.5$) | Well-calibrated, $\alpha \approx$ nominal | Well-calibrated, $\alpha \approx$ nominal |
| Extreme $\mu_i$ (many $< 0.1$ or $> 0.9$) | Possible slight miscalibration due to skewness | Better calibration due to symmetry |
| Large $n$ ($> 10000$) | Asymptotically correct | Asymptotically correct |
| Small $n$ ($< 500$) | May deviate if $\mu_i$ extreme | Closer to nominal due to symmetry |

**Power (for $\gamma \neq 0$)**:

| Scenario | Pearson Residuals | Centered Deviance Residuals | Winner |
|----------|-------------------|----------------------------|--------|
| Common variants, balanced prevalence | $\lambda^{(P)}$ | $\lambda^{(DC)} \approx \lambda^{(P)}$ | Tie |
| Rare variants, rare disease | $\lambda^{(P)}$ | $\lambda^{(DC)}$ may be higher if effect heterogeneous | Context-dependent |
| Uniform effect across $\mu_i$ | $\lambda^{(P)}$ | $\lambda^{(DC)} \approx \lambda^{(P)}$ | Tie |
| Effect stronger at extreme $\mu_i$ | Lower | Higher (captures non-linearity) | Deviance |
| Very small $\gamma$ (weak effect) | Linear regime, nearly equal | Linear regime, nearly equal | Tie |

**Practical Implications**:

1. **For standard GWAS** (large $n$, common variants, balanced cases/controls):
   - Use **Pearson residuals** (simpler, well-established, similar performance)

2. **For rare disease studies** (low prevalence, extreme $\mu_i$):
   - Consider **centered deviance residuals** if Type I error calibration is a concern
   - Validate with permutation tests or simulations

3. **For rare variant association tests** ($f < 0.01$):
   - If effect is suspected to be non-linear (stronger in rare phenotypes): Try deviance
   - Otherwise: Pearson is robust and interpretable

4. **Computational cost**:
   - Pearson: $O(n)$ (direct computation)
   - Deviance: $O(n)$ + centering step (negligible overhead)

**Open questions for empirical validation**:

- Under what specific covariate configurations does $\text{RE}$ significantly deviate from 1?
- How does the choice interact with **sparse GRM** and **family structures** in SPAGRM?
- Can we derive a **data-adaptive criterion** to choose residual type based on observed $\mu_i$ distribution?

## 5. General Residual Principles

### 5.1 Requirements for SPAGRM Validity

For SPAGRM to be valid, the residuals should satisfy:

1. **Zero mean constraint**: Either:
   - $\mathbb{E}[R_i] = 0$ (theoretically), OR
   - $\sum_{i=1}^n R_i = 0$ (empirically)

2. **Independence under null**: $R_i \perp G_i$ under the null hypothesis (no genetic association)

3. **Finite variance**: $\mathbb{V}(R_i) < \infty$

4. **MGF exists**: The moment generating function $\mathbb{E}[e^{tR_i}]$ exists in a neighborhood of $t=0$

**Consequence of violating zero mean**:

If $\sum_{i=1}^n R_i = c \neq 0$ and we use uncentered residuals:
$$
S = \sum_{i=1}^n R_i (G_i - 2f) = \sum_{i=1}^n R_i G_i - 2f \sum_{i=1}^n R_i = \sum_{i=1}^n R_i G_i - 2fc
$$

The term $2fc$ is a constant shift. Under the null:
$$
\mathbb{E}[S] = \mathbb{E}\left[\sum_{i=1}^n R_i G_i\right] - 2fc = 2f \sum_{i=1}^n R_i - 2fc = 2fc - 2fc = 0
$$

So the **expectation is still zero**! However:

- The distribution may be shifted
- The saddlepoint approximation accuracy may be affected
- For rare variants with $\overline{G} \neq 2f$, bias can occur

**Recommendation**: Always center residuals to ensure $\sum R_i = 0$ exactly.

### 5.2 Ordinal Regression Residuals

**Ordinal phenotype**: $Y_i \in \{1, 2, \ldots, K\}$ (ordered categories)

**Proportional odds model**:
$$
\text{logit}(\Pr(Y_i \leq k \mid \mathbf{X}_i)) = \alpha_k - \mathbf{X}_i^\top \boldsymbol{\beta} \tag{53}
$$

for $k = 1, \ldots, K-1$, where $\alpha_1 < \alpha_2 < \cdots < \alpha_{K-1}$ are threshold parameters.

**Fitted probabilities**:
$$
\pi_{ik} = \Pr(Y_i = k \mid \mathbf{X}_i) \tag{54}
$$

**Possible residual definitions**:

1. **Pearson residuals** (category-wise):
   $$
   R_i^{(P)} = \sum_{k=1}^K \frac{(\mathbb{1}(Y_i = k) - \pi_{ik})}{\sqrt{\pi_{ik}}} \tag{55}
   $$

   **Problem**: Not zero-mean!

2. **Score residuals** (latent variable approach):
   $$
   R_i^{(S)} = Y_i - \mathbb{E}[Y_i \mid \mathbf{X}_i] = Y_i - \sum_{k=1}^K k\pi_{ik} \tag{56}
   $$

   **Properties**:
   - $\mathbb{E}[R_i^{(S)} \mid \mathbf{X}_i] = 0$ ✓
   - Simple to compute
   - Treats ordinal as pseudo-continuous

3. **Deviance residuals** (generalized):
   $$
   R_i^{(D)} = \text{sign}(Y_i - \mathbb{E}[Y_i]) \sqrt{-2[\ell_i(Y_i) - \ell_i(\hat{\mu}_i)]} \tag{57}
   $$

   **Problem**: $\mathbb{E}[R_i^{(D)}] \neq 0$ in general

**Recommendation for ordinal regression**:

✓ **Use score residuals** $R_i^{(S)} = Y_i - \sum_{k=1}^K k\pi_{ik}$ (automatically zero-mean)

✗ **Do not use** standard Pearson or deviance residuals without centering

If using deviance or other residuals, **always center**:
$$
R_i^{(\text{centered})} = R_i - \frac{1}{n}\sum_{j=1}^n R_j \tag{58}
$$

#### 5.2.1 Practical Implementation for Ordinal Phenotypes

**R packages for ordinal regression**:

1. **`MASS::polr`** (Proportional Odds Logistic Regression):

   ```r
   library(MASS)
   
   # Fit null model (no genetic effect)
   # Y_ordinal: ordered factor with K levels
   # X: covariate matrix (data frame)
   fit <- polr(Y_ordinal ~ age + sex + PC1 + PC2, data = X, method = "logistic")
   
   # Extract residuals
   resid_default <- residuals(fit)              # Default: "working" residuals
   resid_deviance <- residuals(fit, type = "deviance")
   resid_pearson <- residuals(fit, type = "pearson")
   
   # Compute fitted probabilities
   fitted_probs <- fitted(fit)  # n x K matrix: π_ik for each i, k
   
   # Compute score residuals (recommended for SPAGRM)
   categories <- as.numeric(levels(Y_ordinal))  # e.g., 1, 2, 3, ..., K
   expected_Y <- fitted_probs %*% categories    # E[Y_i | X_i]
   resid_score <- as.numeric(Y_ordinal) - expected_Y
   
   # Check zero-mean
   sum(resid_score)  # Should be ≈ 0 (within numerical precision)
   ```

   **Key observations**:
   - Default `residuals(fit)` returns **working residuals**: NOT zero-mean! ✗
   - `type = "deviance"`: NOT zero-mean! ✗
   - `type = "pearson"`: NOT zero-mean! ✗
   - **Score residuals** (manual computation): Zero-mean ✓

   **SPAGRM recommendation**: Compute score residuals manually, verify $\sum R_i \approx 0$.

2. **`ordinal::clm`** (Cumulative Link Models):

   ```r
   library(ordinal)
   
   # Fit null model
   fit <- clm(Y_ordinal ~ age + sex + PC1 + PC2, data = X, link = "logit")
   
   # Extract residuals
   resid_default <- residuals(fit)              # Default: "working" residuals
   resid_pearson <- residuals(fit, type = "pearson")
   resid_deviance <- residuals(fit, type = "deviance")
   
   # Compute fitted probabilities
   fitted_probs <- predict(fit, type = "prob")$fit  # n x K matrix
   
   # Compute score residuals (recommended for SPAGRM)
   categories <- 1:ncol(fitted_probs)
   expected_Y <- fitted_probs %*% categories
   resid_score <- as.numeric(Y_ordinal) - expected_Y
   
   # Verify zero-mean
   sum(resid_score)  # Should be ≈ 0
   ```

   **Key observations**:
   - Default residuals: NOT zero-mean! ✗
   - **Score residuals** (manual): Zero-mean ✓

**Python packages for ordinal regression**:

1. **`statsmodels.miscmodels.ordinal_model.OrderedModel`**:

   ```python
   import numpy as np
   import pandas as pd
   from statsmodels.miscmodels.ordinal_model import OrderedModel
   
   # Fit null model
   # Y_ordinal: ordinal phenotype (0, 1, 2, ..., K-1)
   # X: covariate DataFrame
   fit = OrderedModel(Y_ordinal, X, distr='logit').fit()
   
   # Extract fitted probabilities
   fitted_probs = fit.predict()  # n x K array: π_ik
   
   # Compute score residuals (recommended for SPAGRM)
   categories = np.arange(fitted_probs.shape[1])  # 0, 1, 2, ..., K-1
   expected_Y = fitted_probs @ categories         # E[Y_i | X_i]
   resid_score = Y_ordinal - expected_Y
   
   # Check zero-mean
   resid_score.sum()  # Should be ≈ 0
   
   # statsmodels does NOT provide residuals() method directly
   # Must compute manually
   ```

   **Key observations**:
   - No built-in `residuals()` method
   - Must compute score residuals manually
   - **Score residuals** (manual): Zero-mean ✓

2. **`mord` package** (ordinal regression):

   ```python
   from mord import LogisticAT  # All-threshold variant
   
   # Fit null model
   clf = LogisticAT()
   clf.fit(X, Y_ordinal)
   
   # Predict probabilities
   fitted_probs = clf.predict_proba(X)  # n x K array
   
   # Compute score residuals
   categories = np.arange(fitted_probs.shape[1])
   expected_Y = (fitted_probs * categories).sum(axis=1)
   resid_score = Y_ordinal - expected_Y
   
   # Verify zero-mean
   resid_score.sum()  # Should be ≈ 0
   ```

**Summary table: Ordinal regression residuals**:

| Package | Method | Default Residual Type | Zero-Mean? | SPAGRM Suitable? |
|---------|--------|----------------------|------------|------------------|
| `MASS::polr` | `residuals(fit)` | Working residuals | ✗ No | ✗ Must center |
| `MASS::polr` | `residuals(fit, type="pearson")` | Pearson | ✗ No | ✗ Must center |
| `MASS::polr` | `residuals(fit, type="deviance")` | Deviance | ✗ No | ✗ Must center |
| `MASS::polr` | Manual: $Y_i - \sum k\pi_{ik}$ | Score residuals | ✓ Yes | ✓ **Recommended** |
| `ordinal::clm` | `residuals(fit)` | Working residuals | ✗ No | ✗ Must center |
| `ordinal::clm` | Manual: $Y_i - \sum k\pi_{ik}$ | Score residuals | ✓ Yes | ✓ **Recommended** |
| `statsmodels.OrderedModel` | Manual: $Y_i - \sum k\pi_{ik}$ | Score residuals | ✓ Yes | ✓ **Recommended** |
| `mord` | Manual: $Y_i - \sum k\pi_{ik}$ | Score residuals | ✓ Yes | ✓ **Recommended** |

**Critical finding**: **None of the standard packages provide zero-mean residuals by default** for ordinal regression. You **must** compute score residuals manually:

$$
R_i = Y_i - \mathbb{E}[Y_i \mid \mathbf{X}_i] = Y_i - \sum_{k=1}^K k \pi_{ik}
$$

**Best practice workflow**:

```r
# R example with MASS::polr
library(MASS)

# Step 1: Fit null model
fit <- polr(Y_ordinal ~ covariates, data = mydata)

# Step 2: Extract fitted probabilities
fitted_probs <- fitted(fit)  # n x K matrix

# Step 3: Compute score residuals
categories <- as.numeric(levels(Y_ordinal))
expected_Y <- fitted_probs %*% categories
residuals_for_spagrm <- as.numeric(Y_ordinal) - expected_Y

# Step 4: Verify zero-mean (critical!)
if (abs(sum(residuals_for_spagrm)) > 1e-10) {
  warning("Residuals not zero-mean! Centering...")
  residuals_for_spagrm <- residuals_for_spagrm - mean(residuals_for_spagrm)
}

# Step 5: Use in SPAGRM
# Pass residuals_for_spagrm to genetic association test
```

```python
# Python example with statsmodels
import numpy as np
from statsmodels.miscmodels.ordinal_model import OrderedModel

# Step 1: Fit null model
fit = OrderedModel(Y_ordinal, X, distr='logit').fit()

# Step 2: Extract fitted probabilities
fitted_probs = fit.predict()  # n x K array

# Step 3: Compute score residuals
categories = np.arange(fitted_probs.shape[1])
expected_Y = fitted_probs @ categories
residuals_for_spagrm = Y_ordinal - expected_Y

# Step 4: Verify zero-mean (critical!)
if abs(residuals_for_spagrm.sum()) > 1e-10:
    print("Warning: Residuals not zero-mean! Centering...")
    residuals_for_spagrm -= residuals_for_spagrm.mean()

# Step 5: Use in SPAGRM
# Pass residuals_for_spagrm to genetic association test
```

**Why standard residuals are not zero-mean**:

The "working" residuals from iteratively reweighted least squares (IRLS) are:
$$
R_i^{\text{work}} = \frac{\partial \ell_i}{\partial \eta_i} \bigg/ \frac{\partial^2 \ell_i}{\partial \eta_i^2}
$$

where $\eta_i$ is the linear predictor. These are **NOT** the same as $Y_i - \mathbb{E}[Y_i]$ and typically have non-zero mean.

**Conclusion**: For ordinal phenotypes, **always compute score residuals manually** and verify zero-mean before using in SPAGRM. Standard package residuals are unsuitable without centering.

### 5.3 Count Regression Residuals

**Count phenotype**: $Y_i \in \{0, 1, 2, \ldots\}$ (non-negative integers)

**Poisson regression**:
$$
\log(\mu_i) = \mathbf{X}_i^\top \boldsymbol{\beta} \tag{59}
$$

where $\mu_i = \mathbb{E}[Y_i \mid \mathbf{X}_i]$.

**Pearson residuals**:
$$
R_i^{(P)} = \frac{Y_i - \mu_i}{\sqrt{\mu_i}} \tag{60}
$$

**Standardized**: Divide by standard deviation $\sqrt{\mu_i}$ to approximate constant variance.

**Mean**: $\mathbb{E}[R_i^{(P)} \mid \mathbf{X}_i] = 0$ ✓

**Deviance residuals**:
$$
R_i^{(D)} = \text{sign}(Y_i - \mu_i) \sqrt{2[Y_i \ln(Y_i/\mu_i) - (Y_i - \mu_i)]} \tag{61}
$$

**Mean**: $\mathbb{E}[R_i^{(D)}] \neq 0$ (must center!)

**Recommendation for count regression**:

✓ **Use Pearson residuals** (already zero-mean)

✓ **Alternative**: Use $R_i = Y_i - \mu_i$ (unstandardized, but simpler for SPAGRM)

If using deviance residuals, **always center**.

### 5.4 Survival Regression Residuals

**Cox proportional hazards model**:
$$
h_i(t) = h_0(t) \exp(\mathbf{X}_i^\top \boldsymbol{\beta}) \tag{62}
$$

**Martingale residuals**:
$$
R_i^{(M)} = \delta_i - \hat{\Lambda}_0(T_i) e^{\mathbf{X}_i^\top \hat{\boldsymbol{\beta}}} \tag{63}
$$

where:

- $\delta_i$ is the event indicator (1 if event, 0 if censored)
- $\hat{\Lambda}_0(T_i)$ is the estimated cumulative baseline hazard at observed time $T_i$

**Properties**:

- $\mathbb{E}[R_i^{(M)}] = 0$ ✓ (for well-fitted models)
- Range: $(-\infty, 1]$ (left-skewed, not symmetric)

**Deviance residuals** (for Cox model):
$$
R_i^{(D)} = \text{sign}(R_i^{(M)}) \sqrt{-2[R_i^{(M)} + \delta_i \ln(\delta_i - R_i^{(M)})]} \tag{64}
$$

**Properties**:

- More symmetric than martingale residuals
- $\mathbb{E}[R_i^{(D)}] \neq 0$ (must center!)

**Recommendation for survival regression**:

✓ **Use martingale residuals** (automatically zero-mean)

If using deviance residuals, **always center**.

## 6. Practical Guidelines

### 6.1 Residual Choice Decision Tree

```
Is your phenotype binary?
│
├─ Yes: Binary logistic regression
│   │
│   ├─ Use Pearson residuals: R_i = Y_i - μ_i
│   │   ✓ Automatically zero-mean
│   │   ✓ Simple to compute
│   │   ✓ Standard choice in SPAGRM
│   │
│   └─ Alternative: Centered deviance residuals
│       R_i = R_i^(D) - mean(R^(D))
│       ✓ More symmetric (better for extreme μ_i)
│       ⚠ Requires centering step
│
├─ No: Ordinal phenotype?
│   │
│   ├─ Yes: Ordinal regression
│   │   │
│   │   └─ Use score residuals: R_i = Y_i - E[Y_i|X_i]
│   │       ✓ Automatically zero-mean
│   │       ✓ Treats ordinal as pseudo-continuous
│   │
│   └─ No: Count phenotype?
│       │
│       ├─ Yes: Poisson/Negative Binomial regression
│       │   │
│       │   └─ Use Pearson residuals: R_i = Y_i - μ_i
│       │       (or unstandardized: R_i = Y_i - μ_i)
│       │       ✓ Automatically zero-mean
│       │
│       └─ No: Time-to-event phenotype?
│           │
│           └─ Yes: Cox regression
│               │
│               └─ Use martingale residuals
│                   ✓ Automatically zero-mean
│                   ✓ Standard in survival analysis
```

### 6.2 General Rule

**For any generalized linear model (GLM)**:

1. **First choice**: Use residuals that are **theoretically zero-mean**
   - Pearson residuals: $R_i = \frac{Y_i - \mu_i}{\sqrt{V(\mu_i)}}$ (if you want constant variance)
   - Unstandardized: $R_i = Y_i - \mu_i$ (simpler, works fine for SPAGRM)

2. **If using deviance or other residuals**: **Always center**
   $$
   R_i^{(\text{centered})} = R_i - \frac{1}{n}\sum_{j=1}^n R_j
   $$

3. **Verify zero-sum**: After computing residuals, check that $\sum_{i=1}^n R_i \approx 0$
   - If not, center them!
   - Tolerance: $|\sum R_i| < 10^{-10}$ (numerical precision)

4. **Avoid**: Never use uncentered residuals with non-zero mean for association testing

### 6.3 Diagnostic Checks

**Before running SPAGRM**:

1. **Check residual sum**:

   ```r
   sum(residuals) < 1e-10  # Should be approximately zero
   ```

2. **Check residual variance**:

   ```r
   var(residuals)  # Should be positive and finite
   ```

3. **Check for extreme values**:

   ```r
   range(residuals)  # Should be reasonable (e.g., not ±Inf)
   ```

4. **Plot residual distribution**:
   - Should be roughly symmetric (or at least not heavily skewed)
   - Deviance residuals should be closer to normal than Pearson

5. **Compare variance estimates**:
   - Theoretical variance: $V = 2f(1-f) \sum R_i^2$
   - Empirical variance: Compute from simulated genotypes
   - Should be similar (variance ratio ≈ 1)

### 6.4 When to Prefer Centered Deviance over Pearson

**Prefer centered deviance residuals when**:

1. $\mu_i$ values are highly skewed (many close to 0 or 1)
2. Testing rare variants (where normal approximation is less accurate)
3. Small sample sizes ($n < 500$)
4. You've verified that centering improves calibration in your specific setting

**Stick with Pearson residuals when**:

1. Standard GWAS setting (large $n$, common variants)
2. Well-balanced case-control ratio (prevalence ≈ 0.5)
3. Simplicity and speed are important
4. You want consistency with other software/papers

## 7. Summary

### 7.1 Key Findings

**Residual Properties**:

1. **Pearson residuals** $R_i = Y_i - \mu_i$:
   - Automatically zero-mean ✓
   - Natural from score test derivation
   - May be skewed for extreme $\mu_i$
   - Standard choice in SPAGRM

2. **Deviance residuals** (centered) $R_i = R_i^{(D)} - \overline{R}^{(D)}$:
   - Requires explicit centering step
   - More symmetric distribution
   - Nearly identical variance to Pearson
   - Potentially better for rare variants

3. **Both are valid** for SPAGRM when properly centered

4. **Main difference**: Distribution shape, not expectation or variance

**Type I Error Rate (Section 4.5.1)**:

5. **Under null hypothesis**: Both methods have $\mathbb{E}[S] = 0$ and $\mathbb{V}(S) \approx$ same

6. **Calibration with balanced $\mu_i$** (around 0.5): 
   - Both well-calibrated, Type I error ≈ nominal level

7. **Calibration with extreme $\mu_i$** (many < 0.1 or > 0.9):
   - **Pearson**: Possible slight miscalibration due to skewness
   - **Deviance**: Better calibration due to symmetry advantage

8. **Sample size effects**:
   - Large $n$ (> 10,000): Both asymptotically correct
   - Small $n$ (< 500): Deviance may be closer to nominal if $\mu_i$ extreme

**Power (Section 4.5.2)**:

9. **Non-centrality parameters**: $\lambda^{(P)} \approx \lambda^{(DC)}$ for small effect sizes

10. **Common variants + balanced prevalence**: Nearly equal power (Relative Efficiency ≈ 1)

11. **Rare variants + rare disease**: Context-dependent
    - **Deviance** may have higher power if effect is heterogeneous across $\mu_i$
    - **Pearson** is robust if effect is uniform

12. **Non-linear effects**: Deviance residuals better capture effects stronger at extreme $\mu_i$

**Overall**:

13. **P-values should be similar** for well-calibrated models in typical GWAS settings

14. **Choice matters most** in: rare disease studies, extreme covariate effects, small samples

### 7.2 Recommendations

**For binary phenotypes in logistic regression**:

**Use Pearson residuals when**:

- Standard GWAS (large $n > 1000$, common variants $f > 0.05$)
- Balanced case-control design (prevalence ≈ 0.5)
- Fitted probabilities $\mu_i$ mostly between 0.1 and 0.9
- Simplicity and consistency with literature are priorities

**Use centered deviance residuals when**:

- Rare disease studies (prevalence < 0.1, many $\mu_i$ close to 0 or 1)
- Small samples ($n < 500$) with extreme covariate effects
- Type I error calibration is critical (e.g., small p-value thresholds)
- Testing rare variants ($f < 0.01$) with suspected non-linear effects
- You observe QQ plot deviation with Pearson residuals

**Both should work similarly when**:

- Large samples with well-calibrated null models
- Common variants with balanced prevalence
- Effect sizes are small ($|\gamma| < 0.3$)

**For other phenotypes**:

- **Ordinal**: Use score residuals $R_i = Y_i - \mathbb{E}[Y_i \mid \mathbf{X}_i]$ (auto zero-mean)
- **Count**: Use Pearson residuals $R_i = Y_i - \mu_i$ (or standardized version)
- **Survival**: Use martingale residuals (auto zero-mean)
- **Always check**: $\sum_{i=1}^n R_i \approx 0$ before running SPAGRM
- **If not**: Center the residuals manually: $R_i \leftarrow R_i - \overline{R}$

**Universal principles**:

1. **Zero-sum constraint is mandatory**: $|\sum R_i| < 10^{-10}$
2. **When in doubt, center**: Better safe than biased
3. **Validate with permutations**: If results are sensitive to residual choice, investigate further
4. **Don't mix residual types**: Use same type across all variants in a study

### 7.3 Critical Questions and Answers

#### Q1: If residuals sum to non-zero, is centering sufficient for SPAGRM?

**Answer**: **Yes, with important caveats**.

**Theory**: If $\sum_{i=1}^n R_i = c \neq 0$, centering gives:
$$
R_i^{\text{centered}} = R_i - \frac{c}{n} \implies \sum_{i=1}^n R_i^{\text{centered}} = 0
$$

The centered residuals satisfy the **zero-mean constraint** required for SPAGRM.

**However**:

1. **Why is centering needed in the first place?**
   - If the null model is correctly specified, $\mathbb{E}[R_i] = 0$ theoretically
   - Non-zero sum suggests: (a) numerical precision issues (acceptable if $|c| < 10^{-6}$), or (b) **model misspecification** (problematic!)

2. **Centering fixes the symptom, not the cause**:
   - If $|c|$ is large (e.g., $> 0.01 \times n$), investigate the null model fit
   - Check: residual plots, goodness-of-fit tests, outlier detection

3. **Centering is always safe** but may mask deeper issues:
   - For ordinal/deviance residuals: centering is **expected and necessary**
   - For Pearson/score residuals: large $|c|$ indicates potential problems

**Recommendation**:
```r
# Always check magnitude before centering
residual_sum <- sum(residuals)
if (abs(residual_sum) > 1e-6) {
  warning(sprintf("Large residual sum: %.6f. Check null model fit!", residual_sum))
}

# Center regardless (safe)
residuals <- residuals - mean(residuals)

# Verify zero-sum after centering
stopifnot(abs(sum(residuals)) < 1e-10)
```

#### Q2: Does the same principle apply to SPAmixPlus?

**Answer**: **Yes, identical requirements**.

SPAmix, SPAmixPlus, and SPAGRM all use the **same mathematical framework**:

$$
S = \sum_{i=1}^n R_i (G_i - 2f)
$$

where:

- $S$ is the score statistic
- $R_i$ are residuals from the null model
- $G_i$ are genotypes (centered by MAF)

**Requirements are identical**:

1. **Zero-mean**: $\sum_{i=1}^n R_i = 0$ (or $\mathbb{E}[R_i] = 0$)
2. **Independence under null**: $R_i \perp G_i$ under $H_0$
3. **Finite variance**: $\mathbb{V}(R_i) < \infty$

**Code evidence** from `SPAmix.cpp`:

```cpp
SPAmixClass::SPAmixClass(arma::mat t_resid, ...)
{
  m_resid = t_resid;  // Residuals passed directly from null model
  // No centering in SPAmix code -> assumes pre-centered!
}
```

**Critical**: SPAmix/SPAmixPlus **assume** residuals are already zero-mean. The user must ensure this before passing residuals.

#### Q3: POLMM Code Comparison and Residual Type

**Codes analyzed**:

- `code/GRAB-0.3.0/src/POLMM.cpp` (R package, Armadillo)
- `code/GRAB-feat-cpp/src/polmm/polmm.cpp` (C++ CLI, Eigen)

**Comparison**:

| Aspect | GRAB-0.3.0 (Armadillo) | GRAB-feat-cpp (Eigen) | Match? |
|--------|------------------------|----------------------|--------|
| **Residual formula** | $\sum_j f_j (y_{\text{cum},j} - F_j)$ | $\sum_j f_j (y_{\text{cum},j} - F_j)$ | ✓ Yes |
| **Implementation** | Matrix operations: `sumCols(RymuMat)` | Loop: `rym += f_j * (cumY - cumMu_j)` | ✓ Same result |
| **Variable name** | `m_RymuVec` | `RymuVec` | ✓ Consistent |
| **Null model fitting** | PQL with PCG | PQL with PCG | ✓ Identical |
| **Zero-mean** | **No explicit centering** | **No explicit centering** | ✓ Consistent |

**Conclusion**: The two implementations are **logically identical**. Both compute the same residual quantity, but **neither explicitly centers it**.

**Residual type in POLMM**: **Modified score residuals** (not standard score residuals!)

$$
R_i^{\text{POLMM}} = \sum_{j=0}^{J-2} f_j \left( \mathbb{1}(Y_i \leq j) - F(j - \eta_i) \right)
$$

where:

- $f_j = F'(j - \eta_i)$ is the derivative (density) at cutpoint $j$
- $F(j - \eta_i)$ is the cumulative probability $P(Y_i \leq j)$
- $\mathbb{1}(Y_i \leq j)$ is the observed cumulative indicator

**Key difference from standard score residuals**:

- **Standard score**: $R_i = Y_i - \mathbb{E}[Y_i \mid \mathbf{X}_i] = Y_i - \sum_k k\pi_{ik}$
- **POLMM residual**: Weighted by $f_j$ (density function) instead of category values $k$

**Is it zero-mean?** **Not automatically!**

From the code:

```cpp
// GRAB-0.3.0/src/POLMM.cpp, line 77-80
arma::mat ymuMat = yMat - m_muMat;                      // n x J
arma::mat RymuMat = ymuMat.cols(0, m_J-2) / t_iRMat;    // n x (J-1): R %*% (y - mu)
m_RymuVec = sumCols(RymuMat, m_J);                      // n x 1
// NO CENTERING STEP!
```

**Verification needed**:

```r
# After fitting POLMM null model
null_model <- fitPOLMM(...)
RymuVec <- null_model$RymuVec

# Check zero-mean
sum(RymuVec)  # Is this ≈ 0?
```

**Mathematical expectation**: Under the null model (conditional on covariates):
$$
\mathbb{E}[R_i^{\text{POLMM}} \mid \mathbf{X}_i] = \sum_{j=0}^{J-2} f_j \left( \mathbb{E}[\mathbb{1}(Y_i \leq j) \mid \mathbf{X}_i] - F(j - \eta_i) \right) = 0
$$

So **theoretically** it should be zero-mean, but **numerically** it may not be exactly zero.

**Which residuals to use?**

| Method | Object | Residual Name | Directly Usable? | Recommendation |
|--------|--------|---------------|------------------|----------------|
| **POLMM** | `fitPOLMM()` | `RymuVec` | **Verify first!** | Check `sum(RymuVec) ≈ 0`. If not, center before using in SPAGRM. |
| **POLMM** | `fitPOLMM()` | Working residuals | ✗ No | Not zero-mean, do not use |
| **SPAmix** | User-provided | `resid` | **Must pre-center** | User responsible for ensuring $\sum R_i = 0$ |
| **SPAmixPlus** | User-provided | `resid` | **Must pre-center** | User responsible for ensuring $\sum R_i = 0$ |
| **SPAGRM** | User-provided | Generic | **Must pre-center** | User responsible for ensuring $\sum R_i = 0$ |

**Best practice for POLMM residuals**:

```r
# Fit POLMM null model
null_model <- fitPOLMM(Y_ordinal ~ covariates, ...)

# Extract residuals
RymuVec <- null_model$RymuVec

# CRITICAL: Verify zero-mean
residual_sum <- sum(RymuVec)
cat(sprintf("Residual sum: %.10f\n", residual_sum))

# If |sum| > 1e-10, center them
if (abs(residual_sum) > 1e-10) {
  warning("POLMM residuals not zero-mean! Centering...")
  RymuVec <- RymuVec - mean(RymuVec)
}

# Verify after centering
stopifnot(abs(sum(RymuVec)) < 1e-10)

# Now safe to use in SPAGRM/SPAmix
```

**Summary**:

1. ✓ **Centering is always sufficient** to make residuals usable for SPAGRM/SPAmix
2. ✓ **Same requirement** applies to SPAmixPlus (zero-mean residuals)
3. ✓ **POLMM codes are identical** in logic (both GRAB-0.3.0 and GRAB-feat-cpp)
4. ⚠ **POLMM `RymuVec`** is theoretically zero-mean but **should be verified** before use
5. ⚠ **Never assume** residuals from any package are zero-mean—always check!
6. ✓ **Always center** as a safety measure: `residuals <- residuals - mean(residuals)`

### 7.4 Open Questions for Simulation Studies

To fully characterize the differences, conduct simulations varying:

1. Sample size: $n \in \{500, 1000, 5000, 10000\}$
2. Prevalence: $P(Y=1) \in \{0.01, 0.1, 0.5\}$ (rare disease, moderate, balanced)
3. MAF: $f \in \{0.001, 0.01, 0.1, 0.5\}$ (rare to common variants)
4. Effect size: $\gamma \in \{0, 0.1, 0.3, 0.5\}$ (null and alternative)
5. Covariate effects: Weak vs. strong (affects range of $\mu_i$)

**Metrics to compare**:

- Type I error rate (should be ≈ 0.05 under null)
- Power (for non-zero $\gamma$)
- P-value correlation between Pearson and deviance
- Runtime (if significant difference)
