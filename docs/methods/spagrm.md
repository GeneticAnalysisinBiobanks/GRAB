# SPAGRM

SPAGRM (Saddlepoint Approximation for Genetic Relationship Matrix) performs score tests for genetic variants accounting for sample relatedness through a sparse GRM and pairwise IBD (Identity-By-Descent) probabilities. SPAGRM extends SPACox to support related samples by modeling the joint genotype distribution within families.

Let:

- $n \times 1$ vector $\mathbf{G} = (G_1, \ldots, G_n)^\top$ represent genotypes (allele counts) for a variant to be tested;
- $n \times 1$ vector $\mathbf{R} = (R_1, \ldots, R_n)^\top$ represent residuals from the null model;
- $n \times n$ sparse matrix $\mathbf{\Phi}$ represent the GRM;
- $\beta$ denote the genetic effect to be tested;
- $f$ denote the minor allele frequency (MAF) for the variant.

The null hypothesis is $H_0: \beta = 0$.

## 1. Score Statistic

**Centering adjustment**: In null models (e.g., Cox regression), residuals $\mathbf{R}$ often have $\sum_{i=1}^n R_i \approx 0$ but not exactly zero. To ensure $\mathbb{E}[S] = 0$ under $H_0$ regardless of $\sum_{i=1}^n R_i$, the score is centered by the observed mean genotype:

$$
S = \mathbf{R}^\top \mathbf{G} - \overline{G} \cdot \sum_{i=1}^n R_i = \sum_{i=1}^n R_i (G_i - \overline{G}) \tag{1}
$$

where $\overline{G} = \frac{1}{n}\sum_{i=1}^n G_i$ is the **observed sample mean genotype**.

Under the null hypothesis with MAF $f$, the expected genotype is $\mathbb{E}[G_i] = 2f$. Since $\mathbb{E}[\overline{G}] = 2f$, we have

$$
\mathbb{E}[S] = \sum_{i=1}^n R_i \mathbb{E}[G_i - \overline{G}] \approx 0 \tag{2}
$$

**Critical Issue: MGF for Uncentered G vs. Observed Centered Score**
The observed score (1) uses the **sample mean** $\overline{G}$ (a random quantity) for centering. However, the MGF and CGF (Section 4) are derived for **uncentered** genotypes $G_i \sim \text{Binomial}(2, f)$. This creates a mathematical inconsistency:

- **Observed score**: $S = \sum R_i (G_i - \overline{G})$ where $\overline{G} = \frac{1}{n}\sum G_i$ is **random**
- **MGF**: $M(t) = \mathbb{E}[e^{t\sum R_i G_i}]$ is for **uncentered** $\sum R_i G_i$

**Why This Approximation Works**: The key assumption is that in a well-fitted null model, $\sum_{i=1}^n R_i \approx 0$. When this holds:

1. The centering term becomes negligible: 
   $$\overline{G} \cdot \sum_{i=1}^n R_i \approx 0$$

2. The observed score is approximately: 
   $$S \approx \sum_{i=1}^n R_i G_i$$

3. For large $n$, the sample mean converges to the expectation: 
   $$\overline{G} \xrightarrow{p} \mathbb{E}[G_i] = 2f$$

4. The CGF-based distribution approximates the distribution of $\sum R_i G_i$, which is close to $S$ when $\sum R_i \approx 0$

**Theoretical Concern**: Strictly speaking, the distribution of $S = \sum R_i (G_i - \overline{G})$ (where $\overline{G}$ is random) is **not identical** to the distribution of $\sum R_i (G_i - 2f)$ (where $2f$ is a constant). The SPAGRM approach is an **approximation** that relies on:

- $\sum R_i \approx 0$ (from well-fitted null model)
- Large sample size (so $\overline{G} \approx 2f$)

In practice, this approximation works well for typical GWAS applications where null models are well-calibrated. However, it may introduce slight inaccuracy when $\sum R_i$ is not sufficiently close to zero.

### 1.1 Genotype Centering with Zero-Mean Residuals

**Assumption**: Throughout this section, we assume:

- **Genotypes are centered by theoretical mean**: $\tilde{G}_i = G_i - 2f$ where $f$ is the minor allele frequency (MAF)
- **Residuals are empirically centered**: $\sum_{i=1}^n R_i = 0$ (exact zero mean)

**Question**: Under these assumptions, does centering the genotypes affect the distribution of the score statistic?

**Answer**: **No**. When residuals have exact zero mean, the score statistic has the same distribution whether we use centered or uncentered genotypes.

**Proof**:

Define two versions of the score statistic:

- **Centered genotypes**:
  $$S_{\text{cent}} = \sum_{i=1}^n R_i \tilde{G}_i = \sum_{i=1}^n R_i (G_i - 2f) = \sum_{i=1}^n R_i G_i - 2f \sum_{i=1}^n R_i \tag{1a}$$

- **Uncentered genotypes**: 
  $$S_{\text{uncent}} = \sum_{i=1}^n R_i G_i \tag{1b}$$

**Relationship**: Since $\sum_{i=1}^n R_i = 0$ (by assumption), we have:
$$
S_{\text{cent}} = \sum_{i=1}^n R_i G_i - 2f \cdot 0 = \sum_{i=1}^n R_i G_i = S_{\text{uncent}} \tag{1c}
$$

The two statistics are **identical**: $S_{\text{cent}} = S_{\text{uncent}}$.

**Implication for MGF**: Since the statistics are the same random variable, their moment generating functions are identical:
$$
M_{\text{cent}}(t) = \mathbb{E}[e^{t S_{\text{cent}}}] = \mathbb{E}[e^{t S_{\text{uncent}}}] = M_{\text{uncent}}(t) \tag{1d}
$$

**MGF derivation**:

For **uncentered genotypes** $G_i \sim \text{Binomial}(2, f)$ with $\mathbb{E}[G_i] = 2f$:
$$
M_{\text{uncent}}(t) = \mathbb{E}\left[\exp\left(t \sum_{i=1}^n R_i G_i\right)\right] = \prod_{i=1}^n \mathbb{E}[e^{t R_i G_i}] = \prod_{i=1}^n M_{G_i}(t R_i) \tag{1e}
$$

where $M_{G_i}(s) = (1-f + f e^s)^2$ is the MGF of $\text{Binomial}(2, f)$.

For **centered genotypes** $\tilde{G}_i = G_i - 2f$ with $\mathbb{E}[\tilde{G}_i] = 0$:
$$
M_{\text{cent}}(t) = \mathbb{E}\left[\exp\left(t \sum_{i=1}^n R_i \tilde{G}_i\right)\right] = \prod_{i=1}^n \mathbb{E}[e^{t R_i (G_i - 2f)}] = \prod_{i=1}^n e^{-t R_i \cdot 2f} M_{G_i}(t R_i) \tag{1f}
$$
$$
= \exp\left(-2ft \sum_{i=1}^n R_i\right) \prod_{i=1}^n M_{G_i}(t R_i) = \prod_{i=1}^n M_{G_i}(t R_i) \tag{1g}
$$

where the last equality uses $\sum_{i=1}^n R_i = 0$.

**Conclusion**: Both formulations give **identical MGFs**:
$$
M_{\text{cent}}(t) = M_{\text{uncent}}(t) = \prod_{i=1}^n (1-f + f e^{t R_i})^2 \tag{1h}
$$

**Cumulant generating function (CGF)**:
$$
K(t) = \ln M(t) = \sum_{i=1}^n 2\ln(1-f + f e^{t R_i}) \tag{1i}
$$

**Key takeaway**: When residuals are exactly centered ($\sum R_i = 0$), genotype centering is **mathematically irrelevant** for the score statistic and its distribution. The MGF, CGF, and saddlepoint approximation are all identical whether we center genotypes or not.

**Practical note**: In real data, $\sum R_i \approx 0$ but may not be exactly zero. SPAGRM uses **sample mean centering** $G_i - \overline{G}$ (equation 1) to ensure exact zero mean of the centered genotypes, which guarantees $\mathbb{E}[S] = 0$ regardless of small deviations in $\sum R_i$.

### 1.2 Impact of Missing Genotypes with Zero Imputation

**Setup**: Assume some genotypes are missing. We consider the imputation strategy of **replacing missing genotypes with zero** (after theoretical mean centering).

**Notation**:

- $\mathcal{O} = \{i : G_i \text{ is observed}\}$ with $|\mathcal{O}| = n_{\text{obs}}$
- $\mathcal{M} = \{i : G_i \text{ is missing}\}$ with $|\mathcal{M}| = n_{\text{miss}}$
- Total sample size: $n = n_{\text{obs}} + n_{\text{miss}}$

**Imputation strategy**: After centering by theoretical mean $2f$, set missing values to zero:
$$
\tilde{G}_i^* = \begin{cases}
G_i - 2f & \text{if } i \in \mathcal{O} \\
0 & \text{if } i \in \mathcal{M}
\end{cases}
\tag{1j}
$$

This is equivalent to **mean imputation** (replacing missing genotypes with $2f$, then centering all genotypes by $2f$).

**Score statistic with imputation**:
$$
S^* = \sum_{i=1}^n R_i \tilde{G}_i^* = \sum_{i \in \mathcal{O}} R_i (G_i - 2f) + \sum_{i \in \mathcal{M}} R_i \cdot 0 = \sum_{i \in \mathcal{O}} R_i (G_i - 2f) \tag{1k}
$$

**Comparison with complete case**: Complete case analysis uses only observed genotypes:
$$
S_{\text{drop}} = \sum_{i \in \mathcal{O}} R_i (G_i - 2f) \tag{1l}
$$

**Key observation**: The two statistics are **identical**:
$$
S^* = S_{\text{drop}} \tag{1m}
$$

**Implication**: With zero imputation (after theoretical mean centering), missing genotypes contribute **nothing** to the score statistic. This is mathematically equivalent to complete case analysis.

**Distribution analysis**:

**Expectation**: Under MCAR (missing completely at random) with $\Pr(i \in \mathcal{M}) = p_{\text{miss}}$ independent of genotypes:
$$
\mathbb{E}[S^*] = \mathbb{E}\left[\sum_{i \in \mathcal{O}} R_i (G_i - 2f)\right] = \sum_{i=1}^n R_i \mathbb{E}[\mathbb{1}(i \in \mathcal{O})] \mathbb{E}[G_i - 2f] = 0 \tag{1n}
$$

since $\mathbb{E}[G_i - 2f] = 0$ by Hardy-Weinberg equilibrium.

**Variance**: For unrelated individuals with $\text{Cov}(G_i, G_j) = 0$ for $i \neq j$:
$$
\mathbb{V}(S^*) = \mathbb{E}\left[\mathbb{V}\left(S^* \mid \mathcal{O}\right)\right] + \mathbb{V}\left(\mathbb{E}[S^* \mid \mathcal{O}]\right) \tag{1o}
$$

Given the observed set $\mathcal{O}$:
$$
\mathbb{V}(S^* \mid \mathcal{O}) = \sum_{i \in \mathcal{O}} R_i^2 \mathbb{V}(G_i) = \sum_{i \in \mathcal{O}} R_i^2 \cdot 2f(1-f) \tag{1p}
$$
$$
\mathbb{E}[S^* \mid \mathcal{O}] = 0 \tag{1q}
$$

Therefore:
$$
\mathbb{V}(S^*) = \mathbb{E}\left[\sum_{i \in \mathcal{O}} R_i^2 \cdot 2f(1-f)\right] = 2f(1-f) \sum_{i=1}^n R_i^2 \mathbb{E}[\mathbb{1}(i \in \mathcal{O})] \tag{1r}
$$
$$
= 2f(1-f) \sum_{i=1}^n R_i^2 (1 - p_{\text{miss}}) = 2f(1-f) (1-p_{\text{miss}}) \sum_{i=1}^n R_i^2 \tag{1s}
$$

**Comparison with complete data**: If all genotypes were observed:
$$
\mathbb{V}(S_{\text{full}}) = 2f(1-f) \sum_{i=1}^n R_i^2 \tag{1t}
$$

**Variance ratio**:
$$
\frac{\mathbb{V}(S^*)}{\mathbb{V}(S_{\text{full}})} = 1 - p_{\text{miss}} \tag{1u}
$$

**Interpretation**: Missing genotypes reduce the variance by a factor of $(1-p_{\text{miss}})$. This is expected because only $(1-p_{\text{miss}}) \times n$ genotypes contribute to the score.

**Impact on MGF and saddlepoint approximation**:

The MGF for the imputed score is:
$$
M^*(t) = \mathbb{E}\left[\exp\left(t \sum_{i \in \mathcal{O}} R_i (G_i - 2f)\right)\right] = \mathbb{E}_{\mathcal{O}}\left[\prod_{i \in \mathcal{O}} M_{\tilde{G}_i}(t R_i)\right] \tag{1v}
$$

where $M_{\tilde{G}_i}(s) = e^{-2fs} M_{G_i}(s) = e^{-2fs} (1-f + f e^s)^2$ is the MGF of the centered genotype.

Under MCAR, for each $i$, the contribution is:
$$
\mathbb{E}[e^{t R_i \tilde{G}_i^*}] = (1-p_{\text{miss}}) M_{\tilde{G}_i}(t R_i) + p_{\text{miss}} \cdot 1 = (1-p_{\text{miss}}) e^{-2f t R_i} (1-f + f e^{t R_i})^2 + p_{\text{miss}} \tag{1w}
$$

**CGF**:
$$
K_i(t) = \ln\left[(1-p_{\text{miss}}) e^{-2f t R_i} (1-f + f e^{t R_i})^2 + p_{\text{miss}}\right] \tag{1x}
$$

This is more complex than the complete-data CGF and does not have a closed-form expression.

**Practical implications**:

1. **Zero imputation = Complete case**: The score statistic is identical, but the interpretation differs:
   - Zero imputation: "Missing genotypes contribute zero"
   - Complete case: "Exclude missing genotypes"

2. **Variance estimation**: With zero imputation, the theoretical variance is:
   $$V^* = 2f(1-f) (1-p_{\text{miss}}) \sum_{i=1}^n R_i^2$$
   
   This must be adjusted for the missingness rate.

3. **Saddlepoint approximation**: The CGF with missing data (equation 1x) is more complex. In practice, SPAGRM can:
   - Use the complete-case CGF (simpler, exact for zero imputation)
   - Use empirical variance $V_{\text{emp}}$ to account for missingness (robust)

4. **Power loss**: The variance reduction by $(1-p_{\text{miss}})$ leads to power loss proportional to the missingness rate:
   - 5% missing → 5% variance reduction → minimal power loss
   - 20% missing → 20% variance reduction → noticeable power loss

**Conclusion**: Zero imputation of missing genotypes (after theoretical mean centering) is mathematically equivalent to complete case analysis. The main impact is a reduction in variance proportional to the fraction of observed genotypes.

### 1.3 Impact of Missing Residuals with Zero Imputation

**Setup**: Assume some residuals are missing. We consider the imputation strategy of **replacing missing residuals with zero** (which is the mean under the assumption $\sum R_i = 0$).

**Notation**:

- $\mathcal{O}_R = \{i : R_i \text{ is observed}\}$ with $|\mathcal{O}_R| = n_R$
- $\mathcal{M}_R = \{i : R_i \text{ is missing}\}$ with $|\mathcal{M}_R| = n - n_R$

**Imputation strategy**: Set missing residuals to zero:
$$
R_i^* = \begin{cases}
R_i & \text{if } i \in \mathcal{O}_R \\
0 & \text{if } i \in \mathcal{M}_R
\end{cases}
\tag{1y}
$$

**Score statistic with imputation**:
$$
S^* = \sum_{i=1}^n R_i^* \tilde{G}_i = \sum_{i \in \mathcal{O}_R} R_i (G_i - 2f) + \sum_{i \in \mathcal{M}_R} 0 \cdot (G_i - 2f) = \sum_{i \in \mathcal{O}_R} R_i (G_i - 2f) \tag{1z}
$$

**Comparison with complete case**: Complete case analysis uses only observed residuals:
$$
S_{\text{drop}} = \sum_{i \in \mathcal{O}_R} R_i (G_i - 2f) \tag{1aa}
$$

**Key observation**: The two statistics are **identical**:
$$
S^* = S_{\text{drop}} \tag{1ab}
$$

**Implication**: With zero imputation of residuals, missing residual positions contribute **nothing** to the score. This is mathematically equivalent to complete case analysis.

**Distribution analysis**:

**Expectation**: Under MCAR with $\Pr(i \in \mathcal{M}_R) = p_{\text{miss}}^R$ independent of genotypes and phenotypes:
$$
\mathbb{E}[S^*] = \mathbb{E}\left[\sum_{i \in \mathcal{O}_R} R_i (G_i - 2f)\right] = \sum_{i=1}^n \mathbb{E}[\mathbb{1}(i \in \mathcal{O}_R)] R_i \mathbb{E}[G_i - 2f] = 0 \tag{1ac}
$$

since $\mathbb{E}[G_i - 2f] = 0$.

**Variance**: For unrelated individuals:
$$
\mathbb{V}(S^*) = \mathbb{E}\left[\mathbb{V}(S^* \mid \mathcal{O}_R)\right] + \mathbb{V}\left(\mathbb{E}[S^* \mid \mathcal{O}_R]\right) \tag{1ad}
$$

Given the observed set $\mathcal{O}_R$:
$$
\mathbb{V}(S^* \mid \mathcal{O}_R) = \sum_{i \in \mathcal{O}_R} R_i^2 \mathbb{V}(G_i) = 2f(1-f) \sum_{i \in \mathcal{O}_R} R_i^2 \tag{1ae}
$$
$$
\mathbb{E}[S^* \mid \mathcal{O}_R] = 0 \tag{1af}
$$

Therefore:
$$
\mathbb{V}(S^*) = \mathbb{E}\left[2f(1-f) \sum_{i \in \mathcal{O}_R} R_i^2\right] = 2f(1-f) \sum_{i=1}^n R_i^2 \mathbb{E}[\mathbb{1}(i \in \mathcal{O}_R)] \tag{1ag}
$$
$$
= 2f(1-f) (1-p_{\text{miss}}^R) \sum_{i=1}^n R_i^2 \tag{1ah}
$$

**Comparison with complete data**:
$$
\mathbb{V}(S_{\text{full}}) = 2f(1-f) \sum_{i=1}^n R_i^2 \tag{1ai}
$$

**Variance ratio**:
$$
\frac{\mathbb{V}(S^*)}{\mathbb{V}(S_{\text{full}})} = 1 - p_{\text{miss}}^R \tag{1aj}
$$

**Interpretation**: Missing residuals reduce the variance by a factor of $(1-p_{\text{miss}}^R)$, similar to missing genotypes.

**Critical difference from missing genotypes**:

**Residuals are fixed quantities**, not random variables. When we compute the variance, we treat $R_i$ as **known constants**. Therefore:

1. **Missing genotypes**: Randomness comes from the Binomial distribution of $G_i$. Missingness affects which $G_i$ are observed, but the distribution of each $G_i$ is known.

2. **Missing residuals**: The residuals $R_i$ are **fixed outputs** from the null model fit. Missingness means we don't know certain $R_i$ values, which is fundamentally different from not observing a random draw.

**Impact on variance estimation**:

With missing residuals, the **true variance** depends on the unknown residuals $\{R_i : i \in \mathcal{M}_R\}$:
$$
\mathbb{V}(S^* \mid \mathbf{R}_{\text{full}}) = 2f(1-f) \sum_{i \in \mathcal{O}_R} R_i^2
$$

This is a **random quantity** because $\mathcal{O}_R$ is random (under MCAR). We can only estimate:
$$
\widehat{\mathbb{V}}(S^*) = 2f(1-f) \sum_{i \in \mathcal{O}_R} R_i^2 \tag{1ak}
$$

**Expected value of the variance estimator**:
$$
\mathbb{E}[\widehat{\mathbb{V}}(S^*)] = 2f(1-f) \sum_{i=1}^n R_i^2 \mathbb{E}[\mathbb{1}(i \in \mathcal{O}_R)] = 2f(1-f) (1-p_{\text{miss}}^R) \sum_{i=1}^n R_i^2 \tag{1al}
$$

This is **biased downward** by a factor of $(1-p_{\text{miss}}^R)$ relative to the complete-data variance!

**Correction**: To obtain an unbiased variance estimator, we should scale up:
$$
\widehat{\mathbb{V}}_{\text{corrected}}(S^*) = \frac{1}{1-p_{\text{miss}}^R} \cdot 2f(1-f) \sum_{i \in \mathcal{O}_R} R_i^2 \tag{1am}
$$

**Practical implications**:

1. **Zero imputation = Complete case**: As with missing genotypes, the score statistic is identical.

2. **Variance downward bias**: Unlike missing genotypes (where variance is inherently reduced), missing residuals create a **biased variance estimator** because we're missing information about the fixed weights $R_i$.

3. **Correction is difficult**: The correction in equation (1am) assumes:
   - MCAR mechanism (known $p_{\text{miss}}^R$)
   - Unrelated individuals (no need to account for GRM structure

   If residuals are missing due to covariate missingness, the mechanism may be MAR or MNAR, making correction unreliable.

4. **Complete case is preferred**: Unlike genotypes (where mean/zero imputation can work well), missing residuals should generally be handled by:
   - Excluding individuals with missing residuals from the analysis
   - Re-fitting the null model on the complete-case sample
   - Computing score statistics only for individuals with both observed residuals and genotypes

**Comparison table**:

| Property | Missing Genotypes (Zero Imputation) | Missing Residuals (Zero Imputation) |
|----------|-------------------------------------|-------------------------------------|
| **Score statistic** | $S^* = \sum_{i \in \mathcal{O}} R_i (G_i - 2f)$ | $S^* = \sum_{i \in \mathcal{O}_R} R_i (G_i - 2f)$ |
| **Equivalence** | Complete case | Complete case |
| **Expectation** | $\mathbb{E}[S^*] = 0$ | $\mathbb{E}[S^*] = 0$ |
| **True variance** | $2f(1-f)(1-p_{\text{miss}}) \sum R_i^2$ | $2f(1-f) \sum_{i \in \mathcal{O}_R} R_i^2$ (random!) |
| **Variance estimator** | Unbiased (after adjustment) | **Biased** (missing $R_i$ information) |
| **Variance reduction** | Expected (fewer data points) | Artifact of missingness |
| **Correction** | Straightforward | Difficult (requires MCAR assumption) |
| **Recommendation** | Acceptable with low missingness | **Avoid**—use complete case |

**Conclusion**:

- **Missing residuals are fundamentally different from missing genotypes** because residuals are fixed quantities, not random variables.
- **Zero imputation of residuals** leads to biased variance estimation unless we know the missingness rate and can assume MCAR.
- **Strong recommendation**: When residuals are missing (e.g., due to covariate missingness), re-fit the null model on the complete-case sample. Do not impute missing residuals.

### 1.4 Residual Choices for Binary Phenotypes in Logistic Regression

**Context**: For binary phenotypes, the null model is typically a logistic regression model. Unlike linear regression where residuals $R_i = Y_i - \mu_i$ naturally have mean zero, logistic regression offers multiple residual definitions with different properties.

**Null model**: Logistic regression with link function $\text{logit}(\mu_i) = \mathbf{X}_i^\top \boldsymbol{\beta}$, where:

- $\mu_i = \mathbb{E}[Y_i \mid \mathbf{X}_i] = \frac{e^{\mathbf{X}_i^\top \boldsymbol{\beta}}}{1 + e^{\mathbf{X}_i^\top \boldsymbol{\beta}}}$ is the fitted probability
- $Y_i \in \{0, 1\}$ is the binary phenotype

We compare three residual definitions and their impact on the score statistic distribution.

#### 1.4.1 Pearson Residuals: $R_i = Y_i - \mu_i$

**Definition**: The **Pearson residual** (also called response residual):
$$
R_i^{(P)} = Y_i - \mu_i \tag{1an}
$$

**Properties**:

- **Mean**: $\mathbb{E}[R_i^{(P)}] = \mathbb{E}[Y_i - \mu_i] = \mu_i - \mu_i = 0$ (automatically zero mean!)
- **Variance**: $\mathbb{V}(R_i^{(P)}) = \mathbb{V}(Y_i) = \mu_i(1-\mu_i)$ (heteroscedastic)
- **Sum**: $\sum_{i=1}^n R_i^{(P)} = \sum_{i=1}^n (Y_i - \mu_i) \approx 0$ for well-fitted models

**Score statistic**:
$$
S^{(P)} = \sum_{i=1}^n R_i^{(P)} (G_i - 2f) = \sum_{i=1}^n (Y_i - \mu_i) (G_i - 2f) \tag{1ao}
$$

**Expectation**: 
$$
\mathbb{E}[S^{(P)}] = \sum_{i=1}^n \mathbb{E}[Y_i - \mu_i] \mathbb{E}[G_i - 2f] = 0 \tag{1ap}
$$

since both $Y_i - \mu_i$ and $G_i - 2f$ have zero mean (under the null hypothesis of no genetic association).

**Variance**:
$$
\mathbb{V}(S^{(P)}) = 2f(1-f) \sum_{i=1}^n (R_i^{(P)})^2 = 2f(1-f) \sum_{i=1}^n (Y_i - \mu_i)^2 \tag{1aq}
$$

for unrelated individuals (ignoring GRM structure for simplicity).

**Advantage**: Simple, automatically has zero mean, widely used in practice.

**Disadvantage**: Residuals are heteroscedastic (variance depends on $\mu_i$), which can affect efficiency.

#### 1.4.2 Deviance Residuals (Likelihood Residuals): Uncentered

**Definition**: The **deviance residual** (or likelihood residual):
$$
R_i^{(D)} = \text{sign}(Y_i - \mu_i) \sqrt{2\left[Y_i \ln\frac{Y_i}{\mu_i} + (1-Y_i) \ln\frac{1-Y_i}{1-\mu_i}\right]} \tag{1ar}
$$

For binary outcomes:

- If $Y_i = 1$: $R_i^{(D)} = \sqrt{-2\ln \mu_i}$
- If $Y_i = 0$: $R_i^{(D)} = -\sqrt{-2\ln(1-\mu_i)}$

**Properties**:

- **Mean**: $\mathbb{E}[R_i^{(D)}] \neq 0$ in general! The deviance residuals are **not** zero-mean.
- **Sum**: $\sum_{i=1}^n R_i^{(D)} \neq 0$ (typically non-zero)
- **Motivation**: Deviance residuals are more symmetric and closer to normality than Pearson residuals

**Empirical mean**: For a specific dataset:
$$
\overline{R}^{(D)} = \frac{1}{n}\sum_{i=1}^n R_i^{(D)} \tag{1as}
$$

**Score statistic (uncentered)**:
$$
S^{(D)} = \sum_{i=1}^n R_i^{(D)} (G_i - 2f) \tag{1at}
$$

**Expectation**: Under the null hypothesis (no genetic association, $G_i \perp Y_i \mid \mathbf{X}_i$):
$$
\mathbb{E}[S^{(D)}] = \sum_{i=1}^n \mathbb{E}[R_i^{(D)}] \mathbb{E}[G_i - 2f] = 0 \tag{1au}
$$

**Key observation**: Even though $\mathbb{E}[R_i^{(D)}] \neq 0$ for individual $i$, the score statistic still has zero expectation because $\mathbb{E}[G_i - 2f] = 0$.

**However**, the **empirical score** using observed genotypes:
$$
S_{\text{obs}}^{(D)} = \sum_{i=1}^n R_i^{(D)} (G_i^{\text{obs}} - 2f)
$$

may have non-zero expected value if $\sum_{i=1}^n R_i^{(D)} \neq 0$ and the sample mean $\overline{G}^{\text{obs}}$ deviates from $2f$.

**Problem**: The non-zero mean of deviance residuals can introduce bias when:

1. The sample size is small
2. The genotypes are unbalanced (e.g., rare variants)
3. The observed sample mean $\overline{G}$ differs from the expected $2f$

**Variance**:
$$
\mathbb{V}(S^{(D)}) = 2f(1-f) \sum_{i=1}^n (R_i^{(D)})^2 + \text{Cov terms} \tag{1av}
$$

The variance is similar to Pearson residuals but deviance residuals may have different scaling.

#### 1.4.3 Deviance Residuals: Manually Centered

**Definition**: To ensure zero mean, manually center the deviance residuals:
$$
R_i^{(DC)} = R_i^{(D)} - \overline{R}^{(D)} \tag{1aw}
$$

where $\overline{R}^{(D)} = \frac{1}{n}\sum_{i=1}^n R_i^{(D)}$ is the sample mean.

**Properties**:

- **Mean**: $\sum_{i=1}^n R_i^{(DC)} = 0$ (exactly zero by construction!)
- **Individual means**: $\mathbb{E}[R_i^{(DC)}] = \mathbb{E}[R_i^{(D)}] - \mathbb{E}[\overline{R}^{(D)}]$ (approximately zero for large $n$)

**Score statistic (centered)**:
$$
S^{(DC)} = \sum_{i=1}^n R_i^{(DC)} (G_i - 2f) = \sum_{i=1}^n (R_i^{(D)} - \overline{R}^{(D)}) (G_i - 2f) \tag{1ax}
$$

**Relationship to uncentered score**:
$$
S^{(DC)} = S^{(D)} - \overline{R}^{(D)} \sum_{i=1}^n (G_i - 2f) \tag{1ay}
$$

Since $\sum_{i=1}^n (G_i - 2f) = \sum_{i=1}^n G_i - n \cdot 2f = n(\overline{G} - 2f)$:
$$
S^{(DC)} = S^{(D)} - n \overline{R}^{(D)} (\overline{G} - 2f) \tag{1az}
$$

**Key observation**:

- If $\overline{G} \approx 2f$ (genotype sample mean close to expectation), then $S^{(DC)} \approx S^{(D)}$
- If $\overline{G} \neq 2f$ (e.g., rare variants with few minor alleles), centering can make a difference

**Expectation**: Under the null:
$$
\mathbb{E}[S^{(DC)}] = \mathbb{E}[S^{(D)}] - n \mathbb{E}[\overline{R}^{(D)}] \mathbb{E}[\overline{G} - 2f] = 0 \tag{1ba}
$$

Both uncentered and centered scores have zero expectation under the null.

**Variance**:
$$
\mathbb{V}(S^{(DC)}) = \mathbb{V}\left(\sum_{i=1}^n R_i^{(D)} (G_i - 2f)\right) + \mathbb{V}\left(n \overline{R}^{(D)} (\overline{G} - 2f)\right) - 2\text{Cov}(\cdots) \tag{1bb}
$$

For large $n$ with independent observations:
$$
\mathbb{V}(S^{(DC)}) \approx \mathbb{V}(S^{(D)}) - \frac{1}{n}\mathbb{V}(R^{(D)}) \mathbb{V}(G) \tag{1bc}
$$

The centering **reduces** the variance slightly (removes a finite-sample correction term).

#### Comparison Summary

| Property | Pearson $R_i = Y_i - \mu_i$ | Deviance (Uncentered) | Deviance (Centered) |
|----------|----------------------------|----------------------|---------------------|
| **Mean** | $\mathbb{E}[R_i] = 0$ | $\mathbb{E}[R_i] \neq 0$ | $\sum R_i = 0$ (empirical) |
| **Score expectation** | $\mathbb{E}[S] = 0$ | $\mathbb{E}[S] = 0$ | $\mathbb{E}[S] = 0$ |
| **Score (empirical sum)** | $\sum R_i \approx 0$ | $\sum R_i \neq 0$ | $\sum R_i = 0$ (exactly) |
| **Bias for rare variants** | Minimal | Potential (if $\overline{G} \neq 2f$) | Minimal |
| **Variance** | Baseline | Similar | Slightly reduced |
| **Symmetry/Normality** | Less symmetric | More symmetric | More symmetric |
| **Recommendation** | **Standard choice** | Avoid | Acceptable alternative |

#### Detailed Impact on SPA

**MGF derivation**: All three residual choices lead to the same MGF formula (equation 1e) because the MGF depends on the **distribution** of $G_i$, not the residuals:
$$
M(t) = \prod_{i=1}^n (1-f + f e^{t R_i})^2
$$

The residuals $R_i$ are **fixed quantities** (computed from the fitted null model), so they only affect the **weights** in the MGF, not the form.

**Saddlepoint equation**: For all three choices, we solve:
$$
K'(\hat{t}) = S_{\text{obs}}
$$

where $S_{\text{obs}}$ is the observed score. The key difference is:

- **Pearson**: $S_{\text{obs}}^{(P)}$ naturally has $\sum R_i^{(P)} \approx 0$, consistent with the MGF derivation assumption
- **Deviance (uncentered)**: $S_{\text{obs}}^{(D)}$ may have $\sum R_i^{(D)} \neq 0$, which can cause slight bias in rare variant tests
- **Deviance (centered)**: $S_{\text{obs}}^{(DC)}$ has $\sum R_i^{(DC)} = 0$ exactly, eliminating the bias

**Practical recommendation**:

1. **Use Pearson residuals** $R_i = Y_i - \mu_i$ as the default choice (standard in SPAGRM and most GWAS software)
2. If deviance residuals are preferred (for better normality), **always center them**: $R_i^{(DC)} = R_i^{(D)} - \overline{R}^{(D)}$
3. Never use uncentered deviance residuals for association testing, especially with rare variants

**Why this matters**: From Section 1.1, we know that when $\sum R_i = 0$, genotype centering is mathematically irrelevant. Using Pearson or centered deviance residuals ensures this condition holds, making the SPA more accurate.

## 2. Score Variance

The variance of the score statistic accounts for genotype correlation induced by the GRM:

$$
\mathbb{V}(S) = \mathbb{V}\left(\sum_{i=1}^n R_i G_i\right) = \sum_{i=1}^n \sum_{j=1}^n R_i R_j \text{Cov}(G_i, G_j) \tag{3}
$$

Under Hardy-Weinberg equilibrium and assuming the GRM $\mathbf{\Phi}$ represents the genotypic correlation structure:

$$
\text{Cov}(G_i, G_j) = 2f(1-f) \cdot \Phi_{ij} \tag{4}
$$

where $\Phi_{ii} \approx 1$ (diagonal elements) and $\Phi_{ij}$ for $i \neq j$ captures the kinship between individuals $i$ and $j$.

Substituting (4) into (3):

$$
\mathbb{V}(S) = 2f(1-f) \sum_{i=1}^n \sum_{j=1}^n R_i \Phi_{ij} R_j = 2f(1-f) \cdot \mathbf{R}^\top \mathbf{\Phi} \mathbf{R} \tag{5}
$$

Define the genotype variance as

$$
\sigma_G^2 = 2f(1-f) \tag{6}
$$

Then the score variance simplifies to

$$
\mathbb{V}(S) = \sigma_G^2 \cdot \mathbf{R}^\top \mathbf{\Phi} \mathbf{R} \tag{7}
$$

The standardized score (z-score) is

$$
Z = \frac{S}{\sqrt{\mathbb{V}(S)}} = \frac{S}{\sqrt{\sigma_G^2 \cdot \mathbf{R}^\top \mathbf{\Phi} \mathbf{R}}} \tag{8}
$$

**Note on centering**: The variance formula (7) applies to both centered and uncentered scores because $\mathbb{V}(S) = \mathbb{V}(S - \text{constant})$. The centering in (1) only affects the mean, not the variance.

## 3. Family Structure and Outlier Partitioning

To improve numerical stability and accuracy of the saddlepoint approximation, SPAGRM partitions the sample into **outliers** and **non-outliers** based on residuals, and further decomposes outliers into family structures based on the sparse GRM connectivity.

### 3.1 Outlier Detection

Let $Q_1, Q_3$ denote the 25th and 75th percentiles of residuals, and $\text{IQR} = Q_3 - Q_1$. Outliers are defined as subjects with residuals outside the range

$$
[\, Q_1 - r \cdot \text{IQR}, \; Q_3 + r \cdot \text{IQR} \,] \tag{10}
$$

where $r = 1.5$ is the outlier ratio (adjustable to control outlier proportion).

### 3.2 Family Decomposition via GRM Connectivity

The sparse GRM $\mathbf{\Phi}$ defines a graph where nodes are subjects and edges exist for non-zero off-diagonal entries. SPAGRM decomposes this graph into connected components:

- **Singletons**: isolated subjects (unrelated outliers) with no GRM edges to other outliers
- **Two-subject families**: outlier pairs $(i,j)$ with $\Phi_{ij} \neq 0$
- **Three-or-more-subject families**: larger connected components with $\geq 3$ outliers

Let:

- $\mathcal{U}$ denote the set of unrelated outlier indices (singletons)
- $\mathcal{F}_2^{(k)}$ denote the $k$-th two-subject family
- $\mathcal{F}_3^{(k)}$ denote the $k$-th family with $\geq 3$ subjects
- $\mathcal{N}$ denote the set of non-outlier indices

The score decomposes as

$$
S = \sum_{i \in \mathcal{U}} R_i G_i + \sum_{k} \sum_{i \in \mathcal{F}_2^{(k)}} R_i G_i + \sum_{k} \sum_{i \in \mathcal{F}_3^{(k)}} R_i G_i + \sum_{i \in \mathcal{N}} R_i G_i \tag{11}
$$

## 4. Moment Generating Function (MGF)

The saddlepoint approximation requires the cumulant generating function (CGF), which is derived from the moment generating function (MGF). SPAGRM computes the MGF separately for each family structure.

### 4.1 Single Genotype MGF

For a single genotype $G_i$ with MAF $f$, following Hardy-Weinberg equilibrium:

$$
G_i \sim \text{Binomial}(2, f), \quad \Pr(G_i = k) = \binom{2}{k} f^k (1-f)^{2-k}, \quad k \in \{0,1,2\} \tag{12}
$$

The MGF of $G_i$ is

$$
M_G(t; f) = \mathbb{E}[e^{tG_i}] = (1-f)^2 + 2f(1-f)e^t + f^2 e^{2t} = (1-f + fe^t)^2 \tag{13}
$$

The first and second derivatives are

$$
M_G'(t; f) = 2fe^t(1-f + fe^t) \tag{14}
$$

$$
M_G''(t; f) = 2fe^t(1-f + 2fe^t) \tag{15}
$$

For a weighted genotype $R_i G_i$, the MGF is $M_G(tR_i; f)$, with derivatives

$$
\frac{d}{dt}M_G(tR_i; f) = R_i M_G'(tR_i; f), \quad \frac{d^2}{dt^2}M_G(tR_i; f) = R_i^2 M_G''(tR_i; f) \tag{16}
$$

The CGF (cumulant generating function) is $K_G^{(i)}(t) = \ln M_G(tR_i; f)$. Expanding in terms of $t, R_i, f$:

$$
K_G^{(i)}(t) = 2\ln(1 - f + fe^{tR_i}) \tag{17}
$$

The first derivative is

$$
K_G^{(i)\prime}(t) = \frac{M_G'(tR_i; f)}{M_G(tR_i; f)} = \frac{2fR_ie^{tR_i}(1 - f + fe^{tR_i})}{(1 - f + fe^{tR_i})^2} = \frac{2fR_ie^{tR_i}}{1 - f + fe^{tR_i}} \tag{18}
$$

The second derivative is

$$
\begin{align}
K_G^{(i)\prime\prime}(t) &= \frac{M_G(tR_i; f)M_G''(tR_i; f) - [M_G'(tR_i; f)]^2}{[M_G(tR_i; f)]^2} \tag{19} \\
&= \frac{2fR_i^2e^{tR_i}(1-f)}{(1 - f + fe^{tR_i})^2} \tag{20}
\end{align}
$$

### 4.2 Unrelated Outliers (Singletons)

For unrelated outliers, genotypes are independent. The joint MGF is the product of individual MGFs:

$$
M_{\mathcal{U}}(t) = \prod_{i \in \mathcal{U}} M_G(tR_i; f) = \prod_{i \in \mathcal{U}} (1-f + fe^{tR_i})^2 \tag{21}
$$

Define $\lambda_i = e^{tR_i}$ and $\alpha_i = 1-f + f\lambda_i = 1-f + fe^{tR_i}$. Then for each outlier $i$:

$$
M_G^{(i)}(t) = \alpha_i^2 \tag{22}
$$

$$
M_G^{(i)\prime}(t) = 2f R_i \lambda_i \alpha_i = 2\alpha_i \cdot f R_i \lambda_i \tag{23}
$$

$$
M_G^{(i)\prime\prime}(t) = 2\left[(fR_i\lambda_i)^2 + fR_i\lambda_i \alpha_i\right] = 2\left[\alpha_i \cdot (fR_i\lambda_i) \cdot R_i + (fR_i\lambda_i)^2\right] \tag{24}
$$

The CGF for all unrelated outliers is

$$
K_{\mathcal{U}}(t) = \sum_{i \in \mathcal{U}} K_G^{(i)}(t) = \sum_{i \in \mathcal{U}} 2\ln(1 - f + fe^{tR_i}) \tag{25}
$$

with derivatives

$$
K_{\mathcal{U}}'(t) = \sum_{i \in \mathcal{U}} \frac{2fR_ie^{tR_i}}{1 - f + fe^{tR_i}} = \sum_{i \in \mathcal{U}} \frac{M_G^{(i)\prime}(t)}{M_G^{(i)}(t)} \tag{26}
$$

$$
K_{\mathcal{U}}''(t) = \sum_{i \in \mathcal{U}} \frac{2fR_i^2e^{tR_i}(1-f)}{(1 - f + fe^{tR_i})^2} = \sum_{i \in \mathcal{U}} \left[\frac{M_G^{(i)\prime\prime}(t)}{M_G^{(i)}(t)} - \left(\frac{M_G^{(i)\prime}(t)}{M_G^{(i)}(t)}\right)^2\right] \tag{27}
$$

### 4.3 Two-Subject Families

For a two-subject family with residuals $(R_1, R_2)$ and pairwise IBD probabilities $(p_a, p_b, p_c)$ where:

- $p_a$ = probability of sharing 2 alleles IBD
- $p_b$ = probability of sharing 1 allele IBD  
- $p_c$ = probability of sharing 0 alleles IBD

satisfying $p_a + p_b + p_c = 1$.

The joint genotype distribution accounts for IBD sharing. Define the kinship coefficient

$$
\rho = p_a + \frac{1}{2}p_b \tag{28}
$$

The correlation between genotypes is modeled through $\rho$. The joint MGF for genotypes $(G_1, G_2)$ is computed by summing over all $3 \times 3 = 9$ genotype configurations weighted by their IBD-adjusted probabilities.

Let $\tau = (1 - \rho) f(1-f)$ be the effective independent variance component. The MGF components are:

$$
\begin{align}
m_1 &= e^{tR_1} \tau \tag{29} \\
m_2 &= e^{tR_2} \tau \tag{30} \\
m_3 &= e^{t(R_1 + R_2)} (f - \tau) \tag{31}
\end{align}
$$

The joint MGF for the two-subject family is

$$
M_{\mathcal{F}_2}(t) = m_1 + m_2 + m_3 - \tau + (1-f) \tag{32}
$$

Expanding in terms of $t, R_1, R_2, f, \rho$:

$$
M_{\mathcal{F}_2}(t) = e^{tR_1}\tau + e^{tR_2}\tau + e^{t(R_1+R_2)}(f - \tau) - \tau + (1-f) \tag{33}
$$

The first and second derivatives are

$$
M_{\mathcal{F}_2}'(t) = R_1 m_1 + R_2 m_2 + (R_1 + R_2) m_3 = R_1e^{tR_1}\tau + R_2e^{tR_2}\tau + (R_1+R_2)e^{t(R_1+R_2)}(f - \tau) \tag{34}
$$

$$
M_{\mathcal{F}_2}''(t) = R_1^2 m_1 + R_2^2 m_2 + (R_1 + R_2)^2 m_3 = R_1^2e^{tR_1}\tau + R_2^2e^{tR_2}\tau + (R_1+R_2)^2e^{t(R_1+R_2)}(f - \tau) \tag{35}
$$

The CGF for the two-subject family is

$$
K_{\mathcal{F}_2}(t) = \ln M_{\mathcal{F}_2}(t) = \ln\left[e^{tR_1}\tau + e^{tR_2}\tau + e^{t(R_1+R_2)}(f - \tau) - \tau + (1-f)\right] \tag{36}
$$

with derivatives

$$
K_{\mathcal{F}_2}'(t) = \frac{M_{\mathcal{F}_2}'(t)}{M_{\mathcal{F}_2}(t)} = \frac{R_1e^{tR_1}\tau + R_2e^{tR_2}\tau + (R_1+R_2)e^{t(R_1+R_2)}(f - \tau)}{e^{tR_1}\tau + e^{tR_2}\tau + e^{t(R_1+R_2)}(f - \tau) - \tau + (1-f)} \tag{37}
$$

$$
K_{\mathcal{F}_2}''(t) = \frac{M_{\mathcal{F}_2}(t)M_{\mathcal{F}_2}''(t) - [M_{\mathcal{F}_2}'(t)]^2}{[M_{\mathcal{F}_2}(t)]^2} \tag{38}
$$

**Interpretation**: The term $\tau = (1 - \rho)f(1-f)$ reduces the effective genotype variance due to correlation induced by relatedness. As $\rho \to 1$ (full siblings), $\tau \to 0$ and genotypes become fully correlated; as $\rho \to 0$ (unrelated), $\tau \to f(1-f)$ and genotypes become independent.

### 4.4 Three-or-More-Subject Families

For families with $n_k \geq 3$ related outliers (where $n_k$ is the number of family members in family $k$), the full joint distribution over $3^{n_k}$ genotype configurations is intractable. SPAGRM uses a **Chow-Liu tree approximation** to construct a tractable factorization of the joint PMF.

#### 4.4.1 Chow-Liu Tree Construction

The Chow-Liu algorithm builds a maximum spanning tree on the family graph by maximizing mutual information between genotypes, given the pairwise IBD probabilities.

**Mutual Information**: For each pair $(i,j)$ in family $k$, compute the mutual information

$$
I(G_i; G_j) = \sum_{a=0}^{2} \sum_{b=0}^{2} \Pr(G_i=a, G_j=b) \log\frac{\Pr(G_i=a, G_j=b)}{\Pr(G_i=a)\Pr(G_j=b)}
$$

where $\Pr(G_i=a, G_j=b)$ is the IBD-adjusted joint probability (computed using IBD probabilities $\phi_{ij}^{(0)}, \phi_{ij}^{(1)}, \phi_{ij}^{(2)}$ and Hardy-Weinberg marginals), and $\Pr(G_i=a) = \binom{2}{a}f^a(1-f)^{2-a}$ is the Hardy-Weinberg marginal.

The Chow-Liu tree is the maximum spanning tree with edge weights equal to $I(G_i; G_j)$.

**Tree Factorization**: Let $\mathbf{G}_k = (G_{i_1}, G_{i_2}, \ldots, G_{i_{n_k}})$ denote the **random vector of genotypes** for all $n_k$ members in family $k$ (each $G_i$ is a random variable taking values in $\{0,1,2\}$). The tree factorization approximates the joint PMF as

$$
\Pr(\mathbf{G}_k = \mathbf{g}) \approx \Pr(G_{r} = g_r) \prod_{(i,j) \in \text{Tree}} \frac{\Pr(G_i = g_i, G_j = g_j)}{\Pr(G_i = g_i)} \tag{39}
$$

where:

- $\mathbf{g} = (g_{i_1}, \ldots, g_{i_{n_k}})$ is a **specific genotype configuration** (a realization of $\mathbf{G}_k$), with each $g_i \in \{0,1,2\}$
- $r$ is the **root node** of the tree (arbitrarily chosen from the $n_k$ family members)
- $(i,j) \in \text{Tree}$ denotes the set of edges in the Chow-Liu tree
- This gives the **joint probability mass** for the specific configuration $\mathbf{g}$ over all $n_k$ genotypes

**Expanded Form with IBD Probabilities**: The pairwise probabilities $\Pr(G_i = g_i, G_j = g_j)$ in (39) are computed using IBD-adjusted formulae. For each pair $(i,j)$ with IBD probabilities $(\phi_{ij}^{(0)}, \phi_{ij}^{(1)}, \phi_{ij}^{(2)})$ (probability of sharing 0, 1, 2 alleles IBD), the joint probability is:

$$
\Pr(G_i = a, G_j = b | f, \phi_{ij}^{(0)}, \phi_{ij}^{(1)}, \phi_{ij}^{(2)}) = \sum_{k=0}^{2} \phi_{ij}^{(k)} \Pr_{\text{IBD}=k}(G_i = a, G_j = b | f) \tag{39a}
$$

where $\Pr_{\text{IBD}=k}(G_i = a, G_j = b | f)$ is the genotype probability conditional on sharing exactly $k$ alleles IBD:

- **IBD = 0** (unrelated): 
  $$\Pr_{\text{IBD}=0}(G_i = a, G_j = b | f) = \Pr(G_i=a|f) \cdot \Pr(G_j=b|f) = \binom{2}{a}f^a(1-f)^{2-a} \cdot \binom{2}{b}f^b(1-f)^{2-b}$$

- **IBD = 1** (share 1 allele): Requires summing over which allele is shared (see SPAGRM paper Appendix)

- **IBD = 2** (identical genotypes): 
  $$\Pr_{\text{IBD}=2}(G_i = a, G_j = b | f) = \begin{cases} \binom{2}{a}f^a(1-f)^{2-a} & \text{if } a = b \\ 0 & \text{if } a \neq b \end{cases}$$

Substituting (39a) into (39) gives the full expansion:

$$
\Pr(\mathbf{G}_k = \mathbf{g} | f, \{\phi_{ij}\}) \approx \binom{2}{g_r}f^{g_r}(1-f)^{2-g_r} \prod_{(i,j) \in \text{Tree}} \frac{\sum_{k=0}^{2} \phi_{ij}^{(k)} \Pr_{\text{IBD}=k}(G_i = g_i, G_j = g_j | f)}{\binom{2}{g_i}f^{g_i}(1-f)^{2-g_i}} \tag{39b}
$$

where $\{\phi_{ij}\}$ denotes the collection of IBD probability triples for all pairs in the family.

#### 4.4.2 Standardized Score and Probability Array

For each family $k$ with $n_k$ members, SPAGRM precomputes:

- **Standardized scores**: For each of the $3^{n_k}$ genotype configurations $\mathbf{g} = (g_1, \ldots, g_{n_k})$ where $g_i \in \{0, 1, 2\}$ is a **specific realization** (observed value) of the random variable $G_i$, compute

  $$
  S_k(\mathbf{g}) = \sum_{i \in \mathcal{F}_3^{(k)}} R_i g_i \tag{40}
  $$
  
  **Notation**: We use lowercase $g_i$ for the **realized genotype value** (0, 1, or 2) and uppercase $G_i$ for the **random variable**. In equation (40), $g_i$ is the specific genotype count in configuration $\mathbf{g}$, not a random variable.
  
  This produces $3^{n_k}$ distinct score values $\{S_k(\mathbf{g}_1), S_k(\mathbf{g}_2), \ldots, S_k(\mathbf{g}_{3^{n_k}})\}$, one for each possible genotype configuration.

- **Probability array**: Define $p(\mathbf{g}; f)$ as the **probability mass function** of the random vector $\mathbf{G}_k$:
  
  $$
  p(\mathbf{g}; f) := \Pr(\mathbf{G}_k = \mathbf{g} | \text{MAF} = f, \{\phi_{ij}\})
  $$
  
  where $\{\phi_{ij}\}$ represents the IBD probabilities for all pairs $(i,j)$ in family $k$. **Note**: The notation $p(\mathbf{g}; f)$ suppresses the dependence on IBD probabilities for brevity, but the probability **does depend on within-family relatedness** through the IBD-adjusted pairwise probabilities in (39a). More precisely:
  
  $$
  p(\mathbf{g}; f) \equiv p(\mathbf{g}; f, \{\phi_{ij}^{(0)}, \phi_{ij}^{(1)}, \phi_{ij}^{(2)}\}_{(i,j) \in \mathcal{F}_3^{(k)}})
  $$
  
  This is computed using the Chow-Liu tree factorization (39b) with Hardy-Weinberg marginals and IBD-adjusted pairwise probabilities (39a). For each configuration $\mathbf{g}$, $p(\mathbf{g}; f)$ gives the probability that the family's genotypes take the specific values $(g_1, \ldots, g_{n_k})$, **accounting for the genetic correlation induced by relatedness**.

**Precomputation on MAF Grid**: Since the probabilities $p(\mathbf{g}; f)$ depend on both **MAF $f$** and **family structure** $\{\phi_{ij}\}$, SPAGRM performs the following precomputation **for each family $k$ separately**:

1. **Build Chow-Liu tree** using the family-specific IBD probabilities $\{\phi_{ij}\}_{(i,j) \in \mathcal{F}_3^{(k)}}$ (this tree structure is fixed for the family)

2. **Compute probability arrays on a MAF grid**:
   $$
   \mathcal{M} = \{M_1, M_2, \ldots, M_{11}\} = \{0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5\} \tag{41}
   $$

   **Alternative geometric grid**: A more regular choice uses ratio 3:
   $$
   \mathcal{M}_{\text{alt}} = \{0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5\}
   $$

   **Comparison**:
   - **Current grid (41)**: 11 points, quasi-logarithmic with varying ratios (2-5)
     - Pros: More points in rare variant region (6 points below MAF=0.01)
     - Cons: Irregular spacing, less intuitive pattern

   - **Geometric grid**: 9 points, strict ratio of 3
     - Pros: Regular spacing, easier to understand and implement, fewer points (faster precomputation)
     - Cons: Only 5 points below MAF=0.01, may need finer grid for ultra-rare variants (MAF < 0.001)

   The current grid provides better resolution for rare variants at the cost of using 2 extra points. For applications focused on common variants (MAF > 0.01), the geometric grid is preferable. For rare variant analysis, the current grid or an even denser grid at low MAF may be better.

3. **For each grid point $M_x \in \mathcal{M}$**, compute $p(\mathbf{g}; M_x)$ for **all $3^{n_k}$ configurations** using equation (39b) with $f = M_x$

**Computational Cost per Family**:

- **Total for family $k$**: $|\mathcal{M}| \times 3^{n_k}$ probability values
- **With current grid**: $11 \times 3^{n_k}$ values
- **With geometric grid**: $9 \times 3^{n_k}$ values
- **Family size $n_k = 3$**: $11 \times 27 = 297$ values (current) vs. $9 \times 27 = 243$ values (geometric)
- **Family size $n_k = 4$**: $11 \times 81 = 891$ values (current) vs. $9 \times 81 = 729$ values (geometric)
- **Family size $n_k = 5$**: $11 \times 243 = 2{,}673$ values (current) vs. $9 \times 243 = 2{,}187$ values (geometric)

**Key Points**:

- Each family $k$ has its **own set of precomputed probability arrays** because each family has a unique IBD structure
- Two families of the same size (e.g., both $n_k = 3$) but with different relatedness patterns (e.g., one with full siblings, another with half-siblings) will have **different** $p(\mathbf{g}; f)$ values
- The precomputation is done **once** per family, and all variants use the same precomputed arrays (with interpolation for intermediate MAF values)

**Linear Interpolation**: For a given variant with MAF $f \in [M_x, M_{x+1}]$, the probability array is obtained by linear interpolation:

$$
p(\mathbf{g}; f) = p(\mathbf{g}; M_x) + \frac{f - M_x}{M_{x+1} - M_x} \left[p(\mathbf{g}; M_{x+1}) - p(\mathbf{g}; M_x)\right] \tag{41a}
$$

This interpolation is applied to **all $3^{n_k}$ configurations** to obtain the full probability array $\{p(\mathbf{g}; f)\}$.

**Computational Advantage**: Without precomputation, one would need to compute $p(\mathbf{g}; f)$ from scratch for each variant by:

1. Computing pairwise IBD-adjusted genotype probabilities $\Pr(G_i, G_j | f, \rho_{ij})$ for all tree edges
2. Applying the tree factorization (39) to get $p(\mathbf{g}; f)$ for all $3^{n_k}$ configurations

With the MAF grid and interpolation, the costly Chow-Liu tree probability computation is done only 11 times (once per grid point), and subsequent variants use fast linear interpolation.

#### 4.4.3 MGF for Three-or-More-Subject Families

**Relationship between PMF and CGF**: The MGF and CGF for family $k$ are constructed from the PMF $p(\mathbf{g}; f)$ as follows.

The MGF for family $k$ is the expectation of $e^{tS_k(\mathbf{G}_k)}$ over the joint distribution:

$$
M_{\mathcal{F}_3^{(k)}}(t; f) = \mathbb{E}[e^{tS_k(\mathbf{G}_k)}] = \sum_{\mathbf{g}} e^{tS_k(\mathbf{g})} p(\mathbf{g}; f) \tag{42}
$$

where the sum is over all $3^{n_k}$ possible configurations $\mathbf{g} = (g_1, \ldots, g_{n_k})$.

The first and second derivatives are

$$
M_{\mathcal{F}_3^{(k)\prime}}(t; f) = \sum_{\mathbf{g}} S_k(\mathbf{g}) \cdot e^{tS_k(\mathbf{g})} p(\mathbf{g}; f) \tag{43}
$$

$$
M_{\mathcal{F}_3^{(k)\prime\prime}}(t; f) = \sum_{\mathbf{g}} S_k(\mathbf{g})^2 \cdot e^{tS_k(\mathbf{g})} p(\mathbf{g}; f) \tag{44}
$$

The CGF for family $k$ is the logarithm of the MGF:

$$
K_{\mathcal{F}_3^{(k)}}(t; f) = \ln M_{\mathcal{F}_3^{(k)}}(t; f) = \ln\left[\sum_{\mathbf{g}} e^{tS_k(\mathbf{g})} p(\mathbf{g}; f)\right] \tag{45}
$$

**Explicit functional relationship**: Given the PMF $\{p(\mathbf{g}; f)\}_{\mathbf{g} \in \{0,1,2\}^{n_k}}$, the CGF is

$$
K_{\mathcal{F}_3^{(k)}}(t; f) = \ln\left[\sum_{g_1=0}^{2} \cdots \sum_{g_{n_k}=0}^{2} \exp\left(t \sum_{i=1}^{n_k} R_i g_i\right) p((g_1,\ldots,g_{n_k}); f)\right]
$$

This makes explicit that the CGF depends on:

- The residuals $\{R_i\}_{i \in \mathcal{F}_3^{(k)}}$ (via the score function)
- The PMF $p(\mathbf{g}; f)$ (via the Chow-Liu tree approximation)
- The saddlepoint parameter $t$
- The MAF $f$

with derivatives

$$
K_{\mathcal{F}_3^{(k)\prime}}(t; f) = \frac{M_{\mathcal{F}_3^{(k)\prime}}(t; f)}{M_{\mathcal{F}_3^{(k)}}(t; f)} = \frac{\sum_{\mathbf{g}} S_k(\mathbf{g}) \cdot e^{tS_k(\mathbf{g})} p(\mathbf{g}; f)}{\sum_{\mathbf{g}} e^{tS_k(\mathbf{g})} p(\mathbf{g}; f)} \tag{46}
$$

$$
K_{\mathcal{F}_3^{(k)\prime\prime}}(t; f) = \frac{M_{\mathcal{F}_3^{(k)}}(t; f) M_{\mathcal{F}_3^{(k)\prime\prime}}(t; f) - [M_{\mathcal{F}_3^{(k)\prime}}(t; f)]^2}{[M_{\mathcal{F}_3^{(k)}}(t; f)]^2} \tag{47}
$$

## 5. Non-Outlier Normal Approximation

For non-outliers $\mathcal{N}$, the score contribution $S_{\mathcal{N}} = \sum_{i \in \mathcal{N}} R_i G_i$ is approximated as normally distributed:

$$
S_{\mathcal{N}} \sim N(\mu_{\mathcal{N}}, \sigma_{\mathcal{N}}^2) \tag{48}
$$

where

$$
\mu_{\mathcal{N}} = 2f \sum_{i \in \mathcal{N}} R_i \tag{49}
$$

$$
\sigma_{\mathcal{N}}^2 = 2f(1-f) \cdot \mathbf{R}_{\mathcal{N}}^\top \mathbf{\Phi}_{\mathcal{N}} \mathbf{R}_{\mathcal{N}} \tag{50}
$$

The MGF of the normal distribution is

$$
M_{\mathcal{N}}(t) = e^{\mu_{\mathcal{N}} t + \frac{1}{2}\sigma_{\mathcal{N}}^2 t^2} \tag{51}
$$

with derivatives

$$
M_{\mathcal{N}}'(t) = (\mu_{\mathcal{N}} + \sigma_{\mathcal{N}}^2 t) M_{\mathcal{N}}(t) \tag{52}
$$

$$
M_{\mathcal{N}}''(t) = \sigma_{\mathcal{N}}^2 M_{\mathcal{N}}(t) + (\mu_{\mathcal{N}} + \sigma_{\mathcal{N}}^2 t)^2 M_{\mathcal{N}}(t) \tag{53}
$$

The CGF is

$$
K_{\mathcal{N}}(t) = \ln M_{\mathcal{N}}(t) = \mu_{\mathcal{N}} t + \frac{1}{2}\sigma_{\mathcal{N}}^2 t^2 \tag{54}
$$

with derivatives

$$
K_{\mathcal{N}}'(t) = \mu_{\mathcal{N}} + \sigma_{\mathcal{N}}^2 t, \quad K_{\mathcal{N}}''(t) = \sigma_{\mathcal{N}}^2 \tag{55}
$$

## 6. Total Cumulant Generating Function (CGF)

The total MGF is the product of MGFs from all components:

$$
M(t) = M_{\mathcal{N}}(t) \prod_{i \in \mathcal{U}} M_G^{(i)}(t) \prod_{k} M_{\mathcal{F}_2^{(k)}}(t) \prod_{k} M_{\mathcal{F}_3^{(k)}}(t) \tag{56}
$$

The CGF is $K(t) = \ln M(t)$. Taking the logarithm:

$$
K(t) = K_{\mathcal{N}}(t) + \sum_{i \in \mathcal{U}} K_G^{(i)}(t) + \sum_{k} K_{\mathcal{F}_2^{(k)}}(t) + \sum_{k} K_{\mathcal{F}_3^{(k)}}(t) \tag{57}
$$

Expanding using (54), (17), (36), and (45):

$$
K(t) = \mu_{\mathcal{N}} t + \frac{1}{2}\sigma_{\mathcal{N}}^2 t^2 + \sum_{i \in \mathcal{U}} 2\ln(1 - f + fe^{tR_i}) + \sum_{k} \ln M_{\mathcal{F}_2^{(k)}}(t) + \sum_{k} \ln M_{\mathcal{F}_3^{(k)}}(t) \tag{58}
$$

The first derivative is

$$
K'(t) = \mu_{\mathcal{N}} + \sigma_{\mathcal{N}}^2 t + \sum_{i \in \mathcal{U}} \frac{2fR_ie^{tR_i}}{1 - f + fe^{tR_i}} + \sum_{k} \frac{M_{\mathcal{F}_2^{(k)\prime}}(t)}{M_{\mathcal{F}_2^{(k)}}(t)} + \sum_{k} \frac{M_{\mathcal{F}_3^{(k)\prime}}(t)}{M_{\mathcal{F}_3^{(k)}}(t)} \tag{59}
$$

The second derivative is

$$
K''(t) = \sigma_{\mathcal{N}}^2 + \sum_{i \in \mathcal{U}} \frac{2fR_i^2e^{tR_i}(1-f)}{(1 - f + fe^{tR_i})^2} + \sum_{k} K_{\mathcal{F}_2^{(k)}}''(t) + \sum_{k} K_{\mathcal{F}_3^{(k)}}''(t) \tag{60}
$$

## 7. Variance Ratio Correction

The variance computed from the CGF (60) **assumes diagonal covariance for outliers** (ignoring GRM correlation among outliers in the same family). To correct for this, SPAGRM computes an **empirical variance** that accounts for the full GRM:

$$
\mathbb{V}_{\text{emp}} = \sigma_G^2 \cdot (\mathbf{R}_{\mathcal{N}\cup\mathcal{U}}^\top \mathbf{\Phi}_{\mathcal{N}\cup\mathcal{U}} \mathbf{R}_{\mathcal{N}\cup\mathcal{U}} + \text{Two-Fam Contribution}) + \mathbb{V}_{\text{Three-Fam}} \tag{61}
$$

where:

- $\mathbf{R}_{\mathcal{N}\cup\mathcal{U}}$ are residuals for non-outliers and unrelated outliers
- $\text{Two-Fam Contribution} = \sum_k \mathbf{R}_{\mathcal{F}_2^{(k)}}^\top \mathbf{\Phi}_{\mathcal{F}_2^{(k)}} \mathbf{R}_{\mathcal{F}_2^{(k)}}$
- $\mathbb{V}_{\text{Three-Fam}}$ is computed from the Chow-Liu tree probabilities:

  $$
  \mathbb{V}_{\text{Three-Fam}} = \sum_k \left[\sum_{\mathbf{g}} S_k(\mathbf{g})^2 p(\mathbf{g}; f) - \left(\sum_{\mathbf{g}} S_k(\mathbf{g}) p(\mathbf{g}; f)\right)^2\right] \tag{62}
  $$

The **variance ratio** (empirical to theoretical) is

$$
\rho = \frac{\mathbb{V}_{\text{emp}}}{\mathbb{V}(S)} \tag{63}
$$

where $\mathbb{V}(S) = \sigma_G^2 \cdot \mathbf{R}^\top \mathbf{\Phi} \mathbf{R}$ is the full GRM-based variance from equation (8).

**Variance ratio scaling**: The score is rescaled before applying SPA to ensure the CGF-based variance matches the GRM-based variance:

$$
\tilde{S} = S \cdot \sqrt{\rho} \tag{64}
$$

**Note on notation**: The scaled score $\tilde{S}$ in (64) is **distinct** from the centered score $S$ in (2):

- Equation (2): $S = \mathbf{R}^\top \mathbf{G} - \overline{G} \sum R_i$ **centers** the score to have $\mathbb{E}[S] = 0$
- Equation (64): $\tilde{S} = S \cdot \sqrt{\rho}$ **scales** the centered score to align CGF variance with GRM variance

Both transformations are applied sequentially: first centering, then scaling.

## 8. Saddlepoint Approximation

The saddlepoint approximation provides an accurate tail probability for the score statistic. Given the scaled score $\tilde{S}$ from (64) and CGF $K(t)$ from (57), the saddlepoint $\hat{t}$ solves

$$
K'(\hat{t}) = \tilde{S} \tag{65}
$$

This is solved using Newton-Raphson iteration:

$$
t_{n+1} = t_n - \frac{K'(t_n) - \tilde{S}}{K''(t_n)} \tag{66}
$$

The Lugannani-Rice formula gives the tail probability:

$$
\Pr(S \geq \tilde{S}) \approx \Phi(u) \tag{67}
$$

where

$$
w = \text{sign}(\hat{t}) \sqrt{2(\hat{t} \tilde{S} - K(\hat{t}))} \tag{68}
$$

$$
v = \hat{t} \sqrt{K''(\hat{t})} \tag{69}
$$

$$
u = w + \frac{1}{w} \ln\left(\frac{v}{w}\right) \tag{70}
$$

and $\Phi$ is the standard normal CDF.

### 8.1 Two-Sided P-Value

For a two-sided test, compute both tails:

$$
\text{p-value} = \Pr(S \geq |\tilde{S}|) + \Pr(S \leq -|\tilde{S}|) \tag{71}
$$

Each tail uses a separate saddlepoint with initial values:

- **Right tail** ($\tilde{S} > 0$): $t_0 = \min(|\tilde{S}|/\mathbb{V}(S), 1.2)$
- **Left tail** ($\tilde{S} < 0$): $t_0 = 0$ (default `zeta` parameter)

## 9. Normal Approximation Cutoff

If the standardized score is small, the normal approximation is sufficiently accurate and SPA is not needed:

$$
|Z| = \left|\frac{S}{\sqrt{\mathbb{V}(S)}}\right| \leq c_{\text{SPA}} \implies \text{p-value} = 2\Phi(-|Z|) \tag{72}
$$

where $c_{\text{SPA}}$ is the SPA cutoff (default 2.0).

## Summary

SPAGRM handles sample relatedness by:

1. **Score centering** (Eq. 2): Compute $S = \mathbf{R}^\top \mathbf{G} - \overline{G} \sum R_i$ to ensure $\mathbb{E}[S] = 0$ regardless of $\sum R_i$

2. **Partitioning** the sample into:
   - Non-outliers $\mathcal{N}$: approximated by normal distribution
   - Unrelated outliers $\mathcal{U}$: independent genotypes
   - Two-subject families $\mathcal{F}_2^{(k)}$: IBD-based joint distribution
   - Three-or-more-subject families $\mathcal{F}_3^{(k)}$: Chow-Liu tree approximation

3. **CGF construction** (Eq. 57-60):
   - For non-outliers: normal CGF $K_{\mathcal{N}}(t) = \mu_{\mathcal{N}} t + \frac{1}{2}\sigma_{\mathcal{N}}^2 t^2$
   - For unrelated outliers: sum of independent CGFs $\sum_i K_G^{(i)}(t)$
   - For two-subject families: IBD-adjusted CGF $K_{\mathcal{F}_2}(t)$
   - For three-or-more families: precomputed PMF-based CGF $K_{\mathcal{F}_3}(t)$

4. **Variance ratio correction** (Eq. 64): Scale score $\tilde{S} = S \sqrt{\rho}$ to align CGF-based variance with full GRM-based variance

5. **Saddlepoint approximation** (Eq. 65-70): Solve $K'(\hat{t}) = \tilde{S}$ and compute tail probability via Lugannani-Rice formula

This approach achieves computational efficiency (avoiding full $3^n$ enumeration) while maintaining accuracy through the tree-structured factorization and empirical variance adjustment. The method correctly handles:

- **Uncentered genotype MGFs** (Eq. 13-47): MGF derived for $G_i \sim \text{Binomial}(2,f)$
- **Centered score matching** (Eq. 2, 65): Saddlepoint equation matches centered observed score
- **Two-stage adjustment**: First centering (Eq. 2) for mean, then scaling (Eq. 64) for variance
