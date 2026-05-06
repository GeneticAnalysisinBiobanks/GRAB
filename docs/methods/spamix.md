# SPAmix / SPAmixPlus

SPAmix (Saddlepoint Approximation for Mixed models) performs score tests for genetic variants accounting for individual-specific allele frequencies and population structure. SPAmixPlus extends SPAmix by incorporating a sparse genetic relationship matrix (GRM) to account for cryptic relatedness.

Let:

- $n \times 1$ vector $\mathbf{G}$ represent genotypes (allele counts) for a variant to be tested;
- $n \times p$ matrix $\mathbf{X}$ represent $p$ principal components (PCs) for ancestry adjustment;
- $n \times 1$ vector $\mathbf{R} = (R_1, \ldots, R_n)^\top$ represent residuals from the null model;
- $n \times 1$ vector $\mathbf{f} = (f_1, \ldots, f_n)^\top$ represent individual-specific allele frequencies;
- $n \times n$ sparse matrix $\mathbf{\Phi}$ represent the GRM (SPAmixPlus only);
- $\beta$ denote the genetic effect to be tested.

The null hypothesis is $H_0: \beta = 0$.

## 1. Score Statistic

The score statistic for testing $\beta$ is

$$
S = \sum_{i=1}^n R_i G_i = \mathbf{R}^\top \mathbf{G} \tag{1}
$$

Under the null hypothesis, assuming Hardy-Weinberg equilibrium, the expected genotype is $\mathbb{E}[G_i | f_i] = 2f_i$. The expected score is

$$
\mathbb{E}[S] = \sum_{i=1}^n R_i \mathbb{E}[G_i] = 2 \sum_{i=1}^n R_i f_i = 2 \mathbf{R}^\top \mathbf{f} \tag{2}
$$

where $f_i$ is the estimated allele frequency for individual $i$. This formula applies to both SPAmix and SPAmixPlus; it does not assume independence between individuals (which is addressed in the variance calculation).

## 2. Individual-Specific Allele Frequency Estimation

Individual-specific allele frequencies are estimated to account for population structure. For each marker, we estimate $f_i$ using one of three methods:

### 2.1 Method Selection Cascade

Let $\text{MAC} = 2n \cdot \overline{f}$ be the minor allele count, where $\overline{f}$ is the overall allele frequency. Let $\tilde{\mathbf{X}} = [\mathbf{1}_n, \mathbf{X}]$ be the design matrix with intercept and PCs.

**Method 0 (Uniform)**: If $\text{MAC} \leq 20$, set

$$
f_i = \overline{f}, \quad \forall i \tag{3}
$$

**Method 1 (Linear Regression)**: If $\text{MAC} > 20$, fit

$$
G_i = \beta_0 + \sum_{j=1}^p \beta_j X_{ij} + \epsilon_i \tag{4}
$$

and compute fitted allele frequency as

$$
\hat{f}_i = \frac{1}{2}\left(\hat{\beta}_0 + \sum_{j=1}^p \hat{\beta}_j X_{ij}\right) \tag{5}
$$

If the proportion of fitted values $\hat{f}_i$ outside $[0, 1]$ is less than 10%, use Method 1. Otherwise, proceed to Method 2.

**Method 2 (Logistic Regression)**: Let $\mathbf{X}_{\text{sig}}$ be the subset of PCs with p-values $< 0.05$ from the linear model (4). Binarize genotypes as $G_i^* = \mathbb{1}(G_i > 0.5)$ and fit

$$
\text{logit}(\mu_i) = \gamma_0 + \sum_{j \in \text{sig}} \gamma_j X_{ij} \tag{6}
$$

The estimated allele frequency is

$$
\hat{f}_i = 1 - \sqrt{1 - \mu_i} \tag{7}
$$

where $\mu_i = 1/(1 + e^{-\eta_i})$ and $\eta_i = \gamma_0 + \sum_{j \in \text{sig}} \gamma_j X_{ij}$.

If no PCs are significant or if the binarized MAC $\leq 20$, fall back to Method 0.

## 3. Score Variance

### 3.1 Diagonal Variance (SPAmix)

When no GRM is provided, the variance is computed as

$$
\widehat{\mathbb{V}}_{\text{diag}}(S) = \sum_{i=1}^n R_i^2 \cdot 2 f_i (1 - f_i) \tag{8}
$$

This assumes independence between individuals conditional on covariates.

### 3.2 GRM-Based Variance (SPAmixPlus)

When a sparse GRM $\mathbf{\Phi}$ is provided, the variance accounts for correlation between individuals' genotypes:

$$
\widehat{\mathbb{V}}_{\text{GRM}}(S) = \mathbb{V}\left(\sum_{i=1}^n R_i G_i\right) = \sum_{i=1}^n \sum_{j=1}^n R_i R_j \text{Cov}(G_i, G_j) \tag{9}
$$

The GRM $\mathbf{\Phi}$ represents the correlation structure of genotypes. Under Hardy-Weinberg equilibrium, the genotype covariance is modeled as:

$$
\text{Cov}(G_i, G_j) = 2\sqrt{f_i(1-f_i) \cdot f_j(1-f_j)} \cdot \Phi_{ij} \tag{10}
$$

where $\Phi_{ii} \approx 1$ (diagonal elements). Substituting (10) into (9):

$$
\begin{align}
\widehat{\mathbb{V}}_{\text{GRM}}(S) &= \sum_{i=1}^n \sum_{j=1}^n R_i R_j \cdot 2\sqrt{f_i(1-f_i) \cdot f_j(1-f_j)} \cdot \Phi_{ij} \tag{11} \\
&= \sum_{i=1}^n \sum_{j=1}^n \left[R_i \sqrt{2f_i(1-f_i)}\right] \Phi_{ij} \left[R_j \sqrt{2f_j(1-f_j)}\right] \tag{12}
\end{align}
$$

This is a quadratic form $\mathbf{w}^\top \mathbf{\Phi} \mathbf{w}$ where $w_i = R_i \sqrt{2f_i(1-f_i)}$. Separating diagonal and off-diagonal contributions:

$$
\widehat{\mathbb{V}}_{\text{GRM}}(S) = \sum_{i=1}^n \Phi_{ii} w_i^2 + 2\sum_{i<j} \Phi_{ij} w_i w_j \tag{13}
$$

In vector notation, define $\mathbf{w} = (w_1, \ldots, w_n)^\top$ where $w_i = R_i\sqrt{2f_i(1-f_i)}$. Then:

$$
\widehat{\mathbb{V}}_{\text{GRM}}(S) = \mathbf{w}^\top \mathbf{\Phi} \mathbf{w} = \mathbf{w}^\top \text{diag}(\mathbf{\Phi}) \mathbf{w} + 2\sum_{i<j} \Phi_{ij} w_i w_j \tag{13a}
$$

In implementation, the sparse GRM is stored in lower-triangular format (including diagonal). The computation iterates over stored entries: diagonal terms contribute once ($\Phi_{ii} w_i^2$), off-diagonal terms contribute twice ($2\Phi_{ij} w_i w_j$ for $i \neq j$).

## 4. Saddlepoint Approximation with Genotype-Dependent CGF

### 4.1 Empirical CGF of Genotype

Consider a single genotype $G_i$ following Hardy-Weinberg equilibrium with allele frequency $f_i$. The genotype follows a binomial distribution:

$$
G_i \sim \text{Binomial}(2, f_i) \tag{14}
$$

with probability mass function:

$$
\Pr(G_i = k) = \binom{2}{k} f_i^k (1-f_i)^{2-k}, \quad k \in \{0, 1, 2\} \tag{15}
$$

The moment generating function (MGF) of a single genotype is:

$$
\begin{align}
M_{G_i}(t; f_i) &= \mathbb{E}[e^{t G_i}] = \sum_{k=0}^2 e^{tk} \Pr(G_i = k) \tag{16} \\
&= (1-f_i)^2 \cdot e^{0} + 2f_i(1-f_i) \cdot e^{t} + f_i^2 \cdot e^{2t} \tag{17} \\
&= (1-f_i)^2 + 2f_i(1-f_i)e^t + f_i^2 e^{2t} \tag{18} \\
&= [(1-f_i) + f_i e^t]^2 \tag{19}
\end{align}
$$

For notational simplicity, we write $M_G(t; f) = (1 - f + f e^t)^2$ for a genotype with allele frequency $f$:

$$
M_G(t; f) = (1 - f + f e^t)^2 \tag{20}
$$

The first and second derivatives are

$$
M_G'(t; f) = 2 f e^t (1 - f + f e^t) \tag{21}
$$

$$
M_G''(t; f) = 2 (f e^t)^2 + 2 f e^t (1 - f + f e^t) = 2f e^t (1 - f + 2f e^t) \tag{22}
$$

The cumulant generating function (CGF) is $K_G(t; f) = \ln M_G(t; f)$. Expanding as a function of $t$ and $f$:

$$
K_G(t; f) = \ln[(1 - f + f e^t)^2] = 2\ln(1 - f + f e^t) \tag{22a}
$$

The first derivative is:

$$
K_G'(t; f) = \frac{M_G'(t; f)}{M_G(t; f)} = \frac{2 f e^t (1 - f + f e^t)}{(1 - f + f e^t)^2} = \frac{2 f e^t}{1 - f + f e^t} \tag{23}
$$

The second derivative, expanded as a function of $t$ and $f$:

$$
\begin{align}
K_G''(t; f) &= \frac{M_G(t; f) M_G''(t; f) - [M_G'(t; f)]^2}{[M_G(t; f)]^2} \tag{24} \\
&= \frac{(1-f+fe^t)^2 \cdot 2fe^t(1-f+2fe^t) - [2fe^t(1-f+fe^t)]^2}{(1-f+fe^t)^4} \tag{24a} \\
&= \frac{2fe^t(1-f+2fe^t) - 4f^2e^{2t}}{(1-f+fe^t)^2} \tag{24b} \\
&= \frac{2fe^t(1-f)}{(1-f+fe^t)^2} \tag{24c}
\end{align}
$$

### 4.2 Outlier-Based Partitioning

To improve numerical stability, the sample is partitioned into **outliers** (extreme residuals) and **non-outliers** (central residuals). Let

- $\mathcal{O}$ be the set of outlier indices (residuals outside $[Q_{25} - r \cdot \text{IQR}, Q_{75} + r \cdot \text{IQR}]$ with $r = 1.5$)
- $\mathcal{N}$ be the set of non-outlier indices

The total score can be decomposed as

$$
S = \sum_{i \in \mathcal{O}} R_i G_i + \sum_{i \in \mathcal{N}} R_i G_i = S_{\mathcal{O}} + S_{\mathcal{N}} \tag{25}
$$

### 4.3 Outlier CGF

For outliers, the exact CGF is computed assuming **diagonal variance** (independent genotypes):

$$
K_{\mathcal{O}}(t) = \sum_{i \in \mathcal{O}} K_G(t R_i; f_i) \tag{26}
$$

$$
K_{\mathcal{O}}'(t) = \sum_{i \in \mathcal{O}} R_i K_G'(t R_i; f_i) \tag{27}
$$

$$
K_{\mathcal{O}}''(t) = \sum_{i \in \mathcal{O}} R_i^2 K_G''(t R_i; f_i) \tag{28}
$$

**Applicability**:

- **SPAmix** (no GRM): These formulas are directly used in the saddlepoint approximation, as the diagonal variance assumption holds.
- **SPAmixPlus** (with GRM): The GRM introduces correlation between genotypes. To use these diagonal-variance-based CGF formulas, a variance ratio correction $\rho$ is applied (Section 5.1).

### 4.4 Non-Outlier Normal Approximation

For non-outliers, the exact CGF would require summing over all individuals in $\mathcal{N}$ (typically thousands). Instead, we use a **normal approximation** based on the Central Limit Theorem: when $|\mathcal{N}|$ is large and individual contributions are small, $S_{\mathcal{N}} = \sum_{i \in \mathcal{N}} R_i G_i$ is approximately normal.

The mean and variance of $S_{\mathcal{N}}$ under the **diagonal variance assumption** (independent genotypes) are:

$$
\mu_{\mathcal{N}} = \mathbb{E}(S_{\mathcal{N}}) = 2 \sum_{i \in \mathcal{N}} R_i f_i \tag{29}
$$

$$
\sigma_{\mathcal{N}}^2 = \mathbb{V}_{\text{diag}}(S_{\mathcal{N}}) = 2 \sum_{i \in \mathcal{N}} R_i^2 f_i (1 - f_i) \tag{30}
$$

The CGF of the normal distribution with mean $\mu_{\mathcal{N}}$ and variance $\sigma_{\mathcal{N}}^2$ is

$$
K_{\mathcal{N}}(t) = \mu_{\mathcal{N}} t + \frac{1}{2} \sigma_{\mathcal{N}}^2 t^2 \tag{31}
$$

with derivatives

$$
K_{\mathcal{N}}'(t) = \mu_{\mathcal{N}} + \sigma_{\mathcal{N}}^2 t, \quad K_{\mathcal{N}}''(t) = \sigma_{\mathcal{N}}^2 \tag{32}
$$

**Applicability**:

- **SPAmix** (no GRM): The diagonal variance assumption $\sigma_{\mathcal{N}}^2$ in (30) is exact, and the CGF (31) is directly used.
- **SPAmixPlus** (with GRM): The true variance of $S_{\mathcal{N}}$ should account for genotype correlation via the GRM (equations 9-13). However, the CGF formulas (26-32) retain the diagonal form. To reconcile this, the entire score $S$ is rescaled by $\sqrt{\rho}$ (Section 5.1), effectively transforming the problem to the diagonal variance scale where these CGF formulas are valid.

### 4.5 Total CGF

The total CGF is

$$
K(t) = K_{\mathcal{O}}(t) + K_{\mathcal{N}}(t) \tag{33}
$$

$$
K'(t) = K_{\mathcal{O}}'(t) + \mu_{\mathcal{N}} + \sigma_{\mathcal{N}}^2 t \tag{34}
$$

$$
K''(t) = K_{\mathcal{O}}''(t) + \sigma_{\mathcal{N}}^2 \tag{35}
$$

## 5. Lugannani-Rice Saddlepoint Formula

### 5.1 Variance Ratio Correction (SPAmixPlus Only)

When using a GRM, the SPA operates on a variance-ratio-corrected score. Let $q$ denote the observed value of $S$. Define the variance ratio using equations (8) and (9):

$$
\rho = \frac{\widehat{\mathbb{V}}_{\text{diag}}(S)}{\widehat{\mathbb{V}}_{\text{GRM}}(S)} \tag{36}
$$

The corrected score when observing $S=q$ and its expectation are:

$$
q\sqrt{\rho} \quad \text{and} \quad \mathbb{E}[S]\sqrt{\rho} \tag{37}
$$

This correction ensures that the SPA, which assumes the diagonal variance structure in the CGF formulas, provides accurate tail probabilities.

### 5.2 Saddlepoint Equation

For the corrected observed score $q\sqrt{\rho}$, solve for the saddlepoint $\zeta$ such that

$$
K'(\zeta) = q\sqrt{\rho} \tag{38}
$$

using Newton-Raphson iteration:

$$
\zeta_{n+1} = \zeta_n - \frac{K'(\zeta_n) - q\sqrt{\rho}}{K''(\zeta_n)} \tag{39}
$$

with convergence tolerance $|\Delta \zeta| < 0.001$ and maximum 100 iterations.

### 5.3 Lugannani-Rice Formula

Compute

$$
w = \text{sgn}(\zeta) \sqrt{2(\zeta \cdot q\sqrt{\rho} - K(\zeta))} \tag{40}
$$

$$
v = \zeta \sqrt{K''(\zeta)} \tag{41}
$$

The tail probability is

$$
\Pr(S\sqrt{\rho} < q\sqrt{\rho}) \simeq \Phi\left(w + \frac{1}{w} \ln\left(\frac{v}{w}\right)\right) \tag{42}
$$

where $\Phi$ is the standard normal CDF.

### 5.4 Two-Sided P-Value

For two-sided testing, compute

$$
\delta = |q\sqrt{\rho} - \mathbb{E}[S]\sqrt{\rho}| = |q - \mathbb{E}[S]|\sqrt{\rho} \tag{43}
$$

$$
p = \Pr(S\sqrt{\rho} > \mathbb{E}[S]\sqrt{\rho} + \delta) + \Pr(S\sqrt{\rho} < \mathbb{E}[S]\sqrt{\rho} - \delta) \tag{44}
$$

Each tail probability is computed using (42) with the saddlepoint equation (38) solved for $\mathbb{E}[S]\sqrt{\rho} \pm \delta$.

### 5.5 Normal Approximation for Small Z-Scores

If the z-score

$$
Z = \frac{q - \mathbb{E}[S]}{\sqrt{\widehat{\mathbb{V}}(S)}} \tag{45}
$$

satisfies $|Z| <$ `spaCutoff` (typically 2), use the normal approximation:

$$
p = 2 \Phi(-|Z|) \tag{46}
$$

where $\widehat{\mathbb{V}}(S)$ is either $\widehat{\mathbb{V}}_{\text{diag}}(S)$ from (8) for SPAmix or $\widehat{\mathbb{V}}_{\text{GRM}}(S)$ from (9) for SPAmixPlus.

This avoids numerical instability in the SPA for statistics close to the mean.
