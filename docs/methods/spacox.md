# SPACox

SPACox (Saddlepoint Approximation for Cox Model) performs score tests for genetic variants in survival analysis.

Let:

- $n \times 1$ vector $\mathbf{G}$ represent genotypes (allele counts) for a variant to be tested;
- $n \times p$ matrix $\mathbf{X}$ represent $p$ covariates (excluding intercept);
- $n \times 1$ vector $\mathbf{R} = (R_1, \ldots, R_n)^\top$ represent martingale residuals from the null Cox model;
- $\beta$ denote the genetic effect to be tested.

The null hypothesis is $H_0: \beta = 0$.

## 1. Score Statistic and Projection

The score statistic for testing $\beta$ is

$$
S = \sum_{i=1}^n R_i G_i = \mathbf{R}^\top \mathbf{G} \tag{1}
$$

Martingale residuals satisfy linear restrictions $\sum_{i=1}^n X_i R_i = 0$ and $\sum_{i=1}^n R_i = 0$. Let $\tilde{\mathbf{X}} = [\mathbf{1}_n, \mathbf{X}]$ be the design matrix including intercept. The restrictions are $\tilde{\mathbf{X}}^\top \mathbf{R} = \mathbf{0}$. Define the projection matrix

$$
\mathbf{Q} = \mathbf{I}_n - \tilde{\mathbf{X}}(\tilde{\mathbf{X}}^\top \tilde{\mathbf{X}})^{-1} \tilde{\mathbf{X}}^\top \tag{2}
$$

We assume $\mathbf{R} = \mathbf{Q}\tilde{\mathbf{R}}$, where $\tilde{\mathbf{R}}$ is a latent random vector without linear restrictions. The score can be rewritten as

$$
S = \mathbf{R}^\top \mathbf{G} = \tilde{\mathbf{R}}^\top \mathbf{Q}\mathbf{G} = \tilde{\mathbf{R}}^\top \tilde{\mathbf{G}} \tag{3}
$$

where

$$
\tilde{\mathbf{G}} = \mathbf{Q}\mathbf{G} = \mathbf{G} - \tilde{\mathbf{X}}(\tilde{\mathbf{X}}^\top \tilde{\mathbf{X}})^{-1} \tilde{\mathbf{X}}^\top \mathbf{G} \tag{4}
$$

is the **covariate-adjusted genotype vector**. Because $\mathbf{R} = \mathbf{Q}\tilde{\mathbf{R}}$, $\mathbf{R}$ is a natural representative of $\tilde{\mathbf{R}}$, and we use the observed martingale residuals $R_i, i \leq n$ to estimate the empirical distribution of $\tilde{\mathbf{R}}$.

## 2. Empirical CGF of Residuals

To construct the CGF of $S$, we first estimate the moment generating function (MGF) of $R_i$. The empirical MGF of $R_i$ is

$$
\hat{M}_0(z) = \mathbb{E}(e^{zR}) \simeq \frac{1}{n} \sum_{i=1}^n e^{z R_i} \tag{5}
$$

and its first and second derivatives are

$$
\hat{M}_0'(z) \simeq \frac{1}{n} \sum_{i=1}^n R_i \cdot e^{z R_i}, \quad \hat{M}_0''(z) \simeq \frac{1}{n} \sum_{i=1}^n R_i^2 \cdot e^{z R_i} \tag{6}
$$

The empirical CGF of $R_i$ is $\hat{K}_0(z) = \ln \hat{M}_0(z)$. Substituting (5) into the CGF definition:

$$
\hat{K}_0(z) = \ln\left[\frac{1}{n}\sum_{i=1}^n e^{z R_i}\right] \tag{7}
$$

The first derivative of the CGF is:

$$
\hat{K}_0'(z) = \frac{\hat{M}_0'(z)}{\hat{M}_0(z)} = \frac{\sum_{i=1}^n R_i e^{z R_i}}{\sum_{i=1}^n e^{z R_i}} \tag{8}
$$

The second derivative is:

$$
\hat{K}_0''(z) = \frac{\hat{M}_0''(z)\hat{M}_0(z) - [\hat{M}_0'(z)]^2}{[\hat{M}_0(z)]^2} = \frac{\left(\sum_{i=1}^n e^{z R_i}\right)\left(\sum_{i=1}^n R_i^2 e^{z R_i}\right) - \left(\sum_{i=1}^n R_i e^{z R_i}\right)^2}{\left(\sum_{i=1}^n e^{z R_i}\right)^2} \tag{9}
$$

## 3. Score Variance and CGF

Considering $G_i, i \leq n$ as constant coefficients, we obtain the empirical variance of the score statistic $S = \sum_{i=1}^n G_i R_i$ as

$$
\hat{\mathbb{V}}(S) = \sum_{i=1}^n \tilde{G}_i^2 \cdot \hat{M}_0''(0) \tag{10}
$$

The estimated CGF of $S$ is obtained by substituting $\tilde{G}_i$ as coefficients into (7):

$$
\hat{K}(z) = \sum_{i=1}^n \hat{K}_0(\tilde{G}_i z) \tag{11}
$$

where each term $\hat{K}_0(\tilde{G}_i z)$ is evaluated using (7). In practice, these values are obtained via interpolation from precomputed table $y_k^{(0)}$ (see Section 5).

The first derivative is obtained by applying (8):

$$
\hat{K}'(z) = \sum_{i=1}^n \tilde{G}_i \hat{K}_0'(\tilde{G}_i z) \tag{12}
$$

where each $\hat{K}_0'(\tilde{G}_i z)$ is interpolated from $y_k^{(1)}$.

The second derivative is obtained by applying (9):

$$
\hat{K}''(z) = \sum_{i=1}^n \tilde{G}_i^2 \hat{K}_0''(\tilde{G}_i z) \tag{13}
$$

where each $\hat{K}_0''(\tilde{G}_i z)$ is interpolated from $y_k^{(2)}$.

## 4. Saddlepoint Approximation

Given an observed score $S = s$, we first calculate the saddlepoint $\zeta$ such that $\hat{K}'(\zeta) = s$. This is solved by Newton-Raphson iteration:

$$
\zeta_{n+1} = \zeta_n - \frac{\hat{K}'(\zeta_n) - s}{\hat{K}''(\zeta_n)} \tag{14}
$$

Then we calculate

$$
w = \text{sgn}(\zeta) \sqrt{2(\zeta s - \hat{K}(\zeta))}, \quad v = \zeta \sqrt{\hat{K}''(\zeta)} \tag{15}
$$

According to the **Lugannani-Rice formula** (Barndorff-Nielsen), the tail probability is

$$
\Pr(S < s) \simeq \Phi\left(w + \frac{1}{w} \cdot \log\left(\frac{v}{w}\right)\right) \tag{16}
$$

where $\Phi$ is the standard normal CDF.

For two-sided testing, we compute tail probabilities at $s = |S|$ and $s = -|S|$:

$$
p = \Pr(|S| > |s|) = \Pr(S > |s|) + \Pr(S < -|s|) \tag{17}
$$

## 5. Interpolation for Efficient CGF Evaluation

Direct evaluation of $\hat{K}_0(t)$ and its derivatives (equations (7)-(9)) requires $O(n)$ exponential operations per call. Since Newton-Raphson iteration (14) requires multiple evaluations over $n$ individuals, the total complexity per variant becomes $O(n^2)$, which is prohibitive for large cohorts ($n \approx 10000$).

**Solution**: Precompute $\hat{K}_0(t)$, $\hat{K}_0'(t)$, $\hat{K}_0''(t)$ on a fixed grid $\{t_1, \ldots, t_L\}$, then use interpolation to approximate values at arbitrary $t = \bar{G}_i \zeta$. This reduces complexity from $O(n^2)$ to $O(n)$ per variant.

### 5.1 Grid Construction

The grid points are generated using **Cauchy-quantile spacing**:

$$
t_k = c \cdot \tan\left[\pi\left(\frac{k}{L+1} - \frac{1}{2}\right)\right], \quad k = 1, \ldots, L \tag{18}
$$

**Parameters**:

- **$L = 10000$**: Number of grid points
- **$c$**: Scaling factor chosen such that the grid spans a desired range (typically $[-100, 100]$)

Specifically, since $\tan[\pi(k/(L+1) - 1/2)]$ ranges from $-\infty$ to $+\infty$, we set:
$$
c = \frac{100}{\max_{k \in \{1,\ldots,L\}} |\tan[\pi(k/(L+1) - 1/2)]|}
$$
to ensure $t_k \in [-100, 100]$ for all $k$.

**Why Cauchy distribution?**

The Cauchy distribution (also known as the Student's t-distribution with 1 degree of freedom) has CDF:
$$
F_{\text{Cauchy}}(x) = \frac{1}{\pi}\arctan(x) + \frac{1}{2}
$$

Inverting this gives the quantile function:
$$
F_{\text{Cauchy}}^{-1}(p) = \tan\left[\pi\left(p - \frac{1}{2}\right)\right]
$$

Setting $p = k/(L+1)$ for $k = 1, \ldots, L$ yields equally spaced quantiles in probability space, which translate to points that are:

- **Sparse in the center** (where $\hat{K}_0$ is smooth and linear interpolation is accurate)
- **Dense in the tails** (where $\hat{K}_0$ has high curvature and SPA is actually used for significant variants)

**Comparison with uniform spacing**:

- **Uniform $t_k \in [-100, 100]$**: Wastes grid points in the center, insufficient resolution in tails → poor accuracy for extreme p-values
- **Cauchy-quantile spacing**: Optimal for tail probability estimation (where $|\zeta|$ is large for significant associations)

**Precision and p-value range**:

The grid range $[-100, 100]$ determines the minimum achievable p-value:

- Maximum $|Z|$ before numerical overflow: $\approx 100 / \max_i |\bar{G}_i| \approx 100$ (since genotypes are standardized with $\sum_i \bar{G}_i^2 \approx 1$)
- Minimum p-value: $\approx 2\Phi(-100) \approx 10^{-2171}$ (far beyond genome-wide significance threshold $5 \times 10^{-8}$)

**Parameter tuning and accuracy**:

1. **Grid size $L$**:
   - **Larger $L$** → higher interpolation accuracy but more precomputation cost
     - $L = 10000$: Interpolation error $\approx 10^{-8}$ for tail probabilities, negligible compared to SPA approximation error ($\approx 10^{-6}$)
     - $L = 1000$: Saves 10× precomputation time, but p-value error $\approx 10^{-6}$ (acceptable for $p > 10^{-10}$, may lose precision for ultra-rare variants)
   - **Trade-off**: Precomputation is one-time per phenotype and amortized over millions of variants, so $L = 10000$ is preferred

2. **Scaling factor $c$** (grid range):
   - **Smaller $c$** (narrow range, e.g., $[-10, 10]$):
     - Better resolution in center
     - Tails truncated → Newton-Raphson fails for extreme $Z$ (significant variants)
     - Minimum p-value limited to $\approx 2\Phi(-10) \approx 10^{-23}$
   - **Larger $c$** (wide range, e.g., $[-1000, 1000]$):
     - Covers ultra-extreme cases
     - Sparser spacing in tails → lower interpolation accuracy
     - May hit numerical limits: $e^{1000 R_i}$ can overflow for large residuals
   - **Standard choice**: $c$ such that range is $[-100, 100]$ balances numerical stability and precision for GWAS (target: $p \geq 10^{-100}$)

### 5.2 Precomputation (Per-Phenotype)

**Input**: Residual vector $\mathbf{R} = (R_1, \ldots, R_n)$ from null Cox model

**Step 1**: Generate grid points $\{t_1, \ldots, t_L\}$ using (18)

**Step 2**: For each grid point $t_k$, evaluate the empirical CGF and its derivatives using (7)-(9):

$$
y_k^{(0)} = \hat{K}_0(t_k) = \ln\left[\frac{1}{n}\sum_{i=1}^n e^{t_k R_i}\right] \tag{19}
$$

$$
y_k^{(1)} = \hat{K}_0'(t_k) = \frac{\sum_{i=1}^n R_i e^{t_k R_i}}{\sum_{i=1}^n e^{t_k R_i}} \tag{20}
$$

$$
y_k^{(2)} = \hat{K}_0''(t_k) = \frac{\left(\sum_{i=1}^n e^{t_k R_i}\right)\left(\sum_{i=1}^n R_i^2 e^{t_k R_i}\right) - \left(\sum_{i=1}^n R_i e^{t_k R_i}\right)^2}{\left(\sum_{i=1}^n e^{t_k R_i}\right)^2} \tag{21}
$$

**Step 3**: Precompute piecewise linear interpolation slopes:

$$
s_k^{(j)} = \frac{y_{k+1}^{(j)} - y_k^{(j)}}{t_{k+1} - t_k}, \quad j = 0, 1, 2 \tag{22}
$$

This avoids repeated division during per-variant evaluation.

**Storage**: Store $(t_k, y_k^{(j)}, s_k^{(j)})$ for $k = 1, \ldots, L$ and $j = 0, 1, 2$ in a `CumulantTable` structure.

**Computational cost**: $O(Ln)$ exponential operations, amortized over all variants (typically millions in GWAS).

### 5.3 Per-Variant Evaluation

**Input**:

- Genotype vector $\mathbf{G}$ (or covariate-adjusted $\tilde{\mathbf{G}}$ from equation (4))
- Standardized score $Z = S / \sqrt{\hat{\mathbb{V}}(S)}$ where $\hat{\mathbb{V}}(S)$ is from (10)
- `CumulantTable` from precomputation (Section 5.2)

**Step 1**: Solve saddlepoint equation $\hat{K}'(\zeta) = Z$ using Newton-Raphson (14)

Initialize: $\zeta_0 = \text{sign}(Z) \cdot 3$

Iterate (up to 100 iterations, convergence tolerance $\epsilon = 0.001$):

1. **Evaluate $\hat{K}'(\zeta_n)$ using (12)**:

   For each individual $i = 1, \ldots, n$:
   - Compute $t_i = \bar{G}_i \zeta_n$ (where $\bar{G}_i$ is normalized genotype from (10))
   - **Find grid interval**: $t_k \leq t_i < t_{k+1}$
     - **O(1) method** (GRAB-feat-cpp): Inverse Cauchy transform
       $$
       k = \left\lfloor \left(\frac{\arctan(t_i / c)}{\pi} + \frac{1}{2}\right)(L+1) - 1 \right\rfloor \tag{23}
       $$
     - **O(log L) method** (GRAB-0.3.0): Binary search
   - **Piecewise linear interpolation**:
     $$
     \hat{K}_0'(t_i) \approx y_k^{(1)} + (t_i - t_k) s_k^{(1)} \tag{24}
     $$
   - Accumulate:
     $$
     \hat{K}'(\zeta_n) = \sum_{i=1}^n \bar{G}_i \cdot \hat{K}_0'(\bar{G}_i \zeta_n) \tag{25}
     $$

2. **Similarly compute $\hat{K}''(\zeta_n)$ using (13)** with $y_k^{(2)}$ and $s_k^{(2)}$

3. **Newton-Raphson update (14)**:
   $$
   \zeta_{n+1} = \zeta_n - \frac{\hat{K}'(\zeta_n) - Z}{\hat{K}''(\zeta_n)}
   $$

4. **Check convergence**: If $|\zeta_{n+1} - \zeta_n| < \epsilon$, break

**Step 2**: After convergence, interpolate $\hat{K}(\zeta)$ using $y_k^{(0)}$, then apply Lugannani-Rice formula (15)-(17):

1. Compute $w$ and $v$ from (15)
2. Evaluate tail probabilities using (16) for both $s = |Z|$ and $s = -|Z|$
3. Sum for two-sided p-value (17)

**Computational cost**: $O(n)$ arithmetic operations per variant (vs. $O(n^2)$ without interpolation)

**Speedup**: For $n = 10000$, interpolation provides $\approx 100$-$1000\times$ speedup in CGF evaluation.

## 6. Residual Centering and Divergence from GRAB 0.2.4 (R)

The SPACox framework treats residuals $R_i$ as fixed and genotypes $G_i$ as random under the empirical permutation null. Under this null, the score statistic $S = \sum_i G_i R_i$ has

$$
\mathbb{E}[S] = \left(\sum_i G_i\right) \bar{R} = 2N \cdot \text{altFreq} \cdot \bar{R}, \qquad
\mathbb{V}(S) = \text{varResid} \cdot \sum_i (G_i - \bar{G})^2 \tag{26}
$$

Two requirements follow:

1. **Unbiasedness of the score test.** The standardized statistic must use the centred score $S - \mathbb{E}[S]$ in the numerator, otherwise $|Z|$ inherits a bias proportional to $\bar{R}$ for every variant. The fastest way to enforce this is to replace $R$ by $R - \bar{R}$ once at the entry of the method; then $\bar{R} = 0$ and $\mathbb{E}[S] = 0$ automatically, so the convenient form $Z = S / \sqrt{\mathbb{V}(S)}$ used by the SPA path is exact.

2. **Invariance under allele recoding.** With $\sum_i R_i = 0$, flipping the genotype to $G' = 2 - G$ yields $S' = 2\sum_i R_i - S = -S$ and $\sum_i (G'_i - \bar{G}')^2 = \sum_i (G_i - \bar{G})^2$, so $|Z|$ is independent of which allele is labelled "ALT". Without centring, $S' \neq -S$ in general, and the implementation must decide on a fixed allele convention before computing $Z$.

### 6.1 Why martingale residuals are the canonical input

If $R$ comes from `survival::coxph()` evaluated at the null (i.e. without the candidate variant), then by construction
$$
\sum_i R_i = 0, \qquad \tilde{\mathbf{X}}^\top \mathbf{R} = \mathbf{0}, \tag{27}
$$
so both requirements above are automatically satisfied and no centring is needed. This is the use case for which SPACox was originally derived (Bi et al. 2020) and is the only configuration in which the GRAB 0.2.4 R implementation is bit-for-bit correct.

### 6.2 The `Residual` traitType path and its pitfalls

GRAB also exposes `traitType = "Residual"`, in which the user supplies an arbitrary numeric vector $R$ and bypasses any null-model fit. In that path there is no guarantee that $\bar{R} = 0$. The GRAB 0.2.4 R implementation handles this path with two coupled choices that both violate the assumptions of (26):

- **No centring of $R$.** `fitNullModel.SPACox()` sets `mresid <- response` for the `Residual` branch without subtracting $\bar{R}$.
- **Genotype flip with unchanged `altFreq`.** `Main.cpp::imputeGenoAndFlip` replaces $G$ by $2 - G$ whenever `altFreq > 0.5`, but `SPACoxClass::getMarkerPval` then computes
  $$
  \text{adjGVec} = G_{\text{flipped}} - 2 \cdot \text{altFreq}_{\text{original}} \tag{28}
  $$
  whose sample mean is $2(1 - \text{altFreq}) - 2 \cdot \text{altFreq} = 2 - 4 \cdot \text{altFreq}$, not zero. The denominator $\sqrt{\sum (\text{adjGVec})^2 \cdot \text{varResid}}$ is therefore inflated by approximately $N \cdot (2 - 4\,\text{altFreq})^2 \cdot \text{varResid}$, and $|Z|$ is biased downward for variants with `altFreq` close to one.

The two issues compound: variants near `altFreq` $= 1$ in the `Residual` path can show $|Z|$ shifts of several units between the flipped and un-flipped implementations of an otherwise identical test. Neither value coincides with the unbiased $Z$ of (26).

### 6.3 This implementation's choice

`src/spacox/spacox.cpp::runSPACox` centres the residual vector once, per phenotype, immediately after extraction:

```cpp
pResid[rc].array() -= pResid[rc].mean();
```

All downstream consumers — the empirical CGF table built by `buildCumulantTable()`, the score `S = GVec.dot(m_resid)` inside `SPACoxMethod::getMarkerPval`, and the SPA root-finder — see a residual vector with $\bar{R} = 0$. As a consequence:

- The score test numerator is unbiased without requiring an explicit subtraction at each marker.
- $|Z|$ is invariant under allele flipping, so no `if (altFreq > 0.5) GVec ← 2 − GVec` block is needed in `getResultVec`.
- For the canonical `time-to-event` path, centring is a numerical no-op (martingale residuals already satisfy (27)), so this implementation produces output that is bit-for-bit identical to GRAB 0.2.4 at six-significant-digit print precision.
- For the `Residual` path with un-centred input the two implementations diverge by exactly the amount of bias $\bar{R}$ introduces. The output of this implementation is the statistically defensible one; reproducing the GRAB 0.2.4 number on the same input would require running the R package directly, or removing the centring line and re-introducing the GVec flip described in (28).

### 6.4 Practical guidance for `Residual` users

If a user supplies residuals from a non-Cox null model (linear regression residuals, deviance residuals from a GLM, etc.), those residuals already satisfy $\sum_i R_i = 0$ by the corresponding orthogonality conditions and no further action is required. If the user supplies a raw phenotype vector (as a numerical convenience, not a model residual), the centring line in `runSPACox` reduces SPACox to a test based on the deviation $R - \bar{R}$; the resulting p-value remains a valid two-sided test under the empirical permutation null, but does not adjust for the original covariates. Users who require covariate adjustment should fit a null model first and pass its residuals to the `Residual` path, or use the `time-to-event` path directly.
