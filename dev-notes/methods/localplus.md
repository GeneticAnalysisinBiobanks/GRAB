# SPAmixLocalPlus

SPAmixLocalPlus performs local-ancestry-specific score tests for admixed
genotype data. For each genetic marker and each ancestry, the method tests the
effect of the allele carried on haplotypes assigned to that ancestry, while
using ancestry-specific sparse Phi matrices to account for relatedness in the
local-ancestry genotype component.

Let:

- $n$ denote the number of used subjects;
- $K$ denote the number of local ancestries;
- $m$ index markers and $k \in \{0,\ldots,K-1\}$ index ancestries;
- $n \times 1$ vector $\mathbf{R} = (R_1,\ldots,R_n)^\top$ represent residuals from the null model;
- $G_{im}^{(k)}$ denote the local-ancestry alternate-allele dosage for subject $i$, marker $m$, and ancestry $k$;
- $H_{im}^{(k)}$ denote the number of haplotypes of subject $i$ at marker $m$ assigned to ancestry $k$;
- $q_m^{(k)}$ denote the local-ancestry allele frequency for marker $m$ and ancestry $k$;
- $\phi_{ij}^{(k,a,b)}$ denote the Phi coefficient for ancestry $k$ and haplotype-count scenario $(H_i,H_j)=(a,b)$, where $(a,b) \in \{(2,2),(2,1),(1,2),(1,1)\}$;
- $\beta_k$ denote the genetic effect for ancestry $k$ at the marker being tested.

The null hypothesis for ancestry $k$ is

$$
H_0: \beta_k = 0. \tag{1}
$$

The following notation suppresses the marker index $m$ and ancestry index $k$
when no ambiguity is possible.

## 1. Local-Ancestry Genotype Representation

For subject $i$, let $A_{i\ell}^{(k)} \in \{0,1\}$ indicate whether haplotype
$\ell \in \{1,2\}$ is assigned to ancestry $k$, and let $B_{i\ell} \in \{0,1\}$
be the alternate-allele indicator on that haplotype. SPAmixLocalPlus stores

$$
H_i^{(k)} = \sum_{\ell=1}^{2} A_{i\ell}^{(k)}, \qquad
G_i^{(k)} = \sum_{\ell=1}^{2} A_{i\ell}^{(k)} B_{i\ell}. \tag{2}
$$

Thus $H_i^{(k)} \in \{0,1,2\}$ and $0 \leq G_i^{(k)} \leq H_i^{(k)}$. Across all
ancestries,

$$
\sum_{k=0}^{K-1} H_i^{(k)} = 2, \qquad
\sum_{k=0}^{K-1} G_i^{(k)} \leq 2. \tag{3}
$$

For a fixed marker and ancestry, the local-ancestry allele frequency is estimated
from the non-missing local haplotypes as

$$
\hat q
= \frac{\sum_{i \in \mathcal{I}_{\mathrm{obs}}} G_i}
       {\sum_{i \in \mathcal{I}_{\mathrm{obs}}} H_i}, \tag{4}
$$

where $\mathcal{I}_{\mathrm{obs}}$ is the set of subjects with finite local
dosage and hapcount. The local minor allele count used for marker filtering is

$$
\mathrm{MAC}
= \min\left\{
\sum_{i \in \mathcal{I}_{\mathrm{obs}}} G_i,\,
\sum_{i \in \mathcal{I}_{\mathrm{obs}}} H_i -
\sum_{i \in \mathcal{I}_{\mathrm{obs}}} G_i
\right\}. \tag{5}
$$

Conditional on the local haplotype count $H_i$, the null genotype model is

$$
G_i \mid H_i, q \sim \mathrm{Binomial}(H_i,q), \tag{6}
$$

so that

$$
\mathbb{E}(G_i \mid H_i,q) = H_i q, \qquad
\mathbb{V}(G_i \mid H_i,q) = H_i q(1-q). \tag{7}
$$

Missing local dosage or hapcount values are set to zero after being counted in
the missing-rate filter. Consequently, missing records make no contribution to
the score, its diagonal mean and variance, or the Phi scenario matching.

## 2. Score Statistic

For a fixed marker and ancestry, the score statistic is

$$
S = \sum_{i=1}^{n} R_i G_i = \mathbf{R}^\top \mathbf{G}. \tag{8}
$$

Using (7), its null mean is

$$
\mu_S
= \mathbb{E}(S \mid \mathbf{H},q)
= q\sum_{i=1}^{n} R_i H_i
= q\,\mathbf{R}^\top \mathbf{H}. \tag{9}
$$

The centered score used in the z-statistic and effect estimate is

$$
U = S - \mu_S. \tag{10}
$$

The reported local-ancestry effect estimate is the score estimate

$$
\hat\beta_k = \frac{U}{\widehat{\mathbb{V}}_{\Phi}(S)}, \tag{11}
$$

where $\widehat{\mathbb{V}}_{\Phi}(S)$ is the Phi-adjusted variance defined in
Section 4.

## 3. Phi Estimation

Phi estimation is performed before GWAS, separately for each ancestry. The
sparse GRM input determines the set of subject pairs on which Phi is estimated.
Let $\mathcal{E}$ denote the off-diagonal sparse-pair support and let

$$
\mathcal{E}^{\rightarrow}
= \{(i,j),(j,i): \{i,j\} \in \mathcal{E},\, i \neq j\} \tag{12}
$$

be the directed pair set used by the implementation. The numerical GRM entry is
not multiplied into the Phi estimator; it is used to define the sparse support.

For ancestry $k$, define the Phi-estimation marker set

$$
\mathcal{M}_{\Phi}^{(k)}
= \{m: c_{\Phi} < \hat q_m^{(k)} < 1-c_{\Phi}\}, \qquad c_{\Phi}=0.01. \tag{13}
$$

For a directed pair $(i,j) \in \mathcal{E}^{\rightarrow}$ and scenario
$(a,b) \in \{(2,2),(2,1),(1,2),(1,1)\}$, define

$$
\mathcal{M}_{ij}^{(k,a,b)}
= \left\{
m \in \mathcal{M}_{\Phi}^{(k)}:
H_{im}^{(k)} = a,\,
H_{jm}^{(k)} = b,\,
G_{im}^{(k)},G_{jm}^{(k)},H_{im}^{(k)},H_{jm}^{(k)} \text{ observed}
\right\}. \tag{14}
$$

Let $M_{ij}^{(k,a,b)} = |\mathcal{M}_{ij}^{(k,a,b)}|$. If
$M_{ij}^{(k,a,b)}>0$, SPAmixLocalPlus estimates

$$
\hat\phi_{ij}^{(k,a,b)}
= \frac{1}{M_{ij}^{(k,a,b)}}
\sum_{m \in \mathcal{M}_{ij}^{(k,a,b)}}
\frac{
\left(G_{im}^{(k)} - H_{im}^{(k)}\hat q_m^{(k)}\right)
\left(G_{jm}^{(k)} - H_{jm}^{(k)}\hat q_m^{(k)}\right)
}{
H_{im}^{(k)}H_{jm}^{(k)}
\hat q_m^{(k)}\left(1-\hat q_m^{(k)}\right)
}. \tag{15}
$$

Because $H_{im}^{(k)}=a$ and $H_{jm}^{(k)}=b$ within the scenario, (15) can be
written equivalently as

$$
\hat\phi_{ij}^{(k,a,b)}
= \frac{1}{M_{ij}^{(k,a,b)}}
\sum_{m \in \mathcal{M}_{ij}^{(k,a,b)}}
\frac{
\left(G_{im}^{(k)} - a\hat q_m^{(k)}\right)
\left(G_{jm}^{(k)} - b\hat q_m^{(k)}\right)
}{
ab\,\hat q_m^{(k)}\left(1-\hat q_m^{(k)}\right)
}. \tag{16}
$$

The four stored scenario labels are

$$
A=(2,2), \qquad B=(2,1), \qquad C=(1,2), \qquad D=(1,1). \tag{17}
$$

Equation (16) implies the working covariance model

$$
\mathrm{Cov}\left(G_i,G_j \mid H_i=a,H_j=b,q\right)
\approx ab\,q(1-q)\,\phi_{ij}^{(k,a,b)}. \tag{18}
$$

Thus $\phi_{ij}^{(k,a,b)}$ is an ancestry-specific, scenario-specific
correlation coefficient for the local-ancestry allele count, normalized by the
number of informative local haplotypes in the pair.

## 4. Phi-Adjusted Score Variance

The diagonal variance under conditional independence is

$$
\widehat{\mathbb{V}}_{\mathrm{diag}}(S)
= q(1-q)\sum_{i=1}^{n} R_i^2 H_i. \tag{19}
$$

The off-diagonal covariance contribution is computed only over stored directed
Phi entries whose scenario matches the current marker's local haplotype counts.
Define

$$
\mathcal{E}_{ab}^{(k)}
= \left\{
(i,j) \in \mathcal{E}^{\rightarrow}:
H_i=a,\ H_j=b,\ \hat\phi_{ij}^{(k,a,b)} \text{ is stored}
\right\}. \tag{20}
$$

Using (18), the Phi-adjusted variance is

$$
\widehat{\mathbb{V}}_{\Phi}(S)
= q(1-q)\sum_{i=1}^{n} R_i^2 H_i
+ q(1-q)
\sum_{(a,b)}
ab \sum_{(i,j)\in \mathcal{E}_{ab}^{(k)}}
\hat\phi_{ij}^{(k,a,b)} R_i R_j, \tag{21}
$$

where the outer sum is over $(a,b)\in\{(2,2),(2,1),(1,2),(1,1)\}$. Expanded by
stored scenario,

$$
\begin{aligned}
\widehat{\mathbb{V}}_{\Phi}(S)
&= q(1-q)\sum_{i=1}^{n} R_i^2 H_i \\
&\quad + 4q(1-q)\sum_{(i,j)\in\mathcal{E}_{22}^{(k)}}
\hat\phi_{ij}^{(k,2,2)}R_iR_j \\
&\quad + 2q(1-q)\sum_{(i,j)\in\mathcal{E}_{21}^{(k)}}
\hat\phi_{ij}^{(k,2,1)}R_iR_j \\
&\quad + 2q(1-q)\sum_{(i,j)\in\mathcal{E}_{12}^{(k)}}
\hat\phi_{ij}^{(k,1,2)}R_iR_j \\
&\quad + q(1-q)\sum_{(i,j)\in\mathcal{E}_{11}^{(k)}}
\hat\phi_{ij}^{(k,1,1)}R_iR_j .
\end{aligned} \tag{22}
$$

No additional factor of two appears in (21) or (22), because
$\mathcal{E}^{\rightarrow}$ is directed. If a symmetric unordered Phi
representation were used instead, the corresponding off-diagonal contribution
would need to be doubled.

The z-statistic reported for ancestry $k$ is

$$
Z_k = \frac{S-\mu_S}{\sqrt{\widehat{\mathbb{V}}_{\Phi}(S)}}. \tag{23}
$$

## 5. Genotype CGF for SPA

For $G_i \mid H_i,q \sim \mathrm{Binomial}(H_i,q)$, the moment generating
function is

$$
M_G(t;q,H_i)
= \mathbb{E}(e^{tG_i}\mid H_i,q)
= (1-q+qe^t)^{H_i}. \tag{24}
$$

The cumulant generating function (CGF) is

$$
K_G(t;q,H_i)
= \log M_G(t;q,H_i)
= H_i\log(1-q+qe^t). \tag{25}
$$

Its first and second derivatives are

$$
K_G'(t;q,H_i)
= \frac{H_i q e^t}{1-q+qe^t}, \tag{26}
$$

and

$$
K_G''(t;q,H_i)
= \frac{H_i q(1-q)e^t}{(1-q+qe^t)^2}. \tag{27}
$$

SPAmixLocalPlus partitions residuals into outliers and non-outliers using the
IQR rule. Let $\mathcal{O}$ denote outlier indices and $\mathcal{N}$ denote
non-outlier indices. For outliers, the diagonal-scale CGF is evaluated exactly:

$$
K_{\mathcal{O}}(t)
= \sum_{i\in\mathcal{O}} K_G(tR_i;q,H_i), \tag{28}
$$

with derivatives

$$
K_{\mathcal{O}}'(t)
= \sum_{i\in\mathcal{O}} R_i K_G'(tR_i;q,H_i), \tag{29}
$$

$$
K_{\mathcal{O}}''(t)
= \sum_{i\in\mathcal{O}} R_i^2 K_G''(tR_i;q,H_i). \tag{30}
$$

For non-outliers, the local score contribution is approximated as normal with

$$
\mu_{\mathcal{N}}
= q\sum_{i\in\mathcal{N}} R_iH_i, \tag{31}
$$

and

$$
\sigma_{\mathcal{N}}^2
= q(1-q)\sum_{i\in\mathcal{N}} R_i^2H_i. \tag{32}
$$

The non-outlier CGF is therefore

$$
K_{\mathcal{N}}(t)
= \mu_{\mathcal{N}}t + \frac{1}{2}\sigma_{\mathcal{N}}^2t^2. \tag{33}
$$

The total diagonal-scale CGF used by the SPA is

$$
K(t)=K_{\mathcal{O}}(t)+K_{\mathcal{N}}(t), \tag{34}
$$

with

$$
K'(t)=K_{\mathcal{O}}'(t)+\mu_{\mathcal{N}}+\sigma_{\mathcal{N}}^2t, \tag{35}
$$

and

$$
K''(t)=K_{\mathcal{O}}''(t)+\sigma_{\mathcal{N}}^2. \tag{36}
$$

## 6. Variance-Ratio Correction and P-Value

The Phi-adjusted variance is used for the normal z-statistic. The SPA tail
calculation uses the diagonal CGF in (34) together with the variance-ratio
correction

$$
\rho
= \frac{\widehat{\mathbb{V}}_{\mathrm{diag}}(S)}
        {\widehat{\mathbb{V}}_{\Phi}(S)}. \tag{37}
$$

The corrected score and corrected mean are

$$
S^\star = S\sqrt{\rho}, \qquad
\mu_S^\star = \mu_S\sqrt{\rho}. \tag{38}
$$

If

$$
|Z_k| < z_{\mathrm{SPA}}, \tag{39}
$$

where $z_{\mathrm{SPA}}$ is the user-specified SPA cutoff, the method uses the
normal approximation

$$
p_{\mathrm{norm}} = 2\Phi_0(-|Z_k|), \tag{40}
$$

where $\Phi_0$ is the standard normal CDF.

Otherwise, define the two-sided tail targets

$$
\delta^\star = |S^\star-\mu_S^\star|, \qquad
s_U = \mu_S^\star+\delta^\star, \qquad
s_L = \mu_S^\star-\delta^\star. \tag{41}
$$

For each target $s \in \{s_U,s_L\}$, solve the saddlepoint equation

$$
K'(\zeta)=s \tag{42}
$$

by Newton-Raphson iteration. Given the solution $\zeta$, compute

$$
w = \mathrm{sgn}(\zeta)
\sqrt{2\{\zeta s-K(\zeta)\}}, \qquad
v = \zeta\sqrt{K''(\zeta)}. \tag{43}
$$

The Lugannani-Rice approximation gives the tail probability for the corrected
statistic:

$$
\Pr(S^\star \leq s)
\simeq
\Phi_0\left(
w + \frac{1}{w}\log\frac{v}{w}
\right). \tag{44}
$$

Thus the upper and lower tail probabilities are approximated by

$$
p_U
\simeq
1-\Phi_0\left(
w_U + \frac{1}{w_U}\log\frac{v_U}{w_U}
\right), \tag{45}
$$

and

$$
p_L
\simeq
\Phi_0\left(
w_L + \frac{1}{w_L}\log\frac{v_L}{w_L}
\right), \tag{46}
$$

where $(w_U,v_U)$ and $(w_L,v_L)$ are evaluated at $s_U$ and $s_L$,
respectively. The two-sided SPA p-value is

$$
p_{\mathrm{SPA}} = p_U+p_L. \tag{47}
$$

If no residual outliers are detected, or if $|Z_k|$ is below the SPA cutoff,
$p_{\mathrm{SPA}}$ is set equal to the normal p-value in (40).

## 7. Multi-Phenotype Form

The implementation supports multiple residual columns. Let
$\mathbf{R}^{(p)}=(R_{1p},\ldots,R_{np})^\top$ denote phenotype-specific
residuals. Equations (8)--(47) apply independently to each phenotype after
replacing $R_i$ by $R_{ip}$.

For computational efficiency, the Phi entries are pre-multiplied as

$$
c_{ab}\,\hat\phi_{ij}^{(k,a,b)} R_{ip}R_{jp},
\qquad c_{ab}=ab, \tag{48}
$$

and stored contiguously across phenotypes. For each marker, a stored entry
contributes to the off-diagonal variance only when the current hapcount pair
matches its scenario $(a,b)$.

## 8. Output Quantities

For each marker and each ancestry, SPAmixLocalPlus reports:

- missing rate, $\#\mathrm{missing}/n$;
- local allele frequency $\hat q$ from (4);
- local MAC from (5);
- p-value $p_{\mathrm{SPA}}$ from (47), or the normal p-value from (40) when SPA is not invoked;
- score effect estimate $\hat\beta_k$ from (11);
- standard error $\widehat{\mathrm{SE}}_k = 1 / \sqrt{\widehat{\mathrm{Var}}(S_k)}$, with $\widehat{\mathrm{Var}}(S_k)$ as in (23).

The z-statistic from (23) is not emitted as a separate column; it is recoverable
as $Z_k = \hat\beta_k / \widehat{\mathrm{SE}}_k$.

Markers failing ancestry-specific filters for missingness, allele frequency, or
MAC are reported as missing for that ancestry.
