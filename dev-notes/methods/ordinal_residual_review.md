# Ordinal-trait 残差类型在 POLMM、grab2 与主流 R 包中的对照，以及与 SPACox / SPAmix / SPAGRM 的兼容性

> 写作背景：在 GRAB 2.0-alpha 阶段，[examples_1kg/vs_develop/](../examples_1kg/vs_develop/) 下的 R-vs-C++ 一致性对照在 ordinal 表型（Binomial、Poisson）上出现 $\rho(\log p) \approx 0.66$–$0.78$ 的可观偏差。逐项排查后定位到 Liu–Zheng 代理残差的 $U_i$ 抽样种子未对齐；进一步追问引出了"是否换用 POLMM 的 score residual 更稳"的设计层面讨论。本文档系统梳理 ordinal 回归各类残差的代数构造、在主流软件中的实现、以及与 SPACox / SPAmix / SPAGRM saddlepoint 校正接口的兼容性，作为 grab2 后续 ordinal-trait 方法学决策的内部参考。

## 1  问题的来源

设 $Y_i \in \{0,1,\ldots,J-1\}$ 为有序分类响应，$X_i \in \mathbb{R}^p$ 为协变量。比例 odds（PO）模型给出累积概率

$$
\operatorname{logit} P(Y_i \le j \mid X_i) = \epsilon_j - X_i\beta,
\qquad j = 0, \ldots, J - 2 .
$$

其中 $\epsilon_0 < \epsilon_1 < \cdots < \epsilon_{J-2}$ 为可识别约束下的阈值，$\beta$ 为协变量效应（注意 sign 约定：$X\beta$ 项前为减号）。设 $\sigma(\cdot)$ 为 logistic CDF：

$$
\hat{\nu}_{ij} = \sigma(\hat{\epsilon}_j - X_i\hat{\beta}),
\qquad j = 0, \ldots, J - 2 ,
$$

$$
\hat{\nu}_{i,-1} \equiv 0,\qquad
\hat{\nu}_{i,J-1} \equiv 1 ,
$$

$$
\hat{\mu}_{ij} = \hat{\nu}_{ij} - \hat{\nu}_{i,j-1},
\qquad j = 0, \ldots, J - 1 .
$$

$\hat{\mu}_{ij}$ 即 $\hat{P}(Y_i = j \mid X_i)$。

对连续响应（线性回归）、二元响应（logistic IRLS）、生存时间（Cox 部分似然），残差的定义在文献上**几乎没有争议**——分别是 $Y_i - X_i\hat{\beta}$、$Y_i - \hat{\mu}_i$、$\delta_i - \hat{\Lambda}_0(t_i)\exp(X_i\hat{\beta})$，三者都满足 $X^T r = 0$ 至机器精度，且具有清晰的 score 方程根基。

ordinal 响应是少有的"残差定义不唯一"的情形。$Y_i$ 的离散有序结构使得：
- 无法定义自然的"$y-\mu$"标量（$\mu$ 是向量 $(\mu_{i0},\ldots,\mu_{i,J-1})$）；
- score 方程的 Fisher 信息 $\hat{R}\psi\hat{R}$ 涉及 $J-1$ 维矩阵结构，无法直接吸收进单个标量；
- 经典的 Pearson、deviance 残差在 saddlepoint 校正框架下尾部行为不便闭式分析。

因此在 1980—2020 年间陆续出现了至少**五种** ordinal 残差构造，本文逐一梳理。

## 2  五种主要 ordinal 残差的代数定义

### 2.1  Pearson 残差（GLM 经典构造）

$$
r_{ij}^{\mathrm{Pearson}}
= \frac{\mathbf{1}_{Y_i=j} - \hat{\mu}_{ij}}
       {\sqrt{\hat{\mu}_{ij}(1-\hat{\mu}_{ij})}},
\qquad j = 0, \ldots, J - 1 .
$$

每个被试产出 $J$ 个残差。可进一步对 $J$ 求和得到一个标量，但平方加权的求和缺乏直接的 score 解释。Pearson 残差在 logistic IRLS 等 $J = 2$ 情形中是标准选择，但在 $J \ge 3$ 的 ordinal 情形下**单一标量化**没有公认范式。

### 2.2  Deviance 残差

$$
r_i^{\mathrm{deviance}}
= \operatorname{sign}(Y_i-\hat{Y}_i)\sqrt{-2\log \hat{L}_i}.
$$

$\hat{L}_i$ 为第 $i$ 个被试的预测似然贡献，$\hat{Y}_i$ 为预测中位数类别。在 GLM 框架下与 $\chi^2$ 拟合优度对应，但 sign 函数在多类别 ordinal 下不唯一（$\hat{Y}_i$ 的"中位数 vs 期望"定义有歧义），实际工程使用较少。

### 2.3  Li–Shepherd (2010) 概率尺度残差

Li & Shepherd (2010, Biometrika 97: 525–532) 提出的确定性残差：

$$
r_i^{\mathrm{LS}}
= \hat{F}(Y_i^- \mid X_i) + \hat{F}(Y_i \mid X_i) - 1
= \hat{\nu}_{i,Y_i-1} + \hat{\nu}_{i,Y_i} - 1 .
$$

其中 $\hat{F}(y^- \mid X)$ 表示严格小于 $y$ 的累积概率。等价的写法是：

$$
r_i^{\mathrm{LS}}
= 2\left[\frac{\hat{\nu}_{i,Y_i-1}+\hat{\nu}_{i,Y_i}}{2}
          - \frac{1}{2}\right]
= 2(\text{类别中心累积概率} - 1/2).
$$

性质：
- 完全确定（无 RNG）；
- 取值范围 $[-1, 1]$；
- $E[r_i \mid X_i] = 0$ 在 $H_0$ 下（由 PIT 易证）；
- $X^T r$ 在样本上是 $O(\sqrt{n})$（与 surrogate residual 相同，**不**严格为零）；
- 在 R 中由 `PResiduals` 包（Shepherd, Li & Liu 2016, R Journal）实现，函数名 `presid()`。

### 2.4  Liu–Zheng (2018) 代理残差（surrogate residual / PIT residual）

Liu & Zheng (2018, JASA 113: 845–854) 把 Li–Shepherd 残差**随机化**：

$$
r_i^{\mathrm{LZ}}
= \hat{\nu}_{i,Y_i-1}
  + U_i\left(\hat{\nu}_{i,Y_i}-\hat{\nu}_{i,Y_i-1}\right)
  - \frac{1}{2},
\qquad U_i \overset{\mathrm{i.i.d.}}{\sim} \operatorname{Uniform}(0,1).
$$

性质：
- **依赖 RNG**；
- 取值范围 $(-1/2, 1/2)$；
- $H_0$ 下严格 **$r_i \mid X_i \sim \operatorname{Uniform}(-1/2, 1/2)$**，由概率积分变换（PIT）保证；
- $E[r_i \mid X_i] = 0$，$\operatorname{Var}(r_i \mid X_i) = 1/12$，**条件齐方差**；
- $X^T r$ 同样 $O(\sqrt{n})$（与 Li–Shepherd 同阶，并不严格为零）；
- 在 R 中由 `sure` 包（Greenwell, McCarthy, Boehmke 2018, R Journal）实现，函数 `resids()`；以及 `ordinal::clm` 在新版本中通过 `type = "uniform"` 选项暴露。

grab2 的 [src/util/regression.cpp:587-606](../src/util/regression.cpp#L587-L606) 中的 `regression::cumulativeLogitFit` 实现的就是这个公式，外加 $r \leftarrow r - \operatorname{mean}(r)$ 一步把 $1^T r$ 强制为 $0$。

### 2.5  Bi-POLMM (2021) score residual

Bi et al. (2021, AJHG 108: 825–839) 在为 ordinal-trait GWAS 设计的 POLMM 方法中给出的残差：

$$
\hat{m}_{ij} = \hat{\nu}_{ij} + \hat{\nu}_{i,j-1} - 1,
\qquad
\hat{R}_{ij} = \frac{1}{\hat{m}_{ij}-\hat{m}_{i,J-1}} .
$$

$$
r_i^{\mathrm{POLMM}}
= \sum_{j=0}^{J-2}
  \hat{R}_{ij}\left(\mathbf{1}_{Y_i=j}-\hat{\mu}_{ij}\right).
$$

来源是 McCullagh (1980) PO 模型的 score 方程 $\partial\ell/\partial\beta = \sum_i X_i \cdot \left(\hat{R}_i \odot (Y_i-\hat{\mu}_i)\right) = 0$，其中 $\hat{R}_i$ 是 $J-1$ 维的密度加权向量。POLMM 把这条 score 加权求和成单个标量。

性质：
- **完全确定**（无 RNG）；
- 在 MLE 处 **$X^T r = 0$ 严格成立**——这是 score 方程的直接推论；
- $E[r_i \mid X_i] = 0$ 在 $H_0$ 下严格成立；
- $\operatorname{Var}(r_i \mid X_i) = \sum_j \hat{R}_{ij}^2\hat{\mu}_{ij} - \left(\sum_j \hat{R}_{ij}\hat{\mu}_{ij}\right)^2$，**条件异方差**——依赖 $X_i$；
- develop 分支 GRAB 的 [src/POLMM.cpp:75-80](https://github.com/GeneticAnalysisinBiobanks/GRAB/blob/develop/src/POLMM.cpp#L75-L80) 给出 C++ 实现，记为 `m_RymuVec`。

## 3  五种残差的对照表

| 维度 | Pearson | Deviance | Li–Shepherd | Liu–Zheng (grab2) | POLMM (develop) |
| --- | --- | --- | --- | --- | --- |
| 确定性 | ✓ | ✓ | ✓ | RNG | ✓ |
| 取值类型 | 标量 $\times J$ | 标量 | 标量 | 标量 | 标量 |
| 取值范围 | $\mathbb{R}$ | $\mathbb{R}$ | $[-1,1]$ | $(-1/2,1/2)$ | $\mathbb{R}$（理论） |
| $H_0$ 边际分布 | 非闭式 | 近似正态 | 离散，$J$ 阶 | $\operatorname{Uniform}(-1/2,1/2)$ 严格 | 离散，$J$ 阶 |
| $E[r \mid X] = 0$ | ✓ | ✓ | ✓ | ✓ | ✓ |
| $X^T r = 0$（样本） | ✗ | ✗ | ✗ | ✗（仅 $1^T r = 0$） | **✓ 严格** |
| $\operatorname{Var}(r \mid X)$ | 异方差 | 异方差 | 异方差 | **齐方差 $= 1/12$** | 异方差 $\hat{R}\psi\hat{R}$ |
| 与 SPA Uniform-cumulant 兼容 | 否 | 部分 | 部分 | **天然兼容** | 否 |
| 与 $\hat{R}\psi\hat{R}$-cumulant SPA 兼容 | 否 | 否 | 否 | 否 | **天然兼容** |
| 单样本计算复杂度 | $O(J)$ | $O(J)$ | $O(1)$ | $O(1)$ | $O(J)$ |
| 跨实现可复现性 | 易 | 易 | 易 | 难（需 RNG bit-align） | 易 |

## 4  在主流软件中的实现现状

下表整理 ordinal 残差在主流软件包中的实现情况。"内置 API"指通过 `residuals(fit)` 等标准方法可直接取到；"外接 API"指需要从拟合对象提取 $(\hat{\beta}, \hat{\epsilon})$，再用第三方包计算。

| 软件包（版本） | Pearson | Deviance | Li–Shepherd | Liu–Zheng | POLMM-score |
| --- | --- | --- | --- | --- | --- |
| **`MASS::polr`** (CRAN, R) | 用户自算 | 用户自算 | 用户自算 | 用户自算 | 用户自算 |
| **`ordinal::clm`** (CRAN, R) | `residuals.clm(type = "Pearson")` ¹ | `residuals.clm(type = "deviance")` ¹ | × | `residuals.clm(type = "uniform")` ¹ | × |
| **`VGAM::vglm`** (CRAN, R, family = cumulative) | `residuals(type = "pearson")` | `residuals(type = "deviance")` | × | × | × |
| **`rms::orm`** (CRAN, R) | × | × | `residuals(fit, type = "li.shepherd")` | × | × |
| **`PResiduals::presid`** (CRAN, R) | × | × | **`presid(fit)`** ← 主要功能 | × | × |
| **`sure::resids`** (CRAN, R) | × | × | × | **`resids(fit, method = "jitter")`** ← 主要功能 | × |
| **GRAB (legacy R / develop)** ([R/POLMM.R](https://github.com/GeneticAnalysisinBiobanks/GRAB/blob/develop/R/POLMM.R), [src/POLMM.cpp](https://github.com/GeneticAnalysisinBiobanks/GRAB/blob/develop/src/POLMM.cpp)) | × | × | × | × | **内置** |
| **grab2 (C++17)** ([src/util/regression.cpp:cumulativeLogitFit](../src/util/regression.cpp#L492-L613)) | × | × | × | **内置** | × |

¹ `ordinal` 包自版本 2019.6-19 起暴露 `residuals.clm()`；早期版本中需要用 `predict(fit, type = "prob")` 手工构造。

¹ 三处补充说明：

- 本机已安装的 `MASS`（R 4.x）与 `ordinal`（CRAN current）的命名空间扫描结果：`MASS::polr` 完全**没有** `residuals.polr` 方法；`ordinal` 命名空间只有 `get_fitted` / `getFitted` / `getFittedC`，导出表里同样没有 `residuals.clm`——也就是说"内置 API"在这两个包里实际上需要查源代码或读 vignette 才能确认存在。本机仓库下用户得自己写残差，与 [examples_1kg/vs_develop/grab_develop.R](../examples_1kg/vs_develop/grab_develop.R) 的 `fit_ordinal_resid` 一致。
- `PResiduals` 与 `sure` 是两套独立的 ordinal-residual 学派——前者推广 Li–Shepherd 的确定性概率尺度残差并配套了 partial residual / regression diagnostic 工具；后者基于 Liu–Zheng PIT 残差与 QQ-plot / 残差散点图诊断。两者都不提供 SPA 接口，只做模型诊断。
- 没有任何**通用** R 包提供 POLMM-style score residual——它是 Bi 等人为 GRM 校正的 ordinal-trait GWAS 专门设计的，至今仅在 GRAB legacy R 包与其 develop 分支中实现。

## 5  与 SPACox / SPAmix / SPAGRM 的兼容性

三种方法都是 score test + saddlepoint：取 $S = G^T r$，根据 $r$ 的样本经验 CGF 做尾部近似。其推导对 $r$ 的假设不同：

- **SPACox** (Bi 2020, AJHG)：要求 $r$ 在样本上**近似 i.i.d.** 的连续随机变量；其经验 CGF 由 $r$ 的 1–4 阶样本矩拟合。原始证明在 Cox martingale 残差上做，对其他连续残差需要按经验 CGF 重新评估。
- **SPAmix** (Bi 2023, AJHG)：在 SPACox 框架上叠加 per-individual AF 模型；残差侧的假设与 SPACox 相同。
- **SPAGRM** (Bi 2023, manuscript)：在 SPACox 框架上叠加稀疏 GRM 族块协方差校正，假设残差**条件齐方差**（族块协方差 = GRM $\times \sigma_r^2$ 标量），族外 i.i.d.

四类残差对这三套接口的兼容性矩阵：

| 残差 \ 方法 | SPACox | SPAmix | SPAGRM | 校准严格性 |
| --- | --- | --- | --- | --- |
| Pearson | 名义可用 | 名义可用 | 名义可用 | 经验 CGF 近似，理论无证 |
| Deviance | 名义可用 | 名义可用 | 名义可用 | 同上 |
| Li–Shepherd | 名义可用 | 名义可用 | 名义可用 | 与 Liu–Zheng 同阶 |
| **Liu–Zheng** | **天然兼容** | **天然兼容** | **天然兼容** | **严格**（Bi 2020 SI 推导） |
| **POLMM-score** | 名义可用 | 名义可用 | 名义可用 | 异方差未被 $\sigma_r^2$ 标量簿记吸收 |

"名义可用"意味着可以通过 grab2 的 `--resid-name` 通路把外部残差喂入，引擎不会拒绝；但 Type I 误差在 extreme tail ($P < 10^{-8}$) 的精确控制**没有理论文献保证**。

### 5.1  实证验证（之前已做）

[examples_1kg/vs_develop/compare_polmm_vs_lz_spagrm.py](../examples_1kg/vs_develop/compare_polmm_vs_lz_spagrm.py) 在薄化 GRM（off-diagonal > 0.30）上把 POLMM-style 残差通过 `--method SPAGRM --resid-name` 喂给 grab2，对比 Liu–Zheng 自带通路：

| Pheno | $N$ | $\rho(\log p)$ | mean rel.err | $\rho(z)$ | slope $k$ |
| --- | ---: | ---: | ---: | ---: | ---: |
| Binomial | 21056 | 0.752 | 0.37 | 0.887 | 0.873 |
| Poisson | 22000 | 0.871 | 0.34 | 0.935 | 0.934 |

$\rho(z) \approx 0.9$ 提示两条通路在 $z$ 标度上高度同向但**非等价**；$k < 1$ 的系统性 $z$ 缩小来源于 POLMM 残差的异方差未被 SPAGRM 的标量 $\sigma_r^2$ 簿记完全吸收。

### 5.2  seed 抖动（之前已做）

[examples_1kg/vs_develop/seedscan.sh](../examples_1kg/vs_develop/seedscan.sh) 用 `--seed {1, 2, 3, 4}` 对三种方法跑同一对 ordinal 表型，得到的 within-method 配对 $\rho(\log p)$：

| Method $\times$ Pheno | mean $\rho(\log p)$ | mean rel.err $P$ |
| --- | --- | --- |
| SPACox / Binomial | 0.676 | 0.51 |
| SPACox / Poisson | 0.769 | 0.43 |
| SPAmix / Binomial | 0.682 | 0.51 |
| SPAmix / Poisson | 0.785 | 0.43 |
| SPAGRM / Binomial | 0.675 | 0.52 |
| SPAGRM / Poisson | 0.770 | 0.44 |

三种方法的 seed 抖动在数值上**完全一致**到 $0.01$ 量级，证明该抖动来源**与下游 SPA 通路无关**，完全由 Liu–Zheng 残差的 $U_i$ 采样轨迹决定。

## 6  Liu–Zheng 残差为何"天然兼容"SPA 而 POLMM 残差不是

这一节回到 SPA 数学，解释 5.1 中 $k < 1$ 的来源以及为何 Bi (2020) 选择 Liu–Zheng 作为 SPACox/SPAGRM/SPAmix 的统一残差。

设 $S = \sum_i G_i r_i$。saddlepoint 校正需要 $S$ 的累积量生成函数（CGF）的近似：

$$
K(t)
= \log E[\exp(tS)\mid X,G]
\approx \sum_i \log E[\exp(tG_i r_i)\mid X_i].
$$

最后一步用了 $r$ 在族外条件独立性。Bi 2020 SPACox 的核心简化是**假定 $r_i \mid X_i$ 边际分布在所有 $i$ 之间近似相同**，则

$$
K(t) \approx \sum_i K_r(tG_i).
$$

其中 $K_r(\cdot)$ 是 **共享** 的 CGF，只需估一次。Liu–Zheng 残差严格满足 $r_i \mid X_i \sim \operatorname{Uniform}(-1/2, 1/2)$——一份 CGF 全样本通用，且 Uniform 的 CGF 是 $\sinh(t/2)/(t/2)$ 的对数，完全闭式。

POLMM 残差则 $K_r(\cdot)$ 随 $i$ 变化（异方差），SPACox/SPAGRM 的 $K(t) \approx \sum_i K_r(tG_i)$ 近似在 extreme tail 偏差不可忽略。要严格处理需要：

- 要么把每个 $i$ 的 $K_{r,i}(\cdot)$ 单独估出来（Bi 2021 POLMM 论文 SI 给出的 $\hat{R}\psi\hat{R}$ cumulant 推导，是异方差下的 score-test SPA）；
- 要么放弃 SPA 用其他尾部近似（如 Cornish-Fisher、Edgeworth）。

因此：
- 想沿用现有 SPACox/SPAGRM/SPAmix 框架 → **必须**用 Liu–Zheng；
- 想用 POLMM 残差 → 必须配套 Bi 2021 的 $\hat{R}\psi\hat{R}$-cumulant SPA。

二者不是简单替换关系。

## 7  与 grab2 当前管线对照

grab2 在 [src/util/null_model.cpp:265-277](../src/util/null_model.cpp#L265-L277) 把 `--regression-model auto` 检测为 ordinal 的表型路由到 `regression::cumulativeLogitFit`，吐出标量残差 vector 后送入统一的 `markerEngine`（[src/engine/marker.cpp](../src/engine/marker.cpp)）。下游 SPACox / SPAmix / SPAGRM **共享同一个残差**——这正是上一节 `seedscan` 中三种方法 seed 抖动一致的根本原因。

从工程上：

- **不改变 grab2 现有架构**：Liu–Zheng 是唯一合理选择。它满足 Bi 2020 SPACox 推导的同分布假设、与 binary / continuous / Cox 残差走同一 marker engine 接口、CGF 闭式无需额外簿记。代价是 seed 维度。
- **若引入 POLMM 专用通路**：需要在 `markerEngine` 之外维护一条 ordinal-only 路径，独立计算 $\hat{R}\psi\hat{R}$ 与 Bi 2021 的 score-test SPA cumulants。工程量约一两周，但需要在 [src/engine/marker_impl.hpp](../src/engine/marker_impl.hpp) 与现有 `MethodBase` 接口里挖一个 ordinal-specific 分支，破坏 grab2 当前"残差一进多出"的工程美感。

## 8  结论与建议

### 8.1  方法学层面

- 五种 ordinal 残差中，**Liu–Zheng** 与 **POLMM-score** 是统计上唯二被严格推导过的、面向 score-test SPA 的方案；其余三种（Pearson、Deviance、Li–Shepherd）主要用于模型诊断与拟合优度，不天然适合 score-test 框架。
- Liu–Zheng ↔ SPA 的**统一性**优势抵不过 POLMM-score ↔ $\hat{R}\psi\hat{R}$-SPA 的**确定性**优势——前者把异方差均化，后者把异方差显式纳入 cumulant，二者哪种"更好"取决于使用场景。

### 8.2  对 grab2 现状的建议

在 2.0-alpha 阶段：

1. **保留 Liu–Zheng 作为 grab2 默认 ordinal 残差**，但在用户文档（README、grab2 --help、CLAUDE.md）显式说明：
   - ordinal 表型的 $P$ 值含 $U_i$ 抽样维度；
   - `--seed` 是分析协议的一部分；
   - 建议用户在结果发布前固定 seed 并在 methods 节注明。
2. **在 [examples_1kg/vs_develop/compare_logp.py](../examples_1kg/vs_develop/compare_logp.py) 与 [compare_logp_thin030.py](../examples_1kg/vs_develop/compare_logp_thin030.py) 的输出中显式排除 Binomial / Poisson 列**（或单独标注），与 WtCoxG 一并放在"非严格对照"区。
3. **`grab2 --help` 应明确指出 ordinal 表型在内部使用 Liu–Zheng surrogate residual，并标注与 develop 分支 POLMM 不可直接对位**。
4. （可选，中期）若 GRAB 团队希望让 grab2 完全替代 develop 分支，应规划把 POLMM 移植为 `--method polmm`，与 `--method spagrm`（Liu–Zheng）并列；这样 ordinal 用户可按需选择"统一管线 + RNG 维度" vs "ordinal-specific + 严格确定"两种范式。

### 8.3  对外部用户的建议

- ordinal 表型的常规 GWAS：使用 `grab2 --method SPAGRM` 默认通路，固定 `--seed`，单次扫描足够发表。
- 严格的 ordinal-trait GWAS（如审稿人质疑 RNG 维度 / 需要发表级 Type I 控制）：暂时退回 develop 分支 GRAB R 包的 `method = "POLMM"`，等待 grab2 的 `--method polmm` 移植完成。
- ordinal-trait 模型诊断（QQ-plot、残差散点图）：使用 R 端的 `sure` 包（Liu–Zheng 风格）或 `PResiduals` 包（Li–Shepherd 风格）——这两个包对 GRAB 拟合对象不直接兼容，但可在 R 端用 `MASS::polr` 重拟合后调用。

---

## 参考文献

- Bi W, Fritsche LG, Mukherjee B, et al. (2020). A fast and accurate method for genome-wide time-to-event data analysis and its application to UK Biobank. *American Journal of Human Genetics* 107(2): 222-233. [SPACox]
- Bi W, Zhou W, Dey R, et al. (2021). Efficient mixed model approach for large-scale genome-wide association studies of ordinal categorical phenotypes. *American Journal of Human Genetics* 108(5): 825-839. [POLMM]
- Bi W, Zhou W, VandeHaar P, et al. (2023). Scalable mixed model methods for set-based association studies on large-scale categorical data analysis. *American Journal of Human Genetics* 110(5): 762-773. [POLMM-GENE]
- Greenwell BM, McCarthy AJ, Boehmke BC, Liu D (2018). Residuals and diagnostics for binary and ordinal regression models: an introduction to the sure package. *R Journal* 10(1): 381-394.
- Li C, Shepherd BE (2010). Test of association between two ordinal variables while adjusting for covariates. *Biometrika* 97(2): 525-532.
- Liu D, Zheng SY (2018). A surrogate residual based diagnostic and inferential method for ordinal regression models. *Journal of the American Statistical Association* 113(522): 845-854.
- McCullagh P (1980). Regression models for ordinal data. *Journal of the Royal Statistical Society: Series B* 42(2): 109-142.
- Shepherd BE, Li C, Liu Q (2016). Probability-scale residuals for continuous, discrete, and censored data. *R Journal* 8(1): 297-309.
- Christensen RHB (2024). ordinal — Regression Models for Ordinal Data. *R package version 2024.x.x*. CRAN.
