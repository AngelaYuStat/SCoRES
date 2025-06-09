
<!-- README.md is generated from README.Rmd. Please edit that file -->

# invSCI

<!-- badges: start -->

[![Codecov test coverage](https://codecov.io/)](https://codecov.io/)
<!-- badges: end -->

### What it does

|  |
|----|
| The identification of domain sets whose outcomes belong to predefined subsets can address fundamental risk assessment challenges in climatology and medicine. A motivating example involves estimating geographical regions where average difference between summer and winter temperatures exceed a certain benchmark, which help policymakers focus on specific areas that are at higher risk for effects of climate change. |
| Mathematically, the target region correspond to the inverse image of $U \subset  \mathbb{R}$ under an unknown function $\mu: \mathcal{S} \to \mathbb{R}$, can be defined as $\mu^{-1}(U) = \{s \in S: \mu(s) \in U\}$, with $U$ a pre-specified subset of a real line $\mathbb{R}$ (e.g., $[c, \infty)$). |
| A point estimator for the inverse set can be constructed as $\hat{\mu}_n^{-1}(U)$, where $\hat{\mu}_n$ is an empirical estimator of $\mu$ based on $n$ observations. To quantify the spatial uncertainty of this estimation, Sommerfeld et al. (2018) introduced Coverage Probability Excursion (CoPE) sets, defined as:  which satisfy:  for a pre-specified confidence level $1-\alpha$ (e.g., $\alpha = 0.05$). |
| Existing approaches require restrictive assumptions, including domain density of $S$ in $R$, continuity of $\hat{\mu}_n$ and $\mu$ near thresholds, and large-sample guarantees, which limit the applicability. Besides, the estimation and coverage depend on setting a fixed threshold level, which is difficult to determine. |
| Ren et al. (2023) proposed a framework that generalizes the estimation of such inverse sets to dense and non-dense domains with protection against inflated Type I error, and constructs multiple upper, lower or interval confidence sets of $\mu^{-1}(U)$ over arbitrary chosen thresholds. The coverage probability is achieved non-asymptotically and simultaneously through inverting simultaneous confidence intervals. For instance, suppose we are interested in inverse set $\mu^{-1}([c,\infty))$ for a single value $c$, the inverse confidence sets (CSs) are constructed by inverting simultaneous confidence intervals (SCIs). Given SCI bounds $\hat{B}_{l}(\boldsymbol{s})$ and $\hat{B}_{u}(\boldsymbol{s})$ satisfying:  the inner and outer CSs for the inverse upper excursion set $\mu^{-1}[c,\infty)$ are defined as:  For details on the estimation of inverse lower excursion set and inverse interval set, see the corresponding package vignette. This package provides useful statistical tools for both the estimation of the inverse set and the corresponding simultaneous confidence intervals. Acceptable forms of input includes both 1D and 2D data for linear regression, logistic regression, functional regression and so on. More details can be found below. |
| \### Installation |

To install from `CRAN`, please use:

``` r
# install.packages("invSCI")
```

To install the latest version directly from Github, please use:

``` r
install.packages("devtools")
devtools::install_github("AngelaYuStat/invSCI")
```

### How to use it

------------------------------------------------------------------------

The first example here is to use ccds data to construct the inverse
confidence sets (CS) from simultaneous confidence bands (SCB).

The ccds dataset contains repeated measures of percent change over time
for multiple subjects under two user categories (use: 1 and no use: 0).
It includes both user and non-user groups, time points, and metadata
related to eye side and frame timing. cleaning process make sure that
the data only includes measurements taken from the right eye at the
post-intervention timepoint (`tp == "post"`).

``` r
library(invSCI)
data(ccds)
```

Before calculating the SCBs, we first process ccds data by fitting a
mean GAM model, extracting residuals and performing FPCA using
`invSCI::prepare_ccds_fpca()`, the function will return an enhanced
dataset includes the FPCA-derived basis scores (Phi1, Phi2, Phi3, Phi4)
for Function-on-Scalar Regression (FoSR) analysis.

Following the FPCA-based data augmentation, we fit a FoSR model using
`mgcv::bam()`, which allows efficient estimation of Generalized Additive
Mixed Models (GAMMs). The model formula is designed to capture both
population-level smooth trends and subject-specific functional
variation.

``` r
data(ccds)
ccds_fpca <- prepare_ccds_fpca(ccds)
fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
            s(seconds, by = use, k=30, bs = "cr") +
            s(subject, by = Phi1, bs="re") +
            s(subject, by = Phi2, bs="re") +
            s(subject, by = Phi3, bs="re") +
            s(subject, by = Phi4, bs="re"),
            method = "fREML", data = ccds_fpca, discrete = TRUE)
```

After obtaining the FoSR object `fosr_mod`, simultaneous confidence
bands (SCB) can be constructed though `invSCI::SCB_functional_outcome`
using pre-specified methods. The `invSCI` package provides two options
for calculating the simultaneous confidence bands (SCB). Use `method` to
specify. Use `groups` to specify the names of grouping variables to
analyze. The input data should have numerical binary 0/1 values for all
scalar group variables. Here, we analyze the user group by specifying
`groups = "use"`.

``` r
# CMA approach
results <- SCB_functional_outcome(object = fosr_mod, method = "cma", est_mean = TRUE, outcome = "percent_change", time = "seconds", groups = "use", subject = "subject")

# Wild bootstrap
results <- SCB_functional_outcome(data = ccds, object = fosr_mod, method = "wild", est_mean = TRUE, outcome = "percent_change", time = "seconds", groups = "use", subject = "subject")
```

The followings are the mathematical details:

#### Correlation and Multiplicity Adjusted (CMA) Confidence Bands Based on Parameter Simulations

Our goal is to obtain the quantile via simulations

$$
q(C_f, 1 - \alpha),
$$

which is the $1 - \alpha$ quantile of

$$
\frac{\|\hat{\mathbf{f}} - \mathbf{f}\|}{\mathbf{D}_f},
$$

where $\hat{\mathbf{f}}$ is the estimated functional effect,
$\mathbf{f}$ is the true effect, and $\mathbf{D}_f$ is a standard
deviation vector. The distribution of this quantity is assumed to be
multivariate normal:

$$
\frac{\hat{\mathbf{f}} - \mathbf{f}}{\mathbf{D}_f} \sim \mathcal{N}(0_n, C_f),
$$

where $C_f$ is a known covariance matrix.

To approximate this, we generate $B$ simulated vectors
$\mathbf{X}_1, \dots, \mathbf{X}_B \sim \mathcal{N}(0_n, C_f)$, and
define:

$$
d_b = \max(|\mathbf{X}_b|), \quad b = 1, \dots, B,
$$

where $|\mathbf{X}_b|$ is the element-wise absolute value. Then the
estimate of $q(C_f, 1 - \alpha)$ is the $100 \times (1 - \alpha)$-th
percentile of $\{d_1, \dots, d_B\}$.

However, directly sampling from a high-dimensional normal distribution
can be computationally expensive. To avoid this, we use the fact that we
can simulate from a lower-dimensional distribution of the parameters:

$$
\boldsymbol{\beta}_1, \dots, \boldsymbol{\beta}_B \sim \mathcal{N}(\hat{\boldsymbol{\beta}}, \hat{V}_{\boldsymbol{\beta}}),
$$

where $\boldsymbol{\beta}$ are the model parameters, and $\mathbf{B}$ is
the matrix mapping $\boldsymbol{\beta}$ to $\mathbf{f}$. We define:

$$
\mathbf{X}_b = \frac{\mathbf{B}(\boldsymbol{\beta}_b - \hat{\boldsymbol{\beta}})}{\mathbf{D}_f},
$$

where the division is element-wise. Then the $\mathbf{X}_b$ are
approximately distributed as $\mathcal{N}(0_n, C_f)$.

#### Non-parametric Wild Bootstrap Procedure for Constructing Confidence Bands

## 2.1 A bootstrap procedure for the mean function

Let $Y_i(s)$ denotes the functional response. The estimation of the mean
function at $s$ is $\hat{\beta}(s)$. Let $q_{\alpha, N}$ be the
$(1 - \alpha)$th quantile for the random variable

$$
M_N = \max_{j \in \Lambda_N} \frac{|\hat{\beta}(s_j) - \beta(s_j)|}{\text{se}(\hat{\beta}(s_j))}
$$

where

$$
\text{se}(\hat{\beta}(s_j)) = \sqrt{ \frac{1}{n} \sum_{i=1}^n \left( \tilde{Y}_i(s_j) - \hat{\beta}(s_j) \right)^2 / (n-1) }.
$$

Since $q_{\alpha, N}$ satisfies

$$
P\left( \frac{|\hat{\beta}(s_j) - \beta(s_j)|}{\text{se}(\hat{\beta}(s_j))} \le q_{\alpha, N},\ \forall j \in \Lambda_N \right) = 1 - \alpha,
$$

then
$\hat{\beta}(s_j) \pm q_{\alpha, N} \times \text{se}(\hat{\beta}(s_j))$
is clearly the $100(1 - \alpha)\%$ simultaneous confidence band (SCB)
for $\{ \beta(s_j), j \in \Lambda_N \}$. Therefore, if we can estimate
the distribution of $M_N$, we can readily obtain an estimate of the SCB.

A wild bootstrap proposed by Chang et al. (2017) to estimate the
distribution of $M_N$ is as follows:

**(1)** Given $\tilde{Y}_1, \ldots, \tilde{Y}_n$, compute the estimated
mean function $\hat{\beta}(s_j)$ for each $j$ and calculate the residual
functions
$e_i(s_j) = \tilde{Y}_i(s_j) - \hat{\beta}(s_j),\ \forall j = 1, \ldots, N$.

**(2) Resampling step:**

- **(2.1)** For each bootstrap sample $b = 1, \ldots, B$, randomly
  generate i.i.d. multipliers $c_1(b), \ldots, c_n(b)$ from a
  distribution with mean zero and variance one.

- **(2.2)** For each $b, i$, generate the sample
  $\tilde{Y}_i^{b}(s_j) = \hat{\beta}(s_j) + c_i(b)e_i(s_j)$, calculate
  the bootstrapped mean function at each $s_j$,
  $\hat{\beta}^b(s_j) = \frac{1}{n} \sum_{i=1}^n \tilde{Y}_i^b(s_j)$,
  and the bootstrapped version of $M_N$,
  $M_N^b = \max_{j=1, \ldots, N} \frac{|\hat{\beta}^b(s_j) - \hat{\beta}(s_j)|}{\text{se}(\hat{\beta}(s_j))}$.

**(3)** Choose the $(1 - \alpha)$ sample quantile,
$\hat{q}_{\alpha, N}$, for $\{ M_N^1, \ldots, M_N^B \}$ as an
approximation of the $(1 - \alpha)$th quantile $q_{\alpha, N}$ of $M_N$.

`invSCI` provides two options for estimating the mean function at $s$,
denoted as $\hat{\beta}(s)$. If `est_mean = TRUE`, the mean function
will be estimated though using the fitted regression object. If
`est_mean = FALSE`, sample mean will be calculated. Default is `FALSE`.

1.  The **sample mean**:  
    $$
    \hat{\beta}(s) = \frac{1}{n} \sum_{i=1}^n \tilde{Y}_i(s),
    $$ where $\tilde{Y}_i(s)$ is a smoothed version of the observed
    functional response.

2.  The **fitted mean value** from a functional regression model (e.g.,
    using `mgcv::bam`).

In the wild bootstrap procedure, `invSCI` supports three types of
multiplier distributions, which is specified by `weights`:

- `"rademacher"`: $c_i \in \{-1, +1\}$ with equal probability
- `"gaussian"`: $c_i \sim \mathcal{N}(0, 1)$
- `"mammen"`: A two-point distribution with mean zero and variance one
  (see Mammen, 1993)

Default is `rademache`.

Two options are available for estimating the standard error
$\text{se}(\hat{\beta}(s_j))$, which is specified by `method_SD`:

- “regular” (empirical standard error based on residuals): $$
  \text{se}(\hat{\beta}(s_j)) = \sqrt{ \frac{1}{n} \sum_{i=1}^n \left( \tilde{Y}_i(s_j) - \hat{\beta}(s_j) \right)^2 / (n-1) }.
  $$

- “t” (bootstrap second moment-based estimator): $$
  \text{se}(\hat{\beta}(s_j)) = \sqrt{ \frac{N}{N-1} \left| \mathbb{E}_b\left[ \tilde{Y}^{b}(s_j)^2 \right] - \left( \mathbb{E}_b\left[ \tilde{Y}^{b}(s_j) \right] \right)^2 \right| },
  $$ where expectations are taken over bootstrap replicates and
  $\tilde{Y}^{b}(s_j)$ is the perturbed sample in bootstrap iteration
  $b$. The absolute value ensures numerical stability when subtracting
  large, nearly equal quantities.

Default is `t`.

The code below visualizes the **inverse confidence sets (CSs)** derived
from SCB results using the `invSCI::plot_cs()` function. The `results`
object is first converted to a tibble for easier manipulation.

The `levels = c(-7, -8, -9, -10)` argument specifies a set of
thresholds, `invSCI::plot_cs()` function estimates multiple inverse
upper excursion sets corresponding to these thresholds, and plot the
estimated inverse set, the inner confidence set, and the outer
confidence set.

``` r
results <- tibble::as_tibble(results)
plot_cs(results,levels = c(-7, -8,-9,-10), x = results$time, mu_hat = results$yhat, xlab = "", ylab = "", level_label = T, min.size = 40, palette = "Spectral", color_level_label = "black")
```

![](README_files/figure-gfm/ccds%20plot%20cs-1.png)<!-- -->

The plot demonstrate how to use SCB to find regions of s where the
estimated mean is greater than or equal to the four levels -7, -8, -9,
-10 for ccds data. The gray shaded area is the 95% SCB, the solid black
line is the estimated mean. The red horizontal line shows the inner
confidence sets (where the lower SCB is greater than the corresponding
level) that are contained in the estimated inverse set represented by
the union of the yellow and red horizontal line (where the estimated
mean is greater than the corresponding levels); the outer confidence
sets are the union of the blue, yellow and red line (where the upper SCB
is greater than the corresponding levels) and contain both the estimated
inverse sets and the inner confidence sets
