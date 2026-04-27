# Optimal Versus Naive Diversification

Replication of DeMiguel, Garlappi, and Uppal (2009, RFS), "Optimal Versus Naive Diversification: How Inefficient is the 1/N Portfolio Strategy?," *The Review of Financial Studies*, 22(5), 1915-1953.

## What is replicated

| Output | Script | Description |
|--------|--------|-------------|
| Table 3 | `table_3.R` | Out-of-sample Sharpe ratios for 14 portfolio strategies across 6 empirical datasets |

Results compare the naive 1/N benchmark against 13 optimizing strategies using risky-only portfolios, a 120-month rolling estimation window, and risk aversion gamma = 1.

## How to run

```r
source("table_3.R")
```

## Methods

### Portfolio strategies (Table 3)

| # | Strategy | Abbreviation | Description |
|---|----------|-------------|-------------|
| 0 | 1/N equal weight | ew | Benchmark: allocate 1/N to each asset |
| 1 | Sample mean-variance | mv | Unconstrained plug-in Markowitz |
| 2 | Bayes-Stein | bs | Jorion (1986) shrinkage of means |
| 3 | Bayesian Data-and-Model | dm | Pastor-Stambaugh (2000) with sigma_alpha = 1% |
| 4 | Minimum-variance | min | Ignore means, minimize variance only |
| 5 | Value-weighted market | vw | 100% in the market index |
| 6 | MacKinlay-Pastor | mp | Missing-factor model (2000) |
| 7 | MV constrained | mv-c | Mean-variance with shortsale constraints |
| 8 | Bayes-Stein constrained | bs-c | Bayes-Stein with shortsale constraints |
| 9 | Min-variance constrained | min-c | Minimum-variance, long-only |
| 10 | Generalized min-variance | g-min-c | Jagannathan-Ma (2003) with lb = 1/(2N) |
| 11 | Kan-Zhou three-fund | mv-min | Optimal combo of tangency + min-variance + risk-free |
| 12 | EW-min mixture | ew-min | Optimal combo of 1/N + min-variance |

### Evaluation metrics

- **Out-of-sample Sharpe ratio** -- mean / std of monthly excess returns from rolling-window portfolios
- **In-sample Sharpe ratio** -- MV-optimal SR computed on the full out-of-sample period (upper bound)
- **Jobson-Korkie p-value** -- tests equality of Sharpe ratios vs. 1/N benchmark, with Memmel (2003) correction

### Key parameters

| Parameter | Value | Meaning |
|-----------|-------|---------|
| Estimation window (M) | 120 months | 10-year rolling window |
| Risk aversion (gamma) | 1 | Mean-variance utility |
| Window type | Rolling | Fixed-length sliding window |
| Portfolio returns | Risky-only | Weights normalized by \|sum(w)\| |
| Sigma (MV, BS) | Unbiased | (T-1)/(T-N-2) * cov() |
| Sigma (min, KZ, PM) | MLE | (T-1)/T * cov() |

## Repository structure

```
data/
  SPSectors.rds         <- S&P 500 sector portfolios (01/1981-12/2002, N=11)
  Industry.rds          <- 10 industry portfolios (07/1963-11/2004, N=11)
  International.rds     <- 8 country indexes (01/1970-07/2001, N=9)
  MKT_SMB_HML.rds       <- SMB and HML portfolios (07/1963-11/2004, N=3)
  FF_1factor.rds        <- 20 size/BM portfolios, 1 factor (07/1963-11/2004, N=21)
  FF_4factor.rds        <- 20 size/BM portfolios, 4 factors (07/1963-11/2004, N=24)
table_3.R               <- main replication script
```

## Data sources

- S&P Sectors: Roberto Wessels
- Industry, FF portfolios: Ken French's website
- International: MSCI
- MKT/SMB/HML: Ken French's website

## References

- DeMiguel, V., Garlappi, L., and Uppal, R. (2009). Optimal versus naive diversification: How inefficient is the 1/N portfolio strategy? *The Review of Financial Studies*, 22(5), 1915-1953.
