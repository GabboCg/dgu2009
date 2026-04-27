# Optimal Versus Naive Diversification: How Inefficient is the 1/N Portfolio Strategy?

## Overview

R replication of **Table 3** (Sharpe ratios for empirical data) from:

> DeMiguel, V., Garlappi, L., & Uppal, R. (2009). Optimal Versus Naive Diversification: How Inefficient is the 1/N Portfolio Strategy? *The Review of Financial Studies*, 22(5), 1915-1953.

The original replication code (MATLAB) is available from the authors. This project translates it to R.

## Requirements

- R (>= 4.0)
- `quadprog` package

```r
install.packages("quadprog")
```

## Usage

```bash
Rscript table_3.R
```

## Project Structure

```
.
├── table_3.R               # Main replication script
├── data/
│   ├── SPSectors.rds       # S&P 500 sector portfolios (01/1981-12/2002, N=11)
│   ├── Industry.rds        # 10 industry portfolios (07/1963-11/2004, N=11)
│   ├── International.rds   # 8 country indexes (01/1970-07/2001, N=9)
│   ├── MKT_SMB_HML.rds     # SMB and HML portfolios (07/1963-11/2004, N=3)
│   ├── FF_1factor.rds      # 20 size/BM portfolios, 1 factor (07/1963-11/2004, N=21)
│   └── FF_4factor.rds      # 20 size/BM portfolios, 4 factors (07/1963-11/2004, N=24)
├── README.md
└── CLAUDE.md
```

## Portfolio Strategies

The script implements 14 asset-allocation models from the paper:

| # | Strategy | Abbreviation |
|---|----------|-------------|
| 0 | 1/N equal weight (benchmark) | ew |
| 1 | Sample-based mean-variance | mv |
| 2 | Bayes-Stein | bs |
| 3 | Bayesian Data-and-Model | dm |
| 4 | Minimum-variance | min |
| 5 | Value-weighted market | vw |
| 6 | MacKinlay-Pastor missing-factor | mp |
| 7 | MV with shortsale constraints | mv-c |
| 8 | Bayes-Stein with shortsale constraints | bs-c |
| 9 | Min-variance with shortsale constraints | min-c |
| 10 | Min-variance with generalized constraints | g-min-c |
| 11 | Kan-Zhou three-fund | mv-min |
| 12 | Mixture of 1/N and min-variance | ew-min |

## Methodology

- **Estimation window**: M = 120 months (rolling)
- **Risk aversion**: gamma = 1
- **Portfolio returns**: Risky-only (weights normalized by |sum(w)|)
- **Performance metric**: Out-of-sample Sharpe ratio
- **Statistical test**: Jobson-Korkie (1981) with Memmel (2003) correction

## Data Sources

- S&P Sectors: Roberto Wessels
- Industry, FF portfolios: Ken French's website
- International: MSCI
- MKT/SMB/HML: Ken French's website
