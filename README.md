# Dominance-model analyses

This repository contains the R scripts used to fit, compare, validate, and diagnose four latent dominance models on dyadic interaction data from `DomArchive`:

- **LBT**: latent Bradley--Terry model
- **LID**: latent intransitive dominance model
- **DABT**: LBT with dyad-specific random effects
- **DALID**: LID with dyad-specific random effects

The models are fitted with `RTMB`, using Laplace approximation for the latent variables and `nlminb()` for optimization.

## Scripts

### `FourModelsComparisonInference.R`
Main analysis script. It:

- filters and prepares count matrices from `DomArchive`;
- fits all four models using multiple starting values;
- extracts parameter estimates, likelihoods, AIC values, fitted probabilities, latent positions, and dyad effects;
- computes gradient and fixed-parameter Hessian diagnostics;
- saves model results and produces tables and figures.

Typical output is written to `DomArchiveResults/`.

### `BootstrappingOnDomArchive.R`
Parametric bootstrap analyses. It:

- simulates data from fitted models;
- computes bootstrap likelihood-ratio tests for the four nested comparisons;
- estimates bootstrap distributions and Monte Carlo p-values;
- performs parameter bootstrapping for selected datasets;
- creates bootstrap summaries, diagnostic plots, and LaTeX tables.

### `LRTGridAnalysis.R`
Simulation study for the null distribution of likelihood-ratio statistics. It:

- simulates complete dyadic systems for grids of group size $n$ and interactions per dyad $s$;
- compares empirical LRT distributions with the reference distribution $0.5\chi^2_0+0.5\chi^2_1$;
- saves simulation results, summaries, histograms, QQ plots, and PP plots.

### `MCMC.R`
Conditional MCMC diagnostics for the latent variables. It:

- samples from the conditional latent distribution given fitted fixed parameters;
- uses adaptive random-walk Metropolis--Hastings;
- computes acceptance rates, effective sample sizes, integrated autocorrelation times, skewness, and excess kurtosis;
- produces marginal-distribution plots and diagnostic tables.

## Main dependencies

```r
install.packages(c(
  "dplyr", "purrr", "tidyr", "ggplot2", "tidyverse",
  "plotly", "patchwork", "scales", "gt", "knitr",
  "kableExtra", "readr", "stringr", "rlang"
))

install.packages("RTMB")
install.packages("DomArchive")
```

## Recommended workflow

1. Run `FourModelsComparisonInference.R` to fit the models and create the main result files.
2. Run `BootstrappingOnDomArchive.R` for bootstrap inference on selected datasets.
3. Run `LRTGridAnalysis.R` for the simulation-based calibration study.
4. Run `MCMC.R` to assess the conditional latent distributions and the adequacy of the Laplace approximation.

## Notes

- Several file paths use `~` or project-relative folders and may need to be changed for a new system.
- The scripts contain both reusable functions and exploratory analysis blocks. They are intended to document the thesis workflow rather than function as a standalone R package.
- Bootstrap and MCMC sections are computationally expensive. Reduce the number of replicates or iterations when testing the code.
- Reproducibility seeds are set in the main simulation sections, commonly to `123`.
