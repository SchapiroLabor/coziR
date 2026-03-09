# coziR

`coziR` is an R implementation of the COZI neighborhood preference method, mirroring the same conceptual steps as [cozipy](https://github.com/SchapiroLabor/COZIpy), 
introduced in [Schiller at al. bioRxiv 2025](https://doi.org/10.1101/2025.03.31.646289).

1. Build a neighborhood graph (`knn`, `radius`, `delaunay`)
2. Quantify interactions (neighbors) across the tissue 
3. Conditionally normalize interaction count
4. Compute permutation-based expected values and z-scores
5. Return conditional cell ratios and z-score matrices

## Installation

```r
# install from GitHub
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("SchapiroLabor/coziR")
```

Or clone the repository and install from source:

```bash
git clone https://github.com/chiaraschiller/coziR.git
cd coziR
```

```r
install.packages(".", repos = NULL, type = "source")
```

## Quick example

```r
library(coziR)

coords <- matrix(c(
  0, 0,
  1, 0,
  0, 1,
  1, 1,
  0.5, 0.5,
  1.5, 0.5
), ncol = 2, byrow = TRUE)

labels <- c("T", "T", "B", "B", "T", "B")

res <- run_cozi(
  coords = coords,
  labels = labels,
  nbh_def = "knn",
  n_neighbors = 3,
  n_permutations = 100,
  random_state = 42
)

res$cond_ratio
res$zscore
```

## Dotplot for KNN results

```r
library(ggplot2)

p <- cozi_dotplot(
  res,
  title = "KNN neighborhood preference"
)

print(p)
```

## Repository guide

Use this order if your goal is to run the method quickly:
1. Install the package (section above)
2. Run the tutorial in `tutorial/coziR_tutorial.Rmd`
3. Adapt the same calls in your own analysis script

Repository structure:

- `DESCRIPTION`, `NAMESPACE`: R package metadata and exported functions.
- `NEWS.md`: package version history and user-facing changes across releases.
- `R/`: core implementation.
- `R/neighbors.R`: neighborhood graph builders (`knn`, `radius`, `delaunay`).
- `R/nep_cozi.R`: main NEP/COZI workflow (`nep_analysis`, `run_cozi`).
- `R/plotting.R`: plotting helper (`cozi_dotplot`).
- `R/internal.R`: internal helper utilities.
- `man/`: function documentation (`.Rd` help pages).
- `tests/test-basic.R`: basic checks/examples for expected behavior.
- `tutorial/coziR_tutorial.Rmd`: step-by-step runnable tutorial.
- `tutorial/tutorial_data/test_sim_data.csv`: example input used by the tutorial.

## Tutorial and interpretation

- For hands-on execution, start with `tutorial/coziR_tutorial.Rmd` and run it top to bottom.
- For biological/statistical interpretation of outputs (for example conditional ratios and z-scores), use the associated manuscript as the primary reference. We are happy to provide this package but ask the user to carefully read the manuscript to be aware of how to propoerly interpret the results.
- For release-by-release package updates, see `NEWS.md`.

## Main API

- `knn_graph(coords, n_neighbors = 6)`
- `radius_graph(coords, radius)`
- `delaunay_graph(coords)`
- `nep_analysis(adj, labels, n_permutations = 1000, random_state = NULL, return_df = TRUE, fixed_type = NULL)`
- `run_cozi(coords, labels, nbh_def = c("knn", "radius", "delaunay"), n_neighbors = 6, radius = 0.2, n_permutations = 100, random_state = NULL, fixed_type = NULL)`
- `cozi_dotplot(results, title = "COZI Dotplot")`

## References

Schiller et al., bioRxiv (2025) "Comparison and Optimization of Cellular Neighbor Preference Methods for Quantitative Tissue Analysis", 
doi: https://doi.org/10.1101/2025.03.31.646289
