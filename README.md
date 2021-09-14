# scDesign2
A transparent simulator that generates high-fidelity single-cell gene expression count data with gene correlations captured

## Installation
```r
if(!require(devtools)) install.packages("devtools")
library(devtools)

devtools::install_github("JSB-UCLA/scDesign2")
```
If errors occur due to the `curl` package, [this discussion](https://stackoverflow.com/questions/20923209/problems-installing-the-devtools-package) on stackoverflow could be helpful.

If dependency R packages are updated when installing the `devtools` R package, try closing and restarting the R session before proceeding.

The scDesign2 R package do not have many dependencies (for a complete list, please refer to the DESCRIPTION file). However, it does depend on the `parallel` package, where we use the `mclapply()` function for the parallel model fitting of multiple cell types. Unfortunately, `mclapply()`'s parallelization can only be executed on Unix-like operating systems but not on Windows. Therefore, the model fitting of multiple cell types will take longer on Windows. We will try to fix this issue in the future.

## Tutorials
This [tutorial](https://htmlpreview.github.io/?https://github.com/JSB-UCLA/scDesign2/blob/master/vignettes/scDesign2.html) introduces how to use the scDesign2 R package, and can be directly used when trustworthy cell type information is available.

If not, then this [tutorial](https://htmlpreview.github.io/?https://github.com/JSB-UCLA/scDesign2/blob/master/vignettes/preprocess_clustering.html) can be read, which introduces how cell clustering can be performed and evaluated before using the scDesign2 method.

## Citation
Sun, T., Song, D., Li, W.V. et al. scDesign2: a transparent simulator that generates high-fidelity single-cell gene expression count data with gene correlations captured. *Genome Biol* **22,** 163 (2021). https://doi.org/10.1186/s13059-021-02367-2
