## NBIMSVG
[![R build status](https://github.com/zhoumeng123456/NBIMSVG/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zhoumeng123456/NBIMSVG/actions)


## Overview
The "NBIMSVG" method propose a novel Bayesian hierarchical framework incorporating non-parametric spatial modeling and across-sample integration. It takes advantage of the non-parametric technique and develops an adaptive spatial process accommodating complex pattern discovery. A novel cross-sample bi-level shrinkage prior is further introduced for robust multi-sample SV gene detection, facilitating more effective information fusion. An efficient variational inference is developed for posterior inference ensuring computational scalability. This architecture synergistically addresses spatial complexity through adaptive pattern learning while maintaining biological interpretability. 

`NBIMSVG` is implemented as an R package.


## Installation
The package can be installed from github as follows, using R version 4.4 or above:

```{r}
install.packages("remotes")
remotes::install_github("zhoumeng123456/NBIMSVG")
```

## Example workflow
A short example workflow is shown below.

**Load packages**
```{r}
library(NBIMSVG)
```

**Generate example dataset**

```{r}
# generate example dataset by sim_create() from Package NBIMSVG
seed <- 123
result1 <- sim_create(gene_size = 20, svgene_size = 0.5, sv_mark = c(0.8, 0.8), inf_size = 0.5, seed = seed)
result2 <- sim_create(gene_size = 20, svgene_size = 0.5, sv_mark = c(0.5, 0.5), inf_size = 0.5, seed = seed + 1)
result3 <- sim_create(gene_size = 20, svgene_size = 0.5, sv_mark = c(0.5, 0.5), inf_size = 0.5, seed = seed + 2)
result4 <- sim_create(gene_size = 20, svgene_size = 0.5, sv_mark = c(0.5, 0.5), inf_size = 0.5, seed = seed + 3)
```

```{r}
dim(result1[[1]])#expression matrix 
dim(result1[[2]])#location matrix
dim(result1[[3]])#covariates
```

```{r}
## [1] 1024   20
## [1] 1024    2
## [1] 1024    6
```

```{r}
# Format input for NBIMSVG()
spelist <- list(list(result1[[1]], result1[[2]]),
                list(result2[[1]], result2[[2]]),
                list(result3[[1]], result3[[2]]),
                list(result4[[1]], result4[[2]]))
c_alpha <- list(result1[[3]],  result2[[3]],result3[[3]],result4[[3]])
```

```{r}
# Run analysis (parallel with 4 cores)
result <- NBIMSVG(spelist = spelist,c_alpha = c_alpha,num_cores = 4)
```

```{r}
# Print key outputs:
# 1. Posterior mean probabilities of each gene
print(result[[2]])
```

```{r}
 [1] 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00
[11] 4.406300e-05 4.260588e-05 4.258319e-05 4.316017e-05 4.405163e-05 4.220867e-05 4.404459e-05 4.407040e-05 4.406212e-05 4.317281e-05
```

```{r}
# 2. Binary classification of SV genes
print(result[[4]])
```

```{r}
## [1] "gene 1"  "gene 2"  "gene 3"  "gene 4"  "gene 5"  "gene 6"  "gene 7"  "gene 8"  "gene 9"  "gene 10"
```
