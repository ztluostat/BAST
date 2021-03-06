# BAST
Code for Bayesian Additive Regression Spanning Trees, a tool for nonparametric regression on complex constrained domains.

Reference:

[Luo, Z. T., Sang, H., & Mallick, B. (2021) BAST: Bayesian Additive Regression Spanning Trees for Complex Constrained Domain. *Advances in Neural Information Processing Systems 34 (NeurIPS 2021)*](https://proceedings.neurips.cc/paper/2021/hash/00b76fddeaaa7d8c2c43d504b2babd8a-Abstract.html)

### Files

`BASTDemo.R`: Demonstration code of BAST

`AralSea.R`: Preprocessing chlorophyll data in Aral Sea for demonstration purposes, following [Niu et al. (2019)](https://github.com/mu2013/Intrinsic-GP-on-complex-constrained-domain)

`BASTFun.R`: Implementation and API documentations for main functions of BAST

`ComplexDomainFun.R`: Utility functions for complex constrained domains

`FEMFun.R`: Utility functions for triangular meshes, modified from `R` package [fdaPDE](https://cran.r-project.org/web/packages/fdaPDE/index.html)

`aral_data.RData`: Preprocessed data from `AralSea.R`

### Dependencies

The code depends on the following `R` packages: `igraph`, `fields`, `FNN`, `mgcv`, `fdaPDE`, `gamair`.
Please make sure they are installed before running the demo code.

### Future development

A more efficient version using `Rcpp` is under active development.
