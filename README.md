# Causal Inference Software Review Summary

- **Traced the evolution** of causal inference from Wright's path diagrams in 1920 to the introduction of robust methods like TMLE by Van Der Laan and Rubin in 2006.
- **Assessed R packages** gfoRmula, ltmle, and Weightit, each implementing G-computation, TMLE, and MSM, for causal effect estimation.
- **Identified barriers** in biostatistics knowledge translation, such as lack of expertise and user-friendly software.
- **Conducted simulations** to compare methods across 1,000 iterations, summarizing performance in terms of bias, standard error, and confidence interval coverage.
- **Found TMLE** to exhibit statistical power and robustness, despite its higher standard error in binary outcome data.
- **Recommended the ltmle package** for longitudinal causal analysis due to its documentation and robust performance.

## Future steps:
- Additional simulation for Bayesian Marginal Structural Model
- Finishing manuscript for publication
