# calipeRmatching

An R package implementing of the caliper matching estimator of the Average Treatment Effect (ATE) and the Average Treatment Effect on the Treated (ATT) of a binary treatment. These effects are estimated by matching on the
- known or 
- estimated parametric
propensity score. See [KPV23](https://arxiv.org/abs/2304.08373) for details.


This package is a wrapper around the C library [`caliper_matching`](https://github.com/kmmate/caliper_matching).

## Installation

From R:

```R
install.packages("remotes")
remotes::install_github("kmmate/calipeRmatching")
```

## References

Kormos, V. d. Pas, V. d. Vaart (2023): Asymptotics of Caliper Matching Estimators for Average Treatment Effects, https://arxiv.org/abs/2304.08373

