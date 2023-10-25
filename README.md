# expo2-stan

`bbr.bayes` is an extension to the `bbr` package for traceable and reproducible Bayesian modeling. Currently, `bbr.bayes` supports Stan models (powered by cmdstanr). Upcoming releases will extend Stan support and add support for NONMEM Bayes models.

This repository contains examples of the various tasks associated with using `bbr.bayes` to fit models in Stan.  It is expected that these
will evolve with time.

# Data

- Model estimation data set 
  - location : `data/derived/exdata3.csv`
  - data specification file : `data/derived/exdata3.yml`
  - exploratory data analysis : `script/info-data.Rmd`


# Scripts
- Initial model creation and submission: `script/initial-model-submission.Rmd`
- Example model diagnostics: `script/initial-model-diagnostics.Rmd`
- Iterative model development: `script/revised-models.Rmd`
- Standalone generated quantities: `script/generated-quantities.Rmd`
- Summarizing and comparing models: `script/modeling-summary.Rmd`


# Helper functions
- Helper functions for model management: `script/functions-model.R`

# Metworx Version
- metworx-22.09.1

# R Version
- 4.1.3


Copied from internal repo at 064060c92fb0f7417157cc34bfa956b814ba44a8

