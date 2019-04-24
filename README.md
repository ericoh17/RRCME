
# raking-RC-ME 1.0.0

This package implements four different estimators to correct 
for bias in the presence of correlated covariate and time-to-event 
measurement error in survival analysis. The four estimators are 
(1) regression calibration (RC), (2) risk set regression calibration (RSRC),
(3) generalized raking regression calibration (GRRC), and 
(4) generalized raking naive (GRN). All of the methods assume the 
existence of a validation subset, which can be selected using
either a simple random sample or case-cohort sampling. 

Here, we demonstrate how to use the method. For more information, 
see our forthcoming paper. 

## Installation

To install and load this package in R from GitHub, run the following commands:
  
```R
install.packages("devtools")
library(devtools)
devtools::install_github("ericoh17/raking-RC-ME")
library(raking-RC-ME)
```  

## Getting Started

In this example, we will read in a simulated dataset and
corresponding validation subset. 

```R
vccc_dat <- read.csv("example_vccc_dat.csv", row.names = 1)
valid_subset <- read.csv("example_valid_subset.csv", row.names = 1)
```

The raking estimators require that the full dataset contains the true
variables for subjects that are in the validation subset.

```R
vccc_dat$time <- vccc_dat$delta <- vccc_dat$x <- NA

vccc_dat$time[vccc_dat$randomized == TRUE] <- valid_subset$time
vccc_dat$delta[vccc_dat$randomized == TRUE] <- valid_subset$delta
vccc_dat$x[vccc_dat$randomized == TRUE] <- valid_subset$x
```

Set the sampling scheme for the validation subset:
('srs' for simple random sampling and
'cc' for case-cohort sampling)
```R
sampling_scheme <- "srs"
#sampling_scheme <- "cc"
```

## Running RC

```R
RC_fit <- FitRCModel(valid_subset, vccc_dat, sampling_scheme)

RC_boot <- boot(vccc_dat, RunRCBootstrap, 
                strata = factor(vccc_dat$randomized), R = num_boot,
                dat_valid = valid_subset, sampling_type = sampling_scheme)
```

## Running RSRC

```R
RSRC_fit <- FitRSRCModel(valid_subset, vccc_dat, sampling_scheme,
                         coef(RC_fit)[1], coef(RC_fit)[2])

RSRC_boot <- boot(vccc_dat, RunRSRCBootstrap,
                  strata = factor(vccc_dat$randomized), R = num_boot,
                  dat_valid = valid_subset, sampling_type = sampling_scheme,
                  beta_x_start = coef(RC_fit)[1], beta_z_start = coef(RC_fit)[2])
```

## Running GRRC

```R
GRRC_fit <- FitRakingModel(valid_subset, vccc_dat, "RC", sampling_scheme)

GRRC_boot <- boot(vccc_dat, RunRakingBootstrap,
                  strata = factor(vccc_dat$randomized), R = num_boot,
                  dat_valid = valid_subset, mod_rake = "RC", 
                  sampling_type = sampling_scheme)
```

## Running GRN

```R
GRN_fit <- FitRakingModel(valid_subset, vccc_dat, "naive")

GRN_boot <- boot(vccc_dat, RunRakingBootstrap,
                 strata = factor(vccc_dat$randomized), R = num_boot,
                 dat_valid = valid_subset, mod_rake = "naive", 
                 sampling_type = sampling_scheme)
```

## Running naive

```R
naive_fit <- FitCoxModel(vccc_dat$time_star, vccc_dat$delta_star, vccc_dat$x_star, vccc_dat$z)

naive_boot <- boot(vccc_dat, RunNaiveBootstrap,
                   strata = factor(vccc_dat$randomized), R = num_boot)

```

## Running complete case

```R
complete_case_fit <- FitCoxModel(valid_subset$time, valid_subset$delta, valid_subset$x, valid_subset$z)

complete_case_boot <- boot(vccc_dat, RunCompleteCaseBootstrap,
                           strata = factor(vccc_dat$randomized), R = num_boot,
                           dat_valid = valid_subset, 
                           sampling_type = sampling_scheme)

```

