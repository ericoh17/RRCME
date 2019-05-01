
# RRCME

This package implements four different estimators to correct 
for bias in time-to-event data in the presence of correlated 
covariate and censored event time measurement error. The four 
estimators are (1) regression calibration (RC), 
(2) risk set regression calibration (RSRC),
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
devtools::install_github("ericoh17/RRCME")
library(RRCME)
```  

## Getting Started

Next we load the boot package to run the bootstrap which
is required to obtain standard errors and set the 
number of bootstrap replicates desired:

```R
library(boot)
num_boot <- 250
``` 

In this example, we will read in a simulated dataset and
corresponding validation subset. 

```R
data(example_full_dat_CC)
data(example_valid_subset_CC)
```

The raking estimators require that the full dataset contains the true
variables for subjects that are in the validation subset. The 
implementation also requires that the full dataset contains a column 
called 'randomized' that is TRUE for those subjects selected in 
the validation subset and FALSE otherwise.

```R
full_dat$time <- full_dat$delta <- full_dat$x <- NA

full_dat$time[full_dat$randomized == TRUE] <- valid_subset$time
full_dat$delta[full_dat$randomized == TRUE] <- valid_subset$delta
full_dat$x[full_dat$randomized == TRUE] <- valid_subset$x
```

Set the sampling scheme for the validation subset
('srs' for simple random sampling and
'cc' for case-cohort sampling):
```R
sampling_scheme <- "cc"
#sampling_scheme <- "srs"
```

## Run RC

The function to run RC is the `FitRCModel` function.
We obtain standard errors by calling the `boot` function with the 
`RunRCBootstrap` function. 

```R
RC_fit <- FitRCModel(valid_subset, full_dat, sampling_scheme)

RC_boot <- boot(full_dat, RunRCBootstrap, 
                strata = factor(full_dat$randomized), R = num_boot,
                dat_valid = valid_subset, sampling_type = sampling_scheme)
```

## Run RSRC

The function to run RSRC is the `FitRSRCModel` function. 
We obtain standard errors by calling the `boot` function with the 
`RunRSRCBootstrap` function. 

```R
RSRC_fit <- FitRSRCModel(valid_subset, full_dat, sampling_scheme,
                         coef(RC_fit)[1], coef(RC_fit)[2])

RSRC_boot <- boot(full_dat, RunRSRCBootstrap,
                  strata = factor(full_dat$randomized), R = num_boot,
                  dat_valid = valid_subset, sampling_type = sampling_scheme,
                  beta_x_start = coef(RC_fit)[1], beta_z_start = coef(RC_fit)[2])
```

## Run GRRC

The function to run GRRC is the `FitRakingModel` function with the 'RC' argument. 
We obtain standard errors by calling the `boot` function with the 
`RunRakingBootstrap` function. 

```R
GRRC_fit <- FitRakingModel(valid_subset, full_dat, "RC", sampling_scheme)

GRRC_boot <- boot(full_dat, RunRakingBootstrap,
                  strata = factor(full_dat$randomized), R = num_boot,
                  dat_valid = valid_subset, mod_rake = "RC", 
                  sampling_type = sampling_scheme)
```

## Run GRN

The function to run GRRC is the `FitRakingModel` function with the 'naive' argument. 
We obtain standard errors by calling the `boot` function with the 
`RunRakingBootstrap` function. 

```R
GRN_fit <- FitRakingModel(valid_subset, full_dat, "naive", sampling_scheme)

GRN_boot <- boot(full_dat, RunRakingBootstrap,
                 strata = factor(full_dat$randomized), R = num_boot,
                 dat_valid = valid_subset, mod_rake = "naive", 
                 sampling_type = sampling_scheme)
```


