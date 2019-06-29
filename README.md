
# RRCME

This package implements four different estimators to correct 
for bias in time-to-event data in the presence of correlated 
covariate and censored event time measurement error. The four 
estimators are 

* Regression Calibration (`RC`), 
* Risk Set Regression Calibration (`RSRC`)
* Generalized Raking Regression Calibration (`GRRC`), and 
* Generalized Raking Naive (`GRN`) 

All of the methods assume the existence of a 
validation subset, which can be selected using
either a simple random sample or case-cohort sampling. 

Here, we demonstrate how to use the method. For more information, 
see our our [paper](https://arxiv.org/abs/1905.08330). 

## Installation

To install and load this package in R from GitHub, run the following commands:
  
```{r}
install.packages("devtools")
library(devtools)
devtools::install_github("ericoh17/RRCME")
```

## Load Package
```R
library(RRCME)
```

## Getting Started

In this example, we will read in a simulated dataset and
corresponding validation subset. The validation subset was
selected retrospectively using a case-cohort design.

```{r}
data(example_full_dat)

# case-cohort validation subset
data(example_valid_subset_cc)
```
The above commands load in two dataframes. 

`full_dat` is a dataframe with the following columns:

* `id`: example identification number for each subject
* `x_star`: error-prone covariate
* `z`: error-free covariate
* `time_star`: error-prone censored event time
* `delta_star`: error-prone event indicator
* `randomized`: boolean indicating whether a subject was selected to be 
in the validation subset or not

`valid_subset` is a dataframe with the following columns:

* `id`: example identification number for each subject
* `x`: true, error-free version of error-prone covariate
* `x_star`: error-prone covariate
* `z`: error-free covariate
* `time`: true, error-free version of error-prone censored event time
* `time_star`: error-prone censored event time
* `delta`: true, error-free version of error-prone event indicator
* `delta_star`: error-prone event indicator
* `total_y_err`: `time_star - time`

The columns `id`, `x_star`, `time_star`, and `delta_star` must be 
contained in both the full dataset as well as the validation subset. 
In addition, the full dataset must contain a column 
called `randomized` that is TRUE for those subjects selected in 
the validation subset and FALSE otherwise.

## Settings 

Set the sampling scheme for the validation subset
('srs' for simple random sampling and
'cc' for case-cohort sampling):

```R
sampling_scheme <- "cc"
```

Set the number of bootstrap replicates to
perform to calculate standard errors:

```R
num_boot <- 250
```

## Run estimators

First, we create a vector with strings representing
the estimators to be run. Any combination of
the `RC`, `RSRC`, `GRRC`, and `GRN` can be run.

```R
estimators <- c("RC", "RSRC", "GRRC", "GRN")
```

Then we call the function `CalcCorrectedlogHR`
with the above arguments.

```R
res <- CalcCorrectedlogHR(full_dat, valid_subset,
                          num_boot, sampling_scheme
                          estimators)
```

`res` is a list containing dataframes for each of the estimators
and contains the estimates and their corresponding
standard errors. Confidence intervals can then be 
calculated using your favorite method. 

