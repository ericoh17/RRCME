#' Runs RRCME procedure
#'
#' Runs some combination of Regression
#' Calibration, Risk set Regression
#' Calibration, Generalized Raking
#' Regression Calibration, and
#' Generalized Raking Naive as specified
#' by the estimators input.
#' Log hazard ratios and corresponding
#' standard errors are returned.
#'
#' @param full_dat Data for full cohort
#' containing error-prone variables
#'
#' @param valid_dat Validation subset
#' from full cohort containing
#' error-free variables
#'
#' @param num_boot Number of bootstrap
#' replicates
#'
#' @param sampling_type Sampling design
#' for validation subset
#'
#' @param estimators Vector of strings
#' indicating which estimators to run
#'
#' @return Matrix of estimates and
#' standard errors
#'
#' @importFrom boot boot
#'
#' @rdname get_estimates
#' @export
CalcCorrectedlogHR <- function(full_dat, valid_dat,
                               num_boot = 250, sampling_type,
                               estimators = c("RC", "RSRC", "GRRC", "GRN")) {

  if (!(sampling_type %in% c("cc", "srs"))) {
    stop("Supported sampling designs are case-cohort (denoted `cc`) and
         simple random sampling (denoted `srs`)")
  }

  if (any(estimators %in% c("RC", "RSRC", "GRRC", "GRN") == FALSE)) {
    stop("Only four estimators are supported:
         Regression Calibration (`RC`),
         Risk Set Regression Calibration (`RSRC`),
         Generalized Raking Regression Calibration (`GRRC`),
         Generalized Raking Naive (`GRN`)")
  }

  if (any(valid_dat$id %in% full_dat$id == FALSE)) {
    stop("Validation data must be a subset of the full data")
  }

  if (num_boot < 100) {
    stop("Number of bootstrap replicates should be at least
         100 to ensure proper confidence intervals")
  }

  full_dat$time[full_dat$randomized == TRUE] <- valid_dat$time
  full_dat$delta[full_dat$randomized == TRUE] <- valid_dat$delta
  full_dat$x[full_dat$randomized == TRUE] <- valid_dat$x
  full_dat$total_y_err[full_dat$randomized == TRUE] <- valid_dat$total_y_err

  if (sampling_type == "cc") {

    full_dat$cc_strata[full_dat$randomized == FALSE] <- 1
    full_dat$cc_strata[full_dat$randomized == TRUE & full_dat$delta_star == 1] <- 2
    full_dat$cc_strata[full_dat$randomized == TRUE & full_dat$delta_star == 0] <- 3

  }

  results_lst <- list()

  if ("RC" %in% estimators) {

    message("Running RC")

    RC_fit <- FitRCModel(valid_dat, full_dat,
                         sampling_type,
                         return_coef = TRUE)

    message("Calculating standard errors for RC")

    if (sampling_type == "cc") {

      # CC RC bootstrap
      RC_boot <- boot(full_dat,
                      RunRCBootstrap,
                      strata = factor(full_dat$cc_strata),
                      R = num_boot,
                      sampling_type = sampling_type)

    } else if (sampling_type == "srs") {

      # SRS RC bootstrap
      RC_boot <- boot(full_dat,
                      RunRCBootstrap,
                      strata = factor(full_dat$randomized),
                      R = num_boot,
                      sampling_type = sampling_type)
    }

    RC_mat <- matrix(c(RC_fit[[1]],
                       sd(RC_boot$t[,1], na.rm = TRUE),
                       RC_fit[[2]],
                       sd(RC_boot$t[,2], na.rm = TRUE)),
                     nrow = 2, byrow = TRUE)

    rownames(RC_mat) <- c("beta_x", "beta_z")
    colnames(RC_mat) <- c("Estimate", "Standard Error")

    results_lst[["RC"]] <- as.data.frame(RC_mat)
  }

  if ("RSRC" %in% estimators) {

    message("Running RSRC")

    RSRC_fit <- FitRSRCModel(valid_dat, full_dat,
                             sampling_type,
                             RC_fit[[1]],
                             RC_fit[[2]])

    message("Calculating standard errors for RSRC")

    if (sampling_type == "cc") {

      # CC RSRC bootstrap
      RSRC_boot <- boot(full_dat,
                        RunRSRCBootstrap,
                        strata = factor(full_dat$cc_strata),
                        R = num_boot,
                        sampling_type = sampling_type,
                        beta_x_start = RC_fit[[1]],
                        beta_z_start = RC_fit[[2]])

    } else if (sampling_type == "srs") {

      # SRS RSRC bootstrap
      RSRC_boot <- boot(full_dat,
                        RunRSRCBootstrap,
                        strata = factor(full_dat$randomized),
                        R = num_boot,
                        sampling_type = sampling_type,
                        beta_x_start = RC_fit[[1]],
                        beta_z_start = RC_fit[[2]])

    }

    RSRC_mat <- matrix(c(RSRC_fit[[1]],
                         sd(RSRC_boot$t[,1], na.rm = TRUE),
                         RSRC_fit[[2]],
                         sd(RSRC_boot$t[,2], na.rm = TRUE)),
                       nrow = 2, byrow = TRUE)

    rownames(RSRC_mat) <- c("beta_x", "beta_z")
    colnames(RSRC_mat) <- c("Estimate", "Standard Error")

    results_lst[["RSRC"]] <- as.data.frame(RSRC_mat)
  }

  if ("GRRC" %in% estimators) {

    message("Running GRRC")

    GRRC_fit <- FitRakingModel(valid_dat, full_dat,
                               mod_rake = "RC",
                               sampling_type)

    message("Calculating standard errors for GRRC")

    if (sampling_type == "cc") {

      # CC GRRC bootstrap
      GRRC_boot <- boot(full_dat,
                        RunRakingBootstrap,
                        strata = factor(full_dat$cc_strata),
                        R = num_boot,
                        mod_rake = "RC",
                        sampling_type = sampling_type)

    } else if (sampling_type == "srs") {

      # SRS GRRC bootstrap
      GRRC_boot <- boot(full_dat,
                        RunRakingBootstrap,
                        strata = factor(full_dat$randomized),
                        R = num_boot,
                        mod_rake = "RC",
                        sampling_type = sampling_type)

    }

    GRRC_mat <- matrix(c(GRRC_fit[[1]],
                         sd(GRRC_boot$t[,1], na.rm = TRUE),
                         GRRC_fit[[2]],
                         sd(GRRC_boot$t[,2], na.rm = TRUE)),
                       nrow = 2, byrow = TRUE)

    rownames(GRRC_mat) <- c("beta_x", "beta_z")
    colnames(GRRC_mat) <- c("Estimate", "Standard Error")

    results_lst[["GRRC"]] <- as.data.frame(GRRC_mat)
  }

  if ("GRN" %in% estimators) {

    message("Running GRN")

    GRN_fit <- FitRakingModel(valid_dat, full_dat,
                              mod_rake = "naive",
                              sampling_type)

    message("Calculating standard errors for GRN")

    if (sampling_type == "cc") {

      # CC GRN bootstrap
      GRN_boot <- boot(full_dat,
                       RunRakingBootstrap,
                       strata = factor(full_dat$cc_strata),
                       R = num_boot,
                       mod_rake = "naive",
                       sampling_type = sampling_type)

    } else if (sampling_type == "srs") {

      # SRS GRN bootstrap
      GRN_boot <- boot(full_dat,
                       RunRakingBootstrap,
                       strata = factor(full_dat$randomized),
                       R = num_boot,
                       mod_rake = "naive",
                       sampling_type = sampling_type)

    }

    GRN_mat <- matrix(c(GRN_fit[[1]],
                        sd(GRN_boot$t[,1], na.rm = TRUE),
                        GRN_fit[[2]],
                        sd(GRN_boot$t[,2], na.rm = TRUE)),
                      nrow = 2, byrow = TRUE)

    rownames(GRN_mat) <- c("beta_x", "beta_z")
    colnames(GRN_mat) <- c("Estimate", "Standard Error")

    results_lst[["GRN"]] <- as.data.frame(GRN_mat)
  }

  message("Done!")
  return(results_lst)

}


