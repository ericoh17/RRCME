
CalcCorrectedlogHR <- function(full_dat, valid_dat, 
                               num_boot = 250, sampling_type,
                               estimators = c("RC", "RSRC", "GRRC", "GRN")) {
  
  ## check time_star, delta_star, x_star exist
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
  
  if ("RC" %in% estimators) {
    
    RC_fit <- FitRCModel(valid_dat, full_dat, 
                         sampling_scheme, 
                         return_coef = TRUE)
    
    if (sampling_type == "cc") {
      
      # CC RC bootstrap
      RC_boot <- boot(full_dat, 
                      RunRCBootstrap, 
                      strata = factor(full_dat$cc_strata), 
                      R = num_boot,
                      sampling_type = sampling_scheme)
      
    } else if (sampling_type == "srs") {
      
      # SRS RC bootstrap
      RC_boot <- boot(full_dat, 
                      RunRCBootstrap, 
                      strata = factor(full_dat$randomized), 
                      R = num_boot,
                      sampling_type = sampling_scheme)
    }
  }

  if ("RSRC" %in% estimators) {
    
    RSRC_fit <- FitRSRCModel(valid_dat, full_dat, 
                             sampling_scheme,
                             RC_fit[[1]], 
                             RC_fit[[2]])
    
    if (sampling_type == "cc") {
      
      # CC RSRC bootstrap
      RSRC_boot <- boot(full_dat, 
                        RunRSRCBootstrap,
                        strata = factor(full_dat$cc_strata), 
                        R = num_boot,
                        sampling_type = sampling_scheme,
                        beta_x_start = RC_fit[[1]], 
                        beta_z_start = RC_fit[[2]])
      
    } else if (sampling_type == "srs") {
      
      # SRS RSRC bootstrap
      RSRC_boot <- boot(full_dat, 
                        RunRSRCBootstrap,
                        strata = factor(full_dat$randomized), 
                        R = num_boot,
                        sampling_type = sampling_scheme,
                        beta_x_start = RC_fit[[1]], 
                        beta_z_start = RC_fit[[2]])
      
    }
  }

  if ("GRRC" %in% estimators) {
    
    GRRC_fit <- FitRakingModel(valid_dat, full_dat, 
                               mod_rake = "RC", 
                               sampling_scheme)
    
    if (sampling_type == "cc") {
      
      # CC GRRC bootstrap
      GRRC_boot <- boot(full_dat, 
                        RunRakingBootstrap,
                        strata = factor(full_dat$cc_strata), 
                        R = num_boot,
                        mod_rake = "RC", 
                        sampling_type = sampling_scheme)
      
    } else if (sampling_type == "srs") {
      
      # SRS GRRC bootstrap
      GRRC_boot <- boot(full_dat, 
                        RunRakingBootstrap,
                        strata = factor(full_dat$randomized), 
                        R = num_boot,
                        mod_rake = "RC", 
                        sampling_type = sampling_scheme)
      
    }
  }

  if ("GRN" %in% estimators) {
    
    GRN_fit <- FitRakingModel(valid_dat, full_dat, 
                              mod_rake = "naive", 
                              sampling_scheme)
    
    if (sampling_type == "cc") {
      
      # CC GRN bootstrap
      GRN_boot <- boot(full_dat, 
                       RunRakingBootstrap,
                       strata = factor(full_dat$cc_strata), 
                       R = num_boot,
                       mod_rake = "naive", 
                       sampling_type = sampling_scheme)
      
    } else if (sampling_type == "srs") {
      
      # SRS GRN bootstrap
      GRN_boot <- boot(full_dat, 
                       RunRakingBootstrap,
                       strata = factor(full_dat$randomized), 
                       R = num_boot,
                       mod_rake = "naive", 
                       sampling_type = sampling_scheme)
      
    }
  }

}


