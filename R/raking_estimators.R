#' Calculates raking estimates
#'
#' These functions implement the raking estimators. First, either
#' a naive model or RC model is fit. Then, influence functions from
#' that model is extracted and then used as auxiliary variables
#' to improve efficiency of Horvitz-Thompson estimates.
#'
#' @param valid_dat Validation subset dataset
#'
#' @param full_dat Full dataset
#'
#' @param mod String indicating which model to get influence
#' functions from
#'
#' @param sampling_type String indicating either simple random
#' sampling or case-cohort sampling
#'
#' @return Raking model fit object
#'
#' @importFrom dplyr bind_cols
#'
#' @importFrom survey twophase calibrate svycoxph
#'
#' @rdname raking_estimators
#' @export
FitRakingModel <- function(valid_dat, full_dat, mod, sampling_type) {

  inf_func <- GetInfluenceFcn(valid_dat, full_dat, mod, sampling_type)

  dat_IF <- bind_cols(full_dat, inf_func)

  phase2_calibration <- CalibrateDesign(dat_IF, sampling_type)

  raking_fit <- FitCalibrationModel(phase2_calibration)

  return(list(raking_fit[[1]], raking_fit[[2]]))

}



GetInfluenceFcn <- function(valid_dat, full_dat, mod, sampling_type) {

  if (mod == "RC") {
    inf_func_mod <- FitRCModel(valid_dat, full_dat, sampling_type,
                               return_coef = FALSE)  # fit RC model
  } else if (mod == "naive") {
    inf_func_mod <- FitCoxModel(full_dat$time_star, full_dat$delta_star, full_dat$x_star, full_dat$z,
                                return_coef = FALSE)
  } else {
    stop("Specify 'RC' or 'naive' for calibration")
  }

  # extract influence functions from model
  inf_func_df <- data.frame(resid(inf_func_mod, "dfbeta"))
  colnames(inf_func_df) <- paste("if", 1:2, sep = "")

  return(inf_func_df)
}



CalibrateDesign <- function(IF_dat, sampling_type) {

  if (sampling_type == "cc") {

    IF_design <- twophase(id = list(~id, ~id), subset = ~randomized,
                          strata = list(NULL, ~delta_star), data = IF_dat)
    IF_raking <- calibrate(IF_design, phase = 2, formula = ~if1+if2+delta_star, calfun = "raking")

  } else if (sampling_type == "srs") {

    IF_design <- twophase(id = list(~id, ~id), subset = ~randomized, data = IF_dat)
    IF_raking <- calibrate(IF_design, phase = 2, formula = ~if1+if2, calfun = "raking")

  }

}



FitCalibrationModel <- function(raking_dat) {

  raking_mod <- svycoxph(Surv(time, delta) ~ x + z, design = raking_dat)

  return(list(raking_mod$coef[1], raking_mod$coef[2]))

}

