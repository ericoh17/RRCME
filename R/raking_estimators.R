#' Calculates raking estimates
#'
#' These functions implement the raking estimators. First, either
#' a naive model or RC model is fit. Then, influence functions from
#' that model is extracted and then used as auxiliary variables
#' to improve efficiency of Horvitz-Thompson estimates.
#'
#' @param valid_dat Validation subset dataset
#'
#' @param dat_sim Full dataset
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
FitRakingModel <- function(valid_dat, dat_sim, mod, sampling_type) {

  inf_func <- GetInfluenceFcn(valid_dat, dat_sim, mod, sampling_type)

  dat_IF <- bind_cols(dat_sim, inf_func)

  phase2_calibration <- CalibrateDesign(dat_IF)

  raking_fit <- FitCalibrationModel(phase2_calibration)

}



GetInfluenceFcn <- function(valid_dat, dat_sim, mod, sampling_type) {

  if (mod == "RC") {
    inf_func_mod <- FitRCModel(valid_dat, dat_sim, sampling_type)  # fit RC model
  } else if (mod == "naive") {
    inf_func_mod <- FitCoxModel(dat_sim$time_star, dat_sim$delta_star, dat_sim$x_star, dat_sim$z)
  } else {
    stop("Specify 'RC' or 'naive' for calibration")
  }

  # extract influence functions from model
  inf_func_df <- data.frame(resid(inf_func_mod, "dfbeta"))
  colnames(inf_func_df) <- paste("if", 1:2, sep = "")

  return(inf_func_df)
}



CalibrateDesign <- function(IF_dat) {

  IF_design <- twophase(id = list(~id, ~id), subset = ~randomized, data = IF_dat)

  IF_raking <- calibrate(IF_design, phase = 2, formula = ~if1+if2, calfun = "raking")

}



FitCalibrationModel <- function(raking_dat) {

  raking_mod <- svycoxph(Surv(time, delta) ~ x + z, design = raking_dat)

  return(list(raking_mod$coef[1], sqrt(raking_mod$var[1,1]),
              raking_mod$coef[2], sqrt(raking_mod$var[2,2])))

}

