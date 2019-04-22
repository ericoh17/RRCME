#' Case-cohort versions of raking estimators
#'
#' These functions implement raking estimators
#' when the validation subset was selected using
#' a case-cohort sampling scheme.
#'
#' @param valid_dat Validation subset dataset
#'
#' @param dat_sim Full dataset
#'
#' @param mod String indicating which model to get influence
#' functions from
#'
#' @return Raking model fit object
#' @rdname caseCohort_functions
#' @export
FitRakingModel_CC <- function(valid_dat, dat_sim, mod) {

  inf_func <- GetInfluenceFcn(valid_dat, dat_sim, mod)

  dat_IF <- dplyr::bind_cols(dat_sim, inf_func)

  phase2_calibration <- CalibrateDesign_CC(dat_IF)

  raking_fit <- FitCalibrationModel(phase2_calibration)

}


CalibrateDesign_CC <- function(IF_dat) {

  IF_design <- twophase(id = list(~id, ~id), subset = ~randomized, strata = list(NULL, ~delta_star), data = IF_dat)

  IF_raking <- calibrate(IF_design, phase = 2, formula = ~if1+if2+delta_star, calfun = "raking")

}




