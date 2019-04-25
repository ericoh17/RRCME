#' Calculates regression calibration estimates
#'
#' These functions implement the regression
#' calibration estimators. Moment estimators
#' are fit using least squares and are used
#' to obtain our best prediction of the
#' true covariate and censored event time.
#'
#' @param cox_time Censored time-to-event
#'
#' @param cox_delta Failure indicator
#'
#' @param cox_x Univariate covariate measured
#' with error
#'
#' @param cox_z Univariate covariate measured
#' without error
#'
#' @param valid_dat Validation subset
#'
#' @param dat_sim Full dataset
#'
#' @param sampling_type String indicating either simple random
#' sampling or case-cohort sampling
#'
#' @return Eiher a list of Cox model coefficients
#' and standard errors of a Cox model fit object
#'
#' @importFrom dplyr mutate
#'
#' @importFrom survival coxph
#'
#' @importFrom survey twophase
#'
#' @rdname RC_estimators
#' @export
FitRCModel <- function(valid_dat, dat_sim, sampling_type, return_coef) {

  if (sampling_type == "cc") {
    valid_dat <- valid_dat %>%
      mutate(wts = 1/twophase(id = list(~id, ~id),
                              subset = ~randomized,
                              strata = list(NULL, ~delta_star),
                              data = dat_sim)$prob)
  }

  x_hat <- CalcExpX(valid_dat, dat_sim, sampling_type)

  w_hat <- CalcExpw(valid_dat, dat_sim, sampling_type)
  time_hat <- dat_sim$time_star - w_hat

  rc_mod <- FitCoxModel(time_hat, dat_sim$delta_star, x_hat, dat_sim$z,
                        return_coef)

  if (return_coef == TRUE) {
    return(list(rc_mod[[1]], rc_mod[[2]]))
  } else {
    return(rc_mod)
  }

}


FitCoxModel <- function(cox_time, cox_delta, cox_x, cox_z, return_coef) {

  cox_mod <- coxph(Surv(cox_time, cox_delta) ~ cox_x + cox_z)

  if (return_coef == TRUE) {
    return(list(cox_mod$coef[1], cox_mod$coef[2]))
  } else {
    return(cox_mod)
  }

}


#' @rdname RC_estimators
#' @export
CalcExpX <- function(valid_dat, dat_sim, sampling_type) {

  if (sampling_type == "cc") {
    x_mod <- lm(x ~ x_star + z, data = valid_dat, weights = wts)
  } else if (sampling_type == "srs") {
    x_mod <- lm(x ~ x_star + z, data = valid_dat)
  }

  x_predict <- predict(x_mod, newdata = dat_sim)

}


#' @rdname RC_estimators
#' @export
CalcExpw <- function(valid_dat, dat_sim, sampling_type) {

  if (sampling_type == "cc") {
    w_mod <- lm(total_y_err ~ x_star + z, data = valid_dat, weights = wts)
  } else if (sampling_type == "srs") {
    w_mod <- lm(total_y_err ~ x_star + z, data = valid_dat)
  }

  w_predict <- predict(w_mod, newdata = dat_sim)

}










