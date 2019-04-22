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
#' @param return_coef Boolean indicating whether
#' to return a list of coefficients
#'
#' @param valid_dat Validation subset
#'
#' @param dat_sim Full dataset
#'
#' @return Eiher a list of Cox model coefficients
#' and standard errors of a Cox model fit object
#'
#' @rdname RC_estimators
#' @export
FitRCModel <- function(valid_dat, dat_sim, return_coef = FALSE) {

  x_hat <- CalcExpX(valid_dat, dat_sim)

  w_hat <- CalcExpw(valid_dat, dat_sim)
  time_hat <- dat_sim$time_star - w_hat

  rc_mod <- FitCoxModel(time_hat, dat_sim$delta_star, x_hat, dat_sim$z)

  if (return_coef == TRUE) {
    return(list(rc_mod$coef[1], summary(rc_mod)$coef[1,3],
                rc_mod$coef[2], summary(rc_mod)$coef[2,3]))
  } else {
    return(rc_mod)
  }

}


FitCoxModel <- function(cox_time, cox_delta, cox_x, cox_z, return_coef = FALSE) {

  cox_mod <- coxph(Surv(cox_time, cox_delta) ~ cox_x + cox_z)

  if (return_coef == TRUE) {
    return(list(cox_mod$coef[1], summary(cox_mod)$coef[1,3],
                cox_mod$coef[2], summary(cox_mod)$coef[2,3]))
  } else {
    return(cox_mod)
  }
}



CalcExpX <- function(valid_dat, dat_sim) {

  x_mod <- lm(x ~ x_star + z, data = valid_dat)
  x_predict <- predict(x_mod, newdata = dat_sim)

}



CalcExpw <- function(valid_dat, dat_sim) {

  w_mod <- lm(total_y_err ~ x_star + z, data = valid_dat)
  w_predict <- predict(w_mod, newdata = dat_sim)

}










