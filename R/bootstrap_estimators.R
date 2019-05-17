#' Calculates bootstrap estimates for all methods
#'
#' These unctions that implement the bootstrap for all considered methods
#' Note: each bootstrap is stratified on inclusion into the
#' validation subset (and on case status for case-cohort sampling)
#'
#' \code{RunRCBootstrap} implements the bootstrap for the RC methods
#'
#' \code{RunRSRCBoostrap} implements the boostrap for the RSRC methods
#'
#' \code{RunRakingBootstrap} implements the bootstrap for the raking
#' methods
#'
#' \code{RunCompleteCaseBootstrap} implements the bootstrap for the
#' complete case analysis
#'
#' \code{RunNaiveBoostrap} implements the bootstrap for the naive
#' analysis
#'
#' @param all_dat Full, original dataset
#'
#' @param inds Vector of indices defining bootstrap sample
#'
#' @param dat_valid Original validation subset
#'
#' @param beta_x_start Initial guess for beta_x in optimization
#'
#' @param beta_z_start Initial guess for beta_z in optimization
#'
#' @param mod_rake String indicating which model to get influence
#' functions from
#'
#' @param sampling_type String indicating either simple random
#' sampling or case-cohort sampling
#'
#' @return A vector of bootstrap beta_x and beta_z estimates
#'
#' @importFrom dplyr filter
#'
#' @importFrom survey twophase svycoxph
#'
#' @importFrom magrittr %>%
#'
#' @rdname bootstrap_estimators
#' @export
RunRCBootstrap <- function(all_dat, inds,
                           sampling_type) {

  full_dat <- all_dat[inds,]
  valid_dat <- full_dat %>% filter(randomized == TRUE)

  rc_boot_mod <- FitRCModel(valid_dat, full_dat,
                            sampling_type,
                            return_coef = TRUE)

  return(c(rc_boot_mod[[1]], rc_boot_mod[[2]]))

}


#' @rdname bootstrap_estimators
#' @export
RunRSRCBootstrap <- function(all_dat, inds,
                             sampling_type,
                             beta_x_start,
                             beta_z_start) {

  full_dat <- all_dat[inds,]
  valid_dat <- full_dat %>% filter(randomized == TRUE)

  RSRC_boot_mod <- FitRSRCModel(valid_dat, full_dat,
                                sampling_type,
                                beta_x_start,
                                beta_z_start)

  return(c(RSRC_boot_mod[[1]], RSRC_boot_mod[[2]]))

}


#' @rdname bootstrap_estimators
#' @export
RunRakingBootstrap <- function(all_dat, inds,
                               mod_rake,
                               sampling_type) {

  full_dat <- all_dat[inds,]
  valid_dat <- full_dat %>% filter(randomized == TRUE)

  raking_boot_mod <- FitRakingModel(valid_dat, full_dat,
                                    mod_rake,
                                    sampling_type)

  return(c(raking_boot_mod[[1]], raking_boot_mod[[2]]))

}


#' @rdname bootstrap_estimators
#' @export
RunCompleteCaseBootstrap <- function(all_dat, inds,
                                     sampling_type) {

  full_dat <- all_dat[inds,]
  valid_dat <- full_dat %>% filter(randomized == TRUE)

  if (sampling_type == "srs") {

    complete_case_boot_mod <- FitCoxModel(valid_dat$time, valid_dat$delta,
                                          valid_dat$x, valid_dat$z,
                                          return_coef = TRUE)

    return(c(complete_case_boot_mod[[1]], complete_case_boot_mod[[2]]))

  } else if (sampling_type == "cc") {

    complete_case_boot_design <- twophase(id = list(~id, ~id), subset = ~randomized,
                                          strata = list(NULL, ~delta_star),
                                          data = full_dat)
    complete_case_boot_mod <- svycoxph(Surv(time, delta) ~ x + z,
                                       design = complete_case_boot_design)

    return(c(coef(complete_case_boot_mod)[1], coef(complete_case_boot_mod)[2]))

  }

}


#' @rdname bootstrap_estimators
#' @export
RunNaiveBootstrap <- function(all_dat, inds) {

  full_dat <- all_dat[inds,]

  naive_boot_mod <- FitCoxModel(full_dat$time_star, full_dat$delta_star,
                                full_dat$x_star, full_dat$z,
                                return_coef = TRUE)

  return(c(naive_boot_mod[[1]], naive_boot_mod[[2]]))

}



