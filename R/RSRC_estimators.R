#' Calculates regression calibration estimates
#'
#' These functions implement the risk set
#' regression calibration estimators.
#' Moment estimators are fit using least
#' squares and are used to obtain our best
#' prediction of the true covariate and
#' censored event time. This is repeated
#' at deciles of the failure times.
#'
#' @param valid_dat Validation subset
#'
#' @param dat_sim Full dataset
#'
#' @param sampling_type String indicating either simple random
#' sampling or case-cohort sampling
#'
#' @param beta_x_start Initial guess for beta_x in optimization
#'
#' @param beta_x_start Initial guess for beta_z in optimization
#'
#' @return List of RSRC beta_x and beta_z estimates
#'
#' @rdname RSRC_estimators
#' @export
FitRSRCModel <- function(valid_dat, dat_sim, sampling_type, beta_x_start, beta_z_start) {

  if (sampling_type == "cc") {
    valid.dat <- valid.dat %>%
      dplyr::mutate(wts = 1/twophase(id = list(~id, ~id),
                                     subset = ~randomized,
                                     strata = list(NULL, ~delta.star),
                                     data = dat_sim)$prob)
  }

  # Recalibrate error-prone covariate and failure time
  RSRC_mats <- CalcRSRCXYhat(valid_dat, dat_sim, sampling_type)

  # Get event times for all subjects
  which_fail <- which(dat_sim$delta_star == 1)
  fail_times <- dat_sim$time_hatXY[which(dat_sim$delta_star == 1)]

  # Get U^{hat} at each event time
  timehat_lst <- lapply(1:length(fail_times),
                        function(i) RSRC_mats[[2]][,(sum(fail_times[i] >= RSRC_mats[[3]]) + 1)])
  timehat_mat <- matrix(unlist(timehat_lst), ncol = length(fail_times))

  # Get risk set at each event time
  risk_set_lst <- lapply(1:length(fail_times),
                         function(i) ifelse(timehat_mat[,i] >= timehat_mat[which_fail[i], i] & !is.na(timehat_mat[,i]), 1, 0))
  risk_set_mat <- matrix(unlist(risk_set_lst), ncol = length(fail_times))

  # Get X^{hat} at each event time
  xhat_lst <- lapply(1:length(fail_times),
                     function(i) RSRC_mats[[1]][, (sum(fail_times[i] >= RSRC_mats[[3]]) + 1)])
  xhat_mat <- matrix(unlist(xhat_lst), ncol = length(fail_times))

  z_mat <- matrix(dat_sim$z, nrow = length(dat_sim$z), ncol = length(fail_times))

  xhat_risk_set <- risk_set_mat * xhat_mat
  z_risk_set <- risk_set_mat * z_mat

  xhat_subj <- unlist(lapply(1:length(fail_times),
                             function(i) xhat_mat[which_fail[i], i]))
  z_subj <- unlist(lapply(1:length(fail_times),
                          function(i) z_mat[which_fail[i], i]))

  # Solve the Cox score equation to find the roots
  score_roots <- tryCatch(multiroot(f = RSRCscore, start = c(beta_x_start, beta_z_start), maxiter = 1000, ctol = 0.001,
                                    xhat_mat = xhat_mat, xhat_subj = xhat_subj, z_mat = z_mat, z_subj = z_subj,
                                    xhat_risk_set = xhat_risk_set, z_risk_set = z_risk_set, risk_set_mat = risk_set_mat),
                          error = function(e) {
                            return(list(NA, NA))
                          }
  )

  if (is.na(score_roots[[1]][1])) {
    return(list(NA, NA))
  } else {
    return(list(score_roots$root[1], score_roots$root[2]))
  }

}



CalcRSRCXYhat <- function(valid_dat, dat_sim, sampling_type) {

  num_fail <- sum(dat_sim$delta_star)

  # sort distinct failure times
  ordered_obs_time <- dat_sim %>%
    dplyr::filter(delta_star == 1) %>%
    dplyr::select(time_hatXY) %>%
    dplyr::arrange(time_hatXY)

  # get step sizes for RC by deciles (can change)
  num_recalibrate <- 10
  step <- floor(num_fail/num_recalibrate)

  # get equally spaced failure times to recalibrate at
  RSRC_times <- ordered_obs_time[seq(from = step, to = (num_recalibrate - 1) * step, by = step),]

  x_hat_mat <- time_hat_mat <- matrix(NA, nrow = dim(dat_sim)[1], ncol = length(RSRC_times) + 1)

  # first column is standard RC using all data
  x_hat_mat[, 1] <- CalcExpX(valid_dat, dat_sim, sampling_type)
  time_hat_mat[, 1] <- dat_sim$time_hatXY - CalcExpw_RSRC(valid_dat, dat_sim, sampling_type)


  for(RSRC_ind in 1:length(RSRC_times)) {

    # get risk set for specific recalibration time
    curr_risk_set <- dat_sim$time_hatXY >= RSRC_times[RSRC_ind]
    valid_risk_set <- curr_risk_set[dat_sim$randomized == TRUE]

    # if there are at least 20 subjects in the risk set, then recalibrate
    if (sum(curr_risk_set) >= 20) {

      x_hat_mat[curr_risk_set, RSRC_ind + 1] <- CalcExpX(valid_dat[valid_risk_set,],
                                                        dat_sim[curr_risk_set,],
                                                        sampling_type)

      time_hat_mat[curr_risk_set, RSRC_ind + 1] <- dat_sim[curr_risk_set,"time_hatXY"] -
        CalcExpw_RSRC(valid_dat[valid_risk_set,], dat_sim[curr_risk_set,], sampling_type)

    } else {
      x_hat_mat[curr_risk_set, RSRC_ind + 1] <- x_hat_mat[curr_risk_set, RSRC_ind]
      time_hat_mat[curr_risk_set, RSRC_ind + 1] <- time_hat_mat[curr_risk_set, RSRC_ind]
    }
  }

  return(list(x_hat_mat, time_hat_mat, RSRC_times))

}



RSRCscore <- function(betas, xhat_mat, xhat_subj, z_mat, z_subj,
                     xhat_risk_set, z_risk_set, risk_set_mat) {

  exp_beta_x <- exp(betas[1] * xhat_mat + betas[2] * z_mat)

  score_x <- sum(xhat_subj - unlist(lapply(1:length(xhat_subj), function(i) sum(xhat_risk_set[,i] * exp_beta_x[,i], na.rm = TRUE) /
                                            sum(risk_set_mat[,i] * exp_beta_x[,i], na.rm = TRUE))))

  score_z <- sum(z_subj - unlist(lapply(1:length(z_subj), function(i) sum(z_risk_set[,i] * exp_beta_x[,i], na.rm = TRUE) /
                                         sum(risk_set_mat[,i] * exp_beta_x[,i], na.rm = TRUE))))

  c(score_x, score_z)

}




CalcExpw_RSRC <- function(valid_dat, dat_sim, sampling_type) {

  valid_dat <- valid_dat %>%
    dplyr::mutate(new_total_y_err = time_hatXY - time)

  if (sampling_type == "cc") {
    w_mod <- lm(new_total_y_err ~ x_star + z, data = valid_dat, weights = wts)
  } else if (sampling_type == "srs") {
    w_mod <- lm(new_total_y_err ~ x_star + z, data = valid_dat)
  }

  w_predict <- predict(w_mod, newdata = dat_sim)

}




