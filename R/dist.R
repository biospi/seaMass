#' summarise using median/mad
#'
#' @export
summarise_normal_robust <- function(value) {
  return(list(m = median(value), s = mad(value), df = Inf))
}


#' fit normal distribution with fitdistrplus
#'
#' @export
summarise_normal_robust_mcmc <- function(chainID, mcmcID, value, plots = FALSE) {
  return(c(list(rhat = rhat(chainID, mcmcID, value, T)), summarise_normal_robust(value)))
}


#' rhat
#'
#' @export
rhat <- function(chainID, mcmcID, value, transform = FALSE) {
  if (length(unique(chainID)) > 1) {
    DT <- data.table::dcast(data.table(mcmcID = mcmcID, chainID = chainID, value = value), mcmcID ~ chainID, value.var = "value")
    DT[, mcmcID := NULL]
    rhat <- coda::gelman.diag(coda::as.mcmc.list(lapply(DT, coda::as.mcmc)), transform, autoburnin = F)$psrf[1]
    return(ifelse(is.nan(rhat), NA_real_, rhat))
  } else {
    return(NA_real_)
  }
}


#' fit normal distribution with fitdistrplus
#'
#' @export
dist_normal <- function(value, plots = FALSE, ...) {
  est <- summarise_normal_robust(value)

  if (est$s > 1e-16) {
    tryCatch({
      ft <- fitdistrplus::fitdist(value, "norm", start = list(mean = est$m, sd = est$s), ...)
      if (plots) plot(ft)
      est <- list(m = ft$estimate[["mean"]], s = ft$estimate[["sd"]], df = Inf)
    }, error = function(e) {
      warning("'dist_normal' fitting failed, falling back to robust approximation.")
    })
  }

  return(est)
}


#' fit normal distribution with fitdistrplus
#'
#' @export
dist_normal_mcmc <- function(chainID, mcmcID, value, plots = FALSE, ...) {
  return(c(list(rhat = rhat(chainID, mcmcID, value, T)), dist_normal(value, plots, method = "mge", gof = "CvM", ...)))
}


#' fit location-scale t distribution with fitdistrplus
#'
#' @export
dist_lst <- function(value, plots = FALSE, ...) {
  est <- summarise_normal_robust(value)
  est$df <- Inf

  if (est$s > 1e-16) {
   d_bayesprot_lst <<- function(x, log_df, mu, log_sigma, log = FALSE) {
     #print(paste(exp(log_df), mu, exp(log_sigma)))
      if (length(x) > 0) {
        return(extraDistr::dlst(x, exp(log_df), mu, exp(log_sigma), log))
      } else {
        return(vector("numeric", 0))
      }
    }

    p_bayesprot_lst <<- function(q, log_df, mu, log_sigma, lower.tail = TRUE, log.p = FALSE) {
      #print(paste(exp(log_df), mu, exp(log_sigma)))
      if (length(q) > 0) {
        return(extraDistr::plst(q, exp(log_df), mu, exp(log_sigma), lower.tail, log.p))
      } else {
        return(vector("numeric", 0))
      }
    }

    q_bayesprot_lst <<- function(p, log_df, mu, log_sigma, lower.tail = TRUE, log.p = FALSE) {
      #print(paste(exp(log_df), mu, exp(log_sigma)))
      if (length(p) > 0) {
        return(extraDistr::qlst(p, exp(log_df), mu, exp(log_sigma), lower.tail, log.p))
      } else {
        return(vector("numeric", 0))
      }
    }

    tryCatch({
      ft <- fitdistrplus::fitdist(value, "_bayesprot_lst", start = list(log_df = 0, mu = est$m, log_sigma = log(est$s)), ...)
      if (plots == T) plot(ft)
      est <- list(m = ft$estimate[["mu"]], s = exp(ft$estimate[["log_sigma"]]), df = exp(ft$estimate[["log_df"]]))
    }, error = function(e) {
      warning("'dist_lst' fitting failed, falling back to robust normal approximation.")
    })
  }

  return(est)
}


#' fit location-scale t distribution
#'
#' @export
dist_lst_mcmc <- function(chainID, mcmcID, value, plots = FALSE, ...) {
  return(c(list(rhat = rhat(chainID, mcmcID, value, T)), dist_lst(value, plots, method = "mge", gof = "CvM", ...)))
}


#' fit scaled inverse chi squared distribution with fitdistrplus
#'
#' @export
dist_invchisq <- function(value, plots = FALSE, ...) {
  est <- summarise_normal_robust(log(value))
  est <- list(v = exp(est$m), df = 2.0 * est$s^-2)

  if (!is.infinite(est$df)) {
    d_bayesprot_invchisq <<- function(x, log_df, log_v, log = F) {
      if (length(x) > 0) {
        return(extraDistr::dinvchisq(x, exp(log_df), exp(log_v), log))
      } else {
        return(vector("numeric", 0))
      }
    }

    p_bayesprot_invchisq <<- function(q, log_df, log_v, lower.tail = T, log.p = F) {
      if (length(q) > 0) {
        return(extraDistr::pinvchisq(q, exp(log_df), exp(log_v), lower.tail, log.p))
      } else {
        return(vector("numeric", 0))
      }
    }

    q_bayesprot_invchisq <<- function(p, log_df, log_v, lower.tail = T, log.p = F) {
      if (length(p) > 0) {
        return(extraDistr::qinvchisq(p, exp(log_df), exp(log_v), lower.tail, log.p))
      } else {
        return(vector("numeric", 0))
      }
    }

    tryCatch({
      ft <- fitdistrplus::fitdist(value, "_bayesprot_invchisq", start = list(log_df = log(est$df), log_v = log(est$v)), ...)
      if (plots == T) plot(ft)
      est <- list(v = exp(ft$estimate[["log_v"]]), df = exp(ft$estimate[["log_df"]]))
    }, error = function(e) {
      warning("'dist_invchisq' fitting failed, falling back to robust approximation.")
    })
  }

  return(est)
}


#' fit scaled inverse chi squared distribution
#'
#' Clarke et al 2012 ("A fast robust method for fitting gamma distributions") are right - CvM is about the only robust way to fit gamma distributions!!
#'
#' @export
dist_invchisq_mcmc <- function(chainID, mcmcID, value, plots = FALSE, ...) {
  return(c(list(rhat = rhat(chainID, mcmcID, value, T)), dist_invchisq(value, plots, method = "mge", gof = "CvM", ...)))
}


#' fit scaled F distribution with per-datapoint fixed df1 with fitdistrplus
#'
#' @export
dist_sf_with_fixed_df1_fitdistrplus <- function(value, df1, plots = FALSE, ...) {
  # first fit a _sample_ of value to an inverse chi squared distribution
  est0 <- dist_invchisq(sapply(1:length(value), function(i) extraDistr::rinvchisq(1, df1[i], value[i])), plots = plots)

  # this then seeds the fit to a scaled F distribution
  d_bayesprot_scaled_f <<- function(x, log_df2, log_scale, log = FALSE) {
    #print(paste("d", exp(log_df2), exp(log_scale)))
    if (length(x) > 0) {
      return(sapply(1:length(x), function(i) extraDistr::dbetapr(2^(x[i]), df1[i]/2, exp(log_df2)/2, exp(log_df2)/df1[i] * exp(log_scale), log)))
    } else {
      return(vector("numeric", 0))
    }
  }

  p_bayesprot_scaled_f <<- function(q, log_df2, log_scale, lower.tail = TRUE, log.p = FALSE) {
    #print(paste("p", exp(log_df2), exp(log_scale)))
    if (length(q) > 0) {
      return(sapply(1:length(q), function(i) extraDistr::pbetapr(q[i], df1[i]/2, exp(log_df2)/2, exp(log_df2)/df1[i] * exp(log_scale), lower.tail, log.p)))
    } else {
      return(vector("numeric", 0))
    }
  }

  q_bayesprot_scaled_f <<- function(p, log_df2, log_scale, lower.tail = TRUE, log.p = FALSE) {
    #print(paste("q", exp(log_df2), exp(log_scale)))
    if (length(p) > 0) {
      return(sapply(1:length(p), function(i) extraDistr::qbetapr(p[i], df1[i]/2, exp(log_df2)/2, exp(log_df2)/df1[i] * exp(log_scale), lower.tail, log.p)))
    } else {
      return(vector("numeric", 0))
    }
  }

  tryCatch({
    ft <- fitdistrplus::fitdist(log2(value), "_bayesprot_scaled_f", start = list(log_df2 = log(est0$df), log_scale = log(est0$v)), ...)
    if (plots == T) plot(ft)
    est <- list(v = exp(ft$estimate[["log_scale"]]), df = exp(ft$estimate[["log_df2"]]))
  }, error = function(e) {
    warning("'dist_sf_with_fixed_df1_fitdistrplus' fitting failed, falling back to 'dist_invchisq' fit to random sample.")
    est <- est0
  })

  return(list(v0 = est0$v, df0 = est0$df, v = est$v, df = est$df))
}

#' Our squeezeVar function - works for small DF unlike Limma squeezeVar, but slower
#'
#' @export
squeeze_var <- function(v, df, use.deconvolution = TRUE, plots = FALSE, ...) {
  if (use.deconvolution) {
    return(dist_sf_with_fixed_df1_fitdistrplus(v, df, plots = plots, ...))
  } else {
    est <- dist_invchisq(sapply(1:length(v), function(i) extraDistr::rinvchisq(1, df[i], v[i])), plots = plots, ...)
    return(list(v0 = est$v, df0 = est$df, v = est$v, df = est$df))
  }
}


#' Limma squeezeVar function
#'
#' @export
squeeze_var_limma <- function(v, df) {
  ft <- limma::fitFDist(v, df)
  return(list(v = ft$scale, df = ft$df2))
}


# dist_invchisq_fitdistrplus <- function(value, plot = F, ...) {
#   est <- summarise_normal_robust(value)
#   s <- est$s
#   est <- list(v = est$m, df = Inf)
#
#   if (s > 0.0000001) {
#     d_bayesprot_invchisq <<- function(value, log_nu, log_tau, log = FALSE) {
#       if (length(value) > 0) {
#         return(extraDistr::dinvchisq(value, exp(log_nu), exp(log_tau), log))
#       } else {
#         return(vector("numeric", 0))
#       }
#     }
#
#     p_bayesprot_invchisq <<- function(q, log_nu, log_tau, lower.tail = TRUE, log.p = FALSE) {
#       if (length(q) > 0) {
#         return(extraDistr::pinvchisq(q, exp(log_nu), exp(log_tau), lower.tail, log.p))
#       } else {
#         return(vector("numeric", 0))
#       }
#     }
#
#     q_bayesprot_invchisq <<- function(p, log_nu, log_tau, lower.tail = TRUE, log.p = FALSE) {
#       if (length(p) > 0) {
#         return(extraDistr::qinvchisq(p, exp(log_nu), exp(log_tau), lower.tail, log.p))
#       } else {
#         return(vector("numeric", 0))
#       }
#     }
#
#     ft <- fitdistrplus::fitdist(value, "_bayesprot_invchisq", start = list(log_nu = 0, log_tau = log(est$v)), ...)
#     if (plot) plot(ft)
#
#     est <- list(v = exp(ft$estimate[["log_tau"]]), df = exp(ft$estimate[["log_nu"]]))
#   }
#
#   return(est)
# }
#
#dist_invchisq <- function(value, plot = F) {
#  return(dist_invchisq_fitdistrplus(value, plot, method = "mge", gof = "KS", optim.method = "BFGS"))
#}
#
# dist_sf_with_fixed_df1_fitdistrplus <- function(value, df1, plot = F, ...) {
#   # first fit a _sample_ of value to an inverse chi squared distribution
#   est <- dist_invchisq_fitdistrplus(sapply(1:length(value), function(i) extraDistr::rinvchisq(1, df1[i], value[i])), ...)
#
#   # this then seeds the fit to a scaled F distribution
#   d_bayesprot_scaled_f <<- function(value, log_df2, log_scale, log = FALSE) {
#     if (length(value) > 0) {
#       return(sapply(1:length(value), function(i) extraDistr::dbetapr(value[i], df1[i]/2, exp(log_df2)/2, exp(log_df2)/df1[i] * exp(log_scale), log)))
#     } else {
#       return(vector("numeric", 0))
#     }
#   }
#
#   p_bayesprot_scaled_f <<- function(q, log_df2, log_scale, lower.tail = TRUE, log.p = FALSE) {
#     if (length(q) > 0) {
#       return(sapply(1:length(q), function(i) extraDistr::pbetapr(q[i], df1[i]/2, exp(log_df2)/2, exp(log_df2)/df1[i] * exp(log_scale), lower.tail, log.p)))
#     } else {
#       return(vector("numeric", 0))
#     }
#   }
#
#   q_bayesprot_scaled_f <<- function(p, log_df2, log_scale, lower.tail = TRUE, log.p = FALSE) {
#     if (length(p) > 0) {
#       return(sapply(1:length(p), function(i) extraDistr::qbetapr(p[i], df1[i]/2, exp(log_df2)/2, exp(log_df2)/df1[i] * exp(log_scale), lower.tail, log.p)))
#     } else {
#       return(vector("numeric", 0))
#     }
#   }
#
#   ft <- fitdistrplus::fitdist(value, "_bayesprot_scaled_f", start = list(log_df2 = log(est$df), log_scale = log(est$v)), ...)
#   if (plot) plot(ft)
#
#   return(list(
#     v0 = est$v, df0 = est$df,
#     v = exp(ft$estimate[["log_scale"]]), df = exp(ft$estimate[["log_df2"]])
#   ))
# }
#
# squeeze_var_limma <- function(v, df) {
#   ft <- limma::fitFDist(v, df)
#   return(list(v = ft$scale, df = ft$df2))
# }

# dist_invchisq_fitdistrplus <- function(value, plots = FALSE, ...) {
#   est <- summarise_normal_robust(value)
#   s <- est$s
#   est <- list(v = est$m, df = Inf)
#
#   if (s > 0.0000001) {
#     ft <- fitdistrplus::fitdist(1/value, "gamma", ...)
#     if (plots == T) plot(ft)
#
#     est <- list(v = ft$estimate[["rate"]] / ft$estimate[["shape"]], df = 2.0 * ft$estimate[["shape"]])
#   }
#
#   return(est)
# }


#' dist_invchisq5 <- function(value, plots = FALSE, ...) {
#'   # work on the log values
#'   log_value <- log(value)
#'
#'   # approximation for seed
#'   est <- summarise_normal_robust(log_value)
#'   est <- list(v = exp(est$m), df = 2*est$s^-2)
#'
#'
#'   return(est)
#' }
#'
#'
#' dist_invchisq6 <- function(value, ...) {
#'   # approximation for seed
#'   est <- summarise_normal_robust(log(value))
#'   est <- list(v = exp(est$m), df = 2*est$s^-2)
#'   #shape = 2.296046
#'   #rate = 38.91579
#'
#'   if (!is.infinite(est$df)) {
#'     ft <- robust::gammaRob(1/value)
#'     est <- list(v = ft$estimate[["scale"]] / ft$estimate[["shape"]], df = 2 * ft$estimate[["shape"]])
#'     #}, error = function(e) {
#'     #  est$df = 0
#'     #})
#'   }
#'
#'   return(est)
#' }
#'
#'
#'
#'
#' dist_invchisq4 <- function(value, plots = FALSE, ...) {
#'   # approximation for seed
#'   est <- summarise_normal_robust(log(value))
#'   est <- list(v = exp(est$m), df = 2*est$s^-2)
#'   #shape = 2.296046
#'   #rate = 38.91579
#'
#'   if (!is.infinite(est$df)) {
#'     d_bayesprot_expgamma <<- function(value, shape, log_scale, log = FALSE) {
#'       return(dgamma(value, shape, scale = exp(log_scale)) * value)
#'     }
#'
#'     #tryCatch({
#'     ft <- fitdistrplus::fitdist(1/value, "_bayesprot_expgamma", start = list(shape = 0.5*est$df, log_scale = log(0.5*est$df*est$v)), lower = c(1e-16, -Inf), control=list(trace=1, REPORT=1))
#'     if (plots == T) plot(ft)
#'     est <- list(v = exp(ft$estimate[["log_scale"]]) / ft$estimate[["shape"]], df = 2 * ft$estimate[["shape"]])
#'     #}, error = function(e) {
#'     #  est$df = 0
#'     #})
#'   }
#'
#'   return(est)
#' }
#'
#'
#'
#' dist_invchisq5 <- function(value, plots = FALSE, ...) {
#'   # approximation for seeding
#'   est <- summarise_normal_robust(log(value))
#'   est <- list(v = exp(est$m), df = 2*est$s^-2)
#'
#'   if (!is.infinite(est$df)) {
#'     d_bayesprot_expsichisq <<- function(x, log_df, log_s, log = F) {
#'       s <- exp(log_s)
#'       df <- exp(log_df)
#'       df2 <- df/2
#'       y <- df2*log(df2) - lgamma(df2) + df*log(s) - x*df2 - s^2*df/(2*exp(x))
#'       if (log) {
#'         return(y)
#'       } else {
#'         return(exp(y))
#'       }
#'     }
#'
#'     #tryCatch({
#'     ft <- fitdistrplus::fitdist(log(value), "_bayesprot_expsichisq", start = list(log_df = log(est$df), log_s = log(sqrt(est$v))), method = "mle", optim.method = "BFGS", control=list(trace=6, REPORT=1))
#'     if (plots == T) plot(ft)
#'     est <- list(v = exp(ft$estimate[["log_s"]])^2, df = exp(ft$estimate[["log_df"]]))
#'     #}, error = function(e) {
#'     #  est$df = 0
#'     #})
#'   }
#'
#'   return(est)
#' }
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' dist_invchisq5 <- function(value, plots = FALSE, ...) {
#'   # work on the log values
#'   log_value <- log(value)
#'
#'   # approximation for seed
#'   est <- summarise_normal_robust(log_value)
#'   est <- list(v = exp(est$m), df = 2*est$s^-2)
#'   #$v
#'   #[1] 0.01166819
#'   #
#'   #$df
#'   #[1] 5.533201
#'
#'   if (!is.infinite(est$df)) {
#'     d_bayesprot_invgamma <<- function(value, shape, scale, log = FALSE) {
#'       logu <- log(scale) - log(x);
#'       return(shape * logu - exp(logu) - exp(value) - lgammafn(shape))
#'     }
#'
#'     #tryCatch({
#'     ft <- fitdistrplus::fitdist(log_value, "_bayesprot_invgamma", start = list(shape = 0.5 * est$v * est$df, scale = 0.5 * est$df), lower = c(1e-16, 1e-16), control=list(trace=1, REPORT=1))
#'     if (plots == T) plot(ft)
#'     est <- list(v = ft$estimate[["shape"]] / ft$estimate[["rate"]], df = 2 * ft$estimate[["shape"]])
#'     #}, error = function(e) {
#'     #  est$df = 0
#'     #})
#'   }
#'
#'   return(est)
#' }
#'
#'
#'
#'
#'
#' #' fit scaled inverse chi squared distribution with fitdistrplus
#' #'
#' #' @export
#' dist_invchisq <- function(value, plots = FALSE, ...) {
#'   est <- summarise_normal_robust(value)
#'   s <- est$s
#'   est <- list(v = est$m, df = Inf)
#'
#'   if (s > 1e-16) {
#'     tryCatch({
#'       ft <- fitdistrplus::fitdist(1/(value*1e8), "gamma", weights = ceiling(1/(value*1e8)), ...)
#'       if (plots == T) plot(ft)
#'       est <- list(v = ft$estimate[["rate"]] / (1e8 * ft$estimate[["shape"]]), df = 2 * ft$estimate[["shape"]])
#'     }, error = function(e) {
#'       est$df = 0
#'     })
#'   }
#'
#'   return(est)
#' }
#'
#'
#' dist_invchisq2 <- function(value, plots = FALSE, ...) {
#'   est <- summarise_normal_robust(log(value))
#'   s <- est$s
#'   est <- list(v = exp(est$m), df = Inf)
#'
#'   if (s > 1e-16) {
#'     #tryCatch({
#'     ft <- fitdistrplus::fitdist(1/value, "gammaAlt", start = list(mean = est$v, cv = s), lower = c(1e-16, 1e-16), ...)
#'     if (plots == T) plot(ft)
#'     est <- list(v = 1 / ft$estimate[["mean"]], df = ft$estimate[["cv"]])
#'     #}, error = function(e) {
#'     #  est$df = 0
#'     #})
#'   }
#'
#'   return(est)
#' }
#'
#' dist_invchisq <- function(value, plots = FALSE, ...) {
# est <- summarise_normal_robust(value)
# s <- est$s
# est <- list(v = est$m, df = Inf)
#
# if (s > 1e-16) {
#   d_bayesprot_invchisq <<- function(value, nu, tau, log = FALSE) {
#     if (length(value) > 0) {
#       return(extraDistr::dinvchisq(value, nu, tau, log))
#     } else {
#       return(vector("numeric", 0))
#     }
#   }
#
#   p_bayesprot_invchisq <<- function(q, nu, tau, lower.tail = TRUE, log.p = FALSE) {
#     if (length(q) > 0) {
#       return(extraDistr::pinvchisq(q, nu, tau, lower.tail, log.p))
#     } else {
#       return(vector("numeric", 0))
#     }
#   }
#
#   q_bayesprot_invchisq <<- function(p, nu, tau, lower.tail = TRUE, log.p = FALSE) {
#     if (length(p) > 0) {
#       return(extraDistr::qinvchisq(p, nu, tau, lower.tail, log.p))
#     } else {
#       return(vector("numeric", 0))
#     }
#   }
#
#   tryCatch({
#     ft <- fitdistrplus::fitdist(value, "_bayesprot_invchisq", start = list(nu = 1, tau = est$v), lower = c(1e-16, 1e-16), ...)
#     if (plots == T) plot(ft)
#     est <- list(v = ft$estimate[["tau"]], df = ft$estimate[["nu"]])
#   }, error = function(e) {
#     est$df = 0
#   })
# }
#
# return(est)
# }

