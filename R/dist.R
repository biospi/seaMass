#' summarise using median/mad
#'
#' @export
summarise_normal_robust <- function(x) {
  return(list(m = median(x), s = mad(x)))
}


#' rhat
#'
#' @export
rhat <- function(mcmcID, chainID, value) {
  if (length(unique(chainID)) > 1) {
    DT <- data.table::dcast(data.table(mcmcID = mcmcID, chainID = chainID, value = value), mcmcID ~ chainID, value.var = "value")
    DT[, mcmcID := NULL]
    return(coda::gelman.diag(coda::as.mcmc.list(lapply(DT, coda::as.mcmc)), autoburnin = F)$psrf[1])
  } else {
    return(NA_real_)
  }
}


#' fit normal distribution with fitdistrplus
#'
#' @export
dist_normal_fitdistrplus <- function(x, plots = FALSE, ...) {
  est <- summarise_normal_robust(x)

  if (est$s > 0.0000001) {
    ft <- fitdistrplus::fitdist(x, "norm", start = list(mean = est$m, sd = est$s), ...)
    if (plots) plot(ft)

    est <- list(m = ft$estimate[["mean"]], s = ft$estimate[["sd"]])
  }

  return(est)
}


#' fit normal distribution with fitdistrplus
#'
#' @export
dist_normal <- function(x, plots = FALSE) {
  return(dist_normal_fitdistrplus(x, plots = plots))
}


#' fit location-scale t distribution with fitdistrplus
#'
#' @export
dist_lst_fitdistrplus <- function(x, plots = FALSE, ...) {
  est <- summarise_normal_robust(x)
  est$df <- Inf

  if (est$s > 0.0000001) {
   d_bayesprot_lst <<- function(x, df, mu, sigma, log = FALSE) {
      if (length(x) > 0) {
        return(extraDistr::dlst(x, df, mu, sigma, log))
      } else {
        return(vector("numeric", 0))
      }
    }

    p_bayesprot_lst <<- function(q, df, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
      if (length(q) > 0) {
        return(extraDistr::plst(q, df, mu, sigma, lower.tail, log.p))
      } else {
        return(vector("numeric", 0))
      }
    }

    q_bayesprot_lst <<- function(p, df, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
      if (length(p) > 0) {
        return(extraDistr::qlst(p, df, mu, sigma, lower.tail, log.p))
      } else {
        return(vector("numeric", 0))
      }
    }

    ft <- fitdistrplus::fitdist(x, "_bayesprot_lst", start = list(df = 1, mu = est$m, sigma = est$s), lower = c(0.00001, -Inf, 0.00001), ...)
    if (plots == T) plot(ft)

    est <- list(m = ft$estimate[["mu"]], s = ft$estimate[["sigma"]], df = ft$estimate[["df"]])
  }

  return(est)
}


#' fit location-scale t distribution with fitdistrplus using quantile matching
#'
#' @export
dist_lst <- function(x, plots = FALSE) {
  return(dist_lst_fitdistrplus(x, plots = plots))
}


#' fit scaled inverse chi squared distribution with fitdistrplus
#'
#' @export
dist_invchisq_fitdistrplus <- function(x, plots = FALSE, ...) {
  est <- summarise_normal_robust(x)
  s <- est$s
  est <- list(v = est$m, df = Inf)

  if (s > 0.0000001) {
    d_bayesprot_invchisq <<- function(x, nu, tau, log = FALSE) {
      if (length(x) > 0) {
        return(extraDistr::dinvchisq(x, nu, tau, log))
      } else {
        return(vector("numeric", 0))
      }
    }

    p_bayesprot_invchisq <<- function(q, nu, tau, lower.tail = TRUE, log.p = FALSE) {
      if (length(q) > 0) {
        return(extraDistr::pinvchisq(q, nu, tau, lower.tail, log.p))
      } else {
        return(vector("numeric", 0))
      }
    }

    q_bayesprot_invchisq <<- function(p, nu, tau, lower.tail = TRUE, log.p = FALSE) {
      if (length(p) > 0) {
        return(extraDistr::qinvchisq(p, nu, tau, lower.tail, log.p))
      } else {
        return(vector("numeric", 0))
      }
    }

    ft <- fitdistrplus::fitdist(x, "_bayesprot_invchisq", start = list(nu = 1, tau = est$v), lower = c(0.00001, 0.00001), ...)
    if (plots == T) plot(ft)

    est <- list(v = ft$estimate[["tau"]], df = ft$estimate[["nu"]])
  }

  return(est)
}


#' fit scaled inverse chi squared distribution with fitdistrplus
#'
#' @export
dist_invchisq <- function(x, plots = FALSE) {
  return(dist_invchisq_fitdistrplus(x, plots = plots))
}


#' fit scaled F distribution with per-datapoint fixed df1 with fitdistrplus
#'
#' @export
dist_sf_with_fixed_df1_fitdistrplus <- function(x, df1, plots = FALSE, ...) {
  # first fit a _sample_ of x to an inverse chi squared distribution
  est <- dist_invchisq_fitdistrplus(sapply(1:length(x), function(i) extraDistr::rinvchisq(1, df1[i], x[i])), plots = plots, ...)

  # this then seeds the fit to a scaled F distribution
  d_bayesprot_scaled_f <<- function(x, df2, scale, log = FALSE) {
    #print(paste("d", df2, scale))
    if (length(x) > 0) {
      return(sapply(1:length(x), function(i) extraDistr::dbetapr(2^(x[i]), df1[i]/2, df2/2, df2/df1[i] * scale, log)))
    } else {
      return(vector("numeric", 0))
    }
  }

  p_bayesprot_scaled_f <<- function(q, df2, scale, lower.tail = TRUE, log.p = FALSE) {
    #print(paste("p", df2, scale))
    if (length(q) > 0) {
      return(sapply(1:length(q), function(i) extraDistr::pbetapr(q[i], df1[i]/2, df2/2, df2/df1[i] * scale, lower.tail, log.p)))
    } else {
      return(vector("numeric", 0))
    }
  }

  q_bayesprot_scaled_f <<- function(p, df2, scale, lower.tail = TRUE, log.p = FALSE) {
    #print(paste("q", df2, scale))
    if (length(p) > 0) {
      return(sapply(1:length(p), function(i) extraDistr::qbetapr(p[i], df1[i]/2, df2/2, df2/df1[i] * scale, lower.tail, log.p)))
    } else {
      return(vector("numeric", 0))
    }
  }

  ft <- fitdistrplus::fitdist(log2(x), "_bayesprot_scaled_f", start = list(df2 = est$df, scale = est$v), lower = c(0.00001, 0.00001), ...)
  if (plots == T) plot(ft)

  return(list(
    v0 = est$v, df0 = est$df,
    v = ft$estimate[["scale"]], df = ft$estimate[["df2"]]
  ))
}


#' Our squeezeVar function - works for small DF unlike Limma squeezeVar, but slower
#'
#' @export
squeeze_var <- function(v, df, use.deconvolution = TRUE, plots = FALSE) {
  if (use.deconvolution) {
    return(dist_sf_with_fixed_df1_fitdistrplus(v, df, plots = plots))
  } else {
    est <- dist_invchisq(sapply(1:length(v), function(i) extraDistr::rinvchisq(1, df[i], v[i])), plots = plots)
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


# dist_invchisq_fitdistrplus <- function(x, plot = F, ...) {
#   est <- summarise_normal_robust(x)
#   s <- est$s
#   est <- list(v = est$m, df = Inf)
#
#   if (s > 0.0000001) {
#     d_bayesprot_invchisq <<- function(x, log_nu, log_tau, log = FALSE) {
#       if (length(x) > 0) {
#         return(extraDistr::dinvchisq(x, exp(log_nu), exp(log_tau), log))
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
#     ft <- fitdistrplus::fitdist(x, "_bayesprot_invchisq", start = list(log_nu = 0, log_tau = log(est$v)), ...)
#     if (plot) plot(ft)
#
#     est <- list(v = exp(ft$estimate[["log_tau"]]), df = exp(ft$estimate[["log_nu"]]))
#   }
#
#   return(est)
# }
#
#dist_invchisq <- function(x, plot = F) {
#  return(dist_invchisq_fitdistrplus(x, plot, method = "mge", gof = "KS", optim.method = "BFGS"))
#}
#
# dist_sf_with_fixed_df1_fitdistrplus <- function(x, df1, plot = F, ...) {
#   # first fit a _sample_ of x to an inverse chi squared distribution
#   est <- dist_invchisq_fitdistrplus(sapply(1:length(x), function(i) extraDistr::rinvchisq(1, df1[i], x[i])), ...)
#
#   # this then seeds the fit to a scaled F distribution
#   d_bayesprot_scaled_f <<- function(x, log_df2, log_scale, log = FALSE) {
#     if (length(x) > 0) {
#       return(sapply(1:length(x), function(i) extraDistr::dbetapr(x[i], df1[i]/2, exp(log_df2)/2, exp(log_df2)/df1[i] * exp(log_scale), log)))
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
#   ft <- fitdistrplus::fitdist(x, "_bayesprot_scaled_f", start = list(log_df2 = log(est$df), log_scale = log(est$v)), ...)
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
