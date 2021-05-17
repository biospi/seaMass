#' Probability density function of a location-scale t distribution
#'
#' @export
dlst <- function(x, m, s, df, log = FALSE) {
  return(extraDistr::dlst(x, df, m, s, log))
}


#' Cumulative distribution function (CDF) of a location-scale t distribution
#'
#' @export
plst <- function(q, m, s, df, lower.tail = TRUE, log.p = FALSE) {
  return(extraDistr::plst(q, df, m, s, lower.tail, log.p))
}


#' Quantile function (inverse CDF) of a location-scale t distribution
#'
#' @export
qlst <- function(p, m, s, df, lower.tail = TRUE, log.p = FALSE) {
  return(extraDistr::qlst(p, df, m, s, lower.tail, log.p))
}


# random number generation for a location-scale t distribution
#'
#' @export
rlst <- function(n, m, s, df) {
  return(extraDistr::rlst(n, df, m, s))
}


#' Probability density function of an inverse nakagami distribution
#'
#' @export
dinaka <- function(x, s, df) {
  return(2 * x^-3 * dgamma(x^-2, 0.5*df, 0.5*df*s*s))
}


#' Cumulative distribution function (CDF) of an inverse nakagami distribution
#'
#' @export
pinaka <- function(x, s, df) {
  return(1 - pgamma(x^-2, 0.5*df, 0.5*df*s*s))
}


#' Quantile function (inverse CDF) of an inverse nakagami distribution
#'
#' @export
qinaka <- function(p, s, df) {
  return(1 / sqrt(qgamma(1-p, 0.5*df, 0.5*df*s*s)))
}


# random number generation for an inverse nakagami distribution
#'
#' @export
rinaka <- function(n, s, df) {
  return(1 / sqrt(rgamma(n, 0.5*df, 0.5*df*s*s)))
}


#' summarise using median/mad
#'
#' @export
dist_normal_robust <- function(value) {
  return(list(m = median(value, na.rm = T), s = mad(value, na.rm = T), df = Inf))
}


#' fit normal distribution with fitdistrplus
#'
#' @export
dist_samples_robust_normal <- function(chain, sample, value) {
  return(c(dist_normal_robust(value), list(rhat = rhat(chain, sample, value, T))))
}


#' rhat
#'
#' @export
rhat <- function(chain, sample, value, transform = FALSE) {
  if (length(unique(chain)) > 1) {
    DT <- data.table::dcast(data.table(sample = sample, chain = chain, value = value), sample ~ chain, value.var = "value")
    DT[, sample := NULL]
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
  est <- dist_normal_robust(value)

  if (!is.na(est$s) && est$s > 1e-16) {
    tryCatch({
      ft <- fitdistrplus::fitdist(as.vector(value), "norm", start = list(mean = est$m, sd = est$s), ...)
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
dist_samples_normal <- function(chain, sample, value, ...) {
  return(c(dist_normal(value, method = "mge", gof = "CvM", ...), list(rhat = rhat(chain, sample, value, T))))
}


#' fit location-scale t distribution with fitdistrplus
#'
#' @export
dist_lst <- function(value, min.df = 0, max.df = Inf, plots = FALSE, ...) {
  est <- dist_normal_robust(value)
  est$df <- Inf

  if (!is.na(est$s) && est$s > 1e-16) {
   d_seaMass_lst <<- function(x, log_df, mu, log_sigma, log = FALSE) {
      if (length(x) > 0) {
        return(extraDistr::dlst(x, exp(log_df), mu, exp(log_sigma), log))
      } else {
        return(vector("numeric", 0))
      }
    }

    p_seaMass_lst <<- function(q, log_df, mu, log_sigma, lower.tail = TRUE, log.p = FALSE) {
      if (length(q) > 0) {
        return(extraDistr::plst(q, exp(log_df), mu, exp(log_sigma), lower.tail, log.p))
      } else {
        return(vector("numeric", 0))
      }
    }

    q_seaMass_lst <<- function(p, log_df, mu, log_sigma, lower.tail = TRUE, log.p = FALSE) {
      if (length(p) > 0) {
        return(extraDistr::qlst(p, exp(log_df), mu, exp(log_sigma), lower.tail, log.p))
      } else {
        return(vector("numeric", 0))
      }
    }

    tryCatch({
      ft <- fitdistrplus::fitdist(as.vector(value), "_seaMass_lst", start = list(log_df = log(min.df + 1), mu = est$m, log_sigma = log(est$s)), lower = c(log(min.df), -Inf, -Inf), upper = c(log(max.df), Inf, Inf), ...)
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
dist_samples_lst <- function(chain, sample, value, ...) {
  return(c(dist_lst(value, method = "mge", gof = "CvM", ...), list(rhat = rhat(chain, sample, value, T))))
}


#' fit location-scale t distribution for ash (minimum df = 2)
#'
#' @export
dist_samples_lst_ash <- function(chain, sample, value, ...) {
  return(c(dist_lst(value, method = "mge", gof = "CvM", min.df = 2, max.df = 100, ...), list(rhat = rhat(chain, sample, value, T))))
}


#' fit inverse nakagami distribution with fitdistrplus, note our input is stdevs not variances!
#'
#' @export
dist_inaka <- function(value, plots = FALSE, ...) {
  est <- dist_normal_robust(log(value^2))
  est <- list(s = sqrt(exp(est$m)), df = 2.0 * est$s^-2)

  if (!is.na(est$df) && !is.infinite(est$df)) {
    d_seaMass_invchisq <<- function(x, log_df, log_v, log = F) {
      if (length(x) > 0) {
        return(extraDistr::dinvchisq(x, exp(log_df), exp(log_v), log))
      } else {
        return(vector("numeric", 0))
      }
    }

    p_seaMass_invchisq <<- function(q, log_df, log_v, lower.tail = T, log.p = F) {
      if (length(q) > 0) {
        return(extraDistr::pinvchisq(q, exp(log_df), exp(log_v), lower.tail, log.p))
      } else {
        return(vector("numeric", 0))
      }
    }

    q_seaMass_invchisq <<- function(p, log_df, log_v, lower.tail = T, log.p = F) {
      if (length(p) > 0) {
        return(extraDistr::qinvchisq(p, exp(log_df), exp(log_v), lower.tail, log.p))
      } else {
        return(vector("numeric", 0))
      }
    }

    tryCatch({
      ft <- fitdistrplus::fitdist(as.vector(value^2), "_seaMass_invchisq", start = list(log_df = log(est$df), log_v = log(est$s^2)), ...)
      if (plots == T) plot(ft)
      est <- list(s = sqrt(exp(ft$estimate[["log_v"]])), df = exp(ft$estimate[["log_df"]]))
    }, error = function(e) {
      warning("'dist_inaka' fitting failed, falling back to robust approximation.")
    })
  }

  return(est)
}


#' fit scaled inverse chi distribution
#'
#' Clarke et al 2012 ("A fast robust method for fitting gamma distributions") are right - CvM is about the only robust way to fit gamma distributions!!
#'
#' @export
dist_samples_inaka <- function(chain, sample, value, ...) {
  return(c(dist_inaka(value, method = "mge", gof = "CvM", ...), list(rhat = rhat(chain, sample, value, T))))
}


#' fit scaled F distribution with per-datapoint fixed df1 with fitdistrplus
#'
#' @export
dist_sf_with_fixed_df1_fitdistrplus <- function(value, df1, plots = FALSE, ...) {
  # first fit a _sample_ of value to an inverse nakagami distribution
  est0 <- dist_inaka(rinaka(length(value), value, df1), plots = plots)

  # this then seeds the fit to a scaled F distribution
  d_seaMass_scaled_f <<- function(x, log_df2, log_scale, log = FALSE) {
    if (length(x) > 0) {
      return(sapply(1:length(x), function(i) extraDistr::dbetapr(2^(x[i]), df1[i]/2, exp(log_df2)/2, exp(log_df2)/df1[i] * exp(log_scale), log)))
    } else {
      return(vector("numeric", 0))
    }
  }

  p_seaMass_scaled_f <<- function(q, log_df2, log_scale, lower.tail = TRUE, log.p = FALSE) {
    if (length(q) > 0) {
      return(sapply(1:length(q), function(i) extraDistr::pbetapr(q[i], df1[i]/2, exp(log_df2)/2, exp(log_df2)/df1[i] * exp(log_scale), lower.tail, log.p)))
    } else {
      return(vector("numeric", 0))
    }
  }

  q_seaMass_scaled_f <<- function(p, log_df2, log_scale, lower.tail = TRUE, log.p = FALSE) {
    if (length(p) > 0) {
      return(sapply(1:length(p), function(i) extraDistr::qbetapr(p[i], df1[i]/2, exp(log_df2)/2, exp(log_df2)/df1[i] * exp(log_scale), lower.tail, log.p)))
    } else {
      return(vector("numeric", 0))
    }
  }

  tryCatch({
    ft <- fitdistrplus::fitdist(log2(as.vector(value)^2), "_seaMass_scaled_f", start = list(log_df2 = log(est0$df), log_scale = log(est0$s^2)), ...)
    if (plots == T) plot(ft)
    est <- list(s = sqrt(exp(ft$estimate[["log_scale"]])), df = exp(ft$estimate[["log_df2"]]))
  }, error = function(e) {
    warning("'dist_sf_with_fixed_df1_fitdistrplus' fitting failed, falling back to 'dist_inaka' fit to random sample.")
    est <- est0
  })

  return(list(s0 = est0$s, df0 = est0$df, s = est$s, df = est$df))
}


#' Our squeezeVar function - works for small DF unlike Limma squeezeVar, but slower
#'
#' @export
squeeze_stdev <- function(s, df, use.deconvolution = TRUE, ...) {
  if (use.deconvolution) {
    return(dist_sf_with_fixed_df1_fitdistrplus(s, df, ...))
  } else {
    est <- dist_inaka(sapply(1:length(s), function(i) rinaka(1, s[i], df[i])), ...)
    return(list(s0 = est$s, df0 = est$df, s = est$s, df = est$df))
  }
}


#' Limma squeezeVar function
#'
#' @export
squeeze_stdev_limma <- function(s, df) {
  ft <- limma::fitFDist(s^2, df)
  return(list(s = sqrt(ft$scale), df = ft$df2))
}

