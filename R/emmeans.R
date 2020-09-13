#' @export
recover_data.MCMCglmm_seaMass = function(object, data, ...) {
  attr(data, "call") <- object$Fixed
  attr(data, "terms") <- delete.response(terms(object$Fixed$formula))
  attr(data, "predictors") <- emmeans::.all.vars(delete.response(attr(data, "terms")))
  return(data)
}

#' @export
emm_basis.MCMCglmm_seaMass = function(object, trms, xlev, grid, ...) {
  nobs.MCMCglmm = function(object, ...) 1   # prevents warning about nobs
  m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
  X = model.matrix(trms, m, contrasts.arg = NULL)
  Sol = as.matrix(object$Sol)[, seq_len(object$Fixed$nfl), drop = F] # toss out random effects if included
  bhat = apply(Sol, 2, mean)
  V = cov(Sol)
  misc = list()
  if (object$family[1] == "poisson" || object$family[1] == "cenpoisson") misc <- emmeans::.std.link.labels(list(link = "log"), misc)
  return(list(
    X = X,
    bhat = bhat,
    nbasis = matrix(NA),
    V = V,
    dffun = function(k, dfargs) Inf,
    dfargs = list(),
    misc = misc, post.beta = Sol)
  )
}
