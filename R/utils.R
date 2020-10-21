#' @import data.table
rbindlists <- function(l) {
  for (j in names(l[[1]])) l[[1]][[j]] <- rbindlist(lapply(1:length(l), function(i) l[[i]][[j]]))
  return(l[[1]])
}


#' @import data.table
batch_split <- function(DT, columns, n, drop = FALSE, keep.by = TRUE) {
  DT[, Batch := Reduce(function(...) paste(..., sep = "."), .SD[, mget(columns)])]
  DT[, Batch := factor(Batch, levels = unique(Batch))]
  levels(DT$Batch) <- rep(1:n, each = ceiling(nlevels(DT$Batch) / n))[1:nlevels(DT$Batch)]
  DT[, Batch := as.integer(Batch)]
  DT <- split(DT, by = "Batch", keep.by = keep.by, drop = T)
  if (drop) DT <- lapply(DT, droplevels)
  return(DT)
}


#' @import foreach
#' @import doRNG
parallel_lapply <- function(items, func, nthread = 0, pred = rep(1, length(items)), .packages = "seaMass") {
  # retreive function arguments from the parent environment
  func.args <- sapply(names(formals(func)), function(arg) as.name(arg), simplify = F, USE.NAMES = T)
  for(n in names(func.args)) if (n != "item") assign(n, get(n, parent.frame(n = 1)))

  if (length(items) > 1) {
    pb <- txtProgressBar(max = sum(pred), style = 3)
    setTxtProgressBar(pb, 0)
  }

  # check if a worker
  if (exists(".master.pid") || nthread == 0) {
    # sequential
    outputs <- lapply(seq_along(items), function(i) {
      item <- items[[i]]
      output <- do.call("func", func.args)
      if (length(items) > 1) setTxtProgressBar(pb, getTxtProgressBar(pb) + pred[i])
      return(output)
    })

    if (length(items) > 1) {
      setTxtProgressBar(pb, sum(pred))
      close(pb)
    }
  } else {
    # restart cluster EVERY TIME just to stop memory leaks (crap GC or other problem?)
    if (!is.null(parallel::getDefaultCluster())) {
        parallel::stopCluster(parallel::getDefaultCluster())
    }

    cl <- parallel::makeCluster(ifelse(length(items) < nthread, length(items), nthread))
    #cl <- parallel::makeCluster(ifelse(length(items) < nthread, length(items), nthread), outfile = "")
    doSNOW::registerDoSNOW(cl)
    parallel::setDefaultCluster(cl)

    # parallel
    if (length(items) > 1) {
      progress <- function(n, i) setTxtProgressBar(pb, getTxtProgressBar(pb) + pred[i])
    } else {
      progress <- NULL
    }

    .master.pid <- Sys.getpid()
    outputs <- foreach::foreach(
      item = iterators::iter(items),
      .errorhandling = "pass",
      .packages = .packages,
      .export = c("func", names(func.args)[func.args != "item"], ".master.pid"),
      .options.snow = list(progress = progress),
      .verbose = F
    ) %dorng% {
      data.table::setDTthreads(1)
      library(data.table)
      options(fst_threads = 1)
      library(fst)

      return(do.call("func", func.args))
    }

    parallel::stopCluster(parallel::getDefaultCluster())

    if (length(items) > 1) {
      setTxtProgressBar(pb, sum(pred))
      close(pb)
    }

    fatal <- F
    foreach::foreach(item = iterators::iter(items), output = iterators::iter(outputs)) %do% {
      if (inherits(output, "simpleError")) {
        message("[", Sys.time(), "] ERROR: parallel_lapply item:")
        print(item)
        message(paste0("[", Sys.time(), "]  message=", output$message))
        message(deparse(output$call))
        fatal <- T
      }
    }
    if (fatal) stop()
  }

  return(outputs)
}


stop_parallel <- function() {
  cl <- parallel::getDefaultCluster()
  if (!is.null(cl)) parallel::stopCluster(cl)
  return(cl)
}


