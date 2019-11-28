get_by_key <- function(obj, key = 1) {
  if (is.null(obj) || is.null(key) || key == 0) {
    return(list(index = 0, key = "NULL", value = NULL))
  } else {
    if (is.character(key)) key <- match(key, names(obj))
    i = (key - 1) %% length(obj) + 1
    return(list(index = i, key = names(obj)[i], value = obj[[i]]))
  }
}


rbindlists <- function(l) {
  for (j in names(l[[1]])) l[[1]][[j]] <- rbindlist(lapply(1:length(l), function(i) l[[i]][[j]]))
  return(l[[1]])
}


batch_split <- function(DT, columns, n) {
  DT.t <- DT[, mget(columns)]
  DT.t[, BatchID := do.call(function(...) paste(..., sep = "."), .SD)]
  DT.t[, BatchID := factor(BatchID, levels = unique(BatchID))]
  nbatch <- ceiling(length(unique(DT.t$BatchID)) / n)
  levels(DT.t$BatchID) <- rep(1:nbatch, each = n)[1:nlevels(DT.t$BatchID)]

  return(split(DT, DT.t$BatchID, drop = T, keep.by = F))
}


parallel_lapply <- function(
  inputs,
  func,
  nthread = 0,
  pred = rep(1, length(inputs)),
  .packages = "seaMass"
) {
  # retreive function arguments from the parent environment
  func.args <- sapply(names(formals(func)), function(arg) as.name(arg), simplify = F, USE.NAMES = T)
  for(n in names(func.args)) if (n != "input") assign(n, get(n, parent.frame(n = 1)))

  pb <- txtProgressBar(max = sum(pred), style = 3)
  setTxtProgressBar(pb, 0)
  if (nthread <= 1) {
    # if nthread is 1 then turn off all multiprocessing
    if (nthread == 1) {
      dt.threads <- data.table::getDTthreads()
      data.table::setDTthreads(1)
      fst.threads <- fst::threads_fst()
      fst::threads_fst(1)
    }

    # sequential
    outputs <- lapply(seq_along(inputs), function(i) {
      input <- inputs[[i]]
      output <- do.call("func", func.args)
      setTxtProgressBar(pb, getTxtProgressBar(pb) + pred[i])
      return(output)
    })

    if (nthread == 1) {
      data.table::setDTthreads(dt.threads)
      fst::threads_fst(fst::threads_fst())
    }
  } else {
    # restart cluster if different number of processes requested
    if (!is.null(parallel::getDefaultCluster()) && length(parallel::getDefaultCluster()) != nthread) {
      parallel::stopCluster(parallel::getDefaultCluster())
    }

    if (is.null(parallel::getDefaultCluster())) {
      cl <- parallel::makeCluster(nthread)
      doSNOW::registerDoSNOW(cl)
      parallel::setDefaultCluster(cl)
    }

    # parallel
    outputs <- foreach::foreach(
      input = iterators::iter(inputs),
      .packages = .packages,
      .export = c("func", names(func.args)[func.args != "input"]),
      .noexport = "output",
      .options.snow = list(progress = function(n, i) setTxtProgressBar(pb, getTxtProgressBar(pb) + pred[i]))
    ) %dorng% {
      data.table::setDTthreads(1)
      library(data.table)
      options(fst_threads = 1)
      library(fst)

      return(do.call("func", func.args))
    }
  }
  setTxtProgressBar(pb, sum(pred))
  close(pb)

  return(outputs)
}


stop_parallel <- function() {
  cl <- parallel::getDefaultCluster()
  if (!is.null(cl)) parallel::stopCluster(cl)
  return(cl)
}

