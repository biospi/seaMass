# Job Schedule Class for SGE, PBS and SLURM HPC systems

require(methods)

# Abstract Class
setClass("ScheduleHPC",
  representation(
    batch = "numeric",
    quantChain = "numeric",
    modelChain = "numeric",
    path = "character"
  ),
  prototype(
    batch = 10,
    quantChain = 10,
    modelChain = 100,
    path = "."
  )
)

# Derived Classes
setClass("SLURM",
  representation(
    cpuNum = "numeric",
    node = "numeric",
    mem = "character",
    himem = "character",
    longQue = "character",
    shortQue = "character",
    totalJobs = "numeric",
    lowCPUNum = "numeric"
  ),
  prototype
  (
    cpuNum = 14,
    node = 1,
    mem = "3G",
    himem = "16G",
    longQue = "cpu",
    shortQue = "serial",
    totalJobs = 1,
    lowCPUNum = 6
  ),
  contains = "ScheduleHPC"
)

setClass("PBS",
  representation(
  ),
  contains = "ScheduleHPC"
)

setClass("SGE",
  representation(
  ),
  contains = "ScheduleHPC"
)

# Class Function defined
setGeneric("quantHPC",
  function(object)
  {
    standardGeneric("quantHPC")
  }
)

setGeneric("exposuresHPC",
  function(object)
  {
    standardGeneric("exposuresHPC")
  }
)

setGeneric("modelHPC",
  function(object)
  {
    standardGeneric("modelHPC")
  }
)

setGeneric("plotsHPC",
  function(object)
  {
    standardGeneric("plotsHPC")
  }
)


setGeneric("outputHPC",
  function(object)
  {
    standardGeneric("outputHPC")
  }
)


setGeneric("genJobFileHPC",
  function(object)
  {
    standardGeneric("genJobFileHPC")
  }
)

setMethod("quantHPC", signature(object = "SLURM"), function(object)
  {
    N <- object@batch * object@quantChain
    batchNum <- vector(mode="numeric",length=N)
    vecChain <- vector(mode="numeric",length=N)

    idx = 1
    for (i in seq (0, object@batch-1, 1))
    {
      for (j in seq(0, object@quantChain-1, 1))
      {
        batchNum[idx] <- i
        vecChain[idx] <- j
        idx <- idx + 1
      }
    }

    idx <- 1
    jobID = 0
    for (i in 1:N)
    {
      if (idx == 1)
      {
        jobID = jobID + 1
        fname=paste(file.path(object@path,"quant","quant-job"),jobID,".sh",sep="",collapse="")
        sink(fname)
        cat("cd quant/results\n")
      }
      cat(sprintf("exec Rscript ../../quant.R %d %d %d &> ../out_std_%d_%d.out &\n",
                    batchNum[i],vecChain[i],object@quantChain,batchNum[i],vecChain[i])
                  )
      if ((idx == object@cpuNum) || i == N)
      {
        cat("wait\n")
        sink()
        idx = 1
      }
      else
      {
        idx = idx + 1
      }
    }
    sink(file.path(object@path,"bayesprot-quant.sh"))

    cat("#!/bin/bash\n")
    cat("#SBATCH -J Quant\n")
    cat("#SBATCH --export=all\n")
    cat("#SBATCH -o quant/out-%A_task-%a.out\n")
    cat("#SBATCH -e quant/error-%A_task-%a.out\n")
    cat(sprintf("#SBATCH --mem-per-cpu=%s\n",object@mem))
    cat(sprintf("#SBATCH -p %s\n",object@longQue))
    cat(sprintf("#SBATCH -N %d\n",object@node))
    cat(sprintf("#SBATCH -c %d\n",object@cpuNum))
    cat(sprintf("#SBATCH --array=1-%d\n",jobID))
    cat("srun sh quant/quant-job$SLURM_ARRAY_TASK_ID.sh\n")

    sink()

    system(paste("chmod u+x",file.path(object@path,"bayesprot-quant.sh")))
  }
)


setMethod("exposuresHPC", signature(object = "SLURM"), function(object)
  {
    singleCPU = 1
    singleNode = 1

    sink(file.path(object@path,"bayesprot-exposures.sh"))

    cat("#!/bin/bash\n")
    cat("#SBATCH -J Exposures\n")
    cat("#SBATCH --export=all\n")
    cat("#SBATCH -o exposures/out-%A_task-%a.out\n")
    cat("#SBATCH -e exposures/error-%A_task-%a.out\n")
    cat(sprintf("#SBATCH --mem=%s\n",object@himem))
    cat(sprintf("#SBATCH -p %s\n",object@shortQue))
    cat(sprintf("#SBATCH -N %d\n",singleNode))
    cat(sprintf("#SBATCH -c %s\n",singleCPU))
    cat("cd exposures/results\n")
    cat("srun Rscript ../../exposures.R\n")

    sink()

    system(paste("chmod u+x",file.path(object@path,"bayesprot-exposures.sh")))
    }
)


setMethod("modelHPC", signature(object = "SLURM"), function(object)
  {
    N <- object@batch * object@modelChain
    batchNum <- vector(mode="numeric",length=N)
    vecChain <- vector(mode="numeric",length=N)

    idx = 1
    for (i in seq (0, object@batch-1, 1))
    {
      for (j in seq(0, object@modelChain-1, 1))
      {
        batchNum[idx] <- i
        vecChain[idx] <- j
        idx <- idx + 1
      }
    } 

    idx <- 1
    jobID = 0
    for (i in 1:N)
    {
      if (idx == 1)
      {
        jobID = jobID + 1
        fname=paste(file.path(object@path,"model","model-job"),jobID,".sh",sep="",collapse="")
        sink(fname)
        cat("cd model/results\n")
      }
      cat(sprintf("exec Rscript ../../model.R %d %d %d &> ../out_std_%d_%d.out &\n",
                    batchNum[i],vecChain[i],object@modelChain,batchNum[i],vecChain[i])
                 )

      if (idx == object@cpuNum || i == N)
      {
        cat("wait\n")
        sink()
        idx = 1
      }
      else
      {
        idx = idx + 1
      }
    }

    if (object@totalJobs == 1)
    {
      sink(file.path(object@path,"bayesprot-model.sh"))
    
      cat("#!/bin/bash\n")
      cat("#SBATCH -J Model\n")
      cat("#SBATCH --export=all\n")
      cat("#SBATCH -o model/out-%A_task-%a.out\n")
      cat("#SBATCH -e model/error-%A_task-%a.out\n")
      cat(sprintf("#SBATCH --mem-per-cpu=%s\n",object@mem))
      cat(sprintf("#SBATCH -p %s\n",object@longQue))
      cat(sprintf("#SBATCH -N %d\n",object@node))
      cat(sprintf("#SBATCH -c %d\n",object@cpuNum))
      cat(sprintf("#SBATCH --array=1-%d\n",jobID))
      cat("srun sh model/model-job$SLURM_ARRAY_TASK_ID.sh\n")

      sink()
    }
    else
    {
      nrun <- round(jobID/object@totalJobs)
      extra <- jobID%%object@totalJobs
      idx <- 0
      for (i in c(1:object@totalJobs))
      {
        fname=paste(file.path(object@path,"model","model-batch-job"),i,".sh",sep="",collapse="")
        sink(fname)
        cat("#!/bin/bash\n")
        for (j in 1:nrun)
        {
          idx = idx + 1
          cat(sprintf("sh model/model-job%d.sh\n",idx))
        }
        if (extra > 0)
        {
          idx = idx + 1
          cat(sprintf("sh model/model-job%d.sh\n",idx))
          extra = extra - 1
        }
        sink()
        if (idx > jobID)
        {
          print(paste("ERROR: JobIB greater than number of subjobs!!! idx: ", idx, " extra: ",extra))
        }
      }

      sink(file.path(object@path,"bayesprot-model.sh"))

      cat("#!/bin/bash\n")
      cat("#SBATCH -J Model\n")
      cat("#SBATCH --export=all\n")
      cat("#SBATCH -o model/out-%A_task-%a.out\n")
      cat("#SBATCH -e model/error-%A_task-%a.out\n")
      cat(sprintf("#SBATCH --mem-per-cpu=%s\n",object@mem))
      cat(sprintf("#SBATCH -p %s\n",object@longQue))
      cat(sprintf("#SBATCH -N %d\n",object@node))
      cat(sprintf("#SBATCH -c %d\n",object@cpuNum))
      cat(sprintf("#SBATCH --array=1-%d\n",object@totalJobs))
      cat("srun sh model/model-job-batch$SLURM_ARRAY_TASK_ID.sh\n")

      sink()
    }
    system(paste("chmod u+x",file.path(object@path,"bayesprot-model.sh")))
  }
)


setMethod("plotsHPC", signature(object = "SLURM"), function(object)
  {
    N <- object@batch 
    batchNum <- seq(0,object@batch-1,1)

    idx <- 1
    jobID = 0
    for (i in 1:N)
    {
      if (idx == 1)
      {
        jobID = jobID + 1
        fname=paste(file.path(object@path,"plots","plots-job"),jobID,".sh",sep="",collapse="")
        sink(fname)
        cat("cd plots/results\n")
      }
      cat(sprintf("exec Rscript ../../plots.R %d &> ../out_std_%d.out &\n",batchNum[i],batchNum[i]))

      if (idx == object@lowCPUNum || i == N)
      {
        cat("wait\n")
        sink()
        idx = 1
      }
      else
      {
        idx = idx + 1
      }
    }

    if (object@totalJobs == 1)
    {
      sink(file.path(object@path,"bayesprot-plots.sh"))

      cat("#!/bin/bash\n")
      cat("#SBATCH -J Plots\n")
      cat("#SBATCH --export=all\n")
      cat("#SBATCH -o plots/out-%A_task-%a.out\n")
      cat("#SBATCH -e plots/error-%A_task-%a.out\n")
      cat(sprintf("#SBATCH --mem=%dG\n",(16*object@lowCPUNum)))
      cat(sprintf("#SBATCH -p %s\n",object@longQue))
      cat(sprintf("#SBATCH -N %d\n",object@node))
      cat(sprintf("#SBATCH -c %d\n",object@cpuNum))
      cat(sprintf("#SBATCH --array=1-%d\n",jobID))
      cat("srun sh plots/plots-job$SLURM_ARRAY_TASK_ID.sh\n")

      sink()
    }
    else
    {
      nrun <- round(jobID/object@totalJobs)
      extra <- jobID%%object@totalJobs
      idx <- 0
      for (i in c(1:object@totalJobs))
      {
        fname=paste(file.path(object@path,"plots","plots-batch-job"),i,".sh",sep="",collapse="")
        sink(fname)
        cat("#!/bin/bash\n")
        for (j in 1:nrun)
        {
          idx = idx + 1
          cat(sprintf("sh plots/plots-job%d.sh\n",idx))
        }
        if (extra > 0)
        {
          idx = idx + 1
          cat(sprintf("sh plots/plot-job%d.sh\n",idx))
          extra = extra - 1
        }
        sink()
        if (idx > jobID)
        {
          print(paste("ERROR: JobIB greater than number of subjobs!!! idx: ", idx, " extra: ",extra))
        }
      }
      sink(file.path(object@path,"bayesprot-plots.sh"))

      cat("#!/bin/bash\n")
      cat("#SBATCH -J Plots\n")
      cat("#SBATCH --export=all\n")
      cat("#SBATCH -o plots/out-%A_task-%a.out\n")
      cat("#SBATCH -e plots/error-%A_task-%a.out\n")
      cat(sprintf("#SBATCH --mem=%dG\n",(16*object@lowCPUNum)))
      cat(sprintf("#SBATCH -p %s\n",object@longQue))
      cat(sprintf("#SBATCH -N %d\n",object@node))
      cat(sprintf("#SBATCH -c %d\n",object@cpuNum))
      cat(sprintf("#SBATCH --array=1-%d\n",object@totalJobs))
      cat("srun sh model/plots-job-batch$SLURM_ARRAY_TASK_ID.sh\n")

      sink()
    }
    system(paste("chmod u+x",file.path(object@path,"bayesprot-plots.sh")))
  }
)


setMethod("outputHPC", signature(object = "SLURM"), function(object)
  {
    singleCPU = 1
    singleNode = 1

    sink(file.path(object@path,"bayesprot-output.sh"))

    cat("#!/bin/bash\n")
    cat("#SBATCH -J Output\n")
    cat("#SBATCH --export=all\n")
    cat("#SBATCH -o output/out-%A_task-%a.out\n")
    cat("#SBATCH -e output/error-%A_task-%a.out\n")
    cat(sprintf("#SBATCH --mem=%s\n",object@himem))
    cat(sprintf("#SBATCH -p %s\n",object@shortQue))
    cat(sprintf("#SBATCH -N %d\n",singleNode))
    cat(sprintf("#SBATCH -c %s\n",singleCPU))
    cat("cd output/results\n")
    cat("srun Rscript ../../output.R\n")

    sink()

    system(paste("chmod u+x",file.path(object@path,"bayesprot-output.sh")))
    }
)


setMethod("genJobFileHPC", signature(object = "SLURM"), function(object)
  {
    sink(file.path(object@path,"jobSubmitScript.sh"))
    
    cat("#!/bin/bash\n")
    cat("# jobSubmitScript.sh\n\n")

    #########################################
    # Quant Array Job
    #########################################
    cat("sbatchReturn=$(sbatch  bayesprot-quant.sh)\n")
    cat("echo sbatch bayesprot-quant.sh\n")
    cat("quantJobID=$sbatchReturn\n")
    cat("quantJobID=${quantJobID//[^0-9]/}\n")
    cat("echo Quant Job Array Submitted: $quantJobID\n\n")

    #########################################
    # Exposures
    #########################################
    cat("sbatchReturn=$(sbatch --dependency=afterok:$quantJobID bayesprot-exposures.sh)\n")
    cat("echo sbatch --dependency=afterok:$quantJobID bayesprot-exposures.sh\n")
    cat("exposuresJobID=$sbatchReturn\n")
    cat("exposuresJobID=${exposuresJobID//[^0-9]/}\n")
    cat("echo Exposures Submitted: $exposuresJobID\n\n")

    #########################################
    # Model Array Job
    #########################################
    cat("sbatchReturn=$(sbatch --dependency=afterok:$exposuresJobID bayesprot-model.sh)\n")
    cat("echo sbatch --dependency=afterok:$exposuresJobID bayesprot-model.sh\n")
    cat("modelJobID=$sbatchReturn\n")
    cat("modelJobID=${modelJobID//[^0-9]/}\n")
    cat("echo Model Array Job Submitted: $modelJobID\n\n")

    #########################################
    # Plots Array Job
    #########################################
    cat("sbatchReturn=$(sbatch --dependency=afterok:$modelJobID bayesprot-plots.sh)\n")
    cat("echo sbatch --dependency=afterok:$modelJobID bayesprot-plots.sh\n")
    cat("plotsJobID=$sbatchReturn\n")
    cat("plotsJobID=${plotsJobID//[^0-9]/}\n")
    cat("echo Plots Array Job Submitted: $plotsJobID\n\n")

    #########################################
    # Output
    #########################################
    cat("sbatchReturn=$(sbatch --dependency=afterok:$plotsJobID bayesprot-output.sh)\n")
    cat("echo sbatch --dependency=afterok:$plotsJobID bayesprot-output.sh\n")
    cat("outputJobID=$sbatchReturn\n")
    cat("outputJobID=${outputJobID//[^0-9]/}\n")
    cat("echo Plots Array Job Submitted: $outputJobID\n\n")

    sink()

    system(paste("chmod u+x",file.path(object@path,"jobSubmitScript.sh")))
  }
)
