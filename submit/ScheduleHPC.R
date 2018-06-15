# Job Schedule Class for SGE, PBS and SLURM HPC systems


# Abstract Class
setClass("ScheduleHPC",
  representation(
    batch = "numeric",
    normChain = "numeric",
    modelChain = "numeric"
  ),
  prototype(
    batch = 10,
    normChain = 10,
    modelChain = 100
  )
)

# Derived Classes
setClass("SLURM",
  representation(
    cpuNum = "numeric",
    node = "numeric",
    mem = "character",
    longQue = "character",
    shortQue = "character",
    totalJobs = "numeric"
  ),
  prototype
  (
    cpuNum = 14,
    node = 1,
    mem = "3G",
    longQue = "cpu",
    shortQue = "serial",
    totalJobs = 1
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
setGeneric("normHPC",
  function(object)
  {
    standardGeneric("normHPC")
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


setGeneric("genJobFileHPC",
  function(object)
  {
    standardGeneric("genJobFileHPC")
  }
)

setMethod("normHPC", signature(object = "SLURM"), function(object)
  {
    N <- object@batch * object@normChain
    batchNum <- vector(mode="numeric",length=N)
    vecChain <- vector(mode="numeric",length=N)

    idx = 1
    for (i in seq (0, object@batch-1, 1))
    {
      for (j in seq(0, object@normChain-1, 1))
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
        fname=paste("norm/norm-job",jobID,".sh",sep="",collapse="")
        sink(fname)
        cat("cd norm/results\n")
      }
      cat(sprintf("exec Rscript ../../norm.R %d %d %d &> ../out/out_std_%d_%d.out &\n",
                    batchNum[i],vecChain[i],object@normChain,batchNum[i],vecChain[i])
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
    sink('bayesprot-norm.sh')

    cat("#!/bin/bash\n")
    cat("#SBATCH -J Norm\n")
    cat("#SBATCH --export=all\n")
    cat("#SBATCH -o norm/out/out-%A_task-%a.out\n")
    cat("#SBATCH -e norm/error/error-%A_task-%a.out\n")
    cat(sprintf("#SBATCH --mem-per-cpu=%s\n",object@mem))
    cat(sprintf("#SBATCH -p %s\n",object@longQue))
    cat(sprintf("#SBATCH -N %d\n",object@node))
    cat(sprintf("#SBATCH -c %d\n",object@cpuNum))
    cat(sprintf("#SBATCH --array=1-%d\n",jobID))
    cat("srun sh norm/norm-job$SLURM_ARRAY_TASK_ID.sh\n")

    sink()

    system("chmod u+x bayesprot-norm.sh")
  }
)


setMethod("exposuresHPC", signature(object = "SLURM"), function(object)
  {
    singleCPU = 1
    singleNode = 1

    sink('bayesprot-exposures.sh')

    cat("#!/bin/bash\n")
    cat("#SBATCH -J Exposures\n")
    cat("#SBATCH --export=all\n")
    cat("#SBATCH -o exposures/out/out-%A_task-%a.out\n")
    cat("#SBATCH -e exposures/error/error-%A_task-%a.out\n")
    cat(sprintf("#SBATCH --mem=16G\n"))
    cat(sprintf("#SBATCH -p %s\n",object@shortQue))
    cat(sprintf("#SBATCH -N %d\n",singleNode))
    cat(sprintf("#SBATCH -c %s\n",singleCPU))
    cat("cd exposures/results\n")
    cat("srun Rscript ../../exposures.R HPC\n")

    sink()

    system("chmod u+x bayesprot-exposures.sh")
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
        fname=paste("model/model-job",jobID,".sh",sep="",collapse="")
        sink(fname)
        cat("cd model/results\n")
      }
      cat(sprintf("exec Rscript ../../model.R %d %d %d &> ../out/out_std_%d_%d.out &\n",
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
      sink('bayesprot-model.sh')
    
      cat("#!/bin/bash\n")
      cat("#SBATCH -J Model\n")
      cat("#SBATCH --export=all\n")
      cat("#SBATCH -o model/out/out-%A_task-%a.out\n")
      cat("#SBATCH -e model/error/error-%A_task-%a.out\n")
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
        fname=paste("model/model-job-batch",i,".sh",sep="",collapse="")
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

      sink('bayesprot-model.sh')

      cat("#!/bin/bash\n")
      cat("#SBATCH -J Model\n")
      cat("#SBATCH --export=all\n")
      cat("#SBATCH -o model/out/out-%A_task-%a.out\n")
      cat("#SBATCH -e model/error/error-%A_task-%a.out\n")
      cat(sprintf("#SBATCH --mem-per-cpu=%s\n",object@mem))
      cat(sprintf("#SBATCH -p %s\n",object@longQue))
      cat(sprintf("#SBATCH -N %d\n",object@node))
      cat(sprintf("#SBATCH -c %d\n",object@cpuNum))
      cat(sprintf("#SBATCH --array=1-%d\n",object@totalJobs))
      cat("srun sh model/model-job-batch$SLURM_ARRAY_TASK_ID.sh\n")

      sink()
    }
    system("chmod u+x bayesprot-model.sh")
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
        fname=paste("plots/plots-job",jobID,".sh",sep="",collapse="")
        sink(fname)
        cat("cd plots/results\n")
      }
      cat(sprintf("exec Rscript ../../plots.R %d &> ../out/out_std_%d.out &\n",batchNum[i],batchNum[i]))

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
      sink('bayesprot-plots.sh')

      cat("#!/bin/bash\n")
      cat("#SBATCH -J Plots\n")
      cat("#SBATCH --export=all\n")
      cat("#SBATCH -o plots/out/out-%A_task-%a.out\n")
      cat("#SBATCH -e plots/error/error-%A_task-%a.out\n")
      cat(sprintf("#SBATCH --mem-per-cpu=%s\n",object@mem))
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
        fname=paste("plots/plots-job-batch",i,".sh",sep="",collapse="")
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
      sink('bayesprot-plots.sh')

      cat("#!/bin/bash\n")
      cat("#SBATCH -J Plots\n")
      cat("#SBATCH --export=all\n")
      cat("#SBATCH -o plots/out/out-%A_task-%a.out\n")
      cat("#SBATCH -e plots/error/error-%A_task-%a.out\n")
      cat(sprintf("#SBATCH --mem-per-cpu=%s\n",object@mem))
      cat(sprintf("#SBATCH -p %s\n",object@longQue))
      cat(sprintf("#SBATCH -N %d\n",object@node))
      cat(sprintf("#SBATCH -c %d\n",object@cpuNum))
      cat(sprintf("#SBATCH --array=1-%d\n",object@totalJobs))
      cat("srun sh model/plots-job-batch$SLURM_ARRAY_TASK_ID.sh\n")

      sink()
    }
    system("chmod u+x bayesprot-plots.sh")
  }
)

setMethod("genJobFileHPC", signature(object = "SLURM"), function(object)
  {
    sink('jobSubmitScript.sh')
    
    cat("#!/bin/bash\n")
    cat("# jobSubmitScript.sh\n\n")

    #########################################
    # Norm Array Job
    #########################################
    cat("sbatchReturn=$(sbatch  bayesprot-norm.sh)\n")
    cat("echo sbatch bayesprot-norm.sh\n")
    cat("normJobID=$sbatchReturn\n")
    cat("echo Norm Job Array Submitted: $normJobID\n\n")

    #########################################
    # Exposures
    #########################################
    cat("sbatchReturn=$(sbatch --dependency=afterok:$normJobID bayesprot-exposures.sh)\n")
    cat("echo sbatch --dependency=afterok:$normJobID bayesprot-exposures.sh\n")
    cat("exposuresJobID=$sbatchReturn\n")
    cat("echo Exposures Submitted: $exposuresJobID\n\n")

    #########################################
    # Model Array Job
    #########################################
    cat("sbatchReturn=$(sbatch --dependency=afterok:$exposuresJobID bayesprot-model.sh)\n")
    cat("echo sbatch --dependency=afterok:$exposuresJobID bayesprot-model.sh\n")
    cat("modelJobID=$sbatchReturn\n")
    cat("echo Model Array Job Submitted: $modelJobID\n\n")

    #########################################
    # Plots Array Job
    #########################################
    cat("sbatchReturn=$(sbatch --dependency=afterok:$modelJobID bayesprot-plots.sh)\n")
    cat("echo sbatch --dependency=afterok:$modelJobID bayesprot-plots.sh\n")
    cat("plotsJobID=$sbatchReturn\n")
    cat("echo Plots Array Job Submitted: $plotsJobID\n\n")

    #########################################
    # Output
    #########################################
    #outputID=$(sbatch --dependency=afterok:$outputSetupJobID bayesprot-output.sh)
    #echo sbatch --dependency=afterok:$outputSetupJobID bayesprot-output.sh
    #echo Output Submitted: $outputID

    sink()

    system("chmod u+x jobSubmitScript.sh")
  }
)

## ****************************************************************************
## ****************************************************************************
## ****************************************************************************

# Default Values batch = 100, norm_chain = 10, num_parallel_jobs = 6
normGen <- function(batch, normChain,numPjobs)
{

  N <- batch * normChain
  batchNum <- vector(mode="numeric",length=N)
  vecChain <- vector(mode="numeric",length=N)

  idx = 1
  for (i in seq (0, batch-1, 1))
  {
    for (j in seq(0, normChain-1, 1))
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
      fname=paste("norm/norm-job",jobID,".sh",sep="",collapse="")
      sink(fname)
      cat("cd norm/results\n")
    }
    cat(sprintf("exec Rscript ../../norm.R %d %d %d &> ../out/out_std_%d_%d.out &\n",batchNum[i],vecChain[i],normChain,batchNum[i],vecChain[i]))

    if (idx == numPjobs || i == N)
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

  sink('bayesprot-norm.sh')
  
  cat("#!/bin/bash\n")
  cat("#SBATCH -J Norm\n")
  cat("#SBATCH --export=all\n")
  cat("#SBATCH -o norm/out/out-%A_task-%a.out\n")
  cat("#SBATCH -e norm/error/error-%A_task-%a.out\n")
  #cat("#SBATCH --mem=125G\n")
  cat("#SBATCH --mem-per-cpu=3G\n")
  cat("#SBATCH -p cpu\n")
  cat("#SBATCH -N 1\n")
  #cat("#SBATCH -c 14\n")
  cat(sprintf("#SBATCH -c %d\n",numPjobs))
  cat(sprintf("#SBATCH --array=1-%d\n",jobID))
  cat("srun sh norm/norm-job$SLURM_ARRAY_TASK_ID.sh\n")
  
  sink()
  
  system("chmod u+x bayesprot-norm.sh")
}


# Exposures Generator
exposuresGen <- function()
{
  sink('bayesprot-exposures.sh')

  cat("#!/bin/bash\n")
  cat("#SBATCH -J Exposures\n")
  cat("#SBATCH --export=all\n")
  cat("#SBATCH -o exposures/out/out-%A_task-%a.out\n")
  cat("#SBATCH -e exposures/error/error-%A_task-%a.out\n")
  cat("#SBATCH --mem=16G\n")
  cat("#SBATCH -p serial\n")
  cat("#SBATCH -N 1\n") 
  cat("#SBATCH -c 1\n")
  cat("cd exposures/results\n")
  cat("srun Rscript ../../exposures.R HPC\n")
  
  sink()
  
  system("chmod u+x bayesprot-exposures.sh")
}

# Default Values batch = 100, model_chain = 100, num_parallel_jobs = 6
modelGen <- function(batch, modelChain,numPjob,totalJobs)
{

  N <- batch * modelChain
  batchNum <- vector(mode="numeric",length=N)
  vecChain <- vector(mode="numeric",length=N)

  idx = 1
  for (i in seq (0, batch-1, 1))
  {
    for (j in seq(0, modelChain-1, 1))
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
      fname=paste("model/model-job",jobID,".sh",sep="",collapse="")
      sink(fname)
      cat("cd model/results\n")
    }
    cat(sprintf("exec Rscript ../../model.R %d %d %d &> ../out/out_std_%d_%d.out &\n",batchNum[i],vecChain[i],modelChain,batchNum[i],vecChain[i]))

    if (idx == numPjobs || i == N)
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

  if (totalJobs == 1)
  {
    sink('bayesprot-model.sh')

    cat("#!/bin/bash\n")
    cat("#SBATCH -J Model\n")
    cat("#SBATCH --export=all\n")
    cat("#SBATCH -o model/out/out-%A_task-%a.out\n")
    cat("#SBATCH -e model/error/error-%A_task-%a.out\n")
    #cat("#SBATCH --mem=125G\n")
    cat("#SBATCH --mem-per-cpu=3G\n")
    cat("#SBATCH -p cpu\n")
    cat("#SBATCH -N 1\n")
    #cat("#SBATCH -c 14\n")
    cat(sprintf("#SBATCH -c %d\n",numPjobs))
    cat(sprintf("#SBATCH --array=1-%d\n",jobID))
    cat("srun sh model/model-job$SLURM_ARRAY_TASK_ID.sh\n")

    sink()
  }
  else
  {
    nrun <- round(jobID/totalJobs)
    extra <- jobID%%totalJobs
    idx <- 0
    for (i in c(1:totalJobs))
    {
      fname=paste("model/model-job-batch",i,".sh",sep="",collapse="")
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

    sink('bayesprot-model.sh')

    cat("#!/bin/bash\n")
    cat("#SBATCH -J Model\n")
    cat("#SBATCH --export=all\n")
    cat("#SBATCH -o model/out/out-%A_task-%a.out\n")
    cat("#SBATCH -e model/error/error-%A_task-%a.out\n")
    #cat("#SBATCH --mem=125G\n")
    cat("#SBATCH --mem-per-cpu=3G\n")
    cat("#SBATCH -p cpu\n")
    cat("#SBATCH -N 1\n")
    #cat("#SBATCH -c 28\n")
    cat(sprintf("#SBATCH -c %d\n",numPjobs))
    cat(sprintf("#SBATCH --array=1-%d\n",totalJobs))
    cat("srun sh model/model-job-batch$SLURM_ARRAY_TASK_ID.sh\n")

    sink()
  }
  system("chmod u+x bayesprot-model.sh")
}

# Default Values batch = 100, num_parallel_jobs = 6
plotsGen <- function(batch, numPjobs,totalJobs)
{
  N <- batch 
  batchNum <- seq(0,batch-1,1)

  idx <- 1
  jobID = 0
  for (i in 1:N)
  {
    if (idx == 1)
    {
      jobID = jobID + 1
      fname=paste("plots/plots-job",jobID,".sh",sep="",collapse="")
      sink(fname)
      cat("cd plots/results\n")
    }
    cat(sprintf("exec Rscript ../../plots.R %d &> ../out/out_std_%d.out &\n",batchNum[i],batchNum[i]))

    if (idx == numPjobs || i == N)
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

  if (totalJobs == 1)
  {
    sink('bayesprot-plots.sh')

    cat("#!/bin/bash\n")
    cat("#SBATCH -J Plots\n")
    cat("#SBATCH --export=all\n")
    cat("#SBATCH -o plots/out/out-%A_task-%a.out\n")
    cat("#SBATCH -e plots/error/error-%A_task-%a.out\n")
    #cat("#SBATCH --mem=125G\n")
    cat("#SBATCH --mem-per-cpu=8G\n")
    cat("#SBATCH -p cpu\n")
    cat("#SBATCH -N 1\n")
    #cat("#SBATCH -c 10\n")
    cat(sprintf("#SBATCH -c %d\n",numPjobs))
    cat(sprintf("#SBATCH --array=1-%d\n",jobID))
    cat("srun sh plots/plots-job$SLURM_ARRAY_TASK_ID.sh\n")

    sink()
  }
  else
  {
    nrun <- round(jobID/totalJobs)
    extra <- jobID%%totalJobs
    idx <- 0
    for (i in c(1:totalJobs))
    {
      fname=paste("plots/plots-job-batch",i,".sh",sep="",collapse="")
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
    sink('bayesprot-plots.sh')

    cat("#!/bin/bash\n")
    cat("#SBATCH -J Plots\n")
    cat("#SBATCH --export=all\n")
    cat("#SBATCH -o plots/out/out-%A_task-%a.out\n")
    cat("#SBATCH -e plots/error/error-%A_task-%a.out\n")
    #cat("#SBATCH --mem=125G\n")
    cat("#SBATCH --mem-per-cpu=8G\n")
    cat("#SBATCH -p cpu\n")
    cat("#SBATCH -N 1\n")
    #cat("#SBATCH -c 14\n")
    cat(sprintf("#SBATCH -c %d\n",numPjobs))
    cat(sprintf("#SBATCH --array=1-%d\n",totalJobs))
    cat("srun sh model/plots-job-batch$SLURM_ARRAY_TASK_ID.sh\n")

    sink()
  }

  system("chmod u+x bayesprot-plots.sh")
}
