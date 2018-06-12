# Job Schedule Class for SGE, PBS and SLURM HPC systems


setClass("ScheduleHPC",
  representation(
    batch = "numeric",
    normChain = "numeric",
    modelChain = "numeric"
  ),
  prototype(
    batch = 100,
    normChain = 10,
    modelChain = 100
  )
)

setClass("SLURM",
  representation(
    check="charactor"
  ),
  contains = "ScheduleHPC"
)

setClass("PBS",
  representation(
    check="charactor"
  ),
  contains = "ScheduleHPC"
)

setClass("SGE",
  representation(
    check="charactor"
  ),
  contains = "ScheduleHPC"
)

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
  cat("#SBATCH --mem=125G\n")
  cat("#SBATCH -p cpu\n")
  cat("#SBATCH -N 1\n")
  cat("#SBATCH -c 14\n")
  cat(sprintf("#SBATCH --array=1-%d\n",jobID))
  cat("srun sh norm/norm-job$SLURM_ARRAY_TASK_ID.sh\n")
  
  sink()
  
  system("chmod u+x bayesprot-norm.sh")
}

exposuresGen <- function()
{
  sink('bayesprot-exposures.sh')

  cat("#!/bin/bash\n")
  cat("#SBATCH -J Exposures\n")
  cat("#SBATCH --export=all\n")
  cat("#SBATCH -o exposures/out/out-%A_task-%a.out\n")
  cat("#SBATCH -e exposures/error/error-%A_task-%a.out\n")
  cat("#SBATCH --mem=64G\n")
  cat("#SBATCH -p cpu\n")
  cat("#SBATCH -N 1\n") 
  cat("#SBATCH -c 14\n")
  cat("cd exposures/results\n")
  cat("srun Rscript ../../exposures.R HPC\n")
  
  sink()
  
  system("chmod u+x bayesprot-exposures.sh")
}

# Default Values batch = 100, model_chain = 100, num_parallel_jobs = 6
modelGen <- function(batch, modelChain,numPjobs,totalJobs)
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
    cat("#SBATCH --mem=125G\n")
    cat("#SBATCH -p cpu\n")
    cat("#SBATCH -N 1\n")
    cat("#SBATCH -c 14\n")
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

#   for (i in c(1:jobID))
#   {
#     if (idx == 1)
#     {
#       jobID = jobID + 1
#       fname=paste("model/model-job-batch",idx,".sh",sep="",collapse="")
#       sink(fname)
#       cat("#!/bin/bash\n")
#     }
#     cat(sprintf("sh model/model-job%d.sh\n",i))
#     if (idx == nrun || i == N)
#     {
#       sink()
#       idx = 1
#     }
#     else
#     {
#       idx = idx + 1
#     }
#   }

    sink('bayesprot-model.sh')

    cat("#!/bin/bash\n")
    cat("#SBATCH -J Model\n")
    cat("#SBATCH --export=all\n")
    cat("#SBATCH -o model/out/out-%A_task-%a.out\n")
    cat("#SBATCH -e model/error/error-%A_task-%a.out\n")
    cat("#SBATCH --mem=125G\n")
    cat("#SBATCH -p cpu\n")
    cat("#SBATCH -N 1\n")
    cat("#SBATCH -c 14\n")
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
    cat("#SBATCH --mem=125G\n")
    cat("#SBATCH -p cpu\n")
    cat("#SBATCH -N 1\n")
    cat("#SBATCH -c 14\n")
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
    cat("#SBATCH --mem=125G\n")
    cat("#SBATCH -p cpu\n")
    cat("#SBATCH -N 1\n")
    cat("#SBATCH -c 14\n")
    cat(sprintf("#SBATCH --array=1-%d\n",totalJobs))
    cat("srun sh model/plots-job-batch$SLURM_ARRAY_TASK_ID.sh\n")

    sink()
  }

  system("chmod u+x bayesprot-plots.sh")
}
