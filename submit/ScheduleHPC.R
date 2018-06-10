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

# Defaul Values batch = 100, norm_chain = 10, num_parallel_jobs = 6
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
    cat(sprintf("exec Rscript ../../norm.R %d %d %d > ../out/out_std_%d_%d.out &\n",batchNum[i],vecChain[i],normChain,batchNum[i],vecChain[i]))

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

# Defaul Values batch = 100, model_chain = 100, num_parallel_jobs = 6
modelGen <- function(batch, modelChain,numPjobs)
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
    cat(sprintf("exec Rscript ../../model.R %d %d %d > ../out/out_std_%d_%d.out &\n",batchNum[i],vecChain[i],modelChain,batchNum[i],vecChain[i]))

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
  
  system("chmod u+x bayesprot-model.sh")
}
