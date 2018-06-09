# Job Schedule Class for SGE, PBS and SLURM HPC systems


setClass("ScheduleHPC",
  representation(
    nbatch = "numeric",
    normChain = "numeric",
    modelChain = "numeric"
  ),
  prototype(
    nbatch = 100,
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


normGen <- function(nbatch, normChain,normNChain,njob)
{
  matbatch=rbind(seq(0,nbatch-1,normNChain))
  for (i in 1:normNChain)
  {
    matbatch=rbind(matbatch,seq(i,nbatch-1,normNChain))
  }

  # Number of parralel jobs to run on NODE
  pjob=6

  rc=dim(matbatch)

  jobID=1
  for (i in 1:rc[1])
  {
    for (k in 1:normChain)
    {
      fname=paste("norm/norm-job",jobID,".sh",sep="",collapse="")
      sink(fname)
      cat("cd norm/results\n")
      idx = 1
      for (j in 1:rc[2])
      {
        cat(sprintf("exec Rscript ../../norm.R %d %d %d > ../out/out_std_%d_%d.out &\n",matbatch[i,j],k,normNChain,matbatch[i,j],k))
        if (idx == pjob)
        {
          cat("wait\n")
          idx = 1
        }
        else
        {
          idx = idx +1
        }
      }
      cat("wait\n")
      sink()
      jobID=jobID+1
    }
  }
  
# for (i in 0:nbatch-1)
# {
#   fname=paste("norm/norm-job",i,".sh",sep="",collapse="")
#   sink(fname)
#   cat("cd norm/results\n")
#   for(j in 1:normChain){
#     cat(sprintf("exec Rscript ../../norm.R %d %d %d\n",i,j,normNChain))
#   }
#   sink()
# }
  
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
  cat(sprintf("#SBATCH --array=1-%d\n",nbatch))
  cat("srun sh norm/norm-job$SLURM_ARRAY_TASK_ID.sh\n")
  
  # Stop writing to the file
  sink()
  
  system("chmod u+x bayesprot-norm.sh")
}
