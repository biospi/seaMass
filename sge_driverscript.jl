#sge_driverscript.jl

#Word of warning: Qsub will send mail for EVERY task in an array job
#Hence why I've not done so for the norm, model and plots jobs
email = "A.M.Phillips@liverpool.ac.uk"

#########################################
# Import
#########################################
print("Import...")
try symlink("../../input/data.Rdata", "import/results/data.Rdata") end
try symlink("../../input/parameters.Rdata", "import/results/parameters.Rdata") end
try symlink("../../input/design.Rdata", "import/results/design.Rdata") end

try mkdir("import/error") end
try mkdir("import/out") end

open("bayesprot-import.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"#\$ -o import/out\n")
  write(f,"#\$ -e import/error\n")
  write(f,"#\$ -M \n")
  write(f,"#\$ -m bes\n")
  write(f,"#\$ -l h_vmem=8G,h_rt=01:00:00\n")
  write(f,"cd import/results\n")
  write(f,"../../../bin/Rscript ../../import.R HPC")
end

run(`qsub bayesprot-import.sh`)
print(" Done.\n")



#########################################
# Norm
#########################################
print("Norm...")
norm_nchain = 10
njobs = 1000

try mkdir("norm/error") end
try mkdir("norm/out") end

open("bayesprot-norm.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"#\$ -o norm/out\n")
  write(f,"#\$ -e norm/error\n")
  write(f,"#\$ -l h_vmem=8G,h_rt=01:00:00\n")
  write(f,"sh norm-job\$SGE_TASK_ID.sh")
end

try
  symlink("../../input/parameters.Rdata", "norm/results/parameters.Rdata")
end

datafiles = filter(x -> ismatch(r"^[0-9]+\.Rdata$", x), readdir("import/results"))

for x in datafiles
  try
    symlink(string("../../import/results/", x), string("norm/results/", x))
  end
end

ids = sort(map(x -> parse(Int, match(r"^([0-9]+)\.Rdata$", x)[1]), datafiles))
all_tasks = vcat(map(y -> map(x -> string(y, ":", x, "/", norm_nchain), 1:norm_nchain), ids)...)

jobs = Array{Array{AbstractString,1},1}(njobs)

for i in 1:njobs
  jobs[i] = all_tasks[collect(i:njobs:length(all_tasks))]
end


for i in 1:njobs
  open("norm-job$i.sh","w") do f
    write(f,"cd norm/results\n")
    write(f,"exec ../../../bin/Rscript ../../norm.R HPC 0")
    for t in jobs[i]
      write(f," $t")
    end
  end
end

run(`qsub -t 1-$njobs bayesprot-norm.sh`)
print("Done\n")

#########################################
# Exposures
#########################################
print("Exposures...")
try symlink("../../import/results/index.Rdata", "exposures/results/index.Rdata") end
try symlink("../../input/parameters.Rdata", "exposures/results/parameters.Rdata") end
try symlink("../../input/design.Rdata", "exposures/results/design.Rdata") end
try symlink("../../norm/results/", "exposures/results/results") end

try mkdir("exposures/error") end
try mkdir("exposures/out") end

open("bayesprot-exposures.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"#\$ -o exposures/out\n")
  write(f,"#\$ -e exposures/error\n")
  write(f,"#\$ -M "*email*"\n")
  write(f,"#\$ -m bes\n")
  write(f,"#\$ -l h_vmem=8G,h_rt=01:00:00\n")
  write(f,"cd exposures/results\n")
  write(f,"../../../bin/Rscript ../../exposures.R HPC")
end

run(`qsub bayesprot-exposures.sh`)
print(" Done.\n")


#########################################
# Model
#########################################
print("Model...")

model_nchain = 100
njobs = 1000

try symlink("../../input/parameters.Rdata", "model/results/parameters.Rdata") end
try symlink("../../exposures/results/exposures.Rdata", "model/results/exposures.Rdata") end
try symlink("../../input/design.Rdata", "model/results/design.Rdata") end

try mkdir("model/error") end
try mkdir("model/out") end

open("bayesprot-model.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"#\$ -o model/out\n")
  write(f,"#\$ -e model/error\n")
  write(f,"#\$ -l h_vmem=8G,h_rt=01:00:00\n")
  write(f,"sh model-job\$SGE_TASK_ID.sh")
end

datafiles = filter(x -> ismatch(r"^[0-9]+\.Rdata$", x), readdir("import/results"))

for x in datafiles
  try
    symlink(string("../../import/results/", x), string("model/results/", x))
  end
end

ids = sort(map(x -> parse(Int, match(r"^([0-9]+)\.Rdata$", x)[1]), datafiles))

all_tasks = vcat(map(y -> map(x -> string(y, ":", x, "/", model_nchain), 1:model_nchain), ids)...)

jobs = Array{Array{AbstractString,1},1}(njobs)

for i in 1:njobs
  jobs[i] = all_tasks[collect(i:njobs:length(all_tasks))]
end

for i in 1:njobs
  open("model-job$i.sh","w") do f
    write(f,"cd model/results\n")
    write(f,"exec ../../../bin/Rscript ../../model.R HPC 0")
    for t in jobs[i]
      write(f," $t")
    end
  end
end

run(`qsub -t 1-$njobs bayesprot-model.sh`)
print("Done\n")

#########################################
# Plots
#########################################
println("Plots...")

try symlink("../../input/parameters.Rdata", "plots/results/parameters.Rdata") end
try symlink("../../input/design.Rdata", "plots/results/design.Rdata") end

try mkdir("plots/error") end
try mkdir("plots/out") end

open("bayesprot-plots.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"#\$ -o plots/out\n")
  write(f,"#\$ -e plots/error\n")
  write(f,"#\$ -l h_vmem=8G,h_rt=01:00:00\n")
  write(f,"sh plots-job\$SGE_TASK_ID.sh")
end

datafiles = filter(x -> ismatch(r"^[0-9]+\.Rdata$", x), readdir("import/results"))

for x in datafiles
  try
    symlink(
      string("../../model/results/", parse(Int, match(r"^([0-9]+)\.Rdata$", x)[1]),"/"),
      string("plots/results/", parse(Int, match(r"^([0-9]+)\.Rdata$", x)[1]))
      )
  end
end

ids = sort(map(x -> parse(Int, match(r"^([0-9]+)\.Rdata$", x)[1]), datafiles))

all_tasks = vcat(map(y -> string(y), ids)...)

jobs = Array{Array{AbstractString,1},1}(njobs)

for i in 1:njobs
  jobs[i] = all_tasks[collect(i:njobs:length(all_tasks))]
end

for i in 1:njobs
  open("plots-job$i.sh","w") do f
    write(f,"cd plots/results\n")
    write(f,"exec ../../../bin/Rscript plots.R HPC 0")
    for t in jobs[i]
      write(f," $t")
    end
  end
end

run(`qsub -t 1-$njobs bayesprot-plots.sh`)
print("Done\n")

#########################################
# Output
#########################################
println("Output...")

try symlink("../../import/results/index.Rdata", "output/results/index.Rdata") end
try symlink("../../input/parameters.Rdata", "output/results/parameters.Rdata") end
try symlink("../../input/design.Rdata", "output/results/design.Rdata") end
try symlink("../../plots/results/stats/", "output/results/stats") end

try mkdir("output/error") end
try mkdir("output/out") end

open("bayesprot-output.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"#\$ -o output/out\n")
  write(f,"#\$ -e output/error\n")
  write(f,"#\$ -M "*email*"\n")
  write(f,"#\$ -m bes\n")
  write(f,"#\$ -l h_vmem=8G,h_rt=01:00:00\n")
  write(f,"cd output/results\n")
  write(f,"../../../bin/Rscript ../../output.R HPC")
end

run(`qsub bayesprot-output.sh`)
print(" Done.\n")
