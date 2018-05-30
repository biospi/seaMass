njobs = parse(ARGS[1])
nCpu = parse(ARGS[2])

N = njobs * nCpu

try symlink("../../input/parameters.Rdata", "plots/results/parameters.Rdata") end
try symlink("../../input/design.Rdata", "plots/results/design.Rdata") end

try mkdir("plots/error") end
try mkdir("plots/out") end

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

#jobs = [all_tasks[collect(i:njobs:length(all_tasks))] for i in 1:njobs]
jobs = [all_tasks[collect(i:N:length(all_tasks))] for i in 1:N]

for i in 1:njobs
  idx=collect(i:njobs:N)
  open("plots/plots-job$i.sh","w") do f
    write(f,"cd plots/results\n")
    for j in idx
      write(f,"exec Rscript ../../plots.R HPC 0")
      for t in jobs[i]
        write(f," $t")
      end
      write(f," &\n")
    end
    write(f,"wait\n")
  end
end
