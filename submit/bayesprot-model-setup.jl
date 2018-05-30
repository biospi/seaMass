model_nchain = parse(ARGS[1])
njobs = parse(ARGS[2])
nCpu = parse(ARGS[3])

N = njobs * nCpu

try symlink("../../input/parameters.Rdata", "model/results/parameters.Rdata") end
try symlink("../../exposures/results/exposures.Rdata", "model/results/exposures.Rdata") end
try symlink("../../input/design.Rdata", "model/results/design.Rdata") end

try mkdir("model/error") end
try mkdir("model/out") end

datafiles = filter(x -> ismatch(r"^[0-9]+\.Rdata$", x), readdir("import/results"))
for x in datafiles
  try
    symlink(string("../../import/results/", x), string("model/results/", x))
  end
end

ids = sort(map(x -> parse(Int, match(r"^([0-9]+)\.Rdata$", x)[1]), datafiles))

all_tasks = vcat(map(y -> map(x -> string(y, ":", x, "/", model_nchain), 1:model_nchain), ids)...)

#jobs = [all_tasks[collect(i:njobs:length(all_tasks))] for i in 1:njobs]
jobs = [all_tasks[collect(i:N:length(all_tasks))] for i in 1:N]

for i in 1:njobs
  idx=collect(i:njobs:N)
  open("model/model-job$i.sh","w") do f
    write(f,"cd model/results\n")
    for j in idx
      write(f,"exec Rscript ../../model.R HPC 0")
      for t in jobs[j]
        write(f," $t")
      end
      write(f," &\n")
    end
    write(f,"wait\n")
  end
end
