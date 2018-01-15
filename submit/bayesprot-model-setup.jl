model_nchain = parse(ARGS[1])
njobs = parse(ARGS[2])

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

jobs = [all_tasks[collect(i:njobs:length(all_tasks))] for i in 1:njobs]

for i in 1:njobs
  open("model/model-job$i.sh","w") do f
    write(f,"cd model/results\n")
    write(f,"module add languages/R-3.4.1-ATLAS\n")
    write(f,"exec Rscript ../../model.R HPC 0")
    for t in jobs[i]
      write(f," $t")
    end
  end
end
