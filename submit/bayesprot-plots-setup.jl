njobs = parse(ARGS[1])

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

jobs = [all_tasks[collect(i:njobs:length(all_tasks))] for i in 1:njobs]

for i in 1:njobs
  open("plots/plots-job$i.sh","w") do f
    write(f,"cd plots/results\n")
    write(f,"module add languages/R-3.4.1-ATLAS\n")
    write(f,"exec Rscript ../../plots.R HPC 0")
    for t in jobs[i]
      write(f," $t")
    end
  end
end
