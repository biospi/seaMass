print("Norm...")
norm_nchain = parse(ARGS[1])
njobs = parse(ARGS[2])

try mkdir("norm/error") end
try mkdir("norm/out") end

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

jobs = [all_tasks[collect(i:njobs:length(all_tasks))] for i in 1:njobs]

for i in 1:njobs
  open("norm/norm-job$i.sh","w") do f
    write(f,"cd norm/results\n")
    write(f,"module add languages/R-3.4.1-ATLAS\n")
    write(f,"exec Rscript ../../norm.R HPC 0")
    for t in jobs[i]
      write(f," $t")
    end
  end
end
