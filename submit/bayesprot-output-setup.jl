try symlink("../../import/results/index.Rdata", "output/results/index.Rdata") end
try symlink("../../input/parameters.Rdata", "output/results/parameters.Rdata") end
try symlink("../../input/design.Rdata", "output/results/design.Rdata") end
try symlink("../../plots/results/stats/", "output/results/stats") end
try symlink("../../plots/results/samplestats/", "output/results/samplestats") end

try mkdir("output/error") end
try mkdir("output/out") end
