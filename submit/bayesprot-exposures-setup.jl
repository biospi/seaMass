try symlink("../../import/results/index.Rdata", "exposures/results/index.Rdata") end
try symlink("../../input/parameters.Rdata", "exposures/results/parameters.Rdata") end
try symlink("../../input/design.Rdata", "exposures/results/design.Rdata") end
try symlink("../../norm/results/", "exposures/results/results") end

try mkdir("exposures/error") end
try mkdir("exposures/out") end
