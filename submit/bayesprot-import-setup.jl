print("Import...\n")
try symlink("../../input/data.Rdata", "import/results/data.Rdata") end
try symlink("../../input/parameters.Rdata", "import/results/parameters.Rdata") end
try symlink("../../input/design.Rdata", "import/results/design.Rdata") end

try mkdir("import/error") end
try mkdir("import/out") end
