load ("DegBoundedIdentification.m2")
load ("experiments/utils.m2")

----------------------------------------
-- Randomly generated graphs on 10 nodes
----------------------------------------
load("experiments/graphs/random-n10-0.2.m2")
n = 10;
degBd = 5;
maxTime = 10;
diago = true;
filePath = "experiments/results/random-n10-0.2-10s.m2";

results = compareIdentificationMethods(adjMatrices, n, degBd, maxTime, diago);
M = saveResults(results, filePath);
M

end
