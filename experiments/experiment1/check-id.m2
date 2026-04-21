load ("DegBoundedIdentification.m2")
load ("experiments/utils.m2")

--------------------------
-- Mixed Graphs on 4 nodes 
--------------------------
load("experiments/graphs/acyclic-mixed-graphs-4-nodes.m2")
n = 4;
degBd = 5;
maxTime = 60;
diago = false;
tian = true;
filePath = "experiments/results/4-nodes.m2";

results = compareIdentificationMethods(adjMatrices, n, degBd, maxTime, diago, tian);
M = saveResults(results, filePath);
M

end
