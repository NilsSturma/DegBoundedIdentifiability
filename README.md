# Degree-Bounded Identification via Adaptive Orderings

This repository contains an implementation of the degree-bounded identfication algroithm developed by Pratik Misra, Benjamin Hollering and Nils Sturma. It can be used to verify rational identifiability of causal effects in linear structural equation models given by mixed graphs. The implementation is in Macaulay2, and the main function is `DegBoundedIdentification(G, degBd, maxTime)`, where `G` is a mixed graph of interest, `degBd` is the maximal degree and `maxTime` is the time limit for the algorithm. An example of how to run the algorithm is given by the following code snippet:
```
load ("DegBoundedIdentification.m2")
M = {
    {
        {0, 1, 0, 1}, 
        {0, 0, 0, 0}, 
        {0, 0, 0, 1}, 
        {0, 0, 0, 0}
    }, 
    {
        {0, 0, 0, 0}, 
        {0, 0, 1, 0}, 
        {0, 1, 0, 1}, 
        {0, 0, 1, 0}
    }
}
n=4
(A1, A2) = toSequence(M / matrix);
D = digraph(toList(1..n), A1);
B = bigraph(toList(1..n), A2);
G = mixedGraph(D, B);
degBd = 2;
maxTime = 10;
DegBoundedIdentification(G, degBd, maxTime)
```

To reproduce the experimental results in the paper, run:
```
bash experiments/experiment1/run.sh
bash experiments/experiment2/run.sh
```