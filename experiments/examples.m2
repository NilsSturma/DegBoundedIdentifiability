load ("DegBoundedIdentification.m2")
load ("experiments/utils.m2")

---------------------------
-- Example for Introduction
---------------------------

M = {
    {
        {0, 1, 0, 1}, 
        {0, 0, 1, 0}, 
        {0, 0, 0, 1}, 
        {0, 0, 0, 0}
    }, 
    {
        {0, 0, 0, 0}, 
        {0, 0, 0, 0}, 
        {0, 0, 0, 1}, 
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
resGP = GarciaPuenteIdentification(G, maxTime)
pairs resGP#1

---------------------------
-- Example for Section 2
---------------------------

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
resGP = GarciaPuenteIdentification(G, maxTime)
pairs resGP#1
