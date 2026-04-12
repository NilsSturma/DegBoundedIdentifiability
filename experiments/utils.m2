GarciaPuenteIdentification = (G, maxTime) -> (

    R := gaussianRing(G, Coefficients => KK);

    try (
        alarm maxTime; 
        H := identifyParameters(R);
    ) else 
        return (false, new HashTable);
    
    lambdas := for e in edges(digraph(G)) list sub(l_(e_0, e_1), R);
    V := sort vertices(digraph(G));
    n := #V;
    w := trekLengthGrading R;
    sigVars := flatten for i from 0 to n-1 list for j from i to n-1 list sub(s_(V_i, V_j),R);

    idHash := new MutableHashTable;

    --print H;

    idParams := for t in lambdas list(

        I := H#t;
        generators := flatten entries gens I;
        f := findLinearPolynomial(t, sigVars, generators);

        if f == 1_ZZ then continue else(

            idHash#t = f;
            sub(t, R)
        )
    );  
    if all(lambdas, t -> member(t, keys idHash)) then return (true, idHash) else return (false, idHash)
);

maxDeg = (idHash, R) -> (
    maxdeg := 0;
    for k in keys idHash do(
        f := sub(idHash#k, R);
        maxdeg = max(maxdeg, first(degree(f)));
    );
    return maxdeg;
);

compareIdentificationMethods = (adjMatrices, n, degBd, maxTime, diago) -> (
    results = {};
    
    for i from 0 to #adjMatrices - 1 do(
        if mod(i,10) == 0 then print("Graph " | toString(i) | "\n");

        (A1, A2) = toSequence(adjMatrices#i / matrix);
        if diago then (
            A2 = A2 - map(ZZ^n);
        );
        D = digraph(toList(1..n), A1);
        B = bigraph(toList(1..n), A2);
        G = mixedGraph(D, B);
        R = gaussianRing(G, Coefficients => KK);

        -- Garcia-Puente
        GP = elapsedTiming(GarciaPuenteIdentification(G, maxTime));
        GPtime = GP#0;
        GPCompleteRes = GP#1;
        GPres = GPCompleteRes#0;
        GPdeg = maxDeg(GPCompleteRes#1, R);
        
        -- Degree-Bounded
        DegBd = elapsedTiming(DegBoundedIdentification(G, degBd, maxTime));
        DegBdtime = DegBd#0;
        DegBdCompleteRes = DegBd#1;
        DegBdres = DegBdCompleteRes#0;
        DegBddeg = maxDeg(DegBdCompleteRes#1, R);

        results = append(results, new HashTable from {
            "objectIndex" => i,
            "GP" => GPres, "timeGP" => GPtime, "GPdeg" => GPdeg,
            "DegBd" => DegBdres, "timeDegBd" => DegBdtime, "DegBddeg" => DegBddeg
        });
    );
    return results;
);

-- Save as matrix
saveResults = (results, filePath) -> (
    toNumber = x -> if x === true then sub(1, RR) else if x === false then sub(0, RR) else x;
    keysOrder = keys results#1;
    keysOrder = {keysOrder#3} | {keysOrder#2} | {keysOrder#5} | {keysOrder#6} | {keysOrder#1} | {keysOrder#4};  -- {GP, DegBd, GPdeg, DegBddeg, timeGP, timeDegBd}
    valuesList = apply(results, ht -> apply(keysOrder, k -> toNumber(ht#k)));
    Mres = matrix valuesList;
    filePath << toString(Mres) << endl << close;
    return(Mres);
);
