needsPackage "GraphicalModels"
needsPackage "Probability"
load("GraphFunctions.m2")

KK = ZZ/32003


-- check if a multigrading refines total degree
hasStdGrading = A -> rank(matrix({toList(numcols(A):1)}) || A) == rank(A)


-- checks if the jacobian is full rank probabilistically
isLocallyIdentifiable = {Coefficients => ZZ/32003, Trials => 1000} >> opts -> R -> (

    L := gaussianParametrization(R);
    S := KK[take(gens R, numParams R)];
    phi := sub(matrix {flatten for i from 0 to numcols(L) - 1 list for j from i to numcols(L) - 1 list L_(i, j)}, S);
    Jac := jacobian(phi);
    KK := opts.Coefficients;
    
    for t from 1 to opts.Trials do if rank(sub(Jac, apply(gens S, i -> i => random(KK)))) == numrows(Jac) then return true;
    
    return false;
)

isIdentifiable = {Coefficients => ZZ/32003, Trials => 10} >> opts -> R -> (

    L := gaussianParametrization(R);
    KK := opts.Coefficients;
    J := semElimIdeal(R);

    degs := for t from 1 to opts.Trials do(

        paramVals := apply(take(gens(R), numParams(R)), i -> i => random(KK));
        L = sub(L, paramVals);
        sigmaVals := flatten for i from 1 to n list for j from i to n list s_(i, j) => L_(i-1, j-1);
        I := sub(J, sigmaVals);
        curDeg := degree I;
        if curDeg != 1 then return false;
    );

    return true
)


-- find a polynomial supported on {t}|oldVars F which is linear in t
findLinearPolynomial = (t, oldVars, F) -> (

    tPos := position(gens ring t, i -> i == t);

    for g in F list(

        supp := support(g);
        if not member(t, supp) then continue; 
        if not all(supp, i -> member(i, oldVars|{t})) then continue;
        
        if all(exponents(g), i -> i_tPos < 2) then return g;
    );

    return 1_ZZ
)


-- check if the set F contains a polynomial supported on {t}|oldVars which is linear in t
containsLinearPolynomial = (t, oldVars, F) -> (

    tPos := position(gens ring t, i -> i == t);

    for g in F list(

        supp := support(g);
        if not member(t, supp) then continue; 
        if not all(supp, i -> member(i, oldVars|{t})) then continue;
        
        if all(exponents(g), i -> i_tPos < 2) then return true;
    );

    return false
)


-- check if the set F contains a polynomial supported on {t}|oldVars 
containsNonLinearPolynomial = (t, oldVars, F) -> (

    if #F == 0 then return false;

    tPos := index t;

    for g in F list(

        supp := support(g);
        if not member(t, supp) then continue; 
        if not all(supp, i -> member(i, oldVars|{t})) then continue;
        
        if any(exponents(g), i -> i_tPos > 1) then return true;
    );

    return false
)


-- makes the elimination ideal of a SEM on a mixed graph
semElimIdeal = R -> ideal(unique flatten entries(covarianceMatrix(R) - gaussianParametrization(R)));


-- computes the total number of model parameters
numParams = R -> #edges(digraph(R.graph)) +  #vertices(digraph(R.graph)) + #edges(bigraph(R.graph));


-- grade the ring R by the trek lengths
trekLengthGrading = R -> (

    L := gaussianParametrization(R);

    n := #vertices(digraph(R.graph));
    
    sigmaWeights := flatten for i from 0 to n-1 list for j from i to n-1 list if L_(i, j) == 0_R then 0 else (degree L_(i, j))_0;

    return hashTable apply(gens R, toList(numParams(R):1)|sigmaWeights, (i,j) -> i => j);
)


-- homogenize the elimination ideal of the SEM
homSemElimIdeal = (w, R) -> (

    newOrder := support(bidirectedEdgesMatrix(R))|support(directedEdgesMatrix(R))|support(covarianceMatrix(R));
    S := KK[newOrder|{h}, Degrees => apply(newOrder, i -> w#i)|{1}, MonomialOrder => Eliminate numParams(R)];
    n := #vertices(R.graph);

    Sigma := sub(covarianceMatrix(R), S);
    Im := sub(gaussianParametrization(R), S);
    M := Sigma-Im;
    homGens := flatten for i from 0 to n-1 list for j from i to n-1 list homogenize(M_(i, j), h);

    return ideal homGens
)


-- elimination but with degree bounding
elimIdeal = (w, elimVars, I) -> (

    R := ring I;
    remainVars := delete(null, apply(gens(R), i -> if not member(i, elimVars) then i));
    newGens := elimVars|remainVars;
    newWeight := apply(newGens, i -> w#i);
    S := KK[newGens, Degrees => newWeight, MonomialOrder => Eliminate(#elimVars)];
    J := sub(I, S);

    return J
)





-- This is an implementation of the main algorithm (Algorithm 2) of our paper
-- d is the degree limit and R is a ring made from gaussianRing
basicDegBoundedId = method(Options => {Verbose => true, Coefficients => QQ, ReturnFullGB => false});
basicDegBoundedId (Number, MixedGraph) := HashTable => opts -> (degBd, G) -> (

    -- graph the graph from the ring and mak the graph ideal
    -- then create associated sets of variables
    KK := opts.Coefficients;
    R := gaussianRing(G, Coefficients => KK);
    V := sort vertices(digraph(G));
    n := #V;
    

    --if opts.Verbose then print("SEM Graph Ideal Computed");
    -- make a hashtable so we can easily get the grading back as we re-order the variables
    w := trekLengthGrading R;
    J := homSemElimIdeal(w, R);
    use R;

    -- create lists of our parameters and covariance variables
    wParams := apply(V, i -> p_(i, i))|(for e in edges(bigraph(G)) list p_(toSequence sort toList e));
    lParams := for e in edges(digraph(G)) list l_(e_0, e_1);
    sigVars := flatten for i from 0 to n-1 list for j from i to n-1 list if w#(s_(V_i, V_j)) != 0 then s_(V_i, V_j) else continue;

    -- create hash tables to store identifying polynomials and gbs as we go
    gbHash := new MutableHashTable;
    idHash := new MutableHashTable;
    idParams := {};
    unidParams := wParams|lParams;

    -- create a ring with an elimination ordering and put the graph ideal in the ring
    S := KK[unidParams|sigVars|{h}, Degrees => apply(unidParams|sigVars, t -> w#t)|{1}, MonomialOrder => Eliminate (#unidParams)];
    I := sub(J, S);

    -- make the weighted deg bound in the trek weighting
    weightedDegBd := max(values w)*degBd;
    --print({degBd, weightedDegBd});

    while #unidParams > 0 do(
        for d from 2 to weightedDegBd do(

            identifiedNewParam := false;

            if opts.Verbose then print(concatenate("computing GB in degree: ", toString d));

            for q in unidParams do(
                
                if opts.Verbose then print(concatenate("Identifying ", toString(q), " in degree: ", toString d));

                -- reorder variables in ring
                rem := select(unidParams, i -> i != q);
                newVarOrder := rem|{q}|toList(idParams)|sigVars;
                S = KK[newVarOrder|{h}, Degrees => apply(newVarOrder, t -> w#t)|{1}, MonomialOrder => Eliminate (#unidParams)];
                I = sub(I, S);
                
                -- compute the GB for the current ordering
                curGB := gb(I, DegreeLimit => d);
                F := flatten entries sub(sub(gens curGB, h => 1), R);

                -- find the new parameters which we identified in this degree and them to our list of identified params
                newIdParams = for t in unidParams list(

                    f := findLinearPolynomial(sub(t, R), toList(idParams)|sigVars|{h}, F);

                    if f == 1_ZZ then continue else(

                        idHash#t = f;
                        sub(t, R)
                    )
                );

                idParams = newIdParams|idParams;

                -- if we identified nothing new in the max degree then we stop
                if #newIdParams == 0 and d == weightedDegBd then return idHash;
                
                -- otherwise we find the remaining parameters
                unidParams = select(unidParams, t -> not member(t, newIdParams));

                -- if we have identified everything then we're done
                if all(lParams, t -> member(t, idParams)) then(

                    if opts.Verbose then print(concatenate("identified: ", toString(newIdParams)));
                    return idHash
                );


                if #newIdParams > 0 then(
                    if opts.Verbose then print(concatenate("identified: ", toString(newIdParams)));
                    if opts.Verbose then print("recomputing GB in lower degree with new ordering");
                    identifiedNewParam = true;   
                    break;
                );
            );

            if identifiedNewParam then break;     
        );
    );
);


-- performs tian decomposition and then runs our main algorithm on each mixed component
tianDegBoundedId = method(Options => {Verbose => true, Coefficients => QQ})
tianDegBoundedId (Number, Ring) := HashTable => opts -> (degBd, R) -> (

    G := R.graph;
    tianComps := tianDecomposition(G);
    if opts.Verbose then print(concatenate("Number of mixed components in graph: ", toString(#tianComps)));
    compCounter := 1;

    idTables := for H in tianComps list(

        if opts.Verbose then print(concatenate("Identifying parameters in component: ", toString compCounter));
        compCounter = compCounter + 1;
        basicDegBoundedId(degBd, H, Verbose => opts.Verbose, Coefficients => opts.Coefficients)
    );

    -- merge the hashtables into one
    idPairs := flatten for H in idTables list apply(keys H, i -> {sub(i, R), sub(H#i, R)});
    lParams := for e in edges(digraph(G)) list sub(l_(e_0, e_1), R);
    idHash := hashTable select(idPairs, i -> member(i_0, lParams));
    finalHash := hashTable idPairs;
    if all(lParams, t -> member(t, keys idHash)) then return (true, finalHash) else return (false, finalHash)
)


tianDegBoundedId (Number, MixedGraph) := HashTable => opts -> (degBd, G) -> (

    KK := opts.Coefficients;
    R := gaussianRing(G, Coefficients => KK);
    return tianDegBoundedId(degBd, R, Verbose => opts.Verbose, Coefficients => KK);
)



-------------------
-- Wrapper function
-------------------


DegBoundedIdentification = (G, degBd, maxTime, tian) -> (
    R := gaussianRing(G, Coefficients => KK);
    if tian then(
        try (
            alarm maxTime; 
            result = tianDegBoundedId(degBd, G, Verbose=>false);
            return(result);
        ) else 
            return (false, new HashTable);
    ) else
        try (
            alarm maxTime; 
            result = basicDegBoundedId(degBd, G, Verbose=>false);
            idPairs := apply(keys result, i -> {sub(i, R), sub(result#i, R)});
            lParams := for e in edges(digraph(G)) list sub(l_(e_0, e_1), R);
            idHash := hashTable select(idPairs, i -> member(i_0, lParams));
            finalHash := hashTable idPairs;
            if all(lParams, t -> member(t, keys idHash)) then return (true, finalHash) else return (false, finalHash)
        ) else 
            return (false, new HashTable);
);

end
