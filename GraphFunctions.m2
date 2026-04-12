needsPackage "StatGraphs"
needsPackage "Probability"

tianDecomposition = method()
tianDecomposition (MixedGraph) := List => (G) -> (    

    Gb = underlyingGraph bigraph G;
    Gb = addVertices(Gb, vertices G);
    Gd = digraph G;
    Gd = addVertices(Gd, vertices G);
    connComp = connectedComponents Gb;
    GGb = bigraph G;
                                

    
    for comp in connComp list (
        gd = inducedSubgraph(Gd, comp);
        ggb = inducedSubgraph(GGb, comp);
        for v in comp list (
            gd = addVertices(gd, toList parents(Gd, v));
        findedg = findPaths(digraphTranspose Gd, v, 1);
        findedg = apply(findedg, e -> reverse e);
        if #findedg > 0 then gd = addEdges'(gd, findedg);
        ggb = addVertices(ggb, toList parents(Gd, v));
        );
        mixedGraph(gd, bigraph edges ggb)
        )
)



treks = method()
treks (Digraph, Thing, Thing) := (G, i, j) -> (
    parents1 = forefathers(G, i);
    parents2 = forefathers(G, j);
    commonParents = intersect(parents1, parents2);
    pathset1 = pathsEndingInSet(G, i, commonParents);
    pathset2 = pathsEndingInSet(G, j, commonParents);
    -- trivial cases
    if #pathset1 == 0 then return pathset2;		   
    if #pathset2 == 0 then return pathset1;
    -- nontrivial case
    if isMember(i, commonParents) then pathset1##pathset1 = {i};
    if isMember(j, commonParents) then pathset2##pathset2 = {j};
    listingParents = toList(commonParents);
    hashingPaths = new MutableHashTable from apply(#listingParents, i -> (listingParents#(i), new MutableList));
    for p in pathset1 do (
	(hashingPaths#(p#-1))##(hashingPaths#(p#-1)) = p;
	);
    listOfResults = for p in pathset2 list (
	(table(toList(hashingPaths#(p#-1)), {p}, (k,l)-> {k,l}))#0#0
	);
    return listOfResults		   
)

pathsEndingInSet = method()
pathsEndingInSet(Digraph, Thing, Set) := (G, i, commonParents) -> (
    listOfPaths = new MutableList;
    queue = new MutableList from {{i}};
    while #queue != 0 do(				    -- kind of bfs for all paths CP->i 
	p = remove(queue, 0);
	extensions = parents(G, p#-1);
	if #extensions != 0 then (
	    for endpoint in toList(extensions) do(
		extPath = new MutableList from p;
		extPath##extPath = endpoint;
		queue##queue = toList(extPath);
		if isMember(endpoint, commonParents) == true then (
		    listOfPaths##listOfPaths = queue#-1;
		    );
		);
	    );
	);
    return listOfPaths;
)


-- maxium in degree in the underyling directed graph of a mixed graph G
maxDiInDegree = G -> (

    D := digraph(G);

    vertices(D) / (i -> degreeIn(D, i)) // max
)


-- typical Erdos-Renyi random undirected graph
erdosRenyiGraph = (n, p) -> (

    P := bernoulliDistribution(p);
    E := for e in subsets(1..n, 2) list if random(P) == 1 then e else continue;

    return graph(toList(1..n), E)
)


-- same Erdos-Renyi model but resulting in a bidirected graph
erdosRenyiBiGraph = (n, p) -> (

    P := bernoulliDistribution(p);
    E := for e in subsets(1..n, 2) list if random(P) == 1 then e else continue;

    return bigraph(toList(1..n), E)
)


-- random digraph obtained by sampling the parents of j with a binomial distribution with probability p*sqrt(n-1/j-1)
-- in this case the expected number of parents for each node is sqrt(j-1)*sqrt(n-1)*p
-- this agrees with the Erdos-Renyi model for the last node in the total order
parentalERDiGraph = {RejectionThreshold => 2} >> opts ->  (n, p) -> (

    E := flatten for j from 2 to n list(

        P := bernoulliDistribution(p*sqrt((n-1)/(j-1)));

        for i from 1 to j-1 list if random(P) == 1 then {i, j} else continue
    );
    
    D := digraph(toList(1..n), E);

    degs := vertices(D) / (i -> degreeIn(D, i));
    k := opts.RejectionThreshold;
    degBounds := {{0,0}}|apply(parentalERExp(n, p), parentalERStdDev(n, p), (i, j) -> {max(0, i-k*j), i+k*j});

    while not isWithinDegreeBounds(degBounds, D) do(

        E = flatten for j from 2 to n list(

            P := bernoulliDistribution(p*sqrt((n-1)/(j-1)));

            for i from 1 to j-1 list if random(P) == 1 then {i, j} else continue

        );
    
        D = digraph(toList(1..n), E);
    );

    return D;
);


-- checks if the in degree of every node in a directed graph is within given bounds
isWithinDegreeBounds = (degBounds, D) -> all(pairs(vertices(D)), (i, j) -> degreeIn(D, j) >= (degBounds_i)_0 and degreeIn(D, j) <= (degBounds_i)_1)


-- Repeatedly samples Erdos-Renyi random graphs until a connected graph is found
connectedErdosRenyiBiGraph = (n, p) -> (

    if binomial(n, 2)*p < n-1 then print("warning: the expected number of edges is less than the required number to connect the graph");

    P := bernoulliDistribution(p);
    E := for e in subsets(1..n, 2) list if random(P) == 1 then e else continue;
    G := graph(toList(1..n), E);

    while not isConnected(G) do(

        E = for e in subsets(1..n, 2) list if random(P) == 1 then e else continue;
        G = graph(toList(1..n), E);
    );

    return bigraph(toList(1..n), E);
)

-- random 
parentalERExp = (n, p) -> apply(toList(2..n), j -> p*sqrt(j-1)*sqrt(n-1))
parentalERStdDev = (n, p) -> apply(toList(2..n), j -> sqrt((j-1)*p*sqrt((n-1)/(j-1))*(1-p*sqrt((n-1)/(j-1)))))

-- independently samples a digraph and bigraph with edge probabilities p1 and p2 to create a mixed graph
erdosRenyiMixedGraph = (n, p1, p2) -> mixedGraph(parentalERDiGraph(n, p1), erdosRenyiBiGraph(n, p))


-- samples a mixed graph whose bidirected part is connected, thus its tian decomposition has a single component
stronglyConnectedMixedGraph = (n, p1, p2) -> mixedGraph(parentalERDiGraph(n, p1), connectedErdosRenyiBiGraph(n, p2))