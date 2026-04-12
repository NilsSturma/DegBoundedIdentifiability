library(rjson)
library(igraph)
source("utils.R")

# Only those with number of edges \leq choose(n, 2)
getAllMixedGraphs <- function(n){
  nLowerDiag <- n * (n-1) / 2
  nOffDiag <- n * (n-1)
  lowerMask <- lower.tri(matrix(0,n,n))
  upperMask <- upper.tri(matrix(0,n,n))

  
  results <- list()

  
  for (m in 0:nLowerDiag){ # since we restrict maximal number of edges
    print(paste("m=", m, sep=""))
    edgeIndicesDir <- t(combn(nOffDiag,m))
    for (i in 1:nrow(edgeIndicesDir)){
      D <- matrix(0,n,n)
      D[upperMask | lowerMask][edgeIndicesDir[i,]] <- 1
      if (is_dag(graph_from_adjacency_matrix(D, mode="directed"))){
        for(k in 0:nLowerDiag){
          if ( (k+m) > nLowerDiag ){
            break # since we restrict maximal number of edges
          }
          edgeIndicesBiDir <- t(combn(nLowerDiag,k))
          for (j in 1:nrow(edgeIndicesBiDir)){
            B <- matrix(0,n,n)
            B[lowerMask][edgeIndicesBiDir[j,]] <- 1
            B <- B + t(B)
            g = list("D"=D, "B"=B)
            #print(g)
              
            if (length(results)==0){
              results = c(results, list(g))
            } else {
              is_new = TRUE
              for (l in 1:length(results)){
                if (is_isomorphic_mixed(results[[l]], g)){
                  is_new = FALSE
                  break
                }
              }
              if (is_new){
                results = c(results, list(g))
              }
            }
          }
        }
      }
    }
  }
  return(results)
}

#######################################
### Acyclic mixed graphs on 4 nodes ###
#######################################
# With at most 6 edges
res <- getAllMixedGraphs(4)
length(res)

# Save to json
jsonData = toJSON(res)  # saved column-wise
write(jsonData, "graphs/acyclic-mixed-graphs-4-nodes.json")

# Save to Macaulay2
graphs_str = getGraphsString(res)
writeLines(graphs_str, "graphs/acyclic-mixed-graphs-4-nodes.m2")
