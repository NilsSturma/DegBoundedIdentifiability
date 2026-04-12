library(SEMID)
library(rjson)

################################
##### Function Definitions #####
################################

rAcyclicDirectedAdjMatrix <- function(n, p) {
  return(1 * (upper.tri(matrix(0, n, n)) & matrix(sample(c(T, F), n^2, replace = T,
                                                         prob = c(p, 1 - p)), ncol = n)))
}

rBidirectedAdjMat <- function(n, p) {
  mat <- matrix(runif(n^2), ncol = n) < p
  mat <- mat * upper.tri(mat)
  return(mat + t(mat) + diag(n))
}

generateGraphs <- function(ngraphs, nobs, pdir, pbidir){
  results = list()
  for (i in 1:ngraphs){
    set.seed(i)
    D = rAcyclicDirectedAdjMatrix(nobs, pdir)
    B = rBidirectedAdjMat(nobs, pbidir)
    results[[i]] = list("D"=D,"B"=B, "ndir"=sum(D), "nbidir"=(sum(B)-nobs)/2)
  }
  return(results)
}

getMatrixString <- function(M){
  n = ncol(M)
  for (i in 1:n){
    col_str = paste(M[,i],collapse=", ")
    col_str = paste("{", col_str, "}", sep="")
    if (i==1){
      M_str = col_str
    }
    else { 
      M_str = paste(M_str, col_str, sep=", ")
    }
  }
  M_str = paste("{", M_str, "}", sep="")
  return(M_str)
}

getSingleGraphString <- function(D,B){
  return(paste("{", getMatrixString(D), ", ", getMatrixString(B), "}", sep=""))
}

getGraphsString <- function(graphs){
  for (i in 1:length(graphs)){
    D = graphs[[i]]$D
    B = graphs[[i]]$B
    if (i==1){
      str = paste("adjMatrices = {", getSingleGraphString(D,B), ", \n", sep="")
    } 
    else if (i==length(graphs)){
      str = paste(str, getSingleGraphString(D,B), "}", sep="")
    }
    else {
      str = paste(str, getSingleGraphString(D,B), ", \n", sep="")
    }
  }
  return(str)
}

####################################
##### Randomly Generate Graphs #####
####################################

pdir = 0.2
pbidir = 0.2
nobs = 10
ngraphs = 1000

set.seed(101)
graphs = generateGraphs(ngraphs, nobs, pdir, pbidir)

# Check HTC-identifiability
for (i in 1:length(graphs)){
  if((i%%10) == 0){print(i)}
  Lvec = graphs[[i]]$D
  L <- matrix(Lvec, nrow=sqrt(length(Lvec)))
  Ovec = graphs[[i]]$B
  O <- matrix(Ovec, nrow=sqrt(length(Ovec))) - diag(nobs)
  g = MixedGraph(L, O)
  
  res <- htcID(g, tianDecompose = TRUE)
  graphs[[i]]$htcID <- (sum(sapply(res$unsolvedParents, length)) == 0)
}

# How many are identifiable by HTC?
countHTC = 0
for (i in 1:length(graphs)){
  if (graphs[[i]]$htcID==TRUE){
    countHTC = countHTC + 1
  }
}
print(countHTC) 

# Save to json
jsonData = toJSON(graphs)  # saved columnwise
write(jsonData, "graphs/random-n10-0.2.json")

# Save to Macaulay2
graphs_str = getGraphsString(graphs)
writeLines(graphs_str, "graphs/random-n10-0.2.m2")



