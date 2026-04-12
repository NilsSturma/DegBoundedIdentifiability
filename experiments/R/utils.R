library(combinat)

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


is_isomorphic_mixed <- function(g1, g2){
  
  if (nrow(g1$D) != nrow(g2$D)){
    return(FALSE)
  }
  
  n <- nrow(g1$D)
  perms <- permn(n)
  
  for (perm in perms) {
    P <- diag(n)[perm, ]  
    test1 <- all(P %*% g1$D %*% t(P) == g2$D)
    test2 <- all(P %*% g1$B %*% t(P) == g2$B)
    if (test1 & test2){
      return(TRUE)
    }
  }
  return(FALSE)
}