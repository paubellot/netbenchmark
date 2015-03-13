zscore.wrap <- function(data){
    Z <- zsc(as.matrix(data))
    Z <- abs(Z)
    colnames(Z) <- colnames(data)
    rownames(Z) <- colnames(data)
    Z <- pmax(Z,t(Z))
    Z
}