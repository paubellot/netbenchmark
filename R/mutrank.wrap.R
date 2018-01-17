mutrank.wrap <- function(data){
    m <- cor(data,method = "pearson")
    r <- apply(m,1,rank)
    diag(r) <- 0
    net <- r*t(r)/2
    colnames(net) <- colnames(data)
    rownames(net) <- colnames(data)
    return(net)
}
