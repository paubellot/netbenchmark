c3net.wrap <- function(data){
    mim <- c3net::makemim(t(data))
    ccc <- mim[upper.tri(mim)]
    ccc <- ccc[ccc != 0] #diagonal removed
    ccc <- abs(ccc)
    x <- sort(ccc, decreasing=T)
    a<- nrow(mim)
    a<- a*(a-1)/2
    RNratio <- round(0.02*a)
    cutoff <- x[RNratio]
    net <- c3net::c3net(t(data),cutoffMI=cutoff)
    return(net);
}
