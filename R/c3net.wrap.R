c3net.wrap <- function(data){
    # cutoff <- 0
    net <- c3net(t(data),cutoffMI=0)
    return(net);
}
