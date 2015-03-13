GeneNet.wrap <- function(data){
    # method <- "static"
    # p <- 0.8
    ngenes <- dim(data)[2]
    pcor <- ggm.estimate.pcor(data,method ="static",verbose=FALSE)
    test.results <- ggm.test.edges(pcor,plot=FALSE,verbose=FALSE)
    idx <- which(test.results$prob > 0.8)
    net <- matrix(0,ngenes,ngenes)
    colnames(net) <- colnames(data)
    rownames(net) <- colnames(data)
    for(i in seq_along(idx)){
        net[test.results[i,2],test.results[i,3]] <- test.results[i,6]
    }
    return(net);
}
