GeneNet.wrap <- function(data){
    # method <- "static"
    # p <- 0.8
    ngenes <- dim(data)[2]
    pcor <- ggm.estimate.pcor(data,method ="static",verbose=FALSE)
    test.results <- network.test.edges(pcor,plot=FALSE,verbose=FALSE)
    idx <- which(test.results$prob > 0.8)
    aux<-network.make.graph(test.results[idx,],node.labels = colnames(data))
    net<-as(aux,"matrix")
    return(net);
}
