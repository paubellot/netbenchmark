GeneNet.wrap <- function(data){
    # method <- "static"
    # p <- 0.8
    ngenes <- dim(data)[2]
    pcor <- GeneNet::ggm.estimate.pcor(data,
                                       method ="static",verbose=FALSE)
    C <- fdrtool::fdrtool(corpcor::sm2vec(pcor),plot=FALSE,
                          statistic="correlation",verbose=FALSE)$param
    test.results <- GeneNet::network.test.edges(pcor,plot=FALSE,verbose=FALSE)
    crop <- GeneNet::extract.network(test.results,cutoff.ggm=C[1,"cutoff"],
                                     verbose=FALSE)
    if(dim(crop)[1]==0){
        d <- round(dim(test.results)[1]/10)
        crop <- GeneNet::extract.network(test.results, 
                                         method.ggm="number", cutoff.ggm=d,
                                         verbose=FALSE)
    }
    aux <- GeneNet::network.make.graph(crop,node.labels = colnames(data))
    net <- as(aux,"matrix")
    return(net)
}
