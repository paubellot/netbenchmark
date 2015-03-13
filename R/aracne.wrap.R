aracne.wrap <- function(data){
    #estimator="spearman"
    #eps=0
    mim <- build.mim(data,estimator="spearman")
    net <- aracne(mim,eps=0)
    return(net);
}
