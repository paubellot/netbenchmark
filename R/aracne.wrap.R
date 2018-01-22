aracne.wrap <- function(data){
    #estimator="spearman"
    #eps=0
    mim <- minet::build.mim(data,estimator="spearman")
    net <- minet::aracne(mim,eps=0)
    return(net);
}
