clr.wrap <- function(data){
    # estimator="spearman"
    mim <- build.mim(data,estimator="spearman")
    net <- clr(mim)
    return(net);
}
