clr.wrap <- function(data){
    # estimator="spearman"
    mim <- minet::build.mim(data,estimator="spearman")
    net <- minet::clr(mim)
    return(net);
}
