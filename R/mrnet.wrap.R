mrnet.wrap <- function(data){
    # estimator <- "spearman"
    mim <- minet::build.mim(data,estimator="spearman")
    net  <- minet::mrnet(mim)
    return(net);
}