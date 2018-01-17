mrnet.wrap <- function(data){
    # estimator <- "spearman"
    mim <- build.mim(data,estimator="spearman")
    net  <- mrnet(mim)
    return(net);
}