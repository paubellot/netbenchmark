mrnetb.wrap <- function(data){
    # estimator <- "spearman"
    mim <- minet::build.mim(data,estimator="spearman")
    net  <- minet::mrnetb(mim)
    return(net);
}
