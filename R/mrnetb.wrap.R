mrnetb.wrap <- function(data){
    # estimator <- "spearman"
    mim <- build.mim(data,estimator="spearman")
    net  <- mrnetb(mim)
    return(net);
}
