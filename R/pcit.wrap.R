pcit.wrap <- function(data){
    m <- cor(data)
    net  <- pcit(m)
    signif <-  idx(net)
    nonsignif <- idxInvert(nrow(m), signif)
    m.new <- m
    m.new[nonsignif] <- 0
    diag(m.new) <- 0
    m.new <- abs(m.new)
    return(m.new);
}
