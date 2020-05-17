.onLoad <- function(libname, pkgname) {
    assign("ntb_globals", new.env(), envir=parent.env(environment()))
    All <- c("aracne.wrap","c3net.wrap","clr.wrap","GeneNet.wrap",
             "Genie3.wrap","mrnet.wrap",
             "mutrank.wrap","mrnetb.wrap","zscore.wrap")
    Fast <- c("aracne.wrap","c3net.wrap","clr.wrap",
             "GeneNet.wrap","mrnetb.wrap","mutrank.wrap",
             "zscore.wrap")
    assign("All",All,ntb_globals)
    assign("Fast",Fast,ntb_globals)
}