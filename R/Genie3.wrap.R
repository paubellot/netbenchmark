Genie3.wrap <- function(data){
    net <- GENIE3::GENIE3(t(data))
    return(net);
}

