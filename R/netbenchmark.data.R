#' Benchmarking of several network inference algorithms for your own data
#' 
#' @param methods A vector of characters containing the names of network i
#' inference algorithms wrappers to be compared (default: "all.fast").
#' @param data data.frame containing the data. Each row should contain a 
#' microarray experiment and each column a gene (default: NULL).
#' @param true.net matrix containg underlying network in the form of adjacency 
#' matrix (default: NULL).
#' @param eval The name of the evaluation metric among the following ones: 
#' "no.truepos", "AUROC" or "AUPR" (default : "AUPR").
#' @param no.topedges Float specifying the percentage number of links to be 
#' considered in the evaluation (default: 20)
#' @param sym Logical specifying if the evaluation is symmetric (default: TRUE)
#'  - see \code{\link{evaluate}}
#' @param plot (default: FALSE)
#' @param verbose Logical specifying if the code should provide a log about 
#' what the function is doing
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
netbenchmark.data <- function(methods="all.fast",data=NULL,true.net=NULL,
                              eval="AUPR",no.topedges=20,sym=TRUE,plot=FALSE,
                              verbose=TRUE)
{
    options(warn=1)
    Fast <- get("Fast", ntb_globals)
    All <- get("All",ntb_globals)
    if(all("all.fast" %in% tolower(methods))) {
        methods <- c(Fast,methods[tolower(methods)!="all.fast"])
    }else if(all("all" %in% tolower(methods))) {
        methods <- c(All,methods[tolower(methods)!="all"])
    }
    if(is.null(data)){
        stop("You should provide the data.")
    }
    if(is.null(true.net)){
        stop("You should provide the underlying network.")
    }
    if(!is.data.frame(data)){
        stop("The provided data shoud be a data.frame.")
    }
    if(!is.matrix(true.net)){
        stop("The provided true.net shoud be a adjacency matrix.")
    }
    if(dim(true.net)[1]!=dim(true.net)[2]){
        stop("The provided true.net shoud be a adjacency matrix.")
    }
    if(all(dim(data)[2]!=dim(true.net)[1])){
        stop("The provided data shoud contain variables in columuns")
    }
    if(all(colnames(data)!=colnames(true.net)[1])){
        stop("The provided data shoud contain variables in columuns")
    }
    nnets <- 1
    ndata <- 1
    nmeths <- length(methods)
    results.table <-  as.data.frame(matrix(0,nrow=nnets,ncol=nmeths+1))
    nlinks.table <- matrix(0,nrow=nnets,ncol=nmeths+1)
    time.table <-  as.data.frame(matrix(0,nrow=nnets,ncol=nmeths))
    names <- as.character(methods)
    colnames(results.table) <- c(names,"rand")
    rownames(results.table) <- as.character(1:nnets)
    colnames(time.table) <- c(names)
    rownames(time.table) <- as.character(1:nnets)
    plots <- list(ndata)
    #########################
    ngenes <- dim(data)[2]
    npos <- sum(true.net) #number of true links in the network
    nlinks <- ngenes^2-ngenes #number of posible links in the network
    if(sym){
        nlinks <- nlinks/2
    }
    no.edges <- round(nlinks*no.topedges/100)
    nneg <- nlinks-npos #number of false links in the network
    conf.mat <-  matrix(0,nrow=nmeths,ncol=4)
    colnames(conf.mat) <- c("TP","FP","TN","FN")
    rownames(conf.mat) <- names
    best <- 0
    if(plot){
        col <- rainbow(nmeths)
    }
    for(j in seq_len(nmeths)){
        if(verbose){
            message(names[j])
        }    
        ptm <- proc.time()
        net <- do.call(names[j],list(data))
        t <- proc.time() - ptm
        if(sum(is.na(net)>0)){
            net[is.na(net)] <- 0
        }
        r <- evaluate(net,true.net,extend=no.edges,sym=sym)
        time.table[1,j]=t[[3]]
        conf.mat[j,] <- r[no.edges,]
        if( tolower(eval)=="no.truepos"){
            results.table[1,j]=mean(r[1:no.edges,"TP"])
        }else if (tolower(eval)== "aupr"){
            results.table[1,j]=aupr(r,no.edges)
        }else if (tolower(eval)== "auroc"){
            results.table[1,j]=auroc(r,no.edges)
        }else stop("unknown evaluation metric")
        if(results.table[1,j]>best){
            best <- results.table[1,j]
            best.net <- net
            M <- j
        }
        if(plot){
            col<-rainbow(nmeths)
            if(j==1){
                pr.plot(r[1:no.edges,],col=col[j],lwd=2)
            }else{
                pr.plot(r[1:no.edges,],col=col[j],lwd=2,device=2)
            }
        }
        table <- as.data.frame(r[1:no.edges,])
        names(table) <- sapply(names(table),tolower)
        plots[[j]] <-  minet::pr(table)
    }
    if(M!=which.max(results.table[1,1:(nmeths)])){
        stop("error")
    }
    rand.net <- matrix(runif(ngenes^2),ngenes,ngenes)
    diag(rand.net) <- 0
    colnames(rand.net) <- colnames(true.net)
    rownames(rand.net) <- colnames(true.net)
    r <- evaluate(rand.net,true.net,extend=no.edges,sym=sym)
    if( tolower(eval)== "no.truepos"){
        results.table[1,nmeths+1]=mean(r[1:no.edges,"TP"])
    }else if (tolower(eval)== "aupr"){
        results.table[1,nmeths+1]=aupr(r,no.edges)
    }else if (tolower(eval)== "auroc"){
        results.table[1,nmeths+1]=auroc(r,no.edges)
    }
    L <- list(results.table,time.table,plots)
    names(L) <- c(paste(eval,"top",as.character(no.topedges),"%",sep=""),
                  "CpuTime","PRcurves")
    return(L)
}