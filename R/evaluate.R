evaluate <- function(inf.net, true.net,sym=TRUE,extend=0) 
{
    if((diff(dim(inf.net))!=0) && dim(inf.net)[2]==3){
        ngenes <- max(as.double(inf.net[,1:2]))
        if((diff(dim(true.net))!=0) && dim(true.net)[2]==3){
            ngenes <-max(as.double(true.net[,1:2]),ngenes)
            true.net <- Matrix::sparseMatrix(i=true.net[,1],
                                             j=true.net[,2],x=true.net[,3],
                                             dims =c(ngenes,ngenes))
            true.net<-as.matrix(true.net)
        }else{
            ngenes <- dim(true.net)[2]
        }
        inf.net <- Matrix::sparseMatrix(i=inf.net[,1],
                                        j=inf.net[,2],x=inf.net[,3],
                                        dims =c(ngenes,ngenes))
        inf.net<-as.matrix(inf.net)
    }else{
        ngenes <- dim(true.net)[2]
    }
    if(sym){
        if(extend>(ngenes*(ngenes+1)/2)){
            extend <- (ngenes*(ngenes-1)/2)
        }
        true.net <- pmax(true.net,t(true.net))
        inf.net <- pmax(inf.net,t(inf.net))
        if(is.null(colnames(true.net))){
            colnames(true.net) <- as.character(1:ngenes)
            rownames(true.net) <- as.character(1:ngenes)
            colnames(inf.net) <- as.character(1:ngenes)
            rownames(inf.net) <- as.character(1:ngenes)
        }else if(is.null(colnames(inf.net))){
            colnames(inf.net) <- colnames(true.net)
            rownames(inf.net) <- rownames(true.net)
        }
        true.net[lower.tri(true.net,diag=TRUE)] <- 0
        inf.net[lower.tri(inf.net,diag=TRUE)] <- 0
        if(sum(inf.net!=0)>0){
            PredEdgeList <- .Adj2Edgelist(inf.net)
        }
        if(sum(inf.net!=0)<extend){
            inv.net <- inf.net*0
            inv.net[PCIT::idxInvert(inf.net,which(inf.net!=0))] <- 1
            diag(inv.net) <- 0
            inv.net[lower.tri(inv.net)] <- 0
            E <- .Adj2Edgelist(inv.net)
            idx <- dim(E)[1]
            E <- E[sample(idx),]
            if(sum(inf.net!=0)==0){
                w <- seq(1,0.01,length.out = extend-sum(inf.net!=0)+1)
                E <- E[1:(extend-sum(inf.net!=0)+1),]
                E[,3] <- w
                PredEdgeList <- E[1:(extend-sum(inf.net!=0)+1),]
            }else{
                w <- seq(min(inf.net[which(inf.net!=0)])-0.01,0.01,length.out = extend-sum(inf.net!=0))
                E <- E[1:(extend-sum(inf.net!=0)),]
                E[,3] <- w
                PredEdgeList <- rbind(PredEdgeList,
                    E[1:(extend-sum(inf.net!=0)),])
            }
        }
        GSEdgeList <- .Adj2Edgelist(true.net)
        r <- rate(PredEdgeList, GSEdgeList, ngenes, 2)
    }else{
        if(extend>(ngenes^2-ngenes)){
            extend <- (ngenes^2-ngenes)
        }
        if(is.null(colnames(true.net))){
            colnames(true.net) <- as.character(1:ngenes)
            rownames(true.net) <- as.character(1:ngenes)
            colnames(inf.net) <- as.character(1:ngenes)
            rownames(inf.net) <- as.character(1:ngenes)
        }else if(is.null(colnames(inf.net))){
            colnames(inf.net) <- colnames(true.net)
            rownames(inf.net) <- rownames(true.net)
        }
        if(sum(inf.net!=0)>0){
            PredEdgeList <- .Adj2Edgelist(inf.net)
        }
        if(sum(inf.net!=0)<extend){
            inv.net <- inf.net*0
            inv.net[idxInvert(inf.net,which(inf.net!=0))] <- 1
            diag(inv.net) <- 0
            E <- .Adj2Edgelist(inv.net)
            idx <- dim(E)[1]
            E <- E[sample(idx),]
            if(sum(inf.net!=0)==0){
                w <- seq(1,0.01,length.out = extend-sum(inf.net!=0))
                E <- E[1:(extend-sum(inf.net!=0)),]
                E[,3] <- w
                PredEdgeList <- E[1:(extend-sum(inf.net!=0)),]
            }else{
                w <- seq(min(inf.net[which(inf.net!=0)])-0.01,0.01,length.out = extend-sum(inf.net!=0))
                E <- E[1:(extend-sum(inf.net!=0)),]
                E[,3] <- w
                PredEdgeList <- rbind(PredEdgeList,
                    E[1:(extend-sum(inf.net!=0)),])
            }
        }
        GSEdgeList <- .Adj2Edgelist(true.net)
        r <- rate(PredEdgeList, GSEdgeList, ngenes, 1)
    }
    colnames(r) <- c("TP","FP","TN","FN")
    return(r)
}

.Adj2Edgelist <- function(Adjmat){
    idx <- which(Adjmat!=0, arr.ind=TRUE)
    r=idx[,1]; c=idx[,2];
    E <- as.character(rep(0,dim(idx)[1]*3))
    l <- dim(idx)[1]
    dim(E) <- c(l,3)
    names <- colnames(Adjmat)
    A <- rep(0,l)
    for(i in seq_len(l)){
        E[i,1] <- names[r[i]]
        E[i,2] <- names[c[i]]
        A[i] <- Adjmat[r[i],c[i]]
    }
    a <- sort(A,decreasing = TRUE,index.return = TRUE)
    E <- E[a$ix,]
    dim(E) <- c(l,3)
    E[,3] <- as.character(a$x)
    return(E)
}
