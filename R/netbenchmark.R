netbenchmark <- function(methods="all.fast",datasources.names="all",
    experiments=150,eval="AUPR",no.topedges=20,datasets.num=5,local.noise=20,
    global.noise=0,noiseType="normal",sym=TRUE,plot=FALSE,seed=NULL,verbose=TRUE)
{
    options(warn=1)
    # set random number generator seed if seed is given
    if (!is.null(seed)) {
        set.seed(seed)
    }else{
        seed <- as.double(Sys.time())
        set.seed(seed)
    }
    if(all("all.fast" %in% tolower(methods))) {
        methods <- c(Fast,methods[tolower(methods)!="all.fast"])
    }else if(all("all" %in% tolower(methods))) {
        methods <- c(All,methods[tolower(methods)!="all"])
    }
    if(length(datasources.names)==1){
        if (tolower(datasources.names)=="all"){
            datasources.names <- c("rogers1000","syntren1000","syntren300",
                "gnw1565","gnw2000")
        }else{
            if(tolower(datasources.names)=="toy"){
                datasources.names <- "toy"
            }
        }
    }
    Availabledata <- eval(parse(text="Availabledata"))
    if(!all(datasources.names %in% Availabledata)){
        stop("The specified datasources are not available")
    }
    seeds <- as.list(round(runif(length(Availabledata),max=1e9)))
    names(seeds) <- Availabledata
    ndata <- length(datasources.names)
    if(!is.na(experiments)){
        for(n in seq_len(ndata)){
            #loading the whole datasource
            datasource <- eval(parse(text=paste(datasources.names[n],
            ".data",sep=""))) 
            if(experiments*datasets.num>dim(datasource)[1]){
                warning(paste("The specified number of experiments and 
                    datasets is bigger than the orginal number of experiments 
                    in the datasource: ",datasources.names[n],", 
                    sampling with replacement will be used",sep=""))
            } 
        }
    }
    nnets <- datasets.num*ndata
    nmeths <- length(methods)
    results.table <-  as.data.frame(matrix(0,nrow=nnets,ncol=nmeths+3))
    nlinks.table <- matrix(0,nrow=nnets,ncol=nmeths+1)
    npos.table <- matrix(0,nrow=nnets,ncol=nmeths+1)
    pval.table <-  as.data.frame(matrix(0,nrow=ndata,ncol=nmeths+2))
    time.table <-  as.data.frame(matrix(0,nrow=nnets,ncol=nmeths+2))
    names <- as.character(methods)
    colnames(results.table) <- c("Origin","experiments",names,"rand")
    rownames(results.table) <- as.character(1:nnets)
    colnames(pval.table) <- c("Origin",names,"rand")
    rownames(pval.table) <- as.character(1:ndata)
    colnames(time.table) <- c("Origin","experiments",names)
    rownames(time.table) <- as.character(1:nnets)
    plots <- list(ndata)
    for(n in seq_len(ndata)){ #for each of the datasources
        if(verbose){
            message(paste("datasource:",datasources.names[n]))
        }
        datasource <- eval(parse(text=paste(datasources.names[n],
            ".data",sep="")))
        true.net <- eval(parse(text=paste(datasources.names[n],".net",
            sep="")))
        #loading the specified datasource
        ngenes <- dim(datasource)[2]
        npos <- sum(true.net) #number of true links in the network
        nlinks <- ngenes^2-ngenes #number of posible links in the network
        if(sym){
            nlinks <- nlinks/2
        }
        no.edges <- round(nlinks*no.topedges/100)
        nneg <- nlinks-npos #number of false links in the network
        l.seed <- eval(parse(text=paste("seeds$",datasources.names[n])))
        set.seed(l.seed)
        data.list <- datasource.subsample(datasource,
            datasets.num=datasets.num,experiments=experiments,
            local.noise=local.noise,global.noise=global.noise,
            noiseType=noiseType)
        nlinks.table[(1:datasets.num)+(n-1)*datasets.num,] <- matrix(no.edges,
            datasets.num,nmeths+1)
        npos.table[(1:datasets.num)+(n-1)*datasets.num,] <- matrix(npos,
            datasets.num,nmeths+1)
        tp.mat <- matrix(0,no.edges,nmeths+1)
        for(i in seq_along(data.list)){
            tp.local.mat <- matrix(0,no.edges,nmeths+1)
            colnames(tp.local.mat) <- c(methods,"rand")
            sd <- data.list[[i]] 
            c.row <- (n-1)*datasets.num+i  
            results.table[c.row,1]=datasources.names[n]
            results.table[c.row,2]=dim(sd)[1]
            if(verbose){
                message(paste(" dataset:",as.character(i)))
            }
            conf.mat <-  matrix(0,nrow=nmeths,ncol=4)
            colnames(conf.mat) <- c("TP","FP","TN","FN")
            rownames(conf.mat) <- names
            best <- 0
            for(j in seq_len(nmeths)){
                if(verbose){
                    message(names[j])
                }    
                ptm <- proc.time()
                net <- do.call(names[j],list(sd))
                t <- proc.time() - ptm
                if(sum(is.na(net)>0)){
                    net[is.na(net)] <- 0
                }
                r <- evaluate(net,true.net,extend=no.edges,sym=sym)
                tp.local.mat[,j] <- r[1:no.edges,"TP"]
                time.table[c.row,j+2]=t[[3]]
                conf.mat[j,] <- r[no.edges,]
                if( tolower(eval)=="no.truepos"){
                    results.table[c.row,j+2]=mean(r[1:no.edges,"TP"])
                }else if (tolower(eval)== "aupr"){
                    results.table[c.row,j+2]=aupr(r,no.edges)
                }else if (tolower(eval)== "auroc"){
                    results.table[c.row,j+2]=auroc(r,no.edges)
                }else stop("unknown evaluation metric")
                if(results.table[c.row,j+2]>best){
                    best <- results.table[c.row,j+2]
                    best.net <- net
                    M <- j
                }
            }
            if(M!=which.max(results.table[c.row,3:(nmeths+2)])){
                stop("error")
            }
            rand.net <- matrix(runif(ngenes^2),ngenes,ngenes)
            diag(rand.net) <- 0
            colnames(rand.net) <- colnames(net)
            rownames(rand.net) <- colnames(net)
            r <- evaluate(rand.net,true.net,extend=no.edges,sym=sym)
            tp.local.mat[,nmeths+1] <- r[1:no.edges,"TP"]
            if( tolower(eval)== "no.truepos"){
                results.table[c.row,nmeths+3]=mean(r[1:no.edges,"TP"])
            }else if (tolower(eval)== "aupr"){
                results.table[c.row,nmeths+3]=aupr(r,no.edges)
            }else if (tolower(eval)== "auroc"){
                results.table[c.row,nmeths+3]=auroc(r,no.edges)
            }
            tp.mat <- tp.mat+tp.local.mat
        }
        res <- results.table[(1:datasets.num)+datasets.num*(n-1),-(1:2)]
        M <- which.max(apply(res,2,mean))
        for(j in seq_len(nmeths)){
            if(j!=M){
                aux <- wilcox.test(res[,j],res[,M])
                pval.table[,j+1]=aux[[3]]
            }else{
                pval.table[,j+1]=1
            }
        }
        colnames(tp.mat) <- c(methods,"rand")
        m.pr <- .get.pr(tp.mat/datasets.num,npos)
        if(plot){
            .get.pr.plot(m.pr,type="l",datasource.name=datasources.names[n],
                lwd=2.5)
        }
        plots[[n]] <- m.pr
    }
    names(plots) <- datasources.names
    pval.table[,1] <- datasources.names
    time.table[,1:2] <- results.table[,1:2]
    aux.table <- results.table[,3:(nmeths+3)]
    mean.table <- results.table[1:2,3:(nmeths+3)]
    mean.table[1,] <- apply(aux.table/nlinks.table,2,mean)
    mean.table[2,] <- apply(aux.table/nlinks.table,2,sd)
    m <- c("Mean","Standard Deviation")
    rownames(mean.table) <- m
    L <- list(results.table,pval.table,mean.table,time.table,plots,seed)
    names(L) <- c(paste(eval,"top",as.character(no.topedges),"%",sep=""),
        "pval","summary","CpuTime","PRcurves","seed")
    return(L)
}
