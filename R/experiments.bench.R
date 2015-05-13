experiments.bench <- function(methods="all.fast",datasources.names="all",
    experiments=c(20,50,150),eval="AUPR",no.topedges=20,datasets.num=3,
    local.noise=20,global.noise=0,noiseType="normal",sym=TRUE,seed=NULL)
{
    options(warn=1)
    if(getRversion() >= "2.15.1")  utils::globalVariables(c("Fast", "All"))
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
    points <- length(experiments)
    nmeths <- length(methods)
    ndata <- length(datasources.names) 
    Availabledata <- eval(parse(text="Availabledata"))
    seeds <- as.list(round(runif(length(Availabledata),max=10000)))
    names(seeds) <- Availabledata
    if (!all(datasources.names %in% Availabledata)){
        stop("unknown dataset")
    }
    results <- as.data.frame(matrix(0,points*ndata,nmeths+3))
    pval <- as.data.frame(matrix(0,points*ndata,nmeths+3))
    rown <- character()
    for(n in seq_len(ndata)){
        message(paste("Dataset:",datasources.names[n],"\n"))
        dataset <- eval(parse(text=paste(datasources.names[n],".data",
            sep="")))
        true.net <- eval(parse(text=paste(datasources.names[n],".net",
            sep="")))
        npos <- sum(true.net)
        ngenes <- dim(dataset)[2] #number of genes in the network
        nlinks <- ngenes^2-ngenes #number of posible links in the network
        if(sym){
            nlinks <- nlinks/2
        }
        no.edges <- round(nlinks*no.topedges/100)
        m <- matrix(0,points,nmeths+1)
        pval.table <- matrix(0,points,nmeths+1)
        tp.local.mat <- matrix(0,no.edges,nmeths+1)
        colnames(tp.local.mat) <- c(methods,"rand")
        l.seed <- eval(parse(text=paste("seeds$",datasources.names[n])))
        set.seed(l.seed)
        for(i in seq_len(points)){
            m.local <- matrix(0,datasets.num,nmeths+1)
            rdata <- datasource.subsample(dataset,experiments=experiments[i],
                datasets.num = datasets.num,local.noise = local.noise,
                global.noise =global.noise,noiseType=noiseType,
                samplevar=FALSE)
            for(j in seq_len(nmeths)){
                message(paste(methods[j],"\n"))
                for(k in seq_len(datasets.num)){
                    net <- do.call(methods[j],list(rdata[[k]]))
                    r <- evaluate(net,true.net,extend=no.edges,sym=sym)
                    tp.local.mat[,j] <- tp.local.mat[,j]+r[1:no.edges,"TP"]
                    if(tolower(eval)=="no.truepos"){
                        m[i,j] <- m[i,j]+mean(r[1:no.edges,"TP"])
                        m.local[k,j] <- mean(r[1:no.edges,"TP"])
                    }else if (tolower(eval)== "aupr"){
                        m[i,j] <- m[i,j]+aupr(r,no.edges)
                        m.local[k,j] <- aupr(r,no.edges)
                    }else if (tolower(eval)== "auroc"){
                        m[i,j] <- m[i,j]+auroc(r,no.edges)
                        m.local[k,j] <- auroc(r,no.edges)
                    }else stop("unknown evaluation metric")
                }
                tp.local.mat[,j] <- tp.local.mat[,j]/datasets.num
            }
            m[i,] <- apply(m.local,2,mean)
            M <- which.max(m[i,])
            precision <- tp.local.mat/matrix(rep(1:no.edges,nmeths+1),
                no.edges)
            for(j in seq_len(nmeths)){
                if(j!=M){
                    aux <- wilcox.test(m.local[,j],m.local[,M])
                    pval.table[i,j] <- aux[[3]]
                }else{
                    pval.table[i,j] <- 1
                }
            }
            rand.net <- matrix(runif(ngenes^2),ngenes,ngenes)
            diag(rand.net) <- 0
            colnames(rand.net) <- colnames(net)
            rownames(rand.net) <- colnames(net)
            r <- evaluate(rand.net,true.net,extend=no.edges,sym=sym)
            tp.local.mat[,nmeths+1] <- r[1:no.edges,"TP"]
            precision <- tp.local.mat/matrix(rep(1:no.edges,nmeths+1),
                no.edges)
            if( tolower(eval)=="no.truepos"){
                m[i,nmeths+1]=mean(r[1:no.edges,"TP"])
            }else if (tolower(eval)== "aupr"){
                m[i,nmeths+1]=aupr(r,no.edges)
            }else if (tolower(eval)== "auroc"){
                m[i,nmeths+1]=auroc(r,no.edges)
            }
            aux <- wilcox.test(precision[,nmeths+1],precision[,M])
            pval.table[i,nmeths+1]=aux[[3]]
        }
        rown <- c(rown,rep(datasources.names[n],points))
        results[(1:points)+(n-1)*points,1] <- rep(datasources.names[n],
            points)
        results[(1:points)+(n-1)*points,2] <- experiments
        results[(1:points)+(n-1)*points,3:(nmeths+3)] <- m
        pval[(1:points)+(n-1)*points,1] <- rep(datasources.names[n],points)
        pval[(1:points)+(n-1)*points,2] <- experiments
        pval[(1:points)+(n-1)*points,3:(nmeths+3)] <- pval.table
    }
    colnames(results) <- c("Dataset","experiments",methods,"rand")
    colnames(pval) <- c("Dataset","experiments",methods,"rand")
    list("results"=results,"pval"=pval,"seed"=seed)
}
