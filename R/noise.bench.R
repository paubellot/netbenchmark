noise.bench <- function(methods = "all.fast", datasources.names = "all", 
                        eval = "AUPR", no.topedges = 20, experiments = 150,
                        datasets.num = 3, local.noise = seq(0, 100, len = 3),
                        global.noise = 0, noiseType = "normal", sym = TRUE,
                        seed = NULL, verbose = TRUE)
{
    options(warn=1)
    Fast <- get("Fast", ntb_globals)
    All <- get("All",ntb_globals)
    # set random number generator seed if seed is given
    if (!is.null(seed)){
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
                experiments<-NULL
            }
        }
    }
    if(length(local.noise)!=length(global.noise)){
        if(length(local.noise)==1){
            points <- length(global.noise)
            local.noise <- rep(local.noise,points)
        }else if(length(global.noise)==1){
            points <- length(local.noise)
            global.noise <- rep(global.noise,points)
        }else{
            stop("Error: mismatch in lengths of local.noise and global.noise")
        }
    }
    nmeths <- length(methods)
    ndata <- length(datasources.names)
    results <- as.data.frame(matrix(0,points*ndata,nmeths+4))
    pval <- as.data.frame(matrix(0,points*ndata,nmeths+4))
    Availabledata <- eval(parse(text="Availabledata"))
    seeds <- as.list(round(runif(length(Availabledata),max=1e9)))
    names(seeds) <- Availabledata
    rown <- character()
    if (!all(datasources.names %in% Availabledata)) stop("unknown datasource")
    for(n in seq_len(ndata)){
        if(verbose){
            message(paste("Datasource:",datasources.names[n]))
        }
        aux <- grndata::getData(datasources.names[n])
        datasource <- aux[[1]]
        true.net <- aux[[2]]
        s <- dim(datasource)
        ngenes <- s[2] #number of genes in the network
        npos <- sum(true.net) #number of true links in the network
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
        if(is.null(experiments)){
            spd <- vector(mode = "list",length = datasets.num)
            for(i in seq_len(datasets.num)){
                spd[[i]] <- datasource
            }
        }else{
            spd <- datasource.subsample(datasource,experiments=experiments,
                                        datasets.num = datasets.num,
                                        local.noise = 0,global.noise = 0)
        }
        for(i in seq_len(points)){
            m.local <- matrix(0,datasets.num,nmeths+1)
            rdata <- vector('list',datasets.num)
            for(k in seq_len(datasets.num)){
                if(local.noise[i]!=0){
                    rdata[[k]] <- apply(spd[[k]],2,.cont,
                                        noise=local.noise[i],noiseType=noiseType)
                }else{
                    rdata[[k]] <- spd[[k]]
                }
                if(global.noise[i]!=0){
                    sds <- apply(spd[[k]], 2, sd)
                    if(noiseType=="normal"){
                        Gnoise <- matrix(rnorm(s[1]*s[2],mean=0,
                                               sd=mean(sds)*global.noise[i]/100),
                                         s[1], s[2])
                    }
                    if(noiseType=="lognormal"){
                        Gnoise <- matrix(rlnorm(s[1]*s[2],meanlog=0,
                                                sdlog=mean(sds)*global.noise[i]/100),
                                         s[1], s[2])
                    }
                    rdata[[k]] <- rdata[[k]]+Gnoise
                }
            }
            for(j in seq_len(nmeths)){
                if(verbose){
                    message(methods[j])
                }
                for(k in seq_len(datasets.num)){
                    net <- do.call(methods[j],list(rdata[[k]]))
                    r <- evaluate(net,true.net,extend=no.edges,sym=sym)
                    tp.local.mat[,j] <- tp.local.mat[,j]+r[1:no.edges,"TP"]
                    if(tolower(eval)=="no.truepos"){
                        m.local[k,j] <- mean(r[1:no.edges,"TP"])
                    }else if (tolower(eval)== "aupr"){
                        m.local[k,j] <- aupr(r,no.edges)
                    }else if (tolower(eval)== "auroc"){
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
                    pval.table[i,j]=aux[[3]]
                }else{
                    pval.table[i,j]=1
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
            if(tolower(eval)=="no.truepos"){
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
        results[(1:points)+(n-1)*points,2] <- local.noise
        results[(1:points)+(n-1)*points,3] <- global.noise
        results[(1:points)+(n-1)*points,4:(nmeths+4)] <- m  
        pval[(1:points)+(n-1)*points,1] <- rep(datasources.names[n],points)
        pval[(1:points)+(n-1)*points,2] <- local.noise
        pval[(1:points)+(n-1)*points,3] <- global.noise
        pval[(1:points)+(n-1)*points,4:(nmeths+4)] <- pval.table
    }
    colnames(results) <- c("Datasource","local.noise","global.noise",
                           methods,"rand")
    colnames(pval) <- c("Datasource","local.noise","global.noise",
                        methods,"rand")
    list("results"=results,"pval"=pval,"seed"=seed)
}

.cont <- function(x,noise=0,noiseType="normal"){
    s.d <- noise*sd(x)/100
    if(noiseType=="normal"){
        n <- rnorm(length(x),mean=0,sd=s.d)
    }
    if(noiseType=="lognormal"){
        n <- rlnorm(length(x), meanlog = 0, sdlog = s.d)
    }
    return(x+n)
}
