RegistrerWrapper<-function(wrapper.name=NULL,all.fast=TRUE){
    tst<-paste("exists('",wrapper.name,"')",sep="")
    if(eval(parse(text=tst))){
        if(all.fast){
            Fast<-get("Fast", ntb_globals)
            assign("Fast", c(Fast,wrapper.name),  ntb_globals)
            message(paste(wrapper.name, "added to all.fast methods"))
            message("The methods registered in all.fast are:")
            message(c(Fast,wrapper.name))
        }else{
            All<-get("All",ntb_globals)
            assign("All",c(All,wrapper.name),  ntb_globals)
            message(paste(wrapper.name, "added to all methods"))
            message("The methods registered in all are:")
            message(c(All,wrapper.name))
        }
    }else{
        stop(paste(wrapper.name,"does not exist!"))
    }
}

UnregistrerWrapper<-function(wrapper.name=NULL,all.fast=TRUE){
    if(all.fast){
        Fast <- get("Fast", ntb_globals)
        idx <- which(wrapper.name==Fast)
        if(idx>0){
            Fast<-Fast[-idx]
            assign("Fast", c(Fast,wrapper.name),  ntb_globals)
            message(paste(wrapper.name, "removed from all.fast methods"))
        }
    }else{
        All <- get("All", ntb_globals)
        idx <- which(wrapper.name==All)
        if(idx>0){
            All<-All[-idx]
            assign("All", c(All,wrapper.name),  ntb_globals)
            message(paste(wrapper.name, "removed from all methods"))
        }
    }
}