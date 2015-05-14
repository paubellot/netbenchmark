RegistrerWrapper<-function(wrapper.name=NULL,all.fast=TRUE){
    tst<-paste("exists('",wrapper.name,"')",sep="")
    if(eval(parse(text=tst))){
        if(all.fast){
            Fast<-get("Fast", ntb_globals)
            assign("Fast", c(Fast,wrapper.name),  ntb_globals)
            message(paste(wrapper.name, "added to all.fast methods"))
        }else{
            All<-get("All",ntb_globals)
            assign("All",c(All,wrapper.name),  ntb_globals)
            message(paste(wrapper.name, "added to all methods"))
        }
    }else{
        stop(paste(wrapper.name,"does not exist!"))
    }
}
