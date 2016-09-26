#' @export
plot_age_hist<-function(model)
{
  #plots age distribution from dat$age_list object which is 
  #output data structure e.g, evomodel$age_list; dat$age_list is
  #initialized in "initialize_evonet_misc" and filled in 
  #"summary_misc")
  
  yearvec=c(0,
            round((model$param[[1]]$n_steps*c(.25,.5,.75,1))/365,1))
  
  age_list <- model$age_list
  for(jj in 1:5){
    out<-NULL
    for(ii in 1:length(age_list)){
      out<- c(out,unlist(age_list[[ii]][[jj]]))
    }
  hist(out,
       breaks=length(model$param[[1]]$min_age:model$param[[1]]$max_age),
       col="blue",main=NA,xlab="age")
  title(paste("age distribution at year",yearvec[jj]))
  }
}