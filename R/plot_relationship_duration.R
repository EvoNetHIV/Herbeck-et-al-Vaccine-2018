#' @export
plot_relationship_duration <- function(model){
  #note: parameter "save_network=T" to run fxn
 #plots only first simulation, ignores censored (begin,end) data
  
if(!isTRUE(model$param[[1]]$save_network)){return(invisible(NULL))}
dframe <- networkDynamic::get.edge.activity(model$nw[[1]], 
                                           as.spellList = TRUE, 
                                           active.default = F)

ix <- dframe$duration<1e6
durations <- dframe$duration[ix]
hist(durations,breaks=100,col="blue",xlab="duration",
     main="relationship duration distribution")
}


