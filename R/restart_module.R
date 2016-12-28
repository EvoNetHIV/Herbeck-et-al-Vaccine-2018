#' @export
restart_module <- function(dat, at)
{
# Description
#   Does one of three things, depending on the value of restart_val 
#     "save" -- saves 'dat' to file 'dat.saved' in getwd() at restart_time
#     "restart" --  replaces 'dat' wt 'dat.saved' at restart_time
#     others -- does nothing
#  
# Inputs
#    param$restart_val
#    param$restart_time
#    tx_type_restart
#    
# Outputs:
#    file 'dat.saved' (if restart_val = "save")
#    Entire contents of dat (if restart_val = "restart") 

 if (at==dat$param$restart_time) {
   if (dat$param$restart_val == "save") {
     cat("About to dat to dat.saved\n")
     save(dat, file= paste(getwd(),"/dat.saved",sep=""))
     cat(paste("Done saving dat as",getwd(),"/dat.saved\n",sep=""))
   }
   if (dat$param$restart_val == "restart") {
     cat("*** Replacing 'dat' with contents of\n ")
     cat(paste(getwd(),"/dat.saved ***\n",sep=""))
     tx_type_save <- dat$param$tx_type_restart
     load(file= paste(getwd(),"/dat.saved",sep=""))
     dat$param$tx_type <- tx_type_save
   }
 }
 return(dat)
}
