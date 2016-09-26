
###################################################################
# 
#' @export
vital_aging_module <- function(dat,at){
  
  ix <- is.element(dat$pop$Status, c(0:1))
  
    dat$pop$age[ix] <- round( dat$pop$age[ix] + (1/365) ,5 )
  
  return(dat)  
}
