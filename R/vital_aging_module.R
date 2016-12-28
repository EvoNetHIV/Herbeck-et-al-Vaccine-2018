
###################################################################
# 
#' @export
vital_aging_module <- function(dat,at){
  
  ix <- which(is.element(dat$pop$Status, c(0:1)))
  dat$pop$age[ix] <- round( dat$pop$age[ix] + (1/365) ,5 )
  mm <- match(ix,dat$attr$id)
  
  #temp qaqc
  if(any(is.na(mm))){browser()}
  
  dat$attr$age[mm]=dat$pop$age[ix]
  if(!is.null(dat[['nw']])){
    network::set.vertex.attribute(x = dat$nw, attr = "age",
                                value = dat$attr$age)
  }
  return(dat)  
}
