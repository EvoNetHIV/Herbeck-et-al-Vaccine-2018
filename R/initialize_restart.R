#' @export
initialize_restart<-function(x,param,init,control,s){
  #------------------------------
  #code directly from EpiModel for sim. restart
  dat <- list()
  dat$nw <- x$network[[s]]
  dat$param <- x$param[[s]] #modified slighly from original EpiModel code (dat$param <- x$param)
  dat$control <- control
  dat$nwparam <- x$nwparam
  dat$epi <- sapply(x$epi, function(var) var[s])
  names(dat$epi) <- names(x$epi)
  dat$attr <- x$attr[[s]]
  dat$stats <- sapply(x$stats, function(var) var[[s]])
  dat$temp <- list()
  #---------------------
  # evonet code for sim. restart
  dat$pop<-x$pop[[s]]
  dat$popsumm <-x$popsumm[[s]]
  dat$coital_acts_list <- x$coital_acts_list[[s]]
  dat$vl_list <- x$vl_list[[s]]
  dat$InfMat <- dat$InfMat[[s]]
  
  #attach fxns to calculte population statistics to "dat" object
  popsumm_fxns <- summary_popsumm_fxns(generic_nodal_att_values= param$generic_nodal_att_values)
  dat$popsumm_fxns <- lapply(1:length(popsumm_fxns),function(x) popsumm_fxns[[x]]$model_value_fxn)
  names(dat$popsumm_fxns)<- names(popsumm_fxns)
  #---------------------
  return(dat)
}
