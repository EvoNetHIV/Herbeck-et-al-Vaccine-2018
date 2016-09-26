
#' @export
setup_epimodel_control_object <- function(evonet_params,module_list)
{ 

  #Description:
  # Helper function to fill in and run EpiModel’s control.net(); 
  # of minor utility, maybe should be removed and control.net just called straight from master_script.
  # Called in master_script
  #calls epimodel::control.net, ?EpiModel::control.net for details
  #input: module_list (modules to be called)
  #input: evoparams "save_network","nsims","n_steps"
  #output: EpiModel "control" object (used in EpiModel netsim fxn)
  
  control_epimodel_params_list <- list(ncores=evoparams$ncores,
                                       type   = "SI",
                                       nsims  = evonet_params$nsims,
                                       nsteps = evonet_params$n_steps,  #set in parameter_list.R
                                        start  = evonet_params$start_timestep,
                                       depend = TRUE,
                                       tea.status    = FALSE,
                                       save.transmat = FALSE,
                                       save.nwstats  = FALSE,
                                       save.network  = evonet_params$save_network,
                                       delete.nodes  = !evonet_params$save_network,
                                       module.order  = names(module_list)[-c(1,length(names(module_list)))],
                                       save.other    = c("attr","pop","param","nw","coital_acts_list",
                                                         "popsumm","vl_list","InfMat","age_list",
                                                         "sessionInfo","evonet_version"))

  control <- do.call(EpiModel::control.net,  
                     c(module_list,control_epimodel_params_list) )
  
  return(control)
}
###############################################
