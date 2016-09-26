                                                      
#' @export
initialize_module <- function(x, param, init, control, s)
{  
  #Description:
  # runs  EpiModel function EpiModel::initialize.net()
  # fills in viral load values for initial infected with initialize_infecteds_vl()
  #	fills in cd4 variable values for initial infected (initial CD4 value based on SPVL)
  # with initialize_infecteds_cd4()
  #	does a few bookkeeping steps with initialize_evonet_misc()
  #	updates vl of initial infecteds (otw, vl jumps from initial value to spvl in one
  # timestep) with viral_update_module_gamma()
  #	creates and inserts initial values for "popsumm" stats (calculated at end of
  # each timestep_ with initialize_popsumm_dynamics()
  
  #sets up basic EpiModel structure
  dat  <-  EpiModel::initialize.net(x, param, init, control,s)
  #sets up agent attributes and initial values 
  dat  <-  initialize_agents(dat, 1)
  #fills in vl variable values for initial infecteds
  dat  <-  initialize_infecteds_vl(dat,1)
  #fills in cd4 variable values for initial infecteds
  #note, this must come after initialize_infecteds_vl because initial cd4 based on spvl 
  dat  <-  initialize_infecteds_cd4(dat,1)
  #does a few random bookkeeping steps
  dat  <-   initialize_evonet_misc(dat)
  #updates vl of initial infecteds (otw, jumps from initial value to spvl)
  #but shouldn't happen for aim3 runs
  if(param$VL_Function != "aim3"){
  dat  <-  viral_update_gamma(dat,1)
  dat  <-  viral_update_cd4_intial_pop(dat)
  }
  #create list, if "save_vl_list" = TRUE to save
  #individual agent's vl/cd4 values for each timestep
  dat<- summary_vl_list(dat,1)
  #creates and fills in initial values for "popsumm" stats (stats calculated at
  #end of each timestep)
  dat <- summary_popsumm(dat,1)
  #keep track of current simulation/replicate
  dat$simulation <- s
  return(dat)
}
