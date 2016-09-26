#' @export
initialize_evonet_misc <- function(dat)
{
  #Description:
  # Miscellaneous internal bookkeeping for model setup:
  
  # 1) #create "dat$attr$status_evo" vector which tracks status of agent on network: 
  # infected(1) ,susceptible(0),dead aids(-2),dead other(-1), aged out (-1.5) 
  # note: dead/aged-out agents removed at end of daily timestep
  # note: dat$attr$status is an EpiModel object but uses different terms than in evonet
  
  # 2) Creates dat$attr$id vector which tracks id of agent on network. Allows to track agent info
  #in "pop" list to network; initial "id" value created in "vital_additions_new()"
  
  #3) Creates "role","sex","att1" vectors which store agent info for three variables on network
   #att1 variable is generic variable that can be used to group agents (eg, risk group, geography)
  # Creates "dat$attr$evo_status” vector for status of agents on network 
  # Sets "agent id" nodal attribute on network
  #Sets sex on network
  # VL list (if flagged): save VL/CD4 data for each agent per timestep
  # QAQC steps (optional): coital_acts_list (all coital acts info per timestep), and 
  # infection matrix (InfMat) (list of newly infected agent and their partner per  timestep)
  
  #------------------------------
  
  #1)
  dat$attr$status_evo <- rep(0,length(dat$attr$status))
  dat$attr$status_evo[which(dat$attr$status=="i" )] <- 1
  
  #2)
  dat$attr$id <- dat$pop$id
  
  #note: values for "att1","role","sex" created in "vital_additions_new()" in "initial" section
  #3)
  if(!is.logical(dat$param$generic_nodal_att_values))
    dat$attr$att1 <- dat$pop$att1
  if(!is.logical(dat$param$role_props) && dat$param$model_sex=="msm")
    dat$attr$role <- dat$pop$role
  if(dat$param$model_sex!="msm")
    dat$attr$sex <- dat$pop$sex
  
  #------------------------------
  
  #create (or not) coital acts list, which is used to save (if flagged)
  #coital acts df for each time step for qaqc purposes
  if(dat$param$save_coital_acts)
    dat$coital_acts_list <- list()
  else
    dat$coital_acts_list  <- NULL
  
  #-----
  #dat$InfMat: data object filled in "transmission_bookkeeping_module"
  #if TRUE, then timestep, id of newly infected agent, and id of infector saved
  #for qaqc purposes
  if(dat$param$save_infection_matrix)
    dat$InfMat <- list()
  else 
    dat$InfMat<-NULL
  
  #---
  #attach fxns to calculte population statistics to "dat" object
  #and plotting methods; see "summary_popsumm_fxns"
  popsumm_fxns <- summary_popsumm_fxns(generic_nodal_att_values=dat$param$generic_nodal_att_values)
  dat$popsumm_fxns <- lapply(1:length(popsumm_fxns),function(x) popsumm_fxns[[x]]$model_value_fxn)
  names(dat$popsumm_fxns)<- names(popsumm_fxns)
  
  #create an empty list based on number of variables, fill in with NAs                                                          
  popsumm_vars     <- names(popsumm_fxns)
  dat$popsumm        <- vector('list', length(popsumm_vars))
  
  #popsumm length is how many times summary stats will be recorded
  #during model run; popsumm_frequency=1 means stats recorded
  #each timestep; setting popsumm_frequency >1 (eg, 30) saves time
  #by skipping calculations of summary stats
  if(dat$param$popsumm_frequency==1)
    popsumm_length <- dat$param$n_steps
  else
    popsumm_length <- (dat$param$n_steps/dat$param$popsumm_frequency)+1
  
  dat$popsumm        <- lapply(dat$popsumm,function(x){ rep(NA_real_,times=popsumm_length)})
  names(dat$popsumm) <- popsumm_vars
  #-----------------------------
  #create age list for plotting age distributions during model
  #run (at start,1/4,1/2,3/4, end of model time period)
  dat$age_list<-vector('list',length=5)
  
  
  
   return(dat)
  
}