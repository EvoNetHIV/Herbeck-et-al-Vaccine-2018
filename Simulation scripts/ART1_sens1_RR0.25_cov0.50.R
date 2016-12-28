hyak=T
hyak_par=T

if(!isTRUE(hyak) & isTRUE(hyak_par)){stop("hyak flags incorrect")}

hyak_path= '/gscratch/csde/kpeebles'
local_path="H://evonet/trunk/scripts/sandbox josh/No treatment/Hyak"

if(hyak)outpath=hyak_path else 
  outpath=local_path

#--------------------------------------------------------------
library(evonet)
library(EpiModelHPC)

#--------------------------------------------------------------
#Load default parameters

primary_parameters  <- input_parameters_primary()
cd4_data            <- input_parameters_cd4_data()

#--- combine individual parameters into single list
evoparams <- c(primary_parameters, cd4_data)

#--------------------------------------------------------------
evoparams$nsims            = 20

if(hyak_par){
  evoparams$ncores=16
}else{evoparams$ncores=1}

# -- ART parameters --------------------------------------------- #
# ART scenario 1
evoparams$tx_type = "cd4_and_time_dist"
evoparams$cd4_trt_guidelines_chgs = list(3:4, 1:4)
evoparams$start_treatment_campaign = c(7*365, 12*365)
evoparams$mean_time_tx = 365*3
evoparams$sd_time_tx = 365
evoparams$prob_eligible_ART = 0.40

# -- Vaccine parameters ----------------------------------------- #
# Vaccine scenario 1
evoparams$start_vacc_campaign = (22*365):(40*365)
evoparams$perc_virus_vaccine_sens <- 1
evoparams$trans_RR_vaccine <- 0.25
evoparams$perc_vaccinated <- (0.50/(1-0.50))/(3*365)
evoparams$vacc_eff_duration = 365*3

# -- Default parameters to override ----------------------------- #
evoparams$initial_pop      = 10000
evoparams$initial_infected = 2000
evoparams$target_stats     = .35*evoparams$initial_pop

evoparams$n_steps          = 365*40
evoparams$popsumm_frequency = 45
evoparams$poisson_birth_lambda = 0.0137*(evoparams$initial_pop/100)
evoparams$trans_RR_age = 1.0
evoparams$fast_edgelist = T
evoparams$vl_peak_agent_flag = T

#add parameters that are functions of other input parameters
evoparams  <- input_parameters_derived(evoparams)
evoparams$age_dist <- seq(50, 10, -10/9)/1110

#convert raw parameter list into EpiModel object
evoparams <- do.call(EpiModel::param.net,evoparams)

# Load network
load(file.path(outpath,"vax_nw_10k.RDATA"))

model_name = paste("ART1_sens", evoparams$perc_virus_vaccine_sens,
                   "_RR", evoparams$trans_RR_vaccine,
                   "_cov0.50",
                   sep = "")

#--------------------------------------------------------------

#-- create initial vector of infection status as an epimodel object
infected_list <- EpiModel::init.net(i.num=evoparams$initial_infected,
                                    status.rand = FALSE)

#--------------------------------------------------------------

#---  Create list of modules to run for input into epimodel_control_fxn() below

# ***   Note: initialize fxn must always be first and verbose fxn last AND death fxn
# ***   must precede birth fxn (these conditions may change in future)
# ***   treatment_fxn must be before update_vl and update_cd4

evo_module_list<- list(
  "initialize.FUN"     = initialize_module,
  "aging.FUN"          = vital_aging_module,
  "testing.FUN"        = social_testing_diagnosis_module,
  "vacc.FUN"           = social_treatment_vaccination,
  "treatment.FUN"      = social_treatment_module,
  "update_vl.FUN"      = viral_update_gamma,
  "update_cd4.FUN"     = viral_update_cd4_daily,
  "coital_acts.FUN"    = social_coital_acts_module,
  "trans.FUN"          = transmission_main_module,
  "trans_book.FUN"     = transmission_bookkeeping_module,
  "trans_cd4.FUN"      = transmission_cd4_module,
  "deaths.FUN"         = vital_deaths_module,
  "births.FUN"         = vital_births_module,
  "summary.FUN"        = summary_module,
  "resim_nets.FUN"     = EpiModel::resim_nets,
  "verbose.FUN"        = NULL)


#--- call epimodel's control fxn (load evonet modules into epimodel)
evocontrol <- setup_epimodel_control_object(evonet_params = evoparams,
                                            module_list   = evo_module_list)

#--------------------------------------------------------------

if(isTRUE(hyak_par)) {
  evomodel  <- EpiModelHPC::netsim_par(x = estimated_nw, 
                                       param = evoparams, 
                                       init = infected_list, 
                                       control = evocontrol)
} else {
  evomodel  <- EpiModel::netsim(x = estimated_nw,
                                param = evoparams,
                                init = infected_list,
                                control = evocontrol)
}

evomodel$epi <- NULL
evomodel$stats <- NULL
evomodel$control <- NULL

plots_popsumm(evomodel,outpath=outpath,
              name=model_name,nw_stats=TRUE,max_points_rep=100,
              evoparams$popsumm_frequency)

assign(model_name,evomodel)
file_name <- paste(model_name,".RData",sep="")
save(list=model_name,
     file = file.path(outpath,file_name) )
remove(evomodel)

#--------------------------------------------------------------
