###################################################################
# Example script to simulate an HIV epidemic over 40 years in     #
# which a partially-effective vaccine is introduced at year 22.   #
# 75% of circulating HIV are sensitive to the vaccine effect. The #
# vaccine reduces risk of transmission by 75% given exposure to a #
# sensitive virus, and has no effect on the risk of transmission  #
# given exposure to a vaccine-resistant virus. Vaccine coverage   #
# is 90%.                                                         #
#                                                                 #
# NOTE: User must input a file path for output files on line 22.  #
###################################################################

if(!require(devtools)) { install.packages("devtools") }
library(devtools)
devtools::install_github("statnet/tergmLite")
devtools::install_github("statnet/EpiModel", ref = "fast_edgelist")

library(evonet)

# Load default parameters
primary_parameters  <- input_parameters_primary()
cd4_data            <- input_parameters_cd4_data()

# Combine individual parameters into single list
evoparams <- c(primary_parameters, cd4_data)

evoparams$output_path = "H:/evonet_runs/Vaccines"

# Override selected defaul parameters
evoparams$nsims            = 3

# -- ART parameters --------------------------------------------- #
# ART scenario 1
evoparams$tx_type = "cd4_and_time_dist"
evoparams$cd4_trt_guidelines_chgs = list(3:4, 1:4)
evoparams$start_treatment_campaign = c(7*365, 12*365)
evoparams$mean_time_tx = 365*3
evoparams$sd_time_tx = 365
evoparams$prob_eligible_ART = 0.40

# -- Vaccine parameters ----------------------------------------- #
evoparams$start_vacc_campaign = (22*365):(40*365)
evoparams$perc_virus_vaccine_sens <- 0.75
evoparams$trans_RR_vaccine <- 0.25
evoparams$perc_vaccinated <- (0.90/(1-0.90))/(3*365)
evoparams$vacc_eff_duration = 365*3

# -- Misc. parameters to override ------------------------------- #
evoparams$initial_pop      = 500
evoparams$initial_infected = 100
evoparams$target_stats     = .35*evoparams$initial_pop

evoparams$n_steps          = 365*40
evoparams$popsumm_frequency = 45
evoparams$poisson_birth_lambda = 0.0137*(evoparams$initial_pop/100)
evoparams$trans_RR_age = 1.0
#evoparams$fast_edgelist = T
evoparams$vl_peak_agent_flag = T

# Add parameters that are functions of other input parameters
evoparams  <- input_parameters_derived(evoparams)
evoparams$age_dist <- seq(50, 10, -10/9)/1110

# Convert raw parameter list into EpiModel object
evoparams <- do.call(EpiModel::param.net,evoparams)

# Create and initialize network
nw <- setup_initialize_network(evoparams)

# Run QAQC on input parameters
input_parameters_qaqc(evoparams)

# Estimate initial network
netest_arg_list <- list(
  nw            =  nw,
  formation     =  as.formula(evoparams$nw_form_terms),
  target.stats  =  evoparams$target_stats,
  coef.form     =  evoparams$nw_coef_form,
  constraints   =  as.formula(evoparams$nw_constraints),
  verbose       =  FALSE,
  coef.diss     =  dissolution_coefs(dissolution =  ~offset(edges),
                                     duration    =  evoparams$relation_dur,
                                     d.rate      =  3e-05))

estimated_nw <- do.call(EpiModel::netest, netest_arg_list)

model_name = "Vax1_ART1_sens0.75_RR0.25_cov0.90_EXAMPLE"

# Create initial vector of infection status as an EpiModel object
infected_list <- EpiModel::init.net(i.num=evoparams$initial_infected,
                                    status.rand = FALSE)

# Create list of modules to run for input into epimodel_control_fxn() below
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

# Call EpiModel's control fxn (load evonet modules into EpiModel)
evocontrol <- setup_epimodel_control_object(evonet_params = evoparams,
                                            module_list   = evo_module_list)

# Run simulation
runtime <- system.time({
  evomodel  <- EpiModel::netsim(x = estimated_nw,
                                param = evoparams,
                                init = infected_list,
                                control = evocontrol)
})

print(runtime)

evomodel$epi <- NULL
evomodel$stats <- NULL
evomodel$control <- NULL

plots_popsumm(evomodel, outpath = evoparams$output_path,
              name = model_name, nw_stats = TRUE, max_points_rep = 100,
              evoparams$popsumm_frequency)

assign(model_name, evomodel)

file_name <- paste(model_name,".RData",sep="")

save(list = model_name,
     file = file.path(evoparams$output_path, file_name))

remove(evomodel)

#--------------------------------------------------------------
