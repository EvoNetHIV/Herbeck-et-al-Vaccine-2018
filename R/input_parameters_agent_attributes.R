
#' @export
input_parameters_agent_attributes <-function(){
  
  #Description:
  # Vector of the individual agent attributes that populate the “pop” list
  
  list( popVariables= c(
    
#--  Virulence model parameters ------------------------------ #

    "s",                               "Status",           
    "NumRecipients",                   "Treated"  , 
    "ViralContribToLogSP0" ,
    "EnvirContribToLogSP0",            "LogSetPoint",           
    "SetPoint",                        "d_acute",  
    "Generation",                      "RandomTimeToAIDS",
    "Time_Inf",                        "V",
     "r0",
    "vl_phase2_trans",                 "rate_phase2",
    "PPP",                              "vl_peak_agent",
    "vl_at_test",                      "cd4_at_test",
    "Donors_V",                        "Donors_treated",  "Donors_treated_2nd_line",
    "Donors_ViralContribToLogSP0",     "Donors_EnvirContribToLogSP0",
    "Donors_Total_Time_Inf_At_Trans",  "Donors_Generation",
    "Donors_Index",                    "Donors_age",                   
    "Donors_LogSetPoint",              "Donors_SetPoint",
    "Donors_d_acute",                  "Time_Death",
    "Donors_diag_status",    
    
    "V_vec",                           "I_vec",
    "M_vec",                           "L_vec",
    "K",                               "Imm_Trig", 
    "CD4count",                        "CD4tot",
    "ChronPhase",                      "OnDrug",
    "Adherence1","Adherence2",         "Adherence3","Adherence4",
    "Drug1", "Drug2",                  "Drug3", "Drug4", 
    "aim3_no_muts",                    "adherence_start",
    "adherence_type",                  "CD4",                                             
    "CD4_time",                        "CD4_initial_value",
    "CD4_treatment_delay_index",       "spvl_cat",
    "CD4_time_death",                  'CD4_time_cat1',
    "CD4_time_cat2",                   "CD4_time_cat3",
    "CD4_time_cat4",                    "CD4_nadir",
    "virus_sens_vacc",                  "virus_sens_drug",
    "virus_part_res_drug",              "virus_3_plus_drug_muts",
    "virus_3_plus_drug_muts",             "Aim3RoundingErrors",               
     "aim3_mutations_long",
    
# -- testing for and treatment of drug resistant viruses (aim 3) --- #    
    "eligible_2nd_line_ART",
    "treated_2nd_line",  
    "diag_resist_status",
    "diag_resist_time", 
    "last_neg_resist_test",
    "time_init_2nd_line",
    
#--  Vital dynamics /social/treatment ------------------------------ #
    "vaccinated",                      "vacc_init_time",           
    "age",                             "arrival_time",
    "last_neg_test",                   "diag_status",
    "diag_time",                       "disclosure_status",
    "id",                              "eligible_care",
    "att1",                            "AI_prob",
    "treated",                         "tx_init_time",
    "circum",                           "tx_schedule",               
    "sti_status",                       "sex",
    "viral_lineage",                    "viral_lineage_init_SPVL",
    "total_acts",                        "role",
    "min_time_tx")
    
    ##############################
  ) #end of  list
  
}
