#' @export
summary_popsumm<-function(dat,at){

  
  #if summary stats calculations doesn't occur at this timesetp,
  #based on popsumm_frequency value, then leave function
  if( (at%%dat$param$popsumm_frequency!=0) & (at!=1)){return(dat)}
  
  #time_index is a time vector based on current value
  #of "at" and parameter value "popsumm_frequency"
  if(at==1)
    time_index <- 1
  else if(at==dat$param$popsumm_frequency)
    time_index <- (at-dat$param$popsumm_frequency+2):at
  else
    time_index<- (at-dat$param$popsumm_frequency+1):at
  
  #"popsumm_index" is an index for the "popsumm" vectors
  #based on value of "at" and paramter "popsumm_frequency"
  if(at==1)
  popsumm_index <- 1
  else
    if(dat$param$popsumm_frequency==1)
      popsumm_index <- at
    else
      popsumm_index <- (at/dat$param$popsumm_frequency)+1
  
  
  #logical vectors and indices helpful to calculate summary stats  
    
  inf_index     <-  dat$pop$Status == 1
  total_inf     <- length(which(inf_index))
  sus_index     <-  dat$pop$Status == 0
  alive_index   <-  inf_index | sus_index
  treated_index <-  dat$pop$treated == 1 & inf_index
  not_treated_index <-  dat$pop$treated == 0 & inf_index
  treated_undetectable <- treated_index & dat$pop$V<dat$param$vl_undetectable
  treated_agents <- which(treated_index)
  not_treated_agents <- which(not_treated_index)
  new_infections <- is.element(dat$pop$Time_Inf, time_index)
  new_infections_count <- length(which(is.element(dat$pop$Time_Inf, time_index)))
  new_infections_virus_vacc_sens_count <- length(which(is.element(dat$pop$Time_Inf, time_index)&
                                                         dat$pop$virus_sens_vacc==1))
  new_infections_virus_vacc_notsens_count <- length(which(is.element(dat$pop$Time_Inf, time_index)&
                                                         dat$pop$virus_sens_vacc==0))
  
  new_infections_virus_drug_sens_count <- length(which(is.element(dat$pop$Time_Inf, time_index)&
                                                         dat$pop$virus_sens_drug==1))
  new_infections_virus_drug_part_res_count <- length(which(is.element(dat$pop$Time_Inf, time_index)&
                                                         dat$pop$virus_part_res_drug==1))
  new_infections_virus_drug_3_plus_res_count <- length(which(is.element(dat$pop$Time_Inf, time_index)&
                                                         dat$pop$virus_3_plus_drug_muts==1))
  
  donor_time_inf  <- ifelse(new_infections_count>0,
                            dat$pop$Donors_Total_Time_Inf_At_Trans[new_infections],
                            NA)
  donor_acute_count <- ifelse(!is.na(donor_time_inf),
                              length(which(donor_time_inf<=dat$param$t_acute)),
                              NA)
  new_births <- is.element(dat$pop$arrival_time, time_index)
  cd4_aids <- dat$pop$CD4 == 4
  new_diagnoses <- dat$pop$diag_status == 1 &  is.element(dat$pop$diag_time,time_index)
  acute_phase_vec <- (at-dat$pop$Time_Inf)<dat$param$t_acute
  acute_phase <- !is.na(acute_phase_vec) & acute_phase_vec==T
  percent_virus_sensitive <- round(100*(length(which(dat$pop$virus_sens_vacc==1 & inf_index))/length(which(inf_index))))
  
  #deaths
  just_died <- is.element(dat$pop$Time_Death,time_index)
  died_aids <- dat$pop$Status == -2 & just_died
  died_aids_mean_age <- mean(dat$pop$age[died_aids])
  died_non_aids <- dat$pop$Status == -1 & just_died
  died_non_aids_inf <- died_non_aids & !is.na(dat$pop$V)
  died_non_aids_sus <- died_non_aids & is.na(dat$pop$V)
  aged_out <- (dat$pop$age>=dat$param$max_age) & just_died
  
#browser()
  #network statistics
  nw_summary    <-  summary(dat$nw~degree(0:1) + concurrent, at = at)
  number_edges <- unname(summary(dat$nw~edges,at=at)[1,1])
  network_size <- network.size(dat$nw)
  total_nodes   <-  sum(nw_summary[1,1]+nw_summary[1,2]+nw_summary[1,3])
  
  #viral load values
  log10_vl_values  <-  log10(dat$pop$V[which(inf_index)]+dat$param$AbsoluteCut)
  spvl_untreated_values <- (
        dat$pop$LogSetPoint[which(inf_index & not_treated_index)])
  edges_by_agent <- unname(summary(dat$nw ~ sociality(base = 0),at=at)[1,]) #use dat$attr$id for index on dat$pop
  edges_untreated <- edges_by_agent[dat$attr$id %in% not_treated_agents ]
  edges_treated <- edges_by_agent[dat$attr$id %in%  treated_agents]

   #aim3 mutations
   inf_undetect_ix <- (dat$pop$Status==1 & dat$pop$V> dat$param$vl_undetectable)
   no_inf_undect <- length(which(inf_undetect_ix))
  mutations0 <- length(which(inf_undetect_ix & dat$pop$aim3_no_muts==0))
  mutations1 <- length(which(inf_undetect_ix & dat$pop$aim3_no_muts>=1))
  mutations2 <- length(which(inf_undetect_ix & dat$pop$aim3_no_muts>=2))
  mutations3 <- length(which(inf_undetect_ix & dat$pop$aim3_no_muts>=3))
  mutations4 <- length(which(inf_undetect_ix & dat$pop$aim3_no_muts>=4))
  mutations5 <- length(which(inf_undetect_ix & dat$pop$aim3_no_muts>=5))
  
  mutations1exact <- length(which(inf_undetect_ix & dat$pop$aim3_no_muts==1))
  mutations2exact <- length(which(inf_undetect_ix & dat$pop$aim3_no_muts==2))
  mutations3exact <- length(which(inf_undetect_ix & dat$pop$aim3_no_muts==3))
  mutations4exact <- length(which(inf_undetect_ix & dat$pop$aim3_no_muts==4))
  
  mutations3plus_long <- length(which(inf_index & dat$pop$aim3_mutations_long>=3))
  mutations4plus_long <- length(which(inf_index & dat$pop$aim3_mutations_long>=4))
  mutations5_long <- length(which(inf_index & dat$pop$aim3_mutations_long==5))
  
  #coital acts
  
  if(!is.null(dat$discord_coital_df)){
  number_coit_acts <- sum(tapply(dat$discord_coital_df$act_id_couple,
                          dat$discord_coital_df$couple_id,
                          max))
  acts_iev <- length(which(dat$discord_coital_df$iev==1))/2
  percent_iev <-  (acts_iev / number_coit_acts)
  transmission_opps_condom_percent <- (length(which(dat$discord_coital_df$condom==1)) / 
                                 nrow(dat$discord_coital_df) )
   trans_probs_mean <- mean(dat$discord_coital_df$trans_probs)
   
  }else{
    number_coit_acts <- 0
    percent_iev <- NA
    percent_condom <- NA
    trans_probs_mean <- NA
    transmission_opps_condom_percent <- NA
  }
    
  
    #actual calculation of summary stats based on indices and vectors from above
    # and functions for each variable in "popsumm_fxns"
    popsumm_vars=names(dat$popsumm)
  for(ii in 1:length(popsumm_vars)){
    temp_var<-popsumm_vars[ii]
    environment(dat$popsumm_fxns[[ii]])<-environment()
    dat$popsumm[[temp_var]][popsumm_index] <- dat$popsumm_fxns[[ii]]()
   }
  
    #calculation of generic attribute stats
    
    #what percent of alive agents are in each category
    #stat: generic_att_percent_cat_xx (xx=1,..,total number of attributes)
    
    #what percent of alive agents are infected in each category
    #stat: generic_att_percent_inf_cat_xx
    
    #stats for generic attribute values need to be treated separately
    #as the number of attributes may vary between model scenarios
    #note: objects below need to be renamed for clarity 
    temp_length <- length(dat$param$generic_nodal_att_values)
    if(temp_length>1){  
       
       #how many alive agents in each category
       temp_table=table(dat$pop$att1[alive_index])
       #how many alive and infected agents in each category
       temp_table2=table(dat$pop$att1[inf_index])
       #total agents
       sum_temp_table=sum(temp_table)
       #this vector makes sure categories from tables above are
       #arranged in ascending order (necessary if zero agents in a particular 
       #category, which would mean they are missing in tables above
       temp_match=match(names(temp_table),1:temp_length)
       
       
       for(zz in 1:length(temp_match)){
         namevec <- paste("generic_att_percent_cat_",temp_match[zz],sep="")
         dat$popsumm[[namevec]][popsumm_index]=temp_table[zz]/sum_temp_table
       }
       for(zz in 1:temp_length){
         namevec2 <- paste("generic_att_percent_inf_cat_",zz,sep="")
         
         ix1<- which(names(temp_table)==zz)
         ix2<- which(names(temp_table2)==zz)
         if(length(ix2)>0){
           val<- temp_table2[ix2]/temp_table[ix1]  
         }else{val<-0}
         dat$popsumm[[namevec2]][popsumm_index] <- val
       }
    }
  # end of calculating summary stats for generic attributes    
return(dat)
}
