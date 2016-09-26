#' @export
social_treatment_module <- function(dat, at)
{
  
  #Description:
  #Assigns treatment status (0/1) to eligible patients based on:
  # is there treatment camagin?
  # is vl > vl threshold for treatment
  # is infection duration > minimum duration for start of treatment
  # is cd4 < minimum threshold
  # is there a postivie diagnosis
  #inputs: tx_type=="VL","CD4","time","vl_and_cd4","vl_and_time","vl_and_cd4_and_time"
          #param$start_treatment_campaign
           #pop$diag_status
           #pop$treated
           #pop$Status
           #pop$CD4
           #pop$Time_Inf
           #param$t_acute
           #param$min_inf_time_for_treat          
  #outputs: pop$treated
            #pop$tx_init_time
  dat$treatment_index <- NULL
  eligible_patients_criteria <- NULL
  
  if(at < dat$param$start_treatment_campaign){return(dat)}
  
  infected <-  which(dat$pop$Status==1)
  if(length(infected)==0){return(dat)}
  
  #eligible_patients: infected,diagnosed,eligible for care,not treated
  eligible_patients <- 
      which(dat$pop$Status == 1 & dat$pop$treated == 0 &
            dat$pop$diag_status == 1 & dat$pop$eligible_care == 1 &
            dat$pop$eligible_ART == 1 &
            ((at - dat$pop$Time_Inf) > dat$param$t_acute)  ) 
  
  if(length(eligible_patients)==0){return(dat)}
  
  #eligible_patients_criteria: agents in eligible_patients not in acute stage and meet various
  #treatment criteria (time,vl,cd4) set in input_parameters_primary
  #Viral load: if VL higher than treatment threshold
  #CD4 count: if CD4 meets treatment threshold
  #Time: if time agent has been infected is greater than the minimun time for treatment
  
  if(dat$param$tx_type=="VL"){
    
    eligible_patients_criteria <- 
      which(dat$pop$V[eligible_patients] > dat$param$treatment_threshold)
  
  }
   
  if(dat$param$tx_type=="CD4") {
    
    eligible_patients_criteria <- 
      which(is.element(dat$pop$CD4[eligible_patients], dat$param$cd4_treatment_threshold))
  }
  if(dat$param$tx_type=="time"){
    
    eligible_patients_criteria <- 
      which((at-dat$pop$Time_Inf[eligible_patients]) > dat$param$min_inf_time_for_treat)
  }
  if(dat$param$tx_type=="vl_and_cd4"){
    
    eligible_patients_criteria <- 
        which(dat$pop$V[eligible_patients] > dat$param$treatment_threshold &
        is.element(dat$pop$CD4[eligible_patients], dat$param$cd4_treatment_threshold) )
  }
  if(dat$param$tx_type=="vl_and_time"){
    
    eligible_patients_criteria <- 
        which(dat$pop$V[eligible_patients] > dat$param$treatment_threshold &
        (at-dat$pop$Time_Inf[eligible_patients]) > dat$param$min_inf_time_for_treat )
  }
  if(dat$param$tx_type=="cd4_and_time"){
    
    eligible_patients_criteria <- 
         which( is.element(dat$pop$CD4[eligible_patients], dat$param$cd4_treatment_threshold) &
         (at-dat$pop$Time_Inf[eligible_patients]) > dat$param$min_inf_time_for_treat)
  }
  
  if(dat$param$tx_type=="vl_and_cd4_and_time"){
    
    eligible_patients_criteria <- 
         which(dat$pop$V[eligible_patients] > dat$param$treatment_threshold &
         is.element(dat$pop$CD4[eligible_patients], dat$param$cd4_treatment_threshold) &
         (at-dat$pop$Time_Inf[eligible_patients]) > dat$param$min_inf_time_for_treat )
  }
  
  if(length(eligible_patients_criteria)==0){return(dat)}

  #subset of eligible agents (eligible_patients) that meet specified tx criteria
  eligible_patients_treated <- eligible_patients[eligible_patients_criteria]
  
  dat$pop$treated[eligible_patients_treated] <- 1
  dat$pop$tx_init_time[eligible_patients_treated] <- at
  dat$treatment_index <- eligible_patients_treated
  
 return(dat)
}
