#' @export
social_2nd_line_treatment_module <- function(dat, at)
{
  
  # Description:
  # Assigns second line treatment status (0/1) to eligible patients based on being 
  #  infected, diagnosed, eligible for care, eligible for ART, and (most importantly)
  #  being diagnosed as being infected with a drug resistant virus 
  #inputs: 
           #pop$diag_resistance_status
  #outputs: pop$treated_2nd_line
            #pop$tx_2nd_line_init_time
  dat$treatment_index <- NULL
  eligible_patients_criteria <- NULL
  
  if(at < dat$param$start_treatment_campaign){return(dat)}
  
  infected <-  which(dat$pop$Status==1)
  if(length(infected)==0){return(dat)}
  
  #eligible_patients: infected, diagnosed with HIV, eligible for care, eligible for 2nd line ART,
  #                   and diagnosed with drug resistant virus
  eligible_patients <- 
      which(dat$pop$Status == 1 & 
            dat$pop$diag_status == 1 & dat$pop$eligible_care == 1 &
            dat$pop$eligible_2nd_line_ART == 1 & 
            dat$pop$diag_resist_status == 1) 
  
  if(length(eligible_patients)==0){return(dat)}
    
  dat$pop$treated_2nd_line[eligible_patients] <- 1
  dat$pop$time_init_second_line[eligible_patients] <- at
  # dat$treatment_2nd_line_index <- eligible_patients
  
 return(dat)
}
