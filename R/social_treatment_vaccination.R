#' @export
social_treatment_vaccination <- function(dat, at)
{
  if(at<dat$param$start_vacc_campaign[1]){return(dat)}

  if(at>dat$param$start_vacc_campaign[1] ){
    vacc_duration <- at-dat$pop$vacc_init_time
    vacc_duration_thresh_index <- which(vacc_duration>dat$param$vacc_eff_duration)
    dat$pop$vaccinated[vacc_duration_thresh_index]= 0
  }
  
  if(!is.element(at,dat$param$start_vacc_campaign)){return(dat)}
  #-----------------------------------------
  #below code is for one time vaccination campaign when at==start_tx_campaign
  #-----------------------------------------

  #eligible_patients: eligible for care,not vaccinated
  eligible_index <- which(dat$pop$Status == 0 & 
                          (dat$pop$vaccinated == 0 | is.na(dat$pop$vaccinated)) &
                          dat$pop$eligible_care == 1) 
  
  if(length(eligible_index)==0){return(dat)}
  
  no_eligible <- length(eligible_index)
  no_treated  <- round(no_eligible*dat$param$perc_vaccinated)
  if(no_treated==0){return(dat)}
  
  treated_index <- sample(eligible_index,no_treated)
  
  dat$pop$vaccinated[treated_index] <- 1
  dat$pop$vacc_init_time[treated_index] <- at
  
  return(dat)
}
