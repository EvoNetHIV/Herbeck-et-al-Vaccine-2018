#' @export
social_treatment_vaccination <- function(dat, at) {
  
  if(at < dat$param$start_vacc_campaign[1]) {return(dat)}

  if(at > dat$param$start_vacc_campaign[1]) {
    vaccinated <- which(dat$pop$vaccinated == 1)
    dat$pop$vaccinated[vaccinated] <- rbinom(length(vaccinated), 1, 1 - (1/dat$param$vacc_eff_duration))
  }
  
  if(!is.element(at,dat$param$start_vacc_campaign)) {return(dat)}

  # Eligible_patients: eligible for care, not vaccinated, not infected
  eligible_index <- which(dat$pop$Status == 0 & 
                          (dat$pop$vaccinated == 0 | is.na(dat$pop$vaccinated)) &
                          dat$pop$eligible_care == 1) 
  
  if(length(eligible_index) == 0) {return(dat)}
  
  no_treated <- sum(rbinom(length(eligible_index), 1, dat$param$perc_vaccinated))
  if(no_treated == 0) {return(dat)}
  
  treated_index <- sample(eligible_index, no_treated)
  
  dat$pop$vaccinated[treated_index] <- 1
  dat$pop$vacc_init_time[treated_index] <- at
  
  return(dat)
}
