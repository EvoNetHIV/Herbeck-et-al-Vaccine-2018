#' @export

social_adherence<-function(dat,at){
  
  #note: adherence_type==1, is default (random), set in vital_new_additions
  
  #type2 adherence: cyclic, asynchronous start
  adh2_ix <- which(dat$pop$adherence_type==2)
  if(length(adh2_ix)>0){

    start <- dat$pop$adherence_start[adh2_ix]
    days_high <- dat$param$adherence_days_high
    days_low <- dat$param$adherence_days_low

    day_within_cycle <- (at + start - 1) %% (days_high + days_low)
    ix1 <- which(day_within_cycle<days_high)
    ix2 <- which(day_within_cycle>=days_high)
    prob<-numeric(length(adh2_ix))
    if(length(ix1)>0){prob[ix1]<- dat$param$aherence_days_high_prob}
    if(length(ix2)>0){prob[ix2]<- dat$param$aherence_days_low_prob}
    
    dat$pop$Adherence1[adh2_ix]=prob
    dat$pop$Adherence2[adh2_ix]=prob
    dat$pop$Adherence3[adh2_ix]=prob
    dat$pop$Adherence4[adh2_ix]=prob
    
  }
  #end type2 adherence
  return(dat)
}