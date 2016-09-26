#' @export
plots_initial_diagnostic<-function(x,param, init,control,s=1)
{
  
  dat <- initialize_module(x,param, init, control, s)

  hist(dat$pop$LogSetPoint,col="blue",freq=TRUE,main="Initial log setpoint",
       ylab="frequency",xlab="log setpoint",breaks=50)
  
  hist(dat$pop$d_acute,col="blue",freq=TRUE,main="Initial d_acute",
       ylab="frequency",xlab="d_acute",breaks=50)
  
  if (param$aids_death_model == "Gamma_Death") {
     hist(dat$pop$RandomTimeToAIDS/365,col="blue",freq=TRUE,main="time to aids (gamma distributed)",
       ylab="frequency",xlab="time (years)",breaks=50)
  }
  
  hist(log10(dat$pop$V),
       col="blue",
       freq=TRUE,
       main="initial viral load for all initial infections",
       ylab="frequency",xlab="log10 viral load",
       breaks=50)
  
  barplot(param$ASMR,col="blue",ylab="ASMR",xlab="ages: 18-70")
  title("ASMR: ages 18 - 70")
  
}
