#' @export
new_additions_fxn_both <- function(input_list,dat,index,type=c("births","initial"),at)
{
  #description:
  #main argument: type="births" or type="initial"
  #fills in agent attribute "pop" list; some variales get same value for
  #initial population and for new additions ("births") while others get
  #different values based on whether its for the initial population or 
  #for new additions. An example is age, ages for initial popn can take a
  #number of distributions/values but ages for births are always the min. age.
  #input: dat$dat$param,dat$pop
  #output: dat$pop
  
  #"pop" attributes with different values/distributions for initial pop vs births:
  #age,last_neg_test,arrival_time,role, att1 (generic attribute),

  
  total_new <- length(index)
  #------------------------------
  #functions for these variables are same initial population and for new births
  
  if(dat$param$model_sex=="msm"){
    input_list$sex[index] <-1
  }else{input_list$sex[index]<-rbinom(total_new,1, 0.5 )}
    
  input_list$s[index]<- rep(dat$param$prog_rate,total_new)

  input_list$id[index] <-  index   
  
  input_list$Status[index] <-  rep(0,total_new)   
  
  input_list$NumRecipients[index] <-  rep(0,total_new)
  
  input_list$sti_status[index] <- rbinom(total_new,1,dat$param$sti_prob)
  
  input_list$circum[index] <- rbinom(total_new,1,dat$param$circum_prob)
  
  input_list$total_acts[index] <- rep(0,total_new)
  
  input_list$Adherence[index] <- runif(total_new, dat$param$min_adherence,dat$param$max_adherence)
  
  input_list$tx_schedule[index] <- sample(names(dat$param$tx_schedule_props),
                                          size=total_new, replace=T,
                                          prob=dat$param$tx_schedule_props)
  
  input_list$CD4_nadir[index] <- 1
  
  input_list$CD4[index] <- 1   # categorical value 1: CD4 > 500,..., 4: CD4: <200
  
  input_list$CD4tot[index] <- 1000  # CD4 T-cell blood count before redistribution
  
  input_list$CD4count[index] <- 1000 # CD4 T-cell blood count w/ redistribution to LN (when V > 0)
  

  
  #-----------------------------
  #these variables need different functions for initial population and births
  if(type=="initial")
  {
    if(!is.logical(dat$param$generic_nodal_att_values))
    {
      input_list$att1[index] <- network::get.vertex.attribute(dat$nw, "att1",unlist = TRUE)
    }
    
    if(!is.logical(dat$param$role_props))
    {
      input_list$role[index] <- network::get.vertex.attribute(dat$nw, "role",unlist = TRUE) 
      
      temp <- index[input_list$role[index]=="I"]
      input_list$insert_quotient[temp] <- 1
      
      temp <- index[input_list$role[index]=="R"]
      input_list$insert_quotient[temp] <- 0
      
      temp <- index[input_list$role[index]=="V"]
      input_list$insert_quotient[temp] <- runif(length(temp))
      
    }else{
      input_list$role[index] <- "V"
      input_list$insert_quotient[index] <- runif(total_new)
    }
    
    
    
    input_list$age <- vital_initial_age_dist( 
                               type      = dat$param$init_age_dist_method,
                               age.range = dat$param$min_age : (dat$param$max_age-1),
                               popsize   = dat$param$initial_pop,
                               user.dist = dat$param$male_age_dist)   
    
    input_list$last_neg_test[index] = sample( - dat$param$mean_test_interval : 0, 
                                                   total_new,
                                                   replace = TRUE )
  }
  
  if(type=="births")
  {
    if(!is.logical(dat$param$generic_nodal_att_values))
    {
      input_list$att1[index] <- sample(dat$param$generic_nodal_att_values ,
                          total_new, replace=TRUE,
                          prob=dat$param$generic_nodal_att_values_props_births)
    }
    if(!is.logical(dat$param$role_props))
    {
      input_list$role[index] <- sample(names(dat$param$role_props) ,
                                            total_new,
                                            replace=TRUE,
                                            prob = dat$param$role_props)
      
      temp <- index[input_list$role[index]=="I"]
      input_list$insert_quotient[temp] <- 1
      
      temp <- index[input_list$role[index]=="R"]
      input_list$insert_quotient[temp] <- 0
      
      temp <- index[input_list$role[index]=="V"]
      input_list$insert_quotient[temp] <- runif(length(temp))
      
    }else{
      input_list$role[index] <- "V"
      input_list$insert_quotient[index] <- runif(total_new)
    }
    
    
    temp_sample_times = (at - dat$param$mean_test_interval): at
    input_list$last_neg_test[index] <- sample(temp_sample_times,  
                                                   total_new,
                                                   replace=TRUE )
    
    temp_sample_times = (at - dat$param$mean_resist_test_interval): at
    input_list$last_neg_resist_test[index] <- sample(temp_sample_times,  
                                                     total_new,
                                                     replace=TRUE )
    input_list$age[index] <- dat$param$min_age
    input_list$arrival_time[index] <- at
    
  }
  #######################################
  return(input_list)
}

