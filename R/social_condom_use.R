
#' @export
social_condom_use <- function(dat,at)
{
  # Assigns condom use (0/1) to each row (analagous to act) in discord_coital_df
  # inputs: param$condom_prob, dat$discord_coital_df
  # output: dat$discord_coital_df$condom
  
  ########################################
  #if no disc. pairs or no sex acts, stop fxn
  if(is.null(dat$discord_coital_df)){return(dat)}
  ###########################################
  
  # Calcuate condom use (0/1) based on "condom_prob"
  ## In MSM model, condom use probabilistically determined by condom_prob
  if(dat$param$model_sex == "msm") {
    dat$discord_coital_df$condom <- rbinom(n = nrow(dat$discord_coital_df),
                                           size = 1,
                                           prob = dat$param$condom_prob)
  }
  
  ## In hetero model, condom use determined by condom_prob and concurrency status of 
  ## male partner
  if(dat$param$model_sex == "hetero") {
    ## 1) Count number of ties in iel and oel lists
    num_ties <- sapply(1:length(dat$nw$iel), function(x) sum(length(dat$nw$iel[[x]]), length(dat$nw$oel[[x]])))
    
    ## 2) Identify the indices of those nodes with more than one tie
    conc_index <- which(num_ties > 1)
    
    ## 3) Identify the ID of each node in conc_index
    id_index <- numeric(0)
    id_index <- sapply(conc_index, function(x) append(id_index, dat$nw$val[[x]]$id))
    
    ## 4) Assign concurrency to corresponding agents in pop dataframe. Note: This value will be updated at each 
    ##    timestep, and so the pop$concurrent vector will not have much meaning at the end of a simulation, as 
    ##    it will only identify those who were concurrent at the last timestep, regardless of concurrency at 
    ##    previous timesteps.
    dat$pop$concurrent[1:length(dat$pop$s)] <- 0
    dat$pop$concurrent[id_index] <- 1
    
    ## 5) Condom use probability is determined by male concurrency, so pull IDs of all males engaging in a 
    ##    discordant sex act
    male_id_vec <- dat$discord_coital_df$sus_id
    male_id_vec[dat$discord_coital_df$inf_sex == "m"] <- dat$discord_coital_df$inf_id[dat$discord_coital_df$inf_sex == "m"]
    
    ## 6) In order to have option of condom use based on female concurrency, also pull IDs of females engaging 
    ##    in discordant sex act
    fem_id_vec <- dat$discord_coital_df$sus_id
    fem_id_vec[dat$discord_coital_df$inf_sex == "f"] <- dat$discord_coital_df$inf_id[dat$discord_coital_df$inf_sex == "f"]
    
    ## 7) Identify which agents engaging in discordant sex act are in concurrent relationships
    dat$discord_coital_df$male_concurrent <- dat$pop$concurrent[male_id_vec]
    dat$discord_coital_df$fem_concurrent <- dat$pop$concurrent[fem_id_vec]
    
    ## 8) Modify condom use probability accordingly.
    ##    Note: Below assumes no interaction in relative risk of condom use when both male and female are concurrent.
    ##    Note: If the product of condom_prob, RR_cond_male_concurrent, and RR_cond_fem_concurrent is greater than 1,
    ##    rbinom will throw an error, as probability must be between 0 and 1, inclusive.
    dat$discord_coital_df$condom <- rbinom(n = nrow(dat$discord_coital_df),
                                           size = 1,
                                           prob = dat$param$condom_prob*((dat$param$RR_cond_male_concurrent^dat$discord_coital_df$male_concurrent)*
                                                                           dat$param$RR_cond_fem_concurrent^dat$discord_coital_df$fem_concurrent))
  }
  
  return(dat)
  
}