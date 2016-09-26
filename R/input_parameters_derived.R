
#' @export
input_parameters_derived  <- function(initial_param)
{
  #Description:
  # Basic parameters calculated from parameter values in input_parameters_primary().  
  # Make changes to parameter values in input_parameters_primary().
 
  #1st, some basic qaqc
  initial_param$n_steps <- round(abs(initial_param$n_steps))
  initial_param$initial_pop <- round(abs(initial_param$initial_pop))
  initial_param$initial_infected  <- round(abs(initial_param$initial_infected))
  
derived_params <- list( 
  poisson_birth_lambda     = (initial_param$initial_pop/100)* initial_param$poisson_birth_lambda_base,
  h     = sqrt(initial_param$Heritability),
  r0    = (log(initial_param$vl_peak_acute / initial_param$V0) /
              initial_param$t_peak),
  pop_growth_rate_timestep  = utilities_annual_prob_conversion(
                                    initial_param$pop_growth_rate_annual,
                                    365),
  male_age_dist = input_parameters_age_distribution(
              data_name="usa_men_18_to_100",
              initial_param$min_age,initial_param$max_age),
  
  mort_per_timestep_male = input_parameters_asmr(
              data_name=initial_param$asmr_data_male ,
              initial_param$min_age,initial_param$max_age), 
  
  mort_per_timestep_female = input_parameters_asmr(
             data_name=initial_param$asmr_data_female,
             initial_param$min_age,initial_param$max_age)) 

final_params <- c(initial_param,derived_params)

 return(final_params)
}
