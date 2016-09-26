#' @export
input_parameters_age_distribution <- function(data_name="usa_men_18_to_100",min_age,max_age)
{
  #description: returns initial age distribution for population
  
  age_dist_data_list <- list()
  ########## Alternative AMSR tables that include elderly persons. #########################
  # Data obtained from the CDC "WONDER" webpage (wonder.cdc.gov) for USA men from 1999-2003
  # downloaded on 8/25/15.  The CDC data only apply to people 85 and under.  To fill in
  # the rest, I (John) made various extrapolations and approximations, as follows... 
  #  (1) In the absence of data, I assumed zero people over 85 years at the beginning of the
  #      simulation.  The model does, however, allow people to age-in to ages 86-100 over time.
  #  (2) Death rate data for those over 86 were obtained from Society of Actuaries' "Social
  #      Security" data set (pretty closely matches the CDC estimate for those 85!).  Obtained from
  #      https://www.soa.org/research/software-tools/research-simple-life-calculator.aspx 
  #      downloaded on 8/25/15.
  #  (3) To keep this from getting totally out of hand, I assumed a 0% chance of living past 100. 
  #
  # While these approximations and extrapolations are imperfect, I figure that the advanced
  # elderly are rare enough (and so inactive sexually) that any imperfections will have
  # little-to-no effect on our conclusions.
  #
  #  Note: Eldery persons can easily be excluded from the simulations by setting max_age to 
  #        something less than 100 (current default, as of 8/25/15, is max_age = 55 years)
  #  
  age_dist_data_list$"usa_men_18_to_100"<- list(
    age_dist = 
      c(0.0205, 0.0206, 0.0204, 0.0202, 0.0200, 0.0196, 0.0194, 0.0192,
        0.0189, 0.0189, 0.0188, 0.0188, 0.0191, 0.0187, 0.0189, 0.0188,
        0.0189, 0.0192, 0.0193, 0.0194, 0.0196, 0.0199, 0.0201, 0.0203,
        0.0204, 0.0204, 0.0203, 0.0202, 0.0202, 0.0201, 0.0199, 0.0197,
        0.0194, 0.0191, 0.0188, 0.0181, 0.0177, 0.0170, 0.0165, 0.0158,
        0.0152, 0.0145, 0.0138, 0.0132, 0.0126, 0.0120, 0.0114, 0.0109,
        0.0102, 0.0096, 0.0092, 0.0088, 0.0084, 0.0080, 0.0076, 0.0073,
        0.0069, 0.0066, 0.0062, 0.0059, 0.0056, 0.0052, 0.0048, 0.0044,
        0.0040, 0.0036, 0.0032, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,
        0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,
        0.0001)) 

  
  age_index        <- min_age : (max_age-1)
  age_dist_index   <- age_index - min_age + 1
  initial_age_dist <- age_dist_data_list[[data_name]]$age_dist[age_dist_index]
  final_age_dist   <- initial_age_dist/sum(initial_age_dist)
  return(final_age_dist)
}

