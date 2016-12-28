#' @export
setup_initialize_network <- function(params){
#create initial network, set nodal attributes "sex","att1","role"
#input: initial parameter list (evoparams)
#output: initial network ("nw")
  
# initialize network
nw <- network::network.initialize(params$initial_pop,
                                  directed = FALSE)

#"id" mainly for qaqc testing
id_vec <- 1:params$initial_pop
network::set.vertex.attribute(x = nw, attr = "id",
                              value = id_vec)


if(params$model_sex!="msm"){
  sex_vec <- sample(c("m","f"),params$initial_pop,prob=c(0.5,0.5),replace=T)
  network::set.vertex.attribute(x = nw, attr = "sex", value =sex_vec)
} else {
  network::set.vertex.attribute(x=nw, attr="sex", value="m")
}
#---------------------------------


#age
age_vec=rep(NA_real_,params$initial_pop)
if(params$model_sex=="msm"){
index_male <- 1:params$initial_pop
age_vec[index_male] <- vital_initial_age_dist( 
  age.range = params$min_age : (params$max_age-1),
  popsize   = length(index_male),
  age_dist = params$male_age_dist)
}else{
  index_male <- which(sex_vec=="m")
  index_female <- which(sex_vec=="f")
  age_vec[index_male] <- vital_initial_age_dist( 
    age.range = params$min_age : (params$max_age-1),
    popsize   = length(index_male),
    age_dist = params$male_age_dist)
  age_vec[index_female] <- vital_initial_age_dist( 
    age.range = params$min_age : (params$max_age-1),
    popsize   = length(index_female),
    age_dist = params$female_age_dist)
}

#temp qaqc
if(any(is.na(age_vec))){browser()}

network::set.vertex.attribute(x=nw, attr="age", value=age_vec)




#--------------------------------
#set generic atts on nw if in use
if(!is.logical(params$generic_nodal_att_values)) 
  nw<- social_set_generic_attribute_on_nw(params, nw)

#set msm roles on nw if in use (default: all versatile)
if(!is.logical(params$role_props) && params$model_sex=="msm")
  nw<-social_set_msm_role_on_nw(params,nw)

return(nw)
}