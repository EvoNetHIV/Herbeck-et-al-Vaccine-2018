#' @export
social_discord_edgelist <- function (nw, status_vec, at) 
{ 
  
  #######################################################
  #-- returns discordant edgelist from network
  #-- called in social_discord_edgelist_df
  #-- individual IDs refer tovnetwork position
  #-- called in social_discord_edgelist_df()
  #-- modification of Sam's EpiModel discord_edgelist()
  
  #input: network, status of agents (inf/sus)
  #output: 2 column matrix of discordant pairs (by agent ids)
  
  #######################################################
  
  status  <- status_vec
  edge_list  <- get.dyads.active(nw, at = at)
  if(length(edge_list)==0){return(NULL)}
  edge_list  <- edge_list[sample(1:nrow(edge_list)), , drop = FALSE]
  stat    <- matrix(status[edge_list], ncol = 2)
  sums    <- rowSums(stat)
  disc    <- which(sums==1)
  
  if(length(disc)==0){return(NULL)}
  
  col2inf <- which( stat[,1] / stat[,2] == 0)
  col1inf <- disc[ which(!is.element(disc,col2inf)) ]
  del_mat <- rbind( edge_list[col1inf,], edge_list[col2inf,2:1] )
  del     <- data.frame(timestep = at, 
                        inf_id = del_mat[,1], 
                        sus_id = del_mat[,2])
  
  return(del)
}
