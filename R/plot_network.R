#' @export
plot_network_fxn <-function(dat,at)
{
  #exit if no network plotting    
  if(dat$param$plot_nw != TRUE){return(dat)}
  #exit if not first timestep or not timestep to plot network
  if( (at != 2) & ((at %% dat$param$network_print_frequency) != 0)){
    return(dat)
  }
  #exit if "save_network==T" (except if first timestep),
  #as plotting doesn't work well on network with complete history
  if(dat$param$save_network== TRUE & at>2){return(dat)} 
    
  #color from low to high values (spvl, connectivity)
  color_vec <- c("white","orange","yellow","red","purple","black")
  
  # Create two plots if we have generic attribute groups
  if(!is.logical(dat$param$generic_nodal_att_values))
     par(mfrow=c(1,2)) 
  else
     par(mfrow=c(1,1)) 
  
  
  
  #get ids of agents on network that are infected
  inf_index <- which(dat$attr$status_evo == 1 & dat$attr$active == 1)  # No longer used now that I try to colorize the nodes
  #position on network
  inf_index_id <- dat$attr$id[inf_index]
  
  spvl_breaks=c(0,1e3,1e4,1e5,1e6,1e10)
  spvl_label <- as.numeric(cut(dat$pop$V[inf_index_id],breaks=spvl_breaks,
                           labels=1:5))
  #sex-specific symbols
  pch_vec <- rep(50,length(dat$attr$status_evo))
  pch_vec[dat$attr$sex=="f"]=3
  
  color_ix <- rep("white",length(dat$attr$status_evo))
  color_ix[inf_index] <- color_vec[spvl_label]
  
  if(length(which(dat$attr$active ==1))<=100)
     plot(dat$nw,vertex.col=color_ix,label="id",label.cex=.85,vertex.sides=pch_vec)  
  else
     plot(dat$nw,vertex.col=color_ix,vertex.sides=pch_vec)
  
    mtext(paste("network at day: ",at," (dark = high VL)"),side=3)
  
    
  # Plot # 1 -- Nodes colored by connection group
  if(!is.logical(dat$param$generic_nodal_att_values)){
    
    inf_index <- which(dat$attr$status_evo >=0 & dat$attr$active ==1) 
    color_ix <- color_vec[dat$pop$att1[inf_index]] 
  
  if(length(which(dat$attr$active ==1))<=100)
    plot(dat$nw, vertex.col = color_ix,label="id",label.cex=.85,vertex.sides=pch_vec)  
  else
    plot(dat$nw, vertex.col = color_ix,vertex.sides=pch_vec)
    
    mtext(paste("network at day ",at," (black to yellow: activity grps 1-5)"),side=3)
 }

    
  return(dat)
}
