#functions:
# 1) output_plot_fxn1
# 2) plot_fxn2
##############################################
#' @export
output_plots_fxn<-function(evomodel,plot_index="all",nsim=1)
{
  #Note: fxn being phased out
  #Due to presence in master scripts, will not be deleted but just
  #not do anything
  #return(invisible(NULL))
  
  evonet_params <- evomodel$param[[nsim]]
  if(evonet_params$n_steps<100){
    cat("\n plotting function requires at least 100 timesteps \n")
    return(invisible(NULL))
  }
  

 if(plot_index=="all" | plot_index %in% c(4:10))
 {
  
   temp_list <- evomodel$pop[nsim][[1]]
   
  pop_df <- as.data.frame(temp_list)
  plot_fxn2(pop_df,evonet_params,type=plot_index,nsim=nsim)
 }

  
}
############################################
#' @export
plot_fxn2 <- function (popdf,params,type=NULL,nsim=1)
{
  
  if(type==4 || type=="all")
  {
    #####################################################
    #plot1
    inf_index <- which(!is.na(popdf$Donors_V))
    if(length(inf_index) >0)
    {
    
    hist(log10(popdf$Donors_V[inf_index]),
         freq=T, xlab="VL at transmission",
         ylab="Frequency",
         main="Plot 4: VL at transmission",
         col="blue")
    }
  }
  ###################################################
  #plot 2
  
  if(type==5 || type=="all")
  {
    logSetPt_breaks=seq(2,7.5,by=.5)
    df<-popdf[which(popdf$Status==-2),]
    if(nrow(df)>0)
    {
      df$totalTimeInf <- df$Time_Death-df$Time_Inf
      df$logSetPt_binned<- round(as.numeric(as.character(
        cut(df$LogSetPoint,
            logSetPt_breaks,
            labels=logSetPt_breaks[-length(logSetPt_breaks)])
      )),2)
      mean_duration_by_VL<- tapply(df$totalTimeInf/365,
                                   df$logSetPt_binned,
                                   mean)
      plot(df$LogSetPoint,
           df$totalTimeInf/365,pch=16,cex=1,col="blue",
           xlim=c(2.5,7),
           xlab="log10 set point viral load", ylab="Duration of infection, Years")
      title("Plot 5")
      lines( as.numeric(names(mean_duration_by_VL)),
             as.numeric(mean_duration_by_VL),
             col='red',type='l',lwd=2)    
      points(as.numeric(names(mean_duration_by_VL)),
             as.numeric(mean_duration_by_VL),col="black",type='p',pch=16,cex=1.3)
      
    }#end of if(nrow(df)>0)
  }
  ###################################################
  
  #plot 5
  if(type==6 || type=="all")
  {
    df<-popdf[which(popdf$Status %in% c(-2,1) & popdf$Generation!=1),]
    temp_sum<-sum(abs(df$Donors_LogSetPoint-df$LogSetPoint))
    if (length(df$LogSetPoint) > 5 & temp_sum >.5) {
      plot(df$Donors_LogSetPoint,df$LogSetPoint,
         main="Plot 6: Heritability of spVL",ylab="Log Set Point",xlab="Donors Log Set Point",
         col="blue",pch=16,cex=1.25)
      h <- (lm(df$LogSetPoint ~ df$Donors_LogSetPoint))
      abline(h, col=2)
      mtext( paste("Regression slope =",
                 format(h$coefficients[2], digits=2)),
           col=2,cex=,outer=F,line=-2,side=3,adj=.2)
    }
  }
  ##########################################################
  #plot 6
  if(type==7 || type=="all")
  {
    df<- popdf[which(popdf$Status %in% c(-2,1) & 
                       popdf$Donors_Total_Time_Inf_At_Trans>0),]
    
    if (length(df$Donors_Total_Time_Inf_At_Trans) > 1) {
      hist(df$Donors_Total_Time_Inf_At_Trans/365,breaks=50,col="blue",
         ylab="No. Infections", xlab="Donor's Total Inf Time at Trans",
         main="")
      title("Plot 7")
    }
  }
  ###############################
  ##########################################################
  #plot 7
  if(type==8 || type=="all")
  {
    df<-popdf[which(popdf$Status %in% c(-2,1) & popdf$Time_Inf>0),]
    if (length(df$LogSetPoint) > 2) {
      plot(df$Time_Inf/365,df$LogSetPoint,
         xlab="Time, in years", ylab="SPVL (log10)",col="blue", cex=1.25)
      lines(lowess(df$Time_Inf/365, df$LogSetPoint, f=0.2), lty=1, lwd=2, col="red")
      title("Plot 8")
    }
  }
  ##########################################################
  ##########################################################
  
  if(type==9 || type=="all")
  {
    aa <- which(is.element(popdf$Status,c(-2,1)))
    if(any(!is.na( popdf$RandomTimeToAIDS[aa])))
    {
      plot(popdf$LogSetPoint[aa], popdf$RandomTimeToAIDS[aa]/365,
           xlab="SPVL", ylab="Time to AIDS, years",
           main="Plot 9: SPVL and time to AIDS", pch=16,col="blue",cex=1.25,
           las=1, xlim=c(2,8), ylim=c(0, 40))
      
      lines(lowess(popdf$LogSetPoint[aa], popdf$RandomTimeToAIDS[aa]/365,
                   f=0.05), col="red", lwd=2)
    }
    

  }
} #end of plot fxn2

