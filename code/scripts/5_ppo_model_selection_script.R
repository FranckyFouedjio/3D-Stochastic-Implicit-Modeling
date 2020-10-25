T1=Sys.time()
##########################################################################################################################################
#                                                                                                                                        #
#                                              SYNTHETIC CASE STUDY: PPO MODEL SELECTION                                                 #
#                                                                                                                                        #
##########################################################################################################################################

## WORKING DIRECTORY SETTING
setwd("/media/fouedjio/095e0241-4724-4a1f-84f8-9ddda0df2da9/fouedjio/3d_stochastic_implicit_modeling/")



## REQUIRING PACKAGES LOADING
library(RGeostats)
library(fields)
library(rgl)
library(plot3D)
library(e1071)
library(png)
library(plot3Drgl)
library(misc3d)
library(data.table)
library(deldir)
library("RColorBrewer")
library(profvis)
library(DescTools)
library(foreach)
library(doParallel)
library(pracma)
library(reticulate)
library(Rvcg)
library(geometry)
library(nat)
library(nloptr)
library(Matrix)
library(compositions)
library(pracma)
library(mco)
library(FNN)
library(data.table)
library(feather)

## SOURCE CODE LOADING
source("./code/functions/global_ppo_computing_function.R")
source("./code/functions/dichotomization_function.R")
source("./code/functions/pseudo_signed_distance_function.R")
source("./code/functions/truncation_rule_domaining_function.R")

##CREATING SOME FOLDERS
if(!file.exists("./outputs/figures/5_ppo_model_selection")){dir.create("./outputs/figures/5_ppo_model_selection")}
if(!file.exists("./outputs/figures/5_ppo_model_selection/variogram_residuals")){dir.create("./outputs/figures/5_ppo_model_selection/variogram_residuals")}
if(!file.exists("./outputs/figures/5_ppo_model_selection/probability_maps")){dir.create("./outputs/figures/5_ppo_model_selection/probability_maps")}
if(!file.exists("./outputs/figures/5_ppo_model_selection/probability_maps/1")){dir.create("./outputs/figures/5_ppo_model_selection/probability_maps/1")}
if(!file.exists("./outputs/figures/5_ppo_model_selection/probability_maps/2")){dir.create("./outputs/figures/5_ppo_model_selection/probability_maps/2")}
if(!file.exists("./outputs/figures/5_ppo_model_selection/probability_maps/3")){dir.create("./outputs/figures/5_ppo_model_selection/probability_maps/3")}
if(!file.exists("./outputs/figures/5_ppo_model_selection/probability_maps/4")){dir.create("./outputs/figures/5_ppo_model_selection/probability_maps/4")}


## DATA LOADING
data_points=as.data.frame(fread("./inputs/data/drillhole_data.csv"))
contact_points=as.data.frame(fread("./inputs/data/contact_points.csv"))
grid_points=as.data.frame(fread("./inputs/data/grid_points.csv"))
contact_points_cat=readRDS('./inputs/data/contact_points_cat.RDS')
data_points_binary=as.data.frame(fread("./outputs/data/drillhole_data_binary.csv",header=TRUE))
psd_data=as.data.frame(fread("./outputs/data/psd_data.csv"))
lithology_names=1:4
nsdf=length(lithology_names)



closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(length(lithology_names))
registerDoParallel(Cl)

full_grid_points=foreach(i=1:nsdf) %dopar%
{
  library(data.table)
  rbindlist(list(grid_points[,1:3],data_points[,2:4],contact_points_cat[[i]][,2:4]))
}
closeAllConnections()
## TARGET GRID DEFINITION
extent_x=range(data_points$X)
extent_y=range(data_points$Y)
extent_z=range(data_points$Z)


cell_size_x=7
cell_size_y=7
cell_size_z=5

xx=seq(min(extent_x),max(extent_x),by=cell_size_x)
yy=seq(min(extent_y),max(extent_y),by=cell_size_y)
zz=seq(min(extent_z),max(extent_z),by=cell_size_z)

D=sqrt((max(extent_x)-min(extent_x))**2+(max(extent_y)-min(extent_y))**2+(max(extent_z)-min(extent_z))**2)

nx=length(xx)
ny=length(yy)
nz=length(zz)


##INITIALIAL MODEL FOR PPO
closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(length(lithology_names))
registerDoParallel(Cl)
phi_g_fields_list=foreach(i=1:nsdf)%dopar%
{
  library(data.table)
  as.data.frame(fread(paste("./outputs/data/rbf/psd_rbf_sdf_",i,".csv",sep="")))[,4]
}
closeAllConnections()



## VOLUME PROPORTION CONSTRAINTS
vol_min=list()
vol_max=list()

vol_limit=readRDS("./inputs/data/volume_proportion_limits.RDS")
vol_limit

for ( i in 1:length(names(table(data_points[,5]))))
{
  vol_min[[i]]=1.001*as.numeric(vol_limit[i,2])
  vol_max[[i]]=0.999*as.numeric(vol_limit[i,3])
}



## VARIOGRAM COMPUTING OF PSEUDO SIGNED DISTANCE VALUES
closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(length(lithology_names))
registerDoParallel(Cl)

db_data=foreach(i=1:nsdf) %dopar%
{
  library(RGeostats)
  db.create(x1=data_points[,2],x2=data_points[,3],x3=data_points[,4],z1=psd_data[,i])
}

vario_exp=foreach(i=1:nsdf) %dopar% 
{
  library(RGeostats)
  nlag=20
  lag=D/(2*nlag)
  vario.calc(db_data[[i]],lag=D/(2*nlag),nlag=nlag)
}

vario_model=foreach(i=1:nsdf) %dopar%
{
  library(RGeostats)
  model.auto(vario_exp[[i]],struct = melem.name(c(3)),upper=c(paste("M1V=",var(psd_data[,i]),sep=""),paste("M1R=",200,sep="")),auth.aniso =FALSE,auth.locksame=FALSE,draw=FALSE,flag.noreduce=TRUE)
}

closeAllConnections()


## STRUCTURAL PARAMETERS (RANGE AND SILL) PRIOR DISTRIBUTION (UNIFORM)
min_range=list()
max_range=list()
min_sill=list()
max_sill=list()


for ( i in 1:nsdf)
{
  min_range[[i]]=0.50*vario_model[[i]]@basics[[1]]@range
  max_range[[i]]=vario_model[[i]]@basics[[1]]@range
  min_sill[[i]]=0.50*vario_model[[i]]@basics[[1]]@sill
  max_sill[[i]]=vario_model[[i]]@basics[[1]]@sill
}




vmin_max=function(x,minx=vol_min[[1]],maxx=vol_max[[1]])
{
  output1=output2=output=0
  if(sum(as.numeric(round(x,digits =4)<minx))>=1){output1=1}
  if(sum(as.numeric(round(x,digits =4)>maxx))>=1){output2=1}
  if((output1+output2)>=1){output=1}
  return(output)
}

## LOADING PPO RESULTS
niter=100
nrun=75
Mean_Error_Select=matrix(NA,nrow=nrun,ncol=nsdf)
Var_Error_Select=matrix(NA,nrow=nrun,ncol=nsdf)
Vol_Prop_Select=matrix(NA,nrow=nrun,ncol=nsdf)
Range_Select=matrix(NA,nrow=nrun,ncol=nsdf)
Sill_Select=matrix(NA,nrow=nrun,ncol=nsdf)
Phi_g_fields_Select=list()
for (l in 1:nsdf)
{
  Phi_g_fields_Select[[l]]=matrix(NA,ncol=nrun,nrow=nrow(full_grid_points[[l]]))
}


path="./outputs/data/ppo/runs/"
path_="./outputs/figures/5_ppo_model_selection/"


## MULTIPLE RUN ANALYSIS ##
for (i in 1:nrun)
{
  Res=readRDS(paste(path,i,"/Result.rds",sep=""))

  
  ## OBJECTIVE FUNCTION AND OPTIMAL PARAMETERS VISUALIZATION
  png(paste(path,i,"/UNIQUE_OBJECTIVE_FUNCTION_",i,".png",sep=""),width=5,height=5,units="in",res=300)
  plot(Res$U_f,xlab="# Iterations",ylab="Objective Function",type="l",lwd=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
  dev.off()
  
  
  for ( l in 1:nsdf)
  {
    png(paste(path,i,"/OPTIMAL_DEFORMATION_PARAMETER_",l,i,".png",sep=""),width=5,height=5,units="in",res=300)
    plot(Res$U_r[,l],xlab="# Iterations",ylab="Optimal Parameter r",type="l",lwd=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
    dev.off()
  }
  
  for ( l in 1:nsdf)
  {
    png(paste(path,i,"/OPTIMAL_SILL_PARAMETER_",l,i,".png",sep=""),width=5,height=5,units="in",res=300)
    plot(Res$sill_opt[,l],xlab="# Iterations",ylab="Optimal Sill",type="l",lwd=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
    dev.off()
  }
  
  for ( l in 1:nsdf)
  {
    png(paste(path,i,"/OPTIMAL_RANGE_PARAMETER_",l,i,".png",sep=""),width=5,height=5,units="in",res=300)
    plot(Res$range_opt[,l],xlab="# Iterations",ylab="Optimal Range",type="l",lwd=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
    dev.off()
  }
  
  
  ## MEAN AND VARIANCE OF SD AT CONTACT POINTS
  Mean_Error=matrix(NA,nrow=niter,ncol=nsdf)
  Var_Error=matrix(NA,nrow=niter,ncol=nsdf)
  
  for ( l in 1:nsdf)
  {
    
    for ( k in 1:niter)
    {  
      Mean_Error[k,l]=mean(Res$phi_current_list[[k]][[l]][-1:-(nrow(grid_points)+nrow(data_points))])
      Var_Error[k,l]=var(Res$phi_current_list[[k]][[l]][-1:-(nrow(grid_points)+nrow(data_points))])
    }

    png(paste(path,i,"/MEAN_SD_CONTACT_POINTS_",l,i,".png",sep=""),width=5,height=5,units="in",res=300)
    plot(Mean_Error[,l],xlab="# Iterations",ylab="Mean of SD at Contact Points",type="l",cex=1.5,lwd=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,ylim=c(-5,5))
    abline(h=0,col=2,lty=2)
    dev.off()
    
    png(paste(path,i,"/VARIANCE_SD_CONTACT_POINTS_",l,i,".png",sep=""),width=5,height=5,units="in",res=300)
    plot(Var_Error[,l],xlab="# Iterations",ylab="Variance of SD at Contact Points",type="l",cex=1.5,lwd=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
    dev.off()
  }
  
  
  
  vp=matrix(NA,nrow=niter,ncol=nsdf)
  domain_grid=matrix(NA,ncol=niter,nrow=nrow(grid_points))
  for (l in 1:niter)
  {
    mm=as.data.frame(cbind(Res$phi_current_list[[l]][[1]][1:nrow(grid_points)],Res$phi_current_list[[l]][[2]][1:nrow(grid_points)],
                           Res$phi_current_list[[l]][[3]][1:nrow(grid_points)],Res$phi_current_list[[l]][[4]][1:nrow(grid_points)]))
    colnames(mm)=1:nsdf
    domain_grid[,l]=R_truncation_rule_domaining(mm)
    vp[l,]=table(domain_grid[,l])/nrow(grid_points)
  }
  
  vp_test=apply(vp,1,vmin_max)
  i_v_test=which(vp_test==1)
  
  
  
  
  mis_classification_rate=rep(NA,niter)

  for (l in 1:niter)
  {
    mm_data=as.data.frame(cbind(Res$phi_current_list[[l]][[1]][-1:-nrow(grid_points)][1:nrow(data_points)],Res$phi_current_list[[l]][[2]][-1:-nrow(grid_points)][1:nrow(data_points)],
                                Res$phi_current_list[[l]][[3]][-1:-nrow(grid_points)][1:nrow(data_points)],Res$phi_current_list[[l]][[4]][-1:-nrow(grid_points)][1:nrow(data_points)]))
    colnames(mm_data)=1:nsdf
    domain=R_truncation_rule_domaining(mm_data)
    
    mis_classification_rate[l]=mean(as.numeric(domain!=as.integer(data_points[,5])))
  }
  png(paste(path,i,"/MIS_CLASSIFICATION_",i,".png",sep=""),width=5,height=5,units="in",res=300)
  plot(mis_classification_rate,type="l",xlab="#Iterations",ylab="Misclassification Rate",cex=1.5,lwd=1.5,cex.lab=1.5,cex.axis=1.5)
  dev.off()

  
  mm_init=as.data.frame(cbind(phi_g_fields_list[[1]][1:nrow(grid_points)],phi_g_fields_list[[2]][1:nrow(grid_points)],
                              phi_g_fields_list[[3]][1:nrow(grid_points)],phi_g_fields_list[[4]][1:nrow(grid_points)]))
  colnames(mm_init)=1:nsdf
  domain_init=R_truncation_rule_domaining(mm_init)
  
  vp_init=table(domain_init)/length(domain_init)
  vp_init
  
  for ( l in 1:(nsdf))
  {
    
    png(paste(path,i,"/VOLUME_PRO_",l,i,".png",sep=""),width=5,height=5,units="in",res=300)
    plot(0:niter,c(vp_init[l],vp[,l]),xlab="# Iterations",ylab="Volume Proportion",type="l",cex=1.5,lwd=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,ylim=c(min(vol_min[[l]],vp[,l]),max(vol_max[[l]],vp[,l])))
    abline(h=vol_min[[l]],col=2,lty=2)
    abline(h=vol_max[[l]],col=2,lty=2)
    dev.off()
  }
  
  
  i_v_test=which(vp_test==1)
  if((length(i_v_test)>=1)&(length(i_v_test)<niter)){i_v_test=i_v_test[-which(i_v_test<30)];Sel_Ind=setdiff(30:niter,i_v_test)}else{if((length(i_v_test)==0)|(length(i_v_test)==niter)){Sel_Ind=30:niter}}
  
  
  ## PPO MODEL SELECTION
  abs_sum_Mean_Error=apply(abs(Mean_Error),1,mean)
  Select=rep(NA,nsdf)
  Select=which(round(abs(abs_sum_Mean_Error[Sel_Ind]),digits=2)==min(round(abs(abs_sum_Mean_Error[Sel_Ind]),digits=2)))[1]
  Sel_Model=Sel_Ind[1]+Select-1
  Sel_Model
  
  for ( l in 1:(nsdf))
  {
    Vol_Prop_Select[i,l]=vp[Sel_Model,l] 
  }
  
  for ( l in 1:nsdf)
  {
    
    Mean_Error_Select[i,l]=Mean_Error[Sel_Model,l]
    Var_Error_Select[i,l]=Var_Error[Sel_Model,l]
    Range_Select[i,l]=Res$range_opt[Sel_Model,l]
    Sill_Select[i,l]=Res$sill_opt[Sel_Model,l]
    Phi_g_fields_Select[[l]][,i]=Res$phi_current_list[[Sel_Model]][[l]]
    
    
    ## VARIOGRAM OF RESIDUALS AT CONTACT POINTS FOR THE SELECTED MODEL
    db_data_e=db.create(x1=contact_points_cat[[l]][,2],x2=contact_points_cat[[l]][,3],x3=contact_points_cat[[l]][,4],z1=Res$phi_current_list[[Sel_Model]][[l]][-1:-nrow(grid_points)][-1:-nrow(psd_data)])
    nlag=15
    lag=D/(3*nlag)
    vario_exp_e=vario.calc(db_data_e,lag=D/(3*nlag),nlag=nlag)

    vario_model_e=model.auto(vario_exp_e,struct = melem.name(c(3)),auth.aniso =FALSE,auth.locksame=FALSE,draw=FALSE,flag.noreduce=FALSE)

    png(paste(path_,"variogram_residuals/VARIO_SD_",l,"_",i,".png",sep=""),width=5,height=5,units="in",res=300)
    vario.plot(vario_exp_e,main=paste("Sim",i,sep=" "),reset=TRUE,flag.norm = FALSE,cex.axis=1.5,cex.lab=1.5,lty=2,lwd=2,npairdw = FALSE,varline =FALSE,xlab="Distance (m)",ylab="Variogram")
    model.plot(vario_model_e,vario=vario_exp_e,add=T,cex.axis=1.5,lwd=c(2),lty=c(1),col=2)
    dev.off()
    
  }
  
  cat("Selection ",i,"Sucessful","\n")
}


## SAVING OPTIMAL SILL AND RANGE 
saveRDS(Sill_Select,"./outputs/data/ppo/PPO_Optimal_Sill.RDS")
saveRDS(Range_Select,"./outputs/data/ppo/PPO_Optimal_Range.RDS")



for ( l in 1:nsdf)
{
  
png(paste("./outputs/figures/5_ppo_model_selection/Final_Mean_sd_contact_points_",lithology_names[l],".png",sep=""),width=5,height=5,units="in",res=300)
plot(Mean_Error_Select[,l],xlab="# Realizations",ylab="Mean of SD at Contact Points",type="l",cex=1.5,lwd=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,ylim=c(-5,5))
abline(h=0,col=2,lty=2)
dev.off()
}

for ( l in 1:nsdf)
{
  
  png(paste("./outputs/figures/5_ppo_model_selection/Final_Vol_Prop_",lithology_names[l],".png",sep=""),width=5,height=5,units="in",res=300)
  plot(Vol_Prop_Select[,l],xlab="# Realizations",ylab="Volume Proportion",type="l",cex=1.5,lwd=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,ylim=c(vol_min[[l]],vol_max[[l]]))
   abline(h=vol_min[[l]],col=2,lty=2)
  abline(h=vol_max[[l]],col=2,lty=2)
  dev.off()
}

for ( l in 1:nsdf)
{
  
  png(paste("./outputs/figures/5_ppo_model_selection/Final_Range_",lithology_names[l],".png",sep=""),width=5,height=5,units="in",res=300)
  hist(Range_Select[,l],xlab="Optimal Range (m)",col="blue",main="",cex=1.5,lwd=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
  dev.off()
}



for ( l in 1:nsdf)
{
  
  png(paste("./outputs/figures/5_ppo_model_selection/Final_Sill_",lithology_names[l],".png",sep=""),width=5,height=5,units="in",res=300)
  hist(Sill_Select[,l],xlab="Optimal Sill (m2)",col="blue",main="",cex=1.5,lwd=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
  dev.off()
}




## DOMAIN TREND MODELS SELECTED
domain_grid_select=matrix(NA,ncol=nrun,nrow=nrow(grid_points))
for (l in 1:nrun)
{
  mm_select=as.data.frame(cbind(Phi_g_fields_Select[[1]][1:nrow(grid_points),l],Phi_g_fields_Select[[2]][1:nrow(grid_points),l],
                                Phi_g_fields_Select[[3]][1:nrow(grid_points),l],Phi_g_fields_Select[[4]][1:nrow(grid_points),l]))
  domain_grid_select[,l]=R_truncation_rule_domaining(mm_select)
}

saveRDS(domain_grid_select,"./outputs/data/ppo/trend_categorical_model.RDS")
saveRDS(Phi_g_fields_Select,"./outputs/data/ppo/trend_sdf_model.RDS")


write_feather(as.data.frame(domain_grid_select),"./outputs/data/ppo/trend_categorical_model.feather")
write_feather(as.data.frame(Phi_g_fields_Select[[1]][1:(nrow(grid_points)+nrow(data_points)),]),"./outputs/data/ppo/trend_sdf_model_1.feather")
write_feather(as.data.frame(Phi_g_fields_Select[[2]][1:(nrow(grid_points)+nrow(data_points)),]),"./outputs/data/ppo/trend_sdf_model_2.feather")
write_feather(as.data.frame(Phi_g_fields_Select[[3]][1:(nrow(grid_points)+nrow(data_points)),]),"./outputs/data/ppo/trend_sdf_model_3.feather")
write_feather(as.data.frame(Phi_g_fields_Select[[4]][1:(nrow(grid_points)+nrow(data_points)),]),"./outputs/data/ppo/trend_sdf_model_4.feather")


write_feather(as.data.frame(Phi_g_fields_Select[[1]][-1:-(nrow(grid_points)+nrow(data_points)),]),"./outputs/data/ppo/trend_sdf_contacts_1.feather")
write_feather(as.data.frame(Phi_g_fields_Select[[2]][-1:-(nrow(grid_points)+nrow(data_points)),]),"./outputs/data/ppo/trend_sdf_contacts_2.feather")
write_feather(as.data.frame(Phi_g_fields_Select[[3]][-1:-(nrow(grid_points)+nrow(data_points)),]),"./outputs/data/ppo/trend_sdf_contacts_3.feather")
write_feather(as.data.frame(Phi_g_fields_Select[[4]][-1:-(nrow(grid_points)+nrow(data_points)),]),"./outputs/data/ppo/trend_sdf_contacts_4.feather")




## DOMAIN DATA MODELS SELECTED
domain_data_select=matrix(NA,ncol=nrun,nrow=nrow(data_points))
for (l in 1:nrun)
{
  mm_data_select=as.data.frame(cbind(Phi_g_fields_Select[[1]][,l][-1:-nrow(grid_points)][1:nrow(psd_data)],Phi_g_fields_Select[[2]][,l][-1:-nrow(grid_points)][1:nrow(psd_data)],
                                     Phi_g_fields_Select[[3]][,l][-1:-nrow(grid_points)][1:nrow(psd_data)],Phi_g_fields_Select[[4]][,l][-1:-nrow(grid_points)][1:nrow(psd_data)]))
  domain_data_select[,l]=R_truncation_rule_domaining(mm_data_select)
                                                     
}
write_feather(as.data.frame(domain_data_select),"./outputs/data/ppo/trend_categorical_data.feather")



## MISCLASSIFICATION RATE 
mis_classification_rate_sel_g=rep(NA,nrun)
for ( i in 1:nrun)
{
  mis_classification_rate_sel_g[i]=sum(as.numeric(domain_data_select[,i]!=data_points[,5]))/nrow(data_points)
  
}

saveRDS(mis_classification_rate_sel_g,"./outputs/data/ppo/mis_classification_rate_ppo_global.RDS")


png("./outputs/figures/5_ppo_model_selection/Mis_Rate_global.png",width=5,height=5,units="in",res=300)
hist(100*mis_classification_rate_sel_g,xlab="Misclassification Rate (%)",col="blue",main=paste("Mean = ",round(mean(100*mis_classification_rate_sel_g),digits=2), sep=" "),cex=1.5,lwd=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
dev.off()

closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(length(lithology_names))
registerDoParallel(Cl)
mis_classification_rate_sel_cat=foreach(k=1:length(lithology_names),.combine = cbind)%dopar%
{
  mis_classification_rate_sel=rep(NA,nrun)
  ii_cat=which(data_points[,5]==k)
  for ( i in 1:nrun)
  {
    mis_classification_rate_sel[i]=sum(as.numeric(domain_data_select[ii_cat,i]!=data_points[ii_cat,5]))/length(ii_cat)
  }
  mis_classification_rate_sel
}
saveRDS(mis_classification_rate_sel_cat,"./outputs/data/ppo/mis_classification_rate_ppo_category.RDS")


for ( l in 1:length(lithology_names))
{
  png(paste("./outputs/figures/5_ppo_model_selection/Mis_Rate_",lithology_names[l],".png",sep=""),width=5,height=5,units="in",res=300)
  hist(100*mis_classification_rate_sel_cat[,l],xlab="Misclassification Rate (%)",col="blue",main=paste("Mean = ",round(mean(100*mis_classification_rate_sel_cat[,l]),digits=2), sep=" "),cex=1.5,lwd=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
  dev.off()
  
}


##PROBABILITY MAPS##
ppo_domains=domain_grid_select


PPO_Prob_Map_Mean=list()
PPO_Prob_Map_Var=list()

for ( j in 1:length(lithology_names))
{
  tempd=1*(ppo_domains==j)
  PPO_Prob_Map_Mean[[j]]=apply(tempd,1,mean)
  PPO_Prob_Map_Var[[j]]=apply(tempd,1,var)
}


PPO_Probability_Mean=cbind(PPO_Prob_Map_Mean[[1]],PPO_Prob_Map_Mean[[2]],PPO_Prob_Map_Mean[[3]],PPO_Prob_Map_Mean[[4]])
colnames(PPO_Probability_Mean)=lithology_names

PPO_Probability_Var=cbind(PPO_Prob_Map_Var[[1]],PPO_Prob_Map_Var[[2]],PPO_Prob_Map_Var[[3]],PPO_Prob_Map_Var[[4]])
colnames(PPO_Probability_Var)=lithology_names


array_Prob_Map_Mean_1=array(PPO_Prob_Map_Mean[[1]], dim=c(nx,ny,nz))
png("./outputs/figures/5_ppo_model_selection/probability_maps/1/PPO_Prob_Map_Mean_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_1, xs=max(data_points$X)-200, ys = min(data_points$Y)+200, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/1/PPO_Prob_Map_Mean_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_1, xs=min(data_points$X)+500, ys = min(data_points$Y)+200, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/1/PPO_Prob_Map_Mean_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_1, xs=min(data_points$X)+500, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()


array_Prob_Map_Var_1=array(PPO_Prob_Map_Var[[1]],dim=c(nx,ny,nz))
png("./outputs/figures/5_ppo_model_selection/probability_maps/1/PPO_Prob_Map_Var_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_1, xs=max(data_points$X)-200, ys = min(data_points$Y)+200, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/1/PPO_Prob_Map_Var_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_1, xs=min(data_points$X)+500, ys = min(data_points$Y)+200, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/1/PPO_Prob_Map_Var_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_1, xs=min(data_points$X)+500, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()


array_Prob_Map_Mean_2=array(PPO_Prob_Map_Mean[[2]], dim=c(nx,ny,nz))
png("./outputs/figures/5_ppo_model_selection/probability_maps/2/PPO_Prob_Map_Mean_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_2, xs=max(data_points$X)-300, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/2/PPO_Prob_Map_Mean_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_2, xs=min(data_points$X)+500, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/2/PPO_Prob_Map_Mean_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_2, xs=min(data_points$X)+200, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()


array_Prob_Map_Var_2=array(PPO_Prob_Map_Var[[2]],dim=c(nx,ny,nz))
png("./outputs/figures/5_ppo_model_selection/probability_maps/2/PPO_Prob_Map_Var_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_2, xs=max(data_points$X)-300, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/2/PPO_Prob_Map_Var_BRXH_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_2, xs=min(data_points$X)+500, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/2/PPO_Prob_Map_Var_BRXH_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_2, xs=min(data_points$X)+200, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

array_Prob_Map_Mean_3=array(PPO_Prob_Map_Mean[[3]], dim=c(nx,ny,nz))
png("./outputs/figures/5_ppo_model_selection/probability_maps/3/PPO_Prob_Map_Mean_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_3, xs=max(data_points$X)-300, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/3/PPO_Prob_Map_Mean_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_3, xs=min(data_points$X)+500, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/3/PPO_Prob_Map_Mean_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_3, xs=min(data_points$X)+200, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()


array_Prob_Map_Var_3=array(PPO_Prob_Map_Var[[3]],dim=c(nx,ny,nz))
png("./outputs/figures/5_ppo_model_selection/probability_maps/3/PPO_Prob_Map_Var_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_3, xs=max(data_points$X)-300, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/3/PPO_Prob_Map_Var_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_3, xs=min(data_points$X)+500, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/3/PPO_Prob_Map_Var_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_3, xs=min(data_points$X)+200, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()


array_Prob_Map_Mean_4=array(PPO_Prob_Map_Mean[[4]], dim=c(nx,ny,nz))
png("./outputs/figures/5_ppo_model_selection/probability_maps/4/PPO_Prob_Map_Mean_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_4, main=" ",xs=max(data_points$X)-300, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1.5,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/4/PPO_Prob_Map_Mean_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_4, main=" ",xs=min(data_points$X)+500, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1.5,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/4/PPO_Prob_Map_Mean_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_4, main=" ",xs=min(data_points$X)+200, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1.5,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()


array_Prob_Map_Var_4=array(PPO_Prob_Map_Var[[4]],dim=c(nx,ny,nz))
png("./outputs/figures/5_ppo_model_selection/probability_maps/4/PPO_Prob_Map_Var_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_4, main=" ",xs=max(data_points$X)-300, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1.5,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/4/PPO_Prob_Map_Var_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_4, main=" ",xs=min(data_points$X)+500, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1.5,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/5_ppo_model_selection/probability_maps/4/PPO_Prob_Map_Var_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_4, main=" ",xs=min(data_points$X)+200, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1.5,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

T2=Sys.time()
difftime(T2,T1)


