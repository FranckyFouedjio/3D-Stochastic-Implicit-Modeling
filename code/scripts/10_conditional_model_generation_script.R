T1=Sys.time()
###############################################################################################
#                                                                                             #
#            SYNTHETIC CASE STUDY: CONDITIONAL MODEL GENERATION (CATEGORICAL MODELS)          #
#                                                                                             #
###############################################################################################

## WORKING DIRECTORY SETTING
setwd("/media/fouedjio/095e0241-4724-4a1f-84f8-9ddda0df2da9/fouedjio/3d_stochastic_implicit_modeling/")


## Packages Loading
library(factoextra)
library(data.table)
library(foreach)
library(doParallel)
library(plot3Drgl)
library(plot3D)
library(pracma)
library(FNN)
library(ggfortify)
library(polycor)
library(psych)
library(philentropy)
library(MASS)
library(ggplot2)
library(misc3d)
library(CCA)
library(RGeostats)
library(rockchalk)
library(Rglpk)
library(DGSA)
library(mvtnorm)
library(feather)

##SOURCE CODE LOADING
source("./code/functions/truncation_rule_domaining_function.R")
source("./code/functions/data_visualization_function.R")

##CREATING SOME FOLDERS
if(!file.exists("./outputs/data/conditional_model_generation")){dir.create("./outputs/data/conditional_model_generation")}
if(!file.exists("./outputs/figures/10_conditional_model_generation")){dir.create("./outputs/figures/10_conditional_model_generation")}
if(!file.exists("./outputs/figures/10_conditional_model_generation/volume_proportion")){dir.create("./outputs/figures/10_conditional_model_generation/volume_proportion")}
if(!file.exists("./outputs/figures/10_conditional_model_generation/probability_maps")){dir.create("./outputs/figures/10_conditional_model_generation/probability_maps")}
if(!file.exists("./outputs/figures/10_conditional_model_generation/probability_maps/1")){dir.create("./outputs/figures/10_conditional_model_generation/probability_maps/1")}
if(!file.exists("./outputs/figures/10_conditional_model_generation/probability_maps/2")){dir.create("./outputs/figures/10_conditional_model_generation/probability_maps/2")}
if(!file.exists("./outputs/figures/10_conditional_model_generation/probability_maps/3")){dir.create("./outputs/figures/10_conditional_model_generation/probability_maps/3")}
if(!file.exists("./outputs/figures/10_conditional_model_generation/probability_maps/4")){dir.create("./outputs/figures/10_conditional_model_generation/probability_maps/4")}
if(!file.exists("./outputs/figures/10_conditional_model_generation/models")){dir.create("./outputs/figures/10_conditional_model_generation/models")}

## DATA LOADING
data_points=as.data.frame(fread("./inputs/data/drillhole_data.csv"))
contact_points=as.data.frame(fread("./inputs/data/contact_points.csv"))
grid_points=as.data.frame(fread("./inputs/data/grid_points.csv"))
contact_points_cat=readRDS('./inputs/data/contact_points_cat.RDS')
data_points_binary=as.data.frame(fread("./outputs/data/drillhole_data_binary.csv",header=TRUE))
psd_data=as.data.frame(fread("./outputs/data/psd_data.csv"))
lithology_names=1:4
ncat=nsdf=length(lithology_names)



closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(length(lithology_names))
registerDoParallel(Cl)

full_grid_points=foreach(i=1:nsdf) %dopar%
{
  library(data.table)
  rbindlist(list(grid_points[,1:3],data_points[,2:4]))
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


full_grid=rbindlist(list(grid_points[,1:3],data_points[,2:4]))



nsimu=1000
n_pos_simu=1000



## RECONSTRUCTION CONDITIONAL LEVEL SETS


signed_distance_fields_posterior_l=list()
for(j in 1:nsdf)
{
  signed_distance_fields_posterior_l[[j]]=read_feather(paste("./outputs/data/conditional_model_generation/conditional_sdf_model_",j,".feather",sep=""),columns=1:250)
}


lithology_domaining_posterior=list()

for ( i in 1:250)
{
  lithology_domaining_posterior[[i]]=R_truncation_rule_domaining(cbind(signed_distance_fields_posterior_l[[1]][1:nrow(full_grid),i],signed_distance_fields_posterior_l[[2]][1:nrow(full_grid),i],signed_distance_fields_posterior_l[[3]][1:nrow(full_grid),i],
                                                                       signed_distance_fields_posterior_l[[4]][1:nrow(full_grid),i]))
}
rm(signed_distance_fields_posterior_l)
gc()


signed_distance_fields_posterior_l=list()
for(j in 1:nsdf)
{
  signed_distance_fields_posterior_l[[j]]=read_feather(paste("./outputs/data/conditional_model_generation/conditional_sdf_model_",j,".feather",sep=""),columns=251:500)
}


for ( i in 1:250)
{
  lithology_domaining_posterior[[i+250]]=R_truncation_rule_domaining(cbind(signed_distance_fields_posterior_l[[1]][1:nrow(full_grid),i],signed_distance_fields_posterior_l[[2]][1:nrow(full_grid),i],signed_distance_fields_posterior_l[[3]][1:nrow(full_grid),i],
                                                                           signed_distance_fields_posterior_l[[4]][1:nrow(full_grid),i]))
}
rm(signed_distance_fields_posterior_l)
gc()



signed_distance_fields_posterior_l=list()
for(j in 1:nsdf)
{
  signed_distance_fields_posterior_l[[j]]=read_feather(paste("./outputs/data/conditional_model_generation/conditional_sdf_model_",j,".feather",sep=""),columns=501:750)
}


for ( i in 1:250)
{
  lithology_domaining_posterior[[i+500]]=R_truncation_rule_domaining(cbind(signed_distance_fields_posterior_l[[1]][1:nrow(full_grid),i],signed_distance_fields_posterior_l[[2]][1:nrow(full_grid),i],signed_distance_fields_posterior_l[[3]][1:nrow(full_grid),i],
                                                                           signed_distance_fields_posterior_l[[4]][1:nrow(full_grid),i]))
}
rm(signed_distance_fields_posterior_l)
gc()


signed_distance_fields_posterior_l=list()
for(j in 1:nsdf)
{
  signed_distance_fields_posterior_l[[j]]=read_feather(paste("./outputs/data/conditional_model_generation/conditional_sdf_model_",j,".feather",sep=""),columns=751:1000)
}


for ( i in 1:250)
{
  lithology_domaining_posterior[[i+750]]=R_truncation_rule_domaining(cbind(signed_distance_fields_posterior_l[[1]][1:nrow(full_grid),i],signed_distance_fields_posterior_l[[2]][1:nrow(full_grid),i],signed_distance_fields_posterior_l[[3]][1:nrow(full_grid),i],
                                                                           signed_distance_fields_posterior_l[[4]][1:nrow(full_grid),i]))
}
rm(signed_distance_fields_posterior_l)
gc()


lithology_domains_posterior=matrix(unlist(lithology_domaining_posterior),ncol=n_pos_simu,nrow=nrow(full_grid),byrow = FALSE)
rm(lithology_domaining_posterior)
gc()

write_feather(as.data.frame(lithology_domains_posterior[1:nrow(grid_points),]), "./outputs/data/conditional_model_generation/conditional_categorical_model.feather")
write_feather(as.data.frame(lithology_domains_posterior[(nrow(grid_points)+1):(nrow(grid_points)+nrow(data_points)),]), "./outputs/data/conditional_model_generation/conditional_categorical_data.feather")



## MISMATCH PERCENTAGE
posterior_LUR_train_domains=lithology_domains_posterior[-1:-nrow(grid_points),][1:nrow(data_points),]
real_train_domains=as.numeric(data_points[,5])
res_match_train=1*(posterior_LUR_train_domains==real_train_domains)
match_percent_train=apply(res_match_train,2,mean)
cat(summary(match_percent_train))



##PROBABILITY MAPS##
lur_domains=lithology_domains_posterior[1:nrow(grid_points),]

LUR_Prob_Map_Mean=list()
LUR_Prob_Map_Var=list()

for ( j in 1:length(lithology_names))
{
  tempd=1*(lur_domains==j)
  LUR_Prob_Map_Mean[[j]]=apply(tempd,1,mean)
  LUR_Prob_Map_Var[[j]]=apply(tempd,1,var)
}


LUR_Probability_Mean=cbind(LUR_Prob_Map_Mean[[1]],LUR_Prob_Map_Mean[[2]],LUR_Prob_Map_Mean[[3]],LUR_Prob_Map_Mean[[4]])
colnames(LUR_Probability_Mean)=lithology_names

LUR_Probability_Var=cbind(LUR_Prob_Map_Var[[1]],LUR_Prob_Map_Var[[2]],LUR_Prob_Map_Var[[3]],LUR_Prob_Map_Var[[4]])
colnames(LUR_Probability_Var)=lithology_names


array_Prob_Map_Mean_1=array(LUR_Prob_Map_Mean[[1]], dim=c(nx,ny,nz))
png("./outputs/figures/10_conditional_model_generation/probability_maps/1/LUR_Prob_Map_Mean_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_1, xs=max(data_points$X)-200, ys = min(data_points$Y)+200, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/1/LUR_Prob_Map_Mean_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_1, xs=min(data_points$X)+500, ys = min(data_points$Y)+200, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/1/LUR_Prob_Map_Mean_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_1, xs=min(data_points$X)+500, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()


array_Prob_Map_Var_1=array(LUR_Prob_Map_Var[[1]],dim=c(nx,ny,nz))
png("./outputs/figures/10_conditional_model_generation/probability_maps/1/LUR_Prob_Map_Var_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_1, xs=max(data_points$X)-200, ys = min(data_points$Y)+200, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/1/LUR_Prob_Map_Var_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_1, xs=min(data_points$X)+500, ys = min(data_points$Y)+200, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/1/LUR_Prob_Map_Var_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_1, xs=min(data_points$X)+500, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()


array_Prob_Map_Mean_2=array(LUR_Prob_Map_Mean[[2]], dim=c(nx,ny,nz))
png("./outputs/figures/10_conditional_model_generation/probability_maps/2/LUR_Prob_Map_Mean_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_2, xs=max(data_points$X)-300, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/2/LUR_Prob_Map_Mean_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_2, xs=min(data_points$X)+500, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/2/LUR_Prob_Map_Mean_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_2, xs=min(data_points$X)+200, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()


array_Prob_Map_Var_2=array(LUR_Prob_Map_Var[[2]],dim=c(nx,ny,nz))
png("./outputs/figures/10_conditional_model_generation/probability_maps/2/LUR_Prob_Map_Var_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_2, xs=max(data_points$X)-300, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/2/LUR_Prob_Map_Var_BRXH_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_2, xs=min(data_points$X)+500, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/2/LUR_Prob_Map_Var_BRXH_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_2, xs=min(data_points$X)+200, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

array_Prob_Map_Mean_3=array(LUR_Prob_Map_Mean[[3]], dim=c(nx,ny,nz))
png("./outputs/figures/10_conditional_model_generation/probability_maps/3/LUR_Prob_Map_Mean_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_3, xs=max(data_points$X)-300, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/3/LUR_Prob_Map_Mean_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_3, xs=min(data_points$X)+500, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/3/LUR_Prob_Map_Mean_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_3, xs=min(data_points$X)+200, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()


array_Prob_Map_Var_3=array(LUR_Prob_Map_Var[[3]],dim=c(nx,ny,nz))
png("./outputs/figures/10_conditional_model_generation/probability_maps/3/LUR_Prob_Map_Var_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_3, xs=max(data_points$X)-300, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/3/LUR_Prob_Map_Var_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_3, xs=min(data_points$X)+500, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/3/LUR_Prob_Map_Var_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_3, xs=min(data_points$X)+200, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()


array_Prob_Map_Mean_4=array(LUR_Prob_Map_Mean[[4]], dim=c(nx,ny,nz))
png("./outputs/figures/10_conditional_model_generation/probability_maps/4/LUR_Prob_Map_Mean_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_4, main=" ",xs=max(data_points$X)-300, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1.5,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/4/LUR_Prob_Map_Mean_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_4, main=" ",xs=min(data_points$X)+500, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1.5,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/4/LUR_Prob_Map_Mean_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Mean_4, main=" ",xs=min(data_points$X)+200, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1.5,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,1))
dev.off()


array_Prob_Map_Var_4=array(LUR_Prob_Map_Var[[4]],dim=c(nx,ny,nz))
png("./outputs/figures/10_conditional_model_generation/probability_maps/4/LUR_Prob_Map_Var_1.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_4, main=" ",xs=max(data_points$X)-300, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1.5,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/4/LUR_Prob_Map_Var_2.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_4, main=" ",xs=min(data_points$X)+500, ys = min(data_points$Y)+100, zs = 1500,theta = -40, phi = 40, scale=FALSE,expand=1.5,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()

png("./outputs/figures/10_conditional_model_generation/probability_maps/4/LUR_Prob_Map_Var_3.png",width=4,height=4,units="in",res=300)
slice3D(xx, yy, zz,colvar=array_Prob_Map_Var_4, main=" ",xs=min(data_points$X)+200, ys = min(data_points$Y)+200, zs = 1500+100,theta = -40, phi = 40, scale=FALSE,expand=1.5,axes=T, pch=20,cex=0.25,colkey= list(length=0.5,width=0.5,cex.axis=1),clim=c(0,0.26))
dev.off()




##VOLUME PROPORTION
volume_posterior_LUR=matrix(NA,nrow=n_pos_simu,ncol=ncat)
vol_limit=readRDS("./inputs/data/volume_proportion_limits.RDS")
vol_limit
low_vol=100*as.numeric(vol_limit[,2])
upp_vol=100*as.numeric(vol_limit[,3])

for ( l in 1:n_pos_simu)
{
  volume_posterior_LUR[l,]=100*(table(lur_domains[,l]))/nrow(grid_points)
}


png(paste("./outputs/figures/10_conditional_model_generation/volume_proportion/Hist_Vol_Proportion_",1,".png",sep=""),width=7,height=5,units="in",res=300)
hist(volume_posterior_LUR[,1],col="blue",xlab="Cat 1 Volume Proportion (%)",main="")
dev.off()

png(paste("./outputs/figures/10_conditional_model_generation/volume_proportion/Hist_Vol_Proportion_",2,".png",sep=""),width=7,height=5,units="in",res=300)
hist(volume_posterior_LUR[,2],col="blue",xlab="Cat 2 Volume Proportion (%)",main="")
dev.off()

png(paste("./outputs/figures/10_conditional_model_generation/volume_proportion/Hist_Vol_Proportion_",3,".png",sep=""),width=7,height=5,units="in",res=300)
hist(volume_posterior_LUR[,3],col="blue",xlab="Cat 3 Volume Proportion (%)",main="")
dev.off()

png(paste("./outputs/figures/10_conditional_model_generation/volume_proportion/Hist_Vol_Proportion_",4,".png",sep=""),width=7,height=5,units="in",res=300)
hist(volume_posterior_LUR[,4],col="blue",xlab="Cat 4 Volume Proportion (%)",main="")
dev.off()

##EXAMPLE OF CATEGORICAL MODELS VISUALIZATION
R_plot_block_model_lithology_1(lur_domains[,1],data_points)
R_plot_block_model_lithology_2(lur_domains[,1],data_points)
T2=Sys.time() 
difftime(T2,T1)




















