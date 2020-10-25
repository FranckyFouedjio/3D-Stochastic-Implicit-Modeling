T1=Sys.time()
#################################################################################################################
#                                                                                                               #
#                         SYNTHETIC CASE STUDY: PSEUDO SIGNED DISTANCE COMPUTING                                #
#                                                                                                               #
#################################################################################################################

## WORKING DIRECTORY SETTING
setwd("/media/fouedjio/095e0241-4724-4a1f-84f8-9ddda0df2da9/fouedjio/3d_stochastic_implicit_modeling/")

## REQUIRING PACKAGES LOADING
library(data.table)
library(png)
library(plot3Drgl)
library(profvis)
library(FNN)
library(doParallel)
library(rgl)
library(scatterplot3d)
library(Rvcg)
library(misc3d)
library(ggplot2)
library(RGeostats)

## LOADING SOURCE CODE
source("./code/functions/pseudo_signed_distance_function.R")
source("./code/functions/dichotomization_function.R")
source("./code/functions/data_visualization_function.R")

##CREATING SOME FOLDERS
if(!file.exists("./outputs")){dir.create("./outputs")}
if(!file.exists("./outputs/figures/")){dir.create("./outputs/figures")}
if(!file.exists("./outputs/figures/1_psd_computing")){dir.create("./outputs/figures/1_psd_computing")}
if(!file.exists("./outputs/data/")){dir.create("./outputs/data")}

## DATA LOADING
data_points=as.data.frame(fread("./inputs/data/drillhole_data.csv"))
contact_points=as.data.frame(fread("./inputs/data/contact_points.csv"))
contact_points_cat=readRDS('./inputs/data/contact_points_cat.RDS')
lithology_names=1:4


##DRILLHOLE DATA VISUALIZATION (WITHOUT CONTACT POINTS)
R_plot_drillhole_lithology(data_points)
writeWebGL(dir="./outputs/figures/1_psd_computing/", filename =paste("./outputs/figures/1_psd_computing/plot_drillhole_data_without_contact_points.html"),width = 2000, reuse = TRUE)


##DRILLHOLE DATA VISUALIZATION (WITH CONTACT POINTS)
R_plot_drillhole_lithology_with_contacts(data_points,contact_points)
writeWebGL(dir="./outputs/figures/1_psd_computing/", filename =paste("./outputs/figures/1_psd_computing/plot_drillhole_data_with_contact_points.html"),width = 2000, reuse = TRUE)


## PSEUDO SIGNED DISTANCE FUNCTION COMPUTING
closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(NbCores)
registerDoParallel(Cl)
psd_data=R_pseudo_signed_distance_p(data_points,contact_points_cat)
psd_data=as.data.frame(psd_data)
colnames(psd_data)=c("psd1","psd2","psd3","psd4")
fwrite(psd_data,"./outputs/data/psd_data.csv",row.names = FALSE)
closeAllConnections()


## PSEUDO SIGNED DISTANCES VISUALIZATION
closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(NbCores)
registerDoParallel(Cl)
foreach(i=1:length(lithology_names)) %dopar%
{
  library(plot3Drgl)
  scatter3Drgl(data_points[,2],data_points[,3],data_points[,4],colvar=psd_data[,i],scale=FALSE,expand=1,pch=20,cex=1,phi=40,theta=0,colkey= list(length=1,width=1,cex.axis=1),zoom=0.85)
  writeWebGL(dir="./outputs/figures/1_psd_computing/", filename =paste("./outputs/figures/1_psd_computing/plot_psd_",lithology_names[i],"_data.html"),width = 2000, reuse = TRUE)
}
closeAllConnections()

closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(NbCores)
registerDoParallel(Cl)
foreach(i=1:length(lithology_names)) %dopar%
{
  library(plot3D)
  png(paste("./outputs/figures/1_psd_computing/plot_psd_",lithology_names[i],"_data.png",sep=""),width=7,height=7,units="in",res=300)
  scatter3D(data_points[,2],data_points[,3],data_points[,4],colvar=psd_data[,i],scale=FALSE,axes=F,expand=1,pch=20,cex=1,phi=40,theta=0,colkey= list(length=1,width=1,cex.axis=1),zoom=0.85)
  dev.off()
}
closeAllConnections()


## DATA PREPARATION FOR RBF SMOOTHING/INTERPOLATING OF PSD
closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(NbCores)
registerDoParallel(Cl)
data=foreach(i=1:length(lithology_names)) %dopar%
{
  cbind(data_points[,2:4],D=psd_data[,i])
}

foreach(i=1:length(lithology_names)) %dopar%
{
  library(data.table)
  fwrite(data[[i]],paste("./outputs/data/data_points_rbf_",lithology_names[i],".csv",sep=""),row.names = FALSE)
}
closeAllConnections()

## INDICATOR COMPUTING
closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(NbCores)
registerDoParallel(Cl)
data_points_binary=R_dichotomization(as.matrix(data_points[,-1]))
closeAllConnections()
fwrite(as.data.frame(data_points_binary),"./outputs/data/drillhole_data_binary.csv",row.names = FALSE)

T2=Sys.time()
difftime(T2,T1)

