T1=Sys.time()
##########################################################################################################################################
#                                                                                                                                        #
#                                     SYNTHETIC CASE STUDY:  PSD RBF INTERPOLATION  (ALL FORMATIONS)                                     #
#                                                                                                                                        #
##########################################################################################################################################

## WORKING DIRECTORY SETTING
setwd("/media/fouedjio/095e0241-4724-4a1f-84f8-9ddda0df2da9/fouedjio/3d_stochastic_implicit_modeling/")

## REQUIRING PACKAGES LOADING
library(fields)
library(rgl)
library(plot3D)
library(foreach)
library(doParallel)
library(pracma)
library(reticulate)
library(e1071)
library(png)
library(plot3Drgl)
library(scatterplot3d)
library(Rvcg)
library(misc3d)
library(OceanView)
library(data.table)
library(FNN)
library(ggplot2)  
library(deldir)
library("RColorBrewer")
library(profvis)
library(DescTools)
library(geometry)




## SOURCE CODE LOADING
source("./code/functions/data_visualization_function.R")
source("./code/functions/dichotomization_function.R")
source("./code/functions/rbf_interpolation_smoothing_function.R")
source("./code/functions/signed_distance_conversion_function.R")
source("./code/functions/truncation_rule_domaining_function.R")


## PYTHON ENVIRONMENT SETTING
library(reticulate)
use_virtualenv("r-reticulate")
numpy=import("numpy")
scipy=import("scipy")
mp=import("multiprocessing")

##CREATING SOME FOLDERS
if(!file.exists("./outputs/data/rbf")){dir.create("./outputs/data/rbf")}
if(!file.exists("./outputs/figures/3_psd_rbf_interpolation")){dir.create("./outputs/figures/3_psd_rbf_interpolation")}


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

##
closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(NbCores)
registerDoParallel(Cl)
sdf=list()
for ( i in 1:nsdf)
{
  sdf[[i]]=as.data.frame(fread(paste("./outputs/data/rbf/","psd_rbf_sdf_",i,".csv",sep="")))
}



## PSEUDO SIGNED DISTANCE VALUES VISUALIZATION

closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(NbCores)
registerDoParallel(Cl)
foreach(i=1:nsdf) %dopar%
{
  library(plot3D)
  png(paste("./outputs/figures/3_psd_rbf_interpolation/psd_",i,"_data.png",sep=""),width=9,height=9,units="in",res=300)
  scatter3D(data_points[,2],data_points[,3],data_points[,4],colvar=as.numeric(psd_data[,i]),scale=FALSE,expand=1,pch=20,cex=1,phi=40,theta=0,axes=FALSE,colkey= list(length=0.7,width=1.5,cex.axis=1),zoom=0.85,clim=c(min(psd_data[,i],sdf[[i]][,4]),max(psd_data[,i],sdf[[i]][,4])))
  dev.off()
}
closeAllConnections()


closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(NbCores)
registerDoParallel(Cl)
foreach(i=1:nsdf) %dopar%
{
  library(plot3D)
  png(paste("./outputs/figures/3_psd_rbf_interpolation/psd_",i,"_interpolated.png",sep=""),width=7,height=7,units="in",res=300)
  scatter3D(grid_points[,1],grid_points[,2],grid_points[,3],colvar=as.numeric((sdf[[i]][1:nrow(grid_points),4])),scale=FALSE,expand=1,pch=20,cex=1,phi=40,theta=0,axes=FALSE,colkey= list(length=0.7,width=1.5,cex.axis=1),zoom=0.85,clim=c(min(psd_data[,i],sdf[[i]][,4]),max(psd_data[,i],sdf[[i]][,4])))
  dev.off()
}
closeAllConnections()


closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(NbCores)
registerDoParallel(Cl)
foreach(i=1:nsdf) %dopar%
{
  library(plot3Drgl)
  scatter3Drgl(data_points[,2],data_points[,3],data_points[,4],colvar=as.numeric(psd_data[,i]),scale=FALSE,expand=1,pch=20,cex=1,phi=40,theta=10,colkey= list(length=0.7,width=1.5,cex.axis=1),zoom=0.85,clim=c(min(psd_data[,i],sdf[[i]][,4]),max(psd_data[,i],sdf[[i]][,4])))
  writeWebGL(dir="./outputs/figures/3_psd_rbf_interpolation/", filename =paste("./outputs/figures/3_psd_rbf_interpolation/psd_",lithology_names[i],"_data.html"),width = 2000, reuse = TRUE)
}
closeAllConnections()


closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(NbCores)
registerDoParallel(Cl)
foreach(i=1:nsdf) %dopar%
{
  library(plot3Drgl)
  scatter3Drgl(grid_points[,1],grid_points[,2],grid_points[,3],colvar=as.numeric((sdf[[i]][1:nrow(grid_points),4])),scale=FALSE,expand=1,pch=20,cex=1,colkey= list(length=0.7,width=1.5,cex.axis=1),zoom=0.85,clim=c(min(psd_data[,i],sdf[[i]][,4]),max(psd_data[,i],sdf[[i]][,4])))
  writeWebGL(dir="./outputs/figures/3_psd_rbf_interpolation/", filename =paste("./figures/3_psd_rbf_interpolation/psd_",lithology_names[i],"_interpolated.html"),width = 2000, reuse = TRUE)
}
closeAllConnections()



## DOMAINING
PSDF=matrix(NA,nrow=nrow(grid_points),ncol=nsdf)
for ( k in 1:nsdf)
{
  PSDF[,k]=sdf[[k]][1:nrow(grid_points),4]
}

PSDF_data=matrix(NA,nrow=nrow(data_points),ncol=nsdf)
for ( k in 1:nsdf)
{
  PSDF_data[,k]=sdf[[k]][-1:-nrow(grid_points),][1:nrow(data_points),4]
}


for ( k in 1:nsdf)
{
  PSDF[,k]=sdf[[k]][1:nrow(grid_points),4]
}


domain=R_truncation_rule_domaining(PSDF)
fwrite(cbind(grid_points,domain),"./outputs/data/rbf/lithology_3D_block_model_rbf.csv",sep=" ",row.names = FALSE)


##RBF BLOCK MODEL VISUALIZATION
R_plot_block_model_lithology_1(domain,data_points)
writeWebGL(dir="./outputs/figures/3_psd_rbf_interpolation/", filename =paste("./outputs/figures/3_psd_rbf_interpolation/lithology_3D_block_model_rbf_1.html"),width = 200, reuse = TRUE)


T2=Sys.time()
difftime(T2,T1)
