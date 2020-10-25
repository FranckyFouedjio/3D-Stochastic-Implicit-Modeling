T1=Sys.time()
##########################################################################################################################################
#                                                                                                                                        #
#                                     SYNTHETIC CASE STUDY:  PSD RBF INTERPOLATION 1 (FORMATION 1)                                       #
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
library(Rvcg)
library(scatterplot3d)
library(misc3d)
library(data.table)
library(FNN)
library(ggplot2)  
library(deldir)
library("RColorBrewer")
library(profvis)
library(DescTools)
library(geometry)



## SOURCE CODE LOADING
source("./code/functions/dichotomization_function.R")
source("./code/functions/rbf_interpolation_smoothing_function.R")
source("./code/functions/signed_distance_conversion_function.R")



## PYTHON ENVIRONMENT SETTING
library(reticulate)
use_virtualenv("r-reticulate")
numpy=import("numpy")
scipy=import("scipy")
mp=import("multiprocessing")


##CREATING SOME FOLDERS
if(!file.exists("./outputs/data/rbf")){dir.create("./outputs/data/rbf")}




## DATA LOADING
data_points=as.data.frame(fread("./inputs/data/drillhole_data.csv"))
contact_points=as.data.frame(fread("./inputs/data/contact_points.csv"))
grid_points=as.data.frame(fread("./inputs/data/grid_points.csv"))
contact_points_cat=list()
contact_points_cat[[1]]=readRDS('./inputs/data/contact_points_cat.RDS')[[1]]
data_points_binary=as.data.frame(fread("./outputs/data/drillhole_data_binary.csv",header=TRUE))
psd_data=as.data.frame(fread("./outputs/data/psd_data.csv")[,1])
lithology_names=1:4
full_grid_points=list()
full_grid_points[[1]]=rbindlist(list(grid_points[,1:3],data_points[,2:4],contact_points_cat[[1]][,2:4]))




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


## RBF SMOOTHING/INTERPOLATING OF PSEUDO SIGNED DISTANCE VALUES 
rbf_result=list()
data_p_=data_p=as.data.frame(fread(paste("./outputs/data/","data_points_rbf_",lithology_names[1],".csv",sep="")))
ii_index=which(data_points_binary[,1]==1)
ii_xyz=get.knnx(grid_points,as.matrix(data_p[ii_index,1:3]),k=1,algorithm = c("kd_tree"))$nn.index
data_p[ii_index,]=cbind(grid_points[ii_xyz,],D=data_p[ii_index,4])
ii_d=which(duplicated(data_p[,1:3]))
if(length(ii_d)>=1){data_p=data_p[-ii_d,]}
grid_p=as.data.frame(fread(paste("./inputs/data/","full_grid_points_",lithology_names[1],".csv",sep="")))
rbf_result[[1]]=R_rbf_smoothing_interpolating(data_p,grid_p,smooth_param=0,type=0)
cat("RBF Interpolation Signed Distance Function Formation 1  Completed","\n")  



## CONVERT THE LEVEL SET INTO SIGNED DISTANCE FUNCTION
sdf=list()
sdf[[1]]=R_signed_distance(grid_points,full_grid_points[[1]],rbf_result[[1]][,4])
fwrite(cbind(full_grid_points[[1]],SD=sdf[[1]]),paste("./outputs/data/rbf/","psd_rbf_sdf_",1,".csv",sep=""),row.names=FALSE)

T2=Sys.time()
difftime(T2,T1)
