T1=Sys.time()
##########################################################################################################################################
#                                                                                                                                        #
#                                       SYNTHETIC CASE STUDY: RESIDUAL MODEL REALIZATIONS (FORMATION 1)                                  #
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
library(dplyr)
library(Hmisc)
library(feather)


##CREATING SOME FOLDERS
if(!file.exists("./outputs/data/unconditional_realizations")){dir.create("./outputs/data/unconditional_realizations")}


## DATA LOADING
data_points=as.data.frame(fread("./inputs/data/drillhole_data.csv"))
contact_points=as.data.frame(fread("./inputs/data/contact_points.csv"))
grid_points=as.data.frame(fread("./inputs/data/grid_points.csv"))


contact_points_cat=list()
contact_points_cat[[1]]=readRDS('./inputs/data/contact_points_cat.RDS')[[1]]

data_points_binary=as.data.frame(fread("./outputs/data/drillhole_data_binary.csv",header=TRUE))
psd_data=as.data.frame(fread("./outputs/data/psd_data.csv")[,1])
lithology_names=1:4
nsdf=length(lithology_names)


full_grid_points=rbindlist(list(grid_points[,1:3],data_points[,2:4]))


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


## LOADING PPO RESULTS##
sill=readRDS("./outputs/data/ppo/PPO_Optimal_Sill.RDS")
range=readRDS("./outputs/data/ppo/PPO_Optimal_Range.RDS")

nb_ppo_simu=75
nb_simu=1000


sample_simu=rep(seq(1,nb_ppo_simu,by=1),ceiling(nb_simu/nb_ppo_simu))
write_feather(as.data.frame(sample_simu),"./outputs/data/unconditional_realizations/sampling_order_1.feather")
neighborhood=neigh.create(ndim=3,type=0)
db_grid=db.create(x1=full_grid_points[,1],x2=full_grid_points[,2],x3=full_grid_points[,3])

closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(NbCores,type="FORK")
registerDoParallel(Cl)

T1=Sys.time()
residual_field=foreach(l=1:nb_simu) %dopar%
{
  library(RGeostats)
  vario_model=model.create(vartype=3,range= range[sample_simu[l],1],sill=sill[sample_simu[l],1],ndim=3)
  seed_set=sample(1:1234567890,1)
  simtub(dbout=db_grid,model=vario_model,neigh=neighborhood,nbsimu=1,mean=0,seed=seed_set)[,5]
}
closeAllConnections()
write_feather(rbindlist(list(residual_field)),"./outputs/data/unconditional_realizations/residual_field_1.feather")
T2=Sys.time()
difftime(T2,T1)
