T1=Sys.time()
##########################################################################################################################################
#                                                                                                                                        #
#                           SYNTHETIC CASE STUDY: PROBABILITY PERTUBATION OPTIMIZATION/GRADUAL DEFORMATION OPTIMIZATION                  #
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
library(pracma)
library(mco)
library(FNN)
library(data.table)
library(parallel)


## SOURCE CODE LOADING
source("./code/functions/global_ppo_computing_function.R")
source("./code/functions/dichotomization_function.R")
source("./code/functions/pseudo_signed_distance_function.R")

##CREATING SOME FOLDERS
if(!file.exists("./outputs/figures/4_ppo")){dir.create("./outputs/figures/4_ppo")}
if(!file.exists("./outputs/data/ppo")){dir.create("./outputs/data/ppo")}
if(!file.exists("./outputs/data/ppo/runs")){dir.create("./outputs/data/ppo/runs")}



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

psd_data_capped=R_psd_capping(psd_data,threshold=50)
db_data=foreach(i=1:nsdf) %dopar%
{
  library(RGeostats)
  db.create(x1=data_points[,2],x2=data_points[,3],x3=data_points[,4],z1=psd_data_capped[,i])
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
  model.auto(vario_exp[[i]],struct = melem.name(c(3)),upper=c(paste("M1V=",var(psd_data_capped[,i]),sep=""),paste("M1R=",200,sep="")),auth.aniso =FALSE,auth.locksame=FALSE,draw=FALSE,flag.noreduce=TRUE)
}


foreach(i=1:nsdf) %dopar% 
{
  png(paste("./outputs/figures/4_ppo/vario_signed_distance_capped_",lithology_names[i],".png",sep=""),width=5,height=5,units="in",res=300)
  vario.plot(vario_exp[[i]],main="",reset=TRUE,flag.norm = FALSE,cex.axis=1.5,cex.lab=1.5,lty=2,lwd=2,npairdw = FALSE,varline =TRUE,xlab="Distance ",ylab="Variogram")
  model.plot(vario_model[[i]],vario=vario_exp[[i]],add=T,cex.axis=1.5,lwd=c(2),lty=c(1),col=2)
  legend("bottomright",legend=c("Experimental Variogram", "Fitted Variogram"),lty=c(1,1),cex=1.2, col=c(1,2), lwd=c(1.5,1.5))
  dev.off()
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





## GLOBAL PPO/GDO PERFORMING 
closeAllConnections()
NbCores=15#detectCores()
Cl=makeCluster(NbCores)
registerDoParallel(Cl)
niter=100
r=1
path="./outputs/data/ppo/runs/"
T1=Sys.time()
Res=R_multiple_run_global_PPO(full_grid_points,grid_points,cell_size_x,cell_size_y,cell_size_z,phi_g_fields_list,psd_data,min_sill,max_sill,min_range,max_range,vol_min,vol_max,niter,path,r)
T2=Sys.time()
difftime(T2,T1)
closeAllConnections()
