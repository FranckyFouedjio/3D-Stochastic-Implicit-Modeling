##############################################################################################################
#                                                                                                            #
#           RBF SMOOTHING/INTERPOLATING of PSEUDO SIGNED DISTANCE VALUES FUNCTION                            #
#                                                                                                            #
##############################################################################################################

## REQUIRING PACKAGES LOADING
library(reticulate)
library(FNN)
library(doParallel)

## PYTHON ENVIRONMENT SETTING
use_virtualenv("r-reticulate")
numpy=import("numpy")
scipy=import("scipy")
mp=import("multiprocessing")
source_python("./code/functions/rbf_prediction.py")


R_rbf_smoothing_interpolating=function(data,grid,smooth_param=0,type)
  
{
  scale=mean(get.knnx(data[,1:3],data[,1:3], k=2,algorithm=c("kd_tree"))$nn.dist[,2])#Minimum distance between adjacent nodes.
  
  x1=numpy$array(data[,1])
  y1=numpy$array(data[,2])
  z1=numpy$array(data[,3])
  t1=numpy$array(data[,4])
  
  x2=numpy$array(grid[,1])
  y2=numpy$array(grid[,2])
  z2=numpy$array(grid[,3])
  
  npart=500
  ilength=nrow(grid)%/%npart
  
  gridl=list()
  for ( k in 1:(npart-1))
  {
    gridl[[k]]=grid[((k-1)*ilength +1):(k*ilength),]
  }
  gridl[[npart]]=grid[(ilength*(npart-1) +1):(nrow(grid)),]
  
  if(type==0){rbf1=scipy$interpolate$Rbf(x1,y1,z1,t1,smooth=smooth_param,epsilon=scale,`function`='cubic')}
  if(type==1){rbf1=scipy$interpolate$Rbf(x1,y1,z1,t1,smooth=smooth_param,epsilon=scale,`function`='linear')}
 
  
  My_List=list()
  for ( i in 1: (npart))
  {
    My_List[[i]]=list()
    My_List[[i]][[1]]=rbf1
    My_List[[i]][[2]]=numpy$array(cbind(gridl[[i]][,1],gridl[[i]][,2],gridl[[i]][,3]))
  }
  
  My_Listp=r_to_py(My_List)
  
  pool=mp$Pool(numpy$int(15))
  d1=pool$starmap(rbf_predict,My_Listp)
  pool$terminate()
  
  return(cbind(grid,unlist(d1)))
}

