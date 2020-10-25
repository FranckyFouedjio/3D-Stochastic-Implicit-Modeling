##############################################################################################################
#                                                                                                            #
#                    GLOBAL PROBABILITY PERTUBATION/GRADUAL DEFORMATION  OPTIMIZATION FUNCTION               #
#                                                                                                            #
##############################################################################################################

## REQUIRING PACKAGES LOADING
library(foreach)
library(doParallel)
library(reticulate)
library(FNN)
library(plot3D)
library(Rvcg)
library(stats)
library(data.table)
library(parallel)


## LOSS/OBJECTIVE FUNCTION DEFINITION
R_loss_function_1=function(y,f)
{
  return(sum(log(1+exp(-y*f)))) #MisMatch at Sampling Points
}

R_loss_function_2=function(f)
{
  return(sum(f**2))    #Mean Squared Error at Contact Points
}

R_loss_function_3=function(f)
{
  return((sum(f))**2) #Squared Mean Error at Contact Points
}

R_loss_function_4=function(vp,a,b) 
{
  return((vp-a)**2 + (vp-b)**2)
}

R_loss_function_5=function(vp,a,b) 
{
  return(as.numeric((vp>b)|(vp<a)))
}

R_loss_function_6=function(vp,a,b) 
{
  l1=R_loss_function_4(vp,a,b)
  l2=R_loss_function_5(vp,a,b)
  return(l1*l2)
}


R_truncation_rule_domaining_single=function(sd_vector)
{
  if(sd_vector[1]<0){domain=1}
  else if (sd_vector[2]<0){domain=2}
  else if (sd_vector[3]<0){domain=3}else{domain=4}
  return(domain)
}

R_truncation_rule_domaining=function(sd_matrix)
{
  return(apply(sd_matrix,1,R_truncation_rule_domaining_single))
}



## WEIGHTS COMPUTING FOR THE OBJECTIVE FUNCTION
R_weight_computing=function(r,target_full_grid_l,target_grid,cellx,celly,cellz,v_current,v,phi_current_l,phi_obs_m,min_sill_l,max_sill_l,min_range_l,max_range_l,vol_min_l,vol_max_l)
{
  
  library(foreach)
  library(doParallel)
  library(reticulate)
  library(FNN)
  library(plot3D)
  library(Rvcg)
  library(stats)
  library(data.table)
  
  
  nx=length(unique(target_grid[,1]))
  ny=length(unique(target_grid[,2]))
  nz=length(unique(target_grid[,3]))
  
  extent_x=range(target_grid[,1])
  extent_y=range(target_grid[,2])
  extent_z=range(target_grid[,3])
  
  n_target_grid=nrow(target_grid)
  n_phi_obs_m=nrow(phi_obs_m)
  n_target_full_grid_l=length(target_full_grid_l)
  
  
  gridx=unique(target_grid[,1])
  gridy=unique(target_grid[,2])
  gridz=unique(target_grid[,3])
  
  iso_phi_current_l=list()
  phi_current_sdf_l=list()
  closest_l=list()
  
  for ( l in 1:n_target_full_grid_l)
  {
    iso_phi_current_l[[l]]=unique(createisosurf(x=gridx, y=gridy, z=gridz,colvar=array(phi_current_l[[l]][1:n_target_grid],dim=c(nx,ny,nz)),level=0))
    closest_l[[l]]=get.knnx(iso_phi_current_l[[l]],target_full_grid_l[[l]], k=1,algorithm=c("kd_tree"))
    phi_current_sdf_l[[l]]=sign(phi_current_l[[l]])*closest_l[[l]]$nn.dist[,1]
  }
  
  
  W=list()
  W$xrange=c(min(target_grid[,1]),max(target_grid[,1]))
  W$yrange=c(min(target_grid[,2]),max(target_grid[,2]))
  W$zrange=c(min(target_grid[,3]),max(target_grid[,3]))
  
  mygrid <- grid.prep(W = W, M = nx%/%3, N = ny%/%3, L=nz%/%3, ext=2)
  cent <- expand.grid(mygrid$mcens, mygrid$ncens,mygrid$lcens)
  Rx <- mygrid$M.ext * mygrid$cell.width
  Ry <- mygrid$N.ext * mygrid$cell.height
  Rz <- mygrid$L.ext * mygrid$cell.length
  m.abs.diff.row1 <- abs(mygrid$mcens.ext[1] - mygrid$mcens.ext)
  m.diff.row1 <- pmin(m.abs.diff.row1, Rx - m.abs.diff.row1)
  n.abs.diff.row1 <- abs(mygrid$ncens.ext[1] - mygrid$ncens.ext)
  n.diff.row1 <- pmin(n.abs.diff.row1, Ry - n.abs.diff.row1)
  l.abs.diff.row1 <- abs(mygrid$lcens.ext[1] - mygrid$lcens.ext)
  l.diff.row1 <- pmin(l.abs.diff.row1, Rz - l.abs.diff.row1)
  cent.ext.row1 <- expand.grid(m.diff.row1, n.diff.row1, l.diff.row1)
  D.ext.row1 <- array(sqrt(cent.ext.row1[, 1]^2 + cent.ext.row1[, 2]^2+ cent.ext.row1[, 3]^2), c(mygrid$M.ext,mygrid$N.ext,mygrid$L.ext))
  
  fftma_closest=list()
  for ( l in 1:n_target_full_grid_l)
  {
    fftma_closest[[l]]=get.knnx(cent[,1:3],iso_phi_current_l[[l]], k=1,algorithm=c("kd_tree"))$nn.index
  }
  
  
  phi_sim=list()
  for(l in 1:n_target_full_grid_l)
  {
    sill_current=(min_sill_l[[l]] + pnorm(v_current[[l]][(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)+1]*cos(r[l]) + v[[l]][(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)+1]*sin(r[l]))*(max_sill_l[[l]]-min_sill_l[[l]]))[1]
    range_current=(min_range_l[[l]]+ pnorm(v_current[[l]][(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)+2]*cos(r[l]) + v[[l]][(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)+2]*sin(r[l]))*(max_range_l[[l]]-min_range_l[[l]]))[1]
    
    v_updated=v_current[[l]][1:(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)]*cos(r[l]) + v[[l]][1:(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)]*sin(r[l])
    r_sim=fft_ma_3d_optimal(D.ext.row1,M=nx%/%3,N=ny%/%3,L=nz%/%3,sill=sill_current,range=range_current,std=v_updated)[fftma_closest[[l]]]
    phi_sim[[l]]=phi_current_sdf_l[[l]] + r_sim[closest_l[[l]]$nn.index]
  } 
  
  phi_sim_grid=matrix(NA,nrow=n_target_grid,ncol=n_target_full_grid_l)
  
  for ( l in 1:n_target_full_grid_l)
  {
    phi_sim_grid[,l]=phi_sim[[l]][1:n_target_grid]
  }
  
  
  domain=R_truncation_rule_domaining(phi_sim_grid)
  
  vp=rep(0,n_target_full_grid_l)
  for ( l in 1:(n_target_full_grid_l))
  {
    vp[l]=length(which(domain==l))/n_target_grid
  }
  
  
  U1=matrix(NA,nrow=n_target_full_grid_l,ncol=3)
  for ( l in 1:n_target_full_grid_l)
  {
    U1[l,1]=R_loss_function_1(sign(phi_obs_m[,l]),phi_sim[[l]][-1:-n_target_grid][1:n_phi_obs_m])
    U1[l,2]=R_loss_function_2(phi_sim[[l]][-1:-n_target_grid][-1:-n_phi_obs_m]) 
    U1[l,3]=R_loss_function_3(phi_sim[[l]][-1:-n_target_grid][-1:-n_phi_obs_m]) 
  }
  

  U2=rep(NA,n_target_full_grid_l)
  for ( l in 1:(n_target_full_grid_l))
  {
    U2[l]=R_loss_function_6(vp[l],vol_min_l[[l]],vol_max_l[[l]]) 
  }
  
  return(list(U1=U1,U2=U2))
}



R_global_weight_computing=function(target_full_grid_l,target_grid,cellx,celly,cellz,v_current,v,phi_current_l,phi_obs_m,min_sill_l,max_sill_l,min_range_l,max_range_l,vol_min_l,vol_max_l)
{
  
  n_target_grid=nrow(target_grid)
  n_phi_obs_m=nrow(phi_obs_m)
  n_target_full_grid_l=length(target_full_grid_l)
  
  
  library(plot3D)
  library(Rvcg)
  library(FNN)
  library(stats)
  library(data.table)
  
  source("./code/functions/fftma_function.R")
  source("./code/functions/global_ppo_computing_function.R")
  
  
  EE=seq(-1,1,length.out=50)
  
  RR1=array(NA,dim=c(n_target_full_grid_l,3,length(EE)))
  RR2=matrix(NA,ncol=length(EE),nrow=n_target_full_grid_l)

  Result=foreach(k=1:length(EE))%dopar%
  {
    source("./code/functions/fftma_function.R")
    source("./code/functions/global_ppo_computing_function.R")
    R_weight_computing(rep(EE[k],n_target_full_grid_l),target_full_grid_l,target_grid,cellx,celly,cellz,v_current,v,phi_current_l,phi_obs_m,min_sill_l,max_sill_l,min_range_l,max_range_l,vol_min_l,vol_max_l)
  }
  
  for ( k in 1:length(EE))
  {
    RR1[,,k]=Result[[k]][[1]]
    RR2[,k]=Result[[k]][[2]]
  }
  WW1=apply(RR1,c(1,2),mean)
  WW2=apply(RR2,1,mean)
  WW1[which(WW1==0)]=1e-40
  WW2[which(WW2==0)]=1e-40
  WWN1=WW1[1,1]/WW1
  WWN2=WW1[1,1]/WW2
  WWN1[,3]=20*WWN1[,3]
  return(list(WWN1=WWN1,WWN2=WWN2))
  
}






## GLOBAL PROBABILITY PERTUBATION/GRADUAL DEFORMATION  OPTIMIZATION 



R_global_PPO=function(target_full_grid_l,target_grid,cellx,celly,cellz,phi_current_l,phi_obs_m,min_sill_l,max_sill_l,min_range_l,max_range_l,vol_min_l,vol_max_l,niter)
{
  
  nx=length(unique(target_grid[,1]))
  ny=length(unique(target_grid[,2]))
  nz=length(unique(target_grid[,3]))
  
  extent_x=range(target_grid[,1])
  extent_y=range(target_grid[,2])
  extent_z=range(target_grid[,3])
  
  n_target_grid=nrow(target_grid)
  n_phi_obs_m=nrow(phi_obs_m)
  n_target_full_grid_l=length(target_full_grid_l)
  
  U_r=matrix(NA,nrow=niter,ncol=n_target_full_grid_l)
  U_f=rep(NA,nrow=niter)
  sill_opt=matrix(NA,nrow=niter,ncol=n_target_full_grid_l)
  range_opt=matrix(NA,nrow=niter,ncol=n_target_full_grid_l)
  
  phi_current_list=list()
  for ( i in 1:niter)
  {
    phi_current_list[[i]]=list()
  }
  
  
  gridx=unique(target_grid[,1])
  gridy=unique(target_grid[,2])
  gridz=unique(target_grid[,3])
  
  
  iso_phi_current_l= foreach(l=1:n_target_full_grid_l)%dopar%
  {
    library(plot3D)
    unique(createisosurf(x=gridx, y=gridy, z=gridz,colvar=array(phi_current_l[[l]][1:n_target_grid],dim=c(nx,ny,nz)),level=0))
  }
  
  closest_l=foreach(l=1:n_target_full_grid_l)%dopar%
  {
    library(FNN)
    get.knnx(iso_phi_current_l[[l]],target_full_grid_l[[l]], k=1,algorithm=c("kd_tree"))
  }
  
  phi_current_sdf_l=foreach(l=1:n_target_full_grid_l)%dopar%
  {
    sign(phi_current_l[[l]])*closest_l[[l]]$nn.dist[,1]
  }
  
  
  
  W=list()
  W$xrange=c(min(target_grid[,1]),max(target_grid[,1]))
  W$yrange=c(min(target_grid[,2]),max(target_grid[,2]))
  W$zrange=c(min(target_grid[,3]),max(target_grid[,3]))
  
  mygrid <- grid.prep(W = W, M = nx%/%3, N = ny%/%3, L=nz%/%3, ext=2)
  cent <- expand.grid(mygrid$mcens, mygrid$ncens,mygrid$lcens)
  Rx <- mygrid$M.ext * mygrid$cell.width
  Ry <- mygrid$N.ext * mygrid$cell.height
  Rz <- mygrid$L.ext * mygrid$cell.length
  m.abs.diff.row1 <- abs(mygrid$mcens.ext[1] - mygrid$mcens.ext)
  m.diff.row1 <- pmin(m.abs.diff.row1, Rx - m.abs.diff.row1)
  n.abs.diff.row1 <- abs(mygrid$ncens.ext[1] - mygrid$ncens.ext)
  n.diff.row1 <- pmin(n.abs.diff.row1, Ry - n.abs.diff.row1)
  l.abs.diff.row1 <- abs(mygrid$lcens.ext[1] - mygrid$lcens.ext)
  l.diff.row1 <- pmin(l.abs.diff.row1, Rz - l.abs.diff.row1)
  cent.ext.row1 <- expand.grid(m.diff.row1, n.diff.row1, l.diff.row1)
  D.ext.row1 <- array(sqrt(cent.ext.row1[, 1]^2 + cent.ext.row1[, 2]^2+ cent.ext.row1[, 3]^2), c(mygrid$M.ext,mygrid$N.ext,mygrid$L.ext))
  
  
  fftma_closest=foreach(l=1:n_target_full_grid_l)%dopar%
  {
    library(FNN)
    get.knnx(cent[,1:3],iso_phi_current_l[[l]], k=1,algorithm=c("kd_tree"))$nn.index
  }
  
  
  v_current_0=rnorm((mygrid$M.ext*mygrid$N.ext*mygrid$L.ext) +2,mean=0,sd=1) 
  v_0=rnorm((mygrid$M.ext*mygrid$N.ext*mygrid$L.ext) +2,mean=0,sd=1) 
  
  
  v_current=list()
  v=list()
  for(l in 1:n_target_full_grid_l)
  {
    v_current[[l]]=v_current_0
    v[[l]]= v_0
  }
  
  
  v_current_l=list()
  v_l=list()
  for ( l in 1:n_target_full_grid_l)
  {
    v_current_l[[l]]=v_current[[l]]
    v_l[[l]]=v[[l]]
  }
  
  source("./code/functions/global_ppo_computing_function.R")
  weights_m=R_global_weight_computing(target_full_grid_l,target_grid,cellx,celly,cellz,v_current,v,phi_current_l,phi_obs_m,min_sill_l,max_sill_l,min_range_l,max_range_l,vol_min_l,vol_max_l)
  
  for ( i in 1:niter)
  {
    
    R_loss_function_min=function(r)
    {
      object=list()
      for(l in 1:n_target_full_grid_l)
      {
        
        object[[l]]=list()
        object[[l]][[1]]=(min_sill_l[[l]] + pnorm(v_current_l[[l]][(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)+1]*cos(r[l]) + v_l[[l]][(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)+1]*sin(r[l]))*(max_sill_l[[l]]-min_sill_l[[l]]))[1]
        object[[l]][[2]]=(min_range_l[[l]]+ pnorm(v_current_l[[l]][(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)+2]*cos(r[l]) + v_l[[l]][(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)+2]*sin(r[l]))*(max_range_l[[l]]-min_range_l[[l]]))[1]
        object[[l]][[3]]=fftma_closest[[l]]
        object[[l]][[4]]=phi_current_sdf_l[[l]]
        object[[l]][[5]]=closest_l[[l]]$nn.index
        object[[l]][[6]]=v_current_l[[l]][1:(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)]*cos(r[l]) + v_l[[l]][1:(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)]*sin(r[l])
      }
      
      fft_ma_3d_p=function(object)
      {
        r_sim=fft_ma_3d_optimal(D.ext.row1,M=nx%/%3,N=ny%/%3,L=nz%/%3,object[[1]],object[[2]],object[[6]])[object[[3]]]
        return(object[[4]]+r_sim[object[[5]]])
      }
      
      
      
      phi_sim=mclapply(object, fft_ma_3d_p, mc.cores = n_target_full_grid_l)
      
      rm(object)
      gc()
      
      phi_sim_grid=matrix(NA,nrow=n_target_grid,ncol=n_target_full_grid_l)
      
      for ( l in 1:n_target_full_grid_l)
      {
        phi_sim_grid[,l]=phi_sim[[l]][1:n_target_grid]
      }
      
      domain=R_truncation_rule_domaining(phi_sim_grid)
      
      vp=rep(0,n_target_full_grid_l)
      for (l in 1:(n_target_full_grid_l))
      {
        vp[l]=length(which(domain==l))/n_target_grid
      }
      
      
      U1=matrix(NA,nrow=n_target_full_grid_l,ncol=3)
      for ( l in 1:n_target_full_grid_l)
      {
        U1[l,1]=weights_m[[1]][l,1]*R_loss_function_1(sign(phi_obs_m[,l]),phi_sim[[l]][-1:-n_target_grid][1:n_phi_obs_m])
        U1[l,2]=weights_m[[1]][l,2]*R_loss_function_2(phi_sim[[l]][-1:-n_target_grid][-1:-n_phi_obs_m]) 
        U1[l,3]=weights_m[[1]][l,3]*R_loss_function_3(phi_sim[[l]][-1:-n_target_grid][-1:-n_phi_obs_m]) 
      }
      
      
      U2=rep(NA,n_target_full_grid_l)
      for ( l in 1:(n_target_full_grid_l))
      {
        U2[l]=weights_m[[2]][l]*R_loss_function_6(vp[l],vol_min_l[[l]],vol_max_l[[l]]) 
      }
      
      return(log(sum(U1)+sum(U2)))
    }
    
  
    r_res=nloptr(x0=rep(0,n_target_full_grid_l), eval_f=R_loss_function_min, lb = rep(-pi/2,n_target_full_grid_l), ub = rep(pi/2,n_target_full_grid_l),opts=list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-20,maxeval=50))
    
    
    ropt=r_res$solution
    U_r[i,]=ropt
    U_f[i]=r_res$objective
    
    v_0=rnorm((mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)+2,mean=0,sd=1)
    
    
    for(l in 1:n_target_full_grid_l)
    {
      v[[l]]=v_0
    }
    
    object_opt=list()
    for(l in 1:n_target_full_grid_l)
    {
      
      object_opt[[l]]=list()
      object_opt[[l]][[1]]=sill_opt[i,l]=(min_sill_l[[l]] + pnorm((v_current_l[[l]][(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)+1]*cos(ropt[l]) + v_l[[l]][(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)+1]*sin(ropt[l]))[1])*(max_sill_l[[l]]-min_sill_l[[l]]))[1]
      object_opt[[l]][[2]]=range_opt[i,l]=(min_range_l[[l]]+ pnorm((v_current_l[[l]][(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)+2]*cos(ropt[l]) + v_l[[l]][(mygrid$M.ext*mygrid$N.ext*mygrid$L.ext)+2]*sin(ropt[l]))[1])*(max_range_l[[l]]-min_range_l[[l]]))[1]
      object_opt[[l]][[3]]=fftma_closest[[l]]
      object_opt[[l]][[4]]=phi_current_sdf_l[[l]]
      object_opt[[l]][[5]]=closest_l[[l]]$nn.index
      object_opt[[l]][[6]]=v_current_l[[l]]=v_current[[l]]=v_current_l[[l]]*cos(ropt[l]) + v_l[[l]]*sin(ropt[l])
      v_l[[l]]=v[[l]]
    }
    
    
    
    fft_ma_3d_p=function(object_opt)
    {
      r_sim=fft_ma_3d_optimal(D.ext.row1,M=nx%/%3,N=ny%/%3,L=nz%/%3,object_opt[[1]],object_opt[[2]],object_opt[[6]])[object_opt[[3]]]
      return(object_opt[[4]]+r_sim[object_opt[[5]]])
    }
    
    phi_current_list[[i]]=mclapply(object_opt, fft_ma_3d_p, mc.cores = n_target_full_grid_l)
    
    rm(object_opt)
    gc()
    
  }
  return(list(U_r=U_r,U_f=U_f,sill_opt=sill_opt,range_opt=range_opt,phi_current_list=phi_current_list))
}



R_multiple_run_global_PPO=function(target_full_grid_l,target_grid,cellx,celly,cellz,phi_current_l,phi_obs_m,min_sill_l,max_sill_l,min_range_l,max_range_l,vol_min_l,vol_max_l,niter,path,r)
{
  library(foreach)
  library(doParallel)
  
  foreach(l=((r-1)*15+1):(r*15)) %dopar%
  { 
    if(!file.exists(paste(path,l,"/",sep=""))){dir.create(paste(path,l,"/",sep=""))}
    library(plot3D)
    library(Rvcg)
    library(nloptr)
    library(FNN)
    library(stats)
    library(foreach)
    library(doParallel)
    source("./code/functions/fftma_function.R")
    source("./code/functions/global_ppo_computing_function.R")
    Res=R_global_PPO(target_full_grid_l,target_grid,cellx,celly,cellz,phi_current_l,phi_obs_m,min_sill_l,max_sill_l,min_range_l,max_range_l,vol_min_l,vol_max_l,niter)
    saveRDS(Res,paste(path,l,"/Result.rds",sep=""))
    rm(Res)
    gc()
    
  }
  
}

