###############################################################
#                                                             #
#       PSEUDO SIGNED DISTANCE COMPUTING FUNCTION             #
#                                                             #
###############################################################

R_pseudo_signed_distance=function(drillhole,contact_points_cat)
{
  
  category=names(table(drillhole[,5]))
  ncat=length(category)
  pseudo_signed_distance=matrix(NA,ncol=ncat,nrow=nrow(drillhole))
  
  for (k in 1:ncat)
  {
    pseudo_signed_distance[,k]=(1-2*as.numeric(drillhole[,5]==category[k]))*get.knnx(contact_points_cat[[k]][,2:4],drillhole[,2:4],k=1,algorithm=c("kd_tree"))$nn.dist[,1]
  }
  return(pseudo_signed_distance)
  
}


R_pseudo_signed_distance_p=function(drillhole,contact_points_cat)
{
  
  category=names(table(drillhole[,5]))
  ncat=length(category)
  
  pseudo_signed_distance=foreach(k=1:ncat,.combine = cbind) %dopar%
  {
    library(FNN)
   (1-2*as.numeric(drillhole[,5]==category[k]))*get.knnx(contact_points_cat[[k]][,2:4],drillhole[,2:4],k=1,algorithm=c("kd_tree"))$nn.dist[,1]
  }
  return(pseudo_signed_distance)
  
}

R_psd_capping=function(psd,threshold)
{
  psd_c=psd
  for ( i in 1:ncol(psd))
  {
    i_sel=which(abs(psd[,i])>threshold)
    psd_c[i_sel,i]=sign(psd[i_sel,i])*threshold
  }
  return(psd_c)
}