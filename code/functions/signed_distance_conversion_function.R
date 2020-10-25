R_signed_distance=function(grid,full_grid_points,full_grid_rbf)
{
  nx=length(unique(grid[,1]))
  ny=length(unique(grid[,2]))
  nz=length(unique(grid[,3]))
  
  gridx=unique(grid[,1])
  gridy=unique(grid[,2])
  gridz=unique(grid[,3])
  
  iso_phi=unique(createisosurf(x=gridx, y=gridy, z=gridz,colvar=array(full_grid_rbf[1:nrow(grid)],dim=c(nx,ny,nz)),level=0))
  full_sdf=sign(full_grid_rbf)*get.knnx(iso_phi[,1:3],as.matrix(full_grid_points[,1:3]), k=1,algorithm=c("kd_tree"))$nn.dist[,1]
  
  return(full_sdf)
  
}

