############################################################################################
#                                                                                          #
#                            DATA VISUALIZATION FUNCTIONS                                  #
#                                                                                          #
############################################################################################

R_plot_drillhole_lithology=function(lith_data,col_set=colors()[c(90,84,43,257)])
{
  categories=names(table(lith_data$FORMATION))
  col_borehole_id=rep(NA,nrow(lith_data))
  
  for ( i in 1: length(categories))
  {
    col_borehole_id[which(lith_data$FORMATION==categories[i])]=col_set[i]
  }
  
  rgl.open()
  par3d(windowRect = 50 + c( 0, 0, 800, 800) )
  rgl.bg(color="white")
  rgl.viewpoint(theta = 0, phi = -40, fov =60,zoom=0.90)

  rgl.points(lith_data[,2],lith_data[,3],lith_data[,4], color = col_borehole_id, size=4,zlim=c(min(extent_z)-10*cell_size_z,max(extent_z)+10*cell_size_z),xlim=c(min(extent_x)-10*cell_size_x,max(extent_x)+10*cell_size_x),ylim=c(min(extent_y)-10*cell_size_y,max(extent_y)+10*cell_size_y))
  box3d(col=1,size=1)
}

R_plot_drillhole_lithology_with_contacts=function(lith_data,contact_points,col_set=colors()[c(90,84,43,257)])
{
  categories=names(table(lith_data$FORMATION))
  col_borehole_id=rep(NA,nrow(lith_data))
  
  for ( i in 1: length(categories))
  {
    col_borehole_id[which(lith_data$FORMATION==categories[i])]=col_set[i]
  }
  
  rgl.open()
  par3d(windowRect = 50 + c( 0, 0, 800, 800) )
  rgl.bg(color="white")
  rgl.viewpoint(theta = 0, phi = -40, fov =60,zoom=0.90)
  rgl.points(lith_data[,2],lith_data[,3],lith_data[,4], color = col_borehole_id, size=4)
  rgl.points(contact_points[,2],contact_points[,3],contact_points[,4], color = "black",size=10)
  box3d(col=1,size=1)
}


R_plot_block_model_lithology_1=function(domain,data_points,col_set=colors()[c(90,84,43,257)])
{
  categories=names(table(data_points$FORMATION))
  col_borehole_id=rep(NA,nrow(data_points))
  
  for ( i in 1: length(categories))
  {
    col_borehole_id[which(data_points$FORMATION==categories[i])]=col_set[i]
  }
  
  
  ii2=which(domain==2)
  ii3=which(domain==3)
  ii4=which(domain==4)
  
  rgl.open()
  par3d(windowRect = 50 + c( 0, 0, 600, 600) )
  rgl.bg(color="white")
  rgl.viewpoint(theta = 0, phi = -40, fov =60,zoom=1)
  rgl.points(grid_points[,1], grid_points[,2],grid_points[,3],color ="white", size=1,xlim=extent_x,ylim=extent_y,zlim=extent_z,alpha=0)
  pch3d(grid_points[ii2,1], grid_points[ii2,2],grid_points[ii2,3],color = col_set[2],bg=col_set[2],add=TRUE,alpha=0.5,pch=15,cex=0.3)
  pch3d(grid_points[ii3,1], grid_points[ii3,2],grid_points[ii3,3],color = col_set[3],bg=col_set[3],add=TRUE,alpha=0.5,pch=15,cex=0.3)
  pch3d(grid_points[ii4,1], grid_points[ii4,2],grid_points[ii4,3],color = col_set[4],bg=col_set[4],add=TRUE,alpha=0.5,pch=15,cex=0.3)
  rgl.points(data_points[,2],data_points[,3],data_points[,4],color=col_borehole_id, size=5,xlim=extent_x,ylim=extent_y,zlim=extent_z,add=TRUE)
  box3d(col=1,size=1)
  
}

R_plot_block_model_lithology_2=function(domain,data_points,col_set=colors()[c(90,84,43,257)])
{
  categories=names(table(data_points$FORMATION))
  col_borehole_id=rep(NA,nrow(data_points))
  
  for ( i in 1: length(categories))
  {
    col_borehole_id[which(data_points$FORMATION==categories[i])]=col_set[i]
  }
  
  
  ii2=which(domain==2)
  ii3=which(domain==3)

  
  rgl.open()
  par3d(windowRect = 50 + c( 0, 0, 600, 600) )
  rgl.bg(color="white")
  rgl.viewpoint(theta = 0, phi = -40, fov =60,zoom=1)
  rgl.points(grid_points[,1], grid_points[,2],grid_points[,3],color ="white", size=1,xlim=extent_x,ylim=extent_y,zlim=extent_z,alpha=0)
  pch3d(grid_points[ii2,1], grid_points[ii2,2],grid_points[ii2,3],color = col_set[2],bg=col_set[2],add=TRUE,alpha=0.5,pch=15,cex=0.3)
  pch3d(grid_points[ii3,1], grid_points[ii3,2],grid_points[ii3,3],color = col_set[3],bg=col_set[3],add=TRUE,alpha=0.5,pch=15,cex=0.3)
  rgl.points(data_points[,2],data_points[,3],data_points[,4],color=col_borehole_id, size=5,xlim=extent_x,ylim=extent_y,zlim=extent_z,add=TRUE)
  box3d(col=1,size=1)
  
}
