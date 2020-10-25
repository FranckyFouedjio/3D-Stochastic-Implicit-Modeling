T1=Sys.time()
###############################################################################################
#                                                                                             #
#           SYNTHETIC CASE STUDY: CONDITIONAL MODEL GENERATION (FORMATION 1)                  #
#                                                                                             #
###############################################################################################


## WORKING DIRECTORY SETTING
setwd("/media/fouedjio/095e0241-4724-4a1f-84f8-9ddda0df2da9/fouedjio/3d_stochastic_implicit_modeling/")

## Packages Loading
library(factoextra)
library(data.table)
library(foreach)
library(doParallel)
library(plot3Drgl)
library(plot3D)
library(pracma)
library(FNN)
library(ggfortify)
library(polycor)
library(psych)
library(philentropy)
library(MASS)
library(ggplot2)
library(misc3d)
library(CCA)
library(RGeostats)
library(rockchalk)
library(Rglpk)
library(mvtnorm)
library(LowRankQP)
library(Rfast)
library(feather)


##CREATING SOME FOLDERS
if(!file.exists("./outputs/data/conditional_model_generation")){dir.create("./outputs/data/conditional_model_generation")}
if(!file.exists("./outputs/figures/9_conditional_model_generation")){dir.create("./outputs/figures/9_conditional_model_generation")}
if(!file.exists("./outputs/figures/9_conditional_model_generation/1")){dir.create("./outputs/figures/9_conditional_model_generation/1")}
if(!file.exists("./outputs/figures/9_conditional_model_generation/1/scores_rqp")){dir.create("./outputs/figures/9_conditional_model_generation/1/scores_rqp")}


## DATA LOADING
data=data_points=as.data.frame(fread("./inputs/data/drillhole_data.csv"))
contact_points=as.data.frame(fread("./inputs/data/contact_points.csv"))
grid=grid_points=as.data.frame(fread("./inputs/data/grid_points.csv"))


contact_points_cat=list()
contact_points_cat[[1]]=readRDS('./inputs/data/contact_points_cat.RDS')[[1]]

data_points_binary=as.data.frame(fread("./outputs/data/drillhole_data_binary.csv",header=TRUE))
psd_data=as.data.frame(fread("./outputs/data/psd_data.csv")[,1])
lithology_names=1:4
nsdf=length(lithology_names)


full_grid_points=rbindlist(list(grid_points[,1:3],data_points[,2:4]))
full_grid=rbindlist(list(grid_points[,1:3],data_points[,2:4]))


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



nsimu=1000



## Prior distribution of PC SCORES (MEAN AND COVARIANCE MATRIX)
pca_prior_sdf=readRDS(paste("./outputs/data/pca_unconditional_realizations/pca_unconditional_sdf_",1,".rds",sep=""))
model_score=pca_prior_sdf$x

eigc=100*pca_prior_sdf$sdev^2/sum(pca_prior_sdf$sdev^2)
eigc1=round(eigc[1],digits=3)
eigc2=round(eigc[2],digits=3)
eigc3=round(eigc[3],digits=3)

alpha_PRIOR=as.matrix(pca_prior_sdf$x)
m_PRIOR=rep(0,nsimu)
c_PRIOR=diag(pca_prior_sdf$sdev^2)


## INEQUALITY CONSTRAINTS FORMULATION
ncat=1
cat_pca_eig_prior_sdf=list()
m_scale=list()
s_scale=list()
for ( j in 1:ncat)
{
  cat_pca_eig_prior_sdf[[j]]=pca_prior_sdf$rotation[((j-1)*nrow(full_grid)+1):(j*nrow(full_grid)),]
  m_scale[[j]]=pca_prior_sdf$center[((j-1)*nrow(full_grid)+1):(j*nrow(full_grid))]
  s_scale[[j]]=pca_prior_sdf$scale[((j-1)*nrow(full_grid)+1):(j*nrow(full_grid))]
}

cat_pca_eig_prior_sdf_data=list()
m_scale_data=list()
s_scale_data=list()
for ( j in 1:ncat)
{
  cat_pca_eig_prior_sdf_data[[j]]=cat_pca_eig_prior_sdf[[j]][(nrow(grid) +1):(nrow(grid)+nrow(data)),]
  m_scale_data[[j]]=m_scale[[j]][(nrow(grid) +1):(nrow(grid)+nrow(data))]
  s_scale_data[[j]]=s_scale[[j]][(nrow(grid) +1):(nrow(grid)+nrow(data))]
}


rhs_l=list()
for ( j in 1:ncat)
{
  rhs_l[[j]]=matrix(rep(m_scale_data[[j]],nsimu),nrow=nsimu,byrow=T)/(matrix(rep(s_scale_data[[j]],nsimu),nrow=nsimu,byrow=T)) 
}


matcoef_l=list()
for ( j in 1:ncat)
{
  matcoef_l[[j]]=cat_pca_eig_prior_sdf_data[[j]]
}
mat_coef=matcoef_l[[1]]

ge=matrix(NA,ncol=ncat,nrow=nrow(data))
for(j in 1:ncat)
{
  ge[,j]=ifelse(as.matrix(psd_data[,j])<0,"<",">")
}



ge=c(ge[,1])


rm(cat_pca_eig_prior_sdf)
gc()

rm(matcoef_l)
gc()


## INITIAL SOLUTION BY LP
ind_select=which(ge==">")
upper_bounds_m=upper_bounds=rhs=-c(rhs_l[[1]][1,])
upper_bounds_m[ind_select]=-upper_bounds[ind_select]

ge_m=ge
ge_m[ind_select]="<"

mat_coef_m=mat_coef
mat_coef_m[ind_select,]=-mat_coef[ind_select,]

bounds=list(lower = list(ind = seq(1,nsimu), val = rep(-Inf,nsimu)), upper = list(ind = seq(1,nsimu), val = rep(Inf,nsimu)))
par_init_model_score=Rglpk_solve_LP(obj = numeric(nsimu), mat = mat_coef_m, dir = ge_m, rhs =upper_bounds_m-1e-6,bounds=bounds,max=TRUE)$solution
ff=mat_coef_m%*%par_init_model_score
table(ff<=(upper_bounds_m))


m_PRIOR=rep(0,nsimu)#colMeans(alpha_PRIOR)
c_PRIOR=diag(pca_prior_sdf$sdev^2)#cov(alpha_PRIOR)

## RANDOMIZED QUADRATIC PROGRAMMING

closeAllConnections()
NbCores=detectCores()
Cl=makeCluster(NbCores)
registerDoParallel(Cl)

m_PRIOR=rep(0,nsimu)
ii=which(pca_prior_sdf$sdev==0)
pca_prior_sdf$sdev[ii]=.Machine$double.eps
c_PRIOR=diag(pca_prior_sdf$sdev^2)
c_PRIOR_=cov(as.matrix(pca_prior_sdf$x))
c_PRIOR_inv=diag(1/pca_prior_sdf$sdev^2)
rm(pca_prior_sdf)
gc()


n_pos_simu=1000
alpha_PRIOR=t(rmvnorm(n_pos_simu,m_PRIOR,c_PRIOR))

LUR_Res=foreach(l=1:n_pos_simu,.combine = rbind)%dopar%
{
  library(quadprog)
  solve.QP(c_PRIOR_inv, alpha_PRIOR[,l]%*%c_PRIOR_inv, -t(mat_coef_m), -upper_bounds_m+1e-6)$solution
}
closeAllConnections()



utest=rep(NA,n_pos_simu)
for (i in 1:n_pos_simu)
{
  ff=mat_coef_m%*%LUR_Res[i,]
  utest[i]=table(ff<=(upper_bounds_m))[1]
}
if(sum(utest)==nrow(mat_coef_m)*n_pos_simu){cat("RQP Verified!")}



for ( l in seq(1,nsimu-1,by=1))
{
  png(paste("./outputs/figures/9_conditional_model_generation/1/scores_rqp/Con_Scores_",l,"_",l+1,".png",sep=""),width=5,height=5,units="in",res=300)
  plot(model_score[,c(l,l+1)],col=1,pch=19,xlab=paste("PC",l," (",round(eigc[l],digits=3),"% )"),ylab=paste("PC",l+1," (",round(eigc[l+1],digits=3),"% )"),xlim=c(min(LUR_Res[,l],model_score[,l]),max(LUR_Res[,l],model_score[,l])),ylim=c(min(LUR_Res[,l+1],model_score[,l+1]),max(LUR_Res[,l+1],model_score[,l+1])))
  points(LUR_Res[,c(l)],LUR_Res[,c(l+1)],col=2,pch=19)
  legend("bottomleft",legend=c("Unconditional","Conditional"),pch=c(19,19),col=c(1,2),cex=c(0.7,0.7))
  dev.off()
}

saveRDS(LUR_Res,"./outputs/data/conditional_model_generation/RQP_1.RDS")


## RECONSTRUCTION LUR POSTERIOR LEVEL SET
n_pos_simu=nrow(LUR_Res)
pca_prior_sdf=readRDS(paste("./outputs/data/pca_unconditional_realizations/pca_unconditional_sdf_",1,".rds",sep=""))
recon_=as.matrix(LUR_Res)%*% t(as.matrix(pca_prior_sdf$rotation))
a=1/pca_prior_sdf$scale
b=-pca_prior_sdf$center
rm(pca_prior_sdf)
gc()
recon_=scale(recon_, center=FALSE, scale=a)
recon_=scale(recon_, center=b, scale=FALSE)
signed_distance_fields_posterior_=t(recon_)
rm(recon_)
gc()

write_feather(as.data.frame(signed_distance_fields_posterior_),paste("./outputs/data/conditional_model_generation/conditional_sdf_model_",1,".feather",sep=""))
T2=Sys.time() 
difftime(T2,T1)

