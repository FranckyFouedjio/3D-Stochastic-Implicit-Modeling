T1=Sys.time()
#############################################################################################################
#                                                                                                           #
#                   SYNTHETIC CASE STUDY: PCA OF UNCONDITIONAL REALIZATIONS (FORMATION 1)                   #
#                                                                                                           #
#############################################################################################################

## WORKING DIRECTORY SETTING
setwd("/media/fouedjio/095e0241-4724-4a1f-84f8-9ddda0df2da9/fouedjio/3d_stochastic_implicit_modeling/")

## REQUIRING PACKAGES LOADING
library(data.table)
library(ggplot2)
library(feather)
library(Rfast)

library(reticulate)
use_virtualenv("r-reticulate")
numpy=import("numpy")
import("gc")


## LOADING SOURCE CODE
source_python("./code/functions/full_evd_fast.py") 

##CREATING SOME FOLDERS
if(!file.exists("./outputs/data/pca_unconditional_realizations")){dir.create("./outputs/data/pca_unconditional_realizations")}
if(!file.exists("./outputs/figures/8_pca_unconditional_realizations")){dir.create("./outputs/figures/8_pca_unconditional_realizations")}
if(!file.exists("./outputs/figures/8_pca_unconditional_realizations/1")){dir.create("./outputs/figures/8_pca_unconditional_realizations/1")}

##
nsimu=1000
lithology_names=1:4


## PCA
pca_unconditional_sdf=full_evd_fast(paste("./outputs/data/unconditional_realizations/unconditional_sdf_model_",1,".feather",sep=""))
names(pca_unconditional_sdf)=c("rotation","sdev","center","scale","x")
pca_unconditional_sdf$sdev=sqrt(abs(pca_unconditional_sdf$sdev)/nsimu)
saveRDS(pca_unconditional_sdf, file= paste("./outputs/data/pca_unconditional_realizations/pca_unconditional_sdf_",1,".rds",sep=""))


eigc=100*pca_unconditional_sdf$sdev^2/sum(pca_unconditional_sdf$sdev^2)
eigc1=round(eigc[1],digits=1)
eigc2=round(eigc[2],digits=1)
eigc3=round(eigc[3],digits=1)

model_score=as.data.frame(pca_unconditional_sdf$x)

l=1
png(paste("./outputs/figures/8_pca_unconditional_realizations/1/Uncon_PC_Scores_",l,"_",l+1,".png",sep=""),width=5,height=5,units="in",res=300)
ggplot(model_score, aes(x=model_score[,l], y=model_score[,l+1])) + geom_point(size=2)  + labs(x = paste("PC",l," (",round(eigc[l],digits=1),"% )"),y =paste("PC",l+1," (",round(eigc[l+1],digits=1),"% )"))+ theme(aspect.ratio=1) #+ coord_fixed()
dev.off()

l=2
png(paste("./outputs/figures/8_pca_unconditional_realizations/1/Uncon_PC_Scores_",l,"_",l+1,".png",sep=""),width=5,height=5,units="in",res=300)
ggplot(model_score, aes(x=model_score[,l], y=model_score[,l+1])) + geom_point(size=2)  + labs(x = paste("PC",l," (",round(eigc[l],digits=1),"% )"),y =paste("PC",l+1," (",round(eigc[l+1],digits=1),"% )")) + theme(aspect.ratio=1) #+ coord_fixed()
dev.off()

l=3
png(paste("./outputs/figures/8_pca_unconditional_realizations/1/Uncon_PC_Scores_",l,"_",l+1,".png",sep=""),width=5,height=5,units="in",res=300)
ggplot(model_score, aes(x=model_score[,l], y=model_score[,l+1])) + geom_point(size=2)  + labs(x = paste("PC",l," (",round(eigc[l],digits=1),"% )"),y =paste("PC",l+1," (",round(eigc[l+1],digits=1),"% )")) +  theme(aspect.ratio=1) #+ coord_fixed()
dev.off()

l=4
png(paste("./outputs/figures/8_pca_unconditional_realizations/1/Uncon_PC_Scores_",l,"_",l+1,".png",sep=""),width=5,height=5,units="in",res=300)
ggplot(model_score, aes(x=model_score[,l], y=model_score[,l+1])) + geom_point(size=2)  + labs(x = paste("PC",l," (",round(eigc[l],digits=1),"% )"),y =paste("PC",l+1," (",round(eigc[l+1],digits=1),"% )")) +  theme(aspect.ratio=1) #+ coord_fixed()
dev.off()


l=5
png(paste("./outputs/figures/8_pca_unconditional_realizations/1/Uncon_PC_Scores_",l,"_",l+1,".png",sep=""),width=5,height=5,units="in",res=300)
ggplot(model_score, aes(x=model_score[,l], y=model_score[,l+1])) + geom_point(size=2)  + labs(x = paste("PC",l," (",round(eigc[l],digits=1),"% )"),y =paste("PC",l+1," (",round(eigc[l+1],digits=1),"% )")) + theme(aspect.ratio=1) #+ coord_fixed()
dev.off()

l=6
png(paste("./outputs/figures/8_pca_unconditional_realizations/1/Uncon_PC_Scores_",l,"_",l+1,".png",sep=""),width=5,height=5,units="in",res=300)
ggplot(model_score, aes(x=model_score[,l], y=model_score[,l+1])) + geom_point(size=2)  + labs(x = paste("PC",l," (",round(eigc[l],digits=1),"% )"),y =paste("PC",l+1," (",round(eigc[l+1],digits=1),"% )"))  +  theme(aspect.ratio=1) #+ coord_fixed()
dev.off()

l=7
png(paste("./outputs/figures/8_pca_unconditional_realizations/1/Uncon_PC_Scores_",l,"_",l+1,".png",sep=""),width=5,height=5,units="in",res=300)
ggplot(model_score, aes(x=model_score[,l], y=model_score[,l+1])) + geom_point(size=2)  + labs(x = paste("PC",l," (",round(eigc[l],digits=1),"% )"),y =paste("PC",l+1," (",round(eigc[l+1],digits=1),"% )")) + theme(aspect.ratio=1) #+ coord_fixed()
dev.off()

l=8
png(paste("./outputs/figures/8_pca_unconditional_realizations/1/Uncon_PC_Scores_",l,"_",l+1,".png",sep=""),width=5,height=5,units="in",res=300)
ggplot(model_score, aes(x=model_score[,l], y=model_score[,l+1])) + geom_point(size=2)  + labs(x = paste("PC",l," (",round(eigc[l],digits=1),"% )"),y =paste("PC",l+1," (",round(eigc[l+1],digits=1),"% )")) +  theme(aspect.ratio=1) #+ coord_fixed()
dev.off()

l=9
png(paste("./outputs/figures/8_pca_unconditional_realizations/1/Uncon_PC_Scores_",l,"_",l+1,".png",sep=""),width=5,height=5,units="in",res=300)
ggplot(model_score, aes(x=model_score[,l], y=model_score[,l+1])) + geom_point(size=2)  + labs(x = paste("PC",l," (",round(eigc[l],digits=1),"% )"),y =paste("PC",l+1," (",round(eigc[l+1],digits=1),"% )")) + theme(aspect.ratio=1) #+ coord_fixed()
dev.off()

l=10
png(paste("./outputs/figures/8_pca_unconditional_realizations/1/Uncon_PC_Scores_",l,"_",l+1,".png",sep=""),width=5,height=5,units="in",res=300)
ggplot(model_score, aes(x=model_score[,l], y=model_score[,l+1])) + geom_point(size=2)  + labs(x = paste("PC",l," (",round(eigc[l],digits=1),"% )"),y =paste("PC",l+1," (",round(eigc[l+1],digits=1),"% )")) +  theme(aspect.ratio=1) #+ coord_fixed()
dev.off()


variance_explained=as.data.frame(cbind(1:nsimu,100*cumsum(pca_unconditional_sdf$sdev^2)/sum(pca_unconditional_sdf$sdev^2)))
colnames(variance_explained)=c("n","v")

ndim_select=variance_explained[which(variance_explained[,2]>=95)[1],1]
ndim_select
rm(pca_unconditional_sdf)
gc()

png("./outputs/figures/8_pca_unconditional_realizations/1/Model_Explained_Variance_PCA.png",width=5,height=5,units="in",res=300)
ggplot(variance_explained, aes(x=n, y=v)) + geom_line(size=1) + labs(x = "Dimension",y=" Explained  Variance (%)") 
dev.off()

T2=Sys.time()
difftime(T2,T1)
