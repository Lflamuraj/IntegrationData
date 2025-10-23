#clears everything on R
rm(list=ls())

#set libraries 
library(StereoMorph) 
library(geomorph)
library(geiger)
library(plotrix)
library(caper)
library(ggplot2)
library(tidyverse)


#loading in data 
#Cl 
CL <- readland.tps(file = "AllCL_400.tps", specID = "imageID")
CL_sub <- CL[-c(13),,]

writeland.tps(CL_sub, file = "CL_land_new.tps", specID = TRUE)


#pca graphs 
Y.gpa_CL<-gpagen(CL_sub)

plotAllSpecimens(Y.gpa_CL$coords)

MYDATA<-gm.prcomp(Y.gpa_CL$coords)
plot(MYDATA)

#textfor plot
#text(MYDATA$x, labels = dimnames(Y.gpa_CL$coords)[[3]])


#Man
Man <- readland.tps(file = "AllrightsideMan.tps", specID = "ID" )

#pca graphs
Y.gpa_Man<-gpagen(Man)

MYDATA2<-gm.prcomp(Y.gpa_Man$coords)
plot(MYDATA2)

#text for plot
#text(MYDATA2$x, labels = dimnames(Y.gpa_Man$coords)[[3]])


#CH
CH <- readland.tps(file = "AllCH_400.tps", specID = "imageID")
#deleting overlapping lms and semis
CH_sub <- CH[-c(8,17,18,27),,]

writeland.tps(CH_sub, file = "HYD_land_new.tps", specID = TRUE)

#visualising the lms
#define.sliders(CH_sub[,,1], write.file = "CH_sliders.csv")

#loading in curves
curves_CH<-as.matrix(read.csv("CH_sliders.csv"))
Y.gpa_CH<-gpagen(CH_sub, curves = curves_CH)
#plots all individuals across all lms
plotAllSpecimens(Y.gpa_CH$coords)



MYDATA3<-gm.prcomp(Y.gpa_CH$coords)
#pc1 and 2
plot(MYDATA3)
#text for plot
#text(MYDATA3$x, labels = dimnames(Y.gpa_CH$coords)[[3]])

#pc2 and 3
#plot(MYDATA3, axis1 = 2, axis2 = 3)


#CB1
CB1 <- readland.tps(file = "AllCB1_400.tps", specID = "imageID")
#deleting overlapping lms and semis
CB1_sub <- CB1[-c(6,15),,]

writeland.tps(CB1_sub, file = "CB1_land_new.tps", specID = TRUE)

#visualising the lms
#define.sliders(CB1_sub[,,1], write.file = "CB1_sliders.csv")

#loading in curves
curves_CB1<-as.matrix(read.csv("CB1_sliders.csv"))
Y.gpa_CB1<-gpagen(CB1_sub, curves = curves_CB1)

write.csv(Y.gpa_CB1$Csize, "CB1Csize.csv")

MYDATA4<-gm.prcomp(Y.gpa_CB1$coords)
plot(MYDATA4)

#text for plot
#text(MYDATA4$x, labels = dimnames(Y.gpa_CB1$coords)[[3]])


#CB5
PJ <- readland.tps(file = "AllrightsidePJ.tps", specID = "ID")

#pca graphs ... reduces variation to fewer dementions 
Y.gpa_PJ<-gpagen(PJ)

MYDATA4<-gm.prcomp(Y.gpa_PJ$coords)
plot(MYDATA4)


#text for plot
#text(MYDATA4$x, labels = dimnames(Y.gpa_PJ$coords)[[3]])

###ALLOMETRY CORRECTED + size corrected + PCA plots 

#Cl 
#Allometry set up 
CL_Shapes<-Y.gpa_CL$coords

gdf_CL <- geomorph.data.frame(coords = CL_Shapes, Csize = as.numeric(Y.gpa_CL$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes) ... finding out if size matters
CLAnova <- procD.lm(coords~Csize, data = gdf_CL) 
summary(CLAnova) #it does?

#plot allometry of the origical not size corrected model
plotAllometry(CLAnova, size = gdf_CL$Csize, method = "RegScore")
#view the residuals of the model 
#CLAnova$residuals


## allometry is having an effect and how we are removing the effect 

# save size-adjusted residuals an make as 3D array
CLResiduals <- arrayspecs(CLAnova$residuals, p=dim(CL_Shapes)[1], k=dim(CL_Shapes)[2])
# make allometry-free shapes
CLResiduals.adj.shape <- CLResiduals+array(Y.gpa_CL$consensus, dim(CLResiduals))

#Run PCA on new  size corrected shape data
CLResiduals_PCA.Allom <- gm.prcomp(CLResiduals.adj.shape)

#gettings pc scores
CLResiduals_PCA.Allom$x
write.csv(CLResiduals_PCA.Allom$x[,c(1,2,3)], "CLResiduals_PC1_3.csv")

#ploting the size corrected specimen
plotAllSpecimens(CLResiduals.adj.shape)

#plotting the PCA pc1 + 2
plot(CLResiduals_PCA.Allom)
#text(psPCA_CL$x, labels = dimnames(Y.gpa_CL$coords)[[3]])
text(CLResiduals_PCA.Allom$x, pos = 4, label = dimnames(CLResiduals.adj.shape)[[3]], cex=.5)

#plotting the PCA pc2 + 3
plot(CLResiduals_PCA.Allom, axis1 = 2, axis2 = 3)
text(CLResiduals_PCA.Allom$x[, 2], CLResiduals_PCA.Allom$x[, 3], pos = 4, label = dimnames(CLResiduals.adj.shape)[[3]], cex=.5)


#new lm SIZE CORRECTED
gdf_CLR <- geomorph.data.frame(coords = CLResiduals.adj.shape, Csize = as.numeric(Y.gpa_CL$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes)##
CLAnovaR <- procD.lm(coords~Csize, data = gdf_CLR) 
summary(CLAnovaR) 


##ball and stick plots

open3d()
#PC1
#plotRefToTarget(M1 = CLResiduals.adj.shape  [,,"LFxTRC5008"], M2 = CLResiduals.adj.shape [,,"LFxTRC5187"], method = "vector", mag = 2)


#Man
#Allometry set up 
Man_Shapes<-Y.gpa_Man$coords

gdf_Man <- geomorph.data.frame(coords = Man_Shapes, Csize = as.numeric(Y.gpa_Man$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes) ... finding out if size matters
ManAnova <- procD.lm(coords~Csize, data = gdf_Man) 
summary(ManAnova) #it does?


#plot allometry of the origical not size corrected model
plotAllometry(ManAnova, size = gdf_Man$Csize, method = "RegScore")
#view the residuals of the model 
#CHAnova$residuals


## allometry is having an effect and how we are removing the effect 

# save size-adjusted residuals an make as 3D array
ManResiduals <- arrayspecs(ManAnova$residuals, p=dim(Man_Shapes)[1], k=dim(Man_Shapes)[2])
# make allometry-free shapes
ManResiduals.adj.shape <- ManResiduals+array(Y.gpa_Man$consensus, dim(ManResiduals))

#Run PCA on new  size corrected shape data
ManResiduals_PCA.Allom <- gm.prcomp(ManResiduals.adj.shape)

#gettings pc scores
ManResiduals_PCA.Allom$x
write.csv(ManResiduals_PCA.Allom$x[,c(1,2,3)], "ManResiduals_PC1_3.csv")

#ploting the size corrected specimen
plotAllSpecimens(ManResiduals.adj.shape)



#plotting the PCA pc1+ 2
plot(ManResiduals_PCA.Allom)
#text(psPCA_CL$x, labels = dimnames(Y.gpa_CL$coords)[[3]])
text(ManResiduals_PCA.Allom$x, pos = 4, label = dimnames(ManResiduals.adj.shape)[[3]], cex=.5)

#plotting the PCA pc2 + 3
plot(ManResiduals_PCA.Allom, axis1 = 2, axis2 = 3)
text(ManResiduals_PCA.Allom$x[, 2], ManResiduals_PCA.Allom$x[, 3], pos = 4, label = dimnames(ManResiduals.adj.shape)[[3]], cex=.5)

#new lm SIZE CORRECTED
gdf_ManR <- geomorph.data.frame(coords = ManResiduals.adj.shape, Csize = as.numeric(Y.gpa_Man$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes)##
ManAnovaR <- procD.lm(coords~Csize, data = gdf_ManR) 
summary(ManAnovaR) 


##ball and stick plots

open3d()
#PC1
plotRefToTarget(M1 = ManResiduals.adj.shape  [,,"LFxTRC5448"], M2 = ManResiduals.adj.shape [,,"LFxTRC5084"],
                method = "vector", mag = 1)



#CH
#Allometry set up 
CH_Shapes<-Y.gpa_CH$coords

gdf_CH <- geomorph.data.frame(coords = CH_Shapes, Csize = as.numeric(Y.gpa_CH$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes) ... finding out if size matters
CHAnova <- procD.lm(coords~Csize, data = gdf_CH) 
summary(CHAnova) #it does?

#plot allometry of the origical not size corrected model
plotAllometry(CHAnova, size = gdf_CH$Csize, method = "RegScore")
#view the residuals of the model 
#CHAnova$residuals


## allometry is having an effect and how we are removing the effect 

# save size-adjusted residuals an make as 3D array
CHResiduals <- arrayspecs(CHAnova$residuals, p=dim(CH_Shapes)[1], k=dim(CH_Shapes)[2])
# make allometry-free shapes
CHResiduals.adj.shape <- CHResiduals+array(Y.gpa_CH$consensus, dim(CHResiduals))

#Run PCA on new  size corrected shape data
CHResiduals_PCA.Allom <- gm.prcomp(CHResiduals.adj.shape)

#gettings pc scores
CHResiduals_PCA.Allom$x
write.csv(CHResiduals_PCA.Allom$x[,c(1,2,3)], "CHResiduals_PC1_3.csv")

#ploting the size corrected specimen
plotAllSpecimens(CHResiduals.adj.shape)

#plotting the PCA pc 1+2 
plot(CHResiduals_PCA.Allom)
#text(psPCA_CL$x, labels = dimnames(Y.gpa_CL$coords)[[3]])
text(CHResiduals_PCA.Allom$x, pos = 4, label = dimnames(CHResiduals.adj.shape)[[3]], cex=.5)

#plotting the PCA pc2 + 3
plot(CHResiduals_PCA.Allom, axis1 = 2, axis2 = 3)
text(CHResiduals_PCA.Allom$x[, 2], CHResiduals_PCA.Allom$x[, 3], pos = 4, label = dimnames(CHResiduals.adj.shape)[[3]], cex=.5)

#new lm SIZE CORRECTED
gdf_CHR <- geomorph.data.frame(coords = CHResiduals.adj.shape, Csize = as.numeric(Y.gpa_CH$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes)##
CHAnovaR <- procD.lm(coords~Csize, data = gdf_CHR) 
summary(CHAnovaR) 


##ball and stick plots

open3d()
#PC1
plotRefToTarget(M1 = CHResiduals.adj.shape  [,,"LFxTRC5448"], M2 = CHResiduals.adj.shape [,,"LFxTRC5084"],
                method = "vector", mag = 2)




#CB1
#Allometry set up 
CB1_Shapes<-Y.gpa_CB1$coords

gdf_CB1 <- geomorph.data.frame(coords = CB1_Shapes, Csize = as.numeric(Y.gpa_CB1$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes) ... finding out if size matters
CB1Anova <- procD.lm(coords~Csize, data = gdf_CB1) 
summary(CB1Anova) #it does?

#plot allometry of the origical not size corrected model
plotAllometry(CB1Anova, size = gdf_CB1$Csize, method = "RegScore")
#view the residuals of the model 
#CB1Anova$residuals


## allometry is having an effect and how we are removing the effect 

# save size-adjusted residuals an make as 3D array
CB1Residuals <- arrayspecs(CB1Anova$residuals, p=dim(CB1_Shapes)[1], k=dim(CB1_Shapes)[2])
# make allometry-free shapes
CB1Residuals.adj.shape <- CB1Residuals+array(Y.gpa_CB1$consensus, dim(CB1Residuals))

#Run PCA on new  size corrected shape data
CB1Residuals_PCA.Allom <- gm.prcomp(CB1Residuals.adj.shape)

#gettings pc scores
CB1Residuals_PCA.Allom$x
write.csv(CB1Residuals_PCA.Allom$x[,c(1,2,3)], "CB1Residuals_PC1_3.csv")

#ploting the size corrected specimen
plotAllSpecimens(CB1Residuals.adj.shape)

#plotting the PCA pc1 + 2
plot(CB1Residuals_PCA.Allom)
#text(psPCA_CL$x, labels = dimnames(Y.gpa_CL$coords)[[3]])
text(CB1Residuals_PCA.Allom$x, pos = 4, label = dimnames(CB1Residuals.adj.shape)[[3]], cex=.5)

#plotting the PCA pc2 + 3
plot(CB1Residuals_PCA.Allom, axis1 = 2, axis2 = 3)
text(CB1Residuals_PCA.Allom$x[, 2], CB1Residuals_PCA.Allom$x[, 3], pos = 4, label = dimnames(CB1Residuals.adj.shape)[[3]], cex=.5)


#new lm SIZE CORRECTED
gdf_CB1R <- geomorph.data.frame(coords = CB1Residuals.adj.shape, Csize = as.numeric(Y.gpa_CB1$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes)##
CB1AnovaR <- procD.lm(coords~Csize, data = gdf_CB1R) 
summary(CB1AnovaR) 

##ball and stick plots

open3d()
#PC1
plotRefToTarget(M1 = CB1Residuals.adj.shape  [,,"LFxTRC5208"], M2 = CB1Residuals.adj.shape [,,"LFxTRC5123"],
                method = "vector", mag = 2)


#CB5
#Allometry set up 
PJ_Shapes<-Y.gpa_PJ$coords

gdf_PJ <- geomorph.data.frame(coords = PJ_Shapes, Csize = as.numeric(Y.gpa_PJ$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes) ... finding out if size matters
PJAnova <- procD.lm(coords~Csize, data = gdf_PJ) 
summary(PJAnova) #it does?


#plot allometry of the origical not size corrected model
plotAllometry(PJAnova, size = gdf_PJ$Csize, method = "RegScore")
#view the residuals of the model 
#CHAnova$residuals


## allometry is having an effect and how we are removing the effect 

# save size-adjusted residuals an make as 3D array
PJResiduals <- arrayspecs(PJAnova$residuals, p=dim(PJ_Shapes)[1], k=dim(PJ_Shapes)[2])
# make allometry-free shapes
PJResiduals.adj.shape <- PJResiduals+array(Y.gpa_PJ$consensus, dim(PJResiduals))

#Run PCA on new  size corrected shape data
PJResiduals_PCA.Allom <- gm.prcomp(PJResiduals.adj.shape)

#gettings pc scores
PJResiduals_PCA.Allom$x
write.csv(PJResiduals_PCA.Allom$x[,c(1,2,3)], "PJResiduals_PC1_3.csv")

#ploting the size corrected specimen
plotAllSpecimens(PJResiduals.adj.shape)

#plotting the PCA pc 1+ 2
plot(PJResiduals_PCA.Allom)
#text(psPCA_CL$x, labels = dimnames(Y.gpa_CL$coords)[[3]])
text(PJResiduals_PCA.Allom$x, pos = 4, label = dimnames(PJResiduals.adj.shape)[[3]], cex=.5)

#plotting the PCA pc2 + 3
plot(PJResiduals_PCA.Allom, axis1 = 2, axis2 = 3)
text(PJResiduals_PCA.Allom$x[, 2], PJResiduals_PCA.Allom$x[, 3], pos = 4, label = dimnames(PJResiduals.adj.shape)[[3]], cex=.5)

#new lm SIZE CORRECTED
gdf_PJR <- geomorph.data.frame(coords = PJResiduals.adj.shape, Csize = as.numeric(Y.gpa_PJ$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes)##
PJAnovaR <- procD.lm(coords~Csize, data = gdf_PJR) 
summary(PJAnovaR) 


##ball and stick plots

open3d()
#PC1
plotRefToTarget(M1 = ManResiduals.adj.shape  [,,"LFxTRC5448"], M2 = ManResiduals.adj.shape [,,"LFxTRC5084"],
                method = "vector", mag = 1)




####extracting lengths
# Distance matrix for lmk distance 
lmks = data.frame(manD = c(6, 1), hyoidL = c(4, 7), CB1L = c(1, 5), CB5L = c(9, 1), CleiL = c(11, 2))
ManD = interlmkdist(gdf_ManR$coords, c(1, 6))
CH_L = interlmkdist(gdf_CH$coords, c(4, 7))
CB1L = interlmkdist(gdf_CB1$coords, c(1, 5))
CB5L = interlmkdist(gdf_PJR$coords, c(1, 9))
CleiL = interlmkdist(gdf_CL$coords, c(2, 11))

write.csv(ManD, file = "Mand_length.csv", row.names = TRUE)
write.csv(CH_L, file = "HYD_length.csv", row.names = TRUE)
write.csv(CB1L, file = "CB1_length.csv", row.names = TRUE)
write.csv(CB5L, file = "CB5_length.csv", row.names = TRUE)
write.csv(CleiL, file = "Clth_length.csv", row.names = TRUE)



#### MAN ON THE X- AXIS CORRS 

## man x ch 
#aligning the common individuals within the data
CommonInds2<-c(dimnames(ManResiduals.adj.shape)[[3]], dimnames(CHResiduals.adj.shape)[[3]])[duplicated(c(dimnames(ManResiduals.adj.shape)[[3]], dimnames(CHResiduals.adj.shape)[[3]]))]

#CH setup
CHnew<-CHResiduals.adj.shape[,,CommonInds2]

#Man setup
ManforCHnew<-ManResiduals.adj.shape[,,CommonInds2]

#running correlations 
MANCHCorrs <- two.b.pls(A1 = ManforCHnew, A2 = CHnew, iter = 9999)
plot(MANCHCorrs, label = CommonInds2)
plot(MANCHCorrs, pch=19) #xlim=c(-0.2,0.2), ylim=c(-0.2,0.2)
#plot(MANCHCorrs, pch=19, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1), col="black")
summary(MANCHCorrs)

MANCHCorrs <- two.b.pls(A1 = ManforCHnew, A2 = CHnew, iter = 9999)
MANCHCorrs.plot <- plot(MANCHCorrs, pch=19, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
#picknplot.shape(MANCHCorrs.plot,  method="vector", mg = 4)
MANCHCorrs.plot_shape = picknplot.shape(MANCHCorrs.plot)
#block 1 man: 
plotRefToTarget(MANCHCorrs.plot_shape$shapes[[1]]$P1, MANCHCorrs.plot_shape$shapes[[2]]$P1, method = "vector", mag=2)

#block 2 ch
plotRefToTarget(MANCHCorrs.plot_shape$shapes[[1]]$P2, MANCHCorrs.plot_shape$shapes[[2]]$P2, method = "vector", mag=2)
#this will insert text to plot 
#text(MANCHCorrs$YScores, labels = dimnames(Y.gpa_Man$coords)[[3]])



## man x cb1
#aligning the common individuals within the data
CommonInds3<-c(dimnames(ManResiduals.adj.shape)[[3]], dimnames(CB1Residuals.adj.shape)[[3]])[duplicated(c(dimnames(ManResiduals.adj.shape)[[3]], dimnames(CB1Residuals.adj.shape)[[3]]))]

#Cb1 setup
CB1new<-CB1Residuals.adj.shape[,,CommonInds3]

#Man setup
ManforCB1new<-ManResiduals.adj.shape[,,CommonInds3]

#running correlations 
MANCB1Corrs <- two.b.pls(A1 = ManforCB1new, A2 = CB1new, iter = 9999)
plot(MANCB1Corrs, label = CommonInds3)
plot(MANCB1Corrs, pch=19)
#plot(MANCB1Corrs, pch=19, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
summary(MANCB1Corrs)

#this will insert text to plot 
#text(MANCB1Corrs$YScores, labels = dimnames(Y.gpa_Man$coords)[[3]])

MANCB1Corrs <- two.b.pls(A1 = ManforCB1new, A2 = CB1new, iter = 9999)
MANCB1Corrs.plot <- plot(MANCB1Corrs, pch=19, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
#picknplot.shape(MANCB1Corrs.plot,  method="vector", mg = 4)
MANCB1Corrs.plot_shape = picknplot.shape(MANCB1Corrs.plot)


#block 1 man: 
plotRefToTarget(MANCB1Corrs.plot_shape$shapes[[1]]$P1, MANCB1Corrs.plot_shape$shapes[[2]]$P1, method = "vector", mag=2)
#block 2 cb1
plotRefToTarget(MANCB1Corrs.plot_shape$shapes[[1]]$P2, MANCB1Corrs.plot_shape$shapes[[2]]$P2, method = "vector", mag=2)




## man x cb5 
#aligning the common individuals within the data
CommonInds4<-c(dimnames(ManResiduals.adj.shape)[[3]], dimnames(PJResiduals.adj.shape)[[3]])[duplicated(c(dimnames(ManResiduals.adj.shape)[[3]], dimnames(PJResiduals.adj.shape)[[3]]))]

#PJ setup
PJnew<-PJResiduals.adj.shape[,,CommonInds4]

#Man setup
MannewforPJ<-ManResiduals.adj.shape[,,CommonInds4]


#running correlations 
MANPJCorrs <- two.b.pls(A1 = MannewforPJ, A2 = PJnew, iter = 9999)
plot(MANPJCorrs, label = CommonInds4)
#plot(MANPJCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
plot(MANPJCorrs, pch=19)
summary(MANPJCorrs)

#this will insert text to plot 
#text(MANPJCorrs$YScores, labels = dimnames(Y.gpa_Man$coords)[[3]])

MANPJCorrs <- two.b.pls(A1 = MannewforPJ, A2 = PJnew, iter = 9999)
MANPJCorrs.plot <- plot(MANPJCorrs, pch=19, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
#picknplot.shape(MANPJCorrs.plot,  method="vector", mg = 4)
MANPJCorrs.plot_shape = picknplot.shape(MANPJCorrs.plot)


#block 1 man: 
plotRefToTarget(MANPJCorrs.plot_shape$shapes[[1]]$P1, MANPJCorrs.plot_shape$shapes[[2]]$P1, method = "vector", mag=2)
#block 2 cb5
plotRefToTarget(MANPJCorrs.plot_shape$shapes[[1]]$P2, MANPJCorrs.plot_shape$shapes[[2]]$P2, method = "vector", mag=2)




## man x cl 
#aligning the common individuals within the data
CommonInds<-c(dimnames(ManResiduals.adj.shape)[[3]], dimnames(CLResiduals.adj.shape)[[3]])[duplicated(c(dimnames(ManResiduals.adj.shape)[[3]], dimnames(CLResiduals.adj.shape)[[3]]))]


#Cl setup
CLnew<-CLResiduals.adj.shape[,,CommonInds]

#Man setup
Mannew<-ManResiduals.adj.shape[,,CommonInds]


#running correlations 
MANCLCorrs <- two.b.pls(A1 = Mannew, A2 = CLnew, iter = 9999)
plot(MANCLCorrs, label = CommonInds)
plot(MANCLCorrs, pch=19)
#plot(MANCLCorrs, pch=19, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
summary(MANCLCorrs)


MANCLCorrs <- two.b.pls(A1 = Mannew, A2 = CLnew, iter = 9999)
MANCLCorrs.plot <- plot(MANCLCorrs, pch=19, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
#picknplot.shape(MANCLCorrs.plot,  method="vector", mg = 4)
MANCLCorrs.plot_shape = picknplot.shape(MANCLCorrs.plot)


#block 1 man: 
plotRefToTarget(MANCLCorrs.plot_shape$shapes[[1]]$P1, MANCLCorrs.plot_shape$shapes[[2]]$P1, method = "vector", mag=2)
#block 2 cl
plotRefToTarget(MANCLCorrs.plot_shape$shapes[[1]]$P2, MANCLCorrs.plot_shape$shapes[[2]]$P2, method = "vector", mag=2)

#this will insert text to plot 
#text(MANCLCorrs$YScores, labels = dimnames(Y.gpa_Man$coords)[[3]])




#### CL ON THE X- AXIS CORRS 

## cl x cb5 
#aligning the common individuals within the data
CommonInds5<-c(dimnames(CLResiduals.adj.shape)[[3]], dimnames(PJResiduals.adj.shape)[[3]])[duplicated(c(dimnames(CLResiduals.adj.shape)[[3]], dimnames(PJResiduals.adj.shape)[[3]]))]

#PJ setup
PJnew_CL<-PJResiduals.adj.shape[,,CommonInds5]

#CL setup
CLnewforPJ<-CLResiduals.adj.shape[,,CommonInds5]


#running correlations 
CLPJCorrs <- two.b.pls(A1 = CLnewforPJ, A2 = PJnew_CL, iter = 9999)
plot(CLPJCorrs, label = CommonInds5)
#plot(CLPJCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
plot(CLPJCorrs, pch=19)
summary(CLPJCorrs)

#this will insert text to plot 
#text(CLPJCorrs$YScores, labels = dimnames(Y.gpa_CL$coords)[[3]])

CLPJCorrs <- two.b.pls(A1 = CLnewforPJ, A2 = PJnew_CL, iter = 9999)
CLPJCorrs.plot <- plot(CLPJCorrs, pch=19, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
#picknplot.shape(CLPJCorrs.plot,  method="vector", mg = 4)
CLPJCorrs.plot_shape = picknplot.shape(CLPJCorrs.plot)


#block 1 cl: 
plotRefToTarget(CLPJCorrs.plot_shape$shapes[[1]]$P1, CLPJCorrs.plot_shape$shapes[[2]]$P1, method = "vector", mag=2)
#block 2 cb5
plotRefToTarget(CLPJCorrs.plot_shape$shapes[[1]]$P2, CLPJCorrs.plot_shape$shapes[[2]]$P2, method = "vector", mag=2)





## cl x cb1 
#aligning the common individuals within the data
CommonInds6<-c(dimnames(CLResiduals.adj.shape)[[3]], dimnames(CB1Residuals.adj.shape)[[3]])[duplicated(c(dimnames(CLResiduals.adj.shape)[[3]], dimnames(CB1Residuals.adj.shape)[[3]]))]

#PJ setup
CB1new_CL<-CB1Residuals.adj.shape[,,CommonInds6]

#CL setup
CLnewforCB1<-CLResiduals.adj.shape[,,CommonInds6]


#running correlations 
CLCB1Corrs <- two.b.pls(A1 = CLnewforCB1, A2 = CB1new_CL, iter = 9999)
plot(CLCB1Corrs, label = CommonInds6)
plot(CLCB1Corrs, pch=19)
#plot(CLCB1Corrs, pch=19, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
summary(CLCB1Corrs)

#this will insert text to plot 
#text(CLCB1Corrs$YScores, labels = dimnames(Y.gpa_CL$coords)[[3]]) 

CLCB1Corrs <- two.b.pls(A1 = CLnewforCB1, A2 = CB1new_CL, iter = 9999)
CLCB1Corrs.plot <- plot(CLCB1Corrs, pch=19, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
#picknplot.shape(CLCB1Corrs.plot,  method="vector", mg = 4)
CLCB1Corrs.plot_shape = picknplot.shape(CLCB1Corrs.plot)


#block 1 cl: 
plotRefToTarget(CLCB1Corrs.plot_shape$shapes[[1]]$P1, CLCB1Corrs.plot_shape$shapes[[2]]$P1, method = "vector", mag=2)
#block 2 cb1
plotRefToTarget(CLCB1Corrs.plot_shape$shapes[[1]]$P2, CLCB1Corrs.plot_shape$shapes[[2]]$P2, method = "vector", mag=2)






## cl x ch 
#aligning the common individuals within the data
CommonInds7<-c(dimnames(CLResiduals.adj.shape)[[3]], dimnames(CHResiduals.adj.shape)[[3]])[duplicated(c(dimnames(CLResiduals.adj.shape)[[3]], dimnames(CHResiduals.adj.shape)[[3]]))]

#PJ setup
CHnew_CL<-CHResiduals.adj.shape[,,CommonInds7]

#CL setup
CLnewforCH<-CLResiduals.adj.shape[,,CommonInds7]


#running correlations 
CLCHCorrs <- two.b.pls(A1 = CLnewforCH, A2 = CHnew_CL, iter = 9999)
plot(CLCHCorrs, label = CommonInds7)
plot(CLCHCorrs, pch=19)
#plot(CLCHCorrs, pch=19, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
summary(CLCHCorrs)

CLCHCorrs <- two.b.pls(A1 = CLnewforCH, A2 = CHnew_CL, iter = 9999)
CLCHCorrs.plot <- plot(CLCHCorrs, pch=19, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
#picknplot.shape(CLCHCorrs.plot,  method="vector", mg = 4)
CLCHCorrs.plot_shape = picknplot.shape(CLCHCorrs.plot)


#block 1 cl: 
plotRefToTarget(CLCHCorrs.plot_shape$shapes[[1]]$P1, CLCHCorrs.plot_shape$shapes[[2]]$P1, method = "vector", mag=2)
#block 2 ch
plotRefToTarget(CLCHCorrs.plot_shape$shapes[[1]]$P2, CLCHCorrs.plot_shape$shapes[[2]]$P2, method = "vector", mag=2)

#this will insert text to plot 
#text(CLCHCorrs$YScores, labels = dimnames(Y.gpa_CL$coords)[[3]])


## cl x man 
#aligning the common individuals within the data
CommonInds8<-c(dimnames(CLResiduals.adj.shape)[[3]], dimnames(ManResiduals.adj.shape)[[3]])[duplicated(c(dimnames(CLResiduals.adj.shape)[[3]], dimnames(ManResiduals.adj.shape)[[3]]))]

#PJ setup
Mannew_CL<-ManResiduals.adj.shape[,,CommonInds8]

#CL setup
CLnewforMan<-CLResiduals.adj.shape[,,CommonInds8]


#running correlations 
CLManCorrs <- two.b.pls(A1 = CLnewforMan, A2 = Mannew_CL, iter = 9999)
plot(CLManCorrs, label = CommonInds8)
plot(CLManCorrs, pch=19)
#plot(CLManCorrs, pch=19, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
summary(CLManCorrs)


CLManCorrs <- two.b.pls(A1 = CLnewforMan, A2 = Mannew_CL, iter = 9999)
CLManCorrs.plot <- plot(CLManCorrs, pch=19, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
#picknplot.shape(CLManCorrs.plot,  method="vector", mg = 4)
CLManCorrs.plot_shape = picknplot.shape(CLManCorrs.plot)


#block 1 cl: 
plotRefToTarget(CLManCorrs.plot_shape$shapes[[1]]$P1, CLManCorrs.plot_shape$shapes[[2]]$P1, method = "vector", mag=2)
#block 2 man
plotRefToTarget(CLManCorrs.plot_shape$shapes[[1]]$P2, CLManCorrs.plot_shape$shapes[[2]]$P2, method = "vector", mag=2)


#this will insert text to plot 
#text(CLManCorrs$YScores, labels = dimnames(Y.gpa_CL$coords)[[3]])


### HYOID on x axis 
####Correlation

# HYD vs Mand 
#aligning the common individuals within the data
CommonInds9<-c(dimnames(CHResiduals.adj.shape)[[3]], dimnames(ManResiduals.adj.shape)[[3]])[duplicated(c(dimnames(CHResiduals.adj.shape)[[3]], dimnames(ManResiduals.adj.shape)[[3]]))]

#MAN setup
Mannew_CH<-ManResiduals.adj.shape[,,CommonInds9]

#CH setup
CHnewforMan<-CHResiduals.adj.shape[,,CommonInds9]


#running correlations 
CHManCorrs <- two.b.pls(A1 = CHnewforMan, A2 = Mannew_CH, iter = 9999)
plot(CHManCorrs, label = CommonInds9)
plot(CHManCorrs, pch=19, xlim=c(-0.2,0.2), ylim=c(-0.2,0.2))
summary(CHManCorrs)


# HYD vs CB1 
#aligning the common individuals within the data
CommonInds10<-c(dimnames(CHResiduals.adj.shape)[[3]], dimnames(CB1Residuals.adj.shape)[[3]])[duplicated(c(dimnames(CHResiduals.adj.shape)[[3]], dimnames(CB1Residuals.adj.shape)[[3]]))]

#CB1 setup
CB1new_CH<-CB1Residuals.adj.shape[,,CommonInds10]

#CH setup
CHnewforCB1<-CHResiduals.adj.shape[,,CommonInds10]


#running correlations 
CHCB1Corrs <- two.b.pls(A1 = CHnewforCB1, A2 = CB1new_CH, iter = 9999)
plot(CHCB1Corrs, label = CommonInds10)
plot(CHCB1Corrs, pch=19, xlim=c(-0.2,0.2), ylim=c(-0.2,0.2))
summary(CHCB1Corrs)


# HYD vs CB5
#aligning the common individuals within the data
CommonInds11<-c(dimnames(CHResiduals.adj.shape)[[3]], dimnames(PJResiduals.adj.shape)[[3]])[duplicated(c(dimnames(CHResiduals.adj.shape)[[3]], dimnames(PJResiduals.adj.shape)[[3]]))]

#CB5 setup
CB5new_CH<-PJResiduals.adj.shape[,,CommonInds11]

#CH setup
CHnewforCB5<-CHResiduals.adj.shape[,,CommonInds11]


#running correlations 
CHCB5Corrs <- two.b.pls(A1 = CHnewforCB5, A2 = CB5new_CH, iter = 9999)
plot(CHCB5Corrs, label = CommonInds11)
plot(CHCB5Corrs, pch=19, xlim=c(-0.2,0.2), ylim=c(-0.2,0.2))
summary(CHCB5Corrs)

# HYD vs CL
#aligning the common individuals within the data
CommonInds12<-c(dimnames(CHResiduals.adj.shape)[[3]], dimnames(CLResiduals.adj.shape)[[3]])[duplicated(c(dimnames(CHResiduals.adj.shape)[[3]], dimnames(CLResiduals.adj.shape)[[3]]))]

#CL setup
CLnew_CH<-CLResiduals.adj.shape[,,CommonInds12]

#CH setup
CHnewforCL<-CHResiduals.adj.shape[,,CommonInds12]


#running correlations 
CHCLCorrs <- two.b.pls(A1 = CHnewforCL, A2 = CLnew_CH, iter = 9999)
plot(CHCLCorrs, label = CommonInds12)
plot(CHCLCorrs, pch=19, xlim=c(-0.2,0.2), ylim=c(-0.2,0.2))
summary(CHCLCorrs)


### CB1 on x axis 
####Correlation

# CB1 vs Hyd 
#aligning the common individuals within the data
CommonInds13<-c(dimnames(CB1Residuals.adj.shape)[[3]], dimnames(CHResiduals.adj.shape)[[3]])[duplicated(c(dimnames(CB1Residuals.adj.shape)[[3]], dimnames(CHResiduals.adj.shape)[[3]]))]

#CH setup
CHnew_CB1<- CHResiduals.adj.shape[,,CommonInds13]

#CB1 setup
CB1newforCH<-CB1Residuals.adj.shape[,,CommonInds13]


#running correlations 
CB1CHCorrs <- two.b.pls(A1 = CB1newforCH, A2 = CHnew_CB1, iter = 9999)
plot(CB1CHCorrs, label = CommonInds13)
plot(CB1CHCorrs, pch=19, xlim=c(-0.2,0.2), ylim=c(-0.2,0.2))
summary(CB1CHCorrs)

# CB1 vs Mand
#aligning the common individuals within the data
CommonInds14<-c(dimnames(CB1Residuals.adj.shape)[[3]], dimnames(ManResiduals.adj.shape)[[3]])[duplicated(c(dimnames(CB1Residuals.adj.shape)[[3]], dimnames(ManResiduals.adj.shape)[[3]]))]

#Man setup
Mannew_CB1<- ManResiduals.adj.shape[,,CommonInds14]

#CB1 setup
CB1newforMan<-CB1Residuals.adj.shape[,,CommonInds14]


#running correlations 
CB1ManCorrs <- two.b.pls(A1 = CB1newforMan, A2 = Mannew_CB1, iter = 9999)
plot(CB1ManCorrs, label = CommonInds14)
plot(CB1ManCorrs, pch=19, xlim=c(-0.2,0.2), ylim=c(-0.2,0.2))
summary(CB1ManCorrs)

# CB1 vs CB5 
#aligning the common individuals within the data
CommonInds15<-c(dimnames(CB1Residuals.adj.shape)[[3]], dimnames(PJResiduals.adj.shape)[[3]])[duplicated(c(dimnames(CB1Residuals.adj.shape)[[3]], dimnames(PJResiduals.adj.shape)[[3]]))]

#CB5 setup
CB5new_CB1<- PJResiduals.adj.shape[,,CommonInds15]

#CB1 setup
CB1newforCB5<-CB1Residuals.adj.shape[,,CommonInds15]


#running correlations 
CB1CB5Corrs <- two.b.pls(A1 = CB1newforCB5, A2 = CB5new_CB1, iter = 9999)
plot(CB1CB5Corrs, label = CommonInds15)
plot(CB1CB5Corrs, pch=19, xlim=c(-0.2,0.2), ylim=c(-0.2,0.2))
summary(CB1CB5Corrs)

# CB1 vs CL 
#aligning the common individuals within the data
CommonInds16<-c(dimnames(CB1Residuals.adj.shape)[[3]], dimnames(CLResiduals.adj.shape)[[3]])[duplicated(c(dimnames(CB1Residuals.adj.shape)[[3]], dimnames(CLResiduals.adj.shape)[[3]]))]

#CL setup
CLnew_CB1<- CLResiduals.adj.shape[,,CommonInds16]

#CB1 setup
CB1newforCL<-CB1Residuals.adj.shape[,,CommonInds16]


#running correlations 
CB1CLCorrs <- two.b.pls(A1 = CB1newforCL, A2 = CLnew_CB1, iter = 9999)
plot(CB1CLCorrs, label = CommonInds16)
plot(CB1CLCorrs, pch=19, xlim=c(-0.2,0.2), ylim=c(-0.2,0.2))
summary(CB1CLCorrs)


### CB5 on x axis 
####Correlation

# CB5 vs Hyd 
#aligning the common individuals within the data
CommonInds17<-c(dimnames(PJResiduals.adj.shape)[[3]], dimnames(CHResiduals.adj.shape)[[3]])[duplicated(c(dimnames(PJResiduals.adj.shape)[[3]], dimnames(CHResiduals.adj.shape)[[3]]))]

#CH setup
CHnew_CB5<- CHResiduals.adj.shape[,,CommonInds17]

#CB5 setup
CB5newforCH<-PJResiduals.adj.shape[,,CommonInds17]


#running correlations 
CB5CHCorrs <- two.b.pls(A1 = CB5newforCH, A2 = CHnew_CB5, iter = 9999)
plot(CB5CHCorrs, label = CommonInds17)
plot(CB5CHCorrs, pch=19, xlim=c(-0.2,0.2), ylim=c(-0.2,0.2))
summary(CB5CHCorrs)

# CB5 vs Mand
#aligning the common individuals within the data
CommonInds18<-c(dimnames(PJResiduals.adj.shape)[[3]], dimnames(ManResiduals.adj.shape)[[3]])[duplicated(c(dimnames(PJResiduals.adj.shape)[[3]], dimnames(ManResiduals.adj.shape)[[3]]))]

#Man setup
Mannew_CB5<- ManResiduals.adj.shape[,,CommonInds18]

#CB5 setup
CB5newforMan<-PJResiduals.adj.shape[,,CommonInds18]


#running correlations 
CB5ManCorrs <- two.b.pls(A1 = CB5newforMan, A2 = Mannew_CB5, iter = 9999)
plot(CB5ManCorrs, label = CommonInds18)
plot(CB5ManCorrs, pch=19, xlim=c(-0.2,0.2), ylim=c(-0.2,0.2))
summary(CB5ManCorrs)

# CB5 vs CB1 
#aligning the common individuals within the data
CommonInds19<-c(dimnames(PJResiduals.adj.shape)[[3]], dimnames(CB1Residuals.adj.shape)[[3]])[duplicated(c(dimnames(PJResiduals.adj.shape)[[3]], dimnames(CB1Residuals.adj.shape)[[3]]))]

#CB1 setup
CB1new_CB5<- CB1Residuals.adj.shape[,,CommonInds19]

#CB5 setup
CB5newforCB1<-PJResiduals.adj.shape[,,CommonInds19]


#running correlations 
CB5CB1Corrs <- two.b.pls(A1 = CB5newforCB1, A2 = CB1new_CB5, iter = 9999)
plot(CB5CB1Corrs, label = CommonInds19)
plot(CB5CB1Corrs, pch=19, xlim=c(-0.2,0.2), ylim=c(-0.2,0.2))
summary(CB5CB1Corrs)

# CB5 vs Clth
#aligning the common individuals within the data
CommonInds20<-c(dimnames(PJResiduals.adj.shape)[[3]], dimnames(CLResiduals.adj.shape)[[3]])[duplicated(c(dimnames(PJResiduals.adj.shape)[[3]], dimnames(CLResiduals.adj.shape)[[3]]))]

#CL setup
CLnew_CB5<- CLResiduals.adj.shape[,,CommonInds20]

#CB5 setup
CB5newforCL<-PJResiduals.adj.shape[,,CommonInds20]


#running correlations 
CB5CLCorrs <- two.b.pls(A1 = CB5newforCL, A2 = CLnew_CB5, iter = 9999)
plot(CB5CLCorrs, label = CommonInds20)
plot(CB5CLCorrs, pch=19, xlim=c(-0.2,0.2), ylim=c(-0.2,0.2))
summary(CB5CLCorrs)


--------
###heat map###
data <- data.frame(
  dims = c("Mandible", "Cleithrum"),
  Mandible = c(0, 4.72116),
  Hyoid = c(4.83999, 7.49431),
  CB1 = c(3.87891, 3.98503),
  CB5 = c(7.0182, 7.05901),
  Cleithrum = c(4.80011, 0)
)



library(reshape2)

# Reshape the data for heatmap
data_melted <- melt(data, id.vars = "dims", variable.name = "variables", value.name = "value")


library(ggplot2)

# Create the heatmap
ggplot(data_melted, aes(x = dims, y = variables, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "cornflowerblue", mid = "khaki1", high = "palevioletred2", midpoint = 4.0,
                       limits = c(0.0, 8.0), oob = scales::squish) +
  labs(title = "Heatmap of Data", fill = "Z Score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45,         # Angle of text (rotate text)
                               hjust = 1,          # Horizontal justification
                               size = 12,          # Font size of the x-axis labels
                               color = "black",     # Color of the x-axis labels
                               family = "Arial"),  # Keep the angle for x-axis text
    axis.text.y = element_text(angle = 45,         # Angle of text (rotate text)
                               hjust = 1,          # Horizontal justification
                               size = 12,          # Font size of the x-axis labels
                               color = "black",     # Color of the x-axis labels
                               family = "Arial"),  # Keep the angle for y-axis text
    axis.title.x = element_blank(),  # Remove x-axis label
    axis.title.y = element_blank(),  # Remove y-axis label
    axis.ticks = element_blank()     # Remove axis ticks
  ) 

