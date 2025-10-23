library(StereoMorph) 
library(geomorph)
library(geiger)
library(plotrix)
library(caper)
library(ggplot2)
library(tidyverse)
library(tools)


#loading in the CL data 
CL <- readland.tps(file = "AllCLshape.tps", specID = "imageID")
CL_sub <- CL[-c(13),,]

writeland.tps(CL_sub, file = "CL_land_new.tps", specID = TRUE)

Y.gpa_CL<-gpagen(CL_sub)
plotAllSpecimens(Y.gpa_CL$coords)

MYDATA<-gm.prcomp(Y.gpa_CL$coords)
plot(MYDATA)
#text(MYDATA$x, labels = dimnames(Y.gpa_CL$coords)[[3]])


#loading in the CH data 
CH <- readland.tps(file = "AllCHshape.tps", specID = "imageID")
#deleting overlapping lms and semis
CH_sub <- CH[-c(8,17,18,27),,]

writeland.tps(CH_sub, file = "HYD_land_new.tps", specID = TRUE)

#visualising the lms
#define.sliders(CH_sub[,,1], write.file = "CH_sliders.csv")

#loading in curves
curves_CH<-as.matrix(read.csv("CH_sliders2.csv"))
Y.gpa_CH<-gpagen(CH_sub, curves = curves_CH)
#plots all individuals across all lms
plotAllSpecimens(Y.gpa_CH$coords)



MYDATA2<-gm.prcomp(Y.gpa_CH$coords)
#pc1 and 2
plot(MYDATA2)
#text for plot
#text(MYDATA3$x, labels = dimnames(Y.gpa_CH$coords)[[3]])

#pc2 and 3
#plot(MYDATA2, axis1 = 2, axis2 = 3)


#loading in CB1 data
CB1 <- readland.tps(file = "AllCB1shape.tps", specID = "imageID")
#deleting overlapping lms and semis
CB1_sub <- CB1[-c(6,15),,]

writeland.tps(CB1_sub, file = "CB1_land_new.tps", specID = TRUE)

#visualising the lms
#define.sliders(CB1_sub[,,1], write.file = "CB1_sliders.csv")


#loading in curves
curves_CB1<-as.matrix(read.csv("CB1_sliders2.csv"))
Y.gpa_CB1<-gpagen(CB1_sub, curves = curves_CB1)

write.csv(Y.gpa_CB1$Csize, "CB1Csize.csv")

#plots all individuals across all lms
plotAllSpecimens(Y.gpa_CB1$coords)


MYDATA3<-gm.prcomp(Y.gpa_CB1$coords)
plot(MYDATA3)

#text for plot
#text(MYDATA4$x, labels = dimnames(Y.gpa_CB1$coords)[[3]])

#loading in man data 
Man <- readland.tps(file = "AllRManshape.tps", specID = "ID" )

#pca graphs
Y.gpa_Man<-gpagen(Man)

MYDATA2<-gm.prcomp(Y.gpa_Man$coords)
plot(MYDATA2)

#text for plot
#text(MYDATA2$x, labels = dimnames(Y.gpa_Man$coords)[[3]])

#loading in pj data 
PJ <- readland.tps(file = "AllRCB5shape.tps", specID = "ID")

#pca graphs ... reduces variation to fewer dementions 
Y.gpa_PJ<-gpagen(PJ)

MYDATA4<-gm.prcomp(Y.gpa_PJ$coords)
plot(MYDATA4)

#text for plot
#text(MYDATA4$x, labels = dimnames(Y.gpa_PJ$coords)[[3]])



### ALLOMETRY AND TREE CORRECTIONS 
#Cl
CL_CichlidMeanShapes<-Y.gpa_CL$coords
CL_CichlidCSizes<-as.data.frame(Y.gpa_CL$Csize)

CL_CichlidCSizes1<-cbind.data.frame(CL_CichlidCSizes, rownames(CL_CichlidCSizes))
colnames(CL_CichlidCSizes1)<-c("Csize","Phylonames")

TaxaIDData<-read.csv("Microevo_uct1.csv")
Tree<-read.tree("CichlidTreeNature2020.tre")

TaxaIDData<-merge(TaxaIDData, CL_CichlidCSizes1, by = "Phylonames")
rownames(TaxaIDData)<-TaxaIDData$Phylonames

Pruning<-treedata(Tree, TaxaIDData)
Ftree<-Pruning$phy
FTree<-ladderize(Ftree)
FData<-Pruning$data

FData <- FData[match(FTree$tip.label, rownames(FData)),]
FData<-as.data.frame(FData)

RM_Cichlid<-dimnames(CL_CichlidMeanShapes)[[3]][dimnames(CL_CichlidMeanShapes)[[3]]%in%rownames(FData)]
CL_CichlidMeanShapes_Tr<-CL_CichlidMeanShapes[,,RM_Cichlid]


CL_CichlidMeanShapes_Tr<-CL_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(CL_CichlidMeanShapes_Tr)[[3]]))]

gdf <- geomorph.data.frame(coords = CL_CichlidMeanShapes_Tr, species = dimnames(CL_CichlidMeanShapes_Tr)[[3]], Csize = as.numeric(FData$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes)##
HybAnova <- procD.pgls(coords~Csize, phy = FTree, data = gdf) 
summary(HybAnova) 

shape.resid.CL <- arrayspecs(HybAnova$pgls.residuals, p=dim(CL_CichlidMeanShapes_Tr)[1], k=dim(CL_CichlidMeanShapes_Tr)[2]) 

# size-adjusted residuals
CL_CichlidMeanShapes_Tr1 <- shape.resid.CL + array(Y.gpa_CL$consensus, dim(shape.resid.CL))

#Run PCA on new  size corrected shape data
CL_CichlidMeanShapes_Tr1_PCA.Allom <- gm.prcomp(CL_CichlidMeanShapes_Tr1)

#ploting the size corrected specimen
plotAllSpecimens(CL_CichlidMeanShapes_Tr1)

#plotting the PCA pc1 + 2
plot(CL_CichlidMeanShapes_Tr1_PCA.Allom)
#text(psPCA_CL$x, labels = dimnames(Y.gpa_CL$coords)[[3]])
text(CL_CichlidMeanShapes_Tr1_PCA.Allom$x, pos = 4, label = dimnames(CL_CichlidMeanShapes_Tr1)[[3]], cex=.5)



#CH 
CH_CichlidMeanShapes<-Y.gpa_CH$coords
CH_CichlidCSizes<-as.data.frame(Y.gpa_CH$Csize)

CH_CichlidCSizes1<-cbind.data.frame(CH_CichlidCSizes, rownames(CH_CichlidCSizes))
colnames(CH_CichlidCSizes1)<-c("Csize","Phylonames")

TaxaIDData<-read.csv("Microevo_uct1.csv")
Tree<-read.tree("CichlidTreeNature2020.tre")

TaxaIDData<-merge(TaxaIDData, CH_CichlidCSizes1, by = "Phylonames")
rownames(TaxaIDData)<-TaxaIDData$Phylonames

Pruning<-treedata(Tree, TaxaIDData)
Ftree<-Pruning$phy
FTree<-ladderize(Ftree)
FData<-Pruning$data

FData <- FData[match(FTree$tip.label, rownames(FData)),]
FData<-as.data.frame(FData)

RM_Cichlid_CH<-dimnames(CH_CichlidMeanShapes)[[3]][dimnames(CH_CichlidMeanShapes)[[3]]%in%rownames(FData)]
CH_CichlidMeanShapes_Tr<-CH_CichlidMeanShapes[,,RM_Cichlid_CH]


CH_CichlidMeanShapes_Tr<-CH_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(CH_CichlidMeanShapes_Tr)[[3]]))]

gdf_CH <- geomorph.data.frame(coords = CH_CichlidMeanShapes_Tr, species = dimnames(CH_CichlidMeanShapes_Tr)[[3]], Csize = as.numeric(FData$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes)##
HybAnova_CH <- procD.pgls(coords~Csize, phy = FTree, data = gdf_CH) 
summary(HybAnova_CH) 

shape.resid.CH <- arrayspecs(HybAnova_CH$pgls.residuals, p=dim(CH_CichlidMeanShapes_Tr)[1], k=dim(CH_CichlidMeanShapes_Tr)[2]) # size-adjusted residuals
CH_CichlidMeanShapes_Tr1 <- shape.resid.CH + array(Y.gpa_CH$consensus, dim(shape.resid.CH))

#Run PCA on new  size corrected shape data
CH_CichlidMeanShapes_Tr1_PCA.Allom <- gm.prcomp(CH_CichlidMeanShapes_Tr1)

#ploting the size corrected specimen
plotAllSpecimens(CH_CichlidMeanShapes_Tr1)

#plotting the PCA pc1 + 2
plot(CH_CichlidMeanShapes_Tr1_PCA.Allom)
#text(psPCA_CL$x, labels = dimnames(Y.gpa_CL$coords)[[3]])
text(CH_CichlidMeanShapes_Tr1_PCA.Allom$x, pos = 4, label = dimnames(CH_CichlidMeanShapes_Tr1)[[3]], cex=.5)



#CB1 
CB1_CichlidMeanShapes<-Y.gpa_CB1$coords
CB1_CichlidCSizes<-as.data.frame(Y.gpa_CB1$Csize)

CB1_CichlidCSizes1<-cbind.data.frame(CB1_CichlidCSizes, rownames(CB1_CichlidCSizes))
colnames(CB1_CichlidCSizes1)<-c("Csize","Phylonames")

TaxaIDData<-read.csv("Microevo_uct1.csv")
Tree<-read.tree("CichlidTreeNature2020.tre")

TaxaIDData<-merge(TaxaIDData, CB1_CichlidCSizes1, by = "Phylonames")
rownames(TaxaIDData)<-TaxaIDData$Phylonames

Pruning<-treedata(Tree, TaxaIDData)
Ftree<-Pruning$phy
FTree<-ladderize(Ftree)
FData<-Pruning$data

FData <- FData[match(FTree$tip.label, rownames(FData)),]
FData<-as.data.frame(FData)

RM_Cichlid_CB1<-dimnames(CB1_CichlidMeanShapes)[[3]][dimnames(CB1_CichlidMeanShapes)[[3]]%in%rownames(FData)]
CB1_CichlidMeanShapes_Tr<-CB1_CichlidMeanShapes[,,RM_Cichlid_CB1]


CB1_CichlidMeanShapes_Tr<-CB1_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(CB1_CichlidMeanShapes_Tr)[[3]]))]

gdf_CB1 <- geomorph.data.frame(coords = CB1_CichlidMeanShapes_Tr, species = dimnames(CB1_CichlidMeanShapes_Tr)[[3]], Csize = as.numeric(FData$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes)##
HybAnova_CB1 <- procD.pgls(coords~Csize, phy = FTree, data = gdf_CB1) 
summary(HybAnova_CB1) 

shape.resid.CB1 <- arrayspecs(HybAnova_CB1$pgls.residuals, p=dim(CB1_CichlidMeanShapes_Tr)[1], k=dim(CB1_CichlidMeanShapes_Tr)[2]) # size-adjusted residuals
CB1_CichlidMeanShapes_Tr1 <- shape.resid.CB1 + array(Y.gpa_CB1$consensus, dim(shape.resid.CB1))

#Run PCA on new  size corrected shape data
CB1_CichlidMeanShapes_Tr1_PCA.Allom <- gm.prcomp(CB1_CichlidMeanShapes_Tr1)

#ploting the size corrected specimen
plotAllSpecimens(CB1_CichlidMeanShapes_Tr1)

#plotting the PCA pc1 + 2
plot(CB1_CichlidMeanShapes_Tr1_PCA.Allom)
#text(psPCA_CL$x, labels = dimnames(Y.gpa_CL$coords)[[3]])
text(CB1_CichlidMeanShapes_Tr1_PCA.Allom$x, pos = 4, label = dimnames(CB1_CichlidMeanShapes_Tr1)[[3]], cex=.5)



#MAN
Man_CichlidMeanShapes<-Y.gpa_Man$coords
Man_CichlidCSizes<-as.data.frame(Y.gpa_Man$Csize)

Man_CichlidCSizes1<-cbind.data.frame(Man_CichlidCSizes, rownames(Man_CichlidCSizes))
colnames(Man_CichlidCSizes1)<-c("Csize","Phylonames")

TaxaIDData<-read.csv("Microevo_uct1.csv")
Tree<-read.tree("CichlidTreeNature2020.tre")

TaxaIDData<-merge(TaxaIDData, Man_CichlidCSizes1, by = "Phylonames")
rownames(TaxaIDData)<-TaxaIDData$Phylonames

Pruning<-treedata(Tree, TaxaIDData)
Ftree<-Pruning$phy
FTree<-ladderize(Ftree)
FData<-Pruning$data

FData <- FData[match(FTree$tip.label, rownames(FData)),]
FData<-as.data.frame(FData)

RM_Cichlid_Man<-dimnames(Man_CichlidMeanShapes)[[3]][dimnames(Man_CichlidMeanShapes)[[3]]%in%rownames(FData)]
Man_CichlidMeanShapes_Tr<-Man_CichlidMeanShapes[,,RM_Cichlid_Man]


Man_CichlidMeanShapes_Tr<-Man_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(Man_CichlidMeanShapes_Tr)[[3]]))]

gdf_Man <- geomorph.data.frame(coords = Man_CichlidMeanShapes_Tr, species = dimnames(Man_CichlidMeanShapes_Tr)[[3]], Csize = as.numeric(FData$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes)##
HybAnova_Man <- procD.pgls(coords~Csize, phy = FTree, data = gdf_Man) 
summary(HybAnova_Man) 

shape.resid.Man <- arrayspecs(HybAnova_Man$pgls.residuals, p=dim(Man_CichlidMeanShapes_Tr)[1], k=dim(Man_CichlidMeanShapes_Tr)[2]) # size-adjusted residuals
Man_CichlidMeanShapes_Tr1 <- shape.resid.Man + array(Y.gpa_Man$consensus, dim(shape.resid.Man))

#Run PCA on new  size corrected shape data
Man_CichlidMeanShapes_Tr1_PCA.Allom <- gm.prcomp(Man_CichlidMeanShapes_Tr1)

#ploting the size corrected specimen
plotAllSpecimens(Man_CichlidMeanShapes_Tr1)

#plotting the PCA pc1 + 2
plot(Man_CichlidMeanShapes_Tr1_PCA.Allom)
#text(psPCA_CL$x, labels = dimnames(Y.gpa_CL$coords)[[3]])
text(Man_CichlidMeanShapes_Tr1_PCA.Allom$x, pos = 4, label = dimnames(Man_CichlidMeanShapes_Tr1)[[3]], cex=.5)



#CB5
PJ_CichlidMeanShapes<-Y.gpa_PJ$coords
PJ_CichlidCSizes<-as.data.frame(Y.gpa_PJ$Csize)

PJ_CichlidCSizes1<-cbind.data.frame(PJ_CichlidCSizes, rownames(PJ_CichlidCSizes))
colnames(PJ_CichlidCSizes1)<-c("Csize","Phylonames")

TaxaIDData<-read.csv("Microevo_uct1.csv")
Tree<-read.tree("CichlidTreeNature2020.tre")

TaxaIDData<-merge(TaxaIDData, PJ_CichlidCSizes1, by = "Phylonames")
rownames(TaxaIDData)<-TaxaIDData$Phylonames

Pruning<-treedata(Tree, TaxaIDData)
Ftree<-Pruning$phy
FTree<-ladderize(Ftree)
FData<-Pruning$data

FData <- FData[match(FTree$tip.label, rownames(FData)),]
FData<-as.data.frame(FData)

RM_Cichlid_PJ<-dimnames(PJ_CichlidMeanShapes)[[3]][dimnames(PJ_CichlidMeanShapes)[[3]]%in%rownames(FData)]
PJ_CichlidMeanShapes_Tr<-PJ_CichlidMeanShapes[,,RM_Cichlid_PJ]


PJ_CichlidMeanShapes_Tr<-PJ_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(PJ_CichlidMeanShapes_Tr)[[3]]))]

gdf_PJ <- geomorph.data.frame(coords = PJ_CichlidMeanShapes_Tr, species = dimnames(PJ_CichlidMeanShapes_Tr)[[3]], Csize = as.numeric(FData$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes)##
HybAnova_PJ <- procD.pgls(coords~Csize, phy = FTree, data = gdf_PJ) 
summary(HybAnova_PJ) 

shape.resid.PJ <- arrayspecs(HybAnova_PJ$pgls.residuals, p=dim(PJ_CichlidMeanShapes_Tr)[1], k=dim(PJ_CichlidMeanShapes_Tr)[2]) # size-adjusted residuals
PJ_CichlidMeanShapes_Tr1 <- shape.resid.PJ + array(Y.gpa_PJ$consensus, dim(shape.resid.PJ))

#Run PCA on new  size corrected shape data
PJ_CichlidMeanShapes_Tr1_PCA.Allom <- gm.prcomp(PJ_CichlidMeanShapes_Tr1)

#ploting the size corrected specimen
plotAllSpecimens(PJ_CichlidMeanShapes_Tr1)

#plotting the PCA pc1 + 2
plot(PJ_CichlidMeanShapes_Tr1_PCA.Allom)
#text(psPCA_CL$x, labels = dimnames(Y.gpa_CL$coords)[[3]])
text(PJ_CichlidMeanShapes_Tr1_PCA.Allom$x, pos = 4, label = dimnames(PJ_CichlidMeanShapes_Tr1)[[3]], cex=.5)




### MAN on x axis 
####Correlation

# man v ch 
CommonIndsCHMAN<-c(gsub( "_.*$", "", dimnames(CH_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(Man_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(CH_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                             file_path_sans_ext(dimnames(Man_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsCHMAN<-c(dimnames(CH_CichlidMeanShapes_Tr1)[[3]], dimnames(Man_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(CH_CichlidMeanShapes_Tr1)[[3]], dimnames(Man_CichlidMeanShapes_Tr)[[3]]))]


#CH setup
CH_CichlidMeanShapes_TrNEW <- CH_CichlidMeanShapes_Tr1[,,(CommonIndsCHMAN)]
dimnames(CH_CichlidMeanShapes_TrNEW)[[3]]<-CommonIndsCHMAN

#man setup
Man_CichlidMeanShapes_TrNEW<-Man_CichlidMeanShapes_Tr[,,CommonIndsCHMAN]
dimnames(Man_CichlidMeanShapes_TrNEW)[[3]]<-CommonIndsCHMAN

FData$PCHKey = FData$Location
FData$PCHKey[FData$PCHKey == 'Malawi'] <- '15'
FData$PCHKey[FData$PCHKey == 'Tanganyika'] <- '16'
FData$PCHKey[FData$PCHKey == 'Victoria'] <- '17'
FData$PCHKey[FData$PCHKey == 'West'] <- '18'

#This function will be the phylo integration one
ManCHCorrs <- two.b.pls(A1 = Man_CichlidMeanShapes_TrNEW, A2 = CH_CichlidMeanShapes_TrNEW, iter = 9999)
plot(ManCHCorrs,label = CommonIndsCHMAN)
# (Levels of Location data) [1] "Malawi"     "Tanganyika" "Victoria"   "West"  
plot(ManCHCorrs, pch = as.numeric(FData$PCHKey), cex = 4) #xlim=c(-0.3,0.3), ylim=c(-0.3,0.3)
summary(ManCHCorrs)


integration_testMandHYD <- modularity.test(Man_CichlidMeanShapes_TrNEW, CH_CichlidMeanShapes_TrNEW,
                                   iter = 9999)
summary(integration_testMandHYD)
plot(integration_testMandHYD)


#pick and plot 
ManCHCorrs <- two.b.pls(A1 = Man_CichlidMeanShapes_TrNEW, A2 = CH_CichlidMeanShapes_TrNEW, iter = 9999)
ManCHCorrs.plot <- plot(ManCHCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
#picknplot.shape(ManCHCorrs.plot,  method="vector", mg = 4)
ManCHCorrs.plot_shape = picknplot.shape(ManCHCorrs.plot)

#magnifications is mag = 1

#block 1 man: 
plotRefToTarget(ManCHCorrs.plot_shape$shapes[[1]]$P1, ManCHCorrs.plot_shape$shapes[[2]]$P1, method = "vector")

#block 2 ch
plotRefToTarget(ManCHCorrs.plot_shape$shapes[[1]]$P2, ManCHCorrs.plot_shape$shapes[[2]]$P2, method = "vector")




# man v cb1 
CommonIndsCB1MAN<-c(gsub( "_.*$", "", dimnames(CB1_CichlidMeanShapes_Tr1)[[3]]),
                    file_path_sans_ext(dimnames(Man_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(CB1_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                              file_path_sans_ext(dimnames(Man_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsCB1MAN<-c(dimnames(CB1_CichlidMeanShapes_Tr1)[[3]], dimnames(Man_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(CB1_CichlidMeanShapes_Tr1)[[3]], dimnames(Man_CichlidMeanShapes_Tr)[[3]]))]



#CB1 setup
CB1_CichlidMeanShapes_TrNEW <- CB1_CichlidMeanShapes_Tr1[,,(CommonIndsCB1MAN)]
dimnames(CB1_CichlidMeanShapes_TrNEW)[[3]]<-CommonIndsCB1MAN

#man setup
Man_CichlidMeanShapes_TrNEW<-Man_CichlidMeanShapes_Tr[,,CommonIndsCB1MAN]
dimnames(Man_CichlidMeanShapes_TrNEW)[[3]]<-CommonIndsCB1MAN

FData$PCHKey = FData$Location
FData$PCHKey[FData$PCHKey == 'Malawi'] <- '15'
FData$PCHKey[FData$PCHKey == 'Tanganyika'] <- '16'
FData$PCHKey[FData$PCHKey == 'Victoria'] <- '17'
FData$PCHKey[FData$PCHKey == 'West'] <- '18'

#This function will be the phylo integration one
ManCB1Corrs <- two.b.pls(A1 = Man_CichlidMeanShapes_TrNEW, A2 = CB1_CichlidMeanShapes_TrNEW, iter = 9999)
plot(ManCB1Corrs, label = CommonIndsCB1MAN)  #pch = c(15, 16, 17, 18)
plot(ManCB1Corrs,pch = as.numeric(FData$PCHKey), cex = 4)
summary(ManCB1Corrs)

integration_testMandCB1 <- integration.test(Man_CichlidMeanShapes_TrNEW, CB1_CichlidMeanShapes_TrNEW,
                                            iter = 9999)
summary(integration_testMandCB1)
plot(integration_testMandCB1)

#pick and plot
ManCB1Corrs <- two.b.pls(A1 = Man_CichlidMeanShapes_TrNEW, A2 = CB1_CichlidMeanShapes_TrNEW, iter = 9999)
ManCB1Corrs.plot <- plot(ManCB1Corrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
#picknplot.shape(ManCB1Corrs.plot,  method="vector", mg = 4)
ManCB1Corrs.plot_shape = picknplot.shape(ManCB1Corrs.plot)


#block 1 man: 
plotRefToTarget(ManCB1Corrs.plot_shape$shapes[[1]]$P1, ManCB1Corrs.plot_shape$shapes[[2]]$P1, method = "vector")

#block 2 cb1
plotRefToTarget(ManCB1Corrs.plot_shape$shapes[[1]]$P2, ManCB1Corrs.plot_shape$shapes[[2]]$P2, method = "vector")




# man v cb5
CommonIndsPJMAN<-c(gsub( "_.*$", "", dimnames(PJ_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(Man_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(PJ_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                             file_path_sans_ext(dimnames(Man_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsPJMAN<-c(dimnames(PJ_CichlidMeanShapes_Tr1)[[3]], dimnames(Man_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(PJ_CichlidMeanShapes_Tr1)[[3]], dimnames(Man_CichlidMeanShapes_Tr)[[3]]))]



#PJ setup
PJ_CichlidMeanShapes_TrNEW <- PJ_CichlidMeanShapes_Tr1[,,(CommonIndsPJMAN)]
dimnames(PJ_CichlidMeanShapes_TrNEW)[[3]]<-CommonIndsPJMAN

#man setup
Man_CichlidMeanShapes_TrNEW<-Man_CichlidMeanShapes_Tr[,,CommonIndsPJMAN]
dimnames(Man_CichlidMeanShapes_TrNEW)[[3]]<-CommonIndsPJMAN

FData$PCHKey = FData$Location
FData$PCHKey[FData$PCHKey == 'Malawi'] <- '15'
FData$PCHKey[FData$PCHKey == 'Tanganyika'] <- '16'
FData$PCHKey[FData$PCHKey == 'Victoria'] <- '17'
FData$PCHKey[FData$PCHKey == 'West'] <- '18'

#This function will be the phylo integration one
ManPJCorrs <- two.b.pls(A1 = Man_CichlidMeanShapes_TrNEW, A2 = PJ_CichlidMeanShapes_TrNEW, iter = 9999)
plot(ManPJCorrs, label =CommonIndsPJMAN)
plot(ManPJCorrs, pch = as.numeric(FData$PCHKey), cex = 4) #, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
summary(ManPJCorrs)


integration_testMandCB5 <- integration.test(Man_CichlidMeanShapes_TrNEW, PJ_CichlidMeanShapes_TrNEW,
                                            iter = 9999)
summary(integration_testMandCB5)
plot(integration_testMandCB5)

ManPJCorrs <- two.b.pls(A1 = Man_CichlidMeanShapes_TrNEW, A2 = PJ_CichlidMeanShapes_TrNEW, iter = 9999)
ManPJCorrs.plot <- plot(ManPJCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
#picknplot.shape(ManPJCorrs.plot,  method="vector", mg = 4)
ManPJCorrs.plot_shape = picknplot.shape(ManPJCorrs.plot)


#block 1 man: 
plotRefToTarget(ManPJCorrs.plot_shape$shapes[[1]]$P1, ManPJCorrs.plot_shape$shapes[[2]]$P1, method = "vector")

#block 2 cb5
plotRefToTarget(ManPJCorrs.plot_shape$shapes[[1]]$P2, ManPJCorrs.plot_shape$shapes[[2]]$P2, method = "vector")




# man v cl 
CommonIndsCLMAN<-c(gsub( "_.*$", "", dimnames(CL_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(Man_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(CL_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                             file_path_sans_ext(dimnames(Man_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsCLMAN<-c(dimnames(CL_CichlidMeanShapes_Tr1)[[3]], dimnames(Man_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(CL_CichlidMeanShapes_Tr1)[[3]], dimnames(Man_CichlidMeanShapes_Tr)[[3]]))]



#CL setup
CL_CichlidMeanShapes_TrNEW <- CL_CichlidMeanShapes_Tr1[,,(CommonIndsCLMAN)]
dimnames(CL_CichlidMeanShapes_TrNEW)[[3]]<-CommonIndsCLMAN

#man setup
Man_CichlidMeanShapes_TrNEW<-Man_CichlidMeanShapes_Tr[,,CommonIndsCLMAN]
dimnames(Man_CichlidMeanShapes_TrNEW)[[3]]<-CommonIndsCLMAN

FData$PCHKey = FData$Location
FData$PCHKey[FData$PCHKey == 'Malawi'] <- '15'
FData$PCHKey[FData$PCHKey == 'Tanganyika'] <- '16'
FData$PCHKey[FData$PCHKey == 'Victoria'] <- '17'
FData$PCHKey[FData$PCHKey == 'West'] <- '18'

#This function will be the phylo integration one
ManCLCorrs <- two.b.pls(A1 = Man_CichlidMeanShapes_TrNEW, A2 = CL_CichlidMeanShapes_TrNEW, iter = 9999)
plot(ManCLCorrs, label = CommonIndsCLMAN)
plot(ManCLCorrs, pch = as.numeric(FData$PCHKey), cex = 4)
summary(ManCLCorrs)

integration_testMandCL <- integration.test(Man_CichlidMeanShapes_TrNEW, CL_CichlidMeanShapes_TrNEW,
                                            iter = 9999)
summary(integration_testMandCL)
plot(integration_testMandCL)


ManCLCorrs <- two.b.pls(A1 = Man_CichlidMeanShapes_TrNEW, A2 = CL_CichlidMeanShapes_TrNEW, iter = 9999)
ManCLCorrs.plot <- plot(ManCLCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
#picknplot.shape(ManCLCorrs.plot,  method="vector", mg = 4)
ManCLCorrs.plot_shape = picknplot.shape(ManCLCorrs.plot)


#block 1 man: 
plotRefToTarget(ManCLCorrs.plot_shape$shapes[[1]]$P1, ManCLCorrs.plot_shape$shapes[[2]]$P1, method = "vector")

#block 2 cl
plotRefToTarget(ManCLCorrs.plot_shape$shapes[[1]]$P2, ManCLCorrs.plot_shape$shapes[[2]]$P2, method = "vector")



### CL on x axis 
####Correlation

# cl v cb5
CommonIndsPJCL<-c(gsub( "_.*$", "", dimnames(PJ_CichlidMeanShapes_Tr1)[[3]]),
                  file_path_sans_ext(dimnames(CL_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(PJ_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                           file_path_sans_ext(dimnames(CL_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsPJCL<-c(dimnames(PJ_CichlidMeanShapes_Tr1)[[3]], dimnames(CL_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(PJ_CichlidMeanShapes_Tr1)[[3]], dimnames(CL_CichlidMeanShapes_Tr)[[3]]))]



#CL setup
CL_CichlidMeanShapes_TrNEW_4CL <- CL_CichlidMeanShapes_Tr1[,,(CommonIndsPJCL)]
dimnames(CL_CichlidMeanShapes_TrNEW_4CL)[[3]]<-CommonIndsPJCL

#PJ setup
PJ_CichlidMeanShapes_TrNEW_4CL<-PJ_CichlidMeanShapes_Tr[,,CommonIndsPJCL]
dimnames(PJ_CichlidMeanShapes_TrNEW_4CL)[[3]]<-CommonIndsPJCL


FData$PCHKey = FData$Location
FData$PCHKey[FData$PCHKey == 'Malawi'] <- '15'
FData$PCHKey[FData$PCHKey == 'Tanganyika'] <- '16'
FData$PCHKey[FData$PCHKey == 'Victoria'] <- '17'
FData$PCHKey[FData$PCHKey == 'West'] <- '18'

#This function will be the phylo integration one
CLPJCorrs <- two.b.pls(A1 = CL_CichlidMeanShapes_TrNEW_4CL, A2 = PJ_CichlidMeanShapes_TrNEW_4CL, iter = 9999)
plot(CLPJCorrs, label = CommonIndsPJCL)
plot(CLPJCorrs, pch = as.numeric(FData$PCHKey), cex = 4)
summary(CLPJCorrs)


integration_testCLCB5 <- integration.test(CL_CichlidMeanShapes_TrNEW_4CL, PJ_CichlidMeanShapes_TrNEW_4CL,
                                           iter = 9999)
summary(integration_testCLCB5)
plot(integration_testCLCB5)


CLPJCorrs <- two.b.pls(A1 = CL_CichlidMeanShapes_TrNEW_4CL, A2 = PJ_CichlidMeanShapes_TrNEW_4CL, iter = 9999)
CLPJCorrs.plot <- plot(CLPJCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
#picknplot.shape(CLPJCorrs.plot,  method="vector", mg = 4)
CLPJCorrs.plot_shape = picknplot.shape(CLPJCorrs.plot)


#block 1 cl: 
plotRefToTarget(CLPJCorrs.plot_shape$shapes[[1]]$P1, CLPJCorrs.plot_shape$shapes[[2]]$P1, method = "vector")

#block 2 cb5
plotRefToTarget(CLPJCorrs.plot_shape$shapes[[1]]$P2, CLPJCorrs.plot_shape$shapes[[2]]$P2, method = "vector")





# cl v cb1
CommonIndsCB1CL<-c(gsub( "_.*$", "", dimnames(CB1_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(CL_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(CB1_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                            file_path_sans_ext(dimnames(CL_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsCB1CL<-c(dimnames(CB1_CichlidMeanShapes_Tr1)[[3]], dimnames(CL_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(CB1_CichlidMeanShapes_Tr1)[[3]], dimnames(CL_CichlidMeanShapes_Tr)[[3]]))]



#CL setup
CL_CichlidMeanShapes_TrNEW_4CL <- CL_CichlidMeanShapes_Tr1[,,(CommonIndsCB1CL)]
dimnames(CL_CichlidMeanShapes_TrNEW_4CL)[[3]]<-CommonIndsCB1CL

#CB1 setup
CB1_CichlidMeanShapes_TrNEW_4CL<-CB1_CichlidMeanShapes_Tr[,,CommonIndsCB1CL]
dimnames(CB1_CichlidMeanShapes_TrNEW_4CL)[[3]]<-CommonIndsCB1CL


FData$PCHKey = FData$Location
FData$PCHKey[FData$PCHKey == 'Malawi'] <- '15'
FData$PCHKey[FData$PCHKey == 'Tanganyika'] <- '16'
FData$PCHKey[FData$PCHKey == 'Victoria'] <- '17'
FData$PCHKey[FData$PCHKey == 'West'] <- '18'

#This function will be the phylo integration one
CLCB1Corrs <- two.b.pls(A1 = CL_CichlidMeanShapes_TrNEW_4CL, A2 = CB1_CichlidMeanShapes_TrNEW_4CL, iter = 9999)
plot(CLCB1Corrs, label = CommonIndsCB1CL)
plot(CLCB1Corrs, pch = as.numeric(FData$PCHKey), cex = 4)
summary(CLCB1Corrs)

integration_testCLCB1 <- integration.test(CL_CichlidMeanShapes_TrNEW_4CL, CB1_CichlidMeanShapes_TrNEW_4CL,
                                          iter = 9999)
summary(integration_testCLCB1)
plot(integration_testCLCB1)



CLCB1Corrs <- two.b.pls(A1 = CL_CichlidMeanShapes_TrNEW_4CL, A2 = CB1_CichlidMeanShapes_TrNEW_4CL, iter = 9999)
CLCB1Corrs.plot <- plot(CLCB1Corrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
#picknplot.shape(CLCB1Corrs.plot,  method="vector", mg = 4)
CLCB1Corrs.plot_shape = picknplot.shape(CLCB1Corrs.plot)


#block 1 cl: 
plotRefToTarget(CLCB1Corrs.plot_shape$shapes[[1]]$P1, CLCB1Corrs.plot_shape$shapes[[2]]$P1, method = "vector")

#block 2 cb1
plotRefToTarget(CLCB1Corrs.plot_shape$shapes[[1]]$P2, CLCB1Corrs.plot_shape$shapes[[2]]$P2, method = "vector")





# cl v ch 
CommonIndsCHCL<-c(gsub( "_.*$", "", dimnames(CH_CichlidMeanShapes_Tr1)[[3]]),
                  file_path_sans_ext(dimnames(CL_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(CH_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                           file_path_sans_ext(dimnames(CL_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsCHCL<-c(dimnames(CH_CichlidMeanShapes_Tr1)[[3]], dimnames(CL_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(CH_CichlidMeanShapes_Tr1)[[3]], dimnames(CL_CichlidMeanShapes_Tr)[[3]]))]



#CL setup
CL_CichlidMeanShapes_TrNEW_4CL <- CL_CichlidMeanShapes_Tr1[,,(CommonIndsCHCL)]
dimnames(CL_CichlidMeanShapes_TrNEW_4CL)[[3]]<-CommonIndsCHCL

#CH setup
CH_CichlidMeanShapes_TrNEW_4CL<-CH_CichlidMeanShapes_Tr[,,CommonIndsCHCL]
dimnames(CH_CichlidMeanShapes_TrNEW_4CL)[[3]]<-CommonIndsCHCL

FData$PCHKey = FData$Location
FData$PCHKey[FData$PCHKey == 'Malawi'] <- '15'
FData$PCHKey[FData$PCHKey == 'Tanganyika'] <- '16'
FData$PCHKey[FData$PCHKey == 'Victoria'] <- '17'
FData$PCHKey[FData$PCHKey == 'West'] <- '18'

#This function will be the phylo integration one
CLCHCorrs <- two.b.pls(A1 = CL_CichlidMeanShapes_TrNEW_4CL, A2 = CH_CichlidMeanShapes_TrNEW_4CL, iter = 9999)
plot(CLCHCorrs, label = CommonIndsCHCL)
plot(CLCHCorrs, pch = as.numeric(FData$PCHKey), cex = 4) #, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
summary(CLCHCorrs)

integration_testCLCH <- integration.test(CL_CichlidMeanShapes_TrNEW_4CL, CH_CichlidMeanShapes_TrNEW_4CL,
                                          iter = 9999)
summary(integration_testCLCH)
plot(integration_testCLCH)

CLCHCorrs <- two.b.pls(A1 = CL_CichlidMeanShapes_TrNEW_4CL, A2 = CH_CichlidMeanShapes_TrNEW_4CL, iter = 9999)
CLCHCorrs.plot <- plot(CLCHCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
#picknplot.shape(CLCHCorrs.plot,  method="vector", mg = 4)
CLCHCorrs.plot_shape = picknplot.shape(CLCHCorrs.plot)


#block 1 cl: 
plotRefToTarget(CLCHCorrs.plot_shape$shapes[[1]]$P1, CLCHCorrs.plot_shape$shapes[[2]]$P1, method = "vector")

#block 2 ch
plotRefToTarget(CLCHCorrs.plot_shape$shapes[[1]]$P2, CLCHCorrs.plot_shape$shapes[[2]]$P2, method = "vector")





# cl v man 
CommonIndsMANCL<-c(gsub( "_.*$", "", dimnames(Man_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(CL_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(Man_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                            file_path_sans_ext(dimnames(CL_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsMANCL<-c(dimnames(Man_CichlidMeanShapes_Tr1)[[3]], dimnames(CL_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(Man_CichlidMeanShapes_Tr1)[[3]], dimnames(CL_CichlidMeanShapes_Tr)[[3]]))]



#CL setup
CL_CichlidMeanShapes_TrNEW_4CL <- CL_CichlidMeanShapes_Tr1[,,(CommonIndsMANCL)]
dimnames(CL_CichlidMeanShapes_TrNEW_4CL)[[3]]<-CommonIndsMANCL

#man setup
Man_CichlidMeanShapes_TrNEW_4CL<-Man_CichlidMeanShapes_Tr[,,CommonIndsMANCL]
dimnames(Man_CichlidMeanShapes_TrNEW_4CL)[[3]]<-CommonIndsMANCL

FData$PCHKey = FData$Location
FData$PCHKey[FData$PCHKey == 'Malawi'] <- '15'
FData$PCHKey[FData$PCHKey == 'Tanganyika'] <- '16'
FData$PCHKey[FData$PCHKey == 'Victoria'] <- '17'
FData$PCHKey[FData$PCHKey == 'West'] <- '18'

#This function will be the phylo integration one
CLManCorrs <- two.b.pls(A1 = CL_CichlidMeanShapes_TrNEW_4CL, A2 = Man_CichlidMeanShapes_TrNEW_4CL, iter = 9999)
plot(CLManCorrs, label = CommonIndsMANCL)
plot(CLManCorrs, pch = as.numeric(FData$PCHKey), cex = 4)
summary(CLManCorrs)

integration_testCLMan <- integration.test(CL_CichlidMeanShapes_TrNEW_4CL, Man_CichlidMeanShapes_TrNEW_4CL,
                                          iter = 9999)
summary(integration_testCLMan)
plot(integration_testCLMan)


CLManCorrs <- two.b.pls(A1 = CL_CichlidMeanShapes_TrNEW_4CL, A2 = Man_CichlidMeanShapes_TrNEW_4CL, iter = 9999)
CLManCorrs.plot <- plot(CLManCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
#picknplot.shape(CLManCorrs.plot,  method="vector", mg = 4)
CLManCorrs.plot_shape = picknplot.shape(CLManCorrs.plot)


#block 1 cl: 
plotRefToTarget(CLManCorrs.plot_shape$shapes[[1]]$P1, CLManCorrs.plot_shape$shapes[[2]]$P1, method = "vector")

#block 2 man
plotRefToTarget(CLManCorrs.plot_shape$shapes[[1]]$P2, CLManCorrs.plot_shape$shapes[[2]]$P2, method = "vector")




### HYOID on x axis 
####Correlation

# HYD vs Mand 
CommonIndsMANCH<-c(gsub( "_.*$", "", dimnames(Man_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(CH_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(Man_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                             file_path_sans_ext(dimnames(CH_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsMANCH<-c(dimnames(Man_CichlidMeanShapes_Tr1)[[3]], dimnames(CH_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(Man_CichlidMeanShapes_Tr1)[[3]], dimnames(CH_CichlidMeanShapes_Tr)[[3]]))]


#CH setup
CH_CichlidMeanShapes_TrNEW_4CH <- CH_CichlidMeanShapes_Tr1[,,(CommonIndsMANCH)]
dimnames(CH_CichlidMeanShapes_TrNEW_4CH)[[3]]<-CommonIndsMANCH

#man setup
Man_CichlidMeanShapes_TrNEW_4CH <-Man_CichlidMeanShapes_Tr[,,CommonIndsMANCH]
dimnames(Man_CichlidMeanShapes_TrNEW_4CH)[[3]]<-CommonIndsMANCH

#This function will be the phylo integration one
CHMANCorrs <- two.b.pls(A1 = CH_CichlidMeanShapes_TrNEW_4CH, A2 = Man_CichlidMeanShapes_TrNEW_4CH, iter = 9999)
plot(CHMANCorrs,label = CommonIndsMANCH)
plot(CHMANCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
summary(CHMANCorrs)

#integration_testMandHYD <- modularity.test(Man_CichlidMeanShapes_TrNEW, CH_CichlidMeanShapes_TrNEW,iter = 9999)
#summary(integration_testMandHYD)
#plot(integration_testMandHYD)


#pick and plot 
#ManCHCorrs <- two.b.pls(A1 = Man_CichlidMeanShapes_TrNEW, A2 = CH_CichlidMeanShapes_TrNEW, iter = 9999)
#ManCHCorrs.plot <- plot(ManCHCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
#picknplot.shape(ManCHCorrs.plot,  method="vector", mg = 4)
#ManCHCorrs.plot_shape = picknplot.shape(ManCHCorrs.plot)
#magnifications is mag = 1

#block 1 man: 
#plotRefToTarget(ManCHCorrs.plot_shape$shapes[[1]]$P1, ManCHCorrs.plot_shape$shapes[[2]]$P1, method = "vector")

#block 2 ch
#plotRefToTarget(ManCHCorrs.plot_shape$shapes[[1]]$P2, ManCHCorrs.plot_shape$shapes[[2]]$P2, method = "vector")



 
# HYD vs CB1 
CommonIndsCB1CH<-c(gsub( "_.*$", "", dimnames(CB1_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(CH_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(CB1_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                            file_path_sans_ext(dimnames(CH_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsCB1CH<-c(dimnames(CB1_CichlidMeanShapes_Tr1)[[3]], dimnames(CH_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(CB1_CichlidMeanShapes_Tr1)[[3]], dimnames(CH_CichlidMeanShapes_Tr)[[3]]))]


#CH setup
CH_CichlidMeanShapes_TrNEW_4CH <- CH_CichlidMeanShapes_Tr1[,,(CommonIndsCB1CH)]
dimnames(CH_CichlidMeanShapes_TrNEW_4CH)[[3]]<-CommonIndsCB1CH

#CB1 setup
CB1_CichlidMeanShapes_TrNEW_4CH <-CB1_CichlidMeanShapes_Tr[,,CommonIndsCB1CH]
dimnames(CB1_CichlidMeanShapes_TrNEW_4CH)[[3]]<-CommonIndsCB1CH

#This function will be the phylo integration one
CHCB1Corrs <- two.b.pls(A1 = CH_CichlidMeanShapes_TrNEW_4CH, A2 = CB1_CichlidMeanShapes_TrNEW_4CH, iter = 9999)
plot(CHCB1Corrs,label = CommonIndsCB1CH)
plot(CHCB1Corrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
summary(CHCB1Corrs)



# HYD vs CB5
CommonIndsCB5CH<-c(gsub( "_.*$", "", dimnames(PJ_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(CH_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(PJ_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                            file_path_sans_ext(dimnames(CH_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsCB5CH<-c(dimnames(PJ_CichlidMeanShapes_Tr1)[[3]], dimnames(CH_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(PJ_CichlidMeanShapes_Tr1)[[3]], dimnames(CH_CichlidMeanShapes_Tr)[[3]]))]


#CH setup
CH_CichlidMeanShapes_TrNEW_4CH <- CH_CichlidMeanShapes_Tr1[,,(CommonIndsCB5CH)]
dimnames(CH_CichlidMeanShapes_TrNEW_4CH)[[3]]<-CommonIndsCB5CH

#CB5 setup
CB5_CichlidMeanShapes_TrNEW_4CH <-PJ_CichlidMeanShapes_Tr[,,CommonIndsCB5CH]
dimnames(CB5_CichlidMeanShapes_TrNEW_4CH)[[3]]<-CommonIndsCB5CH

#This function will be the phylo integration one
CHCB5Corrs <- two.b.pls(A1 = CH_CichlidMeanShapes_TrNEW_4CH, A2 = CB5_CichlidMeanShapes_TrNEW_4CH, iter = 9999)
plot(CHCB5Corrs,label = CommonIndsCB5CH)
plot(CHCB5Corrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
summary(CHCB5Corrs)



# HYD vs CL
CommonIndsCLCH<-c(gsub( "_.*$", "", dimnames(CL_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(CH_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(CL_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                            file_path_sans_ext(dimnames(CH_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsCLCH<-c(dimnames(CL_CichlidMeanShapes_Tr1)[[3]], dimnames(CH_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(CL_CichlidMeanShapes_Tr1)[[3]], dimnames(CH_CichlidMeanShapes_Tr)[[3]]))]


#CH setup
CH_CichlidMeanShapes_TrNEW_4CH <- CH_CichlidMeanShapes_Tr1[,,(CommonIndsCLCH)]
dimnames(CH_CichlidMeanShapes_TrNEW_4CH)[[3]]<-CommonIndsCLCH

#CL setup
CL_CichlidMeanShapes_TrNEW_4CH <-CL_CichlidMeanShapes_Tr[,,CommonIndsCLCH]
dimnames(CL_CichlidMeanShapes_TrNEW_4CH)[[3]]<-CommonIndsCLCH

#This function will be the phylo integration one
CHCLCorrs <- two.b.pls(A1 = CH_CichlidMeanShapes_TrNEW_4CH, A2 = CL_CichlidMeanShapes_TrNEW_4CH, iter = 9999)
plot(CHCLCorrs,label = CommonIndsCLCH)
plot(CHCLCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
summary(CHCLCorrs)


### CB1 on x axis 
####Correlation

# CB1 vs Hyd 
CommonIndsCHCB1<-c(gsub( "_.*$", "", dimnames(CH_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(CB1_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(CH_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                            file_path_sans_ext(dimnames(CB1_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsCHCB1<-c(dimnames(CH_CichlidMeanShapes_Tr1)[[3]], dimnames(CB1_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(CH_CichlidMeanShapes_Tr1)[[3]], dimnames(CB1_CichlidMeanShapes_Tr)[[3]]))]


#CH setup
CH_CichlidMeanShapes_TrNEW_4CB1 <- CH_CichlidMeanShapes_Tr1[,,(CommonIndsCHCB1)]
dimnames(CH_CichlidMeanShapes_TrNEW_4CB1)[[3]]<-CommonIndsCHCB1

#CB1 setup
CB1_CichlidMeanShapes_TrNEW_4CB1 <-CB1_CichlidMeanShapes_Tr[,,CommonIndsCHCB1]
dimnames(CB1_CichlidMeanShapes_TrNEW_4CB1)[[3]]<-CommonIndsCHCB1

#This function will be the phylo integration one
CB1CHCorrs <- two.b.pls(A1 = CB1_CichlidMeanShapes_TrNEW_4CB1, A2 = CH_CichlidMeanShapes_TrNEW_4CB1, iter = 9999)
plot(CB1CHCorrs,label = CommonIndsCHCB1)
plot(CB1CHCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
summary(CB1CHCorrs)


# CB1 vs Mand
CommonIndsManCB1<-c(gsub( "_.*$", "", dimnames(Man_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(CB1_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(Man_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                             file_path_sans_ext(dimnames(CB1_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsManCB1<-c(dimnames(Man_CichlidMeanShapes_Tr1)[[3]], dimnames(CB1_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(Man_CichlidMeanShapes_Tr1)[[3]], dimnames(CB1_CichlidMeanShapes_Tr)[[3]]))]


#Man setup
Man_CichlidMeanShapes_TrNEW_4CB1 <- Man_CichlidMeanShapes_Tr1[,,(CommonIndsManCB1)]
dimnames(Man_CichlidMeanShapes_TrNEW_4CB1)[[3]]<-CommonIndsManCB1

#CB1 setup
CB1_CichlidMeanShapes_TrNEW_4CB1 <-CB1_CichlidMeanShapes_Tr[,,CommonIndsManCB1]
dimnames(CB1_CichlidMeanShapes_TrNEW_4CB1)[[3]]<-CommonIndsManCB1

#This function will be the phylo integration one
CB1ManCorrs <- two.b.pls(A1 = CB1_CichlidMeanShapes_TrNEW_4CB1, A2 = Man_CichlidMeanShapes_TrNEW_4CB1, iter = 9999)
plot(CB1ManCorrs,label = CommonIndsManCB1)
plot(CB1ManCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
summary(CB1ManCorrs)


# CB1 vs CB5 
CommonIndsCB5CB1<-c(gsub( "_.*$", "", dimnames(PJ_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(CB1_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(PJ_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                             file_path_sans_ext(dimnames(CB1_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsCB5CB1<-c(dimnames(PJ_CichlidMeanShapes_Tr1)[[3]], dimnames(CB1_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(PJ_CichlidMeanShapes_Tr1)[[3]], dimnames(CB1_CichlidMeanShapes_Tr)[[3]]))]


#CB5 setup
CB5_CichlidMeanShapes_TrNEW_4CB1 <- PJ_CichlidMeanShapes_Tr1[,,(CommonIndsCB5CB1)]
dimnames(CB5_CichlidMeanShapes_TrNEW_4CB1)[[3]]<-CommonIndsCB5CB1

#CB1 setup
CB1_CichlidMeanShapes_TrNEW_4CB1 <-CB1_CichlidMeanShapes_Tr[,,CommonIndsCHCB1]
dimnames(CB1_CichlidMeanShapes_TrNEW_4CB1)[[3]]<-CommonIndsCB5CB1

#This function will be the phylo integration one
CB1CB5Corrs <- two.b.pls(A1 = CB1_CichlidMeanShapes_TrNEW_4CB1, A2 = CB5_CichlidMeanShapes_TrNEW_4CB1, iter = 9999)
plot(CB1CB5Corrs,label = CommonIndsCB5CB1)
plot(CB1CB5Corrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
summary(CB1CB5Corrs)


# CB1 vs CL 
CommonIndsCLCB1<-c(gsub( "_.*$", "", dimnames(CL_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(CB1_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(CL_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                             file_path_sans_ext(dimnames(CB1_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsCLCB1<-c(dimnames(CL_CichlidMeanShapes_Tr1)[[3]], dimnames(CB1_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(CL_CichlidMeanShapes_Tr1)[[3]], dimnames(CB1_CichlidMeanShapes_Tr)[[3]]))]


#CL setup
CL_CichlidMeanShapes_TrNEW_4CB1 <- CL_CichlidMeanShapes_Tr1[,,(CommonIndsCLCB1)]
dimnames(CL_CichlidMeanShapes_TrNEW_4CB1)[[3]]<-CommonIndsCLCB1

#CB1 setup
CB1_CichlidMeanShapes_TrNEW_4CB1 <-CB1_CichlidMeanShapes_Tr[,,CommonIndsCLCB1]
dimnames(CB1_CichlidMeanShapes_TrNEW_4CB1)[[3]]<-CommonIndsCLCB1

#This function will be the phylo integration one
CB1CLCorrs <- two.b.pls(A1 = CB1_CichlidMeanShapes_TrNEW_4CB1, A2 = CL_CichlidMeanShapes_TrNEW_4CB1, iter = 9999)
plot(CB1CLCorrs,label = CommonIndsCLCB1)
plot(CB1CLCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
summary(CB1CLCorrs)


### CB5 on x axis 
####Correlation

# CB5 vs Hyd 
CommonIndsCHCB5<-c(gsub( "_.*$", "", dimnames(CH_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(PJ_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(CH_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                             file_path_sans_ext(dimnames(PJ_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsCHCB5<-c(dimnames(CH_CichlidMeanShapes_Tr1)[[3]], dimnames(PJ_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(CH_CichlidMeanShapes_Tr1)[[3]], dimnames(PJ_CichlidMeanShapes_Tr)[[3]]))]


#CH setup
CH_CichlidMeanShapes_TrNEW_4CB5 <- CH_CichlidMeanShapes_Tr1[,,(CommonIndsCHCB5)]
dimnames(CH_CichlidMeanShapes_TrNEW_4CB5)[[3]]<-CommonIndsCHCB5

#CB5 setup
CB5_CichlidMeanShapes_TrNEW_4CB5 <-PJ_CichlidMeanShapes_Tr[,,CommonIndsCHCB5]
dimnames(CB5_CichlidMeanShapes_TrNEW_4CB5)[[3]]<-CommonIndsCHCB5

#This function will be the phylo integration one
CB5CHCorrs <- two.b.pls(A1 = CB5_CichlidMeanShapes_TrNEW_4CB5, A2 = CH_CichlidMeanShapes_TrNEW_4CB5, iter = 9999)
plot(CB5CHCorrs,label = CommonIndsCHCB5)
plot(CB5CHCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
summary(CB5CHCorrs)



# CB5 vs Mand
CommonIndsManCB5<-c(gsub( "_.*$", "", dimnames(Man_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(PJ_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(Man_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                            file_path_sans_ext(dimnames(PJ_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsManCB5<-c(dimnames(Man_CichlidMeanShapes_Tr1)[[3]], dimnames(PJ_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(Man_CichlidMeanShapes_Tr1)[[3]], dimnames(PJ_CichlidMeanShapes_Tr)[[3]]))]


#Man setup
Man_CichlidMeanShapes_TrNEW_4CB5 <- Man_CichlidMeanShapes_Tr1[,,(CommonIndsManCB5)]
dimnames(Man_CichlidMeanShapes_TrNEW_4CB5)[[3]]<-CommonIndsManCB5

#CB5 setup
CB5_CichlidMeanShapes_TrNEW_4CB5 <-PJ_CichlidMeanShapes_Tr[,,CommonIndsManCB5]
dimnames(CB5_CichlidMeanShapes_TrNEW_4CB5)[[3]]<-CommonIndsManCB5

#This function will be the phylo integration one
CB5ManCorrs <- two.b.pls(A1 = CB5_CichlidMeanShapes_TrNEW_4CB5, A2 = Man_CichlidMeanShapes_TrNEW_4CB5, iter = 9999)
plot(CB5ManCorrs,label = CommonIndsManCB5)
plot(CB5ManCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
summary(CB5ManCorrs)



# CB5 vs CB1 
CommonIndsCB1CB5<-c(gsub( "_.*$", "", dimnames(CB1_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(PJ_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(CB1_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                            file_path_sans_ext(dimnames(PJ_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsCB1CB5<-c(dimnames(CB1_CichlidMeanShapes_Tr1)[[3]], dimnames(PJ_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(CB1_CichlidMeanShapes_Tr1)[[3]], dimnames(PJ_CichlidMeanShapes_Tr)[[3]]))]


#CB1 setup
CB1_CichlidMeanShapes_TrNEW_4CB5 <- CB1_CichlidMeanShapes_Tr1[,,(CommonIndsCB1CB5)]
dimnames(CB1_CichlidMeanShapes_TrNEW_4CB5)[[3]]<-CommonIndsCB1CB5

#CB5 setup
CB5_CichlidMeanShapes_TrNEW_4CB5 <-PJ_CichlidMeanShapes_Tr[,,CommonIndsCB1CB5]
dimnames(CB5_CichlidMeanShapes_TrNEW_4CB5)[[3]]<-CommonIndsCB1CB5

#This function will be the phylo integration one
CB5CB1Corrs <- two.b.pls(A1 = CB5_CichlidMeanShapes_TrNEW_4CB5, A2 = CB1_CichlidMeanShapes_TrNEW_4CB5, iter = 9999)
plot(CB5CB1Corrs,label = CommonIndsCB1CB5)
plot(CB5CB1Corrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
summary(CB5CB1Corrs)



# CB5 vs Clth
CommonIndsCLCB5<-c(gsub( "_.*$", "", dimnames(CL_CichlidMeanShapes_Tr1)[[3]]),
                   file_path_sans_ext(dimnames(PJ_CichlidMeanShapes_Tr)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(CL_CichlidMeanShapes_Tr1)[[3]]), 
                                                                                            file_path_sans_ext(dimnames(PJ_CichlidMeanShapes_Tr)[[3]])))]

CommonIndsCLCB5<-c(dimnames(CL_CichlidMeanShapes_Tr1)[[3]], dimnames(PJ_CichlidMeanShapes_Tr)[[3]])[duplicated(c(dimnames(CL_CichlidMeanShapes_Tr1)[[3]], dimnames(PJ_CichlidMeanShapes_Tr)[[3]]))]


#CL setup
CL_CichlidMeanShapes_TrNEW_4CB5 <- CL_CichlidMeanShapes_Tr1[,,(CommonIndsCLCB5)]
dimnames(CL_CichlidMeanShapes_TrNEW_4CB5)[[3]]<-CommonIndsCLCB5

#CB5 setup
CB5_CichlidMeanShapes_TrNEW_4CB5 <-PJ_CichlidMeanShapes_Tr[,,CommonIndsCLCB5]
dimnames(CB5_CichlidMeanShapes_TrNEW_4CB5)[[3]]<-CommonIndsCLCB5

#This function will be the phylo integration one
CB5CLCorrs <- two.b.pls(A1 = CB5_CichlidMeanShapes_TrNEW_4CB5, A2 = CL_CichlidMeanShapes_TrNEW_4CB5, iter = 9999)
plot(CB5CLCorrs,label = CommonIndsCLCB5)
plot(CB5CLCorrs, pch=19, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
summary(CB5CLCorrs)



-----
###heat map###
data <- data.frame(
  dims = c("Mandible", "Cleithrum"),
  Mandible = c(0, 2.26159),
  Hyoid = c(2.54064, 3.614),
  CB1 = c(2.8326, 0.16117),
  CB5 = c(3.00013, 2.73384),
  Cleithrum = c(2.26371, 0)
)


install.packages("reshape2")
library(reshape2)

# Reshape the data for heatmap
data_melted <- melt(data, id.vars = "dims", variable.name = "variables", value.name = "value")


install.packages("ggplot2")
library(ggplot2)

# Create the heatmap
ggplot(data_melted, aes(x = dims, y = variables, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "cornflowerblue", mid = "khaki1", high = "palevioletred2", midpoint = 2.0, 
                       limits = c(0.0, 4.0), oob = scales::squish) +
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

ggsave("across_cichlids_heatmap.pdf", width = 5, height = 5, units = "cm", dpi = 300)

