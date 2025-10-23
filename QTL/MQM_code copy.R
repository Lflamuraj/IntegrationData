rm(list = ls(all = TRUE))
install.packages("qtl")
library(qtl)



data<-read.cross(format="csv", file="geno_pheno_MAN_PJ.csv", na.strings=c("NA"), genotypes=c("AA", "AB", "BB", "D", "C"), convertXdata=TRUE)
summary(data)
geno.image(data)

#mqm mapping won't work with missing genotypes, these need to be estimated
#augdata1<-mqmaugment(data, minprob=0.02) #this will take a while, and might not work if too many missing genotypes
#geno.image(augdata1)

simdata<-fill.geno(data,method = c("imp")) #computationally straightforward, and honestly closer to the actual data
geno.image(simdata)

#mqm scan without cofactors, use as a first pass
scan1<-mqmscan(simdata,pheno.col=8,plot=T,model="dominance")  
summary(scan1)
plot(scan1, chr=7) #if you just want to see a single linkage group
write.csv(summary(scan1),"scan1_man_pc1.csv")
scan1 #displays LOD scores at every marker

#choose cofactors for the final qtl model based on scan1 results; note, that this can be an iterative process!
cofactorsindex<-NULL
cofactorsindex<-c(cofactorsindex, find.markerindex(simdata,find.marker(simdata,chr=3,pos=10)))
cofactorsindex

homemadecofactors<-mqmsetcofactors(simdata,cofactors=c(cofactorsindex))
homemadecofactors

#drop without nulling.. add back in change to 1 
#homemadecofactors[63]<-1
#homemadecofactors[686]<-0

#mqm scan with cofactors
homescan1<-mqmscan(cross=simdata,cofactors=homemadecofactors,pheno.col=5,cofactor.significance=0.05, verbose=T, plot=T, model="dominance") #using a dominance model tests for both additive effects and dominance, i.e., AA+AB vs BB and BB+AB vs AA, default is additive only
summary(homescan1)
plot(homescan1,chr=7)
#change back to home scan... have to change to scan1 for clth pc 1
 bayesint(homescan1,7,qtl.index=1,prob=0.95,lodcolumn=1, expandtomarkers=T) #bayes confidence interval; new school
#lodint(homescan1,15,qtl.index=1,drop=1.5,lodcolumn=1, expandtomarkers=T) #confidence interval via lod-drop; old school

#marker effects
find.marker(simdata,chr=19,pos=20)  #finding marker associated with qtl pos
effect<-effectplot(simdata,pheno.col=8,mname1="scaffold_60_1756320",var.flag=c("group")) #effect plot at qtl
data.frame(effect) #give you the AA AB and BB values

effect2<-effectplot(simdata,pheno.col=8,mname1="scaffold_55_1047278", mname2="scaffold_155_1520985",var.flag=c("group")) #effect plot at 2 qtl
data.frame(effect2)



#determine LOD threshold
result<-mqmpermutation(cross=simdata,scanfunction=mqmscan,cofactors=homemadecofactors,pheno.col=8,n.perm=1000,plot=F,verbose=T,model="dominance")
resultqtl<-mqmprocesspermutation(result)
summary(resultqtl)
mqmplot.permutations(result)

#choose cofactors -- two automated procedures, you can also do this by hand, see below
#autocofactors<-mqmautocofactors(simdata,142) #total number across the genome, generally use this one as we have an uneven marker distribution
#mqmplot.cofactors(simdata, autocofactors, justdots=TRUE)
#autocofactors

#setcofactors <- mqmsetcofactors(simdata, 5) #every X number of markers
#mqmplot.cofactors(simdata, setcofactors, justdots=TRUE)

#mqm scan with cofactors
#homescan1 <- mqmscan(simdata, autocofactors, cofactor.significance=0.02)

#homescan1<-mqmscan(cross=simdata,cofactors=autocofactors,pheno.col=5,cofactor.significance=0.02, verbose=T, plot=T, model="dominance") #using a dominance model tests for both additive effects and dominance, i.e., AA+AB vs BB and BB+AB vs AA, default is additive only
#summary(homescan1)
#plot(homescan1,chr=24)
#?bayesint
#bayesint(homescan1,24,qtl.index=1,prob=0.95,lodcolumn=1, expandtomarkers=T)
#?lodint
#lodint(homescan1,15,qtl.index=1,drop=1.5,lodcolumn=1, expandtomarkers=T)
#homescan1


#bayesint(scan1,7,qtl.index=1,prob=0.95,lodcolumn=1, expandtomarkers=T) #interval before cofactor analysis

#find.marker(augdata1,chr=11,pos=30)  #finding markers for to be added to list of cofactors, there is code for automated selection, but using putative qtl from j/qtl here

#cofactorsindex<-NULL
#cofactorsindex<-c(cofactorsindex,find.markerindex(augdata1,find.marker(augdata1,chr=24,pos=20)))
#cofactorsindex

#homemadecofactors<-mqmsetcofactors(augdata1,cofactors=c(cofactorsindex))
#homemadecofactors
#homemadecofactors[919]<-1  #can change as needed, 1=cofactor, 0=no cofactor
#homemadecofactors
#homemadecofactors[437]<-0
#homemadecofactors

#mqmplot.cofactors(simdata,homemadecofactors,justdots=T)



#find.marker(augdata1,chr=24,pos=20)  #finding markers for to be added to list of cofactors, there is code for automated selection, but using putative qtl from j/qtl here
#effect<-effectplot(augdata1,pheno.col=3,mname1="c13.3714129",var.flag=c("group")) #effect plot at qtl
#data.frame(effect)

#effect2<-effectplot(augdata1,pheno.col=3,mname1="c9.6361466", mname2="c6.2501428",var.flag=c("group")) #effect plot at 2 qtl
#data.frame(effect2)

#
PVE<-(1-(10^-((2/190)*(4.397))))
PVE

#result<-mqmpermutation(cross=augdata1,scanfunction=mqmscan,cofactors=homemadecofactors,pheno.col=3,n.perm=1000,plot=F,verbose=T,model="dominance")
#resultqtl<-mqmprocesspermutation(result)
#summary(resultqtl)
#mqmplot.permutations(result)
#write.csv(summary(homescan1),file="/users/rcalbert/documents/cichlids/QTL_LfxTrc/RAD_mapping/color_QTL/MQM/summary_mqm_flankE.csv")
#write.csv(summary(resultqtl),file="/users/rcalbert/documents/cichlids/QTL_LfxTrc/RAD_mapping/color_QTL/MQM/summary_mqm_flankE.csv")