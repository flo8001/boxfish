##### FINAL SCRIPT by TINA for BOXFISH ANALYSES ######

library(ggplot2) 
library(ape)
library(reshape2)
library(phytools)
library(geiger)
library(auteur)
library(nlme)

#####################################################
############### Import and Fix Data #################
#####################################################

# import data
tree <- read.tree("~/Dropbox/Work/boxfish/BOXFISH MANUSCRIPT WORKING DATA/comparative_boxfish_phylogeny.phy")   # this is a pruned tree

PCmeans <- read.csv("~/Dropbox/Work/boxfish/BOXFISH MANUSCRIPT WORKING DATA/species mean data/PCscores_means.csv", row.names=1)

FEAvPC <- read.csv("~/Dropbox/Work/boxfish/BOXFISH MANUSCRIPT WORKING DATA/species mean data/comparativeFEAdata.csv", row.names=1)

CFDvPC <- read.csv("~/Dropbox/Work/boxfish/BOXFISH MANUSCRIPT WORKING DATA/species mean data/comparative_cfd_data.csv", row.names=1)

# add genus column to PCmeans
genusPC <- c("Acanthostracion", "Acanthostracion", "Acanthostracion", "Acanthostracion", "Anoplocapros", "Anoplocapros", "Anoplocapros", "Aracana", "Aracana", "Capropygia", "Kentrocapros", "Kentrocapros", "Kentrocapros", "Lactophrys", "Lactophrys", "Lactophrys", "Lactoria", "Lactoria", "Tetrosomus", "Lactoria", "Ostracion", "Ostracion", "Ostracion", "Ostracion", "Ostracion", "Ostracion", "Ostracion", "Ostracion", "Polyplacapros", "Tetrosomus", "Tetrosomus", "Tetrosomus")
newPC <- PCmeans
newPC$genus <- genusPC

genusFEA <- c("Tetrosomus", "Tetrosomus", "Tetrosomus", "Lactoria", "Lactoria", "Ostracion", "Ostracion", "Ostracion", "Ostracion", "Ostracion", "Ostracion", "Ostracion", "Acanthostracion", "Acanthostracion", "Lactophrys", "Lactophrys", "Lactophrys", "Anoplocapros", "Anoplocapros", "Anoplocapros", "Aracana", "Aracana", "Capropygia", "Kentrocapros", "Kentrocapros")
newFEA <- FEAvPC
newFEA$genus <- genusFEA


genusCFD <- c("Acanthostracion", "Acanthostracion", "Anoplocapros", "Anoplocapros", "Anoplocapros", "Aracana", "Aracana", "Capropygia", "Kentrocapros", "Lactophrys", "Lactophrys", "Lactophrys", "Lactoria", "Lactoria", "Tetrosomus", "Ostracion", "Ostracion", "Ostracion", "Ostracion", "Ostracion", "Ostracion", "Ostracion", "Tetrosomus", "Tetrosomus")
newCFD <- CFDvPC
newCFD$genus <- genusCFD

# correct with new names
rownames(newPC)[19] <- "Tetrosomus_fornasini"
rownames(newFEA)[3] <- "Tetrosomus_fornasini"
rownames(newCFD)[15] <- "Tetrosomus_fornasini"

# have to add this step otherwise treedata converts the data.frame to a list
newPC$species <- rownames(newPC)
newFEA$species <- rownames(newFEA)
newCFD$species <- rownames(newCFD)

# correct tree name
tree$tip.label[4] <- "Tetrosomus_fornasini"

#####################################################
################# Phylomorphospace ##################
#####################################################

# check names
comparison <- name.check(tree, newPC)

# drop the rows
for(i in 1:length(comparison$Data.not.tree)){
	newPC <- newPC[newPC$species != comparison$Data.not.tree[i],]
}

####### get ancestral states ########

# get PC values in same order as phylo tree data
full.data <- data.frame(row.names=newPC$species, genus=newPC$genus, pc1=newPC$PC_1, pc2=newPC$PC_2, pc3=newPC$PC_3)

# match data to tree to make sure in same order and tips dropped, etc.
# TREE DATA WAS NOT USED TO DELETE ROWS OF DATA BECAUSE IF THE FUNCTION DOES, IT RETURNS A LIST. 
td <- treedata(tree, full.data)

pc.scores <- td$data
box.tree <- td$phy

#sort PC data
pc.scores <- pc.scores[box.tree$tip.label,]

## get Anc States
anc.pc1 <- getAncStates(pc.scores$pc1, box.tree)
anc.pc2 <- getAncStates(pc.scores$pc2, box.tree)
anc.pc3 <- getAncStates(pc.scores$pc3, box.tree)

## put together by node number/name
anc.nodes <- data.frame(cbind(anc.pc1, anc.pc2, anc.pc3))


## put together by node number/name
anc.nodes <- data.frame(cbind(anc.pc1, anc.pc2, anc.pc3))

## how to assign shape to genus
genus.identifier <- c("Acanthostracion"=5, "Aracana"=16, "Capropygia"=15, "Kentrocapros"=17, "Lactophrys"=2, "Lactoria"=3, "Ostracion"=1, "Tetrosomus"=6, "Anoplocapros"=18)

# create empty sets
empty <- matrix(data=NA, nrow=length(box.tree$edge[,1]), ncol=1)
lines <- cbind(box.tree$edge, empty, empty, empty, empty, empty, empty)
## order are the edges for col 1 and 2. 3 is daughterX, 4 dY, 5 dZ, 5 nodeX, 6 nY, 7 nZ

#input into line matrix to create lines connecting nodes -- basically reorder data into columns of the same length fro the geom_segment() layer. each column represents a part of the coordinate (x1, x2, x3, y1, y2, y3)
for (i in 1:length(lines[,1])){
	lines[i,3] <- anc.nodes[toString(lines[i,1]), "anc.pc1"]
	lines[i,4] <- anc.nodes[toString(lines[i,1]), "anc.pc2"]
	lines[i,5] <- anc.nodes[toString(lines[i,1]), "anc.pc3"]
				
	if (sum(lines[i,2]==c(1:length(pc.scores$pc1))) > 0){
		lines[i,6] <- pc.scores$pc1[lines[i,2]]
		lines[i,7] <- pc.scores$pc2[lines[i,2]]
		lines[i,8] <- pc.scores$pc3[lines[i,2]]
	}
		
	else {
		lines[i,6] <- anc.nodes[toString(lines[i,2]), "anc.pc1"]
		lines[i,7] <- anc.nodes[toString(lines[i,2]), "anc.pc2"]
		lines[i,8] <- anc.nodes[toString(lines[i,2]), "anc.pc3"]
	}	
}


ggplot() + 
geom_point(data=anc.nodes, aes(x=anc.pc1, y=anc.pc2), colour="azure4") +
geom_segment(aes(x=lines[,3], y=lines[,4], xend=lines[,6], yend=lines[,7]), colour="azure4") +
geom_point(data=pc.scores, aes(x=pc1, y=pc2, shape=factor(genus)), size=3) + 
scale_shape_manual(name="Genera", values=genus.identifier) +
xlab("PC1") + ylab("PC2") 

ggplot() + 
geom_point(data=anc.nodes, aes(x=anc.pc1, y=anc.pc3), colour="azure4") +
geom_segment(aes(x=lines[,3], y=lines[,5], xend=lines[,6], yend=lines[,8]), colour="azure4") +
geom_point(data=pc.scores, aes(x=pc1, y=pc3, shape=factor(genus)), size=3) + 
scale_shape_manual(name="Genera", values=genus.identifier) +
xlab("PC1") + ylab("PC3") 

ggplot() + 
geom_point(data=anc.nodes, aes(x=anc.pc2, y=anc.pc3), colour="azure4") +
geom_segment(aes(x=lines[,4], y=lines[,5], xend=lines[,7], yend=lines[,8]), colour="azure4") +
geom_point(data=pc.scores, aes(x=pc2, y=pc3, shape=factor(genus)), size=3) + 
scale_shape_manual(name="Genera", values=genus.identifier) +
xlab("PC2") + ylab("PC3") 

#####################################################
####################### PGLS ########################
#####################################################

# drop tips
align <- treedata(tree, newCFD)

CFDdata <- align$data
CFDtree <- align$phy

# sort data
CFDdata <- CFDdata[CFDtree$tip.label,]

# use Brownaian model with trait data/tree
bm.boxfish <- corBrownian(phy=CFDtree)

# perform each pgls
pc1.cd_p0 <- gls(cd_p0~PC_1, correlation=bm.boxfish, data=CFDdata)
pc1.cd_p20 <- gls(cd_p20~PC_1, correlation=bm.boxfish, data=CFDdata) 
pc1.cl_p20 <- gls(cl_p20~PC_1, correlation=bm.boxfish, data=CFDdata)
pc1.lift.drag_p20 <- gls(liftdrag_p20~PC_1, correlation=bm.boxfish, data=CFDdata)

pc2.cd_p0 <- gls(cd_p0~PC_2, correlation=bm.boxfish, data=CFDdata)
pc2.cd_p20 <- gls(cd_p20~PC_2, correlation=bm.boxfish, data=CFDdata) 
pc2.cl_p20 <- gls(cl_p20~PC_2, correlation=bm.boxfish, data=CFDdata)
pc2.lift.drag_p20 <- gls(liftdrag_p20~PC_2, correlation=bm.boxfish, data=CFDdata)

pc3.cd_p0 <- gls(cd_p0~PC_3, correlation=bm.boxfish, data=CFDdata)
pc3.cd_p20 <- gls(cd_p20~PC_3, correlation=bm.boxfish, data=CFDdata) 
pc3.cl_p20 <- gls(cl_p20~PC_3, correlation=bm.boxfish, data=CFDdata)
pc3.lift.drag_p20 <- gls(liftdrag_p20~PC_3, correlation=bm.boxfish, data=CFDdata)

# multivariate pgls performed
pc1.pc2.cd_p0 <- gls(cd_p0~PC_1*PC_2, correlation=bm.boxfish, data=CFDdata)
pc1.pc2.cd_p20 <- gls(cd_p20~PC_1*PC_2, correlation=bm.boxfish, data=CFDdata) 
pc1.pc2.cl_p20 <- gls(cl_p20~PC_1*PC_2, correlation=bm.boxfish, data=CFDdata)
pc1.pc2.lift.drag_p20 <- gls(liftdrag_p20~PC_1*PC_2, correlation=bm.boxfish, data=CFDdata)

pc1.pc3.cd_p0 <- gls(cd_p0~PC_1*PC_3, correlation=bm.boxfish, data=CFDdata)
pc1.pc3.cd_p20 <- gls(cd_p20~PC_1*PC_3, correlation=bm.boxfish, data=CFDdata) 
pc1.pc3.cl_p20 <- gls(cl_p20~PC_1*PC_3, correlation=bm.boxfish, data=CFDdata)
pc1.pc3.lift.drag_p20 <- gls(liftdrag_p20~PC_1*PC_3, correlation=bm.boxfish, data=CFDdata)

pc2.pc3.cd_p0 <- gls(cd_p0~PC_2*PC_3, correlation=bm.boxfish, data=CFDdata)
pc2.pc3.cd_p20 <- gls(cd_p20~PC_2*PC_3, correlation=bm.boxfish, data=CFDdata) 
pc2.pc3.cl_p20 <- gls(cl_p20~PC_2*PC_3, correlation=bm.boxfish, data=CFDdata)
pc2.pc3.lift.drag_p20 <- gls(liftdrag_p20~PC_2*PC_3, correlation=bm.boxfish, data=CFDdata)

pc1.pc2.pc3.cd_p0 <- gls(cd_p0~PC_1*PC_2*PC_3, correlation=bm.boxfish, data=CFDdata)
pc1.pc2.pc3.cd_p20 <- gls(cd_p20~PC_1*PC_2*PC_3, correlation=bm.boxfish, data=CFDdata) 
pc1.pc2.pc3.cl_p20 <- gls(cl_p20~PC_1*PC_2*PC_3, correlation=bm.boxfish, data=CFDdata)
pc1.pc2.pc3.lift.drag_p20 <- gls(liftdrag_p20~PC_1*PC_2*PC_3, correlation=bm.boxfish, data=CFDdata)

# create new list to draw from
SuperPGLS <- list(pc1.cd_p0=summary(pc1.cd_p0), pc1.cd_p20=summary(pc1.cd_p20), pc1.cl_p20=summary(pc1.cl_p20), pc1.lift.drag_p20=summary(pc1.lift.drag_p20), pc2.cd_p0=summary(pc2.cd_p0), pc2.cd_p20=summary(pc2.cd_p20), pc2.cl_p20=summary(pc2.cl_p20), pc2.lift.drag_p20=summary(pc2.lift.drag_p20), pc3.cd_p0=summary(pc3.cd_p0), pc3.cd_p20=summary(pc3.cd_p20), pc3.cl_p20=summary(pc3.cl_p20), pc3.lift.drag_p20=summary(pc3.lift.drag_p20))

rows <- c("pc1.cd_p0", "pc1.cd_p20", "pc1.cl_p20", "pc1.lift.drag_p20", "pc2.cd_p0", "pc2.cd_p20", "pc2.cl_p20", "pc2.lift.drag_p20", "pc3.cd_p0", "pc3.cd_p20", "pc3.cl_p20", "pc3.lift.drag_p20")

cols <- c("pValInt", "corr", "slope", "intercept")

# create empty matrix to summarize data -- placeholder
ph <- matrix(data=NA, nrow=12, ncol=4)

# assign row and col names
rownames(ph) <- rows
colnames(ph) <- cols

for (i in 1:length(SuperPGLS)){
	ph[i,1] <- SuperPGLS[[i]]$tTable[2,4]
	ph[i,2] <- SuperPGLS[[i]]$corBeta[2,1]
	ph[i,3] <- SuperPGLS[[i]]$tTable[2,1]
	ph[i,4] <- SuperPGLS[[i]]$tTable[1,1]
}

# write variable to table
write.csv(ph, file="~/Desktop/PGLSdata.csv")

# create another list
multiPGLS <- list(pc1.pc2.cd_p0=summary(pc1.pc2.cd_p0), pc1.pc2.cd_p20=summary(pc1.pc2.cd_p20), pc1.pc2.cl_p20=summary(pc1.pc2.cl_p20), pc1.pc2.lift.drag_p20=summary(pc1.pc2.lift.drag_p20), pc1.pc3.cd_p0=summary(pc1.pc3.cd_p0), pc1.pc3.cd_p20=summary(pc1.pc3.cd_p20), pc1.pc3.cl_p20=summary(pc1.pc3.cl_p20), pc1.pc3.lift.drag_p20=summary(pc1.pc3.lift.drag_p20), pc2.pc3.cd_p0=summary(pc2.pc3.cd_p0), pc2.pc3.cd_p20=summary(pc2.pc3.cd_p20), pc2.pc3.cl_p20=summary(pc2.pc3.cl_p20), pc2.pc3.lift.drag_p20=summary(pc2.pc3.lift.drag_p20), pc1.pc2.pc3.cd_p0=summary(pc1.pc2.pc3.cd_p0), pc1.pc2.pc3.cd_p20=summary(pc1.pc2.pc3.cd_p20), pc1.pc2.pc3.cl_p20=summary(pc1.pc2.pc3.cl_p20), pc1.pc2.pc3.lift.drag_p20=summary(pc1.pc2.pc3.lift.drag_p20))

# create another empty dataset, etc as before
md <- matrix(data=NA, nrow=16, ncol=4)

cols2 <- c("logLik", "BIC", "AIC", "p-value")
rows2 <- c("pc1.pc2.cd_p0", "pc1.pc2.cd_p20", "pc1.pc2.cl_p20", "pc1.pc2.lift.drag_p20","pc1.pc3.cd_p0", "pc1.pc3.cd_p20", "pc1.pc3.cl_p20", "pc1.pc3.lift.drag_p20", "pc2.pc3.cd_p0", "pc2.pc3.cd_p20", "pc2.pc3.cl_p20", "pc2.pc3.lift.drag_p20", "pc1.pc2.pc3.cd_p0", "pc1.pc2.pc3.cd_p20", "pc1.pc2.pc3.cl_p20", "pc1.pc2.pc3.lift.drag_p20")

colnames(md) <- cols2
rownames(md) <- rows2

# forloop as before except modified with conditional to take into consideration multivariate analysis with pcs 1, 2 and 3 which will have a larger table holding the p-values
for (i in 1:length(multiPGLS)){
	md[i,1] <- multiPGLS[[i]]$logLik
	md[i,2] <- multiPGLS[[i]]$BIC
	md[i,3] <- multiPGLS[[i]]$AIC
	if (length(multiPGLS[[i]]$tTable)==16){
		md[i,4] <- multiPGLS[[i]]$tTable[4,4]
	}
	else {
		md[i,4] <- multiPGLS[[i]]$tTable[8,4]
	}
}


write.csv(md, file="~/Desktop/multiPGLSdata.csv")


# saved genus dot, etc
# pch=ifelse(CFDdata$genus=="Acanthostracion",5, ifelse(CFDdata$genus=="Aracana",16, ifelse(CFDdata$genus=="Capropygia",15, ifelse(CFDdata$genus=="Kentrocapros",17, ifelse(CFDdata$genus=="Lactophrys", 2, ifelse(CFDdata$genus=="Lactoria", 3, ifelse(CFDdata$genus=="Ostracion", 1, ifelse(CFDdata$genus=="Tetrosomus", 6, 18))))))))

## PLOT using plot()... set abline(summary(blah))
# mar sets margin sizes (bottom, left, top, right)
pdf("PanelPlot.pdf", width=6.3, height=3.14)
par(mfrow=c(3,4), mar=c(4,4,2,2))
plot(x=CFDdata$PC_1, y=CFDdata$cd_p0, xlab="PC1", ylab="Drag Coefficient at 0 Degrees", pch=16)
plot(x=CFDdata$PC_1, y=CFDdata$cd_p20, xlab="PC1", ylab="Drag Coefficient at 20 Degrees", pch=16)
plot(x=CFDdata$PC_1, y=CFDdata$cl_p20, xlab="PC1", ylab="Lift Coefficient at 20 Degrees", pch=16)
plot(x=CFDdata$PC_1, y=CFDdata$liftdrag_p20, xlab="PC1",  ylab="Lift/Drag at 20 Degrees", pch=16)
plot(x=CFDdata$PC_2, y=CFDdata$cd_p0, xlab="PC2", ylab="Drag Coefficient at 0 Degrees", pch=16)
plot(x=CFDdata$PC_2, y=CFDdata$cd_p20, xlab="PC2", ylab="Drag Coefficient at 20 Degrees", pch=16)
plot(x=CFDdata$PC_2, y=CFDdata$cl_p20, xlab="PC2", ylab="Lift Coefficient at 20 Degrees", pch=16)
abline(summary(pc2.cl_p20))
plot(x=CFDdata$PC_2, y=CFDdata$liftdrag_p20, xlab="PC2", ylab="Lift/Drag at 20 Degrees", pch=16)
abline(summary(pc2.lift.drag_p20))
plot(x=CFDdata$PC_3, y=CFDdata$cd_p0, xlab="PC3", ylab="Drag Coefficient at 0 Degrees", pch=16)
plot(x=CFDdata$PC_3, y=CFDdata$cd_p20, xlab="PC3", ylab="Drag Coefficient at 20 Degrees", pch=16)
plot(x=CFDdata$PC_3, y=CFDdata$cl_p20, xlab="PC3", ylab="Lift Coefficient at 20 Degrees", pch=16)
plot(x=CFDdata$PC_3, y=CFDdata$liftdrag_p20, xlab="PC3", ylab="Lift/Drag at 20 Degrees", pch=16)
dev.off()
#####################################################
###################### Auteur #######################
#####################################################

# uses CFDdata for CFD comparisons and the normal PC dataset for the PCs and newFEA for FEA data

# create new matrix of pertinent FEAdata
fea <- cbind(dv=newFEA$dv_stress, ml=newFEA$ml_stress)

rownames(fea) <- rownames(newFEA)

tc <- treedata(tree, fea)

fea.data <- tc$data
fea.tree <- tc$phy

# sort them
fea.data <- fea.data[fea.tree$tip.label,]

#PC1
rjmcmc.bm(box.tree, pc.scores[,2], SE=0, ngen = 100000, sample.freq = 100, fileBase = "pc1_result")

load("~/BM.pc1_result.parameters/pc1_result.posteriorsamples.rda")

pdf('auteur_pc1.pdf')
shifts.plot(phy=box.tree, base.dir="BM.pc1_result.parameters", burnin=0.5, legend=TRUE, edge.width=2) 
dev.off()

#PC2
rjmcmc.bm(box.tree, pc.scores[,3], SE=0, ngen = 100000, sample.freq = 100, fileBase = "pc2_result")

load("~/BM.pc2_result.parameters/pc2_result.posteriorsamples.rda")

pdf('auteur_pc2.pdf')
shifts.plot(phy=box.tree, base.dir="BM.pc2_result.parameters", burnin=0.5, legend=TRUE, edge.width=2) 
dev.off()

#PC3
rjmcmc.bm(box.tree, pc.scores[,4], SE=0, ngen = 100000, sample.freq = 100, fileBase = "pc3_result")

load("~/BM.pc3_result.parameters/pc3_result.posteriorsamples.rda")

pdf('auteur_pc3.pdf')
shifts.plot(phy=box.tree, base.dir="BM.pc3_result.parameters", burnin=0.5, legend=TRUE, edge.width=2) 
dev.off()

#cd_p0
rjmcmc.bm(CFDtree, CFDdata[,4], SE=0, ngen = 100000, sample.freq = 100, fileBase = "cd_p0_result")

load("~/BM.cd_p0_result.parameters/cd_p0_result.posteriorsamples.rda")

pdf('auteur_cd_p0.pdf')
shifts.plot(phy=CFDtree, base.dir="BM.cd_p0_result.parameters", burnin=0.5, legend=TRUE, edge.width=2) 
dev.off()

#cd_p20
rjmcmc.bm(CFDtree, CFDdata[,5], SE=0, ngen = 100000, sample.freq = 100, fileBase = "cd_p20_result")

load("~/BM.cd_p20_result.parameters/cd_p20_result.posteriorsamples.rda")

pdf('auteur_cd_p20.pdf')
shifts.plot(phy=CFDtree, base.dir="BM.cd_p20_result.parameters", burnin=0.5, legend=TRUE, edge.width=2) 
dev.off()

#cl_p20
rjmcmc.bm(CFDtree, CFDdata[,6], SE=0, ngen = 100000, sample.freq = 100, fileBase = "cl_p20_result")

load("~/BM.cl_p20_result.parameters/cl_p20_result.posteriorsamples.rda")

pdf('auteur_cl_p20.pdf')
shifts.plot(phy=CFDtree, base.dir="BM.cl_p20_result.parameters", burnin=0.5, legend=TRUE, edge.width=2) 
dev.off()

#liftdrag_p20
rjmcmc.bm(CFDtree, CFDdata[,7], SE=0, ngen = 100000, sample.freq = 100, fileBase = "liftdrag_p20_result")

load("~/BM.liftdrag_p20_result.parameters/liftdrag_p20_result.posteriorsamples.rda")

pdf('auteur_liftdrag_p20.pdf')
shifts.plot(phy=CFDtree, base.dir="BM.liftdrag_p20_result.parameters", burnin=0.5, legend=TRUE, edge.width=2) 
dev.off()

#dv stress
rjmcmc.bm(fea.tree, fea.data[,1], SE=0, ngen = 100000, sample.freq = 100, fileBase = "dv_stress_result")

load("~/BM.dv_stress_result.parameters/dv_stress_result.posteriorsamples.rda")

pdf('auteur_dv_stress.pdf')
shifts.plot(phy=fea.tree, base.dir="BM.dv_stress_result.parameters", burnin=0.5, legend=TRUE, edge.width=2) 
dev.off()

#ml stress
rjmcmc.bm(fea.tree, fea.data[,2], SE=0, ngen = 100000, sample.freq = 100, fileBase = "ml_stress_result")

load("~/BM.ml_stress_result.parameters/ml_stress_result.posteriorsamples.rda")

pdf('auteur_ml_stress.pdf')
shifts.plot(phy=fea.tree, base.dir="BM.ml_stress_result.parameters", burnin=0.5, legend=TRUE, edge.width=2) 
dev.off()

#####################################################
################### Model Fitting ###################
#####################################################

# create new matrix of pertinent CFDdata
cfd <- cbind(cd.p0=CFDdata$cd_p0, cd.p20=CFDdata$cd_p20, cl.p20=CFDdata$cl_p20, lift.drag.p20=CFDdata$liftdrag_p20)

rownames(cfd) <- rownames(CFDdata)

#create new matrix of pertinent pc data
pc <- cbind(pc1=pc.scores$pc1, pc2=pc.scores$pc2, pc3=pc.scores$pc3)

rownames(pc) <- rownames(pc.scores)

# model fitting
bmP <- fitContinuous(box.tree, pc, model = "BM")
ouP <- fitContinuous(box.tree, pc, model = "OU")
ebP <- fitContinuous(box.tree, pc, model = "EB", bounds=list(a=c(-0.2,0.2)))

bm <- fitContinuous(CFDtree, cfd, model = "BM")
ou <- fitContinuous(CFDtree, cfd, model = "OU")
eb <- fitContinuous(CFDtree, cfd, model = "EB", bounds=list(a = c(-0.2,0.2)))

bmF <- fitContinuous(fea.tree, fea.data, model = "BM")
ouF <- fitContinuous(fea.tree, fea.data, model = "OU")
ebF <- fitContinuous(fea.tree, fea.data, model = "EB", bounds=list(a=c(-0.2,0.2)))

# source AK weights script, assign name models
source('~/Dropbox/R functions utilities/AKWeights.R', chdir = TRUE)

args(AkaikeWeights)

models <- c("BM", "OU", "ACDC")

##

pc1.aic <- c(bmP$pc1$aicc, ouP$pc1$aicc, ebP$pc1$aicc)
names(pc1.aic) <- models

pc1.aic.res <- AkaikeWeights(pc1.aic)
pc1.aic.res

##

pc2.aic <- c(bmP$pc2$aicc, ouP$pc2$aicc, ebP$pc2$aicc)
names(pc2.aic) <- models

pc2.aic.res <- AkaikeWeights(pc2.aic)
pc2.aic.res

##

pc3.aic <- c(bmP$pc3$aicc, ouP$pc3$aicc, ebP$pc3$aicc)
names(pc3.aic) <- models

pc3.aic.res <- AkaikeWeights(pc3.aic)
pc3.aic.res

##

cd.p0.aic <- c(bm$cd.p0$aicc, ou$cd.p0$aicc, eb$cd.p0$aicc)
names(cd.p0.aic) <- models

cd.p0.aic.res <- AkaikeWeights(cd.p0.aic)
cd.p0.aic.res

##

cl.p20.aic <- c(bm$cl.p20$aicc, ou$cl.p20$aicc, eb$cl.p20$aicc)
names(cl.p20.aic) <- models

cl.p20.aic.res <- AkaikeWeights(cl.p20.aic)
cl.p20.aic.res

##

cd.p20.aic <- c(bm$cd.p20$aicc, ou$cd.p20$aicc, eb$cd.p20$aicc)
names(cd.p20.aic) <- models

cd.p20.aic.res <- AkaikeWeights(cd.p20.aic)
cd.p20.aic.res

##

lift.drag.p20.aic <- c(bm$lift.drag.p20$aicc, ou$lift.drag.p20$aicc, eb$lift.drag.p20$aicc)
names(lift.drag.p20.aic) <- models

lift.drag.p20.aic.res <- AkaikeWeights(lift.drag.p20.aic)
lift.drag.p20.aic.res

##

dv.stress.aic <- c(bmF$dv$aicc, ouF$dv$aicc, ebF$dv$aicc)
names(dv.stress.aic) <- models

dv.stress.aic.res <- AkaikeWeights(dv.stress.aic)
dv.stress.aic.res

##

ml.stress.aic <- c(bmF$ml$aicc, ouF$ml$aicc, ebF$ml$aicc)
names(ml.stress.aic) <- models

ml.stress.aic.res <- AkaikeWeights(ml.stress.aic)
ml.stress.aic.res

##

# make new super matrix for all results
modelfitting <- data.frame(matrix(data=NA, ncol=5, nrow=23))

mcols <- c("tested", "model", "AIC", "dAIC", "Weights")

colnames(modelfitting) <- mcols

#put all above data in modelfitting
items <- list(pc1.aic.res, pc2.aic.res, pc3.aic.res, cd.p0.aic.res, cd.p20.aic.res, cl.p20.aic.res, lift.drag.p20.aic.res, dv.stress.aic.res, ml.stress.aic.res)
item.names <- c("pc1", "pc2", "pc3", "cd.p0", "cd.p20", "cl.p20", "lift.drag.p20", "dv stress", "ml stress")

for (i in 1:9){
	modelfitting[i*3, 2] <- rownames(items[[i]][3,])
	modelfitting[i*3-1, 2] <- rownames(items[[i]][2,])
	modelfitting[i*3-2, 2] <- rownames(items[[i]][1,])
	modelfitting[i*3, 1] <- item.names[i]
	modelfitting[i*3-1, 1] <- item.names[i]
	modelfitting[i*3-2, 1] <- item.names[i]
	modelfitting[i*3-2, 3] <- items[[i]][1,1]
	modelfitting[i*3-2, 4] <- items[[i]][1,2]
	modelfitting[i*3-2, 5] <- items[[i]][1,3]
	modelfitting[i*3-1, 3] <- items[[i]][2,1]
	modelfitting[i*3-1, 4] <- items[[i]][2,2]
	modelfitting[i*3-1, 5] <- items[[i]][2,3]
	modelfitting[i*3, 3] <- items[[i]][3,1]
	modelfitting[i*3, 4] <- items[[i]][3,2]
	modelfitting[i*3, 5] <- items[[i]][3,3]
}

write.csv(modelfitting, file="~/Desktop/modelfittingres.csv")

#####################################################
################### Contour Plots ###################
#####################################################

#####################################################
##################### Scree Plot ####################
#####################################################

#####################################################
################### Variances Plot ##################
#####################################################



