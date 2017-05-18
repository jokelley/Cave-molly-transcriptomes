#Load packages for all analyses
library("limma")
library("DESeq2")
library("edgeR")
library("splines")
library("RColorBrewer")
library("ggplot2")
library("splines")
library("ReactomePA")
library("clusterProfiler")
library("biomaRt")
library("DOSE")
library("Hmisc")
library("WGCNA")
library("org.Hs.eg.db")

################################################################################################
#
#           1. Differential expression analysis: EdgeR                        
#
#
################################################################################################
#Analyses performed with R studio 3.1.2
#Bioconductor version 3.0 [bioC installer 1.16.5]

#set working directory 
setwd("~/Tissue_MS_MolecEcol_Files")

################################################################################################
#
#          1A. Load and subset counts file
#
#
################################################################################################
#Load counts file generated from eXpress (tissue.txt)
data <- read.delim(gzfile("tissue.txt.gz"),header=T, row.names=1) 

#count number of transcripts and number of individuals
dim(data)

#Round data as version of edgeR used for analyses requires whole integers 
data <- round(data)

#trim data: Each row (i.e. transcript) must have greater than 2 counts per million (cpm) and at least 3 of the individuals must have counts data
keep <- rowSums(cpm(data)>2) >= 3     

#Keep data that meet the filtering criteria
datanew <- data[keep,]

#count number of transcripts left after trimming
dim(datanew)

################################################################################################
#
#          1B. Differential expression analysis using EdgeR                    
#          Gill, Liver and Brain were analyzed separately.
#
################################################################################################
###########################################GILL################################################
#subset counts for gill in nonsulfidic surface (control) and nonsulfidic cave (no light)
d <- subset(datanew, select = c(BG1, BG2, BG3, BG4,
                                LG1, LG2, LG4))

#Define group replicates
group <-c(rep("NonSulfurSurf",4),rep("NonSulfurCave",3))

#Count transcripts and subset of individuals
dim(d)

#Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

#Count number of transcripts
dim(d)

#Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

#Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

#Estimate common dispersion and tagwise dispersion 
dg1 <- estimateDisp(dg)

#Calculate differences in the means (i.e. counts) between two groups 
et <- exactTest(dg1)

#Extract the top differentially expressed "tags" for a pair of given groups (based on p-value)
detags <- rownames(topTags(et, n=20))

#Output counts for the top differentially expressed tags 
cpm(dg1)[detags,]

#summarize differentially expresses statistics as up [1], down [-1] or no change [0], note this was based on an adjusted p-value of 0.05 and Benjamin-Hochberg correction "BH"
summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))

#Extract up- (logFC > 0) and downregulated (logFC < 0) out separately based on FDR cut off of 0.05
out <- topTags(et, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.05 
LunaG <- out[keep,]
LunaG <- data.frame(LunaG)
LunaG.up = rownames(LunaG[LunaG$logFC > 0,])
LunaG.down = rownames(LunaG[LunaG$logFC < 0,])

#Subset counts for gill in nonsulfidic surface (control) and sulfidic cave (sulfur, no light)
d <- subset(datanew, select = c(BG1, BG2, BG3, BG4,
                                CG1, CG2, CG3, CG4))
#Group replicates
group <-c(rep("NonSulfurSurf",4),rep("SulfurCave",4))

#Count transcripts and subset of individuals
dim(d)

#Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

#Count number of transcripts
dim(d)

#Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

#Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

#Estimate common dispersion and tagwise dispersion 
dg1 <- estimateDisp(dg)

#Calculate differences in the means (i.e. counts) between two groups 
et <- exactTest(dg1)

#Extract the top differentially expressed "tags" for a pair of given groups (based on p-value)
detags <- rownames(topTags(et, n=20))

#Output counts for the top differentially expressed tags 
cpm(dg1)[detags,]

#summarize differentially expresses statistics as up [1], down [-1] or no change [0], note this was based on a corrected p-value of 0.05 and Benjamin-Hochberg correction "BH"
summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))

#Extract up- (logFC > 0) and downregulated (logFC < 0) out separately based on FDR cut off of 0.05
out <- topTags(et, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.05 
CuevaG <- out[keep,]
CuevaG <- data.frame(CuevaG)
CuevaG.up = rownames(CuevaG[CuevaG$logFC > 0,])
CuevaG.down = rownames(CuevaG[CuevaG$logFC < 0,])

#Subset counts for gill in nonsulfidic surface (control) and sulfidic surface (sulfur)
d <- subset(datanew, select = c(BG1, BG2, BG3, BG4,
                                PSO2G1, PSO2G2, PSO2G3))

#Group replicates
group <-c(rep("NonSulfurSurf",4),rep("SulfurSurface",3))

#Count transcripts and subset of individuals
dim(d)

#Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

#Count number of transcripts
dim(d)

#Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

#Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

#Estimate common dispersion and tagwise dispersion 
dg1 <- estimateDisp(dg)

#Calculate differences in the means (i.e. counts) between two groups 
et <- exactTest(dg1)

#Extract the top differentially expressed "tags" for a pair of given groups (based on p-value)
detags <- rownames(topTags(et, n=20))

#Output counts for the top differentially expressed tags 
cpm(dg1)[detags,]

#summarize differentially expresses statistics as up [1], down [-1] or no change [0], note this was based on a corrected p-value of 0.05 and Benjamin-Hochberg correction "BH"
summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))

#Extract up- (logFC > 0) and downregulated (logFC < 0) out separately based on FDR cut off of 0.05
out <- topTags(et, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.05
PSO2G <- out[keep,]
PSO2G <- data.frame(PSO2G)
PSO2G.up = rownames(PSO2G[PSO2G$logFC > 0,])
PSO2G.down = rownames(PSO2G[PSO2G$logFC < 0,])

####################Venn diagram gill######################### 
####################Upregulated Cave ######################### 
#Create unique variable, and construct venns 
universe.UP <- unique(c(LunaG.up, CuevaG.up))
GroupA <- universe.UP %in% LunaG.up
GroupB <- universe.UP %in% CuevaG.up
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

####################Upregulated Sulfur ######################### 
#Create unique variable, and construct venns 
universe.UP <- unique(c(PSO2G.up, CuevaG.up))
GroupA <- universe.UP %in% PSO2G.up
GroupB <- universe.UP %in% CuevaG.up
input.df <- data.frame(PSO=GroupA, Cueva=GroupB)
colnames(input.df)=c("S Surface","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

####################Downregulated Cave######################### 
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(LunaG.down, CuevaG.down))
GroupA <- universe.DOWN %in% LunaG.down
GroupB <- universe.DOWN %in% CuevaG.down
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

####################Downregulated Sulfur######################### 
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(CuevaG.down, PSO2G.down))
GroupA <- universe.DOWN %in% PSO2G.down
GroupB <- universe.DOWN %in% CuevaG.down
input.df <- data.frame(PSO=GroupA, Cueva=GroupB)
colnames(input.df)=c("S Surface","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

###########################################LIVER################################################
#subset counts for gill in nonsulfidic surface (control) and nonsulfidic cave (no light)
d <- subset(datanew, select = c(BL1, BL2, BL3, BL4,
                                LL1, LL2, LL3, LL4))

#Group replicates
group <-c(rep("NonSulfurSurf",4),rep("NonSulfurCave",4))

#Count transcripts and subset of individuals
dim(d)

#Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

#Count number of transcripts
dim(d)

#Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

#Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

#Estimate common dispersion and tagwise dispersion 
dg1 <- estimateDisp(dg)

#Calculate differences in the means (i.e. counts) between two groups 
et <- exactTest(dg1)

#Extract the top differentially expressed "tags" for a pair of given groups (based on p-value)
detags <- rownames(topTags(et, n=20))

#Output counts for the top differentially expressed tags 
cpm(dg1)[detags,]

#summarize differentially expresses statistics as up [1], down [-1] or no change [0], note this was based on a corrected p-value of 0.05 and Benjamin-Hochberg correction "BH"
summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))

#Extract up- (logFC > 0) and downregulated (logFC < 0) out separately based on FDR cut off of 0.05
out <- topTags(et, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.05 
LunaL <- out[keep,]
LunaL <- data.frame(LunaL)
LunaL.up = rownames(LunaL[LunaL$logFC > 0,])
LunaL.down = rownames(LunaL[LunaL$logFC < 0,])

#subset counts for gill in nonsulfidic surface (control) and nonsulfidic cave (sulfur,no light)
d <- subset(datanew, select = c(BL1, BL2, BL3, BL4,
                                CL1, CL2, CL3, CL4))

#Group replicates
group <-c(rep("NonSulfurSurf",4),rep("SulfurCave",4))

#Count transcripts and subset of individuals
dim(d)

#Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

#Count number of transcripts
dim(d)

#Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

#Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

#Estimate common dispersion and tagwise dispersion 
dg1 <- estimateDisp(dg)

#Calculate differences in the means (i.e. counts) between two groups 
et <- exactTest(dg1)

#Extract the top differentially expressed "tags" for a pair of given groups (based on p-value)
detags <- rownames(topTags(et, n=20))

#Output counts for the top differentially expressed tags 
cpm(dg1)[detags,]

#summarize differentially expresses statistics as up [1], down [-1] or no change [0], note this was based on a corrected p-value of 0.05 and Benjamin-Hochberg correction "BH"
summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))

#Extract up- (logFC > 0) and downregulated (logFC < 0) out separately based on FDR cut off of 0.05
out <- topTags(et, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.05 
CuevaL <- out[keep,]
CuevaL <- data.frame(CuevaL)
CuevaL.up = rownames(CuevaL[CuevaL$logFC > 0,])
CuevaL.down = rownames(CuevaL[CuevaL$logFC < 0,])

#subset counts for gill in nonsulfidic surface (control) and sulfidic surface (sulfur)
d <- subset(datanew, select = c(BL1, BL2, BL3, BL4,
                                PSO2L1, PSO2L2, PSO2L3, PSO2L4))
#Group replicates
group <-c(rep("NonSulfurSurf",4),rep("SulfurSurf",4))

#Count transcripts and subset of individuals
dim(d)

#Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

#Count number of transcripts
dim(d)

#Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

#Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

#Estimate common dispersion and tagwise dispersion 
dg1 <- estimateDisp(dg)

#Calculate differences in the means (i.e. counts) between two groups 
et <- exactTest(dg1)

#Extract the top differentially expressed "tags" for a pair of given groups (based on p-value)
detags <- rownames(topTags(et, n=20))

#Output counts for the top differentially expressed tags 
cpm(dg1)[detags,]

#summarize differentially expresses statistics as up [1], down [-1] or no change [0], note this was based on a corrected p-value of 0.05 and Benjamin-Hochberg correction "BH"
summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))

#Extract up- (logFC > 0) and downregulated (logFC < 0) out separately based on FDR cut off of 0.05
out <- topTags(et, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.05 
PSO2L <- out[keep,]
PSO2L <- data.frame(PSO2L)
PSO2L.up = rownames(PSO2L[PSO2L$logFC > 0,])
PSO2L.down = rownames(PSO2L[PSO2L$logFC < 0,])

####################Venn diagram liver######################### 
####################Upregulated Cave########################### 
#Create unique variable, and construct venns 
universe.UP <- unique(c(LunaL.up, CuevaL.up))
GroupA <- universe.UP %in% LunaL.up
GroupB <- universe.UP %in% CuevaL.up
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

####################Upregulated Sulfur######################### 
#Create unique variable, and construct venns 
universe.UP <- unique(c(PSO2L.up, CuevaL.up))
GroupA <- universe.UP %in% PSO2L.up
GroupB <- universe.UP %in% CuevaL.up
input.df <- data.frame(PSO=GroupA, Cueva=GroupB)
colnames(input.df)=c("S Surface","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

####################Downregulated Cave######################### 
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(LunaL.down, CuevaL.down))
GroupA <- universe.DOWN %in% LunaL.down
GroupB <- universe.DOWN %in% CuevaL.down
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

####################Downregulated Sulfur######################### 
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(PSO2L.down, CuevaL.down))
GroupA <- universe.DOWN %in% PSO2L.down
GroupB <- universe.DOWN %in% CuevaL.down
input.df <- data.frame(PSO=GroupA, Cueva=GroupB)
colnames(input.df)=c("S Surface","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

###########################################BRAIN################################################
#subset counts for gill in nonsulfidic surface (control) and nonsulfidic cave (no light)
d <- subset(datanew, select = c(BB1, BB2, BB3, BB4,
                                LB2, LB3))

#Group replicates
group <-c(rep("NonSulfurSurf",4),rep("NonSulfurCave",2))

#Count transcripts and subset of individuals
dim(d)

#Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

#Count number of transcripts
dim(d)

#Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

#Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

#Estimate common dispersion and tagwise dispersion 
dg1 <- estimateDisp(dg)

#Calculate differences in the means (i.e. counts) between two groups 
et <- exactTest(dg1)

#Extract the top differentially expressed "tags" for a pair of given groups (based on p-value)
detags <- rownames(topTags(et, n=20))

#Output counts for the top differentially expressed tags 
cpm(dg1)[detags,]

#summarize differentially expresses statistics as up [1], down [-1] or no change [0], note this was based on a corrected p-value of 0.05 and Benjamin-Hochberg correction "BH"
summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))

#Extract up- (logFC > 0) and downregulated (logFC < 0) out separately based on FDR cut off of 0.05
out <- topTags(et, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.05 
LunaB <- out[keep,]
LunaB <- data.frame(LunaB)
LunaB.up = rownames(LunaB[LunaB$logFC > 0,])
LunaB.down = rownames(LunaB[LunaB$logFC < 0,])

#subset counts for gill in nonsulfidic surface (control) and sulfidic cave (sulfur, no light)
d <- subset(datanew, select = c(BB1, BB2, BB3, BB4,
                                CB1, CB2, CB3))

#Group replicates
group <-c(rep("NonSulfurSurf",4),rep("SulfurCave",3))

#Count transcripts and subset of individuals
dim(d)

#Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

#Count number of transcripts
dim(d)

#Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

#Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

#Estimate common dispersion and tagwise dispersion 
dg1 <- estimateDisp(dg)

#Calculate differences in the means (i.e. counts) between two groups 
et <- exactTest(dg1)

#Extract the top differentially expressed "tags" for a pair of given groups (based on p-value)
detags <- rownames(topTags(et, n=20))

#Output counts for the top differentially expressed tags 
cpm(dg1)[detags,]

#summarize differentially expresses statistics as up [1], down [-1] or no change [0], note this was based on a corrected p-value of 0.05 and Benjamin-Hochberg correction "BH"
summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))

#Extract up- (logFC > 0) and downregulated (logFC < 0) out separately based on FDR cut off of 0.05based on FDR cut off of 0.05
out <- topTags(et, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.05 
CuevaB <- out[keep,]
CuevaB <- data.frame(CuevaB)
CuevaB.up = rownames(CuevaB[CuevaB$logFC > 0,])
CuevaB.down = rownames(CuevaB[CuevaB$logFC < 0,])

#subset counts for gill in nonsulfidic surface (control) and sulfidic surface (sulfur)
d <- subset(datanew, select = c(BB1, BB2, BB3, BB4,
                                PSO2B1, PSO2B4))

#Group replicates
group <-c(rep("NonSulfurSurf",4),rep("SulfurSurface",2))

#Count transcripts and subset of individuals
dim(d)

#Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

#Count number of transcripts
dim(d)

#Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

#Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

#Estimate common dispersion and tagwise dispersion 
dg1 <- estimateDisp(dg)

#Calculate differences in the means (i.e. counts) between two groups 
et <- exactTest(dg1)

#Extract the top differentially expressed "tags" for a pair of given groups (based on p-value)
detags <- rownames(topTags(et, n=20))

#Output counts for the top differentially expressed tags 
cpm(dg1)[detags,]

#summarize differentially expresses statistics as up [1], down [-1] or no change [0], note this was based on a corrected p-value of 0.05 and Benjamin-Hochberg correction "BH"
summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))

#Extract up- (logFC > 0) and downregulated (logFC < 0) out separately based on FDR cut off of 0.05
out <- topTags(et, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.05 
PSO2B <- out[keep,]
PSO2B <- data.frame(PSO2B)
PSO2B.up = rownames(PSO2B[PSO2B$logFC > 0,])
PSO2B.down = rownames(PSO2B[PSO2B$logFC < 0,])

####################Venn diagram brain######################### 
####################Upregulated Cave########################### 
#Create unique variable, and construct venns 
universe.UP <- unique(c(LunaB.up, CuevaB.up))
GroupA <- universe.UP %in% LunaB.up
GroupB <- universe.UP %in% CuevaB.up
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

####################Upregulated Sulfur########################### 
#Create unique variable, and construct venns 
universe.UP <- unique(c(PSO2B.up, CuevaB.up))
GroupA <- universe.UP %in% PSO2B.up
GroupB <- universe.UP %in% CuevaB.up
input.df <- data.frame(PSO=GroupA, Cueva=GroupB)
colnames(input.df)=c("S Cave","S Surface")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

####################Downregulated Cave########################### 
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(LunaB.down, CuevaB.down))
GroupA <- universe.DOWN %in% LunaB.down
GroupB <- universe.DOWN %in% CuevaB.down
input.df <- data.frame(Cueva=GroupA, PSO=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

####################Downregulated Sulfur########################### 
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(PSO2B.down, CuevaB.down))
GroupA <- universe.DOWN %in% PSO2B.down
GroupB <- universe.DOWN %in% CuevaB.down
input.df <- data.frame(PSO=GroupA, Cueva=GroupB)
colnames(input.df)=c("S Cave","S Surface")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

################################################################################################
#
#          1C. MDS plots 
#
#
################################################################################################
###############################ALL INDIVIDUALS AND ALL ORGANS###################################
#subset counts of individuals used in all analyses (gill, brain, and liver)
d <- subset(datanew, select = c(BB1, BB2, BB3, BB4, BG1, BG2, BG3, BG4, BL1, BL2, BL3, BL4,
                                CB1, CB2, CB3, CG1, CG2, CG3, CG4, CL1, CL2, CL3, CL4,
                                LB2, LB3, LG1, LG2, LG4, LL1, LL2, LL3, LL4,
                                PSO2B1, PSO2B4, PSO2G1, PSO2G2, PSO2G3, PSO2L1, PSO2L2, PSO2L3, PSO2L4))

#Group replicates
group <-c(rep("Bonita",12),rep("Cueva",11),rep("Luna",9),rep("PSO",9))

#Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

#Count number of transcripts and individuals 
dim(d)

#Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

#Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

#Construct MDS plot on top 10,000 expressed transcripts 
mds <- plotMDS(dg, top=10000)

#Convert "names" into symbols and change color to represent the different organs
#Blue corresponds to brains, green corresponds to gills and red corresponds to liver
plot(mds, col=c(rep("blue",4),rep("green",4),rep("red",4),
                rep("blue",3),rep("green",4),rep("red",4),
                rep("blue",2),rep("green",3),rep("red",4),
                rep("blue",2),rep("green",3),rep("red",4)), pch=16)

#Extract residuals
organallmds <- (mds$cmdscale.out)

#Save residuals into .csv file to clean up MDS plots in Excel
write.csv(organallmds, file="all_organ_individual_mds.csv")

#######################################GILL#####################################################
#subset counts for all Gill tissue
d <- subset(datanew, select = c(BG1, BG2, BG3, BG4, 
                                CG1, CG2, CG3, CG4, 
                                LG1, LG2, LG4, 
                                PSO2G1, PSO2G2, PSO2G3))

#Group repliates
group <-c(rep("Bonita",4),rep("Cueva",4),rep("Luna",3),rep("PSO",3))

#Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

#Count number of transcripts and individuals
dim(d)

#Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

#Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

#Construct MDS plot on top 10,000 expressed transcripts 
mds <- plotMDS(dg, top=10000)

#Convert "names" into symbols and change color to represent the different habitats
#Blue corresponds to Nonsulfidic surface, green corresponds to Sulfidic cave 
#Red corresponds to Nonsulfidic cave and orange corresponds to Sulfidic surface
plot(mds, col=c(rep("blue",4),rep("green",4),rep("red",3),rep("orange",3)), pch=16)

#Extract residuals
gillmds <- (mds$cmdscale.out)

#Save residuals into .csv file to clean up MDS plots in Excel
write.csv(gillmds, file="gill_individual_mds.csv")

#######################################LIVER#####################################################
#subset counts for all Liver tissue
d <- subset(datanew, select = c(BL1, BL2, BL3, BL4,
                                CL1, CL2, CL3, CL4,
                                LL1, LL2, LL3, LL4,
                                PSO2L1, PSO2L2, PSO2L3, PSO2L4))

#Group replicates
group <-c(rep("Bonita",4),rep("Cueva",4),rep("Luna",4),rep("PSO",4))

#Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

#Count number of transcripts and individuals
dim(d)

#Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

#Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

#Construct MDS plot on top 10,000 expressed transcripts 
mds <- plotMDS(dg, top=10000)

#Convert "names" into symbols and change color to represent the different habitats
#Blue corresponds to Nonsulfidic surface, green corresponds to Sulfidic cave 
#Red corresponds to Nonsulfidic cave and orange corresponds to Sulfidic surface
plot(mds, col=c(rep("blue",4),rep("green",4),rep("red",4),rep("orange",4)), pch=16)

#Extract residuals
livermds <- (mds$cmdscale.out)

#Save residuals into .csv file to clean up MDS plots in Excel
write.csv(livermds, file="liver_individual_mds.csv")

#######################################BRAIN#####################################################
#Subset counts for all Brain tissue
d <- subset(datanew, select = c(BB1, BB2, BB3, BB4, 
                                CB1, CB2, CB3, 
                                LB2, LB3, 
                                PSO2B1, PSO2B4))

#Group replicates
group <-c(rep("Bonita",4),rep("Cueva",3),rep("Luna",2),rep("PSO",2))

#Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

#Count number of transcripts and individuals
dim(d)

#Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

#Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

#Construct MDS plot on top 10,000 expressed transcripts 
mds <- plotMDS(dg, top=10000)

#Convert "names" into symbols and change color to represent the different habitats
#Blue corresponds to Nonsulfidic surface, green corresponds to Sulfidic cave 
#Red corresponds to Nonsulfidic cave and orange corresponds to Sulfidic surface
plot(mds, col=c(rep("blue",4),rep("green",3),rep("red",2),rep("orange",2)), pch=16)

#Extract residuals
brainmds <- (mds$cmdscale.out)

#Save residuals into .csv file to clean up MDS plots in Excel
write.csv(brainmds, file="brain_individual_mds.csv")

################################################################################################
#
#          2. Genes shared between contrasting environmental conditions
#
#
################################################################################################
## Read in reference 
reference <- read.csv("TableS4.csv") #Table S4 from the manuscript which provides the reference 
# annotations from BLAST to SwissProt

##################################GILL##########################################################
##################################UP REGULATED##################################################
################################NONSULFUR-CAVE##################################################
#Convert transcript list to dataframe
GillNSCave.up <-  data.frame(LunaG.up)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(GillNSCave.up) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
GillNSCAVEmerged.up <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], GillNSCave.up, by = "Query.Sequence.ID")

#subset swissprot accessions 
GillNSCAVEmerged_select.up <- subset(GillNSCAVEmerged.up, select=c("Subject.sequence.ID"))

#remove duplicates
GillNSCAVEmerged_unique.up <- as.character(unique(unlist(GillNSCAVEmerged_select.up)))

###################################SULFUR-CAVE##################################################
#Convert transcript list to dataframe
GillSCave.up <-  data.frame(CuevaG.up)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(GillSCave.up) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
GillSCAVEmerged.up <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], GillSCave.up, by = "Query.Sequence.ID")

#subset swissprot accessions 
GillSCAVEmerged_select.up <- subset(GillSCAVEmerged.up, select=c("Subject.sequence.ID"))

#remove duplicates
GillSCAVEmerged_unique.up <- as.character(unique(unlist(GillSCAVEmerged_select.up)))

###################################SULFUR-SURFACE################################################
#Convert transcript list to dataframe
GillSSurface.up <-  data.frame(PSO2G.up)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(GillSSurface.up) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
GillSSURFACEmerged.up <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], GillSSurface.up, by = "Query.Sequence.ID")

#subset swissprot accessions 
GillSSURFACEmerged_select.up <- subset(GillSSURFACEmerged.up, select=c("Subject.sequence.ID"))

#remove duplicates
GillSSURFACEmerged_unique.up <- as.character(unique(unlist(GillSSURFACEmerged_select.up)))

##################################DOWN REGULATED##################################################
###################################NONSULFUR-CAVE################################################
#Convert transcript list to dataframe
GillNSCave.down <-  data.frame(LunaG.down)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(GillNSCave.down) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
GillNSCAVEmerged.down <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], GillNSCave.down, by = "Query.Sequence.ID")

#subset swissprot accessions 
GillNSCAVEmerged_select.down <- subset(GillNSCAVEmerged.down, select=c("Subject.sequence.ID"))

#remove duplicates
GillNSCAVEmerged_unique.down <- as.character(unique(unlist(GillNSCAVEmerged_select.down)))

###################################SULFUR-CAVE###################################################
#Convert transcript list to dataframe
GillSCave.down <-  data.frame(CuevaG.down)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(GillSCave.down) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
GillSCAVEmerged.down <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], GillSCave.down, by = "Query.Sequence.ID")

#subset swissprot accessions 
GillSCAVEmerged_select.down <- subset(GillSCAVEmerged.down, select=c("Subject.sequence.ID"))

#remove duplicates
GillSCAVEmerged_unique.down <- as.character(unique(unlist(GillSCAVEmerged_select.down)))

###################################SULFUR-SURFACE################################################
#Convert transcript list to dataframe
GillSSurface.down <-  data.frame(PSO2G.down)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(GillSSurface.down) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
GillSSURFACEmerged.down <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], GillSSurface.down, by = "Query.Sequence.ID")

#subset swissprot accessions 
GillSSURFACEmerged_select.down <- subset(GillSSURFACEmerged.down, select=c("Subject.sequence.ID"))

#remove duplicates
GillSSURFACEmerged_unique.down <- as.character(unique(unlist(GillSSURFACEmerged_select.down)))

##################################LIVER#########################################################
##################################UP REGULATED##################################################
################################NONSULFUR-CAVE##################################################
#Convert transcript list to dataframe
LiverNSCave.up <-  data.frame(LunaL.up)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(LiverNSCave.up) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
LIVERNSCAVEmerged.up <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], LiverNSCave.up, by = "Query.Sequence.ID")

#subset swissprot accessions 
LIVERNSCAVEmerged_select.up <- subset(LIVERNSCAVEmerged.up, select=c("Subject.sequence.ID"))

#remove duplicates
LIVERNSCAVEmerged_unique.up <- as.character(unique(unlist(LIVERNSCAVEmerged_select.up)))

###################################SULFUR-CAVE##################################################
#Convert transcript list to dataframe
LiverSCave.up <-  data.frame(CuevaL.up)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(LiverSCave.up) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
LIVERSCAVEmerged.up <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], LiverSCave.up, by = "Query.Sequence.ID")

#subset swissprot accessions
LIVERSCAVEmerged_select.up <- subset(LIVERSCAVEmerged.up, select=c("Subject.sequence.ID"))

#remove duplicates
LIVERSCAVEmerged_unique.up <- as.character(unique(unlist(LIVERSCAVEmerged_select.up)))

###################################SULFUR-SURFACE##################################################
#Convert transcript list to dataframe
LiverSSurface.up <-  data.frame(PSO2L.up)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(LiverSSurface.up) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
LIVERSSURFACEmerged.up <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], LiverSSurface.up, by = "Query.Sequence.ID")

#subset swissprot accessions
LIVERSSURFACEmerged_select.up <- subset(LIVERSSURFACEmerged.up, select=c("Subject.sequence.ID"))

#remove duplicates
LIVERSSURFACEmerged_unique.up <- as.character(unique(unlist(LIVERSSURFACEmerged_select.up)))

##################################DOWN REGULATED##################################################
################################NONSULFUR-CAVE##################################################
#Convert transcript list to dataframe
LiverNSCave.down <-  data.frame(LunaL.down)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(LiverNSCave.down) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
LIVERNSCAVEmerged.down <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], LiverNSCave.down, by = "Query.Sequence.ID")

#subset swissprot accessions
LIVERNSCAVEmerged_select.down <- subset(LIVERNSCAVEmerged.down, select=c("Subject.sequence.ID"))

#remove duplicates
LIVERNSCAVEmerged_unique.down <- as.character(unique(unlist(LIVERNSCAVEmerged_select.down)))

###################################SULFUR-CAVE##################################################
#Convert transcript list to dataframe
LiverSCave.down <-  data.frame(CuevaL.down)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(LiverSCave.down) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
LIVERSCAVEmerged.down <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], LiverSCave.down, by = "Query.Sequence.ID")

#subset swissprot accessions
LIVERSCAVEmerged_select.down <- subset(LIVERSCAVEmerged.down, select=c("Subject.sequence.ID"))

#remove duplicates
LIVERSCAVEmerged_unique.down <- as.character(unique(unlist(LIVERSCAVEmerged_select.down)))

###################################SULFUR-SURFACE##################################################
#Convert transcript list to dataframe
LiverSSurface.down <-  data.frame(PSO2L.down)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(LiverSSurface.down) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
LIVERSSURFACEmerged.down <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], LiverSSurface.down, by = "Query.Sequence.ID")

#subset swissprot accessions
LIVERSSURFACEmerged_select.down <- subset(LIVERSSURFACEmerged.down, select=c("Subject.sequence.ID"))

#remove duplicates
LIVERSSURFACEmerged_unique.down <- as.character(unique(unlist(LIVERSSURFACEmerged_select.down)))

##################################BRAIN#########################################################
##################################UP REGULATED##################################################
################################NONSULFUR-CAVE##################################################
#Convert transcript list to dataframe
BrainNSCave.up <-  data.frame(LunaB.up)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(BrainNSCave.up) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
BrainNSCAVEmerged.up <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], BrainNSCave.up, by = "Query.Sequence.ID")

#subset swissprot accessions
BrainNSCAVEmerged_select.up <- subset(BrainNSCAVEmerged.up, select=c("Subject.sequence.ID"))

#remove duplicates
BrainNSCAVEmerged_unique.up <- as.character(unique(unlist(BrainNSCAVEmerged_select.up)))

################################NONSULFUR-CAVE##################################################
#Convert transcript list to dataframe
BrainSCave.up <-  data.frame(CuevaB.up)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(BrainSCave.up) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
BrainSCAVEmerged.up <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], BrainSCave.up, by = "Query.Sequence.ID")

#subset swissprot accessions
BrainSCAVEmerged_select.up <- subset(BrainSCAVEmerged.up, select=c("Subject.sequence.ID"))

#remove duplicates
BrainSCAVEmerged_unique.up <- as.character(unique(unlist(BrainSCAVEmerged_select.up)))

################################NONSULFUR-SURFACE##################################################
#Convert transcript list to dataframe
BrainSSurface.up <-  data.frame(PSO2B.up)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(BrainSSurface.up ) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
BrainSSURFACEmerged.up <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], BrainSSurface.up, by = "Query.Sequence.ID")

#subset swissprot accessions
BrainSSURFACEmerged_select.up <- subset(BrainSSURFACEmerged.up, select=c("Subject.sequence.ID"))

#remove duplicates
BrainSSURFACEmerged_unique.up <- as.character(unique(unlist(BrainSSURFACEmerged_select.up)))

##################################UP REGULATED##################################################
################################NONSULFUR-CAVE##################################################
#Convert transcript list to dataframe
BrainNSCave.down <-  data.frame(LunaB.down)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(BrainNSCave.down) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
BrainNSCAVEmerged.down <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], BrainNSCave.down, by = "Query.Sequence.ID")

#subset swissprot accessions
BrainNSCAVEmerged_select.down <- subset(BrainNSCAVEmerged.down, select=c("Subject.sequence.ID"))

#remove duplicates
BrainNSCAVEmerged_unique.down <- as.character(unique(unlist(BrainNSCAVEmerged_select.down)))

################################NONSULFUR-CAVE##################################################
#Convert transcript list to dataframe
BrainSCave.down <-  data.frame(CuevaB.down)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(BrainSCave.down) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
BrainSCAVEmerged.down <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], BrainSCave.down, by = "Query.Sequence.ID")

#subset swissprot accessions
BrainSCAVEmerged_select.down <- subset(BrainSCAVEmerged.down, select=c("Subject.sequence.ID"))

#remove duplicates
BrainSCAVEmerged_unique.down <- as.character(unique(unlist(BrainSCAVEmerged_select.down)))

################################NONSULFUR-SURFACE##################################################
#Convert transcript list to dataframe
BrainSSurface.down <-  data.frame(PSO2B.down)

#add column heading that matches the blast output (reference.csv)
col_headings <- c("Query.Sequence.ID")
names(BrainSSurface.down) <- col_headings

##Extract swissprot accessions (subject sequence IDs) from Blast output based on transcript IDs (query sequence IDs)
BrainSSURFACEmerged.down <- merge(reference[, c("Query.Sequence.ID", "Subject.sequence.ID")], BrainSSurface.down, by = "Query.Sequence.ID")

#subset swissprot accessions
BrainSSURFACEmerged_select.down <- subset(BrainSSURFACEmerged.down, select=c("Subject.sequence.ID"))

#remove duplicates
BrainSSURFACEmerged_unique.down <- as.character(unique(unlist(BrainSSURFACEmerged_select.down)))

###################################VENN DIAGRAMS################################################
############Venn Diagram Gill Caves Up-regulated##############
#Create unique variable, and construct venns 
universe.UP <- unique(c(GillNSCAVEmerged_unique.up, GillSCAVEmerged_unique.up))
GroupA <- universe.UP %in% GillNSCAVEmerged_unique.up
GroupB <- universe.UP %in% GillSCAVEmerged_unique.up
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

##############Venn Diagram Gill Sulfur Up-regulated############
#Create unique variable, and construct venns 
universe.UP <- unique(c(GillSSURFACEmerged_unique.up, GillSCAVEmerged_unique.up))
GroupA <- universe.UP %in% GillSSURFACEmerged_unique.up
GroupB <- universe.UP %in% GillSCAVEmerged_unique.up
input.df <- data.frame(PSO=GroupA, Cueva=GroupB)
colnames(input.df)=c("S surface","S cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

##############Venn Diagram Gill Caves Down-regulated#############
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(GillNSCAVEmerged_unique.down, GillSCAVEmerged_unique.down))
GroupA <- universe.DOWN %in% GillNSCAVEmerged_unique.down
GroupB <- universe.DOWN %in% GillSCAVEmerged_unique.down
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

##############Venn Diagram Gill Sulfur Down-regulated#############
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(GillSSURFACEmerged_unique.down, GillSCAVEmerged_unique.down))
GroupA <- universe.DOWN %in% GillSSURFACEmerged_unique.down
GroupB <- universe.DOWN %in% GillSCAVEmerged_unique.down
input.df <- data.frame(PSO=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

######################Venn Diagram Liver Caves Up-regulated########################
#Create unique variable, and construct venns 
universe.UP <- unique(c(LIVERNSCAVEmerged_unique.up, LIVERSCAVEmerged_unique.up))
GroupA <- universe.UP %in% LIVERNSCAVEmerged_unique.up
GroupB <- universe.UP %in% LIVERSCAVEmerged_unique.up
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#############Venn Diagram Liver Sulfur Up-regulated#################
#Create unique variable, and construct venns 
universe.UP <- unique(c(LIVERSSURFACEmerged_unique.up, LIVERSCAVEmerged_unique.up))
GroupA <- universe.UP %in% LIVERSSURFACEmerged_unique.up
GroupB <- universe.UP %in% LIVERSCAVEmerged_unique.up
input.df <- data.frame(PSO=GroupA, Cueva=GroupB)
colnames(input.df)=c("S surface","S cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

##############Venn Diagram Liver Caves Down-regulated##############
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(LIVERNSCAVEmerged_unique.down, LIVERSCAVEmerged_unique.down))
GroupA <- universe.DOWN %in% LIVERNSCAVEmerged_unique.down
GroupB <- universe.DOWN %in% LIVERSCAVEmerged_unique.down
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

######Venn Diagram Liver Sulfur Down-regulated#####################
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(LIVERSSURFACEmerged_unique.down, LIVERSCAVEmerged_unique.down))
GroupA <- universe.DOWN %in% LIVERSSURFACEmerged_unique.down
GroupB <- universe.DOWN %in% LIVERSCAVEmerged_unique.down
input.df <- data.frame(PSO=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

############Venn Diagram Brain Caves Up-regulated###############
#Create unique variable, and construct venns 
universe.UP <- unique(c(BrainNSCAVEmerged_unique.up, BrainSCAVEmerged_unique.up))
GroupA <- universe.UP %in% BrainNSCAVEmerged_unique.up
GroupB <- universe.UP %in% BrainSCAVEmerged_unique.up
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#############Venn Diagram Brain Sulfur Up-regulated###############
#Create unique variable, and construct venns 
universe.UP <- unique(c(BrainSSURFACEmerged_unique.up, BrainSCAVEmerged_unique.up))
GroupA <- universe.UP %in% BrainSSURFACEmerged_unique.up
GroupB <- universe.UP %in% BrainSCAVEmerged_unique.up
input.df <- data.frame(PSO=GroupA, Cueva=GroupB)
colnames(input.df)=c("S surface","S cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#############Venn Diagram Brain Caves Down-regulated##############
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(BrainNSCAVEmerged_unique.down, BrainSCAVEmerged_unique.down))
GroupA <- universe.DOWN %in% BrainNSCAVEmerged_unique.down
GroupB <- universe.DOWN %in% BrainSCAVEmerged_unique.down
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#############Venn Diagram Brain Sulfur Down-regulated##############
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(BrainSSURFACEmerged_unique.down, BrainSCAVEmerged_unique.down))
GroupA <- universe.DOWN %in% BrainSSURFACEmerged_unique.down
GroupB <- universe.DOWN %in% BrainSCAVEmerged_unique.down
input.df <- data.frame(PSO=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

################################################################################################
#
#          3. Physiological pathways (gene function) shared between contrasting environmental conditions
#
#
################################################################################################
###See README for information pertaining to input/ouput files for gene function analyses 

################################################################################################
#
#          3A. Generate input files for Blast2Go
#
#
################################################################################################

#########################GILL###############################
#############SULFIDIC SURFACE UPREGULATED###################
#Convert list of swissprot accesions into a dataframe
GillSSurf.up <- data.frame(GillSSURFACEmerged_unique.up)

#split swissprot accession string based on "|"
GillSSurf_split.up <- data.frame(do.call('rbind', strsplit(as.character(GillSSurf.up$GillSSURFACEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
GillSSurf_access.up <- GillSSurf_split.up$X2

#Write accession to table for Blast2Go input
write.table(GillSSurf_access.up, file = "Ssurfacegillup.txt", col.names = NA) 

#############SULFIDIC SURFACE DOWNREGULATED###################
#Convert list of swissprot accesions into a dataframe
GillSSurf.down <- data.frame(GillSSURFACEmerged_unique.down)

#split swissprot accession string based on "|"
GillSSurf_split.down <- data.frame(do.call('rbind', strsplit(as.character(GillSSurf.down$GillSSURFACEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
GillSSurf_access.down <- GillSSurf_split.down$X2

#Write accession to table for Blast2Go input
write.table(GillSSurf_access.down, file = "SSgilldown.txt", col.names = NA) 

#############SULFIDIC CAVE UPREGULATED###################
#Convert list of swissprot accesions into a dataframe
GillSCcave.up <- data.frame(GillSCAVEmerged_unique.up)

#split swissprot accession string based on "|"
GillSCcave_split.up <- data.frame(do.call('rbind', strsplit(as.character(GillSCcave.up$GillSCAVEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
GillSCcave_access.up <- GillSCcave_split.up$X2

#Write accession to table for Blast2Go input
write.table(GillSCcave_access.up, file = "SCgillup.txt", col.names = NA) 

#############SULFIDIC CAVE DOWNREGULATED#################
#Convert list of swissprot accesions into a dataframe
GillSCcave.down <- data.frame(GillSCAVEmerged_unique.down)

#split swissprot accession string based on "|"
GillSCcave_split.down <- data.frame(do.call('rbind', strsplit(as.character(GillSCcave.down$GillSCAVEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
GillSCcave_access.down <- GillSCcave_split.down$X2

#Write accession to table for Blast2Go input
write.table(GillSCcave_access.down, file = "Scavegilldown.txt", col.names = NA) 

#############NONSULFIDIC CAVE UPREGULATED################
#Convert list of swissprot accesions into a dataframe
GillNSCave.up <- data.frame(GillNSCAVEmerged_unique.up)

#split swissprot accession string based on "|"
GillNSCave_split.up <- data.frame(do.call('rbind', strsplit(as.character(GillNSCave.up$GillNSCAVEmerged_unique.up),'|',fixed=TRUE)))

#Extract out Human accession ID 
GillNSCave_access.up <- GillNSCave_split.up$X2

#Write accession to table for Blast2Go input
write.table(GillNSCave_access.up, file = "NSCgillup.txt", col.names = NA) 

#############NONSULFIDIC CAVE DOWNREGULATED##############
#Convert list of swissprot accesions into a dataframe
GillNSCave.down <- data.frame(GillNSCAVEmerged_unique.down)

#split swissprot accession string based on "|"
GillNSCave_split.down <- data.frame(do.call('rbind', strsplit(as.character(GillNSCave.down$GillNSCAVEmerged_unique.down),'|',fixed=TRUE)))

#Extract out Human accession ID 
GillNSCave_access.down <- GillNSCave_split.down$X2

#Write accession to table for Blast2Go input
write.table(GillNSCave_access.down, file = "NSCgilldown.txt", col.names = NA) 

#########################LIVER##############################
#############SULFIDIC SURFACE UPREGULATED###################
#Convert list of swissprot accesions into a dataframe
LiverSSurf.up <- data.frame(LIVERSSURFACEmerged_unique.up)

#split swissprot accession string based on "|"
LiverSSurf_split.up <- data.frame(do.call('rbind', strsplit(as.character(LiverSSurf.up$LIVERSSURFACEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
LiverSSurf_access.up <- LiverSSurf_split.up$X2

#Write accession to table for Blast2Go input
write.table(LiverSSurf_access.up, file = "SSliverup.txt", col.names = NA) 

#############SULFIDIC SURFACE DOWNREGULATED#################
#Convert list of swissprot accesions into a dataframe
LiverSSurf.down <- data.frame(LIVERSSURFACEmerged_unique.down)

#split swissprot accession string based on "|"
LiverSSurf_split.down <- data.frame(do.call('rbind', strsplit(as.character(LiverSSurf.down$LIVERSSURFACEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
LiverSSurf_access.down <- LiverSSurf_split.down$X2

#Write accession to table for Blast2Go input
write.table(LiverSSurf_access.down, file = "SSliverdown.txt", col.names = NA) 

#############SULFIDIC CAVE UPREGULATED###################
#Convert list of swissprot accesions into a dataframe
LiverSCave.up <- data.frame(LIVERSCAVEmerged_unique.up)

#split swissprot accession string based on "|"
LiverSCave_split.up <- data.frame(do.call('rbind', strsplit(as.character(LiverSCave.up$LIVERSCAVEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
LiverSCave_access.up <- LiverSCave_split.up$X2

#Write accession to table for Blast2Go input
write.table(LiverSCave_access.up, file = "SCliverup.txt", col.names = NA) 

#############SULFIDIC CAVE DOWNREGULATED#################
#Convert list of swissprot accesions into a dataframe
LiverSCave.down <- data.frame(LIVERSCAVEmerged_unique.down)

#split swissprot accession string based on "|"
LiverSCave_split.down <- data.frame(do.call('rbind', strsplit(as.character(LiverSCave.down$LIVERSCAVEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
LiverSCave_access.down <- LiverSCave_split.down$X2

#Write accession to table for Blast2Go input
write.table(LiverSCave_access.down, file = "SCliverdown.txt", col.names = NA) 

#############NONSULFIDIC CAVE UPREGULATED################
#Convert list of swissprot accesions into a dataframe
LiverNSCave.up <- data.frame(LIVERNSCAVEmerged_unique.up)

#split swissprot accession string based on "|"
LiverNSCave_split.up <- data.frame(do.call('rbind', strsplit(as.character(LiverNSCave.up$LIVERNSCAVEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
LiverNSCave_access.up <- LiverNSCave_split.up$X2

#Write accession to table for Blast2Go input
write.table(LiverNSCave_access.up, file = "NSCliverup.txt", col.names = NA) 

#############NONSULFIDIC CAVE DOWNREGULATED##############
#Convert list of swissprot accesions into a dataframe
LiverNSCave.down <- data.frame(LIVERNSCAVEmerged_unique.down)

#split swissprot accession string based on "|"
LiverNSCave_split.down <- data.frame(do.call('rbind', strsplit(as.character(LiverNSCave.down$LIVERNSCAVEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
LiverNSCave_access.down <- LiverNSCave_split.down$X2

#Write accession to table for Blast2Go input
write.table(LiverNSCave_access.down, file = "NSCliverdown.txt", col.names = NA) 

#########################BRAIN##############################
#############SULFIDIC SURFACE UPREGULATED###################
#Convert list of swissprot accesions into a dataframe
BrainSSurf.up <- data.frame(BrainSSURFACEmerged_unique.up)

#split swissprot accession string based on "|"
BrainSSurf_split.up <- data.frame(do.call('rbind', strsplit(as.character(BrainSSurf.up$BrainSSURFACEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
BrainSSurf_access.up <- BrainSSurf_split.up$X2

#Write accession to table for Blast2Go input
write.table(BrainSSurf_access.up, file = "SSbrainup.txt", col.names = NA) 

#############SULFIDIC SURFACE DOWNREGULATED#################
#Convert list of swissprot accesions into a dataframe
BrainSSurf.down <- data.frame(BrainSSURFACEmerged_unique.down)

#split swissprot accession string based on "|"
BrainSSurf_split.down <- data.frame(do.call('rbind', strsplit(as.character(BrainSSurf.down$BrainSSURFACEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
BrainSSurf_access.down <- BrainSSurf_split.down$X2

#Write accession to table for Blast2Go input
write.table(BrainSSurf_access.down, file = "SSbraindown.txt", col.names = NA) 

#############SULFIDIC CAVE UPREGULATED###################
#Convert list of swissprot accesions into a dataframe
BrainSCave.up <- data.frame(BrainSCAVEmerged_unique.up)

#split swissprot accession string based on "|"
BrainSCave_split.up <- data.frame(do.call('rbind', strsplit(as.character(BrainSCave.up$BrainSCAVEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
BrainSCave_access.up <- BrainSCave_split.up$X2

#Write accession to table for Blast2Go input
write.table(BrainSCave_access.up, file = "SCbrainup.txt", col.names = NA) 

#############SULFIDIC CAVE DOWNREGULATED#################
#Convert list of swissprot accesions into a dataframe
BrainSCave.down <- data.frame(BrainSCAVEmerged_unique.down)

#split swissprot accession string based on "|"
BrainSCave_split.down <- data.frame(do.call('rbind', strsplit(as.character(BrainSCave.down$BrainSCAVEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
BrainSCave_access.down <- BrainSCave_split.down$X2

#Write accession to table for Blast2Go input
write.table(BrainSCave_access.down, file = "SCbraindown.txt", col.names = NA) 

#############NONSULFIDIC CAVE UPREGULATED################
#Convert list of swissprot accesions into a dataframe
BrainNSCave.up <- data.frame(BrainNSCAVEmerged_unique.up)

#split swissprot accession string based on "|"
BrainNSCave_split.up <- data.frame(do.call('rbind', strsplit(as.character(BrainNSCave.up$BrainNSCAVEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
BrainNSCave_access.up <- BrainNSCave_split.up$X2

#Write accession to table for Blast2Go input
write.table(BrainNSCave_access.up, file = "NSCbrainup.txt", col.names = NA) 

#############NONSULFIDIC CAVE DOWNREGULATED##############
#Convert list of swissprot accesions into a dataframe
BrainNSCave.down <- data.frame(BrainNSCAVEmerged_unique.down)

#split swissprot accession string based on "|"
BrainNSCave_split.down <- data.frame(do.call('rbind', strsplit(as.character(BrainNSCave.down$BrainNSCAVEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
BrainNSCave_access.down <- BrainNSCave_split.down$X2

#Write accession to table for Blast2Go input
write.table(BrainNSCave_access.down, file = "NSCbraindown.txt", col.names = NA) 

################################################################################################
#
#          3B. Extract GO IDs associated with Biological processes from Blast2GO output
#
#
################################################################################################

#########################GILL###############################
#############NONSULFIDIC CAVE UPREGULATED###################
#read in B2GO output 
GillNScave.up <- read.delim("GO-results/NScavegillup_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
GillNScave_BP.up <- GillNScave.up[GillNScave.up$V9 == "P",]

#Subset only GO IDs and put in dataframe 
GillNScave_BPGO.up <- data.frame(GillNScave_BP.up$V5)

#remove duplicates
GillNScave_BPGO_unique.up  <- as.character(unique(unlist(GillNScave_BPGO.up)))

################SULFIDIC CAVE UPREGULATED###################
#read in B2GO output 
GillScave.up <- read.delim("GO-results/Scavegillup_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
GillScave_BP.up <- GillScave.up[GillScave.up$V9 == "P",]

#Subset only GO IDs and put in dataframe 
GillScave_BPGO.up <- data.frame(GillScave_BP.up$V5)

#remove duplicates
GillScave_BPGO_unique.up  <- as.character(unique(unlist(GillScave_BPGO.up)))

#############SULFIDIC SURFACE UPREGULATED###################
#read in B2GO output 
GillSsurface.up <- read.delim("GO-results/Ssurfacegillup_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
GillSsurface_BP.up <- GillSsurface.up[GillSsurface.up$V9 == "P",]

#Subset only GO IDs and put in dataframe
GillSsurface_BPGO.up <- data.frame(GillSsurface_BP.up$V5)

#remove duplicates
GillSsurface_BPGO_unique.up  <- as.character(unique(unlist(GillSsurface_BPGO.up)))

#############NONSULFIDIC CAVE DOWNREGULATED#################
#read in B2GO output 
GillNScave.down <- read.delim("GO-results/NScavegilldown_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
GillNScave_BP.down <- GillNScave.down[GillNScave.down$V9 == "P",]

#Subset only GO IDs and put in dataframe
GillNScave_BPGO.down <- data.frame(GillNScave_BP.down$V5)

#remove duplicates
GillNScave_BPGO_unique.down  <- as.character(unique(unlist(GillNScave_BPGO.down)))

################SULFIDIC CAVE DOWNREGULATED#################
#read in B2GO output 
GillScave.down <- read.delim("GO-results/Scavegilldown_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
GillScave_BP.down <- GillScave.down[GillScave.down$V9 == "P",]

#Subset only GO IDs and put in dataframe
GillScave_BPGO.down <- data.frame(GillScave_BP.down$V5)

#remove duplicates
GillScave_BPGO_unique.down  <- as.character(unique(unlist(GillScave_BPGO.down)))

################SULFIDIC SURFACE DOWNREGULATED#################
#read in B2GO output 
GillSsurface.down <- read.delim("GO-results/Ssurfacegilldown_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
GillSsurface_BP.down <- GillSsurface.down[GillSsurface.down$V9 == "P",]

#Subset only GO IDs and put in dataframe
GillSsurface_BPGO.down <- data.frame(GillSsurface_BP.down$V5)

#remove duplicates
GillSsurface_BPGO_unique.down  <- as.character(unique(unlist(GillSsurface_BPGO.down)))

#########################LIVER##############################
#############NONSULFIDIC CAVE UPREGULATED###################
#read in B2GO output 
LiverNScave.up <- read.delim("GO-results/NScaveliverup_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
LiverNScave_BP.up <- LiverNScave.up[LiverNScave.up$V9 == "P",]

#Subset only GO IDs and put in dataframe
LiverNScave_BPGO.up <- data.frame(LiverNScave_BP.up$V5)

#remove duplicates
LiverNScave_BPGO_unique.up  <- as.character(unique(unlist(LiverNScave_BPGO.up)))

################SULFIDIC CAVE UPREGULATED###################
#read in B2GO output 
LiverScave.up <- read.delim("GO-results/Scaveliverup_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
LiverScave_BP.up <- LiverScave.up[LiverScave.up$V9 == "P",]

#Subset only GO IDs and put in dataframe
LiverScave_BPGO.up <- data.frame(LiverScave_BP.up$V5)

#remove duplicates
LiverScave_BPGO_unique.up  <- as.character(unique(unlist(LiverScave_BPGO.up)))

#############SULFIDIC SURFACE UPREGULATED###################
#read in B2GO output 
LiverSsurface.up <- read.delim("GO-results/Ssurfaceliverup_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
LiverSsurface_BP.up <- LiverSsurface.up[LiverSsurface.up$V9 == "P",]

#Subset only GO IDs and put in dataframe
LiverSsurface_BPGO.up <- data.frame(LiverSsurface_BP.up$V5)

#remove duplicates
LiverSsurface_BPGO_unique.up  <- as.character(unique(unlist(LiverSsurface_BPGO.up)))

#############NONSULFIDIC CAVE DOWNREGULATED###################
#read in B2GO output 
LiverNScave.down <- read.delim("GO-results/NScaveliverdown_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
LiverNScave_BP.down <- LiverNScave.down[LiverNScave.down$V9 == "P",]

#Subset only GO IDs and put in dataframe
LiverNScave_BPGO.down <- data.frame(LiverNScave_BP.down$V5)

#remove duplicates
LiverNScave_BPGO_unique.down  <- as.character(unique(unlist(LiverNScave_BPGO.down)))

#############NONSULFIDIC CAVE DOWNREGULATED###################
#read in B2GO output 
LiverScave.down <- read.delim("GO-results/Scaveliverdown_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
LiverScave_BP.down <- LiverScave.down[LiverScave.down$V9 == "P",]

#Subset only GO IDs and put in dataframe
LiverScave_BPGO.down <- data.frame(LiverScave_BP.down$V5)

#remove duplicates
LiverScave_BPGO_unique.down  <- as.character(unique(unlist(LiverScave_BPGO.down)))

#############SULFIDIC SURFACE DOWNREGULATED###################
#read in B2GO output 
LiverSsurface.down <- read.delim("GO-results/Ssurfaceliverdown_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
LiverSsurface_BP.down <- LiverSsurface.down[LiverSsurface.down$V9 == "P",]

#Subset only GO IDs and put in dataframe
LiverSsurface_BPGO.down <- data.frame(LiverSsurface_BP.down$V5)

#remove duplicates
LiverSsurface_BPGO_unique.down  <- as.character(unique(unlist(LiverSsurface_BPGO.down)))

########################BRAIN###############################
#############NONSULFIDIC CAVE UPREGULATED###################
#read in B2GO output 
BrainNScave.up <- read.delim("GO-results/NScavebrainup_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
BrainNScave_BP.up <- BrainNScave.up[BrainNScave.up$V9 == "P",]

#Extract out GO terms associated with biological processes 
BrainNScave_BPGO.up <- data.frame(BrainNScave_BP.up$V5)

#Subset only GO IDs and put in dataframe
BrainNScave_BPGO_unique.up <- as.character(unique(unlist(BrainNScave_BPGO.up)))

################SULFIDIC CAVE UPREGULATED###################
#read in B2GO output 
BrainScave.up <- read.delim("GO-results/Scavebrainup_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
BrainScave_BP.up <- BrainScave.up[BrainScave.up$V9 == "P",]

#Extract out GO terms associated with biological processes 
BrainScave_BPGO.up <- data.frame(BrainScave_BP.up$V5)

#Subset only GO IDs and put in dataframe
BrainScave_BPGO_unique.up <- as.character(unique(unlist(BrainScave_BPGO.up)))

################SULFIDIC SURFACE UPREGULATED###################
#read in B2GO output 
BrainSsurface.up <- read.delim("GO-results/Ssurfacebrainup_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
BrainSsurface_BP.up <- BrainSsurface.up[BrainSsurface.up$V9 == "P",]

#Extract out GO terms associated with biological processes 
BrainSsurface_BPGO.up <- data.frame(BrainSsurface_BP.up$V5)

#Subset only GO IDs and put in dataframe
BrainSsurface_BPGO_unique.up <- as.character(unique(unlist(BrainSsurface_BPGO.up)))

#############NONSULFIDIC CAVE DOWNREGULATED#################
#read in B2GO output 
BrainNScave.down <- read.delim("GO-results/NScavebraindown_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
BrainNScave_BP.down <- BrainNScave.down[BrainNScave.down$V9 == "P",]

#Extract out GO terms associated with biological processes 
BrainNScave_BPGO.down <- data.frame(BrainNScave_BP.down$V5)

#Subset only GO IDs and put in dataframe
BrainNScave_BPGO_unique.down  <- as.character(unique(unlist(BrainNScave_BPGO.down)))

################SULFIDIC CAVE DOWNREGULATED#################
#read in B2GO output 
BrainScave.down <- read.delim("GO-results/Scavebraindown_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
BrainScave_BP.down <- BrainScave.down[BrainScave.down$V9 == "P",]

#Extract out GO terms associated with biological processes 
BrainScave_BPGO.down <- data.frame(BrainScave_BP.down$V5)

#Subset only GO IDs and put in dataframe
BrainScave_BPGO_unique.down  <- as.character(unique(unlist(BrainScave_BPGO.down)))

################SULFIDIC SURFACE DOWNREGULATED#################
#read in B2GO output 
BrainSsurface.down <- read.delim("GO-results/Ssurfacebraindown_GO.txt", header=F)

#Extract out GO terms associated with biological processes (denoted as "P")
BrainSsurface_BP.down <- BrainSsurface.down[BrainSsurface.down$V9 == "P",]

#Extract out GO terms associated with biological processes 
BrainSsurface_BPGO.down <- data.frame(BrainSsurface_BP.down$V5)

#Subset only GO IDs and put in dataframe
BrainSsurface_BPGO_unique.down  <- as.character(unique(unlist(BrainSsurface_BPGO.down)))

#######################################Venn Diagram Gill Cave Up-regulated#################################################
#Create unique variable, and construct venns 
universe.UP <- unique(c(GillNScave_BPGO_unique.up, GillScave_BPGO_unique.up))
GroupA <- universe.UP %in% GillNScave_BPGO_unique.up
GroupB <- universe.UP %in% GillScave_BPGO_unique.up
input.df <- data.frame(Cueva=GroupA, PSO=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#####################################Venn Diagram Gill Sulfur Up-regulated#################################################
#Create unique variable, and construct venns 
universe.UP <- unique(c(GillSsurface_BPGO_unique.up, GillScave_BPGO_unique.up))
GroupA <- universe.UP %in% GillSsurface_BPGO_unique.up
GroupB <- universe.UP %in% GillScave_BPGO_unique.up
input.df <- data.frame(Cueva=GroupA, PSO=GroupB)
colnames(input.df)=c("S Surface","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#####################################Venn Diagram Gill Cave Down-regulated#################################################
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(GillNScave_BPGO_unique.down, GillScave_BPGO_unique.down))
GroupA <- universe.DOWN %in% GillNScave_BPGO_unique.down
GroupB <- universe.DOWN %in% GillScave_BPGO_unique.down
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#####################################Venn Diagram Gill Sulfur Down-regulated###############################################
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(GillSsurface_BPGO_unique.down, GillScave_BPGO_unique.down))
GroupA <- universe.DOWN %in% GillSsurface_BPGO_unique.down
GroupB <- universe.DOWN %in% GillScave_BPGO_unique.down
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("S Surface","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#######################################Venn Diagram Liver Cave Up-regulated#################################################
#Create unique variable, and construct venns 
universe.UP <- unique(c(LiverNScave_BPGO_unique.up, LiverScave_BPGO_unique.up))
GroupA <- universe.UP %in% LiverNScave_BPGO_unique.up
GroupB <- universe.UP %in% LiverScave_BPGO_unique.up
input.df <- data.frame(Cueva=GroupA, PSO=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#####################################Venn Diagram Liver Sulfur Up-regulated#################################################
#Create unique variable, and construct venns 
universe.UP <- unique(c(LiverSsurface_BPGO_unique.up, LiverScave_BPGO_unique.up))
GroupA <- universe.UP %in% LiverSsurface_BPGO_unique.up
GroupB <- universe.UP %in% LiverScave_BPGO_unique.up
input.df <- data.frame(Cueva=GroupA, PSO=GroupB)
colnames(input.df)=c("S Surface","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#####################################Venn Diagram Liver Cave Down-regulated#################################################
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(LiverNScave_BPGO_unique.down, LiverScave_BPGO_unique.down))
GroupA <- universe.DOWN %in% LiverNScave_BPGO_unique.down
GroupB <- universe.DOWN %in% LiverScave_BPGO_unique.down
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#####################################Venn Diagram Liver Sulfur Down-regulated###############################################
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(LiverSsurface_BPGO_unique.down, LiverScave_BPGO_unique.down))
GroupA <- universe.DOWN %in% LiverSsurface_BPGO_unique.down
GroupB <- universe.DOWN %in% LiverScave_BPGO_unique.down
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("S Surface","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#######################################Venn Diagram Brain Cave Up-regulated#################################################
#Create unique variable, and construct venns 
universe.UP <- unique(c(BrainNScave_BPGO_unique.up, BrainScave_BPGO_unique.up))
GroupA <- universe.UP %in% BrainNScave_BPGO_unique.up
GroupB <- universe.UP %in% BrainScave_BPGO_unique.up
input.df <- data.frame(Cueva=GroupA, PSO=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#####################################Venn Diagram Brain Sulfur Up-regulated#################################################
#Create unique variable, and construct venns 
universe.UP <- unique(c(BrainSsurface_BPGO_unique.up, BrainScave_BPGO_unique.up))
GroupA <- universe.UP %in% BrainSsurface_BPGO_unique.up
GroupB <- universe.UP %in% BrainScave_BPGO_unique.up
input.df <- data.frame(PSO=GroupA, Cueva=GroupB)
colnames(input.df)=c("S Surface","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#####################################Venn Diagram Brain Cave Down-regulated#################################################
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(BrainNScave_BPGO_unique.down, BrainScave_BPGO_unique.down))
GroupA <- universe.DOWN %in% BrainNScave_BPGO_unique.down
GroupB <- universe.DOWN %in% BrainScave_BPGO_unique.down
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("NS Cave","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

#####################################Venn Diagram Brain Sulfur Down-regulated####################################################
#Create unique variable, and construct venns 
universe.DOWN <- unique(c(BrainSsurface_BPGO_unique.down, BrainScave_BPGO_unique.down))
GroupA <- universe.DOWN %in% BrainSsurface_BPGO_unique.down
GroupB <- universe.DOWN %in% BrainScave_BPGO_unique.down
input.df <- data.frame(Luna=GroupA, Cueva=GroupB)
colnames(input.df)=c("S Surface","S Cave")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)

################################################################################################
#
#          4. Weighted gene co-expression network (WGCNA)
#
#
################################################################################################
###Note this analyses was modified from Oldham et al. 2006 (see references)
#Below is for WGCNA https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/#cranInstall
# R version 3.1.2
#Bioconductor version 3.0 [bioC installer 1.16.5]

# Read table of expression values, row = gene, column = individual 
x <- read.delim(gzfile("tissue.txt.gz"),header=T, row.names=1)

# Remove rows with all zero entries (no expression)
x = x[rowSums(x) != 0,]
dim(x)

# Choose top 10000 expressed genes based on mean expression
select = order(rowMeans(x), decreasing =T)[1:10000]
x.common = x[select,]
dim(x.common)

# Round expression values to integer value
x.round = round(x.common)

# pull out the correct files
x.round <- subset(x.round, select = c(BB1, BB2, BB3, BB4,
                                      CB1, CB2, CB3,
                                      LB2, LB3,
                                      PSO2B1, PSO2B4,
                                      BG1, BG2, BG3, BG4,
                                      CG1, CG2, CG3, CG4, 
                                      LG1, LG2, LG4,
                                      PSO2G1, PSO2G2, PSO2G3,
                                      BL1, BL2, BL3, BL4,
                                      CL1, CL2, CL3, CL4,
                                      LL1, LL2, LL3, LL4,
                                      PSO2L1, PSO2L2, PSO2L3, PSO2L4))

###########################################################################################################################
#
#  4A. Set up Gene Coexpression Network Analysis
#
#
###########################################################################################################################

# Generate information for each of the possible variables, Sulfidic-NonSulfidic (fs.group), Population (pop.group), Light-NoLight (light.group), organ (tissue.group) and individual (sample.group)
fs.group = c(rep(1,4),rep(3,3),rep(1,2),rep(3,2),rep(1,4),rep(3,4),rep(1,3),rep(3,3),rep(1,4),rep(3,4),rep(1,4),rep(3,4))
pop.group = c(rep(1,4),rep(2,3),rep(3,2),rep(4,2),rep(1,4),rep(2,4),rep(3,3),rep(4,3),rep(1,4),rep(2,4),rep(3,4),rep(4,4))
light.group = c(rep(0,4),rep(1,3),rep(1,2),rep(0,2),rep(0,4),rep(1,4),rep(1,3),rep(0,3),rep(0,4),rep(1,4),rep(1,4),rep(0,4))
tissue.group = c(rep(0,11),rep(1,14),rep(2,16))
sample.group = c(1,2,3,4,5,6,7,10,11,13,16,1,2,3,4,5,6,7,8,9,10,12,13,14,15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)

# Create dataframe of variables 
colData.all = data.frame(cbind(fs.group,pop.group,light.group,tissue.group,sample.group))

# Add column names
colnames(colData.all) = c("sulfidic","population","light","tissue","sample")

# Add individual ids to each row
rownames(colData.all) = colnames(x.round)

# Convert pop column to factors because this is required for DESeq dataset design variable for transformation
colData.all[,2] = factor(colData.all[,2])
as.factor(pop.group)

# Set up DEseq data from count matrix, use drainage as design 
ddConsensus = DESeqDataSetFromMatrix(countData = x.round, colData = colData.all, design = ~population)

# Check DEseq matrix
ddConsensus

# Run transformation, blind=T to blind the transformation to the experimental design 
vsdConsensus <- varianceStabilizingTransformation(ddConsensus, blind=T)

# Pull out transformed data as matrix. rows = genes, columns = individuals
transformed.consensus= assay(vsdConsensus)

# Transpose and assign to new variable for transformed expression data
datExpr = t(assay(vsdConsensus))

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, blockSize=5000)

# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
abline(h=0.80,col="blue")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Set up the adjacency matrix with power 3, which is a local max with R^2 greater than 0.8
adjacencyAll = adjacency(datExpr,power=3,type="signed");
diag(adjacencyAll)=0
dissTOMAll   = 1-TOMsimilarity(adjacencyAll, TOMType="signed")
geneTreeAll  = hclust(as.dist(dissTOMAll), method="average")
collectGarbage()

# Set up Trait matrix
datTraits = data.frame(cbind(fs.group,pop.group,light.group,tissue.group,sample.group))
colnames(datTraits) = c("sulfidic","population","light","tissue","sample")
rownames(datTraits) = colnames(x.round)

# Plot gene clustering
plot(geneTreeAll,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (combined dataset)", labels=FALSE,hang=0.04);

# Plot relationship among samples 
sampleTree = hclust(dist(datExpr), method = "average");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Add trait heatmap to sample tree plot 
plotDendroAndColors(sampleTree, datTraits,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

# Set color mapping to NULL
mColorh=NULL

# Run cuttreeHybrid for four different minimum cluster sizes and deep split values 
for (ds in 0:3){
  tree = cutreeHybrid(dendro = geneTreeAll, pamStage=FALSE,
                      minClusterSize = (30-3*ds), cutHeight = 0.99, 
                      deepSplit = ds, distM = dissTOMAll)
  mColorh=cbind(mColorh,labels2colors(tree$labels));
}

# Plot gene tree with different split values and clusters plotted below
plotDendroAndColors(geneTreeAll, mColorh, paste("dpSplt =",0:3), main = "",dendroLabels=FALSE);

# Assign color to each value in trait matrix
traitColors = numbers2colors(datTraits, signed = FALSE);

# Choose cutoff of first split (ds=0)
moduleColors = mColorh[,3]
modulesAll = mColorh[,3]

# Plot gene tree with chosen module 
plotDendroAndColors(geneTreeAll, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")

# Set up calculation of Module Eigengenes 
PCs1A    = moduleEigengenes(datExpr,  colors=modulesAll) 
ME_1A    = PCs1A$eigengenes
distPC1A = 1-abs(cor(ME_1A,use="p"))
distPC1A = ifelse(is.na(distPC1A), 0, distPC1A)
pcTree1A = hclust(as.dist(distPC1A),method="a") 
MDS_1A   = cmdscale(as.dist(distPC1A),2)
colorsA1 = names(table(modulesAll))

# Plot Module Eigengene values, the below is best plotted into a pdf all together
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
plot(pcTree1A, xlab="",ylab="",main="",sub="")
plot(MDS_1A, col= colorsA1,  main="MDS plot", cex=2, pch=19)

ordergenes = geneTreeAll$order
plotMat(scale(log(t(datExpr)[ordergenes,])) , rlabels= modulesAll[ordergenes], clabels= colnames(t(datExpr)), rcols=modulesAll[ordergenes])

for (which.module in names(table(modulesAll))){
  ME = ME_1A[, paste("ME",which.module, sep="")] 
  barplot(ME, col=which.module, main="", cex.main=2, 
          ylab="eigengene expression",xlab="array sample") 
}; 

#########################################################################################################################
#
#         4B. WGCNA: Association with trait
#
#
#########################################################################################################################

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
datTraits = data.frame(cbind(fs.group,pop.group,light.group,tissue.group,sample.group))
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

###########################################################################################################################
#
#          5. ClusterProfiler and ReactomePA dot plots based on enrichment
#
#
###########################################################################################################################
# R version 3.1.2
#Bioconductor version 3.0 [bioC installer 1.16.5]

#Define biomart object
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")

############################GILL############################################################################
######################GILL Sulfidic Surface upregulated##################################################
#Convert list of swissprot accesions into a dataframe
GSS.up <- data.frame(GillSSURFACEmerged_unique.up)

#split swissprot accession string based on "|"
GSS_split.up <- data.frame(do.call('rbind', strsplit(as.character(GSS.up$GillSSURFACEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
GSS_access.up <- data.frame(GSS_split.up$X2)

##Query biomart to find entrez IDs using human accession IDs 
GSS_ENTREZ.up <- getBM(attributes = c("entrezgene"),
                       filters = "uniprotswissprot", values = GSS_access.up,
                       mart = mart)

#remove duplicates
GSS_ENTREZ_unique.up  <- as.character(unique(unlist(GSS_ENTREZ.up)))

######################GILL Sulfidic Cave upregulated##################################################
#Convert list of swissprot accesions into a dataframe
GSC.up <- data.frame(GillSCAVEmerged_unique.up)

#split swissprot accession string based on "|"
GSC_split.up <- data.frame(do.call('rbind', strsplit(as.character(GSC.up$GillSCAVEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
GSC_access.up <- data.frame(GSC_split.up$X2)

##Query biomart to find entrez IDs using human accession IDs 
GSC_ENTREZ.up <- getBM(attributes = c("entrezgene"),
                       filters = "uniprotswissprot", values = GSC_access.up,
                       mart = mart)

#remove duplicates
GSC_ENTREZ_unique.up  <- as.character(unique(unlist(GSC_ENTREZ.up)))

######################GILL Nonsulfidic Cave upregulated##################################################
#Convert list of swissprot accesions into a dataframe
GNSC.up <- data.frame(GillNSCAVEmerged_unique.up)

#split swissprot accession string based on "|"
GNSC_split.up <- data.frame(do.call('rbind', strsplit(as.character(GNSC.up$GillNSCAVEmerged_unique.up),'|',fixed=TRUE)))

#Extract out Human accession ID 
GNSC_access.up <- data.frame(GNSC_split.up$X2)

##Query biomart to find entrez IDs using human accession IDs 
GNSC_ENTREZ.up <- getBM(attributes = c("entrezgene"),
                       filters = "uniprotswissprot", values = GNSC_access.up,
                       mart = mart)

#remove duplicates
GNSC_ENTREZ_unique.up  <- as.character(unique(unlist(GNSC_ENTREZ.up)))

######################GILL Sulfidic Surface downregulated##################################################
#Convert list of swissprot accesions into a dataframe
GSS.down <- data.frame(GillSSURFACEmerged_unique.down)

#split swissprot accession string based on "|"
GSS_split.down <- data.frame(do.call('rbind', strsplit(as.character(GSS.down$GillSSURFACEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
GSS_access.down <- data.frame(GSS_split.down$X2)

##Query biomart to find entrez IDs using human accession IDs 
GSS_ENTREZ.down <- getBM(attributes = c("entrezgene"),
                       filters = "uniprotswissprot", values = GSS_access.down,
                       mart = mart)

#remove duplicates
GSS_ENTREZ_unique.down  <- as.character(unique(unlist(GSS_ENTREZ.down)))

######################GILL Sulfidic Cave downregulated##################################################
#Convert list of swissprot accesions into a dataframe
GSC.down <- data.frame(GillSCAVEmerged_unique.down)

#split swissprot accession string based on "|"
GSC_split.down <- data.frame(do.call('rbind', strsplit(as.character(GSC.down$GillSCAVEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
GSC_access.down <- data.frame(GSC_split.down$X2)

##Query biomart to find entrez IDs using human accession IDs 
GSC_ENTREZ.down <- getBM(attributes = c("entrezgene"),
                       filters = "uniprotswissprot", values = GSC_access.down,
                       mart = mart)

#remove duplicates
GSC_ENTREZ_unique.down  <- as.character(unique(unlist(GSC_ENTREZ.down)))

######################GILL Nonsulfidic Cave downregulated##################################################
#Convert list of swissprot accesions into a dataframe
GNSC.down <- data.frame(GillNSCAVEmerged_unique.down)

#split swissprot accession string based on "|"
GNSC_split.down <- data.frame(do.call('rbind', strsplit(as.character(GNSC.down$GillNSCAVEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
GNSC_access.down <- data.frame(GNSC_split.down$X2)

##Query biomart to find entrez IDs using human accession IDs 
GNSC_ENTREZ.down <- getBM(attributes = c("entrezgene"),
                         filters = "uniprotswissprot", values = GNSC_access.down,
                         mart = mart)

#remove duplicates
GNSC_ENTREZ_unique.down  <- as.character(unique(unlist(GNSC_ENTREZ.down)))

############################LIVER############################################################################
######################Liver Sulfidic Surface upregulated##################################################
#Convert list of swissprot accesions into a dataframe
LSS.up <- data.frame(LIVERSSURFACEmerged_unique.up)

#split swissprot accession string based on "|"
LSS_split.up <- data.frame(do.call('rbind', strsplit(as.character(LSS.up$LIVERSSURFACEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
LSS_access.up <- data.frame(LSS_split.up$X2)

##Query biomart to find entrez IDs using human accession IDs 
LSS_ENTREZ.up <- getBM(attributes = c("entrezgene"),
                       filters = "uniprotswissprot", values = LSS_access.up,
                       mart = mart)

#remove duplicates
LSS_ENTREZ_unique.up  <- as.character(unique(unlist(LSS_ENTREZ.up)))

######################Liver Sulfidic cave upregulated##################################################
#Convert list of swissprot accesions into a dataframe
LSC.up <- data.frame(LIVERSCAVEmerged_unique.up)

#split swissprot accession string based on "|"
LSC_split.up <- data.frame(do.call('rbind', strsplit(as.character(LSC.up$LIVERSCAVEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
LSC_access.up <- data.frame(LSC_split.up$X2)

##Query biomart to find entrez IDs using human accession IDs 
LSC_ENTREZ.up <- getBM(attributes = c("entrezgene"),
                       filters = "uniprotswissprot", values = LSC_access.up,
                       mart = mart)

#remove duplicates
LSC_ENTREZ_unique.up  <- as.character(unique(unlist(LSC_ENTREZ.up)))

######################Liver Nonsulfidic cave upregulated##################################################
#Convert list of swissprot accesions into a dataframe
LNSC.up <- data.frame(LIVERNSCAVEmerged_unique.up)

#split swissprot accession string based on "|"
LNSC_split.up <- data.frame(do.call('rbind', strsplit(as.character(LNSC.up$LIVERNSCAVEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
LNSC_access.up <- data.frame(LNSC_split.up$X2)

##Query biomart to find entrez IDs using human accession IDs 
LNSC_ENTREZ.up <- getBM(attributes = c("entrezgene"),
                       filters = "uniprotswissprot", values = LNSC_access.up,
                       mart = mart)

#remove duplicates
LNSC_ENTREZ_unique.up  <- as.character(unique(unlist(LNSC_ENTREZ.up)))

######################Liver Sulfidic surface downregulated##################################################
#Convert list of swissprot accesions into a dataframe
LSS.down <- data.frame(LIVERSSURFACEmerged_unique.down)

#split swissprot accession string based on "|"
LSS_split.down <- data.frame(do.call('rbind', strsplit(as.character(LSS.down$LIVERSSURFACEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
LSS_access.down <- data.frame(LSS_split.down$X2)

##Query biomart to find entrez IDs using human accession IDs 
LSS_ENTREZ.down <- getBM(attributes = c("entrezgene"),
                       filters = "uniprotswissprot", values = LSS_access.down,
                       mart = mart)

#remove duplicates
LSS_ENTREZ_unique.down  <- as.character(unique(unlist(LSS_ENTREZ.down)))

######################Liver Sulfidic cave downregulated##################################################
#Convert list of swissprot accesions into a dataframe
LSC.down <- data.frame(LIVERSCAVEmerged_unique.down)

#split swissprot accession string based on "|"
LSC_split.down <- data.frame(do.call('rbind', strsplit(as.character(LSC.down$LIVERSCAVEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
LSC_access.down <- data.frame(LSC_split.down$X2)

##Query biomart to find entrez IDs using human accession IDs 
LSC_ENTREZ.down <- getBM(attributes = c("entrezgene"),
                       filters = "uniprotswissprot", values = LSC_access.down,
                       mart = mart)

#remove duplicates
LSC_ENTREZ_unique.down  <- as.character(unique(unlist(LSC_ENTREZ.down)))

######################Liver Nonsulfidic cave upregulated##################################################
#Convert list of swissprot accesions into a dataframe
LNSC.down <- data.frame(LIVERNSCAVEmerged_unique.down)

#split swissprot accession string based on "|"
LNSC_split.down <- data.frame(do.call('rbind', strsplit(as.character(LNSC.down$LIVERNSCAVEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
LNSC_access.down <- data.frame(LNSC_split.down$X2)

##Query biomart to find entrez IDs using human accession IDs 
LNSC_ENTREZ.down <- getBM(attributes = c("entrezgene"),
                         filters = "uniprotswissprot", values = LNSC_access.down,
                         mart = mart)

#remove duplicates
LNSC_ENTREZ_unique.down  <- as.character(unique(unlist(LNSC_ENTREZ.down)))

#########################BRAIN############################################################################
######################Brain Sulfidic Surface upregulated##################################################
#Convert list of swissprot accesions into a dataframe
BSS.up <- data.frame(BrainSSURFACEmerged_unique.up)

#split swissprot accession string based on "|"
BSS_split.up <- data.frame(do.call('rbind', strsplit(as.character(BSS.up$BrainSSURFACEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
BSS_access.up <- data.frame(BSS_split.up$X2)

##Query biomart to find entrez IDs using human accession IDs 
BSS_ENTREZ.up <- getBM(attributes = c("entrezgene"),
                       filters = "uniprotswissprot", values = BSS_access.up,
                       mart = mart)

#remove duplicates
BSS_ENTREZ_unique.up  <- as.character(unique(unlist(BSS_ENTREZ.up)))

##################### Brain Sulfidic Cave upregulated##################################################
#Convert list of swissprot accesions into a dataframe
BSC.up <- data.frame(BrainSCAVEmerged_unique.up)

#split swissprot accession string based on "|"
BSC_split.up <- data.frame(do.call('rbind', strsplit(as.character(BSC.up$BrainSCAVEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
BSC_access.up <- data.frame(BSC_split.up$X2)

##Query biomart to find entrez IDs using human accession IDs 
BSC_ENTREZ.up <- getBM(attributes = c("entrezgene"),
                        filters = "uniprotswissprot", values = BSC_access.up,
                        mart = mart)

#remove duplicates
BSC_ENTREZ_unique.up  <- as.character(unique(unlist(BSC_ENTREZ.up)))

##################### Brain NonSulfidic Cave upregulated##################################################
#Convert list of swissprot accesions into a dataframe
BNSC.up <- data.frame(BrainNSCAVEmerged_unique.up)

#split swissprot accession string based on "|"
BNSC_split.up <- data.frame(do.call('rbind', strsplit(as.character(BNSC.up$BrainNSCAVEmerged_unique.up),'|',fixed=TRUE)))

#Extract out human accession ID 
BNSC_access.up <- data.frame(BNSC_split.up$X2)

##Query biomart to find entrez IDs using human accession IDs 
BNSC_ENTREZ.up <- getBM(attributes = c("entrezgene"),
                          filters = "uniprotswissprot", values = BNSC_access.up,
                          mart = mart)

#remove duplicates
BNSC_ENTREZ_unique.up  <- as.character(unique(unlist(BNSC_ENTREZ.up)))

##################### Brain NonSulfidic Cave downregulated##################################################
#Convert list of swissprot accesions into a dataframe
BNSC.down <- data.frame(BrainNSCAVEmerged_unique.down)

#split swissprot accession string based on "|"
BNSC_split.down <- data.frame(do.call('rbind', strsplit(as.character(BNSC.down$BrainNSCAVEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
BNSC_access.down <- data.frame(BNSC_split.down$X2)

##Query biomart to find entrez IDs using human accession IDs 
BNSC_ENTREZ.down <- getBM(attributes = c("entrezgene"),
                         filters = "uniprotswissprot", values = BNSC_access.down,
                         mart = mart)

#remove duplicates
BNSC_ENTREZ_unique.down  <- as.character(unique(unlist(BNSC_ENTREZ.down)))

##################### Brain Sulfidic cave downregulated##################################################
#Convert list of swissprot accesions into a dataframe
BSC.down <- data.frame(BrainSCAVEmerged_unique.down)

#split swissprot accession string based on "|"
BSC_split.down <- data.frame(do.call('rbind', strsplit(as.character(BSC.down$BrainSCAVEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
BSC_access.down <- data.frame(BSC_split.down$X2)

##Query biomart to find entrez IDs using human accession IDs 
BSC_ENTREZ.down <- getBM(attributes = c("entrezgene"),
                         filters = "uniprotswissprot", values = BSC_access.down,
                         mart = mart)

#remove duplicates
BSC_ENTREZ_unique.down  <- as.character(unique(unlist(BSC_ENTREZ.down)))

##################### Brain Sulfidic surface downregulated##################################################
#Convert list of swissprot accesions into a dataframe
BSS.down <- data.frame(BrainSSURFACEmerged_unique.down)

#split swissprot accession string based on "|"
BSS_split.down <- data.frame(do.call('rbind', strsplit(as.character(BSS.down$BrainSSURFACEmerged_unique.down),'|',fixed=TRUE)))

#Extract out human accession ID 
BSS_access.down <- data.frame(BSS_split.down$X2)

##Query biomart to find entrez IDs using human accession IDs 
BSS_ENTREZ.down <- getBM(attributes = c("entrezgene"),
                           filters = "uniprotswissprot", values = BSS_access.down,
                           mart = mart)

#remove duplicates
BSS_ENTREZ_unique.down  <- as.character(unique(unlist(BSS_ENTREZ.down)))

##################################################################################################
#          
#            5A. EnrichPathway analyses and dot plots
#
#
##################################################################################################

############################NONSULFIDIC CAVE#####################################################
############################UPREGULATED##########################################################
#Create a list containing all three organs for each ecotype   
Luna.up <- list(GNSC_ENTREZ_unique.up = GNSC_ENTREZ_unique.up, LNSC_ENTREZ_unique.up = LNSC_ENTREZ_unique.up, BNSC_ENTREZ_unique.up = BNSC_ENTREZ_unique.up)

#Compare gene clusters (entrez IDs) functional profiles using compareCluster
LunaOrgan.up <- compareCluster(Luna.up, fun = "enrichPathway", organism = "human",  pvalueCutoff = 0.05)

#Summarize compareCluster output
head(summary(LunaOrgan.up))

#construct a dotplot of the enriched terms, colored by the number of genes in the category. 
plot(LunaOrgan.up, type="dot", colorBy = "GeneRatio")

#Extract highest GeneRatio's for all comparisons (note this is based on plot generated above)
LUNAsummary.up <- summary(LunaOrgan.up)
LUNAGeneRatio.up <- subset(LUNAsummary.up, select = c(Cluster, Description, GeneRatio))
LUNAGGeneRatio.up <- subset(LUNAGeneRatio.up, Description == "Metabolism" | Description == "Metabolism of lipids and lipoproteins"
                            | Description == "Transmembrane transport of small molecules" | Description == "Ion channel transport"
                            | Description == "Transport of inorganic cations/anions and amino acids/oligopeptides"
                            | Description == "SLC-mediated transmembrane transport" | Description == "Cholesterol biosynthesis"
                            | Description == "Activation of gene expression by SREBF (SREBP)" | Description == "Regulation of cholesterol biosynthesis by SREBP (SREBF)")

##########################DOWNREGULATED##########################################################
#Create a list containing all three organs for each ecotype   
Luna.down <- list(GNSC_ENTREZ_unique.down = GNSC_ENTREZ_unique.down, LNSC_ENTREZ_unique.down = LNSC_ENTREZ_unique.down, BNSC_ENTREZ_unique.down = BNSC_ENTREZ_unique.down)

#Compare gene clusters (entrez IDs) functional profiles using compareCluster
LunaOrgan.down <- compareCluster(Luna.down, fun = "enrichPathway", organism = "human", pvalueCutoff = 0.05)

#Summarize compareCluster output
head(summary(LunaOrgan.down))

#construct a dotplot of the enriched terms, colored by the number of genes in the category. 
plot(LunaOrgan.down, type="dot", colorBy = "GeneRatio")

#Extract highest GeneRatio's for all comparisons (note this is based on plot generated above)
LUNAsummary.down <- summary(LunaOrgan.down)
LUNAGeneRatio.down <- subset(LUNAsummary.down, select = c(Cluster, Description, GeneRatio))
LUNAGGeneRatio.down <- subset(LUNAGeneRatio.down, Description == "Transmembrane transport of small molecules" | Description == "Hemostasis"
                              | Description == "Neuronal System" | Description == "Axon guidance"
                              | Description == "L1CAM interactions" | Description == "Metabolism" | Description == "Developmental Biology"
                              | Description == "SLC-mediated transmembrane transport" | Description == "Immune System" | Description == "Innate Immune System"
                              | Description == "Adaptive Immune System" | Description == "Cytokine Signaling in Immune system")

###############################SULFIDIC CAVE#####################################################
############################UPREGULATED##########################################################
#Create a list containing all three organs for each ecotype   
Cueva.up <- list(GSC_ENTREZ_unique.up = GSC_ENTREZ_unique.up, LSC_ENTREZ_unique.up = LSC_ENTREZ_unique.up, BSC_ENTREZ_unique.up = BSC_ENTREZ_unique.up)

#Compare gene clusters (entrez IDs) functional profiles using compareCluster
CuevaOrgan.up <- compareCluster(geneCluster = Cueva.up, fun = "enrichPathway", organism = "human", pvalueCutoff = 0.05)

#Summarize compareCluster output
head(summary(CuevaOrgan.up))

#construct a dotplot of the enriched terms, colored by the number of genes in the category. 
plot(CuevaOrgan.up, type="dot", colorBy = "GeneRatio")

#Extract highest GeneRatio's for all comparisons (note this is based on plot generated above)
CUEVAsummary.up <- summary(CuevaOrgan.up)
CUEVAGeneRatio.up <- subset(CUEVAsummary.up, select = c(Cluster, Description, GeneRatio))
CUEVAGGeneRatio.up <- subset(CUEVAGeneRatio.up, Description == "Metabolism" | Description == "Immune System"
                             | Description == "Cell Cycle" | Description == "Cell Cycle, Mitotic"
                             | Description == "Adaptive Immune System" | Description == "Transmembrane transport of small molecules"
                             | Description == "Metabolism of lipids and lipoproteins" | Description == "SLC-mediated transmembrane transport" 
                             | Description == "Glycogen storage diseases")

##########################DOWNREGULATED##########################################################
#Create a list containing all three organs for each ecotype   
Cueva.down <- list(GSC_ENTREZ_unique.down = GSC_ENTREZ_unique.down, LSC_ENTREZ_unique.down = LSC_ENTREZ_unique.down, BSC_ENTREZ_unique.down = BSC_ENTREZ_unique.down)

#Compare gene clusters (entrez IDs) functional profiles using compareCluster
CuevaOrgan.down <- compareCluster(Cueva.down, fun = "enrichPathway", organism = "human", pvalueCutoff = 0.05)

#Summarize compareCluster output
head(summary(CuevaOrgan.down))

#construct a dotplot of the enriched terms, colored by the number of genes in the category. 
plot(CuevaOrgan.down, type="dot", colorBy = "GeneRatio")

#Extract highest GeneRatio's for all comparisons (note this is based on plot generated above)
CUEVAsummary.down <- summary(CuevaOrgan.down)
CUEVAGeneRatio.down <- subset(CUEVAsummary.down, select = c(Cluster, Description, GeneRatio))
CUEVAGGeneRatio.down <- subset(CUEVAGeneRatio.down, Description == "Metabolism" | Description == "Metabolism of lipids and lipoproteins"
                               | Description == "Extracellular matrix organization" | Description == "Glycogen storage diseases"
                               | Description == "Myoclonic epilepsy of Lafora" | Description == "Transmembrane transport of small molecules"
                               | Description == "SLC-mediated transmembrane transport" | Description == "Ion channel transport" 
                               | Description == "Transport of inorganic cations/anions and amino acids/oligopeptides" | Description == "Fatty acid, triacylglycerol, and ketone body metabolism"
                               | Description == "Cholesterol biosynthesis" | Description == "Activation of gene expression by SREBF (SREBP)"
                               | Description == "Regulation of cholesterol biosynthesis by SREBP (SREBF)")

############################SULFIDIC SURFACE#####################################################
############################UPREGULATED##########################################################
#Create a list containing all three organs for each ecotype  
PSO.up <- list(GSS_ENTREZ_unique.up = GSS_ENTREZ_unique.up, LSS_ENTREZ_unique.up = LSS_ENTREZ_unique.up, BSS_ENTREZ_unique.up = BSS_ENTREZ_unique.up)

#Compare gene clusters (entrez IDs) functional profiles using compareCluster
PSOOrgan.up <- compareCluster(PSO.up, fun = "enrichPathway", organism = "human", pvalueCutoff = 0.05)

#Summarize compareCluster output
head(summary(PSOOrgan.up))

#construct a dotplot of the enriched terms, colored by the number of genes in the category. 
plot(PSOOrgan.up, type="dot", colorBy = "GeneRatio")

#Extract highest GeneRatio's for all comparisons (note this is based on plot generated above)
PSOsummary.up <- summary(PSOOrgan.up)
PSOGeneRatio.up <- subset(PSOsummary.up, select = c(Cluster, Description, GeneRatio))
PSOGGeneRatio.up <- subset(PSOGeneRatio.up, Description == "Glucose metabolism" | Description == "Metabolism"
                           | Description == "Glycogen storage diseases" | Description == "Myoclonic epilepsy of Lafora"
                           | Description == "Metabolism of carbohydrates")

##########################DOWNREGULATED##########################################################
#Create a list containing all three organs for each ecotype 
PSO.down <- list(GSS_ENTREZ_unique.down = GSS_ENTREZ_unique.down, LSS_ENTREZ_unique.down = LSS_ENTREZ_unique.down, BSS_ENTREZ_unique.down = BSS_ENTREZ_unique.down)

#Compare gene clusters (entrez IDs) functional profiles using compareCluster
PSOOrgan.down <- compareCluster(PSO.down, fun = "enrichPathway", organism = "human", pvalueCutoff = 0.05)

#Summarize compareCluster output
head(summary(PSOOrgan.down))

#construct a dotplot of the enriched terms, colored by the number of genes in the category. 
plot(PSOOrgan.down, type="dot", colorBy = "GeneRatio")

#Extract highest GeneRatio's for all comparisons (note this is based on plot generated above)
PSOsummary.down <- summary(PSOOrgan.down)
PSOGeneRatio.down <- subset(PSOsummary.down, select = c(Cluster, Description, GeneRatio))
PSOGGeneRatio.down <- subset(PSOGeneRatio.down, Description == "Extracellular matrix organization" | Description == "Transmembrane transport of small molecules"
                             | Description == "SLC-mediated transmembrane transport" | Description == "Ion channel transport"
                             | Description == "Transport of inorganic cations/anions and amino acids/oligopeptides" | Description == "Amino acid transport across the plasma membrane"
                             | Description == "Metabolism" | Description == "Metabolism of lipids and lipoproteins")

#################################################################################################
#          
#            5B. Pearson Correlations-coefficients 
#
#
##################################################################################################
## To measure the strength of association between enrich GO terms, we took the extracted the top GeneRatio's for each comparison (based on the dot plots above)
## and made variables based on compareCluster output. If there was no evidence of enrichment (i.e no value for the enriched GO term)
## we added a "0" to the variable. We then ran the rcorr function, which produces correlations/covariances for pearson correlations 
## between all three organs (gill vs. liver, gill vs. brain and brain vs. liver) and reported the Pearson correlation coefficients (r) 
## and P-values in the supplement which indicates similarities in enrichment among organs.

###############################NONSULFIDIC CAVE#################################################
#Calculate Pearson Correlations based on the top gene Ratio's extracted in compareCluster Summary (see above)
NSCB.up <- c(0.529, 0.235, 0, 0, 0, 0, 0, 0, 0)
NSCG.up <- c(0, 0.254, 0, 0.119, 0.090, 0.134, 0, 0, 0)
NSCL.up <- c(0.551, 0.297, 0, 0, 0, 0, 0.103, 0.097, 0.097)
rcorr(NSCB.up, NSCG.up, type="pearson")
rcorr(NSCB.up, NSCL.up, type="pearson")
rcorr(NSCG.up, NSCL.up, type="pearson")

#Calculate Pearson Correlations based on the top gene Ratio's extracted in compareCluster Summary (see above)
NSCB.down <- c(0.200,0.165,0.153,0.118,0.094,0,0,0,0,0,0,0)
NSCG.down <- c(0.170,0,0,0.094,0.066,0.349,0.123,0.104,0,0,0,0)
NSCL.down <- c(0, 0.144,0,0,0,0,0,0,0.273,0.160,0.134,0.096)
rcorr(NSCB.down, NSCG.down, type="pearson")
rcorr(NSCB.down, NSCL.down, type="pearson")
rcorr(NSCG.down, NSCL.down, type="pearson")

###############################SULFIDIC SURFACE#################################################
#Calculate Pearson Correlations based on the top gene Ratio's extracted in compareCluster Summary (see above)
SSMGU.up <- c(0.059, 0.355,0.091,0.091, 0.091)
SSLU.up <- c(0,0,0,0,0)
SSBU.up <- c(0.086,0.371,0,0,0)
rcorr(SSMGU.up, SSLU.up, type = "pearson")
rcorr(SSMGU.up, SSBU.up, type = "pearson")
rcorr(SSLU.up, SSBU.up, type = 'pearson')

#Calculate Pearson Correlations based on the top gene Ratio's extracted in compareCluster Summary (see above)
SSMG.down <- c(0,0,0,0.190,0.098,0.078,0.059,0.046)
SSL.down <- c(0,0.516,0.226,0,0,0,0,0)
SSB.down <- c(0.167,0,0,0,0,0,0,0)
rcorr(SSMG.down, SSL.down, type = "pearson")
rcorr(SSMG.down, SSB.down, type = "pearson")
rcorr(SSL.down, SSB.down, type = 'pearson')

##################################SULFIDIC CAVE#################################################
#Calculate Pearson Correlations based on the top gene Ratio's extracted in compareCluster Summary (see above)
SSMG.up <- c(0.407, 0, 0, 0, 0, 0.137, 0.121, 0.176, 0.088)
SCL.up <- c(0.359, 0, 0, 0, 0, 0.167, 0.103, 0, 0)
SCB.up <- c(0.323, 0.250, 0.183, 0.177, 0.152, 0, 0, 0, 0)
rcorr(SSMG.up, SCL.up, type="pearson")
rcorr(SSMG.up, SCB.up, type="pearson")
rcorr(SCL.up, SCB.up, type="pearson")

#Calculate Pearson Correlations based on the top gene Ratio's extracted in compareCluster Summary (see above)
SCMG.down <- c(0, 0, 0, 0, 0, 0, 0, 0, 0.208, 0.111, 0.083, 0.069, 0.56)
SCL.down <- c(0.564, 0.302, 0, 0, 0, 0.121, 0.107, 0.107, 0, 0, 0, 0, 0.107)
SCB.down <- c(0.398, 0.136, 0.107, 0.107, 0.107, 0, 0, 0, 0, 0.087, 0, 0, 0.078)
rcorr(SCMG.down, SCL.down, type="pearson")
rcorr(SCMG.down, SCB.down, type="pearson")
rcorr(SCL.down, SCB.down, type="pearson")

