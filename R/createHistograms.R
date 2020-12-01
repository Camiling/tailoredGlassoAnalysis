setwd("~/Documents/Master/Paper")
library(huge)
library(glasso)
library(igraph)
library(Matrix)
library(ggplot2)
library(gridExtra)
source('~/Documents/Master/oppgaven/Alle filer/code/useful.functions.R') 
source('~/Documents/Master/oppgaven/Alle filer/code/WeightTuning/weighttuning_functions.R') 

# Section:CREATE HISTOGRAMS OF THE PARTIAL CORRELATIONS IN SIMULATED DATA -----------------------------------

n=80
p=100
# Case 1: No changed edges, large corr. 
priormat.case1 = create.sample.sf.prior.final(frac.to.mutate=0,v=0.5,n,p)
priorcorrs.case1 = priormat.case1[!diag(rep(1,p))]
priorcorrs.case1 = priorcorrs.case1[which(priorcorrs.case1!=0,arr.ind = T)]
df.case1= data.frame(PartialCor=priorcorrs.case1)

# Case 2: No changed edges, smaller corr. 
priormat.case2 = create.sample.sf.prior.final(frac.to.mutate=0,v=0.03,n,p) 
priorcorrs.case2 = priormat.case2[!diag(rep(1,p))]
priorcorrs.case2 = priorcorrs.case2[which(priorcorrs.case2!=0,arr.ind = T)]
df.case2 = data.frame(PartialCor=priorcorrs.case2)

pdf('~/Documents/Master/Paper/simcorrHist.pdf',8,6)
hist.1 = ggplot(df.case1,aes(x=PartialCor))+labs(title='(a)')+ylab('Frequency')+ xlab('Partial correlations')+geom_histogram(color='deepskyblue2',fill='deepskyblue2', binwidth = 0.005)+theme(legend.position ="none")+xlim(0,0.22)
hist.2 = ggplot(df.case2,aes(x=PartialCor))+labs(title='(b)')+ylab('Frequency')+  xlab('Partial correlations')+ geom_histogram(color='deepskyblue2',fill='deepskyblue2', binwidth = 0.005)+theme(legend.position ="none")+xlim(0,0.22)
grid.arrange(hist.1,hist.2,ncol=2,nrow=1)
dev.off()


# Section:CREATE HISTOGRAMS OF THE PARTIAL CORRELATIONS IN GENOMIC DATA -----------------------------------



# TCGA

load("~/Documents/Master/oppgaven/Alle filer/code/datasetttonje/Newest/RNA_tcga.RData")
# GLASSO ON RNA
rna.fit = huge(RNA.tcga,nlambda=35,method='glasso')
# SELECT BY StARS
set.seed(12345)
rna.select.stars=huge.select(rna.fit, criterion='stars',stars.thresh=0.05)
theta.rna.stars = rna.select.stars$opt.icov # Optimal precision matrix
prior.mat.rna.tcga= as.matrix(abs(cov2cor(forceSymmetric(theta.rna.stars,uplo='U'))+0)) # The prior corr matrix
nonzero.weights.rna = prior.mat.rna.tcga[!diag(rep(1,ncol(RNA.tcga)))]
nonzero.weights.rna=nonzero.weights.rna[which(nonzero.weights.rna!=0,arr.ind = T)]
df.tcga = data.frame(PartialCor=nonzero.weights.rna)


# Oslo 2
load("~/Documents/Master/oppgaven/Alle filer/code/OSLO2/OSLO2RNA.RData")
oslo2.rna.fit = huge(OSLO2.RNA, nlambda=35,method='glasso')
# STARS
set.seed(1235)
oslo2.rna.stars=huge.select(oslo2.rna.fit,criterion='stars',stars.thresh=0.05)
prior.mat.rna.oslo2 = as.matrix(abs(cov2cor(forceSymmetric(oslo2.rna.stars$opt.icov,uplo='U'))+0)) # The prior corr matrix
nonzero.weights.rppa.oslo2 = prior.mat.rna.oslo2[!diag(rep(1,ncol(OSLO2.RNA)))]
df.oslo2 = data.frame(PartialCor=nonzero.weights.rppa.oslo2[nonzero.weights.rppa.oslo2!=0])

pdf('~/Documents/Master/Paper/GenomiccorrHist.pdf',8,6)
hist.3 = ggplot(df.tcga,aes(x=PartialCor))+labs(title='(a)')+ylab('Frequency')+ xlab('Partial correlations')+geom_histogram(color='darkolivegreen3',fill='darkolivegreen3', binwidth = 0.005)+theme(legend.position ="none")+xlim(0,0.27)
hist.4 = ggplot(df.oslo2,aes(x=PartialCor))+labs(title='(b)')+ylab('Frequency')+ xlab('Partial correlations')+geom_histogram(color='darkolivegreen3',fill='darkolivegreen3', binwidth = 0.005)+theme(legend.position ="none")+xlim(0,0.27)
grid.arrange(hist.3,hist.4,ncol=2,nrow=1)
dev.off()

# Also create one single plot for all: 

pdf('~/Documents/Master/Paper/corrHist.pdf',8,8)
hist.3.2 = ggplot(df.tcga,aes(x=PartialCor))+labs(title='(c)')+ylab('Frequency')+ xlab('Partial correlations')+geom_histogram(color='darkolivegreen3',fill='darkolivegreen3', binwidth = 0.005)+theme(legend.position ="none")+xlim(0,0.27)
hist.4.2 = ggplot(df.oslo2,aes(x=PartialCor))+labs(title='(d)')+ylab('Frequency')+ xlab('Partial correlations')+geom_histogram(color='darkolivegreen3',fill='darkolivegreen3', binwidth = 0.005)+theme(legend.position ="none")+xlim(0,0.27)
grid.arrange(hist.1, hist.2, hist.3.2,hist.4.2,ncol=2,nrow=2)
dev.off()

