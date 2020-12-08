# Dowload, preprocess and analyse TCGA data with the tailored graphical lasso

# Neccessary additions
remotes::install_github("Camiling/tailoredGlasso", auth_token = "c19521d5a82ba1a631ffa5de9a569924908c10e4", dependencies = TRUE, build_vignettes = TRUE)
library(tailoredGlasso)
library(huge)
library(glasso)
library(igraph)
library(MASS)
library(ggplot2)
library(gridExtra)
library(GGally)
library(network)

# Section: preprocess data --------------------------------

# Start by downloading data

# Download protein data from TCGA database
utils::download.file("https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/RPPA.gz", "R/RPPA.gz")
RPPA.full <- utils::read.table(gzfile("R/RPPA.gz"))
# Ensure row names and column names are not their own row or column
rownames(RPPA.full) <- (RPPA.full[, 1])
RPPA.full <- RPPA.full[, 2:ncol(RPPA.full)]

# Must extract patient IDs
names.rppa <- unlist(lapply(RPPA.full, FUN = function(l) as.character(l[1])))[1:ncol(RPPA.full)]
names.rppa <- gsub("-", ".", names.rppa)
names(names.rppa) <- NULL
# Fill in values in matrix
RPPA.withIDs <- matrix(0, ncol = ncol(RPPA.full), nrow = nrow(RPPA.full) - 1)
for (i in 1:ncol(RPPA.full)) {
  RPPA.withIDs[, i] <- as.numeric(as.character(RPPA.full[[i]][2:nrow(RPPA.full)]))
}
rownames(RPPA.withIDs) <- rownames(RPPA.full)[2:nrow(RPPA.full)]
colnames(RPPA.withIDs) <- names.rppa
dim(RPPA.withIDs) # 281 x 937
RPPA.full <- RPPA.withIDs

# Only use samples from patients with cancer (indicated by 01 at the end of the sample name)
sampleType <- rep(NA, times = dim(RPPA.full)[2])
for (j in 1:dim(RPPA.full)[2]) {
  sampleType[j] <- substr(colnames(RPPA.full)[j], 14, 15)
}
RPPA <- RPPA.full[, which(sampleType == "01")] # Only patients with cancer

# Map the names of the antibodies used to identify the proteins to the names of the genes that encode them
# Use file showing which genes the proteins detected by the antibodies are encoded from (can be downloaded from TCGA website)
rppa.to.gene <- read.table("annotation_data/RPPA_to_gene.txt", sep = "\t", stringsAsFactors = F)
# We need only the id and the gene name
rppa.to.gene <- rppa.to.gene[, 1:2]
names(rppa.to.gene)[1] <- "id" # Rename
names(rppa.to.gene)[2] <- "gene" # Rename
# Remove last 4 letter of the rownames of the RPPA data frame, to see which antibody we have regardless of animal used
RPPA.geneid <- rep(NA, length(rownames(RPPA)))
for (i in 1:length(rownames(RPPA))) {
  RPPA.geneid[i] <- substr(rownames(RPPA)[i], 1, nchar(rownames(RPPA)[i]) - 4)
}
# Match names in RPPA data frame to the antibody names in the mapping frame we just created
rppa.genes <- rppa.to.gene$gene[match(RPPA.geneid, rppa.to.gene$id)]
# Must check for match if last letters are included as well:
for (i in 1:length(rppa.genes)) {
  if (rppa.genes[i] == "" | is.na(rppa.genes[i])) {
    rppa.genes[i] <- rppa.to.gene$gene[match(rownames(RPPA)[i], rppa.to.gene$id)]
  }
}

# Since not all antibodies were mapped to genes, we also use a mapping downloaded from the University of Stanford
rppa.to.gene2 <- read.table("annotation_data/Allprotein.csv", header = T, sep = ",", stringsAsFactors = F)
# Only check the antibodies that have no match yet.
for (i in 1:length(rppa.genes)) {
  if (rppa.genes[i] == "" | is.na(rppa.genes[i])) {
    rppa.genes[i] <- rppa.to.gene2$Genes[match(rownames(RPPA)[i], rppa.to.gene2$Protein)]
  }
}

# Finally, we have manually created a mapping for the few antibodies/proteins that were not mapped to a gene.
# Map final proteins
rppa.to.gene3 <- read.table("annotation_data/proteinproposed.txt", header = T, stringsAsFactors = F)
for (i in 1:length(rppa.genes)) {
  if (rppa.genes[i] == "" | is.na(rppa.genes[i])) {
    rppa.genes[i] <- rppa.to.gene3$GENE[match(rownames(RPPA)[i], rppa.to.gene3$ID)]
  }
}

# Download TCGA BRCA RNA-seq data from from UCSC Xenabrowser
utils::download.file("https://gdc.xenahubs.net/download/TCGA-BRCA.htseq_fpkm-uq.tsv.gz", "R/RNA_seq.gz")
RNAseq_download <- utils::read.table(gzfile("R/RNA_seq.gz"))
ens.ids <- as.character(RNAseq_download[2:nrow(RNAseq_download), 1])
# Make data readable
names.rna <- unlist(lapply(RNAseq_download, FUN = function(l) as.character(l[1])))[2:ncol(RNAseq_download)]
names.rna <- gsub("-", ".", names.rna)
names(names.rna) <- NULL

# Fill in values in matrix
RNAseq <- matrix(0, ncol = ncol(RNAseq_download) - 1, nrow = nrow(RNAseq_download) - 1)
for (i in 1:ncol(RNAseq)) {
  RNAseq[, i] <- as.numeric(as.character(RNAseq_download[[i + 1]][2:nrow(RNAseq_download)]))
}
rownames(RNAseq) <- ens.ids
colnames(RNAseq) <- names.rna
# Order columns
RNAseq <- RNAseq[, order(colnames(RNAseq))]
dim(RNAseq) # 60483 x 1217


# Download gene annotation for RNAseq from UCSC Xenabrowser: https://xenabrowser.net/datapages/?dataset=TCGA-BRCA.htseq_fpkm.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
utils::download.file("https://gdc.xenahubs.net/download/gencode.v22.annotation.gene.probeMap", "R/gencode.v22.annotation.gene.probeMap")
geneAnnot <- read.table("R/gencode.v22.annotation.gene.probeMap", header = T, sep = "\t")
dim(geneAnnot)

# Check that the annotation corresponds to the names in the RNA-seq data
mean(rownames(RNAseq) %in% geneAnnot$id) # full match
mean(geneAnnot$id %in% rownames(RNAseq)) # full match

# Get it in the same order
geneAnnot <- geneAnnot[match(rownames(RNAseq), geneAnnot$id), c("id", "gene")]

# include genesymbol names in  RNAseq:
RNAseq_withGeneSymbol <- cbind(geneAnnot, RNAseq)



# Select only the unique genes the proteins in RPPA have been mapped to

# First check that all the genes in rppa.genes are present in RNAseq_withGeneSymbol and vice versa.
length(unique(rppa.genes)) # 182
sum(RNAseq_withGeneSymbol$gene %in% rppa.genes) #  182

# Only use the RNA-seq data from the genes present in the RPPA data set.
RNAseq_matchedWithRPPA <- RNAseq_withGeneSymbol[RNAseq_withGeneSymbol$gene %in% rppa.genes, ]
length(unique(RNAseq_matchedWithRPPA$gene)) #  182 unique genes


# Match patients in RNA-seq data with patients in RPPA data set
dim(RNAseq_matchedWithRPPA) # 182 unique genes
RNA.seq <- RNAseq_matchedWithRPPA[, 3:ncol(RNAseq_matchedWithRPPA)] # Skip first two columns, as they contain the ens IDs and gene names
rownames(RNA.seq) <- RNAseq_matchedWithRPPA[, "gene"] # Rownames are now the genes encoding the proteins
dim(RNA.seq) # RNA.seq is now 182 x1217.

# Remove last letter in patient IDs, as it only indicates sample number.
patients <- rep(NA, length(colnames(RNA.seq)))
for (i in 1:length(patients)) {
  patients[i] <- substr(colnames(RNA.seq)[i], 1, nchar(colnames(RNA.seq)[i]) - 1)
}


# Check how many patients that are in both the RPPA and RNA.seq data frame
sum(colnames(RPPA) %in% patients) # 880
sum(patients %in% colnames(RPPA)) # 885
# Not a one-to-one correspondence of patients.

# Get the IDs of the patients with several samples
colnames(RNA.seq)[which(patients %in% patients[duplicated(patients)])]

# Only consider patients in RNA.seq who are also in the RPPA dataset. Disregard second or third samples of patients in RNA.seq.
inboth.indices.rna <- which(patients %in% colnames(RPPA) & !duplicated(patients)) # Use the first sample (A) of each patient
RNA.seq.inboth <- RNA.seq[, inboth.indices.rna]
# Use the patient IDs without the last letter (A/B/C) as there now only is one sample per patient
colnames(RNA.seq.inboth) <- patients[inboth.indices.rna]

# Only patients in RPPA who are also in RNA.seq:
RPPA.inboth <- RPPA[, which(colnames(RPPA) %in% patients)]
dim(RPPA.inboth) # Now 281x880 (but only 181 unique proteins (and hence genes), as this is how many unique proteins the 281 antibodies map to)
# Corresponding associated genes are still in rppa.genes


# We now reduce the proteins we choose to look at

# First, we look at total proteins, instead of phosphorilated ones (which is only a subgroup)
is.total.protein <- rep(0, length(rownames(RPPA.inboth)))
for (i in 1:length(is.total.protein)) {
  # Check if substrings indicating phosphorilation (_pT, _pY or _pS) are present:
  is.total.protein[i] <- !(grepl("_pT", rownames(RPPA.inboth)[i]) | grepl("_pY", rownames(RPPA.inboth)[i]) | grepl("_pS", rownames(RPPA.inboth)[i]))
}
# Only consider total proteins
RPPA.total.proteins <- RPPA.inboth[which(is.total.protein == T), ]

# We then merge the mean value across different animal types used for the probing, if several different were used for one antibody
# Last two letters indicate animal type. We use mean across all animals.
# Remove ending in protein (antibody) names, to see which ones are actually the same
protein.names <- rep(NA, length(rownames(RPPA.total.proteins)))
for (i in 1:length(protein.names)) {
  protein.names[i] <- substr(rownames(RPPA.total.proteins)[i], 1, nchar(rownames(RPPA.total.proteins)[i]) - 4)
}

# Merge proteins across animal types:
several.proteins <- which(duplicated(protein.names)) # Indices of duplicated proteins (the index of the second repetition of a name)
protein.of.interest <- unique(protein.names[several.proteins]) # The unqiue protein names
RPPA.reduced <- RPPA.total.proteins[-several.proteins, ] # New data set, with each name only repeated once
rownames(RPPA.reduced) <- protein.names[-several.proteins] # New rownames, without the animal type-indicating endings
for (i in 1:length(protein.of.interest)) {
  equals.index <- which(protein.names == protein.of.interest[i]) # Indices of duplicates
  # Use mean value across animal types, add to data frame
  RPPA.reduced[which(rownames(RPPA.reduced) == protein.of.interest[i]), ] <- apply(RPPA.total.proteins[equals.index, ], 2, mean, na.rm = T)
}
dim(RPPA.reduced) # Now RPPA.reduced has dimension 170x880

# Create a data frame for mapping proteins to the genes they are encoded by
rppa.gene.mapping <- rppa.genes[which(is.total.protein == T)] # Not the phosphorylated proteins
rppa.gene.mapping <- rppa.gene.mapping[which(!duplicated(protein.names))] # New mapping of current genes corresponding to rownames of RPPA.reduced
mapping.frame <- cbind(rppa.gene.mapping, rownames(RPPA.reduced)) # One-to-one mapping between the genes and the proteins they encode. 170x2.
colnames(mapping.frame) <- c("Gene", "Protein")

# Reduce RNAseq data set down to the genes encoding the 170 proteins we now have
RNA.reduced <- RNA.seq.inboth[which(rownames(RNA.seq.inboth) %in% rppa.gene.mapping), ]
dim(RNA.reduced) # 169 x 880. Note: This actually correponds to 168 genes, as two pairs of proteins map to the same gene.


# Finally, deal with missingness.

# Remove patients that have several of the proteins still missing. Turns out 137 patients miss the same 6 proteins. We remove those proteins.
gene.with.missing <- which(is.na(RPPA.reduced), arr.ind = T)[1] # ID of one of the genes corresponding to a protein with many missing values
patients.with.missing.rppa <- which(is.na(RPPA.reduced[gene.with.missing, ])) # Index in RPPA data set of patients to remove.
patients.with.missing.rna <- which(colnames(RNA.reduced) %in% colnames(RPPA.reduced)[patients.with.missing.rppa]) # Index in RNA data set of patients to remove.
RPPA.reduced <- RPPA.reduced[, -patients.with.missing.rppa] # Now 170x743
RNA.reduced <- RNA.reduced[, -patients.with.missing.rna] # Now 169x743

# Now we remove the 4 proteins where more that half of the patients do not have measurements.
# Start by counting the number of missing measurements for each gene
no.missing <- rep(0, length(rownames(RPPA.reduced)))
for (i in 1:length(rownames(RPPA.reduced))) {
  no.missing[i] <- sum(is.na(RPPA.reduced[i, ]))
}
# We remove the proteins with missing values
proteins.to.remove.index <- which(no.missing != 0)
potential.genes.to.remove.index <- which(mapping.frame[, 2] %in% rownames(RPPA.reduced)[proteins.to.remove.index])
potential.genes.to.remove <- mapping.frame[potential.genes.to.remove.index, 1]
# Do not remove genes if there are other proteins without missing values map to them
genes.to.remove <- potential.genes.to.remove[-which(potential.genes.to.remove %in% mapping.frame[-potential.genes.to.remove.index, 1])]
genes.to.remove.index <- which(rownames(RNA.reduced) %in% genes.to.remove)
# Remove the genes corresponding to the missing proteins in both datasets.
RPPA.reduced <- RPPA.reduced[-proteins.to.remove.index, ] # Now 166x743
RNA.reduced <- RNA.reduced[-genes.to.remove.index, ] # Now 166x743

# Last, remove the gene 'XBP1', and the protein with the same name as it only has zeroes in the RNA-seq dataset.
RPPA.reduced <- RPPA.reduced[-which(rownames(RPPA.reduced) == "XBP1"), ] # Now 165 proteins, 743 patients
RNA.reduced <- RNA.reduced[-which(rownames(RNA.reduced) == "XBP1"), ] # Now 165 genes, 743 patients

# Update the mapping frame telling us which gene encodes which protein to only include those present in the datasets.
new.mapping.indices <- match(rownames(RPPA.reduced), mapping.frame[, 2])
mapping.frame <- mapping.frame[new.mapping.indices, ]


# Scale and transpose, and order genes and patients
rppa.scaled <- scale(t(RPPA.reduced))
rna.scaled <- scale(t(RNA.reduced))
# Rename column names to Gene name for simplicity, and order patients and genes.
colnames(rppa.scaled) <- mapping.frame[match(colnames(rppa.scaled), mapping.frame[, 2]), 1] # Rename
rppa.scaled <- rppa.scaled[, order(colnames(rppa.scaled))] # Order by Gene name
rna.scaled <- rna.scaled[, order(colnames(rna.scaled))] # Order by Gene name
RPPA.tcga <- rppa.scaled[order(rownames(rppa.scaled)), ] # Order patients
RNA.tcga <- rna.scaled[order(rownames(rna.scaled)), ] # Order patients
dim(RPPA.tcga) # 743 x 165
dim(RNA.tcga) # 743 x 165
# Now RPPA.tcga and RNA.tcga are final


# Section: analyse data -----------------------------------

n.tcga <- dim(RNA.tcga)[1]
p.tcga <- dim(RNA.tcga)[2]

# Create a prior weight matrix using the mRNA data
set.seed(1234)
rna.fit <- huge::huge(RNA.tcga, nlambda = 35, method = "glasso")
# Use STARS to construct a precision matrix for RNA data, which will be used to construct the weight matrix
set.seed(12345)
rna.select.stars <- huge::huge.select(rna.fit, criterion = "stars", stars.thresh = 0.05)
a.mat.rna.stars <- igraph::graph.adjacency(rna.select.stars$refit, mode = "undirected", diag = F) # Optimal graph
theta.rna.stars <- rna.select.stars$opt.icov # Optimal precision matrix
# Create the prior weight matrix
prior.mat.rna <- as.matrix(abs(stats::cov2cor(as.matrix(Matrix::forceSymmetric(rna.select.stars$opt.icov, uplo = "U"))) + 0)) # The prior corr matrix

# Use the ordinary graphical lasso on the RPPA data
set.seed(1234)
rppa.fit <- huge::huge(RPPA.tcga, nlambda = 35, method = "glasso")
# Select lambda by StARS
set.seed(12345)
rppa.select.stars <- huge::huge.select(rppa.fit, criterion = "stars", stars.thresh = 0.05)
lambda.rppa.stars <- rppa.select.stars$opt.lambda # Optimal lambda
a.mat.rppa.stars <- igraph::raph.adjacency(rppa.select.stars$refit, mode = "undirected", diag = F) # Optimal graph
theta.rppa.stars <- rppa.select.stars$opt.icov # Optimal precision matrix


# Perform the ordinary graphical lasso without penalised diagonal on RPPA data, to get results comparable to the other methods.
# Use the lambda selected by StARS above
rppa.tcga.glasso <- glasso::glasso(cov(RPPA.tcga), rho = rppa.select.stars$opt.lambda * matrix(1, p.tcga, p.tcga), penalize.diagonal = FALSE)
# Find sparsity of graphical lasso graph with penalty parameter selected by StARS.
# This is the sparsity we will force the estimates of the other methods to
tailoredGlasso::sparsity(rppa.tcga.glasso$wi != 0)
# 0.03370288

# Perform the ordinary weighted graphical lasso
pen.mat.tcga <- matrix(1 - prior.mat.rna, nrow = p.tcga, byrow = T) # the penalty matrix for the weighted graphical lasso
# As in the tailored graphical lasso, we preserve total amount of penalization
lambda.wglasso <- (rppa.select.stars$opt.lambda * p.tcga^2) / (sum(pen.mat.tcga))
fit.w.tcga <- glasso::glasso(cov(RPPA.tcga), rho = lambda.wglasso * pen.mat.tcga, penalize.diagonal = FALSE)
tailoredGlasso::sparsity(fit.w.tcga$wi != 0) # 0.03355506


# Perform the tailored graphical lasso. We already found the common lambda for the unweighted graph above with StARS, do no need to reselect it
set.seed(12345)
res.tcga <- tailoredGlasso::tailoredGlasso(x = RPPA.tcga, prior.matrix = prior.mat.rna, n = n.tcga, ebic.gamma = 0.6, lambda.opt = rppa.select.stars$opt.lambda, k.max = 150)
# Optimal sparsity
res.tcga$opt.sparsity
# optimal k
res.tcga$k.opt # 127


# Get wglasso and tailoredGlasso estimates to the same sparsity level as the graphical lasso estimate, so we can compare the log likelihoods directly
# For logistic method:
res.tcga.final <- tailoredGlasso::tailoredGlasso(RPPA.tcga, lambda.opt = 1.225 * rppa.select.stars$opt.lambda, prior.matrix = prior.mat.rna, n = n.tcga, k = res.tcga$k.opt)
res.tcga.final$opt.sparsity # 0.03370288


# Force wglasso estimate to the same sparsity levels
# Mofidy lambda to get the same sparsity level
lambda.wglasso.new <- (0.996 * rppa.select.stars$opt.lambda * p.tcga^2) / (sum(pen.mat.tcga))
fit.w.tcga.final <- glasso::glasso(cov(RPPA.tcga), rho = lambda.wglasso.new * pen.mat.tcga, penalize.diagonal = FALSE)
tailoredGlasso::sparsity(fit.w.tcga.final$wi != 0) # 0.03370288


# Finally: compare the log likelihoods
tailoredGlasso::gaussianloglik(cov(RPPA.tcga), fit.w.tcga.final$wi, n.tcga) # weighted graphical lasso
tailoredGlasso::gaussianloglik(cov(RPPA.tcga), rppa.tcga.glasso$wi, n.tcga) # graphical lasso
tailoredGlasso::gaussianloglik(cov(RPPA.tcga), res.tcga.final$theta.opt, n.tcga) # tailored graphical lasso
# Output:
# > tailoredGlasso::gaussianloglik(cov(RPPA.tcga),fit.w.tcga.final$wi,n.tcga) # weighted graphical lasso
# [1] -161928.4
# > tailoredGlasso::gaussianloglik(cov(RPPA.tcga),rppa.tcga.glasso$wi,n.tcga) # graphical lasso
# [1] -162226.8
# tailoredGlasso::gaussianloglik(cov(RPPA.tcga),res.tcga.final$theta.opt,n.tcga) # tailored graphical lasso
# [1] -160697.1


# The tailored graphical lasso estimate has the highest log likelihood



# Section: create histogram of data --------------------------------


# Only show the nonzero off-diagonal elements of the prior weight matrix
nonzero.weights.rna <- prior.mat.rna[!diag(rep(1, ncol(RNA.tcga)))]
nonzero.weights.rna <- nonzero.weights.rna[which(nonzero.weights.rna != 0, arr.ind = T)]
df.tcga <- data.frame(PartialCor = nonzero.weights.rna)
ggplot2::ggplot(df.tcga, aes(x = PartialCor)) +
  labs(title = "(c)") +
  ylab("Frequency") +
  xlab("Partial correlations") +
  geom_histogram(color = "darkolivegreen3", fill = "darkolivegreen3", breaks = seq(0, 0.27, by = 0.006)) +
  theme(legend.position = "none")


# Section: plot graph ------------------------------------

# Plot graphs
graph.tcga <- igraph::graph.adjacency(res.tcga.final$theta.opt != 0, mode = "undirected", diag = F)
V(graph.tcga)$name <- colnames(RPPA.tcga)
net.tcga <- network::network(igraph::as_adj(graph.tcga))
network.vertex.names(net.tcga) <- V(graph.tcga)$name
set.seed(123456)
GGally::ggnet2(net.tcga,
  node.size = 10, label.size = 2, edge.size = 0.3, node.label = V(graph.tcga)$name, alpha = 0.6,
  mode = "fruchtermanreingold", color = "maroon2"
)


# Section: check results with STRING database --------------------------------------------------------

# Function for making a data frame with the edges that are in one graph (mat1), but not the other (mat2)
make.df.of.changed.edges <- function(mat1, mat2, colnam) {
  changes.ind <- which(mat1 != 0 & mat2 == 0, arr.ind = T)
  changes <- data.frame(as.character(colnam[changes.ind[, 1]]), as.character(colnam[changes.ind[, 2]]))
  df <- data.frame(t(apply(changes, 1, sort)))
  df <- unique(df)
  rownames(df) <- 1:nrow(df)
  colnames(df) <- c("Gene1", "Gene2")
  return(df)
}


# Which edges are in tailoredGlasso but not in wglasso rppa networks for TCGA? Make data frame of them
df.changed.rppa.tcga <- make.df.of.changed.edges(res.tcga.final$theta.opt, fit.w.tcga.final$wi, colnames(RPPA.tcga))
list.unique.rppa.tcga <- unique(c(as.character(df.changed.rppa.tcga$Gene1), as.character(df.changed.rppa.tcga$Gene2)))
# Print a list of genes adjacent to any of these edges that can be inserted into STRING
lapply(list.unique.rppa.tcga, function(m) cat(m, "\n")) # Print results than can be checked in STRING
# In STRING: we enter the list into into STRING and get this graph: https://version-11-0b.string-db.org/cgi/network?networkId=bDw7GLqfzIFM
# A list of the edges in the STRING PPI network can be downloaded as a tsv:
STRING.output.tcga <- "https://string-db.org/cgi/generatetaskspecificdownloadfile?taskId=bDw7GLqfzIFM&downloadDataFormat=tsv&downloadFileName=string_interactions.tsv"
# STRING tells us:
# STRING DATABASE" PPI enrichment p-value: < 1.0e-16, your network has significantly more interactions than expected"


# Which edges are in wglasso rppa network but not in tailoredGlasso for TCGA:
df.lost.rppa.tcga <- make.df.of.changed.edges(fit.w.tcga.final$wi, res.tcga.final$theta.opt, colnames(RPPA.tcga))
list.unique.lost.rppa.tcga <- unique(c(as.character(df.lost.rppa.tcga$Gene1), as.character(df.lost.rppa.tcga$Gene2)))
# Print a list of genes that can be inserted into STRING
lapply(list.unique.lost.rppa.tcga, function(m) cat(m, "\n")) # Print results than can be checked in the STRING
# In STRING: we enter the list into into STRING and get this graph:
# A list of the edges in the STRING PPI network can be downloaded as a tsv: https://version-11-0b.string-db.org/cgi/network?networkId=bxBKZPl260eV
STRING.output.tcga.lost <- "https://string-db.org/cgi/generatetaskspecificdownloadfile?taskId=bxBKZPl260eV&downloadDataFormat=tsv&downloadFileName=string_interactions.tsv"
# STRING tells us:
# STRING DATABASE" PPI enrichment p-value: < 1.0e-16, your network has significantly more interactions than expected"



# Check agreement with STRING results for tailored graphical lasso
STRING.rppa.tcga <- utils::read.csv(STRING.output.tcga, sep = "\t")
# Sort so that gene pairs are given in alphabetical order.
STRING.rppa.tcga[, 1:2] <- t(apply(STRING.rppa.tcga[, 1:2], 1, sort)) # First gene in alphabet is always in col 1
colnames(STRING.rppa.tcga)[1:2] <- c("Gene1", "Gene2")
# Note that gene pairs in df.changed.rppa already are ordered this way.
# Combine edges in STRING and in tailoredGlasso into one table, so that we can check how many edges occur twice in the
#     table - this means they are present in both the STRING and the tailoredGlasso graph
df.all.edges.tcga <- rbind(df.changed.rppa.tcga, STRING.rppa.tcga[, 1:2])
ids.in.STRING <- which(duplicated(df.all.edges.tcga)) - nrow(df.changed.rppa.tcga) # Trues the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Get the mean number of edges that tailoredGlasso found but not wglasso, which are present in the STRING network.
sum(duplicated(df.all.edges.tcga)) / nrow(df.changed.rppa.tcga)
#  0.3217993

# New data frame with information about edges with evidence in STRING
df.with.evidence.tcga <- STRING.rppa.tcga[ids.in.STRING, ]
# Find average score of different evidence types
colMeans(df.with.evidence.tcga[, 5:13])
# neighborhood_on_chromosome       gene_fusion                              phylogenetic_cooccurrence           homology
# 0.0004516129                     0.0000000000                             0.0104946237                        0.0608494624
# coexpression                     experimentally_determined_interaction    database_annotated                  automated_textmining
# 0.1856021505                     0.2106559140                             0.2174193548                        0.5742903226
# combined_score
# 0.7301075269


# Check agreement with STRING results for weighted graphical lasso
STRING.rppa.tcga.removed <- utils::read.csv(STRING.output.tcga.lost, sep = "\t")
# Sort so that gene pairs are given in alphabetical order. Note that gene pairs in df.changed.rppa already are ordered this way.
STRING.rppa.tcga.removed[, 1:2] <- t(apply(STRING.rppa.tcga.removed[, 1:2], 1, sort)) # First gene in alphabet is always in col 1
colnames(STRING.rppa.tcga.removed)[1:2] <- c("Gene1", "Gene2")
# Combine edges in STRING and in wglasso into one table, so that we can check how many edges occur twice in the
#     table - this means they are present in both the STRING and the wglasso graph
df.all.edges.tcga.removed <- rbind(df.lost.rppa.tcga, STRING.rppa.tcga.removed[, 1:2])
ids.in.STRING.removed <- which(duplicated(df.all.edges.tcga.removed)) - nrow(df.lost.rppa.tcga) # Trues the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Get the mean number of edges that tailoredGlasso found but not wglasso, which are present in the STRING network.
sum(duplicated(df.all.edges.tcga.removed)) / nrow(df.lost.rppa.tcga)
# 0.1799308

# New data frame with information about edges with evidence in STRING
df.with.evidence.tcga.removed <- STRING.rppa.tcga.removed[ids.in.STRING.removed, ]
# Find average score of different evidence types
colMeans(df.with.evidence.tcga.removed[, 5:13])
# neighborhood_on_chromosome               gene_fusion                                 phylogenetic_cooccurrence             homology
# 0.001557692                              0.000000000                                 0.010461538                           0.042153846
# coexpression                             experimentally_determined_interaction       database_annotated                    automated_textmining
# 0.078846154                              0.235596154                                 0.321153846                           0.457326923
# combined_score
# 0.740230769

# Write list of edges to file: -------------------------------------------

# We now write a list of all edges indentified by the tailored graphical lasso to a .csv file.
# We also write a list of all edges that the tailored graphical lasso was able to find, but not wglasso, to file

get_and_print_edges <- function(a.mat, col.names) {
  # Function for printing all the edges in a graph with a layout that can be inserted into a latex table.
  # Also returns a data frame containing the edges
  # a.mat:          the adjacency matrix
  # col.names:      the names of the nodes in the graph
  a.mat[which(diag(rep(1, ncol(a.mat))) == 1, arr.ind = T)] <- 0
  pairs <- which(a.mat[, ] == 1, arr.ind = T)
  df <- data.frame(t(apply(pairs, 1, sort))) # Sort so that the node in the pair whose name is first in the alphabet is first.
  df <- unique(df)
  names <- cbind(col.names[df[, 1]], col.names[df[, 2]])
  for (i in 1:nrow(names)) {
    cat(names[i, 1], " & ", names[i, 2], " \\\\ \n")
  }
  return(names)
}


# Print all the edges that the tailored graphical lasso found, and write to file
edges.tcga <- get_and_print_edges(res.tcga.final$theta.opt != 0, colnames(RPPA.tcga))
colnames(edges.tcga) <- c("Gene1", "Gene2")
write.csv(edges.tcga, file = "Edge_lists/edgesTCGA", row.names = F)

# Write the list of edges that the tailored graphical lasso was able to find, but not wglasso, to file
write.csv(df.changed.rppa.tcga, file = "Edge_lists/edgesTCGA_unique", row.names = F)
