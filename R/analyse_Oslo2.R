# Download, preprocess and analyse the Oslo 2 data with the tailored graphical lasso

# Neccessary additions
remotes::install_github("Camiling/tailoredGlasso", auth_token = "c19521d5a82ba1a631ffa5de9a569924908c10e4", dependencies = TRUE, build_vignettes = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("GEOquery")
library(GEOquery)
library(Biobase)
library(tailoredGlasso)
library(huge)
library(glasso)
library(igraph)
library(MASS)
library(ggplot2)
library(gridExtra)
library(GGally)
library(network)
library(readxl)


# Section: download and preprocess data --------------------------------


# Load mRNA Data set using Bioconductor

tmp <- GEOquery::getGEO("GSE58212", GSEMatrix = TRUE, getGPL = TRUE)
if (length(tmp) > 1) {
  idx <- grep("GPL14550", attr(tmp, "names"))
} else {
  idx <- 1
}
tmp2 <- tmp[[idx]]
# Load gene names
genes.oslo2.rna <- Biobase::fData(tmp2)$GENE_SYMBOL
OSLO2.RNA.temp <- Biobase::exprs(tmp2) # pxn (42405x283) mRNA expression matrix. Row names are the probe names.
# Name mRNA measures after their associated genes.
rownames(OSLO2.RNA.temp) <- genes.oslo2.rna
# Check dimension
dim(OSLO2.RNA.temp) # 42405x283


# Load RPPA Data set from supplementary material to Oslo 2 paper

# Download supplementary material from supplementary material of Oslo 2 paper (https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-015-0135-5#Sec25)
RPPA.link <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-015-0135-5/MediaObjects/13073_2015_135_MOESM4_ESM.xlsx"
utils::download.file(RPPA.link, "Oslo2_all_data.xlsx")
# RPPA data is found in sheet Q
tmp3 <- readxl::read_xlsx("Oslo2_all_data.xlsx", sheet = "Q. Oslo2 RPPA data", skip = 1)
tmp3 <- as.data.frame(tmp3)
# We get a 102x285 data frame. First column name is Antibody, second Genes, rest are patient IDs.
dim(tmp3)
# Matrix of only the expression measurements
OSLO2.RPPA.temp <- tmp3[, 3:ncol(tmp3)]
# Let the rownames be the genes known to encode the protein. Dimension pxn (102x283)
rownames(OSLO2.RPPA.temp) <- tmp3[, 2]
colnames(OSLO2.RPPA.temp) <- gsub("-", ".", colnames(OSLO2.RPPA.temp)) # Use . instead of - in IDs


# Convert patient IDs of mRNA and RPPA data so they correspond to each other
# Need the mapping between the sample names and geo accession ID
mapping.temp <- Biobase::pData(tmp2)[, c(1, 2)]
# Create data frame with sample names in first column and geo accession ID in second
patient.mapping.oslo2 <- as.data.frame(cbind(mapping.temp[, 2], as.character(mapping.temp[, 1])))
patient.mapping.oslo2[, 2] <- sapply(patient.mapping.oslo2[, 2], function(s) substring(s, 13)) # Remove first sentence in names
patient.mapping.oslo2[, 2] <- gsub("-", ".", patient.mapping.oslo2[, 2]) # Use . instead of - in IDs
# Change patient IDs in mRNA data set to correspond to IDs in the RPPA data set
colnames(OSLO2.RNA.temp) <- patient.mapping.oslo2[match(colnames(OSLO2.RNA.temp), patient.mapping.oslo2[, 1]), 2]
# Now patients have same ID in both data sets.

# Check if all the gene names in the RPPA data set are in mRNA data set.
rppa.in.rna.oslo2 <- rownames(OSLO2.RPPA.temp) %in% rownames(OSLO2.RNA.temp)
mean(rppa.in.rna.oslo2) # Not all are. We remove the ones that are not in the mRNA data set.
rownames(OSLO2.RPPA.temp)[which(!rppa.in.rna.oslo2)] # There are two antibodies that map to several genes.
# We remove antibodies mapping to several genes:
OSLO2.RPPA.temp <- OSLO2.RPPA.temp[which(rppa.in.rna.oslo2), ]
dim(OSLO2.RPPA.temp) # New dimension 100x283

# Check if some RPPA measurements are from the same genes
unique.genes.in.rppa.oslo2 <- unique(rownames(OSLO2.RPPA.temp))
length(unique.genes.in.rppa.oslo2) # All are unique, so no issue

# Only choose genes in RNA present in RPPA as well
OSLO2.RNA.newtemp <- OSLO2.RNA.temp[which(rownames(OSLO2.RNA.temp) %in% rownames(OSLO2.RPPA.temp)), ]
dim(OSLO2.RNA.newtemp) # New dimension 155x283

# Check if some genes are duplicated in the RNA set
unique.genes.in.rna.oslo2 <- unique(rownames(OSLO2.RNA.newtemp))
length(unique.genes.in.rna.oslo2) # 100 are unique, which are not all.
# Get the genes that several probes map to
several.probes.rna <- unique(rownames(OSLO2.RNA.newtemp)[duplicated(rownames(OSLO2.RNA.newtemp))])
several.probes.rna
# To get only one measure per gene, we merge the expression measures of these
OSLO2.RNA.merged <- OSLO2.RNA.newtemp
for (i in 1:length(several.probes.rna)) {
  # Indices of all probes that map to the gene several.probes.rna[i]
  ind <- which(rownames(OSLO2.RNA.merged) == several.probes.rna[i])
  # Scale before merging
  scaled.rows <- t(apply(OSLO2.RNA.merged[ind, ], 1, scale))
  # Merge values
  merged.row <- t(apply(scaled.rows, 2, mean))
  # Put the merged values in the place of the probe occuring first in the data frame
  OSLO2.RNA.merged[ind[1], ] <- merged.row
  # Remove other rows from dataset
  OSLO2.RNA.merged <- OSLO2.RNA.merged[-ind[2:length(ind)], ]
}

# If some patients have two samples, remove the second
# First part of IDs tell us which patient it is, last part which sample.
patients.oslo2 <- substr(colnames(OSLO2.RNA.merged), 1, 10) # Get patient ID, not including sample number
# The patients with several samples
relapse.patients <- colnames(OSLO2.RNA.merged)[which(patients.oslo2 %in% patients.oslo2[duplicated(patients.oslo2)])]
# Remove the second sample of these patients. To get those to remove, sort by name and choose every other sample name (with an even index)
remove.patients <- sort(relapse.patients)[which(rep(c(0, 1), length(relapse.patients) / 2) == 1)]
to.remove.rna <- which(colnames(OSLO2.RNA.merged) %in% remove.patients)
OSLO2.RNA.unscaled <- OSLO2.RNA.merged[, -to.remove.rna] # Dimension 100x280
to.remove.rppa <- which(colnames(OSLO2.RPPA.temp) %in% remove.patients)
OSLO2.RPPA.unscaled <- OSLO2.RPPA.temp[, -to.remove.rppa] # Dimension 100x280


# Finally: scale variables
OSLO2.RNA.scaled <- t(apply(OSLO2.RNA.unscaled, 1, scale))
OSLO2.RPPA.scaled <- t(apply(OSLO2.RPPA.unscaled, 1, scale))
colnames(OSLO2.RNA.scaled) <- colnames(OSLO2.RNA.unscaled)
colnames(OSLO2.RPPA.scaled) <- colnames(OSLO2.RPPA.unscaled)
# Sort rows and columns so that the orders are the same for the RPPA and mRNA data sets
OSLO2.RNA <- OSLO2.RNA.scaled[, order(colnames(OSLO2.RNA.scaled))]
OSLO2.RNA <- OSLO2.RNA[order(rownames(OSLO2.RNA)), ]
OSLO2.RPPA <- OSLO2.RPPA.scaled[order(rownames(OSLO2.RPPA.scaled)), ]
# Transpose to make matrices nxp
OSLO2.RNA <- t(OSLO2.RNA) # 280x100
OSLO2.RPPA <- t(OSLO2.RPPA) # 280x100
n.oslo2 <- nrow(OSLO2.RPPA)
p.oslo2 <- ncol(OSLO2.RPPA)

# Now both data frames are 100x280. Save as RData objects:
save(OSLO2.RNA, file = "OSLO2RNA.RData")
save(OSLO2.RPPA, file = "OSLO2RPPA.RData")


# Section: analyse data --------------------------------------

# Create a prior weight matrix using the RNA data
set.seed(1234)
oslo2.rna.fit <- huge::huge(OSLO2.RNA, nlambda = 35, method = "glasso")
# Use STARS to construct a precision matrix for RNA data, which will be used to construct the weight matrix
set.seed(1235)
oslo2.rna.stars <- huge::huge.select(oslo2.rna.fit, criterion = "stars", stars.thresh = 0.05)
oslo2.rna.stars$opt.sparsity
a.mat.oslo2.rna.stars <- igraph::graph.adjacency(oslo2.rna.stars$path[oslo2.rna.stars$opt.index][[1]], mode = "undirected", diag = F)
# Create the prior weight matrix
prior.mat.rna.oslo2 <- as.matrix(abs(stats::cov2cor(as.matrix(Matrix::forceSymmetric(oslo2.rna.stars$opt.icov, uplo = "U"))) + 0)) # The prior corr matrix

# Use the ordinary graphical lasso on the RPPA data
set.seed(1234)
oslo2.rppa.fit <- huge::huge(OSLO2.RPPA, nlambda = 35, method = "glasso")
# Select lambda by StARS
oslo2.rppa.stars <- huge::huge.select(oslo2.rppa.fit, criterion = "stars", stars.thresh = 0.05)
oslo2.rppa.stars$opt.sparsity
a.mat.oslo2.rppa.stars <- igraph::graph.adjacency(oslo2.rppa.stars$path[oslo2.rppa.stars$opt.index][[1]], mode = "undirected", diag = F)

# Perform the ordinary graphical lasso without penalised diagonal, to get results comparable to the other methods.
# Use the lambda selected by StARS above
oslo2.glasso <- glasso::glasso(cov(OSLO2.RPPA), rho = oslo2.rppa.stars$opt.lambda * matrix(1, p.oslo2, p.oslo2), penalize.diagonal = FALSE)
# Find sparsity of graphical lasso graph with penalty parameter selected by StARS.
# This is the sparsity we will force the estimates of the other methods to
tailoredGlasso::sparsity(oslo2.glasso$wi != 0)
# 0.02646465

# Perform the ordinary weighted graphical lasso
pen.mat.oslo2 <- matrix(1 - prior.mat.rna.oslo2, nrow = p.oslo2, byrow = T) # the penalty matrix for the weighted graphical lasso
# As in the tailored graphical lasso, we preserve total amount of penalization
lambda.wglasso <- (oslo2.rppa.stars$opt.lambda * p.oslo2^2) / (sum(pen.mat.oslo2))
fit.w.oslo2 <- glasso::glasso(cov(OSLO2.RPPA), rho = lambda.wglasso * pen.mat.oslo2, penalize.diagonal = FALSE)
tailoredGlasso::sparsity(fit.w.oslo2$wi != 0)


# Perform the tailored graphical lasso. We already found the common lambda for the unweighted graph above with StARS, do no need to reselect it
set.seed(1234)
res.oslo2 <- tailoredGlasso::tailoredGlasso(x = OSLO2.RPPA, prior.matrix = prior.mat.rna.oslo2, n = n.oslo2, ebic.gamma = 0.6, lambda.opt = oslo2.rppa.stars$opt.lambda, k.max = 151)
# Optimal sparsity
res.oslo2$opt.sparsity # 0.04
# optimal k
res.oslo2$k.opt # 149


# Get wglasso and tailoredGlasso estimates to the same sparsity level as the graphical lasso estimate, so we can compare the log likelihood directly
# For logistic method:
res.oslo2.final <- tailoredGlasso::tailoredGlasso(OSLO2.RPPA,
  lambda.opt = 3.9 * oslo2.rppa.stars$opt.lambda, prior.matrix = prior.mat.rna.oslo2, n = n.oslo2, ebic.gamma = 0,
  k = res.oslo2$k.opt
)
res.oslo2.final$opt.sparsity # 0.02646465

# Force wglasso estimate to the same sparsity levels
# Mofidy lambda to get the same sparsity level
lambda.wglasso.new <- (0.99 * oslo2.rppa.stars$opt.lambda * p.oslo2^2) / (sum(pen.mat.oslo2))
fit.w.oslo2.final <- glasso::glasso(cov(OSLO2.RPPA), rho = lambda.wglasso.new * pen.mat.oslo2, penalize.diagonal = FALSE)
tailoredGlasso::sparsity(fit.w.oslo2.final$wi != 0) # 0.02646465


# Finally: compare the log likelihoods
tailoredGlasso::gaussianloglik(cov(OSLO2.RPPA), fit.w.oslo2.final$wi, n.oslo2) # weighted graphical lasso
tailoredGlasso::gaussianloglik(cov(OSLO2.RPPA), oslo2.glasso$wi, n.oslo2) # graphical lasso
tailoredGlasso::gaussianloglik(cov(OSLO2.RPPA), res.oslo2.final$theta.opt, n.oslo2) # tailored graphical lasso
# Output:
# > gaussianloglik(cov(OSLO2.RPPA),fit.w.oslo2.final$wi,n.oslo2) # weighted graphical lasso
# [1] -38419.44
# > gaussianloglik(cov(OSLO2.RPPA),oslo2.glasso$wi,n.oslo2) # graphical lasso
# [1] -38489.65
# > gaussianloglik(cov(OSLO2.RPPA),res.oslo2.final$theta.opt,n.oslo2) # tailored graphical lasso
# [1] -38132.04

# The tailored graphical lasso estimate has the highest log likelihood



# Section: plot resulting tailored graphical lasso graph -----------------------------------------------------

# Plot graph with ggnet2. Must then convert graph to network object
graph.oslo2 <- igraph::graph.adjacency(res.oslo2.final$theta.opt != 0, mode = "undirected", diag = F)
V(graph.oslo2)$name <- colnames(OSLO2.RPPA)
net.oslo2 <- network::network(igraph::as_adj(graph.oslo2))
# Ensure node names will be shown in the plot
network.vertex.names(net.oslo2) <- V(graph.oslo2)$name
# Plot with Fruchterman-Reingold layout for a nice visualisation
set.seed(123456)
GGally::ggnet2(net.oslo2,
  node.size = 10, label.size = 2, edge.size = 0.3, node.label = V(graph.oslo2)$name, alpha = 0.6,
  mode = "fruchtermanreingold", color = "maroon2"
)



# Section: create histogram of non-zero prior weights -----------------------------------------------------

# Only show the nonzero off-diagonal elements of the prior weight matrix
nonzero.weights.rppa.oslo2 <- prior.mat.rna.oslo2[!diag(rep(1, ncol(OSLO2.RNA)))]
df.oslo2 <- data.frame(PartialCor = nonzero.weights.rppa.oslo2[nonzero.weights.rppa.oslo2 != 0])
ggplot2::ggplot(df.oslo2, aes(x = PartialCor)) +
  labs(title = "(b)") +
  ylab("Frequency") +
  xlab("Partial correlations") +
  geom_histogram(color = "darkolivegreen3", fill = "darkolivegreen3", breaks = seq(0, 0.27, by = 0.006)) +
  theme(legend.position = "none")


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

# Which edges are in tailoredGlasso but not in wglasso rppa networks for OSLO2? Make data frame of them
df.changed.rppa.oslo2 <- make.df.of.changed.edges(res.oslo2.final$theta.opt, fit.w.oslo2.final$wi, colnames(OSLO2.RPPA))
list.unique.rppa.oslo2 <- unique(c(as.character(df.changed.rppa.oslo2$Gene1), as.character(df.changed.rppa.oslo2$Gene2)))
# Print a list of genes adjacent to any of these edges that can be inserted into STRING
lapply(list.unique.rppa.oslo2, function(m) cat(m, "\n"))
# In STRING: we enter the list into into STRING and get this graph: https://version-11-0b.string-db.org/cgi/network?networkId=bkagSA6agzsn
# A list of the edges in the STRING PPI network can be downloaded as a tsv:
STRING.output.oslo2 <- "https://string-db.org/cgi/generatetaskspecificdownloadfile?taskId=bkagSA6agzsn&downloadDataFormat=tsv&downloadFileName=string_interactions.tsv"
# STRING tells us:
# STRING DATABASE" PPI enrichment p-value: < 1.0e-16, your network has significantly more interactions than expected"


# Which edges are in wglasso rppa network but not in tailoredGlasso for OSLO2:
df.lost.rppa.oslo2 <- make.df.of.changed.edges(fit.w.oslo2.final$wi, res.oslo2.final$theta.opt, colnames(OSLO2.RPPA))
list.unique.lost.rppa.oslo2 <- unique(c(as.character(df.lost.rppa.oslo2$Gene1), as.character(df.lost.rppa.oslo2$Gene2)))
# Print a list of genes that can be inserted into STRING
lapply(list.unique.lost.rppa.oslo2, function(m) cat(m, "\n"))
# In STRING: we enter the list into into STRING and get this graph: https://version-11-0b.string-db.org/cgi/network?networkId=bH5g3qrRnmlA
# A list of the edges in the STRING PPI network can be downloaded as a tsv:
STRING.output.oslo2.lost <- "https://string-db.org/cgi/generatetaskspecificdownloadfile?taskId=bH5g3qrRnmlA&downloadDataFormat=tsv&downloadFileName=string_interactions.tsv"
# STRING tells us:
# STRING DATABASE" PPI enrichment p-value: < 1.0e-16, your network has significantly more interactions than expected"

# The below uses the tsv file that STRING gives as output, which we saved above.

# Check agreement with STRING results for tailored graphical lasso
STRING.rppa.oslo2 <- utils::read.csv(STRING.output.oslo2, sep = "\t")
# Sort so that gene pairs are given in alphabetical order.
STRING.rppa.oslo2[, 1:2] <- t(apply(STRING.rppa.oslo2[, 1:2], 1, sort)) # First gene in alphabet is always in col 1
colnames(STRING.rppa.oslo2)[1:2] <- c("Gene1", "Gene2")
# Note that gene pairs in df.changed.rppa already are ordered this way.
# Combine edges in STRING and in tailoredGlasso into one table, so that we can check how many edges occur twice in the
#     table - this means they are present in both the STRING and the tailoredGlasso graph
df.all.edges.oslo2 <- rbind(df.changed.rppa.oslo2, STRING.rppa.oslo2[, 1:2])
ids.in.STRING <- which(duplicated(df.all.edges.oslo2)) - nrow(df.changed.rppa.oslo2) # Trues the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Get the mean number of edges that tailoredGlasso found but not wglasso, which are present in the STRING network.
sum(duplicated(df.all.edges.oslo2)) / nrow(df.changed.rppa.oslo2)
# 0.4368932

# New data frame with information about edges with evidence in STRING
df.with.evidence.oslo2 <- STRING.rppa.oslo2[ids.in.STRING, ]
# Find average score of different evidence types
colMeans(df.with.evidence.oslo2[, 5:13])
# neighborhood_on_chromosome            gene_fusion                           phylogenetic_cooccurrence
# 0.002444444                           0.000000000                           0.018911111
# homology                              coexpression                          experimentally_determined_interaction
# 0.072333333                           0.291933333                           0.341200000
# database_annotated                    automated_textmining                  combined_score
# 0.303555556                           0.627333333                           0.801511111


# Check agreement with STRING results for weighted graphical lasso
STRING.rppa.oslo2.removed <- utils::read.csv(STRING.output.oslo2.lost, sep = "\t")
# Sort so that gene pairs are given in alphabetical order. Note that gene pairs in df.changed.rppa already are ordered this way.
STRING.rppa.oslo2.removed[, 1:2] <- t(apply(STRING.rppa.oslo2.removed[, 1:2], 1, sort)) # First gene in alphabet is always in col 1
colnames(STRING.rppa.oslo2.removed)[1:2] <- c("Gene1", "Gene2")
# Combine edges in STRING and in wglasso into one table, so that we can check how many edges occur twice in the
#     table - this means they are present in both the STRING and the wglasso graph
df.all.edges.oslo2.removed <- rbind(df.lost.rppa.oslo2, STRING.rppa.oslo2.removed[, 1:2])
ids.in.STRING.removed <- which(duplicated(df.all.edges.oslo2.removed)) - nrow(df.lost.rppa.oslo2) # Trues the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Get the mean number of edges that tailoredGlasso found but not wglasso, which are present in the STRING network.
sum(duplicated(df.all.edges.oslo2.removed)) / nrow(df.lost.rppa.oslo2)
#  0.2815534

# New data frame with information about edges with evidence in STRING
df.with.evidence.oslo2.removed <- STRING.rppa.oslo2.removed[ids.in.STRING.removed, ]
# Find average score of different evidence types
colMeans(df.with.evidence.oslo2.removed[, 5:13])
# neighborhood_on_chromosome            gene_fusion                           phylogenetic_cooccurrence
# 0.00000000                            0.00000000                            0.00000000
# homology                              coexpression                          experimentally_determined_interaction
# 0.07848276                            0.05175862                            0.22886207
# database_annotated                    automated_textmining                  combined_score
# 0.33103448                            0.53500000                            0.67006897

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
edges.oslo2 <- get_and_print_edges(res.oslo2.final$theta.opt != 0, colnames(OSLO2.RPPA))
colnames(edges.oslo2) <- c("Gene1", "Gene2")
write.csv(edges.oslo2, file = "Edge_lists/edgesOslo2", row.names = F)

# Write the list of edges that the tailored graphical lasso was able to find, but not wglasso, to file
write.csv(df.changed.rppa.oslo2, file = "Edge_lists/edgesOslo2_unique", row.names = F)
