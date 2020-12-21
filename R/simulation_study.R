# Simulate realistic data.  7 different cases. Also plot an example graph and histograms of the prior weights. 

# Neccessary additions
remotes::install_github("Camiling/tailoredGlasso", dependencies = TRUE,build_vignettes = TRUE)
library(tailoredGlasso)
library(huge)
library(glasso)
library(igraph)
library(MASS)
library(ggplot2)
library(gridExtra)
library(GGally)
library(network)

# The below code requires functions found in 'simulation_study.R'

# Section: Case 1 ----------------------------------
# No edges changed. Large cor (0.2) in true graph, same in prior. 
n=80
p=100
frac_to_mutate_sf = 0
true.cor= 0.2
prior.cor = 0.2
N=100
# Perform simulations
simulations.sf.realistic.priors.1 = tailoredGlasso_simulation(N,frac.to.mutate=frac_to_mutate_sf,true.cor=true.cor,prior.cor=prior.cor,n,p,seed=1)
# Print only the results that are relevant for the paper
print_paper_results(simulations.sf.realistic.priors.1,frac_to_mutate_sf,partial.cor=true.cor,prior.partial.cor=prior.cor)
#\text{Edge disagreement} \% & \text{Partial cor} & \text{Prior partial cor} & \text{Method} & k_{\text{opt}}, & \text{Sparsity} & \text{Precision} & \text{Recall} \\ 
#0 & 0.2 & 0.2 &Glasso &-&  0.035  &  0.283 & 0.493  \\ 
#& & & Wglasso &-& 0.032  &  0.312 & 0.503  \\ 
#& & & TailoredGlasso & 49.64  &  0.031  &  0.389 & 0.606 \\ 

# Section: Case 2 ----------------------------------
# No edges changed. Large cor (0.2) in true graph, smaller cor (0.1) in prior. 
frac_to_mutate_sf = 0
true.cor= 0.2
prior.cor = 0.1
N=100
# Perform simulations
simulations.sf.realistic.priors.2 = tailoredGlasso_simulation(N,frac.to.mutate=frac_to_mutate_sf,true.cor=true.cor,prior.cor=prior.cor,n,p)
# Print only the results that are relevant for the paper
print_paper_results(simulations.sf.realistic.priors.2,frac_to_mutate_sf,partial.cor=true.cor,prior.partial.cor=prior.cor)
#\text{Edge disagreement} \% & \text{Partial cor} & \text{Prior partial cor} & \text{Method} & k_{\text{opt}}, & \text{Sparsity} & \text{Precision} & \text{Recall} \\ 
#0 & 0.2 & 0.1 &Glasso &-&  0.035  &  0.285 & 0.499  \\ 
#& & & Wglasso &-& 0.034  &  0.293 & 0.493  \\ 
#& & & TailoredGlasso & 13.39  &  0.034  &  0.295 & 0.496 \\ 

# Section: Case 3 ----------------------------------
# No edges changed. Small cor (0.1) in true graph, large cor (0.2) in prior. 
frac_to_mutate_sf = 0
true.cor= 0.1
prior.cor = 0.2
N=100
# Perform simulations
simulations.sf.realistic.priors.3 = tailoredGlasso_simulation(N,frac.to.mutate=frac_to_mutate_sf,true.cor=true.cor,prior.cor=prior.cor,n,p)
# Print only the results that are relevant for the paper
print_paper_results(simulations.sf.realistic.priors.3,frac_to_mutate_sf,partial.cor=0.1,prior.partial.cor=0.2)
##\text{Edge disagreement} \% & \text{Partial cor} & \text{Prior partial cor} & \text{Method} & k_{\text{opt}}, & \text{Sparsity} & \text{Precision} & \text{Recall} \\ 
#0 & 0.1 & 0.2 &Glasso &-&  0.022  &  0.079 & 0.085  \\ 
#& & & Wglasso &-& 0.021  &  0.096 & 0.099  \\ 
#& & & TailoredGlasso & 3.34  &  0.021  &  0.1 & 0.103 \\ 


# Section: Case 4 ----------------------------------
# No edges changed. Small cor (0.1) in true graph, same in prior. 
frac_to_mutate_sf = 0
true.cor= 0.1
prior.cor = 0.1
N=100
# Perform simulations
simulations.sf.realistic.priors.4 = tailoredGlasso_simulation(N,frac.to.mutate=frac_to_mutate_sf,true.cor=true.cor,prior.cor=prior.cor,n,p)
# Print only the results that are relevant for the paper
print_paper_results(simulations.sf.realistic.priors.4,frac_to_mutate_sf,partial.cor=0.1,prior.partial.cor=0.1)
#\text{Edge disagreement} \% & \text{Partial cor} & \text{Prior partial cor} & \text{Method} & k_{\text{opt}}, & \text{Sparsity} & \text{Precision} & \text{Recall} \\ 
#0 & 0.1 & 0.1 &Glasso &-&  0.022  &  0.079 & 0.085  \\ 
#& & & Wglasso &-& 0.02  &  0.082 & 0.083  \\ 
#& & & TailoredGlasso & 5.63  &  0.02  &  0.083 & 0.084 \\ 


# Section: Case 5 ----------------------------------
# 10% of edges changed. Large cor (0.2) in true graph, same in prior. 
frac_to_mutate_sf = 0.1
true.cor= 0.2
prior.cor = 0.2
N=100
# Perform simulations
simulations.sf.realistic.priors.5 = tailoredGlasso_simulation(N,frac.to.mutate=frac_to_mutate_sf,true.cor=true.cor,prior.cor=prior.cor,n,p,seed=1)
# Print only the results that are relevant for the paper
print_paper_results(simulations.sf.realistic.priors.5,frac_to_mutate_sf,partial.cor=0.2,prior.partial.cor=0.2)
#\text{Edge disagreement} \% & \text{Partial cor} & \text{Prior partial cor} & \text{Method} & k_{\text{opt}}, & \text{Sparsity} & \text{Precision} & \text{Recall} \\ 
#10 & 0.2 & 0.2 &Glasso &-&  0.035  &  0.283 & 0.493  \\ 
#& & & Wglasso &-& 0.033  &  0.305 & 0.493  \\ 
#& & & TailoredGlasso & 41.2  &  0.032  &  0.335 & 0.532 \\ 

# Section: Case 6 ----------------------------------
# 20% of edges changed. Large cor (0.2) in true graph, same in prior. 
frac_to_mutate_sf = 0.2
true.cor= 0.2
prior.cor = 0.2
N=100
# Perform simulations
simulations.sf.realistic.priors.6 = tailoredGlasso_simulation(N,frac.to.mutate=frac_to_mutate_sf,true.cor=true.cor,prior.cor=prior.cor,n,p,seed=1)
# Print only the results that are relevant for the paper
print_paper_results(simulations.sf.realistic.priors.6,frac_to_mutate_sf,partial.cor=0.2,prior.partial.cor=0.2)
#\text{Edge disagreement} \% & \text{Partial cor} & \text{Prior partial cor} & \text{Method} & k_{\text{opt}}, & \text{Sparsity} & \text{Precision} & \text{Recall} \\ 
#20 & 0.2 & 0.2 &Glasso &-&  0.035  &  0.283 & 0.493  \\ 
#& & & Wglasso &-& 0.034  &  0.291 & 0.491  \\ 
#& & & TailoredGlasso & 4.31  &  0.034  &  0.291 & 0.493 \\ 

# Section: Case 7 ----------------------------------
# 100% of edges changed. Large cor (0.2) in true graph, same in prior. 
frac_to_mutate_sf = 1
true.cor= 0.2
prior.cor = 0.2
N=100
# Perform simulations
simulations.sf.realistic.priors.7 = tailoredGlasso_simulation(N,frac.to.mutate=frac_to_mutate_sf,true.cor=true.cor,prior.cor=prior.cor,n,p,seed=1)
# Print only the results that are relevant for the paper
print_paper_results(simulations.sf.realistic.priors.7,frac_to_mutate_sf,partial.cor=0.2,prior.partial.cor=0.2)
#\text{Edge disagreement} \% & \text{Partial cor} & \text{Prior partial cor} & \text{Method} & k_{\text{opt}}, & \text{Sparsity} & \text{Precision} & \text{Recall} \\ 
#100 & 0.2 & 0.2 &Glasso &-&  0.035  &  0.283 & 0.493  \\ 
#& & & Wglasso &-& 0.034  &  0.289 & 0.485  \\ 
#& & & TailoredGlasso & 2.94  &  0.034  &  0.29 & 0.485 \\ 

# Section: Create example histograms of the non-zero prior weights -----------------------------

# Case 1: No changed edges, large corr. 
priormat.case1 = create_sample_prior(frac.to.mutate=0,prior.cor=0.2,n,p)
priorcorrs.case1 = priormat.case1[!diag(rep(1,p))]
priorcorrs.case1 = priorcorrs.case1[which(priorcorrs.case1!=0,arr.ind = T)]
df.case1= data.frame(PartialCor=priorcorrs.case1)

# Case 2: No changed edges, smaller corr. 
priormat.case2 = create_sample_prior(frac.to.mutate=0,prior.cor=0.1,n,p) 
priorcorrs.case2 = priormat.case2[!diag(rep(1,p))]
priorcorrs.case2 = priorcorrs.case2[which(priorcorrs.case2!=0,arr.ind = T)]
df.case2 = data.frame(PartialCor=priorcorrs.case2)

# Create combined plot of prior weights
hist.1 = ggplot2::ggplot(df.case1,aes(x=PartialCor))+labs(title='(a)')+ylab('Frequency')+ xlab('Partial correlations')+
          geom_histogram(color='deepskyblue2',fill='deepskyblue2',breaks=seq(0,0.22,by=0.006))+theme(legend.position ="none")
hist.2 = ggplot2::ggplot(df.case2,aes(x=PartialCor))+labs(title='(b)')+ylab('Frequency')+  xlab('Partial correlations')+ 
          geom_histogram(color='deepskyblue2',fill='deepskyblue2',breaks=seq(0,0.22,by=0.006))+theme(legend.position ="none")
gridExtra::grid.arrange(hist.1,hist.2,ncol=2,nrow=1)



# Section: create illustratory graph -----------------------------------

n=80
p=100
# Use same seed as was used to create the graph in the simulations
set.seed(12345)
data.scalefree = huge::huge.generator(n=n, d=p,graph = 'scale-free',v=0.5,u=0.05) 
# Create graph object
g.scalefree =igraph::graph.adjacency(data.scalefree$theta,mode="undirected",diag=F) 

# Function for assigning colors to each vertex according to their degree. Returns vector of color codes. 
color.by.degree = function(adj.mat){
  # adj.mat:   the igraph adjacency matrix object
  # Assign luminance in [0,100], where a low node degree corresponds to high luminance
  luminance= 100-(igraph::degree(adj.mat)-min(igraph::degree(adj.mat)))/(max(igraph::degree(adj.mat))-min(igraph::degree(adj.mat)))*100 
  # All nodes are given a pink color, where a darker pink corresponds to higher degree
  return( hcl(c=80,l=luminance,h=0) ) 
}

# Plot graph
net.sf=network::network(igraph::as_adj(g.scalefree))
col.sf = color.by.degree(g.scalefree)
set.seed(123456)
GGally::ggnet2(net.sf,node.size = 10, edge.size = 0.3,alpha=0.9,mode = "fruchtermanreingold",color = col.sf)


