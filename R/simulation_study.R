# Simulate realistic data.  7 different cases. 


remotes::install_github("Camiling/tailoredGlasso", auth_token='c19521d5a82ba1a631ffa5de9a569924908c10e4', dependencies = TRUE,build_vignettes = TRUE)
library(tailoredGlasso)

# Section: Generate data with scale-free property-----------------------------------

n=80
p=100
set.seed(12345)
data.scalefree = huge.generator(n=n, d=p,graph = 'scale-free',v=0.5,u=0.05) 
g.true.scalefree = data.scalefree$theta # True adjacency matrix
g.scalefree =graph.adjacency(data.scalefree$theta,mode="undirected",diag=F) # true igraph object
x.scalefree = data.scalefree$data # Observed attributes. nxp matrix.
x.scalefree.scaled= scale(x.scalefree) # Scale columns/variables.
s.scalefree.scaled = cov(x.scalefree.scaled) # Empirical covariance matrix
data.scalefree$sparsity # True sparsity

# Section: Use glasso on scale-free graph -----------------------------------

# Use data to select lambdas, estimate graph with huge. Note that huge scales internally. 
fit.huge.scalefree = huge(x.scalefree,method="glasso",nlambda=35) # Look at 35 different lambda values


# Section: Case 1 ----------------------------------
# No edges changed. Large cor (0.2) in true graph, same in prior. 
w0_vals_sf = c(0,0.01,0.1,0.25,0.5)
frac_to_mutate_sf = 0
v_sf = 0.5
N=100
simulations.sf.realistic.priors.1 = weighttuning.sf.with.realistic.priors(N,data.scalefree,frac.to.mutate=frac_to_mutate_sf,v=v_sf,w0.quantile=w0_vals_sf,n,p)
print.results.sf.realistic.priors(simulations.sf.realistic.priors.1,frac_to_mutate_sf,partial.cor=0.2,prior.partial.cor=0.2)


# Section: Case 2 ----------------------------------
# No edges changed. Large cor (0.2) in true graph, smaller cor (0.1) in prior. 
w0_vals_sf = c(0,0.01,0.1,0.25,0.5)
frac_to_mutate_sf = 0
v_sf = 0.03
N=100
simulations.sf.realistic.priors.2 = weighttuning.sf.with.realistic.priors(N,data.scalefree,frac.to.mutate=frac_to_mutate_sf,v=v_sf,w0.quantile=w0_vals_sf,n,p)
print.results.sf.realistic.priors(simulations.sf.realistic.priors.2,frac_to_mutate_sf,partial.cor=0.2,prior.partial.cor=0.1)

# Section: Case 3 ----------------------------------
# No edges changed. Small cor (0.1) in true graph, large cor (0.2) in prior. 
# Must first create a true graph with same edges, but smaller partial correlations
set.seed(12345) # Same seed as before
data.scalefree.0.1 = huge.generator(n=n, d=p,graph = 'scale-free',v=0.03,u=0.05) # Omega is the prec matrix
w0_vals_sf = c(0,0.01,0.1,0.25,0.5)
frac_to_mutate_sf = 0
v_sf = 0.5
N=100
simulations.sf.realistic.priors.3 = weighttuning.sf.with.realistic.priors(N,data.scalefree.0.1,frac.to.mutate=frac_to_mutate_sf,v=v_sf,w0.quantile=w0_vals_sf,n,p)
print.results.sf.realistic.priors(simulations.sf.realistic.priors.3,frac_to_mutate_sf,partial.cor=0.1,prior.partial.cor=0.2)

# Section: Case 4 ----------------------------------
# No edges changed. Small cor (0.1) in true graph, same in prior. 
w0_vals_sf = c(0,0.01,0.1,0.25,0.5)
frac_to_mutate_sf = 0
v_sf = 0.03
N=100
simulations.sf.realistic.priors.4 = weighttuning.sf.with.realistic.priors(N,data.scalefree.0.1,frac.to.mutate=frac_to_mutate_sf,v=v_sf,w0.quantile=w0_vals_sf,n,p)
print.results.sf.realistic.priors(simulations.sf.realistic.priors.4,frac_to_mutate_sf,partial.cor=0.1,prior.partial.cor=0.1)

# Section: Case 5 ----------------------------------
# 10% of edges changed. Large cor (0.2) in true graph, same in prior. 
w0_vals_sf = c(0,0.01,0.1,0.25,0.5)
frac_to_mutate_sf = 0.1
v_sf = 0.5
N=100
simulations.sf.realistic.priors.5 = weighttuning.sf.with.realistic.priors(N,data.scalefree,frac.to.mutate=frac_to_mutate_sf,v=v_sf,w0.quantile=w0_vals_sf,n,p)
print.results.sf.realistic.priors(simulations.sf.realistic.priors.5,frac_to_mutate_sf,partial.cor=0.2,prior.partial.cor=0.2)

# Section: Case 6 ----------------------------------
# 20% of edges changed. Large cor (0.2) in true graph, same in prior. 
w0_vals_sf = c(0,0.01,0.1,0.25,0.5)
frac_to_mutate_sf = 0.2
v_sf = 0.5
N=100
simulations.sf.realistic.priors.6 = weighttuning.sf.with.realistic.priors(N,data.scalefree,frac.to.mutate=frac_to_mutate_sf,v=v_sf,w0.quantile=w0_vals_sf,n,p)
print.results.sf.realistic.priors(simulations.sf.realistic.priors.6,frac_to_mutate_sf,partial.cor=0.2,prior.partial.cor=0.2)

# Section: Case 7 ----------------------------------
# 100% of edges changed. Large cor (0.2) in true graph, same in prior. 
w0_vals_sf = c(0,0.01,0.1,0.25,0.5)
frac_to_mutate_sf = 1
v_sf = 0.5
N=100
simulations.sf.realistic.priors.7 = weighttuning.sf.with.realistic.priors(N,data.scalefree,frac.to.mutate=frac_to_mutate_sf,v=v_sf,w0.quantile=w0_vals_sf,n,p)
print.results.sf.realistic.priors(simulations.sf.realistic.priors.7,frac_to_mutate_sf,partial.cor=0.2,prior.partial.cor=0.2)

# Section: Create some histograms of the off-diagonal precision matrix elements -----------------------------

# Case 1: No changed edges, large corr. 
priormat.case1 = create.sample.sf.prior(frac.to.mutate=0,v=0.5,n,p)
priorcorrs.case1 = priormat.case1[!diag(rep(1,p))]
priorcorrs.case1 = priorcorrs.case1[which(priorcorrs.case1!=0,arr.ind = T)]

# Case 2: No changed edges, smaller corr. 
priormat.case2 = create.sample.sf.prior(frac.to.mutate=0,v=0.03,n,p)
priorcorrs.case2 = priormat.case2[!diag(rep(1,p))]
priorcorrs.case2 = priorcorrs.case2[which(priorcorrs.case2!=0,arr.ind = T)]

par(mfrow=c(1,2))
hist(priorcorrs.case1,breaks=100,xlab='Resulting prior weights',main = 'True prior partial correlations of 0.2')
hist(priorcorrs.case2,breaks=100,xlab='Resulting prior weights',main = 'True prior partial correlations of 0.1')


