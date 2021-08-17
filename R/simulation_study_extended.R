# Extended simulation study with more methods included. 
# Simulate realistic data.  7 different cases.

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
library(space)
library(GeneNet)
library(MBESS)
source('CMI2NI_R.r')
source('simulation_functions.R')
source('simulation_study_extended_functions.R')
devtools::build('espace') # Must build espace package from repository
devtools::install('espace')
library(espace)

# Section: Case 1 ----------------------------------
# No edges changed. Large cor (0.2) in true graph, same in prior. 
n=80
p=100
frac_to_mutate_sf = 0
true.cor= 0.2
prior.cor = 0.2
N=100
# Perform simulations
simulations.sf.realistic.priors.1.extended = tailoredGlasso_simulation_extended(N,frac.to.mutate=frac_to_mutate_sf,true.cor=true.cor,prior.cor=prior.cor,n,p,seed=1)
#save(simulations.sf.realistic.priors.1.extended, file='extended1.RData')
load('extended1.RData')
# Print only the results that are relevant for the paper
print_paper_results_extended(simulations.sf.realistic.priors.1.extended,frac_to_mutate_sf,partial.cor=true.cor,prior.partial.cor=prior.cor)
#& & & & Space &-& 0.028  &  0.365 & 0.471  \\ 
#& & & & Espace &-& 0.018  &  0.452 & 0.39  \\ 
#& & & & NS &-& 0.001  &  0.929 & 0.052  \\ 
#& & & & GeneNet &-& 0  & - & - \\ 
#& & & & CMI2NI &-& 0.04  &  0.236 & 0.472  \\ 

# Orig result: ebic.gamma=0.1, lambda.min = 20 for space, lambda.min=40 for espace. sparsity zero tested moved out of loop for alpha. alpha.min=0.001

# Section: Case 2 ----------------------------------
# No edges changed. Large cor (0.2) in true graph, smaller cor (0.1) in prior. 
frac_to_mutate_sf = 0
true.cor= 0.2
prior.cor = 0.1
N=100
# Perform simulations
simulations.sf.realistic.priors.2.extended = tailoredGlasso_simulation_extended(N,frac.to.mutate=frac_to_mutate_sf,true.cor=true.cor,
                                                                                prior.cor=prior.cor,n,p, lambda.min.espace = 30, n.lambda.espace=20)
save(simulations.sf.realistic.priors.2.extended, file='extended2.RData')
load('extended2.RData')
# Print only the results that are relevant for the paper
print_paper_results_extended(simulations.sf.realistic.priors.2.extended,frac_to_mutate_sf,partial.cor=true.cor,prior.partial.cor=prior.cor)
#& & & & Space &-& 0.029  &  0.361 & 0.486  \\ 
#& & & & Espace &-& 0.025  &  0.398 & 0.455  \\ 
#& & & & NS &-& 0.001  &  0.948 & 0.053  \\ 
#& & & & GeneNet &-& 0  & - & - \\ 
#& & & & CMI2NI &-& 0.04  &  0.24 & 0.484  \\ 

# Section: Case 3 ----------------------------------
# No edges changed. Small cor (0.1) in true graph, large cor (0.2) in prior. 
frac_to_mutate_sf = 0
true.cor= 0.1
prior.cor = 0.2
N=100
# Perform simulations
simulations.sf.realistic.priors.3.extended = tailoredGlasso_simulation_extended(N,frac.to.mutate=frac_to_mutate_sf,true.cor=true.cor,prior.cor=prior.cor,n,p)
#save(simulations.sf.realistic.priors.3.extended, file='extended3.RData')
load('extended3.RData')
# Print only the results that are relevant for the paper
print_paper_results_extended(simulations.sf.realistic.priors.3.extended,frac_to_mutate_sf,partial.cor=0.1,prior.partial.cor=0.2)
#& & & & Space &-& 0.001  &  0.585 & 0.007  \\ 
#& & & & Espace &-& 0.002  &  0.409 & 0.019  \\ 
#& & & & NS &-& 0  & - & - \\ 
#& & & & GeneNet &-& 0  & - & - \\ 
#& & & & CMI2NI &-& 0.031  &  0.07 & 0.107  \\ 

# Section: Case 4 ----------------------------------
# No edges changed. Small cor (0.1) in true graph, same in prior. 
frac_to_mutate_sf = 0
true.cor= 0.1
prior.cor = 0.1
N=100
# Perform simulations
simulations.sf.realistic.priors.4.extended = tailoredGlasso_simulation_extended(N,frac.to.mutate=frac_to_mutate_sf,true.cor=true.cor,prior.cor=prior.cor,n,p, 
                                                                                lambda.min.espace = 30,n.lambda.espace=20)
save(simulations.sf.realistic.priors.4.extended, file='extended4.RData')
load('extended4.RData')
# Print only the results that are relevant for the paper
print_paper_results_extended(simulations.sf.realistic.priors.4.extended,frac_to_mutate_sf,partial.cor=0.1,prior.partial.cor=0.1)
#& & & & Space &-& 0.001  &  0.585 & 0.007  \\ 
#& & & & Espace &-& 0.001  &  0.218 & 0.01  \\ 
#& & & & NS &-& 0  & - & - \\ 
#& & & & GeneNet &-& 0  & - & - \\ 
#& & & & CMI2NI &-& 0.031  &  0.07 & 0.107  \\ 


# Section: Case 5 ----------------------------------
# 10% of edges changed. Large cor (0.2) in true graph, same in prior. 
frac_to_mutate_sf = 0.1
true.cor= 0.2
prior.cor = 0.2
N=100
# Perform simulations
simulations.sf.realistic.priors.5.extended = tailoredGlasso_simulation_extended(N,frac.to.mutate=frac_to_mutate_sf,true.cor=true.cor,prior.cor=prior.cor,n,p,seed=1)
#save(simulations.sf.realistic.priors.5.extended, file='extended5.RData')
load('extended5.RData')
# Print only the results that are relevant for the paper
print_paper_results_extended(simulations.sf.realistic.priors.5.extended,frac_to_mutate_sf,partial.cor=0.2,prior.partial.cor=0.2)
#& & & & Space &-& 0.028  &  0.365 & 0.471  \\ 
#& & & & Espace &-& 0.018  &  0.452 & 0.394  \\ 
#& & & & NS &-& 0.001  &  0.929 & 0.052  \\ 
#& & & & GeneNet &-& 0  & - & - \\ 
#& & & & CMI2NI &-& 0.04  &  0.236 & 0.472  \\ 

# Section: Case 6 ----------------------------------
# 20% of edges changed. Large cor (0.2) in true graph, same in prior. 
frac_to_mutate_sf = 0.2
true.cor= 0.2
prior.cor = 0.2
N=100
# Perform simulations
simulations.sf.realistic.priors.6.extended = tailoredGlasso_simulation_extended(N,frac.to.mutate=frac_to_mutate_sf,true.cor=true.cor,prior.cor=prior.cor,n,p,seed=1)
#save(simulations.sf.realistic.priors.6.extended, file='extended6.RData')
load('extended6.RData')
# Print only the results that are relevant for the paper
print_paper_results_extended(simulations.sf.realistic.priors.6.extended,frac_to_mutate_sf,partial.cor=0.2,prior.partial.cor=0.2)
#& & & & Space &-& 0.028  &  0.365 & 0.471  \\ 
#& & & & Espace &-& 0.026  &  0.382 & 0.456  \\ 
#& & & & NS &-& 0.001  &  0.929 & 0.052  \\ 
#& & & & GeneNet &-& 0  & - & - \\ 
#& & & & CMI2NI &-& 0.04  &  0.236 & 0.472  \\ 

# Section: Case 7 ----------------------------------
# 100% of edges changed. Large cor (0.2) in true graph, same in prior. 
frac_to_mutate_sf = 1
true.cor= 0.2
prior.cor = 0.2
N=100
# Perform simulations
simulations.sf.realistic.priors.7.extended = tailoredGlasso_simulation_extended(N,frac.to.mutate=frac_to_mutate_sf,true.cor=true.cor,prior.cor=prior.cor,n,p,seed=1)
#save(simulations.sf.realistic.priors.7.extended, file='extended7.RData')
load('extended7.RData')
# Print only the results that are relevant for the paper
print_paper_results_extended(simulations.sf.realistic.priors.7.extended,frac_to_mutate_sf,partial.cor=0.2,prior.partial.cor=0.2)
#& & & & Space &-& 0.028  &  0.365 & 0.471  \\ 
#& & & & Espace &-& 0.025  &  0.384 & 0.456  \\ 
#& & & & NS &-& 0.001  &  0.929 & 0.052  \\ 
#& & & & GeneNet &-& 0  & - & - \\ 
#& & & & CMI2NI &-& 0.04  &  0.236 & 0.472  \\ 

