tailoredGlasso_simulation_extended = function(N,frac.to.mutate,true.cor,prior.cor,n,p, include.tailored=F, stars.thresh=0.05,ebic.gamma=0, ebic.gamma.space=0, 
                                              ebic.gamma.espace=0, alpha.genenet=0.2, alpha.ns = 0.2, alpha.cmi2ni = 0.03, lambda.min.espace=20,seed=NULL){
  # Simulate N scale-free graphs with multivariate Gaussian data, and N sets of prior graphs and prior Gaussian graphical data,
  #       and reconstruct graphs using the tailored graphical lasso, the graphical lasso and the weighted graphical lasso. 
  #
  # N:                  the number of simulations to perform 
  # frac.to.mutate:     the fraction of nodes to swap when constructing the prior graph
  # true.cor:           the size of the partial correlations of the data of the 'true' graph
  # prior.cor:          the size of the partial correlations of the data of the prior graph
  # n:                  the number of observations to generate in each simulated data set
  # p:                  the number of variables/nodes
  # include.tailored:   should the tailored graphical lasso (and the ordinary unweighted and weighted graphical lasso) be included? Logical
  # stars.thresh:       the variability threshold to use in StARS when selecting the sparsity for the unweighted graph
  # ebic.gamma:         the value of gamma to use in the eBIC selection of the optimal model in tailoredGlasso
  # ebic.gamma.space:   the the value of gamma to use in the eBIC selection of the optimal model in Space
  # ebic.gamma.espace:  the the value of gamma to use in the eBIC selection of the optimal model in Espace 
  # alpha.genenet:      the the value of the FDR-controllingalpha to use in GeneNet
  # alpha.ns:           the the value of the FDR-controlling alpha to use in the neighbourhood selection method
  # alpha.cmi2ni:       the the value of alpha to use in CMI2NI
  # seed:               the seed to use in the simulations, to ensure reproducible results
  # 
  # Make list for storing results for each simulation. 
  res = list(space.sparsities=rep(0,N), space.precisions=rep(0,N), space.recalls=rep(0,N), espace.sparsities=rep(0,N), espace.precisions=rep(0,N), 
             espace.recalls=rep(0,N), genenet.sparsities=rep(0,N), genenet.precisions=rep(0,N), genenet.recalls=rep(0,N), ns.sparsities=rep(0,N), 
             ns.precisions=rep(0,N), ns.recalls=rep(0,N), cmi2ni.sparsities=rep(0,N), cmi2ni.precisions=rep(0,N), cmi2ni.recalls=rep(0,N))
  if(include.tailored){
    res = c(res, list(tailored.sparsities=rep(0,N),k.opts=rep(0,N),tailored.precisions=rep(0,N),tailored.recalls=rep(0,N),
               w.sparsities=rep(0,N),w.precisions=rep(0,N),w.recalls=rep(0,N),
               glasso.sparsities=rep(0,N),glasso.precisions=rep(0,N),glasso.recalls=rep(0,N),
               w0_ruleofthumbs=rep(0,N)))
    # The quantile used to determine the sigmoid midpoint is equal to stars.thresh
    res$w0.quantile=stars.thresh 
  }
  # Generate 'true' graph. Choose v so that the partial correlations are of the desired size. 
  if(true.cor==0.2) v=0.5
  else v=0.03
  set.seed(12345)
  graph = huge::huge.generator(n=n, d=p,graph = 'scale-free',v=v,u=0.05,verbose = F) 
  # True precision matrix
  theta.true = graph$omega 
  # Omit almost-zero elements, as huge.generator often has rounding errors 
  theta.true[which(abs(theta.true)<10^(-4),arr.ind = T)] = 0 
  # True adjacency matrix
  adj.true = as.matrix(graph$theta!=0)+0 
  # The true covariance matrix to sample from
  cov.mat.true = stats::cov2cor(graph$sigma) 
  # Generate prior graph. Choose v so that the partial correlations are of the desired size. 
  if(prior.cor==0.2) v.prior=0.5
  else v.prior=0.03
  # Use same seed as for 'true' graph to get the same graph structure
  set.seed(12345) 
  data.prior.tmp = huge::huge.generator(n=n, d=p,graph = 'scale-free',v=v.prior,u=0.05,verbose = F) 
  # Create a prior graph structure with the desired fraction of nodes swapped
  prior.graph = mutate.graph(data.prior.tmp,frac.to.mutate)
  # The covariance matrix we will sample prior data from
  cov.mat.prior = prior.graph$cov.mat 
  # Set seed, if provided
  if(!is.null(seed)) set.seed(seed)
  # Sample N new true and prior data sets, use prior data to create a prior weight matrix, and use graph reconstruction methods. 
  for(i in 1:N){
    # Create an nxp matrix of simulated data. Scale for fair penalisation
    X = scale(MASS::mvrnorm(n,mu=rep(0,p),Sigma = cov.mat.true)) 
    X.prior = scale(MASS::mvrnorm(n,mu=rep(0,p),Sigma = cov.mat.prior))
    X.cov.mat = stats::cov(X)
    
    # Create prior weight matrix
    fit.huge.prior = huge::huge(X.prior,method="glasso",nlambda=35,verbose = F) 
    fit.stars.prior = huge::huge.select(fit.huge.prior,criterion='stars',stars.thresh,verbose = F) 
    prior.theta.est = stats::cov2cor(as.matrix(Matrix::forceSymmetric(fit.stars.prior$opt.icov,uplo='U')))+0
    prior.mat = as.matrix(abs(prior.theta.est))
    
    if(include.tailored){
      # Select optimal common lambda for the unweighted graph with StARS.
      fit.huge = huge::huge(X,method="glasso",nlambda=35,verbose = F) 
      fit.stars= huge::huge.select(fit.huge,criterion='stars',stars.thresh,verbose = F) 
      lambda.opt=fit.stars$opt.lambda
  
      # Tailored graphical lasso 
      prior.mat.off.diag = prior.mat[!diag(rep(1,p))]
      # Select w0 according to our rule of thumb. 
      w0_ruleofthumb = stats::quantile(prior.mat.off.diag[prior.mat.off.diag!=0],stars.thresh) # Times two bc all edges in prec matrix are duplicated
      res$w0_ruleofthumbs[i] = w0_ruleofthumb
      # Use the tailored graphical lasso with w0 and lambda.opt set
      tuned.w.ruleofthumb = tailoredGlasso::tailoredGlasso(lambda.opt=lambda.opt, x=X.cov.mat, prior.matrix=prior.mat,n=n,ebic.gamma=ebic.gamma,k.max=81,w0=w0_ruleofthumb,verbose = F)
      # Save sparsity, chosen steepness parameter k, precision and recall of estimate
      res$tailored.sparsities[i] = tuned.w.ruleofthumb$opt.sparsity
      res$k.opts[i] = tuned.w.ruleofthumb$k.opt
      res$tailored.precisions[i] = tailoredGlasso::precision(adj.true,tuned.w.ruleofthumb$theta.opt!=0)
      res$tailored.recalls[i] = tailoredGlasso::recall(adj.true,tuned.w.ruleofthumb$theta.opt!=0)
      
      # Ordinary weighted graphical lasso
      # The weight matrix for the weighted graphical lasso 
      weights = matrix(1-prior.mat,nrow=p,byrow=T)
      # Preserve total amount of penalization
      lambda.wglasso = (lambda.opt*p^2)/(sum(weights)) 
      fit.w = glasso::glasso(X.cov.mat, rho=lambda.wglasso*weights,penalize.diagonal = FALSE)
      res$w.sparsities[i] = tailoredGlasso::sparsity(fit.w$wi!=0) 
      res$w.precisions[i] = tailoredGlasso::precision(adj.true,fit.w$wi!=0)
      res$w.recalls[i]  = tailoredGlasso::recall(adj.true,fit.w$wi!=0)
      
      # Ordinary graphical lasso
      # Use the lambda selected by StARS
      lambda.glasso = lambda.opt 
      fit.glasso = glasso::glasso(X.cov.mat, rho=lambda.glasso*matrix(1,p,p),penalize.diagonal = FALSE)
      res$glasso.sparsities[i] = tailoredGlasso::sparsity(fit.glasso$wi!=0) 
      res$glasso.precisions[i] = tailoredGlasso::precision(adj.true,fit.glasso$wi!=0)
      res$glasso.recalls[i]  = tailoredGlasso::recall(adj.true,fit.glasso$wi!=0)
    }
    
    # Space
    fit.space = space.select.eBIC(X, ebic.gamma = ebic.gamma.space)
    res$space.sparsities[i] = fit.space$spars.opt
    res$space.precisions[i] = tailoredGlasso::precision(adj.true, fit.space$corr.mat!=0) 
    res$space.recalls[i] = tailoredGlasso::recall(adj.true, fit.space$corr.mat!=0) 
    
    # Espace
    fit.espace = espace.select.eBIC(X, graph.prior = prior.theta.est, ebic.gamma= ebic.gamma.espace, lambda.min=lambda.min.espace)
    res$espace.sparsities[i] = fit.espace$spars.opt
    res$espace.precisions[i] = tailoredGlasso::precision(adj.true, fit.espace$corr.mat!=0) 
    res$espace.recalls[i] = tailoredGlasso::recall(adj.true, fit.espace$corr.mat!=0) 
    
    # GeneNet
    fit.genenet = get_GeneNet(X, alpha=alpha.genenet)
    res$genenet.sparsities[i] = fit.genenet$sparsity
    res$genenet.precisions[i] = tailoredGlasso::precision(adj.true, fit.genenet$theta!=0) 
    res$genenet.recalls[i] = tailoredGlasso::recall(adj.true, fit.genenet$theta!=0) 
    
    # Neighbourhood selection by Meinshausen & Buhlmann
    # Select lambda according to their proposed rule-of-thumb to control the FDR
    lambdas.ns = sqrt(colSums(X^2))*(-qnorm(alpha.ns/(2*p^2)))/n # One lambda per node, but because we have standardized variables they are equal
    ns.list = huge::huge(X, method='mb', lambda=lambdas.ns[1], sym='and') # Because the lambda is equal for all vars, we provide the first value only
    res$ns.sparsities[i] = ns.list$sparsity
    res$ns.precisions[i] = tailoredGlasso::precision(adj.true, as.matrix(ns.list$path[[1]])!=0) 
    res$ns.recalls[i] = tailoredGlasso::recall(adj.true, as.matrix(ns.list$path[[1]])!=0) 
    
    # CMI2NI
    res.cmi2ni = CMI2NI(X,lamda=alpha.cmi2ni)
    res$cmi2ni.sparsities[i] = tailoredGlasso::sparsity(res.cmi2ni$G!=0)
    res$cmi2ni.precisions[i] = tailoredGlasso::precision(adj.true, res.cmi2ni$G!=0) 
    res$cmi2ni.recalls[i] = tailoredGlasso::recall(adj.true, res.cmi2ni$G!=0) 
    
    # To track progress, inform every time the completion percentage is dividable by 5. 
    done <- round(100 * i / N)
    if(done %% 5 == 0) {
      cat(done, " % done \n")
    }
  }
  return(res)
}


space.select.eBIC = function(X, ebic.gamma=0, lambda.min=20, lambda.max=100){
  # Select lambda by ebic, using 20 as the minimum value as proposed in space paper
  n = nrow(X)
  lambda.vals = seq(lambda.min,lambda.max, length.out = 100)
  ebic.vals = rep(NA,length(lambda.vals)) # rows for alpha vals, cols for lambda vals
  for(i in 1:length(lambda.vals)){
    res.space = space::space.joint(X, lambda.vals[i])
    # Find precision matrix from the partial correlation matrix
    res.space$ParCor[which(abs(res.space$ParCor) < 1e-6, arr.ind=T)] = 0 
    partial.cor = -res.space$ParCor
    diag(partial.cor) = 1 # Only off-diagonal elements get a change of sign
    icov.space = MBESS::cor2cov(partial.cor, sqrt(res.space$sig.fit)) # invert this to get the covariance matrix (sigma). w_d is the diagonal elements of the inverse cov matrix
    ebic.vals[i] = tailoredGlasso::eBIC(cov(X), icov.space, n, ebic.gamma)
    # If sparsity zero is reached, terminate
    if(sparsity(res.space$ParCor!=0) == 0){
      ebic.vals = ebic.vals[1:i]
      lambda.vals = lambda.vals[1:i]
      break()
    }
  }
  ind.opt = which.min(ebic.vals)
  lambda.opt = lambda.vals[ind.opt[1]]
  mod.opt = space::space.joint(X,lambda.opt)
  mod.opt$ParCor[which(abs(mod.opt$ParCor)<1e-6)] = 0 
  spars.opt = sparsity(mod.opt$ParCor)
  return(list(mod.opt = mod.opt, corr.mat=mod.opt$ParCor, spars.opt = spars.opt,lambda.opt=lambda.opt))
} 

espace.select.eBIC = function(X, hubs.ind=NULL, graph.prior=NULL, ebic.gamma=0, alpha.min=0.001, alpha.max=1, lambda.min=20, lambda.max=100){
  # Select alpha and lambda by eBIC
  # hubs.ind is the index of the hubs. Alternatively, a prior graph graph.prior can be provided. If so, we identify hubs as in espace paper
  # Use a grid as this is proposed in the paper, despite computational limitations
  # lambda from 20 and up, as in ordinary space 
  n = nrow(X)
  alpha.vals = seq(alpha.min, alpha.max, length.out = 20)
  lambda.vals = seq(lambda.min,lambda.max, length.out = 50)
  ebic.vals = matrix(NA,length(alpha.vals), length(lambda.vals)) # rows for alpha vals, cols for lambda vals
  if(is.null(hubs.ind)){
    espace.degree = igraph::degree(igraph::graph.adjacency(graph.prior!=0, mode='undirected', diag=F))
    hubs.ind = which(espace.degree > 7 & espace.degree > quantile(espace.degree, 0.95))
  }
  done = F
  for(i in 1:length(lambda.vals)){
    for(j in 1:length(alpha.vals)){
      res.espace = espace::espace(X,hubs.ind, alpha=alpha.vals[j], lambda=lambda.vals[i])
      res.espace$rho[which(abs(res.espace$rho) < 1e-6, arr.ind=T)] = 0 
      # Find precision matrix from the partial correlation matrix
      partial.cor = -res.espace$rho
      diag(partial.cor) = 1 # Only off-diagonal elements get a change of sign
      icov.espace = MBESS::cor2cov(partial.cor, sqrt(res.espace$w_d)) # invert this to get the covariance matrix (sigma). w_d is the diagonal elements of the inverse cov matrix
      # If resulting estimate is not of full rank, estimate is not valid and more shrinkage is needed:
      if(!is.positive.definite(icov.espace)){
        ebic.vals[j,i] = Inf
      }
      else{
        ebic.vals[j,i] = tailoredGlasso::eBIC(cov(X), icov.espace, n, ebic.gamma)
      }
      # If sparsity zero for smallest alpha is reached, terminate
      if(j == 1 & sparsity(res.espace$rho!=0) == 0){
        ebic.vals = ebic.vals[,1:i]
        lambda.vals = lambda.vals[1:i]
        done = T
        break()
      }
      if(done) break()
    }
    if(done) break()
  }
  ind.opt = which(ebic.vals==min(ebic.vals), arr.ind=T)
  alpha.opt = alpha.vals[ind.opt[,1]]
  lambda.opt = lambda.vals[ind.opt[,2]]
  mod.opt = espace::espace(X,hubs.ind, alpha=alpha.opt, lambda=lambda.opt)
  spars.opt = sparsity(mod.opt$rho)
  return(list(mod.opt = mod.opt, corr.mat=mod.opt$rho, spars.opt = spars.opt, alpha.opt=alpha.opt, lambda.opt=lambda.opt, ebic.vals=ebic.vals))
}  


GeneNetToGraph = function(edge.list, p){
  # Convert data frame of values returned from extract.network to adjacency matrix
  g = matrix(rep(0,p^2),p)
  diag(g) = 1
  for (i in 1:nrow(edge.list)){
    g[edge.list[i,'node1'], edge.list[i,'node2']] = edge.list[i,'pcor']
    g[edge.list[i,'node2'], edge.list[i,'node1']] = edge.list[i,'pcor']
  }
  return(g)
}

get_GeneNet = function(X, alpha=0.2){
  genenet.est = GeneNet::ggm.estimate.pcor(X)
  res.genenet<- GeneNet::network.test.edges(genenet.est)
  # Extract network containing edges with prob > 1-alpha (i.e. local fdr < alpha)
  net.genenet <- GeneNet::extract.network(res.genenet, cutoff.ggm=1-alpha)
  theta.genenet = GeneNetToGraph(net.genenet,ncol(X))
  return(list(theta=theta.genenet, sparsity = tailoredGlasso::sparsity(theta.genenet!=0)))
}

print_paper_results_extended = function(obj,frac.mutated,partial.cor,prior.partial.cor){
  # Print simulation results from extended study that can be inserted directly into latex table
  # Print mean sparsities, precisions and recalls when several data sets were generated
  # Results from the tailored graphical lasso and the ordinary weighted and unweighted graphical lasso are printed separately with print_paper_results.
  # obj:                    the object returned by tailoredGlasso_simulation
  # frac.mutated:           the fraction of the nodes that ahd their edges swapped in the illustratory graph
  # partial.cor:            the size of the partial correlations in the 'true' graph
  # prior.partial.cor:      the size of the partial correlations in the prior graph
  # 
  # Find average results across simulations
  # Space
  space.spars = mean(obj$space.sparsities)
  space.prec = mean(obj$space.precisions)
  space.recall = mean(obj$space.recalls)
  # Espace
  espace.spars = mean(obj$espace.sparsities)
  espace.prec = mean(obj$espace.precisions)
  espace.recall = mean(obj$espace.recalls)
  # NS
  ns.spars = mean(obj$ns.sparsities)
  ns.prec = mean(obj$ns.precisions)
  ns.recall = mean(obj$ns.recalls)
  # GeneNet
  genenet.spars = mean(obj$genenet.sparsities)
  genenet.prec = mean(obj$genenet.precisions)
  genenet.recall = mean(obj$genenet.recalls)
  # CMI2NI
  cmi2ni.spars = mean(obj$cmi2ni.sparsities)
  cmi2ni.prec = mean(obj$cmi2ni.precisions)
  cmi2ni.recall = mean(obj$cmi2ni.recalls)
  
  # Print results that can be inserted into latex tabular enviroment
  cat(paste0('\\text{Edge disagreement} \\% & \\text{Partial cor} & \\text{Prior partial cor} & \\text{Method} & k_{\\text{opt}}, & \\text{Sparsity} & \\text{Precision} & \\text{Recall} \\\\ \n'))
  #cat(frac.mutated*100, '&', partial.cor,  '&', prior.partial.cor,  '&')
  if(round(space.spars,3)==0){
    cat(paste0('& & & & Space &-&'),  round(space.spars,3), ' & - & - \\\\ \n')
  }
  else{
    cat(paste0('& & & & Space &-&'),  round(space.spars,3),' & ', round(space.prec,3),'&',round(space.recall,3),' \\\\ \n')
  }
  if(round(espace.spars,3)==0){
    cat(paste0('& & & & Espace &-&'),  round(espace.spars,3),' & - & - \\\\ \n')
  }
  else{
    cat(paste0('& & & & Espace &-&'),  round(espace.spars,3),' & ', round(espace.prec,3),'&',round(espace.recall,3),' \\\\ \n')
  }
  if(round(ns.spars,3)==0){
    cat(paste0('& & & & NS &-&'),  round(ns.spars,3),' & - & - \\\\ \n')
  }
  else{
    cat(paste0('& & & & NS &-&'),  round(ns.spars,3),' & ', round(ns.prec,3),'&',round(ns.recall,3),' \\\\ \n')
  }
  if(round(genenet.spars,3)==0){
    cat(paste0('& & & & GeneNet &-&'),  round(genenet.spars,3),' & - & - \\\\ \n')
  }
  else{
    cat(paste0('& & & & GeneNet &-&'),  round(genenet.spars,3),' & ', round(genenet.prec,3),'&',round(genenet.recall,3),' \\\\ \n')
  }
  if(round(cmi2ni.spars,3)==0){
    cat(paste0('& & & & CMI2NI &-&'),  round(cmi2ni.spars,3), ' & - & - \\\\ \n')
  }
  else{
    cat(paste0('& & & & CMI2NI &-&'),  round(cmi2ni.spars,3),' & ', round(cmi2ni.prec,3),'&',round(cmi2ni.recall,3),' \\\\ \n')
  }
}

