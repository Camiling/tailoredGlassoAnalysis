tailoredGlasso_simulation = function(N,frac.to.mutate,true.cor,prior.cor,n,p,stars.thresh=0.05,ebic.gamma=0, seed=NULL){
  # Simulate N scale-free graphs with multivariate Gaussian data, and N sets of prior graphs and prior Gaussian graphical data,
  #       and reconstruct graphs using the tailored graphical lasso, the graphical lasso and the weighted graphical lasso. 
  #
  # N:                  the number of simulations to perform 
  # frac.to.mutate:     the fraction of nodes to swap when constructing the prior graph
  # true.cor:           the size of the partial correlations of the data of the 'true' graph
  # prior.cor:          the size of the partial correlations of the data of the prior graph
  # n:                  the number of observations to generate in each simulated data set
  # p:                  the number of variables/nodes
  # stars.thresh:       the variability threshold to use in StARS when selecting the sparsity for the unweighted graph. 
  # ebic.gamma:         the value of gamma to use in the eBIC selection of the optimal model. 
  # 
  # Make list for storing results for each simulation. 
  res = list(tailored.sparsities=rep(0,N),k.opts=rep(0,N),tailored.precisions=rep(0,N),tailored.recalls=rep(0,N),
             w.sparsities=rep(0,N),w.precisions=rep(0,N),w.recalls=rep(0,N),
             glasso.sparsities=rep(0,N),glasso.precisions=rep(0,N),glasso.recalls=rep(0,N),
             w0_ruleofthumbs=rep(0,N))
  # Generate 'true' graph. Choose v so that the partial correlations are of the desired size. 
  if(true.cor==0.2) v=0.5
  else v=0.03
  set.seed(12345)
  graph = huge::huge.generator(n=n, d=p,graph = 'scale-free',v=v,u=0.05,verbose = F) 
  # The quantile used to determine the sigmoid midpoint is equal to stars.thresh
  res$w0.quantile=stars.thresh 
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
    # To avoid having to perform StARS twice, we select the optimal common lambda for the unweighted graph 
    #     outside the tailoredGlasso function. 
    #
    # Select optimal common lambda for the unweighted graph with StaRS.
    fit.huge = huge::huge(X,method="glasso",nlambda=35,verbose = F) 
    fit.stars= huge::huge.select(fit.huge,criterion='stars',stars.thresh,verbose = F) 
    lambda.opt=fit.stars$opt.lambda
    # Create prior weight matrix
    fit.huge.prior = huge::huge(X.prior,method="glasso",nlambda=35,verbose = F) 
    fit.stars.prior = huge::huge.select(fit.huge.prior,criterion='stars',stars.thresh,verbose = F) 
    prior.theta.est = stats::cov2cor(as.matrix(Matrix::forceSymmetric(fit.stars.prior$opt.icov,uplo='U')))+0
    prior.mat = as.matrix(abs(prior.theta.est))

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
    
    # To track progress, inform every time the completion percentage is dividable by 5. 
    done <- round(100 * i / N)
    if(done %% 5 == 0) {
      cat(done, " % done \n")
    }
  }
  return(res)
}

create_sample_prior = function(frac.to.mutate,prior.cor,n,p,stars.thresh=0.05){
  # Create a sample prior matrix. Used for illustratory purposes
  # frac.to.mutate:   the fraction of the nodes to swap edges for in the illustratory graph
  # prior.cor:        the size of the partial correlations in the prior graph
  # n:                the number of observations to generate
  # p:                the number of variables/nodes
  # stars.thresh:     the variability threshold to use in StARS when selecting the sparsity for the unweighted graph. 
  #
  # Use same seed as before to get the same graph structure
  set.seed(12345) 
  if(prior.cor==0.2) v.prior=0.5
  else v.prior=0.03
  # Generate data
  data.prior.tmp = huge::huge.generator(n=n, d=p,graph = 'scale-free',v=v.prior,u=0.05,verbose=F) # Omega is the prec matrix
  # Change desired amount of edges
  prior.graph = mutate.graph(data.prior.tmp,frac.to.mutate)
  cov.mat.prior = prior.graph$cov.mat
  # Scale prior data
  X.prior = scale(MASS::mvrnorm(n,mu=rep(0,p),Sigma = cov.mat.prior))
  # Create prior precision matrix like we did in the simulations
  fit.huge.prior = huge::huge(X.prior,method="glasso",nlambda=35,verbose=F) 
  fit.stars.prior = huge::huge.select(fit.huge.prior,criterion='stars',stars.thresh,verbose=F) 
  # Create prior weight matrix
  prior.theta.est = stats::cov2cor(as.matrix(Matrix::forceSymmetric(fit.stars.prior$opt.icov,uplo='U')))+0 
  prior.mat = as.matrix(abs(prior.theta.est))
  return(prior.mat)
}

print_paper_results = function(obj,frac.mutated,partial.cor,prior.partial.cor){
  # Print simulation results that can be inserted directly into latex table
  # Print mean sparsities, precisions and recalls when several data sets were generated
  # obj:                    the object returned by tailoredGlasso_simulation
  # frac.mutated:   the fraction of the nodes that ahd their edges swapped in the illustratory graph
  # partial.cor:            the size of the partial correlations in the 'true' graph
  # prior.partial.cor:      the size of the partial correlations in the prior graph
  # 
  # Find average results across simulations
  # Tailored graphical lasso
  tailored.spars = mean(obj$tailored.sparsities) 
  tailored.prec = mean(obj$tailored.precisions)
  tailored.recall = mean(obj$tailored.recalls)
  k.opt = mean(obj$k.opts)
  # Weighted graphical lasso
  w.spars = mean(obj$w.sparsities)
  w.recall = mean(obj$w.recalls)
  w.prec = mean(obj$w.precisions)
  # Graphical lasso
  g.spars = mean(obj$glasso.sparsities) 
  g.prec = mean(obj$glasso.precisions)
  g.recall = mean(obj$glasso.recalls)
  # Print results that can be inserted into latex tabular enviroment
  cat(paste0('\\text{Edge disagreement} \\% & \\text{Partial cor} & \\text{Prior partial cor} & \\text{Method} & k_{\\text{opt}}, & \\text{Sparsity} & \\text{Precision} & \\text{Recall} \\\\ \n'))
  cat(frac.mutated*100, '&', partial.cor,  '&', prior.partial.cor,  '&')
  cat('Glasso &-& ',  round(g.spars,3),' & ', round(g.prec,3),'&',round(g.recall,3),' \\\\ \n')
  cat(paste0('& & & Wglasso &-&'),  round(w.spars,3),' & ', round(w.prec,3),'&',round(w.recall,3),' \\\\ \n')
  # Print results for logistic weight tuning method with rule-of-thumb choice of w_0:
  cat(paste0('& & & TailoredGlasso'),'&',round(k.opt,2),' & ',  round(tailored.spars,3),' & ', round(tailored.prec,3),'&',round(tailored.recall,3),'\\\\ \n')
}

mutate.graph = function(graph,fraction){
  # Swap a certain fraction of nodes in order to change graph
  # graph:      the huge.generate() object to mutate
  # fraction:   the fraction of edges to change
  prec.mat = graph$omega
  # Avoid rounding errors, which can happen with huge.generator
  prec.mat[which(abs(prec.mat)<10^(-4),arr.ind=T)]=0
  cov.mat = graph$sigma
  adj.mat = as.matrix(graph$theta)+0
  data=graph$data
  p = ncol(graph$omega)
  # If no nodes are to be swapped, simply return graph as is
  if(fraction==0){
    ans = list()
    ans$cov.mat=cov.mat
    ans$prec.mat = prec.mat
    ans$adj.mat = adj.mat
    ans$data = data
    return(ans)
  }
  # We basically 'swap pairs of nodes' by switching their cols and rows
  # Indices of edge pairs
  edges = which(adj.mat==1,arr.ind=T) 
  # Avoid doubling it up
  edges=edges[1:(nrow(edges)/2),] 
  n.mutations = floor(nrow(edges)*fraction)
  # Sample nodes to give new edges
  nodes.add = sample(1:p,n.mutations) 
  # Sample nodes to remove edges from
  nodes.remove = edges[sample(1:nrow(edges),n.mutations),1] 
  # Must swap edges fpr adjacency matrix, precision matrix, covariance matrix and data matrix. 
  for(i in 1:n.mutations){
    tmp.prec=prec.mat
    tmp.adj = adj.mat
    tmp.dat = data
    tmp.cov.mat = cov.mat
    # swap precision matrix rows. Then cols, the order does not matter
    prec.mat[nodes.remove[i],] = tmp.prec[nodes.add[i],]
    prec.mat[nodes.add[i],] = tmp.prec[nodes.remove[i],]
    tmp.prec=prec.mat
    # Swap precision matrix cols
    prec.mat[,nodes.remove[i]] = tmp.prec[,nodes.add[i]]
    prec.mat[,nodes.add[i]] = tmp.prec[,nodes.remove[i]]
    # Swap adjacency matrix rows
    adj.mat[nodes.remove[i],] = tmp.adj[nodes.add[i],]
    adj.mat[nodes.add[i],] = tmp.adj[nodes.remove[i],]
    tmp.adj = adj.mat
    # Swap adjacency matrix cols
    adj.mat[,nodes.remove[i]] = tmp.adj[,nodes.add[i]]
    adj.mat[,nodes.add[i]] = tmp.adj[,nodes.remove[i]]
    # Swap covariance mat rows
    cov.mat[nodes.remove[i],] = tmp.cov.mat[nodes.add[i],]
    cov.mat[nodes.add[i],] = tmp.cov.mat[nodes.remove[i],]
    tmp.cov.mat = cov.mat
    # Swap covariance mat cols
    cov.mat[,nodes.remove[i]] = tmp.cov.mat[,nodes.add[i]]
    cov.mat[,nodes.add[i]] = tmp.cov.mat[,nodes.remove[i]]
    # Swap data matrix cols (only the cols since it is nxp)
    data[,nodes.remove[i]] = tmp.dat[,nodes.add[i]]
    data[,nodes.add[i]] = tmp.dat[,nodes.remove[i]]
  }
  ans = list()
  ans$cov.mat=cov.mat
  ans$prec.mat = prec.mat
  ans$adj.mat = adj.mat
  ans$data = data
  return(ans)
}

