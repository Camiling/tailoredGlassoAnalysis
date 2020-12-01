weighttuning.sf.with.realistic.priors = function(N,graph,frac.to.mutate,v,w0.quantile,n,p,stars.thresh=0.05){
  # Simulate N scale-free prior graphs
  # N is the number of simulations to perform 
  # graph is the huge.generate() object to mutate, fraction is the fraction of edges to change.
  # lambda.opt is the optimal lambda for 
  # frac.to.mutate is the fraction of nodes to swap
  # v is the number to be added to the off-diagonal elements of the precision matrix
  # w0.quantile is the (array of) lower quantiles to use to determine the sigmoid midpoint w0.  
  res = list(k.sparsities=rep(0,N),k.opts=rep(0,N),k.dists=rep(0,N),k.precisions=rep(0,N),k.recalls=rep(0,N),
             k.rot.sparsities=rep(0,N),k.rot.opts=rep(0,N),k.rot.dists=rep(0,N),k.rot.precisions=rep(0,N),k.rot.recalls=rep(0,N),
             w.sparsities=rep(0,N),w.dists=rep(0,N),w.precisions=rep(0,N),w.recalls=rep(0,N),
             glasso.sparsities=rep(0,N),glasso.dists=rep(0,N),glasso.precisions=rep(0,N),glasso.recalls=rep(0,N),
             estimated.prior.dists=rep(0,N),estimated.prior.precisions=rep(0,N),w0_ruleofthumbs=rep(0,N))
  if(length(w0.quantile)>1){ # If several w0 are to be considered, store results in Nxlength(w0.quantile) matrix. 
    res$k.sparsities=res$k.opts=res$k.dists=res$k.precisions=res$k.recalls=matrix(rep(0,N*length(w0.quantile)),nrow=N)
  }
  res$w0=w0.quantile
  theta.true = graph$omega # omit almost-zero elements
  theta.true[which(abs(theta.true)<10^(-4),arr.ind = T)] = 0
  adj.true = graph$theta+0 # True adjacency matrix. 
  # The true covariance matrix to sample from:
  cov.mat.true = cov2cor(graph$sigma) # This is actually already scaled for the scale-free data. 
  # Now we create the prior matrix
  set.seed(12345) # Same seed as before
  data.prior.tmp = huge.generator(n=n, d=p,graph = 'scale-free',v=v,u=0.05) # Omega is the prec matrix
  prior.graph = mutate.graph(data.prior.tmp,frac.to.mutate)
  cov.mat.prior = prior.graph$cov.mat # The matrix we will sample from. 
  # Now: sample N new true and prior data sets, estimate prior prec matrix, use methods. 
  for(i in 1:N){
    X = scale(mvrnorm(n,mu=rep(0,p),Sigma = cov.mat.true)) # An nxp matrix
    X.prior = scale(mvrnorm(n,mu=rep(0,p),Sigma = cov.mat.prior))
    X.cov.mat = cov(X)
    # Find optimal lambda and orginary glasso estimated precision matrix:
    fit.huge = huge(X,method="glasso",nlambda=35) 
    fit.stars= huge.select(fit.huge,criterion='stars',stars.thresh) 
    lambda.opt=fit.stars$opt.lambda
    # Create prior precision matrix:
    fit.huge.prior = huge(X.prior,method="glasso",nlambda=35) 
    fit.stars.prior = huge.select(fit.huge.prior,criterion='stars',stars.thresh) 
    prior.theta.est = cov2cor(forceSymmetric(fit.stars.prior$opt.icov,uplo='U'))+0 # The prior corr matrix
    prior.mat = as.matrix(abs(prior.theta.est))
    res$estimated.prior.dists[i] = matrix.distance.simple(prior.theta.est,fit.stars$opt.icov) # Both are sparse dgCMatrix
    res$estimated.prior.precisions[i] = precision(fit.stars$refit+0, (prior.theta.est!=0)+0)
    # Tune linear transformation 1-q*w (included since wglasso is included here as well for q=1)
    tuned.q = linear.weight.transformation.select(lambda.opt,cov.mat=X.cov.mat,prior.mat,n,ebic.gamma = 0,q.max=1)[[1]]
    
    # Rule-of-thumb value for logistic tuning
    prior.mat.off.diag = prior.mat[!diag(rep(1,p))]
    w0_ruleofthumb = quantile(prior.mat.off.diag[prior.mat.off.diag!=0],stars.thresh) # Times two bc all edges in prec matrix are duplicated
    res$w0_ruleofthumbs[i] = w0_ruleofthumb
    tuned.w.ruleofthumb = tailoredGlasso(lambda.opt=lambda.opt, x=X.cov.mat, prior.matrix=prior.mat,n=n,ebic.gamma=0,k.max=81,w0=w0_ruleofthumb)
    res$k.rot.sparsities[i] = tuned.w.ruleofthumb$opt.sparsity
    res$k.rot.opts[i] = tuned.w.ruleofthumb$k.opt
    res$k.rot.dists[i] = matrix.distance.simple(theta.true,tuned.w.ruleofthumb$theta.opt)
    res$k.rot.precisions[i] = precision(adj.true,tuned.w.ruleofthumb$theta.opt!=0)
    res$k.rot.recalls[i] = recall(adj.true,tuned.w.ruleofthumb$theta.opt!=0)
    # Logistic tuning
    if(length(w0.quantile)>1){
      for(j in 1:length(w0.quantile)){
        w0_val = quantile(prior.mat.off.diag[prior.mat.off.diag!=0],w0.quantile[j]) 
        tuned.w = tailoredGlasso(lambda.opt=lambda.opt, x=X.cov.mat, prior.matrix=prior.mat,n=n,ebic.gamma=0,k.max=81,w0=w0_val)
        res$k.sparsities[i,j] = tuned.w$opt.sparsity
        res$k.opts[i,j] = tuned.w$k.opt
        res$k.dists[i,j] = matrix.distance.simple(theta.true,tuned.w$theta.opt)
        res$k.precisions[i,j] = precision(adj.true,tuned.w$theta.opt!=0)
        res$k.recalls[i,j] = recall(adj.true,tuned.w$theta.opt!=0)
      }
    }
    else{
      w0_val = quantile(prior.mat.off.diag[prior.mat.off.diag!=0],w0.quantile) 
      tuned.w =tailoredGlasso(lambda.opt=lambda.opt, x=X.cov.mat, prior.matrix=prior.mat,n=n,ebic.gamma=0,k.max=81,w0=w0_val)
      # Results from logistic method
      res$k.sparsities[i] = tuned.w$opt.sparsity
      res$k.opts[i] = tuned.w$k.opt
      res$k.dists[i] = matrix.distance.simple(theta.true,tuned.w$theta.opt)
      res$k.precisions[i] = precision(adj.true,tuned.w$theta.opt!=0)
      res$k.recalls[i] = recall(adj.true,tuned.w$theta.opt!=0)
    }
    # Save results from orginary wglasso
    n_q = length(tuned.q$q)
    res$w.sparsities[i] = tuned.q$sparsity[n_q]
    theta.est.w= tuned.q$thetas[,,n_q]
    res$w.dists[i]= matrix.distance.simple(theta.true,theta.est.w)
    res$w.precisions[i] = precision(adj.true,theta.est.w!=0)
    res$w.recalls[i]  = recall(adj.true,theta.est.w!=0)
    # Results from ordinary glasso
    res$glasso.sparsities[i] = tuned.q$sparsity[1]
    theta.est.g = tuned.q$thetas[,,1]
    res$glasso.dists[i] = matrix.distance.simple(theta.true,theta.est.g)
    res$glasso.precisions[i] = precision(adj.true,theta.est.g !=0)
    res$glasso.recalls[i] = recall(adj.true,theta.est.g!=0)
    cat('Index: ', i, '\n ')
  }
  return(res)
}

# Used to create corr. priors for scale-free network 
create.sample.sf.prior = function(frac.to.mutate,v,n,p,stars.thresh=0.05){
  # Function for creating a sample prior matrix, for a given frac.to.mutate and v. Used to make example histogram of weights
  set.seed(12345) # Same seed as before
  data.prior.tmp = huge.generator(n=n, d=p,graph = 'scale-free',v=v,u=0.05) # Omega is the prec matrix
  prior.graph = mutate.graph(data.prior.tmp,frac.to.mutate)
  cov.mat.prior = prior.graph$cov.mat # The matrix we will sample from. 
  X.prior = scale(mvrnorm(n,mu=rep(0,p),Sigma = cov.mat.prior))
  # Create prior precision matrix:
  fit.huge.prior = huge(X.prior,method="glasso",nlambda=35) 
  fit.stars.prior = huge.select(fit.huge.prior,criterion='stars',stars.thresh) 
  prior.theta.est = cov2cor(forceSymmetric(fit.stars.prior$opt.icov,uplo='U'))+0 # The prior corr matrix
  prior.mat = as.matrix(abs(prior.theta.est))
  return(prior.mat)
}

print.results.sf.realistic.priors.final = function(obj,frac.mutated,partial.cor,prior.partial.cor){
  # obj is an object returned by weighttuning.sf.with.realistic.priors.final().
  # Function for print mean sparsities, precisions, recalls and matrix distances when several data sets were generated
  # obj is returned by perform.realistic.simulations(). 
  
  # Results from rule-of-thumb:
  k.rot.spars = mean(obj$k.rot.sparsities) # Logistic method, mean sparsity
  k.rot.prec = mean(obj$k.rot.precisions)
  k.rot.recall = mean(obj$k.rot.recalls)
  k.rot.dist=mean(obj$k.rot.dists)
  k.rot.opt = mean(obj$k.rot.opts)
  w.spars = mean(obj$w.sparsities) # Wglasso, mean sparsity
  w.recall = mean(obj$w.recalls)
  w.prec = mean(obj$w.precisions)
  w.dist=mean(obj$w.dists)
  g.spars = mean(obj$glasso.sparsities) # Glasso, mean sparsity
  g.prec = mean(obj$glasso.precisions)
  g.recall = mean(obj$glasso.recalls)
  g.dist=mean(obj$glasso.dists)
  
  # Print results that can be inserted into latex tabular enviroment
  cat('\text{Edge disagreement} \% & \text{Partial cor} & \text{Prior partial cor} & \text{Method} & k_{\text{opt}}, & \text{Sparsity} & \text{Precision} & \text{Recall}')
  cat(frac.mutated*100, '&', partial.cor,  '&', prior partial cor,  '&')
  cat('Glasso &-& ',  round(g.spars,3),' & ', round(g.prec,3),'&',round(g.recall,3),' \\\\ \n')
  cat('& & & Wglasso $\\boldmath{P}_1$ &-&',  round(w.spars,3),' & ', round(w.prec,3),'&',round(w.recall,3),' \\\\ \n')
  # Print results for logistic weight tuning method with rule-of-thumb choice of w_0:
  cat(paste0('& & & Wglasso $\\boldmath{P}_{k_{\\text{opt}}}$&'),'&',round(k.rot.opt,2),' & ',  round(k.rot.spars,3),' & ', round(k.rot.prec,3),'&',round(k.rot.recall,3),'\\\\ \n')
  
}

mutate.graph = function(graph,fraction){
  # graph is the huge.generate() object to mutate, fraction is the fraction of edges to change. 
  # We basically 'swap pairs of nodes' by switching their cols and rows. 
  prec.mat = graph$omega
  prec.mat[which(abs(prec.mat)<10^(-4),arr.ind=T)]=0
  cov.mat = graph$sigma
  adj.mat = graph$theta
  data=graph$data
  p = ncol(graph$omega)
  if(fraction==0){
    ans = list()
    ans$cov.mat=cov.mat
    ans$prec.mat = prec.mat
    ans$adj.mat = adj.mat
    ans$data = data
    return(ans)
  }
  edges = which(adj.mat==1,arr.ind=T) # Edge pairs.
  edges=edges[1:(nrow(edges)/2),] # Avoid doubling it up. Now dim n_edges x 2
  n.mutations = floor(nrow(edges)*fraction)
  nodes.add = sample(1:p,n.mutations) # id of nodes to give the edges
  nodes.remove = edges[sample(1:nrow(edges),n.mutations),1] # The nodes to 'leave out'
  for(i in 1:n.mutations){
    tmp.prec=prec.mat
    tmp.adj = adj.mat
    tmp.dat = data
    tmp.cov.mat = cov.mat
    # swap prec mat rows. Then cols, the order does not matter!
    prec.mat[nodes.remove[i],] = tmp.prec[nodes.add[i],]
    prec.mat[nodes.add[i],] = tmp.prec[nodes.remove[i],]
    tmp.prec=prec.mat
    # swap prec mat cols
    prec.mat[,nodes.remove[i]] = tmp.prec[,nodes.add[i]]
    prec.mat[,nodes.add[i]] = tmp.prec[,nodes.remove[i]]
    # swap adj mat rows
    adj.mat[nodes.remove[i],] = tmp.adj[nodes.add[i],]
    adj.mat[nodes.add[i],] = tmp.adj[nodes.remove[i],]
    tmp.adj = adj.mat
    # swap adj mat cols
    adj.mat[,nodes.remove[i]] = tmp.adj[,nodes.add[i]]
    adj.mat[,nodes.add[i]] = tmp.adj[,nodes.remove[i]]
    # swap cov mat rows
    cov.mat[nodes.remove[i],] = tmp.cov.mat[nodes.add[i],]
    cov.mat[nodes.add[i],] = tmp.cov.mat[nodes.remove[i],]
    tmp.cov.mat = cov.mat
    # swap cov mat cols
    cov.mat[,nodes.remove[i]] = tmp.cov.mat[,nodes.add[i]]
    cov.mat[,nodes.add[i]] = tmp.cov.mat[,nodes.remove[i]]
    # swap data mat cols (only cols since it is nxp)
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

linear.weight.transformation.select = function(lambda.opt, cov.mat, prior.matrix,n, ebic.gamma,q.max){
  # Function for tuning the parameter q, where q determines the differentiation between groups: weights = 1-q*priorinfo
  # priorinfo may be estimated correlations, or 0/1 if no edge/edge. They are in the prior.matrix.
  # lambda.opt is the optimal lambda tuned by ordinary glasso. 
  # ebic.gamma is the eBIC parameter. It may be an array of several values. 
  # q.max is the largest q to consider.
  p=nrow(cov.mat)
  q_vals = seq(0,q.max,by=0.01) # Range of q-values
  likelihoods= rep(0,length(q_vals))
  theta.hats= array(dim=c(p,p,length(q_vals)))
  lambdas = rep(0,length(q_vals))
  # First: fit models, one for each q
  for (i in 1:length(q_vals)){
    weights = matrix(1-q_vals[i]*prior.matrix,nrow=p,byrow=T)
    #diag(weights) = rep(0,p)
    lambdas[i] = (lambda.opt*p^2)/(sum(weights)) # Same amount of total penalization
    #weighted = q_vals[i] != 0 # Should weights be used? Only false if q=0 (ordinary glasso)
    weighted=T
    fit.w = glasso(cov.mat, rho=lambdas[i]*weights,penalize.diagonal = !(weighted))
    theta.hats[,,i] = fit.w$wi # Save precision matrix
    likelihoods[i] = gaussianloglik(cov.mat,fit.w$wi,n)
  }
  aic.scores = apply(theta.hats,3,gaussianAIC,sample.cov=cov.mat,n=n)
  # Find eBIC for different choices of gamma, create list of objects
  ans = list()
  for (i in 1:length(ebic.gamma)){
    ebic=apply(theta.hats,3,eBIC,sample.cov=cov.mat, n=n,gamma=ebic.gamma[i])
    obj = list() # Object for current value of gamma
    obj$gamma = ebic.gamma[i] # Current value of gamma
    obj$ebic = ebic # ebic scores for all q
    obj$q = q_vals # All q values
    obj$sparsity=apply(theta.hats,3,sparsity) # sparsities for all q
    obj$opt.sparsity = sparsity(theta.hats[,,which.min(ebic)]) # Sparsity of optimal model
    obj$q.opt = q_vals[which.min(ebic)] # Optimal q
    obj$opt.ebic = min(ebic) # eBIC of optimal model
    obj$lambda.common = lambdas[which.min(ebic)] # Common lambda for optimal model
    obj$loglikelihoods= likelihoods # Likelihoods for all q
    obj$thetas = theta.hats # Precision matrices for all q
    obj$theta.opt = theta.hats[,,which.min(ebic)] # Optimal precision matrix
    obj$AIC = aic.scores# Include AIC as well
    ans[[i]] = obj
  }
  # now ans[[i]] is the optimal model wrt the eBIC when using gamma = ebic.gamma[i]
  return(ans) # Return list of values of optimal q and common penalty param, for each possible choice of ebic.gamma 
}