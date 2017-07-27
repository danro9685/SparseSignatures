# perform the discovery of K somatic mutational signatures given the observations x
"nmfLasso" <- function( x, K, genome_freq = NULL, iterations = 20, lambda_rate = 0.01, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)
    
    if(verbose) {
        cat("Computing the initial beta values by standard NMF...","\n")
    }
    
    # compute the initial values for beta
    beta = t(nmfDecomposition(x=t(x),r=K)$w)
    
    # normalize beta so that each row sums to 1
    beta = beta / rowSums(beta)
    
    # add a signature to beta (leading to K+1 signatures in total) to explicitly model noise
    if(is.null(genome_freq)) {
        genome_freq = svd(x)$d
        genome_freq = genome_freq/sum(genome_freq)
    }
    beta = rbind(genome_freq,beta)
    
    if(verbose) {
        cat("Performing the discovery of the signatures by NMF with Lasso...","\n")
    }
    
    # perform the discovery
    results = nmfLassoDecomposition(x,beta,iterations,lambda_rate)
    
    # normalize the signatures to sum to 1
    beta_rate = results$beta/rowSums(results$beta)
    results[["beta_rate"]] = beta_rate
    
    
    return(results)
    
}

# perform de novo discovery of somatic mutational signatures using NMF with Lasso to ensure sparsity
"nmfLassoDecomposition" <- function( x, beta, iterations = 20, lambda_rate = 0.01, verbose = TRUE ) {
    
    # n is the number of observations in x, i.e., the patients
    n = dim(x)[1]
    
    # J is the number of trinucleotides, i.e., 96 categories
    J = dim(x)[2]
    
    # K is the number of signatures to be called
    K = dim(beta)[1]
    
    # initialize alpha
    alpha = matrix(0,nrow=n,ncol=K)
    
    # structure where to save the log-likelihood at each iteration 
    loglik = rep(NA,iterations)
    
    if(verbose) {
        cat("Performing",iterations,"iterations...","\n")
    }
    
    # repeat a 2 step algorithm iteratively, where first alpha is estimated by Non-Negative Linear Least Squares 
    # and, then, beta is estimated by Non-Negative Lasso
    for(i in 1:iterations) {
        
        # update alpha independently for each patient by Non-Negative Linear Least Squares
        for(j in 1:n) {
            alpha[j,] = nnls(t(beta),as.vector(x[j,]))$x
        }
        
        # compute independently for each trinucleotide the error between the observed counts, i.e., x,
        # and the ones predicted by the first signature (i.e., which represents the noise model)
        error = rep(list(),J)
        max_lambda_value = rep(NA,J)
        for(k in 1:J) {
            error[[k]] = x[,k] - alpha[,1] * beta[1,k]
            max_lambda_value[k] = max(abs(t(alpha[,2:K]) %*% error[[k]]))
        }
        mean_max_lambda_value = mean(max_lambda_value)
        lambda_value = max(mean_max_lambda_value*lambda_rate,1e-4)
        
        # estimate beta for the remaining signatues by Non-Negative Lasso to ensure sparsity
        for(k in 1:J) {
            beta[2:K,k] = as.vector(nnlasso(x=alpha[,2:K],y=error[[k]],lambda=lambda_value,intercept=FALSE,normalize=FALSE,path=FALSE)$coef[2,])
        }
        
        # compute the log-likelihood for the current iteration
        loglik[i] = -sum((x - alpha %*% beta)^2) - lambda_value * sum(beta[2:K,])
    
        if(verbose) {
            cat("Progress",paste0((i/iterations)*100,"%..."),"\n")
        }
        
    }
    
    return(list(alpha=alpha,beta=beta,loglik=loglik))
    
}
