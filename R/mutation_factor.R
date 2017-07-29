# perform the discovery of K somatic mutational signatures given the observations x
"nmfLasso" <- function( x, K, background_signature = NULL, iterations = 20, lambda_rate = 0.01, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)
    
    if(verbose) {
        cat("Computing the initial values of beta by standard NMF...","\n")
    }
    
    # compute the initial values of beta
    beta = t(nmfDecomposition(x=t(x),r=K)$w)
    
    # add a signature to beta (leading to K+1 signatures in total) to explicitly model noise
    if(is.null(background_signature)) {
        background_signature = svd(x)$d
        if(min(background_signature)<0) {
            background_signature = background_signature + min(background_signature)
        }
    }
    beta = rbind(background_signature,beta)
    
    # normalize beta so that each row sums to 1
    beta = beta / rowSums(beta)
    
    if(verbose) {
        cat("Performing the discovery of the signatures by NMF with Lasso...","\n")
    }
    
    # perform the discovery
    results = nmfLassoDecomposition(x,beta,iterations,lambda_rate)
    
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
        cat("Performing a maximum of",iterations,"EM iterations...","\n")
    }
    
    # repeat a 2 step algorithm iteratively, where first alpha is estimated by Non-Negative Linear Least Squares 
    # and, then, beta is estimated by Non-Negative Lasso
    best_loglik = -Inf
    best_alpha = NA
    best_beta = NA
    for(i in 1:iterations) {
        
        # initialize the value of the log-likelihood for the current iteration
        loglik[i] = 0
            
        # normalize the rate of the current signature to sum to 1
        beta = beta / rowSums(beta)
        
        # update alpha independently for each patient by Non-Negative Linear Least Squares
        for(j in 1:n) {
            alpha[j,] = nnls(t(beta),as.vector(x[j,]))$x
        }
        
        # update beta by Non-Negative Lasso
        for(k in 1:J) {
            
            # compute independently for each trinucleotide the error between the observed counts, i.e., x, 
            # and the ones predicted by the first signature (i.e., which represents the noise model)
            error = x[,k] - alpha[,1] * beta[1,k]
            
            # estimate beta for the remaining signatues by Non-Negative Lasso to ensure sparsity
            max_lambda_value = max(abs(t(alpha[,2:K]) %*% error))
            lambda_value = max(max_lambda_value*lambda_rate,1e-4)
            beta[2:K,k] = as.vector(nnlasso(x=alpha[,2:K],y=error,lambda=lambda_value,intercept=FALSE,normalize=FALSE,path=FALSE)$coef[2,])
            
            # update the log-likelihood for the current iteration
            curr_loglik = -sum((x[,k] - alpha %*% beta[,k])^2) - lambda_value * sum(beta[2:K,k])
            loglik[i] = loglik[i] + curr_loglik
            
        }
        
        # count how many EM steps do not lead to a better log-likelihood
        if(loglik[i]>best_loglik) {
            best_loglik = loglik[i]
            best_alpha = alpha
            best_beta = beta
        }
        
        if(verbose) {
            cat("Progress",paste0((i/iterations)*100,"%..."),"\n")
        }
        
    }
    alpha = best_alpha
    beta = best_beta
    
    # normalize the rate of the current signature to sum to 1
    beta = beta / rowSums(beta)
    
    return(list(alpha=alpha,beta=beta,best_loglik=best_loglik,loglik=loglik))
    
}
