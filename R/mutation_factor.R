# perform the discovery of K (unknown) somatic mutational signatures given a set of observations x
"nmfLasso" <- function( x, K = 2:15, background_signature = NULL, lambda_values = seq(0.3, 0.7, by = 0.1), iterations = 20, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)
    
    # set 10% of the entries to 0 in order to perform cross validation
    x_cv = x
    valid_entries = which(x_cv>0,arr.ind=TRUE)
    cv_entries = valid_entries[sample(1:nrow(valid_entries),size=round(nrow(valid_entries)*0.10),replace=FALSE),]
    x_cv[cv_entries] = 0
    
    # perform a grid search to estimate the best values of K and lambda
    cont = 0
    if(verbose) {
        cat("Performing a grid search to estimate the best values of K and lambda...","\n")
    }
    
    # structure to save the results of the grid search
    grid_search = array(list(),c(length(K),length(lambda_values)))
    rownames(grid_search) = paste0(as.character(K),"_signatures")
    colnames(grid_search) = paste0(as.character(lambda_values),"_lambda")
    
    # consider all the values for K
    pos_k = 0
    for(k in K) {
            
        # get the first k signatures to be used for the current configuration
        curr_beta = t(nmfDecomposition(x=t(x),r=k)$w)
        pos_k = pos_k + 1
        
        # consider all the values for lambda
        pos_l = 0
        for(l in lambda_values) {
            
            # perform the inference
            curr_results = nmfLassoK(x=x_cv,K=k,beta=curr_beta,background_signature=background_signature,lambda_rate=l,iterations=iterations,seed=round(runif(1)*100000),verbose=FALSE)
            
            # save the results for the current configuration
            pos_l = pos_l + 1
            grid_search[[pos_k,pos_l]] = curr_results
            
            if(verbose) {
                cont = cont + 1
                cat("Progress",paste0(round((cont/(length(K)*length(lambda_values)))*100,digits=3),"%..."),"\n")
            }
            
        }
        
    }
    
    if(verbose) {
        cat("Evaluating the results of the cross validation in terms of mean squared error...","\n")
    }
    
    # structure to save the mean squared errors
    mean_squared_error = array(NA,c(length(K),length(lambda_values)))
    rownames(mean_squared_error) = paste0(as.character(K),"_signatures")
    colnames(mean_squared_error) = paste0(as.character(lambda_values),"_lambda")
    
    # assess the results of the cross validation by mean squared error
    J = dim(x)[2]
    pos_k = 0
    for(k in K) {
        pos_k = pos_k + 1
        pos_l = 0
        for(l in lambda_values) {
            
            # consider the current configuration
            pos_l = pos_l + 1
            curr_alpha = grid_search[[pos_k,pos_l]][["alpha"]]
            curr_beta = grid_search[[pos_k,pos_l]][["beta"]]
            
            # compute the mean squared error
            error = 0
            for(i in 1:J) {
                # compute for each trinucleotide the error between the observed counts, i.e., x, and the predicted ones
                curr_error = mean((x[,i] - as.vector(curr_alpha %*% curr_beta[,i]))^2)
                error = error + curr_error
            }
            error = error / J
            mean_squared_error[pos_k,pos_l] = error
            
        }
    }
    
    # save the results
    results = list(grid_search=grid_search,mean_squared_error=mean_squared_error)
    
    return(results)
    
}

# perform the discovery of K somatic mutational signatures given a set of observations x
"nmfLassoK" <- function( x, K, beta = NULL, background_signature = NULL, lambda_rate = 0.5, iterations = 20, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)
    
    # compute the initial values of beta
    if(is.null(beta)) {
        if(verbose) {
            cat("Computing the initial values of beta by standard NMF...","\n")
        }
        beta = t(nmfDecomposition(x=t(x),r=K)$w)
    }
    
    # add a signature to beta (leading to K+1 signatures in total) to explicitly model noise
    if(is.null(background_signature)) {
        background_signature = svd(x)$d
        if(min(background_signature)<0) {
            background_signature = background_signature + abs(min(background_signature))
        }
    }
    beta = rbind(background_signature,beta)
    
    if(verbose) {
        cat("Performing the discovery of the signatures by NMF with Lasso...","\n")
    }
    
    # perform the discovery of the signatures
    results = nmfLassoDecomposition(x,beta,lambda_rate,iterations,verbose)
    
    return(results)
    
}

# perform de novo discovery of somatic mutational signatures using NMF with Lasso to ensure sparsity
"nmfLassoDecomposition" <- function( x, beta, lambda_rate = 0.5, iterations = 20, verbose = TRUE ) {
    
    # n is the number of observations in x, i.e., the patients
    n = dim(x)[1]
    
    # J is the number of trinucleotides, i.e., 96 categories
    J = dim(x)[2]
    
    # K is the number of signatures to be called
    K = dim(beta)[1]
    
    # initialize alpha
    alpha = matrix(0,nrow=n,ncol=K)
    rownames(alpha) = 1:nrow(alpha)
    colnames(alpha) = rownames(beta)
    
    # structure where to save the log-likelihood at each iteration 
    loglik = rep(NA,iterations)
    
    # structure where to save the lambda values for each J
    lambda_values = rep(NA,J)
    
    if(verbose) {
        cat("Performing a total of",iterations,"EM iterations...","\n")
    }
            
    # normalize the rate of the current signature to sum to 1
    beta = (beta / rowSums(beta))
    
    # repeat a 2 step algorithm iteratively, where first alpha is estimated by Non-Negative Linear Least Squares 
    # and, then, beta is estimated by Non-Negative Lasso
    best_loglik = -Inf
    best_alpha = NA
    best_beta = NA
    for(i in 1:iterations) {
        
        # initialize the value of the log-likelihood for the current iteration
        loglik[i] = 0
        
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
            if(is.na(lambda_values[k])) {
                max_lambda_value = max(abs(t(alpha[,2:K]) %*% error)) / n
                lambda_values[k] = max(max_lambda_value*lambda_rate,1e-4)
            }
            beta[2:K,k] = as.vector(nnlasso(x=alpha[,2:K],y=error,family="normal",lambda=lambda_values[k],intercept=FALSE,normalize=FALSE,path=FALSE)$coef[2,])
            
            # update the log-likelihood for the current iteration
            curr_loglik = -sum((x[,k] - alpha %*% beta[,k])^2) - lambda_values[k] * sum(beta[2:K,k])
            loglik[i] = loglik[i] + curr_loglik
            
        }
        
        # save the grid_search at maximum log-likelihood
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
    
    # normalize the rate of the signatures to sum to 1
    beta = beta / rowSums(beta)
        
    # final computation of alpha by Non-Negative Linear Least Squares using the normalized version o the best beta
    for(j in 1:n) {
        alpha[j,] = nnls(t(beta),as.vector(x[j,]))$x
    }
    
    # save the results
    results = list(alpha=alpha,beta=beta,best_loglik=best_loglik,loglik_progression=loglik)
    
    return(results)
    
}
