# perform the discovery of K (unknown) somatic mutational signatures given a set of observations x
"nmfLasso" <- function( x, K = 2:15, background_signature = NULL, lambda_values = c(0.01, 0.05, 0.10, 0.15), cross_validation_entries = 0.15, iterations = 20, max_iterations_lasso = 10000, num_processes = Inf, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)
    
    # set a percentage of cross_validation_entries entries to 0 in order to perform cross validation
    x_cv = x
    valid_entries = which(x_cv>0,arr.ind=TRUE)
    cv_entries = valid_entries[sample(1:nrow(valid_entries),size=round(nrow(valid_entries)*cross_validation_entries),replace=FALSE),]
    x_cv[cv_entries] = 0
    
    # perform a grid search to estimate the best values of K and lambda
    cont = 0
    if(verbose) {
        cat("Performing a grid search to estimate the best values of K and lambda...","\n")
    }
    
    # setting up parallel execution
    if(is.na(num_processes) || is.null(num_processes)) {
        parallel = NULL
    }
    else if(num_processes==Inf) {
        cores = as.integer((detectCores()-1))
        if(cores < 2) {
            parallel = NULL
        }
        else {
            num_processes = cores
            parallel = makeCluster(num_processes,outfile="")
            clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
        }
    }
    else {
        parallel = makeCluster(num_processes,outfile="")
        clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
    }
    
    if(verbose && !is.null(parallel)) {
        cat("Executing",num_processes,"processes via parallel...","\n")
    }
    
    # structure to save the results of the grid search
    grid_search = array(list(),c(length(K),length(lambda_values)))
    rownames(grid_search) = paste0(as.character(K),"_signatures")
    colnames(grid_search) = paste0(as.character(lambda_values),"_lambda")
    
    # structure to save the starting values of beta for each K
    starting_beta = array(list(),c(length(K),1))
    rownames(starting_beta) = paste0(as.character(K),"_signatures")
    colnames(starting_beta) = "Value"
    
    # consider all the values for K
    pos_k = 0
    for(k in K) {
            
        # get the first k signatures to be used for the current configuration
        curr_beta = basis(nmf(t(x),rank=k))
        curr_beta = t(curr_beta)
        pos_k = pos_k + 1
        starting_beta[[pos_k,1]] = curr_beta
        
        # consider all the values for lambda
        pos_l = 0
        for(l in lambda_values) {
            
            # perform the inference
            curr_results = nmfLassoK(x = x_cv, 
                                     K = k, 
                                     beta = curr_beta, 
                                     background_signature = background_signature, 
                                     lambda_rate = l, 
                                     iterations = iterations, 
                                     max_iterations_lasso = max_iterations_lasso, 
                                     num_processes = NULL, 
                                     parallel = parallel, 
                                     seed = round(runif(1)*100000), 
                                     verbose = FALSE)
            
            # save the results for the current configuration
            pos_l = pos_l + 1
            grid_search[[pos_k,pos_l]] = curr_results
            
            if(verbose) {
                cont = cont + 1
                cat("Progress",paste0(round((cont/(length(K)*length(lambda_values)))*100,digits=3),"%..."),"\n")
            }
            
        }
        
    }
    
    # close parallel
    stopCluster(parallel)
    
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
    results = list(grid_search=grid_search,mean_squared_error=mean_squared_error,starting_beta=starting_beta)
    
    return(results)
    
}

# estimate the range of lambda values to be considered in the signature inference
"evaluateLambdaRange" <- function( x, K = 8, background_signature = NULL, lambda_values = c(0.01, 0.05, 0.10, 0.15), iterations = 20, max_iterations_lasso = 10000, num_processes = Inf, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)
    
    # compute the initial values of beta
    if(verbose) {
        cat("Computing the initial values of beta by standard NMF...","\n")
    }
    beta = basis(nmf(t(x),rank=K))
    beta = t(beta)
    
    if(verbose) {
        cat("Estimating the signatures for the different values of lambda...","\n")
    }
    
    # setting up parallel execution
    if(is.na(num_processes) || is.null(num_processes)) {
        parallel = NULL
    }
    else if(num_processes==Inf) {
        cores = as.integer((detectCores()-1))
        if(cores < 2) {
            parallel = NULL
        }
        else {
            num_processes = cores
            parallel = makeCluster(num_processes,outfile="")
            clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
        }
    }
    else {
        parallel = makeCluster(num_processes,outfile="")
        clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
    }
    
    if(verbose && !is.null(parallel)) {
        cat("Executing",num_processes,"processes via parallel...","\n")
    }
    
    # structure to save the estimated signatures
    lambda_results = array(list(),c(length(K),length(lambda_values)))
    rownames(lambda_results) = paste0(as.character(K),"_signatures")
    colnames(lambda_results) = paste0(as.character(lambda_values),"_lambda")
    
    # perform signature discovery for all the values of lambda
    cont = 0
    for(l in lambda_values) {
            
        # perform the inference
        curr_results = nmfLassoK(x = x, 
                                 K = K, 
                                 beta = beta, 
                                 background_signature = background_signature, 
                                 lambda_rate = l, 
                                 iterations = iterations, 
                                 max_iterations_lasso = max_iterations_lasso, 
                                 num_processes = NULL, 
                                 parallel = parallel, 
                                 seed = round(runif(1)*100000), 
                                 verbose = FALSE)
                                 
        # save the results for the current configuration
        cont = cont + 1
        lambda_results[[1,cont]] = curr_results
            
        if(verbose) {
            cat("Progress",paste0(round((cont/(length(K)*length(lambda_values)))*100,digits=3),"%..."),"\n")
        }
        
    }
    
    # close parallel
    stopCluster(parallel)
    
    # save the results
    results = lambda_results
    
    return(results)
    
}

# perform the discovery of K somatic mutational signatures given a set of observations x
"nmfLassoK" <- function( x, K, beta = NULL, background_signature = NULL, lambda_rate = 0.05, iterations = 20, max_iterations_lasso = 10000, num_processes = Inf, parallel = NULL, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)
    
    # compute the initial values of beta
    if(is.null(beta)) {
        if(verbose) {
            cat("Computing the initial values of beta by standard NMF...","\n")
        }
        beta = basis(nmf(t(x),rank=K))
        beta = t(beta)
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
    
    # setting up parallel execution
    close_parallel = FALSE
    if(is.null(parallel)) {
        if(is.na(num_processes) || is.null(num_processes)) {
            parallel = NULL
        }
        else if(num_processes==Inf) {
            cores = as.integer((detectCores()-1))
            if(cores < 2) {
                parallel = NULL
            }
            else {
                num_processes = cores
                parallel = makeCluster(num_processes,outfile="")
                clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
                close_parallel = TRUE
            }
        }
        else {
            parallel = makeCluster(num_processes,outfile="")
            clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
            close_parallel = TRUE
        }
        
        if(verbose && !is.null(parallel)) {
            cat("Executing",num_processes,"processes via parallel...","\n")
        }
    }
    
    # perform the discovery of the signatures
    result = tryCatch({
        nmfLassoDecomposition(x,beta,lambda_rate,iterations,max_iterations_lasso,parallel,verbose)
    }, error = function(e) {
        warning("Lasso did not converge, you should try a lower value of lambda! Current settings: K = ",K,", lambda_rate = ",lambda_rate,".")
        return(list(alpha=NA,beta=NA,best_loglik=NA,loglik_progression=rep(NA,iterations)))
    })
    
    results = 
    
    # close parallel
    if(close_parallel) {
        stopCluster(parallel)
    }
    
    return(results)
    
}

# perform de novo discovery of somatic mutational signatures using NMF with Lasso to ensure sparsity
"nmfLassoDecomposition" <- function( x, beta, lambda_rate = 0.05, iterations = 20, max_iterations_lasso = 10000, parallel = NULL, verbose = TRUE ) {
    
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
    
    # repeat a 2 step algorithm iteratively, where first alpha is estimated by Non-Negative Linear Least Squares 
    # and, then, beta is estimated by Non-Negative Lasso
    best_loglik = -Inf
    best_alpha = NA
    best_beta = NA
    for(i in 1:iterations) {
        
        # initialize the value of the log-likelihood for the current iteration
        loglik[i] = 0
        
        # update alpha independently for each patient by Non-Negative Linear Least Squares
        if(is.null(parallel)) {
            for(j in 1:n) {
                alpha[j,] = nnls(t(beta),as.vector(x[j,]))$x
            }
        }
        else {
            
            # compute alpha in parallel
            j = 1:n
            res_clusterEvalQ = clusterEvalQ(parallel,library("nnls"))
            clusterExport(parallel,varlist=c("x","beta"),envir=environment())
            alpha_res = parLapply(parallel,j,function(j) {
                return(nnls(t(beta),as.vector(x[j,]))$x)
            })
            
            # reduce the results
            for(j in 1:n) {
                alpha[j,] = alpha_res[[j]]
            }
            
        }
        
        # update beta by Non-Negative Lasso
        if(is.null(parallel)) {
            
            for(k in 1:J) {
                
                # compute independently for each trinucleotide the error between the observed counts, i.e., x, 
                # and the ones predicted by the first signature (i.e., which represents the noise model)
                error = x[,k] - alpha[,1] * beta[1,k]
                
                # estimate beta for the remaining signatues by Non-Negative Lasso to ensure sparsity
                if(is.na(lambda_values[k])) {
                    max_lambda_value = max(abs(t(alpha[,2:K]) %*% error))
                    lambda_values[k] = max_lambda_value * lambda_rate
                }
                beta[2:K,k] = as.vector(nnlasso(x = alpha[,2:K], 
                                                y = error, 
                                                family = "normal", 
                                                lambda = lambda_values[k], 
                                                intercept = FALSE, 
                                                normalize = FALSE, 
                                                maxiter = max_iterations_lasso, 
                                                path = FALSE)$coef[2,])
                
                # update the log-likelihood for the current iteration
                curr_loglik = -sum((x[,k] - alpha %*% beta[,k])^2) - lambda_values[k] * sum(beta[2:K,k])
                loglik[i] = loglik[i] + curr_loglik
                
            }
            
        }
        else {
            
            # if this is the first iteration, compute the values of lambda sequentially
            if(is.na(lambda_values[1])) {
                
                for(k in 1:J) {
                    error = x[,k] - alpha[,1] * beta[1,k]
                    max_lambda_value = max(abs(t(alpha[,2:K]) %*% error))
                    lambda_values[k] = max_lambda_value * lambda_rate
                }
                
            }
            
            # perform the computations of beta in parallel
            k = 1:J
            res_clusterEvalQ = clusterEvalQ(parallel,library("nnlasso"))
            clusterExport(parallel,varlist=c("x","alpha","beta","lambda_values","K","max_iterations_lasso"),envir=environment())
            beta_res = parLapply(parallel,1:J,function(k) {
                
                # compute independently for each trinucleotide the error between the observed counts, i.e., x, 
                # and the ones predicted by the first signature (i.e., which represents the noise model)
                error = x[,k] - alpha[,1] * beta[1,k]
                
                # estimate beta for the remaining signatues by Non-Negative Lasso to ensure sparsity
                beta[2:K,k] = as.vector(nnlasso(x = alpha[,2:K], 
                                                y = error, 
                                                family = "normal", 
                                                lambda = lambda_values[k], 
                                                intercept = FALSE, 
                                                normalize = FALSE, 
                                                maxiter = max_iterations_lasso, 
                                                path = FALSE)$coef[2,])
                
                # compute the log-likelihood for the current iteration
                curr_loglik = -sum((x[,k] - alpha %*% beta[,k])^2) - lambda_values[k] * sum(beta[2:K,k])
                
                return(list(beta=beta[2:K,k],loglik=curr_loglik))
                
            })
            
            # reduce the results
            for(k in 1:J) {
                
                beta[2:K,k] = beta_res[[k]][["beta"]]
                loglik[i] = loglik[i] + beta_res[[k]][["loglik"]]
                
            }
            
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
    
    # check if the likelihood is increasing
    if(length(loglik)>1) {
        cont = 1
        for(i in 2:length(loglik)) {
            if(loglik[i]>loglik[(i-1)]) {
                cont = cont + 1
            }
        }
        if(cont/iterations<0.5) {
            warning("The likelihood is not increasing during the EM, you should try a lower value of lambda! Current settings: K = ",K,", lambda_rate = ",lambda_rate,".")
        }
    }
    
    # normalize the rate of the signatures to sum to 1
    beta = beta / rowSums(beta)
        
    # final computation of alpha by Non-Negative Linear Least Squares using the normalized version o the best beta
    if(is.null(parallel)) {
        for(j in 1:n) {
            alpha[j,] = nnls(t(beta),as.vector(x[j,]))$x
        }
    }
    else {
        
        # compute alpha in parallel
        j = 1:n
        res_clusterEvalQ = clusterEvalQ(parallel,library("nnls"))
        clusterExport(parallel,varlist=c("x","beta"),envir=environment())
        alpha_res = parLapply(parallel,j,function(j) {
            return(nnls(t(beta),as.vector(x[j,]))$x)
        })
        
        # reduce the results
        for(j in 1:n) {
            alpha[j,] = alpha_res[[j]]
        }
        
    }
    
    # save the results
    results = list(alpha=alpha,beta=beta,best_loglik=best_loglik,loglik_progression=loglik)
    
    return(results)
    
}
