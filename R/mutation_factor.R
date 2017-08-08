# perform the discovery by cross validation of K (unknown) somatic mutational signatures given a set of observations x
"nmfLasso" <- function( x, K = 2:15, starting_beta = NULL, background_signature = NULL, lambda_values = c(0.01, 0.05, 0.10, 0.20, 0.30), cross_validation_entries = 0.15, cross_validation_iterations = 10, iterations = 20, max_iterations_lasso = 10000, num_processes = Inf, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)
    
    # perform a grid search to estimate the best values of K and lambda
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
    
    # structure to save the results of the grid search for all the cross_validation_iterations
    grid_search_iterations = list()
    
    # structure to save the starting values of beta for each K for all the cross_validation_iterations
    starting_beta_iterations = list()

    # repeat the estimation for a number of cross_validation_iterations
    for(cv_iteration in 1:cross_validation_iterations) {
                
        if(verbose) {
            cat(paste0("Performing cross validation iteration ",cv_iteration," out of ",cross_validation_iterations,"..."),"\n")
        }

        # set a percentage of cross_validation_entries entries to 0 in order to perform cross validation
        x_cv = x
        valid_entries = which(x_cv>0,arr.ind=TRUE)
        cv_entries = valid_entries[sample(1:nrow(valid_entries),size=round(nrow(valid_entries)*cross_validation_entries),replace=FALSE),]
        x_cv[cv_entries] = 0
        
        # structure to save the results of the grid search
        grid_search = array(list(),c(length(K),length(lambda_values)))
        rownames(grid_search) = paste0(as.character(K),"_signatures")
        colnames(grid_search) = paste0(as.character(lambda_values),"_lambda")
        
        # structure to save the starting values of beta for each K
        if(is.null(starting_beta)) {
            starting_beta = array(list(),c(length(K),1))
            rownames(starting_beta) = paste0(as.character(K),"_signatures")
            colnames(starting_beta) = "Value"
        }

        # consider all the values for K
        cont = 0
        pos_k = 0
        for(k in K) {
                
            # get the first k signatures to be used for the current configuration
            pos_k = pos_k + 1
            if(is.null(starting_beta[[pos_k,1]])) {
                curr_beta = basis(nmf(t(x),rank=k))
                curr_beta = t(curr_beta)
                starting_beta[[pos_k,1]] = curr_beta
            }
            else {
                curr_beta = starting_beta[[pos_k,1]]
            }
            
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
        
        # save the results for the current iteration
        starting_beta_iterations[[cv_iteration]] = starting_beta
        grid_search_iterations[[cv_iteration]] = grid_search

    }
    
    if(verbose) {
        cat("Evaluating the results of the cross validation in terms of mean squared error...","\n")
    }

    # structure to save the mean squared errors for all the cross_validation_iterations
    mean_squared_error_iterations = list()

    # repeat the estimation for a number of cross_validation_iterations
    for(cv_iteration in 1:cross_validation_iterations) {
    
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
                curr_alpha = grid_search_iterations[[cv_iteration]][[pos_k,pos_l]][["alpha"]]
                curr_beta = grid_search_iterations[[cv_iteration]][[pos_k,pos_l]][["beta"]]
                
                # compute the mean squared error
                if(!is.na(curr_alpha)&&!is.na(curr_beta)) {
                    error = 0
                    for(i in 1:J) {
                        # compute for each trinucleotide the error between the observed counts, i.e., x, and the predicted ones
                        curr_error = mean((x[,i] - as.vector(curr_alpha %*% curr_beta[,i]))^2)
                        error = error + curr_error
                    }
                    error = error / J
                }
                else {
                    error = NA
                }

                mean_squared_error[pos_k,pos_l] = error
                
            }
        }

        # save the results for the current iteration
        mean_squared_error_iterations[[cv_iteration]] = mean_squared_error

    }
    
    if(verbose) {
        cat("Estimating the best configuration...","\n")
    }
    
    # structure to save the mean squared errors averaged over the cross_validation_iterations
    mean_squared_error = array(list(),c(length(K),length(lambda_values)))
    rownames(mean_squared_error) = paste0(as.character(K),"_signatures")
    colnames(mean_squared_error) = paste0(as.character(lambda_values),"_lambda")
    
    # average the results over the different cross_validation_iterations
    for(i in 1:length(mean_squared_error_iterations)) {
        curr_mean_squared_error = mean_squared_error_iterations[[i]]
        for(j in 1:nrow(mean_squared_error)) {
            for(k in 1:ncol(mean_squared_error)) {
                mean_squared_error[j,k] = list(c(unlist(mean_squared_error[j,k]),curr_mean_squared_error[j,k]))
            }
        }
    }
    for(j in 1:nrow(mean_squared_error)) {
        for(k in 1:ncol(mean_squared_error)) {
            mean_squared_error[j,k] = mean(unlist(mean_squared_error[j,k]),na.rm=TRUE)
        }
    }
    
    # find the best configuration
    best_j = NA
    best_k = NA
    best_result = NA
    for(j in 1:nrow(mean_squared_error)) {
        for(k in 1:ncol(mean_squared_error)) {
            if(is.na(best_result)&&!is.nan(as.numeric(mean_squared_error[j,k]))) {
                best_result = as.numeric(mean_squared_error[j,k])
                best_j = j
                best_k = k
            }
            else if(!is.nan(as.numeric(mean_squared_error[j,k]))) {
                if(as.numeric(mean_squared_error[j,k])<best_result) {
                    best_result = as.numeric(mean_squared_error[j,k])
                    best_j = j
                    best_k = k
                }
            }
        }
    }
    
    # compute the signatures for the best configuration
    if(!is.na(best_j)&&!is.na(best_k)) {

        # set the starting beta values
        if(is.null(starting_beta)) {
            curr_beta = NULL
        }
        else {
            curr_beta = starting_beta[[K[best_j],1]]
        }

        # perform the inference
        curr_results = nmfLassoK(x = x, 
                                 K = K[best_j], 
                                 beta = curr_beta, 
                                 background_signature = background_signature, 
                                 lambda_rate = lambda_values[best_k], 
                                 iterations = iterations, 
                                 max_iterations_lasso = max_iterations_lasso, 
                                 num_processes = NULL, 
                                 parallel = parallel, 
                                 seed = round(runif(1)*100000), 
                                 verbose = TRUE)

        # save the results
        best_configuration = curr_results
        best_configuration[["starting_beta"]] = curr_beta
        best_configuration[["K"]] = K[best_j]
        best_configuration[["lambda_rate"]] = lambda_values[best_k]

    }
    else {
        best_configuration = NA
    }
    
    # close parallel
    stopCluster(parallel)
    
    # save the results
    results = list(grid_search=grid_search_iterations,mean_squared_error=mean_squared_error_iterations,starting_beta=starting_beta_iterations,best_configuration=best_configuration)
    
    return(results)
    
}

# estimate the range of lambda values to be considered in the signature inference
"evaluateLambdaRange" <- function( x, K = 8, beta = NULL, background_signature = NULL, lambda_values = c(0.01, 0.05, 0.10, 0.20, 0.30), iterations = 20, max_iterations_lasso = 10000, num_processes = Inf, seed = NULL, verbose = TRUE ) {
    
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
"nmfLassoK" <- function( x, K, beta = NULL, background_signature = NULL, lambda_rate = 0.10, iterations = 20, max_iterations_lasso = 10000, num_processes = Inf, parallel = NULL, seed = NULL, verbose = TRUE ) {
    
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
    results = tryCatch({
        nmfLassoDecomposition(x,beta,lambda_rate,iterations,max_iterations_lasso,parallel,verbose)
    }, error = function(e) {
        warning("Lasso did not converge, you should try a lower value of lambda! Current settings: K = ",K,", lambda_rate = ",lambda_rate,".")
        return(list(alpha=NA,beta=NA,best_loglik=NA,loglik_progression=rep(NA,iterations)))
    })
    
    # close parallel
    if(close_parallel) {
        stopCluster(parallel)
    }
    
    return(results)
    
}

# perform de novo discovery of somatic mutational signatures using NMF with Lasso to ensure sparsity
"nmfLassoDecomposition" <- function( x, beta, lambda_rate = 0.10, iterations = 20, max_iterations_lasso = 10000, parallel = NULL, verbose = TRUE ) {
    
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
