#' Perform the discovery by cross validation of K (unknown) somatic mutational signatures given a set of observations x. The estimation can slow down because of 
#' memory usage, when I high number of cross validation repetitions is asked and when the grid search is performed for a lot of configurations. In this case, 
#' we advice to split the computation into multiple smaller sets. 
#' @title nmf.LassoCV
#' @param x count matrix.
#' @param K a range of numeric value (each of them greater than 1) indicating the number of signatures to be discovered.
#' @param starting_beta a list of starting beta value for each configuration of K. If it is NULL, starting betas are estimated by 
#' NMF.
#' @param background_signature background signature to be used. If not provided, a warning is thrown.
#' @param nmf_runs number of iteration of NMF to be performed for a robust estimation of starting beta. If beta is not NULL, 
#' this parameter is ignored.
#' @param lambda_values range of values of LASSO to be used between 0 and 1. This value should be greater than 0. 1 is the value of 
#' LASSO 
#' that would shrink all the signatures to 0 within one step. The higher lambda_rate is, the sparser are the resulting signatures, 
#' but too large values result in a poor fit of the counts.
#' @param cross_validation_entries Percentage of cells in the count matrix to be replaced by 0s.
#' @param cross_validation_iterations For each configuration, the first time the signatures are discovered form a matrix with a 
#' ercentage of values replaced by 0s. This may result in a poor. This parameter is the number of restarts to be performed to 
#' improve this estimate.
#' @param cross_validation_repetitions Number of time cross-validation should be repeated. Higher values result in better estimate, but are computationally expensive.
#' @param iterations Number of iterations to be performed. Each iteration correspond to a first step where the counts are fitted 
#' and a second step where sparsity is enhanced.
#' @param max_iterations_lasso Number of maximum iterations to be performed during the sparsification.
#' @param num_processes Number of processes to be used during parallel execution. If executing in single process mode, 
#' this is ignored.
#' @param seed Seed for reproducibility.
#' @param verbose boolean; Shall I print all messages?
#' @return A list corresponding with 3 elements: grid_search, starting_beta and mean_squared_error. Here, grid_search provides all the results of the executions within 
#' the grid search; starting_beta is the set of initial values of beta used for each configuration and mean_squared_error is the mean squared error between the 
#' observed counts and the predicted ones for each configuration.
#' @export nmf.LassoCV
#' @import NMF
#' @import nnlasso
#' @import nnls
#' @import parallel
#'
"nmf.LassoCV" <- function( x, K = 3:10, starting_beta = NULL, background_signature = NULL, nmf_runs = 10, lambda_values = c(0.10, 0.20, 0.30), cross_validation_entries = 0.05, cross_validation_iterations = 5, cross_validation_repetitions = 10, iterations = 20, max_iterations_lasso = 10000, num_processes = Inf, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)

    nmf_method <- "nmf_standard"
    
    # perform a grid search to estimate the best values of K and lambda
    if(verbose) {
        cat("Performing a grid search to estimate the best values of K and lambda with a total of",cross_validation_repetitions,"cross validation repetitions...","\n")
    }
    
    # setting up parallel execution
    if(is.na(num_processes) || is.null(num_processes)) {
        parallel <- NULL
    }
    else if(num_processes==Inf) {
        cores <- as.integer((detectCores()-1))
        if(cores < 2) {
            parallel <- NULL
        }
        else {
            num_processes <- cores
            parallel <- makeCluster(num_processes,outfile="")
            clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
        }
    }
    else {
        parallel <- makeCluster(num_processes,outfile="")
        clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
    }
    
    if(verbose && !is.null(parallel)) {
        cat("Executing",num_processes,"processes via parallel...","\n")
    }

    if(verbose) {
        cat("Starting cross validation with a total of",cross_validation_repetitions,"repetitions...","\n")
    }

    # now starting cross validations
    results <- list()
    for(cv_repetitions in 1:cross_validation_repetitions) {

        if(verbose) {
            cat(paste0("Performing repetition ",cv_repetitions," out of ",cross_validation_repetitions,"..."),"\n")
        }
    
        # structure to save the results of the grid search for all the cross_validation_iterations
        grid_search_iterations <- list()

        # repeat the estimation for a number of cross_validation_iterations
        for(cv_iteration in 1:cross_validation_iterations) {
                    
            if(verbose) {
                cat(paste0("Performing cross validation iteration ",cv_iteration," out of ",cross_validation_iterations,"..."),"\n")
            }

            # set a percentage of cross_validation_entries entries to 0 in order to perform cross validation
            if(cv_iteration==1) {
                x_cv <- x
                valid_entries <- which(x_cv>0,arr.ind=TRUE)
                cv_entries <- valid_entries[sample(1:nrow(valid_entries),size=round(nrow(valid_entries)*cross_validation_entries),replace=FALSE),]
                x_cv[cv_entries] <- 0
            }
            
            # structure to save the results of the grid search
            grid_search <- array(list(),c(length(K),length(lambda_values)))
            rownames(grid_search) <- paste0(as.character(K),"_signatures")
            colnames(grid_search) <- paste0(as.character(lambda_values),"_lambda")
            
            # structure to save the starting values of beta for each K
            if(is.null(starting_beta)) {
                starting_beta <- array(list(),c(length(K),1))
                rownames(starting_beta) <- paste0(as.character(K),"_signatures")
                colnames(starting_beta) <- "Value"
            }

            # consider all the values for K
            cont <- 0
            pos_k <- 0
            for(k in K) {
                    
                # get the first k signatures to be used for the current configuration
                pos_k <- pos_k + 1
                if(is.null(starting_beta[[pos_k,1]])) {
                    if(nmf_method=="nmf_lasso") {
                        
                        if(verbose) {
                            cat("Computing the initial values of beta by NMF with Lasso...","\n")
                        }
                        
                        # add a signature to beta (leading to K+1 signatures in total) to explicitly model noise
                        if(is.null(background_signature)) {
                            warning("No background signature has been specified...")
                        }
                        
                        # compute the starting points nmf_runs times
                        beta_estimation <- list()
                        beta_mse <- NULL
                        for(i in 1:nmf_runs) {
                            
                            # set the initial random values for beta
                            if(is.null(background_signature)) {
                                curr_beta <- matrix(0,nrow=k,ncol=dim(x)[2])
                                for(i in 1:k) {
                                    curr_beta[i,] <- runif(dim(x)[2])
                                }
                                colnames(curr_beta) <- colnames(x)
                            }
                            else {
                                curr_beta <- matrix(0,nrow=(k+1),ncol=dim(x)[2])
                                curr_beta[1,] <- background_signature
                                for(i in 2:(k+1)) {
                                    curr_beta[i,] <- runif(dim(x)[2])
                                }
                                rownames(curr_beta) <- c("background_signature",rep("",k))
                                colnames(curr_beta) <- colnames(x)
                            }
                            
                            # compute the starting beta given these initial values
                            curr_starting_beta_estimation = tryCatch({
                                res <- nmfLassoDecomposition(x,curr_beta,lambda_rate=0.01,iterations=20,max_iterations_lasso=10000,parallel=parallel,verbose=FALSE)
                                mse <- sum((x-round(res$alpha%*%res$curr_beta))^2)/nrow(x)
                                list(curr_beta=res$curr_beta,mse=mse)
                            }, error = function(e) {
                                list(curr_beta=NA,mse=NA)
                            })
                            
                            # save the results at the current step if not NA
                            if(!is.na(curr_starting_beta_estimation$mse)) {
                                beta_estimation[[(length(beta_estimation)+1)]] <- curr_starting_beta_estimation$curr_beta
                                beta_mse <- c(beta_mse,curr_starting_beta_estimation$mse)
                            }
                            
                        }
                        # estimate the best starting betas
                        if(is.null(beta_mse)) {
                            stopCluster(parallel)
                            stop("Something went wrong while estimating the starting beta values, you may try again or consider to use nmf_standard initialization...")
                        }
                        else {
                            curr_beta <- beta_estimation[[which.min(beta_mse)]]
                        }
                        
                    }
                    else if(nmf_method=="nmf_standard") {
                        
                        if(verbose) {
                            cat("Computing the initial values of beta by standard NMF...","\n")
                        }
                        curr_beta <- basis(nmf(t(x),rank=k,nrun=nmf_runs))
                        curr_beta <- t(curr_beta)
                        
                        # add a signature to beta (leading to K+1 signatures in total) to explicitly model noise
                        if(is.null(background_signature)) {
                            warning("No background signature has been specified...")
                        }
                        else {
                            curr_beta <- rbind(background_signature,curr_beta)
                        }
                        
                    }
                    curr_beta <- curr_beta / rowSums(curr_beta)
                    starting_beta[[pos_k,1]] <- curr_beta
                }
                else {
                    curr_beta <- starting_beta[[pos_k,1]]
                }
                
                # consider all the values for lambda
                pos_l <- 0
                for(l in lambda_values) {
                    
                    # set the predicted values for the cross validation entries
                    pos_l <- pos_l + 1
                    if(cv_iteration>1 && !is.na(grid_search_iterations[[(cv_iteration-1)]][[pos_k,pos_l]])) {
                        best_alpha <- grid_search_iterations[[(cv_iteration-1)]][[pos_k,pos_l]][["alpha"]]
                        best_beta <- grid_search_iterations[[(cv_iteration-1)]][[pos_k,pos_l]][["beta"]]
                        predicted_counts <- best_alpha %*% best_beta
                        x_cv[cv_entries] <- predicted_counts[cv_entries]
                    }
                    
                    # perform the inference
                    curr_results <- nmf.LassoK(x = x_cv, 
                                             K = k, 
                                             beta = curr_beta, 
                                             background_signature = background_signature, 
                                             nmf_runs = 10, 
                                             lambda_rate = l, 
                                             iterations = iterations, 
                                             max_iterations_lasso = max_iterations_lasso, 
                                             num_processes = NULL, 
                                             parallel = parallel, 
                                             seed = round(runif(1)*100000), 
                                             verbose = FALSE)
                    
                    # save the results for the current configuration
                    grid_search[[pos_k,pos_l]] <- curr_results
                    
                    if(verbose) {
                        cont <- cont + 1
                        cat("Progress",paste0(round((cont/(length(K)*length(lambda_values)))*100,digits=3),"%..."),"\n")
                    }
                    
                }
                
            }
            
            # save the results for the current iteration
            grid_search_iterations[[cv_iteration]] <- grid_search

        }
        
        if(verbose) {
            cat("Evaluating the results of the cross validation in terms of mean squared error...","\n")
        }

        # structure to save the mean squared errors for all the cross_validation_iterations
        mean_squared_error_iterations <- list()

        # repeat the estimation for a number of cross_validation_iterations
        for(cv_iteration in 1:cross_validation_iterations) {
        
            # structure to save the mean squared errors
            mean_squared_error <- array(NA,c(length(K),length(lambda_values)))
            rownames(mean_squared_error) <- paste0(as.character(K),"_signatures")
            colnames(mean_squared_error) <- paste0(as.character(lambda_values),"_lambda")
            
            # assess the results of the cross validation by mean squared error
            J <- dim(x)[2]
            pos_k <- 0
            for(k in K) {
                pos_k <- pos_k + 1
                pos_l <- 0
                for(l in lambda_values) {
                    
                    # consider the current configuration
                    pos_l <- pos_l + 1
                    curr_alpha <- grid_search_iterations[[cv_iteration]][[pos_k,pos_l]][["alpha"]]
                    curr_beta <- grid_search_iterations[[cv_iteration]][[pos_k,pos_l]][["beta"]]
                    
                    # compute the mean squared error
                    if(!is.na(curr_alpha)&&!is.na(curr_beta)) {
                        curr_predicted_counts = round(curr_alpha%*%curr_beta)
                        curr_true_considered_counts = as.vector(x[cv_entries])
                        curr_predicted_considered_counts = as.vector(curr_predicted_counts[cv_entries])
                        error = mean((curr_true_considered_counts-curr_predicted_considered_counts)^2)
                    }
                    else {
                        error <- NA
                    }

                    mean_squared_error[pos_k,pos_l] <- error
                    
                }
            }

            # save the results for the current iteration
            mean_squared_error_iterations[[cv_iteration]] <- mean_squared_error

        }
        
        if(verbose) {
            cat("Estimating the best configuration...","\n")
        }

        # save the results
        curr_results <- list(grid_search=grid_search_iterations,starting_beta=starting_beta,mean_squared_error=mean_squared_error_iterations)
        results[[cv_repetitions]] <- curr_results

    }
    
    # close parallel
    if(!is.null(parallel)) {
        stopCluster(parallel)
    }
    
    return(results)
    
}

#' Perform a robust estimation of the starting beta for the nmfLasso method 
#' @title starting.betas.estimation
#' @param x count matrix.
#' @param K range of numeric values (each of them greater than 1) indicating the number of signatures to be discovered.
#' @param background_signature background signature to be used. If not provided, a warning is thrown.
#' @param nmf_runs number of iteration of NMF to be performed for a robust estimation of starting beta. If beta is not NULL, 
#' this parameter is ignored.
#' @param num_processes Number of processes to be used during parallel execution. If executing in single process mode, 
#' this is ignored.
#' @param seed Seed for reproducibility.
#' @param verbose boolean; Shall I print all messages?
#' @return A list of starting beta values for each configuration of K.
#' @export starting.betas.estimation
#' @importFrom NMF basis nmf
#' @import nnls
#' @import parallel
#'
"starting.betas.estimation" <- function( x, K = 3:10, background_signature = NULL, nmf_runs = 10, num_processes = Inf, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)

    nmf_method <- "nmf_standard"
    
    # setting up parallel execution
    close_parallel <- FALSE
    if(nmf_method=="nmf_lasso") {
        if(is.na(num_processes) || is.null(num_processes)) {
            parallel <- NULL
        }
        else if(num_processes==Inf) {
            cores <- as.integer((detectCores()-1))
            if(cores < 2) {
                parallel <- NULL
            }
            else {
                num_processes <- cores
                parallel <- makeCluster(num_processes,outfile="")
                clusterSetRNGStream(parallel,iseed<-round(runif(1)*100000))
                close_parallel <- TRUE
            }
        }
        else {
            parallel <- makeCluster(num_processes,outfile="")
            clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
            close_parallel <- TRUE
        }
        
        if(verbose && !is.null(parallel)) {
            cat("Executing",num_processes,"processes via parallel...","\n")
        }
    }
    
    # perform a robust estimation of the starting beta for the nmfLasso method
    if(verbose) {
        cat("Performing a robust estimation of the starting betas for the nmfLasso method...","\n")
    }

    # structure to save the starting values of beta for each K
    starting_beta <- array(list(),c(length(K),1))
    rownames(starting_beta) <- paste0(as.character(K),"_signatures")
    colnames(starting_beta) <- "Value"

    # consider all the values for K
    pos_k <- 0
    for(k in K) {
    
        # compute the initial values of beta
        pos_k <- pos_k + 1
        beta <- NULL
        if(is.null(beta)) {
            
            if(nmf_method=="nmf_lasso") {
                
                if(verbose) {
                    cat("Computing the initial values of beta by NMF with Lasso...","\n")
                }
                
                # add a signature to beta (leading to K+1 signatures in total) to explicitly model noise
                if(is.null(background_signature)) {
                    warning("No background signature has been specified...")
                }
                
                # compute the starting points nmf_runs times
                beta_estimation <- list()
                beta_mse <- NULL
                for(i in 1:nmf_runs) {
                    
                    # set the initial random values for beta
                    if(is.null(background_signature)) {
                        beta <- matrix(0,nrow=k,ncol=dim(x)[2])
                        for(i in 1:k) {
                            beta[i,] <- runif(dim(x)[2])
                        }
                        colnames(beta) <- colnames(x)
                    }
                    else {
                        beta <- matrix(0,nrow=(k+1),ncol=dim(x)[2])
                        beta[1,] <- background_signature
                        for(i in 2:(k+1)) {
                            beta[i,] <- runif(dim(x)[2])
                        }
                        rownames(beta) <- c("background_signature",rep("",k))
                        colnames(beta) <- colnames(x)
                    }
                    
                    # compute the starting beta given these initial values
                    curr_starting_beta_estimation = tryCatch({
                        res <- nmfLassoDecomposition(x,beta,lambda_rate=0.01,iterations=20,max_iterations_lasso=10000,parallel=parallel,verbose=FALSE)
                        mse <- sum((x-round(res$alpha%*%res$beta))^2)/nrow(x)
                        list(beta=res$beta,mse=mse)
                    }, error = function(e) {
                        list(beta=NA,mse=NA)
                    })
                    
                    # save the results at the current step if not NA
                    if(!is.na(curr_starting_beta_estimation$mse)) {
                        beta_estimation[[(length(beta_estimation)+1)]] <- curr_starting_beta_estimation$beta
                        beta_mse <- c(beta_mse,curr_starting_beta_estimation$mse)
                    }
                    
                }
                # estimate the best starting betas
                if(is.null(beta_mse)) {
                    if(close_parallel) {
                        stopCluster(parallel)
                    }
                    stop("Something went wrong while estimating the starting beta values, you may try again or consider to use nmf_standard initialization...")
                }
                else {
                    beta <- beta_estimation[[which.min(beta_mse)]]
                }
                
            }
            else if(nmf_method=="nmf_standard") {
                
                if(verbose) {
                    cat("Computing the initial values of beta by standard NMF...","\n")
                }
                beta <- basis(nmf(t(x),rank=k,nrun=nmf_runs))
                beta <- t(beta)
                
                # add a signature to beta (leading to K+1 signatures in total) to explicitly model noise
                if(is.null(background_signature)) {
                    warning("No background signature has been specified...")
                }
                else {
                    beta <- rbind(background_signature,beta)
                }
                
            }

        }
        beta <- beta / rowSums(beta)
        starting_beta[[pos_k,1]] <- beta

        if(verbose) {
            cat("Progress",paste0(round((pos_k/length(K))*100,digits=3),"%..."),"\n")
        }

    }
    
    # close parallel
    if(close_parallel) {
        stopCluster(parallel)
    }

    return(starting_beta)
    
}

#' Estimate the range of lambda values to be considered in the signature inference. Note that too small values of lambda 
#' result in dense signatures, but too large values lead to bad fit of the counts.
#' @title evaluate.lambda.range
#' @param x count matrix.
#' @param K numeric value (greater than 1) indicating the number of signatures to be discovered.
#' @param beta starting beta for the estimation. If it is NULL, starting beta is estimated by NMF.
#' @param background_signature background signature to be used. If not provided, a warning is thrown.
#' @param nmf_runs number of iteration of NMF to be performed for a robust estimation of starting beta. If beta is not NULL, 
#' this parameter is ignored.
#' @param lambda_values range of values of LASSO to be used between 0 and 1. This value should be greater than 0. 1 is the value of LASSO 
#' that would shrink all the signatures to 0 within one step. The higher lambda_rate is, the sparser are the resulting signatures, 
#' but too large values result in a poor fit of the counts.
#' @param iterations Number of iterations to be performed. Each iteration correspond to a first step where the counts are fitted 
#' and a second step where sparsity is enhanced.
#' @param max_iterations_lasso Number of maximum iterations to be performed during the sparsification.
#' @param num_processes Number of processes to be used during parallel execution. If executing in single process mode, 
#' this is ignored.
#' @param seed Seed for reproducibility.
#' @param verbose boolean; Shall I print all messages?
#' @return A list corresponding to results of the function nmf.LassoK for each value of lambda to be tested. This function allows 
#' to test a good range of lambda values to be considered. One should keep in mind that too small values generate dense solution, 
#' while too high ones leads to poor fit. This behavior is resampled in the values of loglik_progression, which should be increasing: 
#' too small values of lamda results in unstable log-likelihood through the iterations, while too large values make log-likelihood 
#' drop. 
#' @export evaluate.lambda.range
#' @import NMF
#' @import nnlasso
#' @import nnls
#' @import parallel
#'
"evaluate.lambda.range" <- function( x, K = 6, beta = NULL, background_signature = NULL, nmf_runs = 10, lambda_values = c(0.10, 0.20, 0.30, 0.40, 0.50), iterations = 20, max_iterations_lasso = 10000, num_processes = Inf, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)

    nmf_method <- "nmf_standard"
    
    # setting up parallel execution
    parallel <- NULL
    close_parallel <- FALSE
    if(is.null(parallel)) {
        if(is.na(num_processes) || is.null(num_processes)) {
            parallel <- NULL
        }
        else if(num_processes==Inf) {
            cores <- as.integer((detectCores()-1))
            if(cores < 2) {
                parallel <- NULL
            }
            else {
                num_processes <- cores
                parallel <- makeCluster(num_processes,outfile="")
                clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
                close_parallel <- TRUE
            }
        }
        else {
            parallel <- makeCluster(num_processes,outfile="")
            clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
            close_parallel <- TRUE
        }
        
        if(verbose && !is.null(parallel)) {
            cat("Executing",num_processes,"processes via parallel...","\n")
        }
    }
    
    # compute the initial values of beta
    if(is.null(beta)) {
        if(nmf_method=="nmf_lasso") {
            
            if(verbose) {
                cat("Computing the initial values of beta by NMF with Lasso...","\n")
            }
            
            # add a signature to beta (leading to K+1 signatures in total) to explicitly model noise
            if(is.null(background_signature)) {
                warning("No background signature has been specified...")
            }
            
            # compute the starting points nmf_runs times
            beta_estimation <- list()
            beta_mse <- NULL
            for(i in 1:nmf_runs) {
                
                # set the initial random values for beta
                if(is.null(background_signature)) {
                    beta <- matrix(0,nrow=K,ncol=dim(x)[2])
                    for(i in 1:K) {
                        beta[i,] <- runif(dim(x)[2])
                    }
                    colnames(beta) <- colnames(x)
                }
                else {
                    beta <- matrix(0,nrow=(K+1),ncol=dim(x)[2])
                    beta[1,] <- background_signature
                    for(i in 2:(K+1)) {
                        beta[i,] <- runif(dim(x)[2])
                    }
                    rownames(beta) <- c("background_signature",rep("",K))
                    colnames(beta) <- colnames(x)
                }
                
                # compute the starting beta given these initial values
                curr_starting_beta_estimation = tryCatch({
                    res <- nmfLassoDecomposition(x,beta,lambda_rate=0.01,iterations=20,max_iterations_lasso=10000,parallel=parallel,verbose=FALSE)
                    mse <- sum((x-round(res$alpha%*%res$beta))^2)/nrow(x)
                    list(beta=res$beta,mse=mse)
                }, error = function(e) {
                    list(beta=NA,mse=NA)
                })
                
                # save the results at the current step if not NA
                if(!is.na(curr_starting_beta_estimation$mse)) {
                    beta_estimation[[(length(beta_estimation)+1)]] <- curr_starting_beta_estimation$beta
                    beta_mse <- c(beta_mse,curr_starting_beta_estimation$mse)
                }
                
            }
            # estimate the best starting betas
            if(is.null(beta_mse)) {
                if(close_parallel) {
                    stopCluster(parallel)
                }
                stop("Something went wrong while estimating the starting beta values, you may try again or consider to use nmf_standard initialization...")
            }
            else {
                beta <- beta_estimation[[which.min(beta_mse)]]
            }
            
        }
        else if(nmf_method=="nmf_standard") {
            
            if(verbose) {
                cat("Computing the initial values of beta by standard NMF...","\n")
            }
            beta <- basis(nmf(t(x),rank=K,nrun=nmf_runs))
            beta <- t(beta)
            
            # add a signature to beta (leading to K+1 signatures in total) to explicitly model noise
            if(is.null(background_signature)) {
                warning("No background signature has been specified...")
            }
            else {
                beta <- rbind(background_signature,beta)
            }
            
        }
    }
    
    # structure to save the estimated signatures
    lambda_results <- array(list(),c(length(K),length(lambda_values)))
    rownames(lambda_results) <- paste0(as.character(K),"_signatures")
    colnames(lambda_results) <- paste0(as.character(lambda_values),"_lambda")
    
    # perform signature discovery for all the values of lambda
    cont <- 0
    for(l in lambda_values) {
            
        # perform the inference
        curr_results <- nmf.LassoK(x = x, 
                                 K = K, 
                                 beta = beta, 
                                 background_signature = background_signature, 
                                 nmf_runs = 10, 
                                 lambda_rate = l, 
                                 iterations = iterations, 
                                 max_iterations_lasso = max_iterations_lasso, 
                                 num_processes = NULL, 
                                 parallel = parallel, 
                                 seed = round(runif(1)*100000), 
                                 verbose = FALSE)
                                 
        # save the results for the current configuration
        cont <- cont + 1
        lambda_results[[1,cont]] <- curr_results
            
        if(verbose) {
            cat("Progress",paste0(round((cont/(length(K)*length(lambda_values)))*100,digits=3),"%..."),"\n")
        }
        
    }
    
    # close parallel
    if(close_parallel) {
        stopCluster(parallel)
    }
    
    # save the results
    results <- lambda_results
    
    return(results)
    
}

#' Perform the discovery of K somatic mutational signatures given a set of observed counts x.
#'
#' @examples
#' data(starting_betas_example)
#' beta = starting_betas_example[["5_signatures","Value"]]
#' res = nmf.LassoK(x=patients,K=5,beta=beta,background=background,lambda_rate=0.10,iterations=5,num_processes=NA)
#'
#' @title nmf.LassoK
#' @param x count matrix.
#' @param K numeric value (greater than 1) indicating the number of signatures to be discovered.
#' @param beta starting beta for the estimation. If it is NULL, starting beta is estimated by NMF.
#' @param background_signature background signature to be used. If not provided, a warning is thrown.
#' @param nmf_runs number of iteration of NMF to be performed for a robust estimation of starting beta. If beta is not NULL, 
#' this parameter is ignored.
#' @param lambda_rate value of LASSO to be used between 0 and 1. This value should be greater than 0. 1 is the value of LASSO 
#' that would shrink all the signatures to 0 within one step. The higher lambda_rate is, the sparser are the resulting signatures, 
#' but too large values result in a poor fit of the counts.
#' @param iterations Number of iterations to be performed. Each iteration correspond to a first step where the counts are fitted 
#' and a second step where sparsity is enhanced.
#' @param max_iterations_lasso Number of maximum iterations to be performed during the sparsification.
#' @param num_processes Number of processes to be used during parallel execution. If executing in single process mode, 
#' this is ignored.
#' @param parallel Cluster object for parallel execution.
#' @param seed Seed for reproducibility.
#' @param verbose boolean; Shall I print all messages?
#' @return A list with the discovered signatures. It includes 5 elements: 
#'              alpha: matrix of the discovered alpha values
#'              beta: matrix of the discovered signatures
#'              starting_beta: initial signatures on which the method has been apploid
#'              best_loglik: log-likelihood of the best signatures configuration
#'              loglik_progression: log-likelihood values during the iterations. This values should be increasing, if not the selected value of lambda is too high
#' @export nmf.LassoK
#' @import NMF
#' @import nnlasso
#' @import nnls
#' @import parallel
#'
"nmf.LassoK" <- function( x, K, beta = NULL, background_signature = NULL, nmf_runs = 10, lambda_rate = 0.20, iterations = 20, max_iterations_lasso = 10000, num_processes = Inf, parallel = NULL, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)

    nmf_method <- "nmf_standard"
    
    # setting up parallel execution
    close_parallel <- FALSE
    if(is.null(parallel)) {
        if(is.na(num_processes) || is.null(num_processes)) {
            parallel <- NULL
        }
        else if(num_processes==Inf) {
            cores <- as.integer((detectCores()-1))
            if(cores < 2) {
                parallel <- NULL
            }
            else {
                num_processes <- cores
                parallel <- makeCluster(num_processes,outfile="")
                clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
                close_parallel <- TRUE
            }
        }
        else {
            parallel <- makeCluster(num_processes,outfile="")
            clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
            close_parallel <- TRUE
        }
        
        if(verbose && !is.null(parallel)) {
            cat("Executing",num_processes,"processes via parallel...","\n")
        }
    }
    
    # compute the initial values of beta
    if(is.null(beta)) {
        if(nmf_method=="nmf_lasso") {
            
            if(verbose) {
                cat("Computing the initial values of beta by NMF with Lasso...","\n")
            }
            
            # add a signature to beta (leading to K+1 signatures in total) to explicitly model noise
            if(is.null(background_signature)) {
                warning("No background signature has been specified...")
            }
            
            # compute the starting points nmf_runs times
            beta_estimation <- list()
            beta_mse <- NULL
            for(i in 1:nmf_runs) {
                
                # set the initial random values for beta
                if(is.null(background_signature)) {
                    beta <- matrix(0,nrow=K,ncol=dim(x)[2])
                    for(i in 1:K) {
                        beta[i,] <- runif(dim(x)[2])
                    }
                    colnames(beta) <- colnames(x)
                }
                else {
                    beta <- matrix(0,nrow=(K+1),ncol=dim(x)[2])
                    beta[1,] <- background_signature
                    for(i in 2:(K+1)) {
                        beta[i,] <- runif(dim(x)[2])
                    }
                    rownames(beta) <- c("background_signature",rep("",K))
                    colnames(beta) <- colnames(x)
                }
                
                # compute the starting beta given these initial values
                curr_starting_beta_estimation = tryCatch({
                    res <- nmfLassoDecomposition(x,beta,lambda_rate=0.01,iterations=20,max_iterations_lasso=10000,parallel=parallel,verbose=FALSE)
                    mse <- sum((x-round(res$alpha%*%res$beta))^2)/nrow(x)
                    list(beta=res$beta,mse=mse)
                }, error = function(e) {
                    list(beta=NA,mse=NA)
                })
                
                # save the results at the current step if not NA
                if(!is.na(curr_starting_beta_estimation$mse)) {
                    beta_estimation[[(length(beta_estimation)+1)]] <- curr_starting_beta_estimation$beta
                    beta_mse <- c(beta_mse,curr_starting_beta_estimation$mse)
                }
                
            }
            # estimate the best starting betas
            if(is.null(beta_mse)) {
                if(close_parallel) {
                    stopCluster(parallel)
                }
                stop("Something went wrong while estimating the starting beta values, you may try again or consider to use nmf_standard initialization...")
            }
            else {
                beta <- beta_estimation[[which.min(beta_mse)]]
            }
            
        }
        else if(nmf_method=="nmf_standard") {
            
            if(verbose) {
                cat("Computing the initial values of beta by standard NMF...","\n")
            }
            beta <- basis(nmf(t(x),rank=K,nrun=nmf_runs))
            beta <- t(beta)
            
            # add a signature to beta (leading to K+1 signatures in total) to explicitly model noise
            if(is.null(background_signature)) {
                warning("No background signature has been specified...")
            }
            else {
                beta <- rbind(background_signature,beta)
            }
            
        }
    }
    
    if(verbose) {
        cat("Performing the discovery of the signatures by NMF with Lasso...","\n")
    }
    
    # perform the discovery of the signatures
    results = tryCatch({
        nmfLassoDecomposition(x,beta,lambda_rate,iterations,max_iterations_lasso,parallel,verbose)
    }, error = function(e) {
        warning("Lasso did not converge, you should try a lower value of lambda! Current settings: K = ",K,", lambda_rate = ",lambda_rate,"...")
        list(alpha=NA,beta=NA,starting_beta=NA,best_loglik=NA,loglik_progression=rep(NA,iterations))
    })
    
    # close parallel
    if(close_parallel) {
        stopCluster(parallel)
    }
    
    return(results)
    
}

# perform de novo discovery of somatic mutational signatures using NMF with Lasso to ensure sparsity
"nmfLassoDecomposition" <- function( x, beta, lambda_rate = 0.20, iterations = 20, max_iterations_lasso = 10000, parallel = NULL, verbose = TRUE ) {
    
    # n is the number of observations in x, i.e., the patients
    n <- dim(x)[1]
    
    # J is the number of trinucleotides, i.e., 96 categories
    J <- dim(x)[2]
    
    # K is the number of signatures to be called
    K <- dim(beta)[1]
    
    # initialize alpha
    alpha <- matrix(0,nrow=n,ncol=K)
    rownames(alpha) <- 1:nrow(alpha)
    colnames(alpha) <- rownames(beta)
    
    # starting beta
    starting_beta <- beta / rowSums(beta)
    beta <- starting_beta * 1000
    
    # structure where to save the log-likelihood at each iteration 
    loglik <- rep(NA,iterations)
    
    # structure where to save the lambda values for each J
    lambda_values <- rep(NA,J)
    
    if(verbose) {
        cat("Performing a total of",iterations,"iterations...","\n")
    }
    
    # repeat a 2 step algorithm iteratively, where first alpha is estimated by Non-Negative Linear Least Squares 
    # and, then, beta is estimated by Non-Negative Lasso
    best_loglik <- -Inf
    best_alpha <- NA
    best_beta <- NA
    for(i in 1:iterations) {
        
        # initialize the value of the log-likelihood for the current iteration
        loglik[i] <- 0
        
        # update alpha independently for each patient by Non-Negative Linear Least Squares
        if(is.null(parallel)) {
            for(j in 1:n) {
                alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
            }
        }
        else {
            
            # compute alpha in parallel
            j <- 1:n
            res_clusterEvalQ <- clusterEvalQ(parallel,library("nnls"))
            clusterExport(parallel,varlist=c("x","beta"),envir=environment())
            alpha_res <- parLapply(parallel,j,function(j) {
                return(nnls(t(beta),as.vector(x[j,]))$x)
            })
            
            # reduce the results
            for(j in 1:n) {
                alpha[j,] <- alpha_res[[j]]
            }
            
        }
        
        # update beta by Non-Negative Lasso
        if(is.null(parallel)) {
            
            for(k in 1:J) {
                
                # compute independently for each trinucleotide the error between the observed counts, i.e., x, 
                # and the ones predicted by the first signature (i.e., which represents the noise model)
                error <- x[,k] - alpha[,1] * beta[1,k]
                
                # estimate beta for the remaining signatues by Non-Negative Lasso to ensure sparsity
                if(is.na(lambda_values[k])) {
                    max_lambda_value <- max(abs(t(alpha[,2:K]) %*% error))
                    lambda_values[k] <- max_lambda_value * lambda_rate
                }
                beta[2:K,k] <- as.vector(nnlasso(x = alpha[,2:K], 
                                                y = error, 
                                                family = "normal", 
                                                lambda = lambda_values[k], 
                                                intercept = FALSE, 
                                                normalize = FALSE, 
                                                maxiter = max_iterations_lasso, 
                                                path = FALSE)$coef[2,])
                
                # update the log-likelihood for the current iteration
                curr_loglik <- -sum((x[,k] - alpha %*% beta[,k])^2) - lambda_values[k] * sum(beta[2:K,k])
                loglik[i] <- loglik[i] + curr_loglik
                
            }
            
        }
        else {
            
            # if this is the first iteration, compute the values of lambda sequentially
            if(is.na(lambda_values[1])) {
                
                for(k in 1:J) {
                    error <- x[,k] - alpha[,1] * beta[1,k]
                    max_lambda_value <- max(abs(t(alpha[,2:K]) %*% error))
                    lambda_values[k] <- max_lambda_value * lambda_rate
                }
                
            }
            
            # perform the computations of beta in parallel
            k <- 1:J
            res_clusterEvalQ <- clusterEvalQ(parallel,library("nnlasso"))
            clusterExport(parallel,varlist=c("x","alpha","beta","lambda_values","K","max_iterations_lasso"),envir<-environment())
            beta_res <- parLapply(parallel,1:J,function(k) {
                
                # compute independently for each trinucleotide the error between the observed counts, i.e., x, 
                # and the ones predicted by the first signature (i.e., which represents the noise model)
                error <- x[,k] - alpha[,1] * beta[1,k]
                
                # estimate beta for the remaining signatues by Non-Negative Lasso to ensure sparsity
                beta[2:K,k] <- as.vector(nnlasso(x = alpha[,2:K], 
                                                y = error, 
                                                family = "normal", 
                                                lambda = lambda_values[k], 
                                                intercept = FALSE, 
                                                normalize = FALSE, 
                                                maxiter = max_iterations_lasso, 
                                                path = FALSE)$coef[2,])
                
                # compute the log-likelihood for the current iteration
                curr_loglik <- -sum((x[,k] - alpha %*% beta[,k])^2) - lambda_values[k] * sum(beta[2:K,k])
                
                return(list(beta=beta[2:K,k],loglik=curr_loglik))
                
            })
            
            # reduce the results
            for(k in 1:J) {
                
                beta[2:K,k] <- beta_res[[k]][["beta"]]
                loglik[i] <- loglik[i] + beta_res[[k]][["loglik"]]
                
            }
            
        }
        
        # save the grid_search at maximum log-likelihood
        if(loglik[i]>best_loglik) {
            best_loglik <- loglik[i]
            best_alpha <- alpha
            best_beta <- beta
        }
        
        if(verbose) {
            cat("Progress",paste0((i/iterations)*100,"%..."),"\n")
        }
        
    }
    alpha <- best_alpha
    beta <- best_beta
    
    # check if the likelihood is increasing
    if(length(loglik)>1) {
        cont <- 1
        for(i in 2:length(loglik)) {
            if(loglik[i]>loglik[(i-1)]) {
                cont <- cont + 1
            }
        }
        if(cont/iterations<0.5) {
            warning("The likelihood is not increasing, you should try a lower value of lambda! Current settings: K = ",K,", lambda_rate = ",lambda_rate,"...")
        }
    }
    
    # normalize the rate of the signatures to sum to 1
    beta <- beta / rowSums(beta)
        
    # final computation of alpha by Non-Negative Linear Least Squares using the normalized version o the best beta
    if(is.null(parallel)) {
        for(j in 1:n) {
            alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
        }
    }
    else {
        
        # compute alpha in parallel
        j <- 1:n
        res_clusterEvalQ <- clusterEvalQ(parallel,library("nnls"))
        clusterExport(parallel,varlist=c("x","beta"),envir=environment())
        alpha_res <- parLapply(parallel,j,function(j) {
            return(nnls(t(beta),as.vector(x[j,]))$x)
        })
        
        # reduce the results
        for(j in 1:n) {
            alpha[j,] <- alpha_res[[j]]
        }
        
    }
    
    # save the results
    results <- list(alpha=alpha,beta=beta,starting_beta=starting_beta,best_loglik=best_loglik,loglik_progression=loglik)
    
    return(results)
    
}
