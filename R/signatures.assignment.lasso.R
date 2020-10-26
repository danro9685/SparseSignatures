#' Perform the assessment of different sigAssignmentLasso solutions by cross validation for a given set of somatic mutational signatures and observations x. The estimation can slow down because of 
#' memory usage and intensive computations, when a big number of cross validation repetitions is asked and when the grid search is performed for a lot of configurations. In this case, 
#' an advice may be to split the computation into multiple smaller sets.
#'
#' @examples
#' data(patients)
#' data(starting_betas_example)
#' beta = starting_betas_example[["5_signatures","Value"]]
#' res = sigAssignmentCV(x=patients[1:100,],
#'      beta=beta,
#'      lambda_values_alpha=c(0.00),
#'      cross_validation_repetitions=2,
#'      num_processes=NA,
#'      seed=12345)
#'
#' @title sigAssignmentCV
#' @param x count matrix for a set of n patients and 96 trinucleotides.
#' @param beta beta to be fixed during the estimation of alpha.
#' @param normalize_counts if true, the input count matrix x is normalize such that the patients have the same number of mutation.
#' @param lambda_values_alpha value of LASSO to be used for alpha between 0 and 1. This value should be greater than 0. 1 is the value of LASSO 
#' that would shrink all the signatures to 0 within one step. The higher lambda_rate_alpha is, the sparser are the resulting exposures, 
#' but too large values may result in a reduced fit of the observed counts.
#' @param cross_validation_entries Percentage of cells in the count matrix to be replaced by 0s during cross validation.
#' @param cross_validation_iterations For each configuration, the first time the signatures are discovered form a matrix with a 
#' percentage of values replaced by 0s. This may result in poor fit/results. Then, we perform predictions of these entries and replace them with 
#' such predicted values. This parameter is the number of restarts to be performed to improve this estimate and obtain more stable solutions.
#' @param cross_validation_repetitions Number of time cross-validation should be repeated. Higher values result in better estimate, but are computationally more expensive.
#' @param max_iterations_lasso Number of maximum iterations to be performed during the sparsification via Lasso.
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode, 
#' this parameter needs to be set to either NA or NULL.
#' @param seed Seed for reproducibility.
#' @param verbose boolean; Shall I print all messages?
#' @param log_file log file where to print outputs when using parallel. If parallel execution is disabled, this parameter is ignored.
#' @return A list of 2 elements: grid_search_mse and and grid_search_loglik. Here, grid_search_mse reports the mean squared error for each configuration of performed 
#' cross validation; grid_search_loglik reports for each configuration the number of times the algorithm converged.
#' @export sigAssignmentCV
#' @import nnlasso
#' @import nnls
#' @import parallel
#'
"sigAssignmentCV" <- function( x, beta, normalize_counts = TRUE, lambda_values_alpha = c(0.00, 0.01, 0.05, 0.10), cross_validation_entries = 0.01, cross_validation_iterations = 5, cross_validation_repetitions = 50, max_iterations_lasso = 10000, num_processes = Inf, seed = NULL, verbose = TRUE, log_file = "" ) {
    
    # set the seed
    set.seed(seed)

    # check the input parameters
    x <- as.matrix(x)
    if(normalize_counts) {
        x <- (x/rowSums(x))*2500
    }
    lambda_values_alpha <- sort(unique(lambda_values_alpha))

    # initialize beta
    beta <- beta / rowSums(beta)
    rownames(beta) <- c("Background",paste0("S",1:(nrow(beta)-1)))
    colnames(beta) <- colnames(x)

    # perform a grid search to estimate the best values of K and lambda
    if(verbose) {
        cat("Performing a search to estimate the best values of lambda with a total of",cross_validation_repetitions,"cross validation repetitions...","\n")
    }

    # setting up parallel execution
    parallel <- NULL
    close_parallel <- FALSE
    if(is.null(parallel)&&cross_validation_repetitions>1) {
        if(is.na(num_processes) || is.null(num_processes)) {
            parallel <- NULL
        }
        else if(num_processes==Inf) {
            cores <- as.integer((detectCores()-1))
            if(cores < 2) {
                parallel <- NULL
            }
            else {
                num_processes <- min(num_processes,cross_validation_repetitions)
                parallel <- makeCluster(num_processes,outfile=log_file)
                close_parallel <- TRUE
            }
        }
        else {
            num_processes <- min(num_processes,cross_validation_repetitions)
            parallel <- makeCluster(num_processes,outfile=log_file)
            close_parallel <- TRUE
        }
        
        if(verbose && !is.null(parallel)) {
            cat("Executing",num_processes,"processes via parallel...","\n")
        }
    }

    # structure to save the results of the grid search
    grid_search_mse <- array(list(),c(length(lambda_values_alpha),1,1),dimnames=list(paste0(as.character(lambda_values_alpha),"_Lambda_Alpha"),"Fixed_Lambda_Beta",paste0(as.character(dim(beta)[1]),"_Signatures")))
    grid_search_loglik <- array(list(),c(length(lambda_values_alpha),1,1),dimnames=list(paste0(as.character(lambda_values_alpha),"_Lambda_Alpha"),"Fixed_Lambda_Beta",paste0(as.character(dim(beta)[1]),"_Signatures")))

    # now starting cross validation
    if(is.null(parallel)) {

        # perform a total of cross_validation_repetitions repetitions of cross validation
        set.seed(round(runif(1)*100000))
        for(cv_repetitions in 1:cross_validation_repetitions) {

            if(verbose) {
                cat(paste0("Performing repetition ",cv_repetitions," out of ",cross_validation_repetitions,"..."),"\n")
            }

            # randomly set the cross validation entries for the current iteration
            valid_entries <- which(x>0,arr.ind=TRUE)
            cv_entries <- valid_entries[sample(1:nrow(valid_entries),size=round(nrow(valid_entries)*cross_validation_entries),replace=FALSE),]

            # consider all the values for K
            cont <- 0
            pos_k <- 0
            for(k in 1:1) {
                # consider the first k signatures to be used for the current configuration
                pos_k <- pos_k + 1

                # consider all the values for lambda
                pos_l_alpha <- 0
                for(l_alpha in lambda_values_alpha) {
                    pos_l_alpha <- pos_l_alpha + 1

                    pos_l_beta <- 0
                    for(l_beta in 1:1) {
                        pos_l_beta <- pos_l_beta + 1

                        # repeat the estimation for a number of cross_validation_iterations
                        x_cv <- NULL
                        cont <- cont + 1
                        for(cv_iteration in 1:cross_validation_iterations) {
                    
                            if(verbose) {
                                cat(paste0("Performing cross validation iteration ",cv_iteration," out of ",cross_validation_iterations,"..."),"\n")
                            }

                            # set a percentage of cross_validation_entries entries to 0 in order to perform cross validation
                            if(cv_iteration==1) {
                                x_cv <- x
                                x_cv[cv_entries] <- 0
                            }
                            else {
                                if(!is.na(curr_iter_alpha[1])&&!is.na(curr_iter_beta[1])) {
                                    predicted_counts <- curr_iter_alpha %*% curr_iter_beta
                                    x_cv[cv_entries] <- predicted_counts[cv_entries]
                                }
                            }

                            # perform the inference
                            curr_results <- sigAssignmentLasso(x = x_cv, 
                                                                beta = beta, 
                                                                normalize_counts = normalize_counts, 
                                                                lambda_rate_alpha = l_alpha, 
                                                                max_iterations_lasso = max_iterations_lasso, 
                                                                seed = round(runif(1)*100000), 
                                                                verbose = FALSE)

                            curr_iter_alpha <- curr_results[["alpha"]]
                            curr_iter_beta <- curr_results[["beta"]]

                        }
                        
                        # save an estimate of mean squared error and stability of the solution
                        if(!is.na(curr_iter_alpha[1])&&!is.na(curr_iter_beta[1])) {
                            curr_predicted_counts <- round(curr_iter_alpha%*%curr_iter_beta)
                            curr_true_considered_counts <- as.vector(x[cv_entries])
                            curr_predicted_considered_counts <- as.vector(curr_predicted_counts[cv_entries])
                            error <- mean((curr_true_considered_counts-curr_predicted_considered_counts)^2)
                            grid_search_mse[pos_l_alpha,pos_l_beta,pos_k][[1]] <- c(grid_search_mse[pos_l_alpha,pos_l_beta,pos_k][[1]],error)
                            grid_search_loglik[pos_l_alpha,pos_l_beta,pos_k][[1]] <- c(grid_search_loglik[pos_l_alpha,pos_l_beta,pos_k][[1]],"increasing")
                        }
                        else {
                            grid_search_loglik[pos_l_alpha,pos_l_beta,pos_k][[1]] <- c(grid_search_loglik[pos_l_alpha,pos_l_beta,pos_k][[1]],"decreasing")
                        }
                        
                        if(verbose) {
                            cat("Progress",paste0(round((cont/(1*length(lambda_values_alpha)*1))*100,digits=3),"%..."),"\n")
                        }

                    }

                }

            }

        }

        # save the results
        results <- list(grid_search_mse=grid_search_mse,grid_search_loglik=grid_search_loglik)

    }
    else {

        # perform the inference
        res_clusterEvalQ <- clusterEvalQ(parallel,library("nnls"))
        res_clusterEvalQ <- clusterEvalQ(parallel,library("nnlasso"))
        clusterExport(parallel,varlist=c("x","beta","normalize_counts","cross_validation_entries"),envir=environment())
        clusterExport(parallel,varlist=c("lambda_values_alpha","cross_validation_iterations","max_iterations_lasso"),envir=environment())
        clusterExport(parallel,varlist=c("grid_search_mse","grid_search_loglik","cross_validation_repetitions","verbose"),envir=environment())
        clusterExport(parallel,c('sigAssignmentLasso','nmfLassoDeconstruction'),envir=environment())
        clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
        curr_results <- parLapply(parallel,1:cross_validation_repetitions,function(cv_rep) {

            if(verbose) {
                cat(paste0("Performing repetition ",cv_rep," out of ",cross_validation_repetitions,"..."),"\n")
            }

            # randomly set the cross validation entries for the current iteration
            valid_entries <- which(x>0,arr.ind=TRUE)
            cv_entries <- valid_entries[sample(1:nrow(valid_entries),size=round(nrow(valid_entries)*cross_validation_entries),replace=FALSE),]

            # consider all the values for K
            cont <- 0
            pos_k <- 0
            for(k in 1:1) {
                # consider the first k signatures to be used for the current configuration
                pos_k <- pos_k + 1

                # consider all the values for lambda
                pos_l_alpha <- 0
                for(l_alpha in lambda_values_alpha) {
                    pos_l_alpha <- pos_l_alpha + 1

                    pos_l_beta <- 0
                    for(l_beta in 1:1) {
                        pos_l_beta <- pos_l_beta + 1

                        # repeat the estimation for a number of cross_validation_iterations
                        x_cv <- NULL
                        cont <- cont + 1
                        for(cv_iteration in 1:cross_validation_iterations) {
                    
                            # set a percentage of cross_validation_entries entries to 0 in order to perform cross validation
                            if(cv_iteration==1) {
                                x_cv <- x
                                x_cv[cv_entries] <- 0
                            }
                            else {
                                if(!is.na(curr_iter_alpha[1])&&!is.na(curr_iter_beta[1])) {
                                    predicted_counts <- curr_iter_alpha %*% curr_iter_beta
                                    x_cv[cv_entries] <- predicted_counts[cv_entries]
                                }
                            }

                            # perform the inference
                            curr_results <- sigAssignmentLasso(x = x_cv, 
                                                                beta = beta, 
                                                                normalize_counts = normalize_counts, 
                                                                lambda_rate_alpha = l_alpha, 
                                                                max_iterations_lasso = max_iterations_lasso, 
                                                                seed = round(runif(1)*100000), 
                                                                verbose = FALSE)

                            curr_iter_alpha <- curr_results[["alpha"]]
                            curr_iter_beta <- curr_results[["beta"]]

                        }
                        
                        # save an estimate of mean squared error and stability of the solution
                        if(!is.na(curr_iter_alpha[1])&&!is.na(curr_iter_beta[1])) {
                            curr_predicted_counts <- round(curr_iter_alpha%*%curr_iter_beta)
                            curr_true_considered_counts <- as.vector(x[cv_entries])
                            curr_predicted_considered_counts <- as.vector(curr_predicted_counts[cv_entries])
                            error <- mean((curr_true_considered_counts-curr_predicted_considered_counts)^2)
                            grid_search_mse[pos_l_alpha,pos_l_beta,pos_k][[1]] <- c(grid_search_mse[pos_l_alpha,pos_l_beta,pos_k][[1]],error)
                            grid_search_loglik[pos_l_alpha,pos_l_beta,pos_k][[1]] <- c(grid_search_loglik[pos_l_alpha,pos_l_beta,pos_k][[1]],"increasing")
                        }
                        else {
                            grid_search_loglik[pos_l_alpha,pos_l_beta,pos_k][[1]] <- c(grid_search_loglik[pos_l_alpha,pos_l_beta,pos_k][[1]],"decreasing")
                        }

                    }

                }

            }

            # save the results
            curr_results <- list(grid_search_mse=grid_search_mse,grid_search_loglik=grid_search_loglik)

            return(curr_results)

        })

        # save the results
        results <- list(grid_search_mse=grid_search_mse,grid_search_loglik=grid_search_loglik)
        for(par_res in 1:length(curr_results)) {
            curr_par_res <- curr_results[[par_res]]
            curr_par_res_grid_search_mse <- curr_par_res$grid_search_mse
            curr_par_res_grid_search_loglik <- curr_par_res$grid_search_loglik
            for(a in 1:dim(curr_par_res_grid_search_mse)[1]) {
                for(b in 1:dim(curr_par_res_grid_search_mse)[2]) {
                    for(c in 1:dim(curr_par_res_grid_search_mse)[3]) {
                        if(!is.null(unlist(curr_par_res_grid_search_mse[a,b,c]))) {
                            results$grid_search_mse[[a,b,c]] <- c(unlist(results$grid_search_mse[a,b,c]),unlist(curr_par_res_grid_search_mse[a,b,c]))
                            results$grid_search_loglik[[a,b,c]] <- c(unlist(results$grid_search_loglik[a,b,c]),unlist(curr_par_res_grid_search_loglik[a,b,c]))
                        }
                    }
                }
            }
        }

    }
    
    # close parallel
    if(close_parallel) {
        stopCluster(parallel)
    }
    
    return(results)
    
}

#' Estimate the range of lambda values for alpha to be considered in the signature assignment. Note that too small values of lambda 
#' result in dense exposures, but too large values lead to bad fit of the counts.
#'
#' @examples
#' data(patients)
#' data(starting_betas_example)
#' beta = starting_betas_example[["5_signatures","Value"]]
#' res = sigAssignmentEvaluation(x=patients[1:100,],
#'      beta=beta,
#'      lambda_values=c(0.01,0.05),
#'      num_processes=NA,
#'      seed=12345)
#'
#' @title sigAssignmentEvaluation
#' @param x count matrix for a set of n patients and 96 trinucleotides.
#' @param beta beta to be fixed during the estimation of alpha.
#' @param normalize_counts if true, the input count matrix x is normalize such that the patients have the same number of mutation.
#' @param lambda_values value of LASSO to be used for alpha between 0 and 1. This value should be greater than 0. 1 is the value of LASSO 
#' that would shrink all the signatures to 0 within one step. The higher lambda_values is, the sparser are the resulting exposures, 
#' but too large values may result in a reduced fit of the observed counts.
#' @param max_iterations_lasso Number of maximum iterations to be performed during the sparsification via Lasso.
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode, 
#' this parameter needs to be set to either NA or NULL.
#' @param seed Seed for reproducibility.
#' @param verbose boolean; Shall I print all messages?
#' @param log_file log file where to print outputs when using parallel. If parallel execution is disabled, this parameter is ignored.
#' @return A list corresponding to results of the function sigAssignmentLasso for each value of lambda to be tested. This function allows 
#' to test a good range of lambda values for alpha to be considered. One should keep in mind that too small values generate dense solution, 
#' while too high ones leads to poor fit. This behavior is resampled in the values of loglik_progression, which should be increasing: 
#' too small values of lambda results in unstable log-likelihood through the iterations, while too large values make log-likelihood 
#' drop. 
#' @export sigAssignmentEvaluation
#' @import nnlasso
#' @import nnls
#' @import parallel
#'
"sigAssignmentEvaluation" <- function( x, beta, normalize_counts = TRUE, lambda_values = c(0.01, 0.05, 0.10, 0.20), max_iterations_lasso = 10000, num_processes = Inf, seed = NULL, verbose = TRUE, log_file = "" ) {
    
    # set the seed
    set.seed(seed)

    # check the input parameters
    x <- as.matrix(x)
    if(normalize_counts) {
        x_not_normalized <- x
        x <- (x/rowSums(x))*2500
    }
    lambda_values <- sort(unique(lambda_values))

    # initialize beta
    beta <- beta / rowSums(beta)
    rownames(beta) <- c("Background",paste0("S",1:(nrow(beta)-1)))
    colnames(beta) <- colnames(x)
    
    # setting up parallel execution
    parallel <- NULL
    close_parallel <- FALSE
    if(is.null(parallel)&&length(lambda_values)>1) {
        if(is.na(num_processes) || is.null(num_processes)) {
            parallel <- NULL
        }
        else if(num_processes==Inf) {
            cores <- as.integer((detectCores()-1))
            if(cores < 2) {
                parallel <- NULL
            }
            else {
                num_processes <- min(num_processes,length(lambda_values))
                parallel <- makeCluster(num_processes,outfile=log_file)
                close_parallel <- TRUE
            }
        }
        else {
            num_processes <- min(num_processes,length(lambda_values))
            parallel <- makeCluster(num_processes,outfile=log_file)
            close_parallel <- TRUE
        }
        
        if(verbose && !is.null(parallel)) {
            cat("Executing",num_processes,"processes via parallel...","\n")
        }
    }
    
    # structure to save the estimated signatures
    lambda_results <- array(list(),c(1,length(lambda_values)))
    rownames(lambda_results) <- paste0(as.character(dim(beta)[1]),"_signatures")
    colnames(lambda_results) <- paste0(as.character(lambda_values),"_lambda")
    
    if(verbose) {
        cat("Performing estimation of lambda range for alpha...","\n")
    }
    
    # perform signature discovery for all the values of lambda
    if(is.null(parallel)) {
        set.seed(round(runif(1)*100000))
        cont <- 0
        for(l in lambda_values) {

            # perform the inference
            curr_results <- sigAssignmentLasso(x = x, 
                                                beta = beta, 
                                                normalize_counts = normalize_counts, 
                                                lambda_rate_alpha = l, 
                                                max_iterations_lasso = max_iterations_lasso, 
                                                seed = round(runif(1)*100000), 
                                                verbose = FALSE)
                                     
            # save the results for the current configuration
            cont <- cont + 1

            # rescale alpha to the original magnitude
            if(normalize_counts) {
                curr_results$alpha <- curr_results$alpha * (rowSums(x_not_normalized)/2500)
            }

            lambda_results[[1,cont]] <- curr_results
                
            if(verbose) {
                cat("Progress",paste0(round((cont/(1*length(lambda_values)))*100,digits=3),"%..."),"\n")
            }
            
        }
    }
    else {
        res_clusterEvalQ <- clusterEvalQ(parallel,library("nnls"))
        res_clusterEvalQ <- clusterEvalQ(parallel,library("nnlasso"))
        clusterExport(parallel,varlist=c("x","beta","normalize_counts","max_iterations_lasso"),envir=environment())
        clusterExport(parallel,c('sigAssignmentLasso','nmfLassoDeconstruction'),envir=environment())
        clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
        lambda_results_res <- parLapply(parallel,lambda_values,function(l) {

            # perform the inference
            curr_results <- sigAssignmentLasso(x = x, 
                                                beta = beta, 
                                                normalize_counts = normalize_counts, 
                                                lambda_rate_alpha = l, 
                                                max_iterations_lasso = max_iterations_lasso, 
                                                seed = round(runif(1)*100000), 
                                                verbose = FALSE)

            return(curr_results)

        })
        for(i in 1:ncol(lambda_results)) {
            # rescale alpha to the original magnitude
            if(normalize_counts) {
                lambda_results_res[[i]]$alpha <- lambda_results_res[[i]]$alpha * (rowSums(x_not_normalized)/2500)
            }
            lambda_results[[1,i]] <- lambda_results_res[[i]]
        }
    }
    
    # close parallel
    if(close_parallel) {
        stopCluster(parallel)
    }
    
    return(lambda_results)
    
}

#' Perform the assignment of somatic mutational signatures to patients given a set of observed counts x and signatures beta.
#'
#' @examples
#' data(patients)
#' data(starting_betas_example)
#' beta = starting_betas_example[["5_signatures","Value"]]
#' res = sigAssignmentLasso(x=patients[1:100,],beta=beta,lambda_rate_alpha=0.05,seed=12345)
#'
#' @title sigAssignmentLasso
#' @param x count matrix for a set of n patients and 96 trinucleotides.
#' @param beta beta to be fixed during the estimation of alpha.
#' @param normalize_counts if true, the input count matrix x is normalize such that the patients have the same number of mutation.
#' @param lambda_rate_alpha value of LASSO to be used for alpha between 0 and 1. This value should be greater than 0. 1 is the value of LASSO 
#' that would shrink all the exposure values to 0 within one step. The higher lambda_rate_alpha is, the sparser are the resulting exposure values, 
#' but too large values may result in a reduced fit of the observed counts.
#' @param max_iterations_lasso Number of maximum iterations to be performed during the sparsification via Lasso.
#' @param seed Seed for reproducibility.
#' @param verbose boolean; Shall I print all messages?
#' @return A list with the discovered signatures and their assignment to patients. It includes 2 elements: 
#'              alpha: matrix of the assigned exposure values
#'              beta: matrix of the discovered signatures
#' @export sigAssignmentLasso
#' @import nnls
#' @import nnlasso
#'
"sigAssignmentLasso" <- function( x, beta, normalize_counts = TRUE, lambda_rate_alpha = 0.05, max_iterations_lasso = 10000, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)

    # check the input parameters
    if(lambda_rate_alpha<0) {
        if(verbose) {
            cat("The minimum value of lambda_rate_alpha is 0...","\n")
        }
        lambda_rate_alpha <- 0
    }
    if(lambda_rate_alpha>1) {
        if(verbose) {
            cat("The maximum value of lambda_rate_alpha is 1...","\n")
        }
        lambda_rate_alpha <- 1
    }
    x <- as.matrix(x)
    if(normalize_counts) {
        x_not_normalized <- x
        x <- (x/rowSums(x))*2500
    }

    # initialize beta
    beta <- beta / rowSums(beta)
    rownames(beta) <- c("Background",paste0("S",1:(nrow(beta)-1)))
    colnames(beta) <- colnames(x)
    
    if(verbose) {
        cat("Performing the assignment of signatures to patients with Lasso...","\n")
    }
    
    # perform the assignment of signatures to samples
    results = tryCatch({
        results = nmfLassoDeconstruction(x,beta,lambda_rate_alpha,max_iterations_lasso,round(runif(1)*100000))
        # rescale alpha to the original magnitude
        if(normalize_counts) {
            results$alpha <- results$alpha * (rowSums(x_not_normalized)/2500)
        }
        return(results)
    }, error = function(e) {
        warning("Lasso did not converge, you should try a lower value of lambda! Current settings: lambda_rate_alpha = ",lambda_rate_alpha,"...")
        return(list(alpha=NA,beta=NA))
    })
    
    return(results)
    
}

# perform assignment of somatic mutational signatures to samples using Lasso to ensure sparsity
"nmfLassoDeconstruction" <- function( x, beta, lambda_rate_alpha = 0.05, max_iterations_lasso = 10000, seed = NULL ) {
    
    # set the seed
    set.seed(seed)
    
    # n is the number of observations in x, i.e., the patients
    n <- dim(x)[1]
    
    # K is the number of signatures to be assigned
    K <- dim(beta)[1]
    
    # initialize beta
    beta <- beta / rowSums(beta)
    rownames(beta) <- c("Background",paste0("S",1:(nrow(beta)-1)))
    colnames(beta) <- colnames(x)
    
    # initialize alpha
    alpha <- matrix(0,nrow=n,ncol=K)
    rownames(alpha) <- 1:nrow(alpha)
    colnames(alpha) <- rownames(beta)
    
    # configuration 1 (turn off sparsity): alpha is estimated by Non-Negative Linear Least Squares 
    if(lambda_rate_alpha==0) {

        # estimate alpha by Non-Negative Linear Least Squares
        for(j in 1:n) {
            alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
        }

    }

    # configuration 2 (sparsity on alpha): alpha is estimated by Non-Negative Lasso 
    if(lambda_rate_alpha>0) {

        # estimate alpha by Non-Negative Lasso
        for(j in 1:n) {

            # compute independently for each patient the counts to be approximted
            target <- x[j,]
            
            # estimate alpha by Non-Negative Lasso to ensure sparsity
            max_lambda_value <- max(abs(beta %*% target))
            lambda_values_alpha <- max_lambda_value * lambda_rate_alpha
            alpha[j,] <- as.vector(nnlasso(x = t(beta), 
                                            y = target, 
                                            family = "normal", 
                                            lambda = lambda_values_alpha, 
                                            intercept = FALSE, 
                                            normalize = FALSE, 
                                            tol = 1e-05, 
                                            maxiter = max_iterations_lasso, 
                                            path = FALSE)$coef[2,])
            
        }

    }
    
    # save the results
    results <- list(alpha=alpha,beta=beta)
    
    return(results)
    
}
