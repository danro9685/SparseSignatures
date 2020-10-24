#' Perform the assessment of different nmfLasso solutions by cross validation for K (unknown) somatic mutational signatures given a set of observations x. The estimation can slow down because of 
#' memory usage and intensive computations, when a big number of cross validation repetitions is asked and when the grid search is performed for a lot of configurations. In this case, 
#' an advice may be to split the computation into multiple smaller sets.
#'
#' @examples
#' data(background)
#' data(patients)
#' res = nmfLassoCV(x=patients[1:100,],
#'      K=3:5,
#'      background_signature=background,
#'      nmf_runs=1,
#'      lambda_values_alpha=c(0.00),
#'      lambda_values_beta=c(0.00),
#'      cross_validation_repetitions=2,
#'      num_processes=NA,
#'      seed=12345)
#'
#' @title nmfLassoCV
#' @param x count matrix for a set of n patients and 96 trinucleotides.
#' @param K a range of numeric value (each of them greater than 1) indicating the number of signatures to be discovered.
#' @param starting_beta a list of starting beta value for each configuration of K. If it is NULL, starting betas are estimated by 
#' NMF.
#' @param background_signature background signature to be used. If not provided, a warning is thrown and an initial value for it is 
#' estimated by NMF. If beta is not NULL, this parameter is ignored.
#' @param normalize_counts if true, the input count matrix x is normalize such that the patients have the same number of mutation.
#' @param nmf_runs number of iteration (minimum 1) of NMF to be performed for a robust estimation of starting beta. If beta is not NULL, 
#' this parameter is ignored.
#' @param lambda_values_alpha value of LASSO to be used for alpha between 0 and 1. This value should be greater than 0. 1 is the value of LASSO 
#' that would shrink all the signatures to 0 within one step. The higher lambda_rate_alpha is, the sparser are the resulting exposures, 
#' but too large values may result in a reduced fit of the observed counts.
#' @param lambda_values_beta value of LASSO to be used for beta between 0 and 1. This value should be greater than 0. 1 is the value of LASSO 
#' that would shrink all the signatures to 0 within one step. The higher lambda_rate_beta is, the sparser are the resulting exposures, 
#' but too large values may result in a reduced fit of the observed counts.
#' @param cross_validation_entries Percentage of cells in the count matrix to be replaced by 0s during cross validation.
#' @param cross_validation_iterations For each configuration, the first time the signatures are discovered form a matrix with a 
#' percentage of values replaced by 0s. This may result in poor fit/results. Then, we perform predictions of these entries and replace them with 
#' such predicted values. This parameter is the number of restarts to be performed to improve this estimate and obtain more stable solutions.
#' @param cross_validation_repetitions Number of time cross-validation should be repeated. Higher values result in better estimate, but are computationally more expensive.
#' @param iterations Number of iterations to be performed. Each iteration corresponds to a first step where beta is fitted 
#' and a second step where alpha is fitted.
#' @param max_iterations_lasso Number of maximum iterations to be performed during the sparsification via Lasso.
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode, 
#' this parameter needs to be set to either NA or NULL.
#' @param seed Seed for reproducibility.
#' @param verbose boolean; Shall I print all messages?
#' @param log_file log file where to print outputs when using parallel. If parallel execution is disabled, this parameter is ignored.
#' @return A list of 2 elements: grid_search_mse and and grid_search_loglik. Here, grid_search_mse reports the mean squared error for each configuration of performed 
#' cross validation; grid_search_loglik reports for each configuration the number of times the algorithm converged.
#' @export nmfLassoCV
#' @import NMF
#' @import nnlasso
#' @import nnls
#' @import parallel
#'
"nmfLassoCV" <- function( x, K = 3:10, starting_beta = NULL, background_signature = NULL, normalize_counts = TRUE, nmf_runs = 10, lambda_values_alpha = c(0.00, 0.01, 0.05, 0.10), lambda_values_beta = c(0.00, 0.01, 0.05, 0.10), cross_validation_entries = 0.01, cross_validation_iterations = 5, cross_validation_repetitions = 50, iterations = 30, max_iterations_lasso = 10000, num_processes = Inf, seed = NULL, verbose = TRUE, log_file = "" ) {
    
    # set the seed
    set.seed(seed)

    # check the input parameters
    if(min(K)<2) {
        if(verbose) {
            cat("The minimum value of K is 2...","\n")
        }
        K <- K[K>=2]
        if(length(K)==0) {
            K <- 2
            warning("No value of K greater than 2 is provided: setting K to 2...")
        }
    }
    if(nmf_runs<1) {
        if(verbose) {
            cat("The minimum value of nmf_runs is 1...","\n")
        }
        nmf_runs <- 1
    }
    x <- as.matrix(x)
    if(normalize_counts) {
        x <- (x/rowSums(x))*2500
    }
    lambda_values_alpha <- sort(unique(lambda_values_alpha))
    lambda_values_beta <- sort(unique(lambda_values_beta))

    # perform a grid search to estimate the best values of K and lambda
    if(verbose) {
        cat("Performing a grid search to estimate the best values of K and lambda with a total of",cross_validation_repetitions,"cross validation repetitions...","\n")
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

    # estimate the background signature
    if(is.null(background_signature)&&is.null(starting_beta)) {
        warning("No background signature has been specified...")
        background_signature_beta <- basis(nmf(t(x),rank=1,nrun=nmf_runs))
        background_signature_beta <- t(background_signature_beta)
        background_signature_beta <- background_signature_beta / rowSums(background_signature_beta)
        background_signature <- as.numeric(background_signature_beta[1,])
    }
    
    # structure to save the starting values of beta for each K
    if(is.null(starting_beta)) {
        
        background_signature <- as.numeric(background_signature)
        background_signature <- background_signature / sum(background_signature)

        starting_beta <- array(list(),c(length(K),1))
        rownames(starting_beta) <- paste0(as.character(K),"_Signatures")
        colnames(starting_beta) <- "Value"

        # estimate the contribution of the background to the observed counts
        background_signature_beta <- t(background_signature)
        background_signature_alpha <- matrix(0,nrow=dim(x)[1],ncol=1)
        for(i in 1:dim(x)[1]) {
            background_signature_alpha[i,] <- nnls(t(background_signature_beta),as.vector(x[i,]))$x
        }
        curr_x <- x - (background_signature_alpha %*% background_signature_beta)
        curr_x[curr_x<0] <- 0
        if(length(which(rowSums(t(curr_x))==0))>0) {
            invalid_cols <- as.numeric(which(rowSums(t(curr_x))==0))
            for(inv_cols in invalid_cols) {
                curr_x[sample(1:length(curr_x[,inv_cols]),size=1),inv_cols] <- 1e-05
            }
        }

        # compute starting values for beta
        pos_k <- 0
        for(k in K) {
            # compute the initial values of beta
            if(verbose) {
                cat(paste0("Computing the initial values of beta by standard NMF for K equals to ",k,"..."),"\n")
            }
            pos_k <- pos_k + 1
            beta <- NULL
            beta <- basis(nmf(t(curr_x),rank=(k-1),nrun=nmf_runs))
            beta <- t(beta)
            beta <- rbind(background_signature,beta)
            
            # normalize and save beta for the current value of K
            rownames(beta) <- c("Background",paste0("S",1:(nrow(beta)-1)))
            colnames(beta) <- colnames(x)
            beta <- beta / rowSums(beta)
            starting_beta[[pos_k,1]] <- beta
        }

    }

    if(verbose) {
        cat("Starting cross validation with a total of",cross_validation_repetitions,"repetitions...","\n")
    }

    # structure to save the results of the grid search
    grid_search_mse <- array(list(),c(length(lambda_values_alpha),length(lambda_values_beta),length(K)),dimnames=list(paste0(as.character(lambda_values_alpha),"_Lambda_Alpha"),paste0(as.character(lambda_values_beta),"_Lambda_Beta"),paste0(as.character(K),"_Signatures")))
    grid_search_loglik <- array(list(),c(length(lambda_values_alpha),length(lambda_values_beta),length(K)),dimnames=list(paste0(as.character(lambda_values_alpha),"_Lambda_Alpha"),paste0(as.character(lambda_values_beta),"_Lambda_Beta"),paste0(as.character(K),"_Signatures")))

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
            for(k in K) {
                # consider the first k signatures to be used for the current configuration
                pos_k <- pos_k + 1

                # consider all the values for lambda
                pos_l_alpha <- 0
                for(l_alpha in lambda_values_alpha) {
                    pos_l_alpha <- pos_l_alpha + 1

                    pos_l_beta <- 0
                    for(l_beta in lambda_values_beta) {
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
                            curr_results <- nmfLasso(x = x_cv, 
                                                     K = k, 
                                                     beta = starting_beta[[pos_k,1]], 
                                                     background_signature = background_signature, 
                                                     normalize_counts = normalize_counts, 
                                                     nmf_runs = 10, 
                                                     lambda_rate_alpha = l_alpha, 
                                                     lambda_rate_beta = l_beta, 
                                                     iterations = iterations, 
                                                     max_iterations_lasso = max_iterations_lasso, 
                                                     seed = round(runif(1)*100000), 
                                                     verbose = FALSE)

                            curr_iter_alpha <- curr_results[["alpha"]]
                            curr_iter_beta <- curr_results[["beta"]]
                            curr_iter_loglik <- curr_results[["loglik_progression"]]

                        }
                        
                        # save an estimate of mean squared error and stability of the solution
                        if(!is.na(curr_iter_alpha[1])&&!is.na(curr_iter_beta[1])) {
                            curr_predicted_counts <- round(curr_iter_alpha%*%curr_iter_beta)
                            curr_true_considered_counts <- as.vector(x[cv_entries])
                            curr_predicted_considered_counts <- as.vector(curr_predicted_counts[cv_entries])
                            error <- mean((curr_true_considered_counts-curr_predicted_considered_counts)^2)
                            grid_search_mse[pos_l_alpha,pos_l_beta,pos_k][[1]] <- c(grid_search_mse[pos_l_alpha,pos_l_beta,pos_k][[1]],error)
                        }
                        if((!is.na(curr_iter_loglik)[1])&&(curr_iter_loglik[1]<curr_iter_loglik[length(curr_iter_loglik)])) {
                            grid_search_loglik[pos_l_alpha,pos_l_beta,pos_k][[1]] <- c(grid_search_loglik[pos_l_alpha,pos_l_beta,pos_k][[1]],"increasing")
                        }
                        else if((!is.na(curr_iter_loglik)[1])&&(curr_iter_loglik[1]>=curr_iter_loglik[length(curr_iter_loglik)])) {
                            grid_search_loglik[pos_l_alpha,pos_l_beta,pos_k][[1]] <- c(grid_search_loglik[pos_l_alpha,pos_l_beta,pos_k][[1]],"decreasing")
                        }
                        
                        if(verbose) {
                            cat("Progress",paste0(round((cont/(length(K)*length(lambda_values_alpha)*length(lambda_values_beta)))*100,digits=3),"%..."),"\n")
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
        clusterExport(parallel,varlist=c("x","K","starting_beta","background_signature","normalize_counts","cross_validation_entries"),envir=environment())
        clusterExport(parallel,varlist=c("lambda_values_alpha","lambda_values_beta","cross_validation_iterations","iterations","max_iterations_lasso"),envir=environment())
        clusterExport(parallel,varlist=c("grid_search_mse","grid_search_loglik","cross_validation_repetitions","verbose"),envir=environment())
        clusterExport(parallel,c('nmfLasso','nmfLassoDecomposition','basis','nmf'),envir=environment())
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
            for(k in K) {
                # consider the first k signatures to be used for the current configuration
                pos_k <- pos_k + 1

                # consider all the values for lambda
                pos_l_alpha <- 0
                for(l_alpha in lambda_values_alpha) {
                    pos_l_alpha <- pos_l_alpha + 1

                    pos_l_beta <- 0
                    for(l_beta in lambda_values_beta) {
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
                            curr_results <- nmfLasso(x = x_cv, 
                                                     K = k, 
                                                     beta = starting_beta[[pos_k,1]], 
                                                     background_signature = background_signature, 
                                                     normalize_counts = normalize_counts, 
                                                     nmf_runs = 10, 
                                                     lambda_rate_alpha = l_alpha, 
                                                     lambda_rate_beta = l_beta, 
                                                     iterations = iterations, 
                                                     max_iterations_lasso = max_iterations_lasso, 
                                                     seed = round(runif(1)*100000), 
                                                     verbose = FALSE)

                            curr_iter_alpha <- curr_results[["alpha"]]
                            curr_iter_beta <- curr_results[["beta"]]
                            curr_iter_loglik <- curr_results[["loglik_progression"]]

                        }
                        
                        # save an estimate of mean squared error and stability of the solution
                        if(!is.na(curr_iter_alpha[1])&&!is.na(curr_iter_beta[1])) {
                            curr_predicted_counts <- round(curr_iter_alpha%*%curr_iter_beta)
                            curr_true_considered_counts <- as.vector(x[cv_entries])
                            curr_predicted_considered_counts <- as.vector(curr_predicted_counts[cv_entries])
                            error <- mean((curr_true_considered_counts-curr_predicted_considered_counts)^2)
                            grid_search_mse[pos_l_alpha,pos_l_beta,pos_k][[1]] <- c(grid_search_mse[pos_l_alpha,pos_l_beta,pos_k][[1]],error)
                        }
                        if((!is.na(curr_iter_loglik)[1])&&(curr_iter_loglik[1]<curr_iter_loglik[length(curr_iter_loglik)])) {
                            grid_search_loglik[pos_l_alpha,pos_l_beta,pos_k][[1]] <- c(grid_search_loglik[pos_l_alpha,pos_l_beta,pos_k][[1]],"increasing")
                        }
                        else if((!is.na(curr_iter_loglik)[1])&&(curr_iter_loglik[1]>=curr_iter_loglik[length(curr_iter_loglik)])) {
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

#' Perform the evaluation of different nmfLasso solutions by bootstrap for K (unknown) somatic mutational signatures given a set of observations x. The estimation can slow down because of 
#' memory usage and intensive computations, when a big number of bootstrap repetitions is asked and when the analysis is performed for a big range of signatures (K). In this case, 
#' an advice may be to split the computation into multiple smaller sets.
#'
#' @examples
#' data(background)
#' data(patients)
#' res = nmfLassoBootstrap(x=patients[1:100,],
#'      K=3:5,
#'      background_signature=background,
#'      nmf_runs=1,
#'      bootstrap_repetitions=2,
#'      num_processes=NA,
#'      seed=12345)
#'
#' @title nmfLassoBootstrap
#' @param x count matrix for a set of n patients and 96 trinucleotides.
#' @param K a range of numeric value (each of them greater than 1) indicating the number of signatures to be discovered.
#' @param starting_beta a list of starting beta value for each configuration of K. If it is NULL, starting betas are estimated by 
#' NMF.
#' @param background_signature background signature to be used. If not provided, a warning is thrown and an initial value for it is 
#' estimated by NMF. If beta is not NULL, this parameter is ignored.
#' @param normalize_counts if true, the input count matrix x is normalize such that the patients have the same number of mutation.
#' @param nmf_runs number of iteration (minimum 1) of NMF to be performed for a robust estimation of starting beta. If beta is not NULL, 
#' this parameter is ignored.
#' @param bootstrap_repetitions Number of time bootstrap should be repeated. Higher values result in better estimate, but are computationally more expensive.
#' @param iterations Number of iterations to be performed. Each iteration corresponds to a first step where beta is fitted 
#' and a second step where alpha is fitted.
#' @param max_iterations_lasso Number of maximum iterations to be performed during the sparsification via Lasso.
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode, 
#' this parameter needs to be set to either NA or NULL.
#' @param seed Seed for reproducibility.
#' @param verbose boolean; Shall I print all messages?
#' @param log_file log file where to print outputs when using parallel. If parallel execution is disabled, this parameter is ignored.
#' @return A list of 3 elements: stability, RSS and evar. Here, stability reports the estimared cosine similarity for alpha and beta at each bootstrap 
#' repetition; RSS reports for each configuration the estimated residual sum of squares; finally, evar reports the explained variance.
#' @export nmfLassoBootstrap
#' @import NMF
#' @import nnlasso
#' @import nnls
#' @import parallel
#'
"nmfLassoBootstrap" <- function( x, K = 3:10, starting_beta = NULL, background_signature = NULL, normalize_counts = TRUE, nmf_runs = 10, bootstrap_repetitions = 50, iterations = 30, max_iterations_lasso = 10000, num_processes = Inf, seed = NULL, verbose = TRUE, log_file = "" ) {
    
    # set the seed
    set.seed(seed)

    # check the input parameters
    if(min(K)<2) {
        if(verbose) {
            cat("The minimum value of K is 2...","\n")
        }
        K <- K[K>=2]
        if(length(K)==0) {
            K <- 2
            warning("No value of K greater than 2 is provided: setting K to 2...")
        }
    }
    if(nmf_runs<1) {
        if(verbose) {
            cat("The minimum value of nmf_runs is 1...","\n")
        }
        nmf_runs <- 1
    }
    x <- as.matrix(x)
    if(normalize_counts) {
        x <- (x/rowSums(x))*2500
    }

    # no Lasso is considered here to speed up computations, i.e., this is an heuristic/approximated result
    lambda_values_alpha <- 0.00
    lambda_values_beta <- 0.00

    # perform a grid search to estimate the best values of K and lambda
    if(verbose) {
        cat("Performing a total of",bootstrap_repetitions,"bootstrap repetitions to assess nmfLasso solutions for different K ranks...","\n")
    }

    # setting up parallel execution
    parallel <- NULL
    close_parallel <- FALSE
    if(is.null(parallel)&&bootstrap_repetitions>1) {
        if(is.na(num_processes) || is.null(num_processes)) {
            parallel <- NULL
        }
        else if(num_processes==Inf) {
            cores <- as.integer((detectCores()-1))
            if(cores < 2) {
                parallel <- NULL
            }
            else {
                num_processes <- min(num_processes,bootstrap_repetitions)
                parallel <- makeCluster(num_processes,outfile=log_file)
                close_parallel <- TRUE
            }
        }
        else {
            num_processes <- min(num_processes,bootstrap_repetitions)
            parallel <- makeCluster(num_processes,outfile=log_file)
            close_parallel <- TRUE
        }
        
        if(verbose && !is.null(parallel)) {
            cat("Executing",num_processes,"processes via parallel...","\n")
        }
    }

    # estimate the background signature
    if(is.null(background_signature)&&is.null(starting_beta)) {
        warning("No background signature has been specified...")
        background_signature_beta <- basis(nmf(t(x),rank=1,nrun=nmf_runs))
        background_signature_beta <- t(background_signature_beta)
        background_signature_beta <- background_signature_beta / rowSums(background_signature_beta)
        background_signature <- as.numeric(background_signature_beta[1,])
    }
    
    # structure to save the starting values of beta for each K
    if(is.null(starting_beta)) {
        
        background_signature <- as.numeric(background_signature)
        background_signature <- background_signature / sum(background_signature)

        starting_beta <- array(list(),c(length(K),1))
        rownames(starting_beta) <- paste0(as.character(K),"_Signatures")
        colnames(starting_beta) <- "Value"

        # estimate the contribution of the background to the observed counts
        background_signature_beta <- t(background_signature)
        background_signature_alpha <- matrix(0,nrow=dim(x)[1],ncol=1)
        for(i in 1:dim(x)[1]) {
            background_signature_alpha[i,] <- nnls(t(background_signature_beta),as.vector(x[i,]))$x
        }
        curr_x <- x - (background_signature_alpha %*% background_signature_beta)
        curr_x[curr_x<0] <- 0
        if(length(which(rowSums(t(curr_x))==0))>0) {
            invalid_cols <- as.numeric(which(rowSums(t(curr_x))==0))
            for(inv_cols in invalid_cols) {
                curr_x[sample(1:length(curr_x[,inv_cols]),size=1),inv_cols] <- 1e-05
            }
        }

        # compute starting values for beta
        pos_k <- 0
        for(k in K) {
            # compute the initial values of beta
            if(verbose) {
                cat(paste0("Computing the initial values of beta by standard NMF for K equals to ",k,"..."),"\n")
            }
            pos_k <- pos_k + 1
            beta <- NULL
            beta <- basis(nmf(t(curr_x),rank=(k-1),nrun=nmf_runs))
            beta <- t(beta)
            beta <- rbind(background_signature,beta)
            
            # normalize and save beta for the current value of K
            rownames(beta) <- c("Background",paste0("S",1:(nrow(beta)-1)))
            colnames(beta) <- colnames(x)
            beta <- beta / rowSums(beta)
            starting_beta[[pos_k,1]] <- beta
        }

    }

    if(verbose) {
        cat("Starting bootstrap for a total of",bootstrap_repetitions,"repetitions...","\n")
    }

    # structure to save the results of the grid search
    bootstrap_rss <- array(list(),c(length(K),2))
    rownames(bootstrap_rss) <- paste0(as.character(K),"_Signatures")
    colnames(bootstrap_rss) <- c("RSS bootstrapped data","RSS randomized data")
    bootstrap_evar <- bootstrap_rss
    colnames(bootstrap_evar) <- c("Explained Variance bootstrapped data","Explained Variance randomized data")
    bootstrap_stability <- array(list(),c(length(K),4))
    rownames(bootstrap_stability) <- paste0(as.character(K),"_Signatures")
    colnames(bootstrap_stability) <- c("Similarity alpha bootstrapped data","Similarity beta bootstrapped data","Similarity alpha randomized data","Similarity beta randomized data")
    starting_values_bootstrap = bootstrap_stability

    # generate a set of randomized data from x
    set.seed(round(runif(1)*100000))
    randomized_x <- x
    for(i in 1:ncol(randomized_x)) {
        randomized_x[,i] <- randomized_x[sample(1:dim(randomized_x)[1],size=dim(randomized_x)[1],replace=FALSE),i]
    }
    for(j in 1:nrow(randomized_x)) {
        randomized_x[j,] <- randomized_x[j,sample(1:dim(randomized_x)[2],size=dim(randomized_x)[2],replace=FALSE)]
    }

    # compute starting values of alpha and beta in order to estimate stability of the considered solutions (and ranks)
    pos_k <- 0
    for(k in K) {

        pos_k <- pos_k + 1

        is.valid <- FALSE
        while(!is.valid) {
            curr_results <- nmfLasso(x = x, 
                                     K = k, 
                                     beta = starting_beta[[pos_k,1]], 
                                     background_signature = background_signature, 
                                     normalize_counts = normalize_counts, 
                                     nmf_runs = 10, 
                                     lambda_rate_alpha = lambda_values_alpha, 
                                     lambda_rate_beta = lambda_values_beta, 
                                     iterations = iterations, 
                                     max_iterations_lasso = max_iterations_lasso, 
                                     seed = round(runif(1)*100000), 
                                     verbose = FALSE)
            if(!is.na(curr_results$best_loglik)) {
                is.valid <- TRUE
            }
        }

        starting_values_bootstrap[[pos_k,1]] <- curr_results$alpha
        starting_values_bootstrap[[pos_k,2]] <- curr_results$beta

        is.valid <- FALSE
        while(!is.valid) {
            curr_results <- nmfLasso(x = randomized_x, 
                                     K = k, 
                                     beta = starting_beta[[pos_k,1]], 
                                     background_signature = background_signature, 
                                     normalize_counts = normalize_counts, 
                                     nmf_runs = 10, 
                                     lambda_rate_alpha = lambda_values_alpha, 
                                     lambda_rate_beta = lambda_values_beta, 
                                     iterations = iterations, 
                                     max_iterations_lasso = max_iterations_lasso, 
                                     seed = round(runif(1)*100000), 
                                     verbose = FALSE)
            if(!is.na(curr_results$best_loglik)) {
                is.valid <- TRUE
            }
        }

        starting_values_bootstrap[[pos_k,3]] <- curr_results$alpha
        starting_values_bootstrap[[pos_k,4]] <- curr_results$beta

    }

    # now starting cross validation
    if(is.null(parallel)) {

        # perform a total of bootstrap_repetitions repetitions of bootstrap
        for(nboot in 1:bootstrap_repetitions) {

            if(verbose) {
                cat(paste0("Performing repetition ",nboot," out of ",bootstrap_repetitions,"..."),"\n")
            }

            # generate bootstrapped dataset from data
            bootstrapped_data <- x
            for(i in 1:nrow(x)) {
                num_mutations <- sum(x[i,])
                mut_distribution <- as.numeric(x[i,]/sum(x[i,]))
                sampling <- sample(1:96,size=num_mutations,replace=TRUE,prob=mut_distribution)
                bootstrapped_data[i,] <- 0
                for(j in sampling) {
                    bootstrapped_data[i,j] <- bootstrapped_data[i,j] + 1
                }
            }
            
            # generate randomized data from the bootstrapped data
            randomized_data <- bootstrapped_data
            for(i in 1:ncol(randomized_data)) {
                randomized_data[,i] <- randomized_data[sample(1:dim(randomized_data)[1],size=dim(randomized_data)[1],replace=FALSE),i]
            }
            for(j in 1:nrow(randomized_data)) {
                randomized_data[j,] <- randomized_data[j,sample(1:dim(randomized_data)[2],size=dim(randomized_data)[2],replace=FALSE)]
            }
    
            # perform inference on bootstrapped and randomized data
            pos_k <- 0
            for(k in K) {
                
                pos_k <- pos_k + 1

                # perform the analysis for bootstrapped data
                res <- nmfLasso(x=bootstrapped_data,K=k,beta=starting_beta[[pos_k,1]],lambda_rate_alpha=lambda_values_alpha,lambda_rate_beta=lambda_values_beta,verbose=FALSE)

                if(!is.na(res$best_loglik)) {

                    # compute cosine similarity of alpha for bootstrapped data
                    similarity <- NULL
                    for(i in 1:ncol(res$alpha)) {
                        curr_similarity <- ((sum(starting_values_bootstrap[[pos_k,1]][,i]*res$alpha[,i]))/sqrt(sum(starting_values_bootstrap[[pos_k,1]][,i]^2)*sum(res$alpha[,i]^2)))
                        similarity <- c(similarity,curr_similarity)
                    }
                    bootstrap_stability[[pos_k,"Similarity alpha bootstrapped data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity alpha bootstrapped data"]]),mean(similarity))

                    # compute cosine similarity of beta for bootstrapped data
                    similarity <- NULL
                    for(i in 1:nrow(res$beta)) {
                        curr_similarity <- ((sum(starting_values_bootstrap[[pos_k,2]][i,]*res$beta[i,]))/sqrt(sum(starting_values_bootstrap[[pos_k,2]][i,]^2)*sum(res$beta[i,]^2)))
                        similarity <- c(similarity,curr_similarity)
                    }
                    bootstrap_stability[[pos_k,"Similarity beta bootstrapped data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity beta bootstrapped data"]]),mean(similarity))

                    # compute RSS for bootstrapped data
                    predictions <- res$alpha%*%res$beta
                    mse <- NULL
                    for(i in 1:nrow(predictions)) {
                        for(j in 1:ncol(predictions)) {
                            mse <- c(mse,((bootstrapped_data[i,j]-predictions[i,j])^2))
                        }
                    }
                    mse <- sum(mse)
                    bootstrap_rss[[pos_k,"RSS bootstrapped data"]] <- c(unlist(bootstrap_rss[[pos_k,"RSS bootstrapped data"]]),mse)

                    # compute Explained Variance for bootstrapped data
                    evar <- (1 - (mse/sum(bootstrapped_data^2)))
                    bootstrap_evar[[pos_k,"Explained Variance bootstrapped data"]] <- c(unlist(bootstrap_evar[[pos_k,"Explained Variance bootstrapped data"]]),evar)

                }
                else {
                    bootstrap_stability[[pos_k,"Similarity alpha bootstrapped data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity alpha bootstrapped data"]]),NA)
                    bootstrap_stability[[pos_k,"Similarity beta bootstrapped data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity beta bootstrapped data"]]),NA)
                    bootstrap_rss[[pos_k,"RSS bootstrapped data"]] <- c(unlist(bootstrap_rss[[pos_k,"RSS bootstrapped data"]]),NA)
                    bootstrap_evar[[pos_k,"Explained Variance bootstrapped data"]] <- c(unlist(bootstrap_evar[[pos_k,"Explained Variance bootstrapped data"]]),NA)
                }

                # perform the analysis for randomized data
                res <- nmfLasso(x=randomized_data,K=k,beta=starting_beta[[pos_k,1]],lambda_rate_alpha=lambda_values_alpha,lambda_rate_beta=lambda_values_beta,verbose=FALSE)

                if(!is.na(res$best_loglik)) {

                    # compute cosine similarity of alpha for randomized data
                    similarity <- NULL
                    for(i in 1:ncol(res$alpha)) {
                        curr_similarity <- ((sum(starting_values_bootstrap[[pos_k,3]][,i]*res$alpha[,i]))/sqrt(sum(starting_values_bootstrap[[pos_k,3]][,i]^2)*sum(res$alpha[,i]^2)))
                        similarity <- c(similarity,curr_similarity)
                    }
                    bootstrap_stability[[pos_k,"Similarity alpha randomized data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity alpha randomized data"]]),mean(similarity))

                    # compute cosine similarity of beta for randomized data
                    similarity <- NULL
                    for(i in 1:nrow(res$beta)) {
                        curr_similarity <- ((sum(starting_values_bootstrap[[pos_k,4]][i,]*res$beta[i,]))/sqrt(sum(starting_values_bootstrap[[pos_k,4]][i,]^2)*sum(res$beta[i,]^2)))
                        similarity <- c(similarity,curr_similarity)
                    }
                    bootstrap_stability[[pos_k,"Similarity beta randomized data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity beta randomized data"]]),mean(similarity))

                    # compute RSS for randomized data
                    predictions <- res$alpha%*%res$beta
                    mse <- NULL
                    for(i in 1:nrow(predictions)) {
                        for(j in 1:ncol(predictions)) {
                            mse <- c(mse,((randomized_data[i,j]-predictions[i,j])^2))
                        }
                    }
                    mse <- sum(mse)
                    bootstrap_rss[[pos_k,"RSS randomized data"]] <- c(unlist(bootstrap_rss[[pos_k,"RSS randomized data"]]),mse)

                    # compute Explained Variance for randomized data
                    evar <- (1 - (mse/sum(randomized_data^2)))
                    bootstrap_evar[[pos_k,"Explained Variance randomized data"]] <- c(unlist(bootstrap_evar[[pos_k,"Explained Variance randomized data"]]),evar)

                }
                else {
                    bootstrap_stability[[pos_k,"Similarity alpha randomized data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity alpha randomized data"]]),NA)
                    bootstrap_stability[[pos_k,"Similarity beta randomized data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity beta randomized data"]]),NA)
                    bootstrap_rss[[pos_k,"RSS randomized data"]] <- c(unlist(bootstrap_rss[[pos_k,"RSS randomized data"]]),NA)
                    bootstrap_evar[[pos_k,"Explained Variance randomized data"]] <- c(unlist(bootstrap_evar[[pos_k,"Explained Variance randomized data"]]),NA)
                }
                
            }

            if(verbose) {
                cat("Progress",paste0(round((nboot/bootstrap_repetitions)*100,digits=3),"%..."),"\n")
            }

        }

        # save the results
        results <- list(stability=bootstrap_stability,RSS=bootstrap_rss,evar=bootstrap_evar)

    }
    else {

        # perform the inference
        res_clusterEvalQ <- clusterEvalQ(parallel,library("nnls"))
        res_clusterEvalQ <- clusterEvalQ(parallel,library("nnlasso"))
        clusterExport(parallel,varlist=c("x","K","starting_beta","background_signature","normalize_counts"),envir=environment())
        clusterExport(parallel,varlist=c("lambda_values_alpha","lambda_values_beta","iterations","max_iterations_lasso"),envir=environment())
        clusterExport(parallel,varlist=c("bootstrap_stability","bootstrap_rss","bootstrap_evar","bootstrap_repetitions","verbose"),envir=environment())
        clusterExport(parallel,c('nmfLasso','nmfLassoDecomposition','basis','nmf'),envir=environment())
        clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
        curr_results <- parLapply(parallel,1:bootstrap_repetitions,function(nboot) {

            if(verbose) {
                cat(paste0("Performing repetition ",nboot," out of ",bootstrap_repetitions,"..."),"\n")
            }

            # generate bootstrapped dataset from data
            bootstrapped_data <- x
            for(i in 1:nrow(x)) {
                num_mutations <- sum(x[i,])
                mut_distribution <- as.numeric(x[i,]/sum(x[i,]))
                sampling <- sample(1:96,size=num_mutations,replace=TRUE,prob=mut_distribution)
                bootstrapped_data[i,] <- 0
                for(j in sampling) {
                    bootstrapped_data[i,j] <- bootstrapped_data[i,j] + 1
                }
            }
            
            # generate randomized data from the bootstrapped data
            randomized_data <- bootstrapped_data
            for(i in 1:ncol(randomized_data)) {
                randomized_data[,i] <- randomized_data[sample(1:dim(randomized_data)[1],size=dim(randomized_data)[1],replace=FALSE),i]
            }
            for(j in 1:nrow(randomized_data)) {
                randomized_data[j,] <- randomized_data[j,sample(1:dim(randomized_data)[2],size=dim(randomized_data)[2],replace=FALSE)]
            }
    
            # perform inference on bootstrapped and randomized data
            pos_k <- 0
            for(k in K) {
                
                pos_k <- pos_k + 1

                # perform the analysis for bootstrapped data
                res <- nmfLasso(x=bootstrapped_data,K=k,beta=starting_beta[[pos_k,1]],lambda_rate_alpha=lambda_values_alpha,lambda_rate_beta=lambda_values_beta,verbose=FALSE)

                if(!is.na(res$best_loglik)) {

                    # compute cosine similarity of alpha for bootstrapped data
                    similarity <- NULL
                    for(i in 1:ncol(res$alpha)) {
                        curr_similarity <- ((sum(starting_values_bootstrap[[pos_k,1]][,i]*res$alpha[,i]))/sqrt(sum(starting_values_bootstrap[[pos_k,1]][,i]^2)*sum(res$alpha[,i]^2)))
                        similarity <- c(similarity,curr_similarity)
                    }
                    bootstrap_stability[[pos_k,"Similarity alpha bootstrapped data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity alpha bootstrapped data"]]),mean(similarity))

                    # compute cosine similarity of beta for bootstrapped data
                    similarity <- NULL
                    for(i in 1:nrow(res$beta)) {
                        curr_similarity <- ((sum(starting_values_bootstrap[[pos_k,2]][i,]*res$beta[i,]))/sqrt(sum(starting_values_bootstrap[[pos_k,2]][i,]^2)*sum(res$beta[i,]^2)))
                        similarity <- c(similarity,curr_similarity)
                    }
                    bootstrap_stability[[pos_k,"Similarity beta bootstrapped data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity beta bootstrapped data"]]),mean(similarity))

                    # compute RSS for bootstrapped data
                    predictions <- res$alpha%*%res$beta
                    mse <- NULL
                    for(i in 1:nrow(predictions)) {
                        for(j in 1:ncol(predictions)) {
                            mse <- c(mse,((bootstrapped_data[i,j]-predictions[i,j])^2))
                        }
                    }
                    mse <- sum(mse)
                    bootstrap_rss[[pos_k,"RSS bootstrapped data"]] <- c(unlist(bootstrap_rss[[pos_k,"RSS bootstrapped data"]]),mse)

                    # compute Explained Variance for bootstrapped data
                    evar <- (1 - (mse/sum(bootstrapped_data^2)))
                    bootstrap_evar[[pos_k,"Explained Variance bootstrapped data"]] <- c(unlist(bootstrap_evar[[pos_k,"Explained Variance bootstrapped data"]]),evar)

                }
                else {
                    bootstrap_stability[[pos_k,"Similarity alpha bootstrapped data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity alpha bootstrapped data"]]),NA)
                    bootstrap_stability[[pos_k,"Similarity beta bootstrapped data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity beta bootstrapped data"]]),NA)
                    bootstrap_rss[[pos_k,"RSS bootstrapped data"]] <- c(unlist(bootstrap_rss[[pos_k,"RSS bootstrapped data"]]),NA)
                    bootstrap_evar[[pos_k,"Explained Variance bootstrapped data"]] <- c(unlist(bootstrap_evar[[pos_k,"Explained Variance bootstrapped data"]]),NA)
                }

                # perform the analysis for randomized data
                res <- nmfLasso(x=randomized_data,K=k,beta=starting_beta[[pos_k,1]],lambda_rate_alpha=lambda_values_alpha,lambda_rate_beta=lambda_values_beta,verbose=FALSE)

                if(!is.na(res$best_loglik)) {

                    # compute cosine similarity of alpha for randomized data
                    similarity <- NULL
                    for(i in 1:ncol(res$alpha)) {
                        curr_similarity <- ((sum(starting_values_bootstrap[[pos_k,3]][,i]*res$alpha[,i]))/sqrt(sum(starting_values_bootstrap[[pos_k,3]][,i]^2)*sum(res$alpha[,i]^2)))
                        similarity <- c(similarity,curr_similarity)
                    }
                    bootstrap_stability[[pos_k,"Similarity alpha randomized data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity alpha randomized data"]]),mean(similarity))

                    # compute cosine similarity of beta for randomized data
                    similarity <- NULL
                    for(i in 1:nrow(res$beta)) {
                        curr_similarity <- ((sum(starting_values_bootstrap[[pos_k,4]][i,]*res$beta[i,]))/sqrt(sum(starting_values_bootstrap[[pos_k,4]][i,]^2)*sum(res$beta[i,]^2)))
                        similarity <- c(similarity,curr_similarity)
                    }
                    bootstrap_stability[[pos_k,"Similarity beta randomized data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity beta randomized data"]]),mean(similarity))

                    # compute RSS for randomized data
                    predictions <- res$alpha%*%res$beta
                    mse <- NULL
                    for(i in 1:nrow(predictions)) {
                        for(j in 1:ncol(predictions)) {
                            mse <- c(mse,((randomized_data[i,j]-predictions[i,j])^2))
                        }
                    }
                    mse <- sum(mse)
                    bootstrap_rss[[pos_k,"RSS randomized data"]] <- c(unlist(bootstrap_rss[[pos_k,"RSS randomized data"]]),mse)

                    # compute Explained Variance for randomized data
                    evar <- (1 - (mse/sum(randomized_data^2)))
                    bootstrap_evar[[pos_k,"Explained Variance randomized data"]] <- c(unlist(bootstrap_evar[[pos_k,"Explained Variance randomized data"]]),evar)

                }
                else {
                    bootstrap_stability[[pos_k,"Similarity alpha randomized data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity alpha randomized data"]]),NA)
                    bootstrap_stability[[pos_k,"Similarity beta randomized data"]] <- c(unlist(bootstrap_stability[[pos_k,"Similarity beta randomized data"]]),NA)
                    bootstrap_rss[[pos_k,"RSS randomized data"]] <- c(unlist(bootstrap_rss[[pos_k,"RSS randomized data"]]),NA)
                    bootstrap_evar[[pos_k,"Explained Variance randomized data"]] <- c(unlist(bootstrap_evar[[pos_k,"Explained Variance randomized data"]]),NA)
                }
                
            }

            # save the results
            curr_results <- list(stability=bootstrap_stability,RSS=bootstrap_rss,evar=bootstrap_evar)

            return(curr_results)

        })

        # save the results
        results <- list(stability=bootstrap_stability,RSS=bootstrap_rss,evar=bootstrap_evar)
        for(par_res in 1:length(curr_results)) {
            curr_par_res <- curr_results[[par_res]]
            curr_par_res_stability <- curr_par_res$stability
            curr_par_res_rss <- curr_par_res$RSS
            curr_par_res_evar <- curr_par_res$evar
            for(a in 1:dim(curr_par_res_stability)[1]) {
                for(b in 1:dim(curr_par_res_stability)[2]) {
                    results$stability[[a,b]] <- c(unlist(results$stability[a,b]),unlist(curr_par_res_stability[a,b]))
                }
            }
            for(a in 1:dim(curr_par_res_rss)[1]) {
                for(b in 1:dim(curr_par_res_rss)[2]) {
                    results$RSS[[a,b]] <- c(unlist(results$RSS[a,b]),unlist(curr_par_res_rss[a,b]))
                    results$evar[[a,b]] <- c(unlist(results$evar[a,b]),unlist(curr_par_res_evar[a,b]))
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

#' Perform a robust estimation of the starting betas for the nmfLasso method
#'
#' @examples
#' data(background)
#' data(patients)
#' res = startingBetaEstimation(x=patients[1:100,],
#'      K=3:5,
#'      background_signature=background,
#'      nmf_runs=1,
#'      seed=12345)
#'
#' @title startingBetaEstimation
#' @param x count matrix for a set of n patients and 96 trinucleotides.
#' @param K numeric value (minimum 2) indicating the number of signatures to be discovered.
#' @param background_signature background signature to be used. If not provided, a warning is thrown and an initial value for it is 
#' estimated by NMF. If beta is not NULL, this parameter is ignored.
#' @param normalize_counts if true, the input count matrix x is normalize such that the patients have the same number of mutation.
#' @param nmf_runs number of iteration (minimum 1) of NMF to be performed for a robust estimation of starting beta. If beta is not NULL, 
#' this parameter is ignored.
#' @param seed Seed for reproducibility.
#' @param verbose boolean; Shall I print all messages?
#' @return A list containing the starting beta values for each configuration of K.
#' @export startingBetaEstimation
#' @import NMF
#' @import nnls
#'
"startingBetaEstimation" <- function( x, K = 3:10, background_signature = NULL, normalize_counts = TRUE, nmf_runs = 10, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)

    # check the input parameters
    if(min(K)<2) {
        if(verbose) {
            cat("The minimum value of K is 2...","\n")
        }
        K <- K[K>=2]
        if(length(K)==0) {
            K <- 2
            warning("No value of K greater than 2 is provided: setting K to 2...")
        }
    }
    if(nmf_runs<1) {
        if(verbose) {
            cat("The minimum value of nmf_runs is 1...","\n")
        }
        nmf_runs <- 1
    }
    x <- as.matrix(x)
    if(normalize_counts) {
        x <- (x/rowSums(x))*2500
    }

    # estimate the background signature
    if(is.null(background_signature)) {
        warning("No background signature has been specified...")
        background_signature_beta <- basis(nmf(t(x),rank=1,nrun=nmf_runs))
        background_signature_beta <- t(background_signature_beta)
        background_signature_beta <- background_signature_beta / rowSums(background_signature_beta)
        background_signature <- as.numeric(background_signature_beta[1,])
    }
    background_signature <- as.numeric(background_signature)
    background_signature <- background_signature / sum(background_signature)

    # estimate the contribution of the background to the observed counts
    background_signature_beta <- t(background_signature)
    background_signature_alpha <- matrix(0,nrow=dim(x)[1],ncol=1)
    for(i in 1:dim(x)[1]) {
        background_signature_alpha[i,] <- nnls(t(background_signature_beta),as.vector(x[i,]))$x
    }
    curr_x <- x - (background_signature_alpha %*% background_signature_beta)
    curr_x[curr_x<0] <- 0
    if(length(which(rowSums(t(curr_x))==0))>0) {
        invalid_cols <- as.numeric(which(rowSums(t(curr_x))==0))
        for(inv_cols in invalid_cols) {
            curr_x[sample(1:length(curr_x[,inv_cols]),size=1),inv_cols] <- 1e-05
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

    # consider all the values of K
    pos_k <- 0
    for(k in K) {
    
        pos_k <- pos_k + 1
        beta <- NULL

        # compute the initial values of beta
        if(verbose) {
            cat("Computing the initial values of beta by standard NMF...","\n")
        }
        beta <- basis(nmf(t(curr_x),rank=(k-1),nrun=nmf_runs))
        beta <- t(beta)
        beta <- rbind(background_signature,beta)
        
        # normalize and save beta for the current value of K
        beta <- beta / rowSums(beta)
        rownames(beta) <- c("Background",paste0("S",1:(nrow(beta)-1)))
        colnames(beta) <- colnames(x)
        starting_beta[[pos_k,1]] <- beta

        if(verbose) {
            cat("Progress",paste0(round((pos_k/length(K))*100,digits=3),"%..."),"\n")
        }

    }

    return(starting_beta)
    
}

#' Estimate the range of lambda values for alpha to be considered in the signature inference. Note that too small values of lambda 
#' result in dense exposures, but too large values lead to bad fit of the counts.
#'
#' @examples
#' data(background)
#' data(patients)
#' res = lambdaRangeAlphaEvaluation(x=patients[1:100,],
#'      K=5,background_signature=background,
#'      nmf_runs=1,lambda_values=c(0.01,0.05),
#'      num_processes=NA,
#'      seed=12345)
#'
#' @title lambdaRangeAlphaEvaluation
#' @param x count matrix for a set of n patients and 96 trinucleotides.
#' @param K numeric value (minimum 2) indicating the number of signatures to be discovered.
#' @param beta starting beta for the estimation. If it is NULL, starting beta is estimated by NMF.
#' @param background_signature background signature to be used. If not provided, a warning is thrown and an initial value for it is 
#' estimated by NMF. If beta is not NULL, this parameter is ignored.
#' @param normalize_counts if true, the input count matrix x is normalize such that the patients have the same number of mutation.
#' @param nmf_runs number of iteration (minimum 1) of NMF to be performed for a robust estimation of starting beta. If beta is not NULL, 
#' this parameter is ignored.
#' @param lambda_values value of LASSO to be used for alpha between 0 and 1. This value should be greater than 0. 1 is the value of LASSO 
#' that would shrink all the signatures to 0 within one step. The higher lambda_values is, the sparser are the resulting exposures, 
#' but too large values may result in a reduced fit of the observed counts.
#' @param iterations Number of iterations to be performed. Each iteration corresponds to a first step where beta is fitted 
#' and a second step where alpha is fitted.
#' @param max_iterations_lasso Number of maximum iterations to be performed during the sparsification via Lasso.
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode, 
#' this parameter needs to be set to either NA or NULL.
#' @param seed Seed for reproducibility.
#' @param verbose boolean; Shall I print all messages?
#' @param log_file log file where to print outputs when using parallel. If parallel execution is disabled, this parameter is ignored.
#' @return A list corresponding to results of the function nmfLasso for each value of lambda to be tested. This function allows 
#' to test a good range of lambda values for alpha to be considered. One should keep in mind that too small values generate dense solution, 
#' while too high ones leads to poor fit. This behavior is resampled in the values of loglik_progression, which should be increasing: 
#' too small values of lambda results in unstable log-likelihood through the iterations, while too large values make log-likelihood 
#' drop. 
#' @export lambdaRangeAlphaEvaluation
#' @import NMF
#' @import nnlasso
#' @import nnls
#' @import parallel
#'
"lambdaRangeAlphaEvaluation" <- function( x, K = 5, beta = NULL, background_signature = NULL, normalize_counts = TRUE, nmf_runs = 10, lambda_values = c(0.01, 0.05, 0.10, 0.20), iterations = 30, max_iterations_lasso = 10000, num_processes = Inf, seed = NULL, verbose = TRUE, log_file = "" ) {
    
    # set the seed
    set.seed(seed)

    # check the input parameters
    if(K<2) {
        if(verbose) {
            cat("The minimum value of K is 2...","\n")
        }
        K <- 2
    }
    if(nmf_runs<1) {
        if(verbose) {
            cat("The minimum value of nmf_runs is 1...","\n")
        }
        nmf_runs <- 1
    }
    x <- as.matrix(x)
    if(normalize_counts) {
        x_not_normalized <- x
        x <- (x/rowSums(x))*2500
    }
    lambda_values <- sort(unique(lambda_values))

    # compute the initial values of beta
    if(is.null(beta)) {

        if(verbose) {
            cat("Computing the initial values of beta by standard NMF...","\n")
        }
        
        # add a signature to beta (leading to K signatures in total) to explicitly model background noise
        if(is.null(background_signature)) {
            warning("No background signature has been specified...")
            background_signature_beta <- basis(nmf(t(x),rank=1,nrun=nmf_runs))
            background_signature_beta <- t(background_signature_beta)
            background_signature_beta <- background_signature_beta / rowSums(background_signature_beta)
            background_signature <- as.numeric(background_signature_beta[1,])
        }
        else {
            background_signature <- as.numeric(background_signature)
            background_signature <- background_signature / sum(background_signature)
            background_signature_beta <- t(background_signature)
        }

        background_signature_alpha <- matrix(0,nrow=dim(x)[1],ncol=1)
        for(i in 1:dim(x)[1]) {
            background_signature_alpha[i,] <- nnls(t(background_signature_beta),as.vector(x[i,]))$x
        }
        curr_x <- x - (background_signature_alpha %*% background_signature_beta)
        curr_x[curr_x<0] <- 0
        if(length(which(rowSums(t(curr_x))==0))>0) {
            invalid_cols <- as.numeric(which(rowSums(t(curr_x))==0))
            for(inv_cols in invalid_cols) {
                curr_x[sample(1:length(curr_x[,inv_cols]),size=1),inv_cols] <- 1e-05
            }
        }
        beta <- basis(nmf(t(curr_x),rank=(K-1),nrun=nmf_runs))
        beta <- t(beta)
        beta <- rbind(background_signature,beta)

    }
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
    lambda_results <- array(list(),c(length(K),length(lambda_values)))
    rownames(lambda_results) <- paste0(as.character(K),"_signatures")
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
            curr_results <- nmfLasso(x = x, 
                                     K = K, 
                                     beta = beta, 
                                     normalize_counts = normalize_counts, 
                                     lambda_rate_alpha = l, 
                                     lambda_rate_beta = 0.0, 
                                     iterations = iterations, 
                                     max_iterations_lasso = max_iterations_lasso, 
                                     seed = round(runif(1)*100000), 
                                     verbose = FALSE)
                                     
            # save the results for the current configuration
            cont <- cont + 1

            # rescale alpha to the original magnitude
            if(normalize_counts) {
                curr_results$starting_alpha <- curr_results$starting_alpha * (rowSums(x_not_normalized)/2500)
                curr_results$alpha <- curr_results$alpha * (rowSums(x_not_normalized)/2500)
            }

            lambda_results[[1,cont]] <- curr_results
                
            if(verbose) {
                cat("Progress",paste0(round((cont/(length(K)*length(lambda_values)))*100,digits=3),"%..."),"\n")
            }
            
        }
    }
    else {
        res_clusterEvalQ <- clusterEvalQ(parallel,library("nnls"))
        res_clusterEvalQ <- clusterEvalQ(parallel,library("nnlasso"))
        clusterExport(parallel,varlist=c("x","K","beta","normalize_counts","iterations","max_iterations_lasso"),envir=environment())
        clusterExport(parallel,c('nmfLasso','nmfLassoDecomposition','basis','nmf'),envir=environment())
        clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
        lambda_results_res <- parLapply(parallel,lambda_values,function(l) {

            # perform the inference
            curr_results <- nmfLasso(x = x, 
                                     K = K, 
                                     beta = beta, 
                                     normalize_counts = normalize_counts, 
                                     lambda_rate_alpha = l, 
                                     lambda_rate_beta = 0.0, 
                                     iterations = iterations, 
                                     max_iterations_lasso = max_iterations_lasso, 
                                     seed = round(runif(1)*100000), 
                                     verbose = FALSE)

            return(curr_results)

        })
        for(i in 1:ncol(lambda_results)) {
            # rescale alpha to the original magnitude
            if(normalize_counts) {
                lambda_results_res[[i]]$starting_alpha <- lambda_results_res[[i]]$starting_alpha * (rowSums(x_not_normalized)/2500)
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

#' Estimate the range of lambda values for beta to be considered in the signature inference. Note that too small values of lambda 
#' result in dense signatures, but too large values lead to bad fit of the counts.
#'
#' @examples
#' data(background)
#' data(patients)
#' res = lambdaRangeBetaEvaluation(x=patients[1:100,],
#'      K=5,
#'      background_signature=background,
#'      nmf_runs=1,
#'      lambda_values=c(0.01,0.05),
#'      num_processes=NA,
#'      seed=12345)
#'
#' @title lambdaRangeBetaEvaluation
#' @param x count matrix for a set of n patients and 96 trinucleotides.
#' @param K numeric value (minimum 2) indicating the number of signatures to be discovered.
#' @param beta starting beta for the estimation. If it is NULL, starting beta is estimated by NMF.
#' @param background_signature background signature to be used. If not provided, a warning is thrown and an initial value for it is 
#' estimated by NMF. If beta is not NULL, this parameter is ignored.
#' @param normalize_counts if true, the input count matrix x is normalize such that the patients have the same number of mutation.
#' @param nmf_runs number of iteration (minimum 1) of NMF to be performed for a robust estimation of starting beta. If beta is not NULL, 
#' this parameter is ignored.
#' @param lambda_values value of LASSO to be used for beta between 0 and 1. This value should be greater than 0. 1 is the value of LASSO 
#' that would shrink all the signatures to 0 within one step. The higher lambda_values is, the sparser are the resulting signatures, 
#' but too large values may result in a reduced fit of the observed counts.
#' @param iterations Number of iterations to be performed. Each iteration corresponds to a first step where beta is fitted 
#' and a second step where alpha is fitted.
#' @param max_iterations_lasso Number of maximum iterations to be performed during the sparsification via Lasso.
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode, 
#' this parameter needs to be set to either NA or NULL.
#' @param seed Seed for reproducibility.
#' @param verbose boolean; Shall I print all messages?
#' @param log_file log file where to print outputs when using parallel. If parallel execution is disabled, this parameter is ignored.
#' @return A list corresponding to results of the function nmfLasso for each value of lambda to be tested. This function allows 
#' to test a good range of lambda values for beta to be considered. One should keep in mind that too small values generate dense solution, 
#' while too high ones leads to poor fit. This behavior is resampled in the values of loglik_progression, which should be increasing: 
#' too small values of lambda results in unstable log-likelihood through the iterations, while too large values make log-likelihood 
#' drop. 
#' @export lambdaRangeBetaEvaluation
#' @import NMF
#' @import nnlasso
#' @import nnls
#' @import parallel
#'
"lambdaRangeBetaEvaluation" <- function( x, K = 5, beta = NULL, background_signature = NULL, normalize_counts = TRUE, nmf_runs = 10, lambda_values = c(0.01, 0.05, 0.10, 0.20), iterations = 30, max_iterations_lasso = 10000, num_processes = Inf, seed = NULL, verbose = TRUE, log_file = "" ) {
    
    # set the seed
    set.seed(seed)

    # check the input parameters
    if(K<2) {
        if(verbose) {
            cat("The minimum value of K is 2...","\n")
        }
        K <- 2
    }
    if(nmf_runs<1) {
        if(verbose) {
            cat("The minimum value of nmf_runs is 1...","\n")
        }
        nmf_runs <- 1
    }
    x <- as.matrix(x)
    if(normalize_counts) {
        x_not_normalized <- x
        x <- (x/rowSums(x))*2500
    }
    lambda_values <- sort(unique(lambda_values))

    # compute the initial values of beta
    if(is.null(beta)) {

        if(verbose) {
            cat("Computing the initial values of beta by standard NMF...","\n")
        }
        
        # add a signature to beta (leading to K signatures in total) to explicitly model background noise
        if(is.null(background_signature)) {
            warning("No background signature has been specified...")
            background_signature_beta <- basis(nmf(t(x),rank=1,nrun=nmf_runs))
            background_signature_beta <- t(background_signature_beta)
            background_signature_beta <- background_signature_beta / rowSums(background_signature_beta)
            background_signature <- as.numeric(background_signature_beta[1,])
        }
        else {
            background_signature <- as.numeric(background_signature)
            background_signature <- background_signature / sum(background_signature)
            background_signature_beta <- t(background_signature)
        }

        background_signature_alpha <- matrix(0,nrow=dim(x)[1],ncol=1)
        for(i in 1:dim(x)[1]) {
            background_signature_alpha[i,] <- nnls(t(background_signature_beta),as.vector(x[i,]))$x
        }
        curr_x <- x - (background_signature_alpha %*% background_signature_beta)
        curr_x[curr_x<0] <- 0
        if(length(which(rowSums(t(curr_x))==0))>0) {
            invalid_cols <- as.numeric(which(rowSums(t(curr_x))==0))
            for(inv_cols in invalid_cols) {
                curr_x[sample(1:length(curr_x[,inv_cols]),size=1),inv_cols] <- 1e-05
            }
        }
        beta <- basis(nmf(t(curr_x),rank=(K-1),nrun=nmf_runs))
        beta <- t(beta)
        beta <- rbind(background_signature,beta)

    }
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
    lambda_results <- array(list(),c(length(K),length(lambda_values)))
    rownames(lambda_results) <- paste0(as.character(K),"_signatures")
    colnames(lambda_results) <- paste0(as.character(lambda_values),"_lambda")
    
    if(verbose) {
        cat("Performing estimation of lambda range for beta...","\n")
    }
    
    # perform signature discovery for all the values of lambda
    if(is.null(parallel)) {
        set.seed(round(runif(1)*100000))
        cont <- 0
        for(l in lambda_values) {

            # perform the inference
            curr_results <- nmfLasso(x = x, 
                                     K = K, 
                                     beta = beta, 
                                     normalize_counts = normalize_counts, 
                                     lambda_rate_alpha = 0.0, 
                                     lambda_rate_beta = l, 
                                     iterations = iterations, 
                                     max_iterations_lasso = max_iterations_lasso, 
                                     seed = round(runif(1)*100000), 
                                     verbose = FALSE)
                                     
            # save the results for the current configuration
            cont <- cont + 1

            # rescale alpha to the original magnitude
            if(normalize_counts) {
                curr_results$starting_alpha <- curr_results$starting_alpha * (rowSums(x_not_normalized)/2500)
                curr_results$alpha <- curr_results$alpha * (rowSums(x_not_normalized)/2500)
            }

            lambda_results[[1,cont]] <- curr_results
                
            if(verbose) {
                cat("Progress",paste0(round((cont/(length(K)*length(lambda_values)))*100,digits=3),"%..."),"\n")
            }
            
        }
    }
    else {
        res_clusterEvalQ <- clusterEvalQ(parallel,library("nnls"))
        res_clusterEvalQ <- clusterEvalQ(parallel,library("nnlasso"))
        clusterExport(parallel,varlist=c("x","K","beta","normalize_counts","iterations","max_iterations_lasso"),envir=environment())
        clusterExport(parallel,c('nmfLasso','nmfLassoDecomposition','basis','nmf'),envir=environment())
        clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
        lambda_results_res <- parLapply(parallel,lambda_values,function(l) {

            # perform the inference
            curr_results <- nmfLasso(x = x, 
                                     K = K, 
                                     beta = beta, 
                                     normalize_counts = normalize_counts, 
                                     lambda_rate_alpha = 0.0, 
                                     lambda_rate_beta = l, 
                                     iterations = iterations, 
                                     max_iterations_lasso = max_iterations_lasso, 
                                     seed = round(runif(1)*100000), 
                                     verbose = FALSE)

            return(curr_results)

        })
        for(i in 1:ncol(lambda_results)) {
            # rescale alpha to the original magnitude
            if(normalize_counts) {
                lambda_results_res[[i]]$starting_alpha <- lambda_results_res[[i]]$starting_alpha * (rowSums(x_not_normalized)/2500)
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

#' Perform the discovery of K somatic mutational signatures given a set of observed counts x.
#'
#' @examples
#' data(patients)
#' data(starting_betas_example)
#' beta = starting_betas_example[["5_signatures","Value"]]
#' res = nmfLasso(x=patients[1:100,],
#'      K=5,
#'      beta=beta,
#'      lambda_rate_alpha=0.05,
#'      lambda_rate_beta=0.05,
#'      iterations=5,
#'      seed=12345)
#'
#' @title nmfLasso
#' @param x count matrix for a set of n patients and 96 trinucleotides.
#' @param K numeric value (minimum 2) indicating the number of signatures to be discovered.
#' @param beta starting beta for the estimation. If it is NULL, starting beta is estimated by NMF.
#' @param background_signature background signature to be used. If not provided, a warning is thrown and an initial value for it is 
#' estimated by NMF. If beta is not NULL, this parameter is ignored.
#' @param normalize_counts if true, the input count matrix x is normalize such that the patients have the same number of mutation.
#' @param nmf_runs number of iteration (minimum 1) of NMF to be performed for a robust estimation of starting beta. If beta is not NULL, 
#' this parameter is ignored.
#' @param lambda_rate_alpha value of LASSO to be used for alpha between 0 and 1. This value should be greater than 0. 1 is the value of LASSO 
#' that would shrink all the exposure values to 0 within one step. The higher lambda_rate_alpha is, the sparser are the resulting exposure values, 
#' but too large values may result in a reduced fit of the observed counts.
#' @param lambda_rate_beta value of LASSO to be used for beta between 0 and 1. This value should be greater than 0. 1 is the value of LASSO 
#' that would shrink all the signatures to 0 within one step. The higher lambda_rate_beta is, the sparser are the resulting signatures, 
#' but too large values may result in a reduced fit of the observed counts.
#' @param iterations Number of iterations to be performed. Each iteration corresponds to a first step where beta is fitted 
#' and a second step where alpha is fitted.
#' @param max_iterations_lasso Number of maximum iterations to be performed during the sparsification via Lasso.
#' @param seed Seed for reproducibility.
#' @param verbose boolean; Shall I print all messages?
#' @return A list with the discovered signatures. It includes 6 elements: 
#'              alpha: matrix of the discovered exposure values
#'              beta: matrix of the discovered signatures
#'              starting_alpha: initial alpha on which the method has been applied
#'              starting_beta: initial beta on which the method has been applied
#'              loglik_progression: log-likelihood values during the iterations. This values should be increasing, if not the selected value of lambda is too high
#'              best_loglik: log-likelihood of the best signatures configuration
#' @export nmfLasso
#' @import NMF
#' @import nnls
#' @import nnlasso
#'
"nmfLasso" <- function( x, K, beta = NULL, background_signature = NULL, normalize_counts = TRUE, nmf_runs = 10, lambda_rate_alpha = 0.05, lambda_rate_beta = 0.05, iterations = 30, max_iterations_lasso = 10000, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)

    # check the input parameters
    if(K<2) {
        if(verbose) {
            cat("The minimum value of K is 2...","\n")
        }
        K <- 2
    }
    if(nmf_runs<1) {
        if(verbose) {
            cat("The minimum value of nmf_runs is 1...","\n")
        }
        nmf_runs <- 1
    }
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
    if(lambda_rate_beta<0) {
        if(verbose) {
            cat("The minimum value of lambda_rate_beta is 0...","\n")
        }
        lambda_rate_beta <- 0
    }
    if(lambda_rate_beta>1) {
        if(verbose) {
            cat("The maximum value of lambda_rate_beta is 1...","\n")
        }
        lambda_rate_beta <- 1
    }
    x <- as.matrix(x)
    if(normalize_counts) {
        x_not_normalized <- x
        x <- (x/rowSums(x))*2500
    }

    # compute the initial values of beta
    if(is.null(beta)) {

        if(verbose) {
            cat("Computing the initial values of beta by standard NMF...","\n")
        }
        
        # add a signature to beta (leading to K signatures in total) to explicitly model background noise
        if(is.null(background_signature)) {
            warning("No background signature has been specified...")
            background_signature_beta <- basis(nmf(t(x),rank=1,nrun=nmf_runs))
            background_signature_beta <- t(background_signature_beta)
            background_signature_beta <- background_signature_beta / rowSums(background_signature_beta)
            background_signature <- as.numeric(background_signature_beta[1,])
        }
        else {
            background_signature <- as.numeric(background_signature)
            background_signature <- background_signature / sum(background_signature)
            background_signature_beta <- t(background_signature)
        }

        background_signature_alpha <- matrix(0,nrow=dim(x)[1],ncol=1)
        for(i in 1:dim(x)[1]) {
            background_signature_alpha[i,] <- nnls(t(background_signature_beta),as.vector(x[i,]))$x
        }
        curr_x <- x - (background_signature_alpha %*% background_signature_beta)
        curr_x[curr_x<0] <- 0
        if(length(which(rowSums(t(curr_x))==0))>0) {
            invalid_cols <- as.numeric(which(rowSums(t(curr_x))==0))
            for(inv_cols in invalid_cols) {
                curr_x[sample(1:length(curr_x[,inv_cols]),size=1),inv_cols] <- 1e-05
            }
        }
        beta <- basis(nmf(t(curr_x),rank=(K-1),nrun=nmf_runs))
        beta <- t(beta)
        beta <- rbind(background_signature,beta)

    }
    beta <- beta / rowSums(beta)
    rownames(beta) <- c("Background",paste0("S",1:(nrow(beta)-1)))
    colnames(beta) <- colnames(x)
    
    if(verbose) {
        cat("Performing the discovery of the signatures by NMF with Lasso...","\n")
    }
    
    # perform the discovery of the signatures
    results = tryCatch({
        results = nmfLassoDecomposition(x,beta,lambda_rate_alpha,lambda_rate_beta,iterations,max_iterations_lasso,round(runif(1)*100000),verbose)
        # rescale alpha to the original magnitude
        if(normalize_counts) {
            results$starting_alpha <- results$starting_alpha * (rowSums(x_not_normalized)/2500)
            results$alpha <- results$alpha * (rowSums(x_not_normalized)/2500)
        }
        return(results)
    }, error = function(e) {
        warning("Lasso did not converge, you should try a lower value of lambda! Current settings: K = ",K,", lambda_rate_alpha = ",lambda_rate_alpha,", lambda_rate_beta = ",lambda_rate_beta,"...")
        return(list(alpha=NA,beta=NA,starting_alpha=NA,starting_beta=NA,loglik_progression=rep(NA,iterations),best_loglik=NA))
    })
    
    return(results)
    
}

# perform de novo discovery of somatic mutational signatures using NMF with Lasso to ensure sparsity
"nmfLassoDecomposition" <- function( x, starting_beta, lambda_rate_alpha = 0.05, lambda_rate_beta = 0.05, iterations = 30, max_iterations_lasso = 10000, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)
    
    # n is the number of observations in x, i.e., the patients
    n <- dim(x)[1]
    
    # J is the number of trinucleotides, i.e., 96 categories
    J <- dim(x)[2]
    
    # K is the number of signatures to be called
    K <- dim(starting_beta)[1]
    
    # initialize beta
    starting_beta <- starting_beta / rowSums(starting_beta)
    rownames(starting_beta) <- c("Background",paste0("S",1:(nrow(starting_beta)-1)))
    colnames(starting_beta) <- colnames(x)
    beta <- starting_beta
    
    # initialize alpha
    alpha <- matrix(0,nrow=n,ncol=K)
    rownames(alpha) <- 1:nrow(alpha)
    colnames(alpha) <- rownames(beta)
    for(i in 1:n) {
        alpha[i,] <- nnls(t(beta),as.vector(x[i,]))$x
    }
    starting_alpha <- alpha

    # rescale the two matrices of a factor of 50 for numerical reasons
    beta <- beta * 50
    for(i in 1:n) {
        alpha[i,] <- nnls(t(beta),as.vector(x[i,]))$x
    }
    
    # structure where to save the log-likelihood at each iteration
    loglik <- rep(NA,iterations)

    # structure where to save the lambda values for alpha and beta
    lambda_values_alpha <- rep(NA,n) # sparsity over the patients (reduce the number of signatures per patient)
    lambda_values_beta <- rep(NA,J) # sparsity over the signatures (produce sparse signatures)

    if(verbose) {
        cat("Performing a total of",iterations,"iterations...","\n")
    }

    # start the learning procedure
    best_loglik <- NA
    best_alpha <- NA
    best_beta <- NA
    for(i in 1:iterations) {
        
        # initialize the value of the log-likelihood for the current iteration
        loglik[i] <- 0

        # configuration 1 (turn off sparsity): repeat a 2 step algorithm iteratively, where 
        # first beta is estimated by Non-Negative Linear Least Squares and, then, alpha is also estimated by Non-Negative Linear Least Squares 
        if(lambda_rate_alpha==0&&lambda_rate_beta==0) {

            # update beta by Non-Negative Linear Least Squares
            for(k in 1:J) {
                 # the first signature represents the background model, thus it is not changed
                beta[2:K,k] <- nnls(alpha[,2:K,drop=FALSE],as.vector(x[,k]-(alpha[,1]*beta[1,k])))$x
            }

            # update alpha by Non-Negative Linear Least Squares
            for(j in 1:n) {
                alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
            }

            # update the log-likelihood for the current iteration
            loglik[i] <- -sum((x - alpha %*% beta)^2)

        }

        # configuration 2 (sparsity only on beta): repeat a 2 step algorithm iteratively, where 
        # first beta is estimated by Non-Negative Lasso and, then, alpha is estimated by Non-Negative Linear Least Squares 
        if(lambda_rate_alpha==0&&lambda_rate_beta>0) {
            
            # update beta by Non-Negative Lasso
            for(k in 1:J) {
                
                # compute independently for each trinucleotide the difference between the observed counts, i.e., x, 
                # and the counts predicted by considering only the first signature (i.e., which represents the background model) 
                target <- x[,k] - (alpha[,1] * beta[1,k])
                
                # estimate beta for the remaining signatues by Non-Negative Lasso to ensure sparsity
                if(is.na(lambda_values_beta[k])) {
                    max_lambda_value <- max(abs(t(alpha[,2:K,drop=FALSE]) %*% target))
                    lambda_values_beta[k] <- max_lambda_value * lambda_rate_beta
                }
                beta[2:K,k] <- as.vector(nnlasso(x = alpha[,2:K,drop=FALSE], 
                                                y = target, 
                                                family = "normal", 
                                                lambda = lambda_values_beta[k], 
                                                intercept = FALSE, 
                                                normalize = FALSE, 
                                                tol = 1e-05, 
                                                maxiter = max_iterations_lasso, 
                                                path = FALSE)$coef[2,])
                
            }

            # update alpha by Non-Negative Linear Least Squares
            for(j in 1:n) {
                alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
            }
            
            # update the log-likelihood for the current iteration
            for(k in 1:J) {
                curr_loglik <- -sum((x[,k] - alpha %*% beta[,k])^2) - (lambda_values_beta[k] * sum(beta[2:K,k]))
                loglik[i] <- loglik[i] + curr_loglik
            }

        }

        # configuration 3 (sparsity only on alpha): repeat a 2 step algorithm iteratively, where 
        # first beta is estimated by Non-Negative Linear Least Squares and, then, alpha is estimated by Non-Negative Lasso 
        if(lambda_rate_alpha>0&&lambda_rate_beta==0) {

            # update beta by Non-Negative Linear Least Squares
            for(k in 1:J) {
                # the first signature represents the background model, thus it is not changed
                beta[2:K,k] <- nnls(alpha[,2:K,drop=FALSE],as.vector(x[,k]-(alpha[,1]*beta[1,k])))$x
            }

            # update alpha by Non-Negative Lasso
            for(j in 1:n) {

                # compute independently for each patient the counts to be approximted
                target <- x[j,]
                
                # estimate alpha by Non-Negative Lasso to ensure sparsity
                if(is.na(lambda_values_alpha[j])) {
                    max_lambda_value <- max(abs(beta %*% target))
                    lambda_values_alpha[j] <- max_lambda_value * lambda_rate_alpha
                }
                alpha[j,] <- as.vector(nnlasso(x = t(beta), 
                                                y = target, 
                                                family = "normal", 
                                                lambda = lambda_values_alpha[j], 
                                                intercept = FALSE, 
                                                normalize = FALSE, 
                                                tol = 1e-05, 
                                                maxiter = max_iterations_lasso, 
                                                path = FALSE)$coef[2,])
                
            }
            
            # update the log-likelihood for the current iteration
            for(j in 1:n) {
                curr_loglik <- -sum((x[j,] - alpha[j,] %*% beta)^2) - (lambda_values_alpha[j] * sum(alpha[j,]))
                loglik[i] <- loglik[i] + curr_loglik
            }

        }

        # configuration 4 (sparsity on both alpha and beta): repeat a 2 step algorithm iteratively, where 
        # first beta is estimated by Non-Negative Lasso and, then, also alpha is estimated by Non-Negative Lasso 
        if(lambda_rate_alpha>0&&lambda_rate_beta>0) {
            
            # update beta by Non-Negative Lasso
            for(k in 1:J) {
                
                # compute independently for each trinucleotide the difference between the observed counts, i.e., x, 
                # and the counts predicted by considering only the first signature (i.e., which represents the background model) 
                target <- x[,k] - (alpha[,1] * beta[1,k])
                
                # estimate beta for the remaining signatues by Non-Negative Lasso to ensure sparsity
                if(is.na(lambda_values_beta[k])) {
                    max_lambda_value <- max(abs(t(alpha[,2:K,drop=FALSE]) %*% target))
                    lambda_values_beta[k] <- max_lambda_value * lambda_rate_beta
                }
                beta[2:K,k] <- as.vector(nnlasso(x = alpha[,2:K,drop=FALSE], 
                                                y = target, 
                                                family = "normal", 
                                                lambda = lambda_values_beta[k], 
                                                intercept = FALSE, 
                                                normalize = FALSE, 
                                                tol = 1e-05, 
                                                maxiter = max_iterations_lasso, 
                                                path = FALSE)$coef[2,])
                
            }
            
            # update alpha by Non-Negative Lasso
            for(j in 1:n) {

                # compute independently for each patient the counts to be approximted
                target <- x[j,]
                
                # estimate alpha by Non-Negative Lasso to ensure sparsity
                if(is.na(lambda_values_alpha[j])) {
                    max_lambda_value <- max(abs(beta %*% target))
                    lambda_values_alpha[j] <- max_lambda_value * lambda_rate_alpha
                }
                alpha[j,] <- as.vector(nnlasso(x = t(beta), 
                                                y = target, 
                                                family = "normal", 
                                                lambda = lambda_values_alpha[j], 
                                                intercept = FALSE, 
                                                normalize = FALSE, 
                                                tol = 1e-05, 
                                                maxiter = max_iterations_lasso, 
                                                path = FALSE)$coef[2,])
                
            }
            
            # update the log-likelihood for the current iteration
            for(k in 1:J) {
                for(j in 1:n) {
                    curr_loglik <- -((x[j,k] - alpha[j,] %*% beta[,k])^2) - (lambda_values_alpha[j] * sum(alpha[j,])) - (lambda_values_beta[k] * sum(beta[2:K,k]))
                    loglik[i] <- loglik[i] + curr_loglik
                }
            }

        }
        
        # save the results at maximum log-likelihood
        if(is.na(best_loglik)||loglik[i]>best_loglik) {
            best_alpha <- alpha
            best_beta <- beta
            best_loglik <- loglik[i]
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
            warning("The likelihood is not increasing, you should try a lower value of lambda! Current settings: K = ",K,", lambda_rate_alpha = ",lambda_rate_alpha,", lambda_rate_beta = ",lambda_rate_beta,"...")
        }
    }
    
    # normalize the rate of the signatures to sum to 1
    beta <- beta / rowSums(beta)
        
    # final computation of alpha using the normalized version of beta
    if(lambda_rate_alpha>0) {

        # update alpha by Non-Negative Lasso
        for(j in 1:n) {

            # compute independently for each patient the counts to be approximted
            target <- x[j,]
            
            # estimate alpha by Non-Negative Lasso to ensure sparsity
            alpha[j,] <- as.vector(nnlasso(x = t(beta), 
                                            y = target, 
                                            family = "normal", 
                                            lambda = (max(abs(beta %*% target)) * lambda_rate_alpha), 
                                            intercept = FALSE, 
                                            normalize = FALSE, 
                                            tol = 1e-05, 
                                            maxiter = max_iterations_lasso, 
                                            path = FALSE)$coef[2,])
            
        }
        
    }
    else {

        # update alpha by Non-Negative Linear Least Squares
        for(j in 1:n) {
            alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
        }

    }
    
    # save the results
    results <- list(alpha=alpha,beta=beta,starting_alpha=starting_alpha,starting_beta=starting_beta,loglik_progression=loglik,best_loglik=best_loglik)
    
    return(results)
    
}
