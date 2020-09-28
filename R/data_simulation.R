#' Simulate a count matrix for a single cell type based on a copula model
#'
#' @param copula_result A list that contains the parameters of a copula model.
#' @param n             An integer value that indicates the number of cells to generate.
#' @param marginal      A character string that indicates whether the generated values should
#'                      stay as discrete or switch to continuous. Default value is 'nb', which
#'                      should be used for generating a count marix. The alternative 'Gamma' is
#'                      only needed when this function is being called by other functions that
#'                      generate data with a user-specified sequencing depth. Normally, users
#'                      do not need to change this value.
#' @return A matrix of shape p by n that contains the simulated count values. p is derived from
#' \code{copula_result}
#'
# @export
simulate_count_copula <- function(copula_result, n = 100,
                                  marginal = c('nb', 'Gamma')){
  marginal <- match.arg(marginal)

  p1 <- length(copula_result$gene_sel1)
  if(p1 > 0){
    result1 <- mvrnorm(n = n, mu = rep(0.0, p1), Sigma = copula_result$cov_mat)
    result1 <- matrix(result1, nrow = n)
    result2 <- apply(result1, 2, pnorm)
    result2 <- matrix(result2, nrow = n)
  }
  p2 <- length(copula_result$gene_sel2)
  if(marginal == 'nb'){
    if(p1 > 0){
      result31 <- t(sapply(1:p1, function(iter){
        param <- copula_result$marginal_param1[iter, ]
        qnbinom(pmax(0.0, result2[, iter] - param[1]) / (1-param[1]),
                size = param[2], mu = param[3])
      }))
    }
    if(p2 > 0){
      result32 <- t(sapply(1:p2, function(iter){
        param <- copula_result$marginal_param2[iter, ]
        rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
      }))
    }
  }else if(marginal == 'Gamma'){
    if(p1 > 0){
      result31 <- t(sapply(1:p1, function(iter){
        param <- copula_result$marginal_param1[iter, ]
        qgamma(max(0.0, result2[, iter] - param[1]), shape = param[2], scale = param[3] / param[2])
      }))
    }
    if(p2 > 0){
      result32 <- t(sapply(1:p2, function(iter){
        param <- copula_result$marginal_param2[iter, ]
        rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
      }))
    }
  }

  result <- matrix(0, nrow = p1 + p2 + length(copula_result$gene_sel3), ncol = n)
  if(p1 > 0){
    result[copula_result$gene_sel1, ] <- result31
  }
  if(p2 > 0){
    result[copula_result$gene_sel2, ] <- result32
  }
  result
}



#' Simulate a count matrix for a single cell type based on a (w/o copula model)
#'
#' @param model_params A list that contains the model parameters (can be either the copula
#'                     model or the (w/o copula) model).
#' @inheritParams simulate_count_copula
#' @return A matrix of shape p by n that contains the simulated count values. p is derived from
#' \code{model_params}.
#'
# @export
simulate_count_ind <- function(model_params, n = 100,
                               marginal = c('nb', 'Gamma')){
  marginal <- match.arg(marginal)

  if(model_params$sim_method == 'copula' || 'gene_sel3' %in% names(model_params)){
    p1 <- length(model_params$gene_sel1)
    p2 <- length(model_params$gene_sel2)
    if(marginal == 'nb'){
      if(p1 > 0){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
        }))
      }
      if(p2 > 0){
        result32 <- t(sapply(1:p2, function(iter){
          param <- model_params$marginal_param2[iter, ]
          rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
        }))
      }
    }else if(marginal == 'Gamma'){
      if(p1 > 0){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
        }))
      }
      if(p2 > 0){
        result32 <- t(sapply(1:p2, function(iter){
          param <- model_params$marginal_param2[iter, ]
          rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
        }))
      }
    }
    result <- matrix(0, nrow = p1 + p2 + length(model_params$gene_sel3), ncol = n)
    if(p1 > 0){
      result[model_params$gene_sel1, ] <- result31
    }
    if(p2 > 0){
      result[model_params$gene_sel2, ] <- result32
    }
  }else{
    p1 <- length(model_params$gene_sel1)
    p2 <- length(model_params$gene_sel2)
    result <- matrix(0, nrow = p1 + p2, ncol = n)
    if(p1 > 0){
      if(marginal == 'nb'){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
        }))
      }else if(marginal == 'Gamma'){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
        }))
      }
      result[model_params$gene_sel1, ] <- result31
    }
  }
  result
}



#' Simulate a count matrix for multiple cell types
#'
#' @param model_params    A list with the same length as \code{cell_type_prop} that contains
#'                        the fitted model as each of its element (can be either the copula
#'                        model or the (w/o copula) model).
#' @param total_count_new The (expected) total number of reads or UMIs in the simulated count
#'                        matrix.
#' @param n_cell_new      The total number of cells in the simulated count matrix.
#' @param cell_type_prop  The cell type proportion in the simulated count matrix.
#' @param total_count_old The total number of reads or UMIs in the original count matrix where
#'                        \code{model_params} was fitted.
#' @param n_cell_old      The The total number of cells in the original count matrix where
#'                        \code{model_params} was fitted.
#' @param sim_method      Specification of the type of model for data simulation.
#'                        Default value is 'copula', which selects the copula model.
#'                        'ind' will select the (w/o copula) model.
#' @param reseq_method    Specification of how the new count matrix should be derived under the
#'                        new sequencing depth.
#'                        Default is 'mean_scale', which scales the original parameters and
#'                        then simulate new data.
#'                        'multinomial' will do a resampling. It ensures that the simulated
#'                        count matrix has the exact total number of reads as specified in
#'                        \code{total_count_new}.
#' @param cell_sample     Logical, whether cells for each cell type should be sampled from a
#'                        multinomial distribution or follows the exact same proportion as
#'                        specified in \code{cell_type_prop}.
#' @return A matrix of shape p by n that contains the simulated count values. p is derived from
#' \code{model_params}.
#'
#' @export
simulate_count_scDesign2 <- function(model_params, total_count_new = NULL, n_cell_new = NULL,
                                     cell_type_prop = NULL, total_count_old = NULL,
                                     n_cell_old = NULL, sim_method = c('copula', 'ind'),
                                     reseq_method = c('mean_scale', 'multinomial'),
                                     cell_sample = TRUE){
  sim_method <- match.arg(sim_method)
  reseq_method <- match.arg(reseq_method)

  n_cell_vec <- sapply(model_params, function(x) x$n_cell)
  n_read_vec <- sapply(model_params, function(x) x$n_read)
  if(is.null(total_count_new)) total_count_new <- sum(n_read_vec)
  if(is.null(n_cell_new))      n_cell_new      <- sum(n_cell_vec)
  if(is.null(cell_type_prop))  cell_type_prop  <- n_cell_vec
  if(is.null(total_count_old)) total_count_old <- sum(n_read_vec)
  if(is.null(n_cell_old))      n_cell_old      <- sum(n_cell_vec)

  if(length(model_params)!=cell_type_prop){
    stop('Cell type proportion should have the same length as the number of models.')
  }

  n_cell_type <- length(cell_type_prop)
  if(cell_sample == TRUE){
    n_cell_each <- as.numeric(rmultinom(1, size = n_cell_new, prob = cell_type_prop))
  }else{
    cell_type_prop <- cell_type_prop / sum(cell_type_prop)
    n_cell_each <- round(cell_type_prop * n_cell_new)
    if(sum(n_cell_each) != n_cell_new){
      idx <- sample(n_cell_type, size = 1)
      n_cell_each[idx] <- n_cell_each[idx] + n_cell_new - sum(n_cell_each)
    }
  }

  p <- length(model_params[[1]]$gene_sel1) + length(model_params[[1]]$gene_sel2) +
    length(model_params[[1]]$gene_sel3)
  new_count <- matrix(0, nrow = p, ncol = n_cell_new)
  if(reseq_method == 'mean_scale'){
    n_cell_each
    if(length(total_count_new) == 1){
      r <- rep(total_count_new / sum((total_count_old / n_cell_old) * n_cell_each),
               n_cell_type)
    }else{
      r <- (total_count_new / n_cell_new) / (total_count_old / n_cell_old)
    }
    for(iter in 1:n_cell_type)
      if(n_cell_each[iter] > 0){
        ulim <- sum(n_cell_each[1:iter])
        llim <- ulim - n_cell_each[iter] + 1
        params_new <- model_params[[iter]]
        params_new$marginal_param1[, 3] <- params_new$marginal_param1[, 3] * r[iter]
        if(sim_method == 'copula'){
          params_new$marginal_param2[, 3] <- params_new$marginal_param2[, 3] * r[iter]
          new_count[, llim:ulim] <- simulate_count_copula(params_new, n = n_cell_each[iter],
                                                          marginal = 'nb')
        }else if(sim_method == 'ind'){
          new_count[, llim:ulim] <- simulate_count_ind(params_new, n = n_cell_each[iter],
                                                       marginal = 'nb')
        }
      }
    if(is.null(names(model_params))){
      colnames(new_count) <- unlist(lapply(1:n_cell_type, function(x){rep(x, n_cell_each[x])}))
    }else{
      colnames(new_count) <- unlist(lapply(1:n_cell_type, function(x){
        rep(names(model_params)[x], n_cell_each[x])}))
    }
    return(new_count)
  }else if(reseq_method == 'multinomial'){
    for(iter in 1:n_cell_type){
      ulim <- sum(n_cell_each[1:iter])
      llim <- ulim - n_cell_each[iter] + 1
      if(sim_method == 'copula'){
        new_count[, llim:ulim] <- simulate_count_copula(model_params[[iter]],
                                                        n = n_cell_each[iter], marginal = 'Gamma')
      }else if(sim_method == 'ind'){
        new_count[, llim:ulim] <- simulate_count_ind(model_params[[iter]],
                                                     n = n_cell_each[iter], marginal = 'Gamma')
      }
    }

    new_count[which(is.infinite(new_count))] <- 0
    new_count[which(is.na(new_count))] <- 0

    bam_file <- sample(x = p*n_cell_new, size = total_count_new,
                       replace = TRUE, prob = as.vector(new_count))
    hist_result <- hist(bam_file, breaks = 0:(n_cell_new*p), plot = FALSE)
    result <- matrix(hist_result$counts, nrow = nrow(new_count))
    if(is.null(names(model_params))){
      colnames(result) <- unlist(lapply(1:n_cell_type, function(x){rep(x, n_cell_each[x])}))
    }else{
      colnames(result) <- unlist(lapply(1:n_cell_type, function(x){
        rep(names(model_params)[x], n_cell_each[x])}))
    }
    return(result)
  }
}


