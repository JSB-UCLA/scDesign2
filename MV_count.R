library(MASS)
library(pscl)
fit_marginals <- function(x, marginal = NULL, pval_cutoff = 0.05, epsilon = 1e-5,
                          jitter = TRUE, DT = TRUE)
{
  p <- nrow(x)
  n <- ncol(x)
  
  if(is.null(marginal))
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v)
      {
        mle_Poisson <- glm(gene ~ 1, family = poisson)
        tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }
      else
      {
        mle_NB <- glm.nb(gene ~ 1)
        if(min(gene) > 0)
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        else
          tryCatch({
            mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
            chisq_val <- 2 * (logLik(mle_ZINB) - logLik(mle_NB))
            pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
            if(pvalue < pval_cutoff)
              c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
            else
              c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          },
          error = function(cond){
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  else if(marginal == 'nb')
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v)
      {
        mle_Poisson <- glm(gene ~ 1, family = poisson)
        tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }
      else
      {
        mle_NB <- glm.nb(gene ~ 1)
        c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
      }
    }))
  else if(marginal == 'zinb')
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v)
      {
        mle_Poisson <- glm(gene ~ 1, family = poisson)
        tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }
      else
      {      
        if(min(gene) > 0)
        {
          mle_NB <- glm.nb(gene ~ 1)
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        }
        else
          tryCatch({
            mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
            c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
          },
          error = function(cond){
            mle_NB <- glm.nb(gene ~ 1)
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  else
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v)
      {
        mle_Poisson <- glm(gene ~ 1, family = poisson)
        tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }
      else
      {
        mle_NB <- glm.nb(gene ~ 1)
        if(min(gene) > 0)
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        else
          tryCatch({
            mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
            chisq_val <- 2 * (logLik(mle_ZINB) - logLik(mle_NB))
            pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
            if(pvalue < pval_cutoff)
              c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
            else
              c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          },
          error = function(cond){
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  # params <- t(apply(x, 1, function(gene){
  #   m <- mean(gene)
  #   v <- var(gene)
  #   if(m < v)
  #     c(m^2/(v-m), m)
  #   else
  #     c(0.001, m)
  # }))
  
  if(DT)
      u <- t(sapply(1:p, function(iter){
        param <- params[iter, ]
        gene <- x[iter, ]
        prob0 <- param[1]
        u1 <- prob0 + (1 - prob0) * pnbinom(gene, size = param[2], mu = param[3])
        u2 <- (prob0 + (1 - prob0) * pnbinom(gene - 1, size = param[2], mu = param[3])) *
          as.integer(gene > 0)
        if(jitter)
          v <- runif(n)
        else
          v <- rep(0.5, n)
        r <- u1 * v + u2 * (1 - v)
        idx_adjust <- which(1-r < epsilon)
        # print(idx_adjust)
        r[idx_adjust] <- r[idx_adjust] - epsilon
        # print(r[idx_adjust])
        
        idx_adjust <- which(r < epsilon)
        r[idx_adjust] <- r[idx_adjust] + epsilon
        
        r
      }))
  else
    u <- NULL
  
  return(list(params = params, u = u))
}



fit_Gaussian_copula <- function(x, marginal = NULL, jitter = TRUE, zp_cutoff = 0.8,
                                min_nonzero_num = 2)
{
  n <- ncol(x)
  p <- nrow(x)
  
  gene_zero_prop <- apply(x, 1, function(y){
    sum(y < 1e-5) / n
  })
  
  gene_sel1 <- which(gene_zero_prop < zp_cutoff)
  gene_sel2 <- which(gene_zero_prop < 1.0 - min_nonzero_num/n &
                       gene_zero_prop >= zp_cutoff)
  gene_sel3 <- (1:p)[-c(gene_sel1, gene_sel2)]
  
  if(length(gene_sel1) > 0)
  {
    marginal_result1 <- fit_marginals(x[gene_sel1, ], jitter = jitter, DT = TRUE)
    quantile_normal <- qnorm(marginal_result1$u)
    cov_mat <- cor(t(quantile_normal))
  }
  else
  {
    cov_mat = NULL
    marginal_result1 = NULL
  }
  
  if(length(gene_sel2) > 0)
    marginal_result2 <- fit_marginals(x[gene_sel2, ], DT = FALSE)
  else
    marginal_result2 = NULL
  return(list(cov_mat = cov_mat, marginal_param1 = marginal_result1$params,
              marginal_param2 = marginal_result2$params,
              gene_sel1 = gene_sel1, gene_sel2 = gene_sel2, gene_sel3 = gene_sel3,
              zp_cutoff = zp_cutoff))#,u = marginal_result$u))
}

simulate_new_count <- function(copula_result, n = 100,
                               marginal = c('negative binomial', 'Gamma'))
{
  p1 <- length(copula_result$gene_sel1)
  result1 <- mvrnorm(n = n, mu = rep(0.0, p1), Sigma = copula_result$cov_mat)
  result1 <- matrix(result1, nrow = n)
  result2 <- apply(result1, 2, pnorm)
  result2 <- matrix(result2, nrow = n)
  # str(result2)
  # print(head(result2))
  p2 <- length(copula_result$gene_sel2)
  
  if(marginal == 'negative binomial')
  {
    result31 <- t(sapply(1:p1, function(iter){
      param <- copula_result$marginal_param1[iter, ]
      qnbinom(pmax(0.0, result2[, iter] - param[1]) / (1-param[1]),
              size = param[2], mu = param[3])
    }))
    if(p2 > 0)
      result32 <- t(sapply(1:p2, function(iter){
        param <- copula_result$marginal_param2[iter, ]
        rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
      }))
  }
  
  else if(marginal == 'Gamma')
  {
    result31 <- t(sapply(1:p1, function(iter){
      param <- copula_result$marginal_param1[iter, ]
      qgamma(max(0.0, result2[, iter] - param[1]), shape = param[2], scale = param[3] / param[2])
    }))
    if(p2 > 0)
      result32 <- t(sapply(1:p2, function(iter){
        param <- copula_result$marginal_param2[iter, ]
        rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
      }))
  }
  
  result <- matrix(0, nrow = p1 + p2 + length(copula_result$gene_sel3), ncol = n)
  result[copula_result$gene_sel1, ] <- result31
  if(p2 > 0)
    result[copula_result$gene_sel2, ] <- result32
  result
}


simulate_new_count_ind <- function(copula_result, n = 100)
{
  p1 <- length(copula_result$gene_sel1)
  p2 <- length(copula_result$gene_sel2)
  
  result31 <- t(sapply(1:p1, function(iter){
    param <- copula_result$marginal_param1[iter, ]
    rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
    }))
  if(p2 > 0)
    result32 <- t(sapply(1:p2, function(iter){
      param <- copula_result$marginal_param2[iter, ]
      rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
    }))

  result <- matrix(0, nrow = p1 + p2 + length(copula_result$gene_sel3), ncol = n)
  result[copula_result$gene_sel1, ] <- result31
  if(p2 > 0)
    result[copula_result$gene_sel2, ] <- result32
  result
}



### copula_result is a list that contains one set of parameter for each cell type
simulate_new_count_exp_design <- function(copula_result, total_count_new, n_cell_new,
                                          cell_type_prop = 1,
                                          total_count_old = NULL, n_cell_old = NULL,
                                          method = c('mean_scale', 'multinomial'),
                                          cell_sample = TRUE)
{
  n_cell_type <- length(cell_type_prop)
  if(cell_sample == TRUE)
    n_cell_each <- as.numeric(rmultinom(1, size = n_cell_new, prob = cell_type_prop))
  else
  {
    cell_type_prop <- cell_type_prop / sum(cell_type_prop)
    n_cell_each <- round(cell_type_prop * n_cell_new)
    if(sum(n_cell_each) != n_cell_new)
    {
      idx <- sample(n_cell_type, size = 1)
      n_cell_each[idx] <- n_cell_each[idx] + n_cell_new - sum(n_cell_each)
    }
  }

  p <- length(copula_result[[1]]$gene_sel1) + length(copula_result[[1]]$gene_sel2) +
    length(copula_result[[1]]$gene_sel3)
  new_count <- matrix(0, nrow = p, ncol = n_cell_new)
  if(method == 'multinomial')
  {
    for(iter in 1:n_cell_type)
    {
      ulim <- sum(n_cell_each[1:iter])
      llim <- ulim - n_cell_each[iter] + 1
      new_count[, llim:ulim] <- simulate_new_count(copula_result[[iter]],
                                                   n = n_cell_each[iter], marginal = 'Gamma')
    }
    
    new_count[which(is.infinite(new_count))] <- 0
    new_count[which(is.na(new_count))] <- 0
    
    bam_file <- sample(x = p*n_cell_new, size = total_count_new,
                       replace = TRUE, prob = as.vector(new_count))
    hist_result <- hist(bam_file, breaks = 0:(n_cell_new*p), plot = FALSE)
    result <- matrix(hist_result$counts, nrow = nrow(new_count))
    colnames(result) <- unlist(lapply(1:n_cell_type, function(x){rep(x, n_cell_each[x])}))
    result
  }
  else if(method == 'mean_scale')
  {
    n_cell_each
    if(length(total_count_new) == 1)
      r <- rep(total_count_new / sum((total_count_old / n_cell_old) * n_cell_each),
               n_cell_type)
    else
      r <- (total_count_new / n_cell_new) / (total_count_old / n_cell_old)
    for(iter in 1:n_cell_type)
      if(n_cell_each[iter] > 0)
      {
        ulim <- sum(n_cell_each[1:iter])
        llim <- ulim - n_cell_each[iter] + 1
        params_new <- copula_result[[iter]]
        params_new$marginal_param1[, 3] <- params_new$marginal_param1[, 3] * r[iter]
        params_new$marginal_param2[, 3] <- params_new$marginal_param2[, 3] * r[iter]
        new_count[, llim:ulim] <- simulate_new_count(params_new,
                                                     n = n_cell_each[iter], marginal = 'negative binomial')
      }
    colnames(new_count) <- unlist(lapply(1:n_cell_type, function(x){rep(x, n_cell_each[x])}))
    new_count
  }
}



simulate_new_count_multiple_types <-
  function(copula_result, n_cell_new,
           cell_type_prop = 1, cell_sample = FALSE,
           sim_method = 'ind')
{
  n_cell_type <- length(cell_type_prop)
  if(cell_sample == TRUE)
    n_cell_each <- as.numeric(rmultinom(1, size = n_cell_new, prob = cell_type_prop))
  else
  {
    cell_type_prop <- cell_type_prop / sum(cell_type_prop)
    n_cell_each <- round(cell_type_prop * n_cell_new)
    if(sum(n_cell_each) != n_cell_new)
    {
      idx <- sample(n_cell_type, size = 1)
      n_cell_each[idx] <- n_cell_each[idx] + n_cell_new - sum(n_cell_each)
    }
  }
  
  p <- length(copula_result[[1]]$gene_sel1) + length(copula_result[[1]]$gene_sel2) +
    length(copula_result[[1]]$gene_sel3)
  new_count <- matrix(0, nrow = p, ncol = n_cell_new)

  n_cell_each
  for(iter in 1:n_cell_type)
    if(n_cell_each[iter] > 0)
    {
      ulim <- sum(n_cell_each[1:iter])
      llim <- ulim - n_cell_each[iter] + 1
      params_new <- copula_result[[iter]]
      if(sim_method == 'copula')
        new_count[, llim:ulim] <- 
        simulate_new_count(params_new, n = n_cell_each[iter], marginal = 'negative binomial')
      else if(sim_method == 'ind')
        new_count[, llim:ulim] <-
        simulate_new_count_ind(params_new, n = n_cell_each[iter])
    }
  colnames(new_count) <- unlist(lapply(1:n_cell_type, function(x){rep(x, n_cell_each[x])}))
  new_count
}

simulate_new_count_multiple_types_count <-
  function(traincount, n_cell_new, cell_type_sel,
           cell_type_prop = 1, cell_sample = FALSE,
           sim_method = 'simple')
{
  n_cell_type <- length(cell_type_prop)
  if(cell_sample == TRUE)
    n_cell_each <- as.numeric(rmultinom(1, size = n_cell_new, prob = cell_type_prop))
  else
  {
    cell_type_prop <- cell_type_prop / sum(cell_type_prop)
    n_cell_each <- round(cell_type_prop * n_cell_new)
    if(sum(n_cell_each) != n_cell_new)
    {
      idx <- sample(n_cell_type, size = 1)
      n_cell_each[idx] <- n_cell_each[idx] + n_cell_new - sum(n_cell_each)
    }
  }
  
  n_cell_each
  new_count <- mclapply(1:n_cell_type, function(iter){
    if(n_cell_each[iter] > 0)
    {
      realmat <- traincount[, colnames(traincount) == cell_type_sel[iter]]
      colnames(realmat) <- 1:ncol(realmat)
      if(sim_method == "SPARSim")
        return(get_data_simulated(realmat, sim_method)[[1]][[1]])
      else if(sim_method == "zinb-wave")
      {
        temp <- get_data_simulated(realmat, sim_method)[[1]]
        result_mat <- matrix(0, nrow = nrow(realmat), ncol = ncol(realmat))
        result_mat[temp$gene_sel, ] <- temp[[1]]
        result_mat
      }
      else
        get_data_simulated(realmat, sim_method)[[1]]
    }
    
  }, mc.cores = n_cell_type)
  new_count <- Reduce(cbind, new_count)
  
  n_cell_each <- table(colnames(traincount))[cell_type_sel]
  colnames(new_count) <- unlist(lapply(1:n_cell_type, function(x){rep(x, n_cell_each[x])}))
  new_count
}
    
### copula_result contains the estimated parameter of one cell type
get_diff_path_param <- function(copula_result, ngroup = 3,
                                pUp = 0.05, pDown = 0.05, fU = 5, fL = 1.5)
{
  result <- list(copula_result)
  for(iter in 2:ngroup)
  {
    result[[iter]] <- result[[iter-1]]
    n_sel1 <- length(result[[iter]]$gene_sel1)
    # n_sel2 <- length(result[[iter]]$gene_sel2)
    
    
    d1 <- sample(c(-1, 0, 1), size = n_sel1, replace = TRUE,
                 prob = c(pDown, 1-pUp-pDown, pUp))
    # d2 <- sample(c(-1, 0, 1), size = n_sel2, replace = TRUE,
    #              prob = c(pDown, 1-pUp-pDown, pUp))
    
    f1 <- runif(n_sel1, fL, fU)
    # f2 <- runif(n_sel2, fL, fU)
    
    result[[iter]]$marginal_param1[, 3] <- result[[iter]]$marginal_param1[, 3] *
      f1 ^ d1
    # result[[iter]]$marginal_param2[, 3] <- result[[iter]]$marginal_param2[, 3] *
    #   f2 ^ d2
  }
  result
}



