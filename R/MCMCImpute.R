#' Impute single cell HiC data
#' 
#' This function imputes the single cell data by borrowing information from beighborbood, similar cells, and bulk data.
#' 
#' @param niter Number of MCMC iteration.
#' @param burnin Number of MCMC burnin.
#' @param single Single cell matrix with each column being the upper triangular matrix of a single cell.
#' @param bulk Bulk data. If NULL, bulk is set to be the sum of all single cells.
#' @param startval starting value of MCMC chain.
#' @param n Dimension of single cell matrix(original matrix is square and symmetrical).
#' @param epsilon1 Range size of delta, default is 0.5.
#' @param epsilon2 Range size of B, default is 5.
#' @param mc.cores Number of cores to be used in mclapply function, default is 1.
#' @param cutoff The threshold of \eqn{\pi_{ij}} that is used to define structural zeros.
#' 
#' @import parallel
#' @import Rcpp
#' @import RcppArmadillo
#' 
#' @return A list of posterior mean of probability, imputed data without defining SZ, and imputed data with SZ.
#' @export
#'
#' @examples
#' data("K562_T1_7k")
#' data("K562_bulk")
#' single=K562_T1_7k
#' T1_7k_res=MCMCImpute(niter=100000,burnin=5000,single=K562_T1_7k,bulk=K562_bulk,
#' startval=c(100,100,10,8,10,0.1,900,0.2,0,replicate(dim(single)[2],8)),n=61,mc.cores = 1,cutoff=0.5)

MCMCImpute <- function(niter=30000,burnin=15000,single,bulk=bulk,startval=c(100,100,10,8,10,0.1,900,0.2,0,replicate(dim(single)[2],8)),n,epsilon1=0.5,epsilon2=5,mc.cores = 1,cutoff=0.5){
  ## normalize single cells to get same depth as max
  single_sum <- apply(single, 2, sum)
  max_single <- max(single_sum)
  lambda <- single_sum/max_single
  single_norm <- t(t(single)/lambda)
  
  ## variance for each contact pair
  mu_sigma = neivar(single_norm, nei=5, n)
  mu_sigma_vec = mu_sigma[upper.tri(mu_sigma)]
  
  ## single weighted sum 
  n_single = dim(single)[2]
  weight <- single_sum/sum(single_sum)
  type_bulk_vec <- apply(t(t(single_norm)*weight),1,sum)
  type_bulk <- matrix(0, n, n) 
  type_bulk[upper.tri(type_bulk, diag=FALSE)] <- type_bulk_vec
  
  ## correction matrix
  correct = correctfac(type_bulk,nei = 5)
  if(is.null(bulk)){
    bulk =  apply(single, 1, sum)
  }
  
  ## combine information into one matrix
  m <- single
  rowmax <- apply(m, 1, max) 
  prop0 <- rowMeans(m==0, na.rm = FALSE, dims = 1)
  correctvec <- correct[upper.tri(correct,diag = FALSE)]
  B <- bulk*correctvec*rowmax/sum(bulk) 
  m <- as.matrix(cbind(m,rowmax,prop0,correctvec,B))
  
  ## apply oneimpute function to each row of matrix
  result=mclapply(0:(nrow(single)-1), oneimpute, niter=niter, burnin=burnin, n_single=n_single, m=as.matrix(m), lambda=lambda, mu_sigma_vec=mu_sigma_vec, startval=startval, n=n, epsilon1=epsilon1, epsilon2=epsilon2, mc.cores = mc.cores)
  result=do.call(rbind,result)
  
  ## organize MCMC result
  pii_mean <- result[,1]
  pii <- matrix(0, n, n)
  pii[upper.tri(pii, diag=FALSE)] <- pii_mean
  post_mu_k <- result[,-1]
  
  ## true 0 positions figured out by MCMC
  # posi <- which(pii>=cutoff,TRUE)
  # posirows <- NULL
  # for (i in 1:(dim(posi)[1])) {
  #   posirows <- c(posirows,matrow(posi[i,1],posi[i,2]))
  # }
  
  ## Imputation
  IMP1 <- NULL
  for (i in 1:ncol(single)) {
    imp <- post_mu_k[,i]*lambda[i]
    IMP1 <- cbind(IMP1,imp)
  }
  
  IMP2=IMP1
  IMP2[pii_mean>0.5,]<-0
  
  ## output
  output=list(pii,IMP1, IMP2)
  names(output) = c("pii","IMP1","IMP2")
  return(output)
}



