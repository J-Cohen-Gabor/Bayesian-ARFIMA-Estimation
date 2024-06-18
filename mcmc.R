if (!require(nsarfima))install.packages('nsarfima',lib = "~/R_lib",repos="https://cran.csiro.au/")
if (!require(arfima))install.packages('arfima',lib = "~/R_lib",repos="https://cran.csiro.au/")
if (!require(MASS))install.packages('MASS',lib = "~/R_lib",repos="https://cran.csiro.au/")
if (!require(DataCombine))install.packages('DataCombine',lib = "~/R_lib",repos="https://cran.csiro.au/")
if (!require(FastGP))install.packages('FastGP',lib = "~/R_lib",repos="https://cran.csiro.au/")
if (!require(spam))install.packages('spam',lib = "~/R_lib",repos="https://cran.csiro.au/")
if (!require(mvtnorm))install.packages('mvtnorm',lib = "~/R_lib",repos="https://cran.csiro.au/")
if (!require(tmvtnorm))install.packages('tmvtnorm',lib = "~/R_lib",repos="https://cran.csiro.au/")
if (!require(fracdiff))install.packages('fracdiff',lib = "~/R_lib",repos="https://cran.csiro.au/")
#if (!require(sarima))install.packages('sarima',lib = "~/R_lib",repos="https://cran.csiro.au/")
#if (!require(gsl))install.packages('gsl',lib = "~/R_lib",repos="https://cran.csiro.au/")

library(fracdiff)
library(arfima)
library(nsarfima)
library(MASS)
library(FastGP)
library(extraDistr)
library(mvtnorm)
library(spam)
library(truncnorm)
#library(gsl)
library(DataCombine)
library(tmvtnorm)
library(sarima)
library(coda)

#ensure location has subfolders 'par_sims' and 'Y'
location <- 'a:'

start_time <- Sys.time()


args = commandArgs(trailingOnly=TRUE)
d <- as.double(args[1])
phi <- as.double(args[2])
theta <- as.double(args[3])

max_sims <- as.integer(args[4])

n_chain_start <- as.integer(args[5])
n_chain_end <- as.integer(args[6])

alpha <- 1
beta <- 1

d_sig <- 0.025
arma_sigm_reduction_factor <- 1

#order selection
p <- 0
q <- 0
if (phi != 0) p = 1
if (theta != 0) q = 1



true_pars <- list(d=d, mu=0, sig2=1, phi=c(), theta=c())
if (p==1){
  true_pars$phi = phi
}
if (q==1){
  true_pars$theta = theta
}
#if (READ.IN.Y){
#  Y <- read.csv('H:/2024/Honours Revisited/SNP_study/gnp_data.csv')$nomean_logchange
#  Y <- Y * 100
#  
#  #Y <- nsarfima::arfima.sim(360,d=0.3,ar=0.2)
#  
#  #Y <- Y - mean(Y)
#  n <- length(Y)
#}else{
#  n <- 1000
#}

n <- 1000

save_details <- function(loc){
  write.csv(par_sims, paste(loc,'/par_sims/par_sims_convolution_d',d,'_phi',phi,'_theta',theta,'_rep',chain,'.csv',sep=''))
  write.csv(Y, paste(loc,'/Y/Y',d,'_phi',phi,'_theta',theta,'_rep',chain,'.csv',sep=''))
  
  #png(
  #  paste(loc,'/plots/densities',d,'_phi',phi,'_theta',theta,'_rep',chain,'.png',sep=''), 
  #  width=8,height=8, units='in',res=1200,
  #  type="cairo")
  #gen_plots(par_sims[,1:3], par_lims=matrix(c(0,0.5,-1,1,-1,1),ncol=2, byrow=T))
  
  #dev.off()
  
  write.csv(par_sims_noconv, paste(loc,'/par_sims/par_sims_noconvolution_d',d,'_phi',phi,'_theta',theta,'_rep',chain,'.csv',sep=''))
  #write.csv(Y, paste(loc,'/Y/Y',d,'_phi',phi,'_theta',theta,'_rep',chain,'.csv',sep=''))
  
  #png(
  #  paste(loc,'/plots/densities',d,'_phi',phi,'_theta',theta,'_rep',chain,'.png',sep=''), 
  #  width=8,height=8, units='in',res=1200,
  #  type="cairo")
  #gen_plots(par_sims[,1:3], par_lims=matrix(c(0,0.5,-1,1,-1,1),ncol=2, byrow=T))
  
  #dev.off()
  #dev.print(png,paste(loc,'/densities.png',sep=''), width=6,height=6)
}
gen_plots <- function(mcmc, par_lims=c(), rev.order=F){
  n_pars <- ncol(mcmc)
  par(mfrow=c(n_pars,2))
  if (rev.order){
    par(mfrow=c(2,n_pars))
    for (i in 1:n_pars){
      if (length(par_lims) == 0){
        plot(as.vector(mcmc[,i]), type='l', main=paste('Trace',colnames(mcmc)[i]))
      }
      else{
        plot(as.vector(mcmc[,i]), type='l', main=paste('Trace',colnames(mcmc)[i]),
             ylim=par_lims[i,])
      }
    }
    for (i in 1:n_pars){
      if (length(par_lims) == 0){
        plot(density(mcmc[,i]), main=paste('Density',colnames(mcmc)[i]))
        abline(v=mean(mcmc[,i]))
      }
      else{
        cat(par_lims[i,],'\n')
        plot(density(as.vector(mcmc[,i])), main=paste('Density',colnames(mcmc)[i]),
             xlim=par_lims[i,])
        abline(v=mean(mcmc[,i]))
      }
    }
  }
  else{
    for (i in 1:n_pars){
      cat(i)
      if (length(par_lims) == 0){
        plot(as.vector(mcmc[,i]), type='l', main=paste('Trace',colnames(mcmc)[i]))
        plot(density(mcmc[,i]), main=paste('Density',colnames(mcmc)[i]))
        abline(v=mean(mcmc[,i]))
      }
      else{
        cat(par_lims[i,],'\n')
        plot(as.vector(mcmc[,i]), type='l', main=paste('Trace',colnames(mcmc)[i]),
             ylim=par_lims[i,])
        plot(density(as.vector(mcmc[,i])), main=paste('Density',colnames(mcmc)[i]),
             xlim=par_lims[i,])
        abline(v=mean(mcmc[,i]))
      }
    }
  }
  par(mfrow=c(1,1))
}

#source('K:/PhD/2024/Honours Revisited/arfima_library.R')
gamma_ratio_est1 <- function(x, a, b=0){
  x^a
  #(x+(a+b-1)/2)^(a-b)
}
gamma_ratio_est2 <- function(x, a, b=0){
  #x^a
  (x+(a+b-1)/2)^(a-b)
}
indic_vec <- function(n, val=1){
  matrix(rep(val, n), ncol=1, nrow=n)
}
gauss_hypergeom_2nd <- function(a,b,c,x){
  1 + x * a*b/c + x^2 * a*(a+1)*b*(b+1)/(c*(c+1)*2)
}
acf_arma.1.1 <- function(n,phi,theta){
  acf <- rep(phi,n)
  acf[1] <- (1+2*phi*theta+theta^2)/(1-phi*phi)
  acf[2] <- (phi+theta)*(1+phi*theta)/(1-phi*phi)
  acf[2:n] <- cumprod(acf[2:n])
  return(acf)
}
acf_arfima.0.d.0 <- function(n, d){
  h <- seq(0,n-1)
  ret_arr <- numeric(length(h))
  ret_arr[1] <- gamma(1-2*d)/ (gamma(1-d) * gamma(d)) * gamma(h[1] + d) / gamma(1+h[1] - d)
  ret_arr[2:n] = 
    gamma(1-2*d)/ (gamma(1-d) * gamma(d)) * 1/gamma_ratio_est2(h[-1]+d, 1-2*d)
  return(ret_arr)
}
acf_arfima.0.d.0.fast <- function(n, d){
  h <- seq(0,n-1)
  ret_arr <- numeric(length(h))
  ret_arr[1] <- gamma(1-2*d)/ (gamma(1-d) * gamma(d)) * exp(lgamma(h[-1]+d) - lgamma(1+h[-1]-d))
  ret_arr[2:n] = 
    gamma(1-2*d)/ (gamma(1-d) * gamma(d)) * 
    exp(lgamma(h[-1]+d) - lgamma(1+h[-1]-d))
  return(ret_arr)
}
acf_arfima.0.d.0.true <- function(n, d){
  h <- seq(0,n-1)
  ret_arr <- numeric(length(h))
  ret_arr[1] <- gamma(1-2*d)/ (gamma(1-d) * gamma(d)) * gamma(h[1] + d) / gamma(1+h[1] - d)
  ret_arr[2:n] = 
    gamma(1-2*d)/ (gamma(1-d) * gamma(d)) * 
    exp(lgamma(h[-1]+d) - lgamma(1+h[-1]-d))
  
  return(ret_arr)
}

acf_c_func <- function(lambda, d,h,rho,p=1){
  lambda[abs(h)+1] * 
    (rho^(2*p)*gauss_hypergeom_2nd(d+h,1,1-d+h,rho) + gauss_hypergeom_2nd(d-h,1,1-d-h,rho)-1)
}
acf_arfima.1.d.1 <- function(n, d, phi, theta){
  acf0.d.0_arr <- acf_arfima.0.d.0(n,d)
  h <- seq(0,n-1)
  c1 <- acf_c_func(acf0.d.0_arr, d,-h,-phi)
  c2 <- acf_c_func(acf0.d.0_arr, d,1-h,-phi)
  c3 <- acf_c_func(acf0.d.0_arr, d,2-h,-phi)
  return(
    (theta*c1+(1+theta^2)*c2+theta*c3)/(phi*(phi^2-1))
  )
}

#cant run gsl on artemis.... (ie hyperg_2F1)
acf_c_func.true <- function(lambda, d, h, rho, p=1, sig2=1){
  #lambda[abs(h)+1]/sig2 * 
  #  (rho^(2*p)*Re(hypergeo(d+h,1,1-d+h,rho, tol=1e-10)) + Re(hypergeo(d-h,1,1-d-h,rho, tol=1e-10))-1)
  lambda[abs(h)+1]/sig2 * 
    (rho^(2*p)*Re(hyperg_2F1(d+h,1,1-d+h,rho)) + Re(hyperg_2F1(d-h,1,1-d-h,rho))-1)
}
acf_arfima.1.d.1.true <- function(n, d, phi, theta){
  theta = -theta #palma flips notation! we use normal positive notation
  acf0.d.0_arr <- acf_arfima.0.d.0.true(n,d)
  h <- seq(0,n-1)
  c1 <- acf_c_func.true(acf0.d.0_arr, d,-h,-phi)
  c2 <- acf_c_func.true(acf0.d.0_arr, d,1-h,-phi)
  c3 <- acf_c_func.true(acf0.d.0_arr, d,2-h,-phi)
  return(
    (theta*c1+(1+theta^2)*c2+theta*c3)/(phi*(phi^2-1))
  )
}
acf_arfima.1.d.1.2ndorder <- function(n, d, phi, theta){
  theta = -theta #palma flips notation! we use normal positive notation
  acf0.d.0_arr <- acf_arfima.0.d.0.fast(n,d)
  h <- seq(0,n-1)
  c1 <- acf_c_func(acf0.d.0_arr, d,-h,-phi)
  c2 <- acf_c_func(acf0.d.0_arr, d,1-h,-phi)
  c3 <- acf_c_func(acf0.d.0_arr, d,2-h,-phi)
  return(
    (theta*c1+(1+theta^2)*c2+theta*c3)/(phi*(phi^2-1))
  )
}

post_psi0_arfima0d0_cond_y.V2 <- function(Y, d, alpha, beta){
  n <- length(Y)
  acf_vec <- acf_arfima.0.d.0.true(n,d)
  
  cov_mat <- toeplitz(acf_vec)
  #cov_mat[cov_mat < 1e-10] <- 0
  
  Sig_inv <- tinv(cov_mat)
  log_det_covmat <- determinant(cov_mat, logarithm=T)$modulus
  
  Ev <- t(Y)%*%Sig_inv%*%Y
  
  if (Ev < 0){
    return(-Inf)
  }
  return(
    -1/2 * log_det_covmat - 
      (n/2 + alpha)*log(Ev + 2*beta)
  )
}
post_psi0_cond_y.V2 <- function(Y, d, phi, theta, alpha, beta, ret.Sigma=F){
  n <- length(Y)
  #acf_vec <- acf_arfima.1.d.1.true(n,d,-phi,-theta)
  acf_vec <- tacvfARFIMA(phi=phi,theta=-theta,dfrac=d,maxlag=n-1)
  cov_mat <- toeplitz(acf_vec)
  
  Sig_inv <- tinv(cov_mat)
  log_det_covmat <- determinant(cov_mat, logarithm=T)$modulus
  
  
  Ev <- t(Y)%*%Sig_inv%*%Y
  
  if (ret.Sigma){
    if (Ev < 0){
      return(
        list(posterior=-Inf,
             cov_mat=-Inf)
      )
    }
    return(
      list(posterior=
             -1/2 * log_det_covmat - 
             (n/2 + alpha)*log(Ev + 2*beta),
           cov_mat=cov_mat
      )
    )
  }
  
  if (Ev < 0){
    return(-Inf)
  }
  return(
    -1/2 * log_det_covmat - 
      (n/2 + alpha)*log(Ev + 2*beta)
  )
}
post_psi0_arma_cond_y.V2 <- function(Y, phi, theta, alpha, beta){
  n <- length(Y)
  acf_vec <- acf_arma.1.1(n, phi,theta) #ARMAacf(ar=phi,ma=theta,lag.max=length(Y)-1)
  
  cov_mat <- toeplitz(acf_vec)
  #note: R struggles with verrrry small numbers, which occur here (e-150+)
  #round these to zero beyond roughly 8b precision
  cov_mat[abs(cov_mat) < 1e-8] <- 0
  Sig_inv <- tinv(cov_mat)
  log_det_covmat <- determinant(cov_mat, logarithm=T)$modulus
  
  Ev <- t(Y)%*%Sig_inv%*%Y
  
  if (Ev < 0){
    return(-Inf)
  }
  return(
    -1/2 * log_det_covmat - 
      (n/2 + alpha)*log(Ev + 2*beta)
  )
}
get_arfima_mle <- function(Y, p=1,q=1, d.range=c(0,0.5), incl.covmat=F){
  while(T){
    flag = F
    par0 = tryCatch(
      {
        par0 <- mle.arfima(Y, p=p, q=q, d.range=c(0,0.5), incl.mean=F)
        par0
      }, 
      error=function(e){flag=T},
      finally={}
    )
    if (!is.logical(par0)){break}
  }
  if (incl.covmat){
    return(list(pars=par0$pars, covmat=par0$cov.mat))
  }
  return(par0$pars)
}

sig2.proposal <- function(Y, alpha, beta, phi, theta, d){
  n <- length(Y)
  acf_vec <- tacvfARFIMA(phi=phi,theta=-theta,dfrac=d,maxlag=n-1)
  cov_mat <- toeplitz(acf_vec)
  
  rinvgamma(1, n/2+alpha, (t(Y)%*%tinv(cov_mat)%*%Y + 2*beta)/2)
}
sig2.proposal.covmat <- function(alpha, beta, cov_mat){
  rinvgamma(1, n/2+alpha, (t(Y)%*%tinv(cov_mat)%*%Y + 2*beta)/2)
}

#set.seed(1)

#YDF <- readRDS(
#  paste('K:/PhD/2024/Honours Revisited/SIMULATIONS/generated_Y/phi_eq_theta_eq_0.5_sig1_mu0/',
#        'd_',d,'.rda',sep='')
#)

start_time <- Sys.time()

for (chain in n_chain_start:n_chain_end){
  set.seed(chain)
  Y <- nsarfima::arfima.sim(n=n, d=true_pars$d, ar=true_pars$phi, ma=true_pars$theta,
                            mu=true_pars$mu, sig2=true_pars$sig2)

  #fix_Y_idx <- 1
  #fix_Y_idx <- 3
  #Y <- YDF$Y[(1+(fix_Y_idx-1)*n):(fix_Y_idx*n)]
  
  par_sims <- matrix(rep(0,(length(true_pars)+6)*max_sims), nrow=max_sims,
                     dimnames=list(rep(0,max_sims),c( 
                       "d","phi","theta",
                       "prop_d","prop_phi","prop_theta",
                       "posterior_d","alpha_d","posterior_arma","alpha_arma",
                       "sig2")))
  
  par_sims_noconv <- matrix(rep(0,(length(true_pars)+4)*max_sims), nrow=max_sims,
                            dimnames=list(rep(0,max_sims),c( 
                              "d","phi","theta",
                              "prop_d","prop_phi","prop_theta",
                              "posterior","alpha",
                              "sig2")))
  
  psi0_pars <- c('d')
  psi0_pars2 <- c('d')
  if (p == 1){
    psi0_pars <- c(psi0_pars, 'ar.1')
    psi0_pars2 <- c(psi0_pars2, 'phi')
  }
  if (q == 1){
    psi0_pars <- c(psi0_pars, 'ma.1')
    psi0_pars2 <- c(psi0_pars2, 'theta')
  }
  
  par_count <- 0
  par0_arfima <- get_arfima_mle(Y,p=p,q=q, incl.covmat=T)
  par0 <- par0_arfima$pars
  par0_covmat <- par0_arfima$covmat[psi0_pars,psi0_pars] / arma_sigm_reduction_factor
  
  sig2.mle <- par0['sig2']
  
  mle_arfima1d1 <- (par0[psi0_pars])
  names(mle_arfima1d1) <- psi0_pars2
  if (is.na(mle_arfima1d1['phi'])) mle_arfima1d1['phi'] <- 0
  if (is.na(mle_arfima1d1['theta'])) mle_arfima1d1['theta'] <- 0

  arma_pars <- c()
  arma_pars2 <- c()
  if (p==1) {
    arma_pars <- c(arma_pars, 'ar1')
    arma_pars2 <- c(arma_pars2, 'phi')
  }
  if (q==1){
    arma_pars <- c(arma_pars, 'ma1') 
    arma_pars2 <- c(arma_pars2, 'theta')
  }
  if (p > 0 | q > 0){
    #two chain analysis - start with the arma chain
    U <- diffseries(Y, d=mle_arfima1d1['d'])[-1] #diffseries just includes the first point
    arma_fit <- arima(U, order=c(p,0,q))
  
    arma_mle <- (arma_fit$coef[arma_pars])
    names(arma_mle) <- c(arma_pars2)
    
    arma_sigm <- arma_fit$var.coef[arma_pars,arma_pars]/arma_sigm_reduction_factor
    
    par_sims[1, arma_pars2] <- rtmvnorm(1, mean=unlist(arma_mle), sigma=arma_sigm,
                                             lower=rep(-1,p+q),upper=rep(1,p+q))
    par_sims[1, 'posterior_arma'] <- 
      post_psi0_arma_cond_y.V2(U, par_sims[1,'phi'], par_sims[1,'theta'], alpha, beta)
    prev_accepted_post_arma <- par_sims[1, 'posterior_arma']
    
    num_accept_arma <- 1
    
    #setup the d chain
    Z <- xarmaFilter(x=Y, model=list(ar=mle_arfima1d1['phi'], ma=mle_arfima1d1['theta']), whiten=T)
    d_mle <- list(d=mle_arfima1d1['d'])
    
    #run 1 had 0.03
    #d_sig2 <- 0.03 #tweak parameter?
    
    
    par_sims[1,'d'] <- rtruncnorm(1, mean=d_mle, sd=d_sig, a=-0.5, b=0.5)
    par_sims[1, 'posterior_d'] <- 
      post_psi0_arfima0d0_cond_y.V2(Z, par_sims[1,'d'], alpha, beta)
    prev_accepted_post_d <- par_sims[1, 'posterior_d']
    
    num_accept_d <- 1
    prev_num_accept_d <- 0
  }
  ###########################
  #setup for non-conv version
  first_val <- rtmvnorm(1, mean=par0[psi0_pars],sigma=par0_covmat,
                        lower=c(-0.5,rep(-1,p+q)), upper=c(0.5,rep(1,p+q)))
  par_sims_noconv[1,psi0_pars2] <- first_val
  par_sims_noconv[1,'posterior'] <- 
    post_psi0_cond_y.V2(Y, par_sims_noconv[1,'d'],par_sims_noconv[1,'phi'],par_sims_noconv[1,'theta'],
                        alpha, beta)
  
  par_sims[1,'sig2'] <- sig2.proposal(Y, alpha, beta, 
                                      par_sims[1,'phi'], par_sims[1,'theta'], par_sims[1,'d'])
  par_sims_noconv[1,'sig2'] <- sig2.proposal(Y, alpha, beta, 
                                             par_sims_noconv[1,'phi'], par_sims_noconv[1,'theta'], par_sims_noconv[1,'d'])
  
  num_accept_arfima <- 1
  prev_num_accept_arfima <- 0
  prev_prop_val_arfima <- 999
  
  curr_sim <- 2
  while(curr_sim <= max_sims){
    if (curr_sim %% 10 == 0){
      cat('Sim =',curr_sim, '-- ARFIMA Chain - Num accept (d,ar,ma) =',num_accept_arfima, ' Filtered ARFIMA Chain - Num accept d =',
          num_accept_d, 'Num accept arma =',num_accept_arma, '   \r')
    }    
    
    #proposal for non-filter
    if (T){
      prop_psi0_mu <- par_sims_noconv[curr_sim-1,psi0_pars2]  
      #c(par_sims_noconv[curr_sim-1,'d'], par_sims_noconv[curr_sim-1,'phi'], par_sims_noconv[curr_sim-1,'theta'])
      names(prop_psi0_mu) <- psi0_pars2
      prop_psi0 <- as.vector(rtmvnorm(1, mean=prop_psi0_mu, sigma=par0_covmat,
                                      lower=c(-0.5,rep(-1,p+q)),upper=c(0.5,rep(1,p+q))))
      names(prop_psi0) <- psi0_pars2
      if (is.na(prop_psi0['phi'])) prop_psi0['phi'] <- 0
      if (is.na(prop_psi0['theta'])) prop_psi0['theta'] <- 0
      
      post_ret <- post_psi0_cond_y.V2(Y, prop_psi0['d'], prop_psi0['phi'], prop_psi0['theta'], alpha, beta, ret.Sigma = T)
      post_psi_prop <- post_ret$posterior
      psi0_covmat <- post_ret$cov_mat
      post_psi_prop2 <- prev_prop_val_arfima
      if (num_accept_arfima != prev_num_accept_arfima){
        post_ret2 <- post_psi0_cond_y.V2(Y, par_sims_noconv[curr_sim-1,'d'], par_sims_noconv[curr_sim-1,'phi'], 
                                         par_sims_noconv[curr_sim-1,'theta'], alpha, beta, ret.Sigma=T)
        post_psi_prop2 <- post_ret2$posterior
        prev_num_accept_arfima <- prev_num_accept_arfima + 1
        prev_prop_val_arfima <- post_psi_prop2
      }
      
      alpha_psi0 <- 0
      dens_prop1 <- 0
      dens_prop2 <- 0
      if (!is.na(post_psi_prop) && !is.infinite(post_psi_prop)){
        #calc proposal densities
        dens_prop1 <- dtmvnorm(par_sims_noconv[curr_sim-1, psi0_pars2],mean=as.vector(prop_psi0[psi0_pars2]), sigma=par0_covmat,
                               lower=c(-0.5,rep(-1,p+q)),upper=c(0.5,rep(1,p+q)), log=T)
        dens_prop2 <- dtmvnorm(as.vector(prop_psi0[psi0_pars2]),mean=par_sims_noconv[curr_sim-1, psi0_pars2], sigma=par0_covmat,
                               lower=c(-0.5,rep(-1,p+q)),upper=c(0.5,rep(1,p+q)), log=T)
        
        alpha_s <- post_psi_prop + dens_prop1 - post_psi_prop2 - dens_prop2
        alpha_psi0 <- min(1, exp(alpha_s))
      }
      w <- rbinom(1, 1, alpha_psi0)
      if (w == 1){
        par_sims_noconv[curr_sim, psi0_pars2] <- prop_psi0[psi0_pars2]
        num_accept_arfima <- num_accept_arfima + 1
        
        prop_psi0_mu <- par_sims_noconv[curr_sim, psi0_pars2]
      }
      else{
        psi0_covmat <- post_ret2$cov_mat
        par_sims_noconv[curr_sim, psi0_pars2] <- par_sims_noconv[curr_sim-1, psi0_pars2]
      }
      par_sims_noconv[curr_sim, 'prop_d'] <- prop_psi0['d']
      par_sims_noconv[curr_sim, 'prop_phi'] <- prop_psi0['phi']
      par_sims_noconv[curr_sim, "prop_theta"] <- prop_psi0['theta']
      par_sims_noconv[curr_sim, 'posterior'] <- post_psi_prop
      par_sims_noconv[curr_sim, 'alpha'] <- alpha_psi0
      par_sims_noconv[curr_sim, 'sig2'] <- 
        sig2.proposal.covmat(alpha, beta, cov_mat=psi0_covmat)
    }
    if (T & (p > 0 | q > 0)){
      prop_arma_mu <- par_sims[curr_sim-1,arma_pars2]
      names(prop_arma_mu) <- arma_pars2
      prop_arma <- rtmvnorm(1, mean=unlist(prop_arma_mu), sigma=arma_sigm,
                            lower=rep(-1,p+q),upper=rep(1,p+q))
      names(prop_arma) <- arma_pars2
      prop_d_mu <- par_sims[curr_sim-1,'d']
      prop_d <- rtruncnorm(1, mean=prop_d_mu, sd=d_sig, a=-0.5,b=0.5)
      
      #######################################
      #NON-CONVOLUTION FILTERED CODE FINISHED 
      #CONVOLUTION COMPONENT CALCULATIONS HERE
      #generate z~arfima(0,d,0) from the last accepted (phi,theta)
      Z <- xarmaFilter(x=Y, model=list(ar=par_sims[curr_sim-1,'phi'],ma=par_sims[curr_sim-1,'theta']), whiten=T)
      
      #determine the acceptance of the long-memory parameter d
      post_d_prop <- post_psi0_arfima0d0_cond_y.V2(Z, prop_d, alpha, beta)
      post_d_prop2 <- post_psi0_arfima0d0_cond_y.V2(Z, par_sims[curr_sim-1,'d'], alpha, beta)
      
      alpha_d <- 0
      if (!is.na(post_d_prop) && !is.infinite(post_d_prop)){
        dens_prop1 <- dtruncnorm(prop_d, mean=prop_d_mu, sd=d_sig, a=-0.5,b=0.5)
        dens_prop2 <- dtruncnorm(prop_d_mu, mean=prop_d, sd=d_sig, a=-0.5,b=0.5)
        
        alpha_s <- post_d_prop + dens_prop1 - post_d_prop2 - dens_prop2
        alpha_d <- min(1, exp(alpha_s))
      }
      
      w_d <- rbinom(1,1,alpha_d)
      if (w_d == 1){
        par_sims[curr_sim, 'd'] <- prop_d
        prev_accepted_post_d <- post_d_prop
        num_accept_d <- num_accept_d + 1
      }
      else{
        par_sims[curr_sim, 'd'] <- par_sims[curr_sim-1, 'd']
      }
      
      #generate u~arma(1,1) from the last accepted d
      U <- diffseries(Y, d=par_sims[curr_sim,'d'])[-1] #diffseries just includes the first point
      
      if (is.na(prop_arma['phi'])) prop_arma['phi'] <- 0
      if (is.na(prop_arma['theta'])) prop_arma['theta'] <- 0
      
      post_arma_prop <- post_psi0_arma_cond_y.V2(U, prop_arma['phi'],prop_arma['theta'], alpha, beta)
      post_arma_prop2 <- post_psi0_arma_cond_y.V2(U, par_sims[curr_sim-1,'phi'],par_sims[curr_sim-1,'theta'], alpha, beta)
      
      #determine the acceptance of the ARMA parameters (phi,theta)
      alpha_arma <- 0
      if (!is.na(post_arma_prop) && !is.infinite(post_arma_prop)){
        dens_prop1 <- dtmvnorm(prop_arma[arma_pars2], mean=prop_arma_mu[arma_pars2], sigma=arma_sigm,
                               lower=rep(-1,p+q),upper=rep(1,p+q))
        
        dens_prop2 <- dtmvnorm(prop_arma_mu[arma_pars2], mean=as.vector(prop_arma[arma_pars2]), sigma=arma_sigm,
                               lower=rep(-1,p+q),upper=rep(1,p+q))
        alpha_s <- post_arma_prop + dens_prop1 - post_arma_prop2 - dens_prop2
        alpha_arma <- min(1, exp(alpha_s))
      }
      w_arma <- rbinom(1,1,alpha_arma)
      if (w_arma == 1){
        par_sims[curr_sim, arma_pars2] <- prop_arma[arma_pars2]
        prev_accepted_post_arma <- post_arma_prop
        num_accept_arma <- num_accept_arma + 1
      }
      else{
        par_sims[curr_sim, arma_pars2] <- par_sims[curr_sim-1, arma_pars2]
      }
      
      par_sims[curr_sim, 'prop_d'] <- prop_d
      par_sims[curr_sim, 'posterior_d'] <- post_d_prop
      par_sims[curr_sim, 'alpha_d'] <- alpha_d
      
      par_sims[curr_sim, c("prop_phi","prop_theta")] <- prop_arma
      par_sims[curr_sim, 'posterior_arma'] <- post_arma_prop
      par_sims[curr_sim, 'alpha_arma'] <- alpha_arma
      
      par_sims[curr_sim, "sig2"] <- sig2.proposal(Y, 
                                                  alpha, beta, par_sims[curr_sim,'phi'],
                                                  par_sims[curr_sim,'theta'], par_sims[curr_sim,'d'])
    }
    
    curr_sim <- curr_sim + 1
  }
  save_details(location)
}

end_time <- Sys.time()
end_time-start_time


cat('\nDone\n')

#gen_plots(par_sims[,1:3], par_lims=matrix(c(0,0.5,-1,1,-1,1),ncol=2, byrow=T))

