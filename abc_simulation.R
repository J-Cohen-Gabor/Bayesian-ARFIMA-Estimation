#install.packages("nsarfima")
#library(nsarfima,lib.loc="~/R_libs")
library(tseries)
library(truncnorm)
library(nsarfima)
library(fracdiff)

location <- 'A:'
start_time <- Sys.time()

#arg[1] must be the index of the d value to run from 1-101 (running 100 d's from 0-1)
args = commandArgs(trailingOnly=TRUE)

keep.quantile <- c(0.001,0.005,0.01)
n <- 1000

d <- as.double(args[1])
phi.real <- as.double(args[2])
theta.real <- as.double(args[3])

p <- 0
q <- 0

if (phi.real != 0){
  p <- 1
}
if (theta.real != 0){
  q <- 1
}

n <- 1000

sig2.quant <- 0.5

m <- as.integer(args[4])

n_rep_start <- as.integer(args[5])
n_rep_end <- as.integer(args[6])

calc_logperiog <- F
J <- 2 #as.integer(args[7])
l <- 0 #as.integer(args[8])

#n_rep <- 25
num_periog <- 20

#d.raw <- 5

cat('running ',d,'\n')

get_logperiog_regressor <- function(Y.periog, J=5, l=1){
  m <- length(Y.periog)
  
  K <- seq(l+J, m, by=J)
  Ykj <- sapply(K, function(k){
    kj_seq <- k+J-1:J
    log(sum(Y.periog[kj_seq]))
  })
  return(Ykj[-length(Ykj)])
}
get_arima_mle <- function(Y, p=1,q=1){
  par_count <- 1
  while(T){
    if (par_count > 20){
      par_count <- 0
      #cat('\nPar count error\n')
      return(NULL)
    }
    flag = F
    par0 = tryCatch(
      {
        par0 <- arima(Y, order=c(p,0,q), include.mean=F)
        par0
      }, 
      error=function(e){flag=T},
      finally={}
    )
    if (!is.logical(par0)){break}
    par_count <- par_count + 1
  }
  return(par0$coef)
}
get_arfima_mle <- function(Y, p=1,q=1, d.range=c(0,0.5)){
  while(T){
    #if (par_count > 20){
    #  par_count <- 0
    #  #cat('new Y\n')
    #  Y <- nsarfima::arfima.sim(n=n, d=true_pars$d, ar=true_pars$theta, ma=true_pars$phi,
    #                            mu=true_pars$mu, sig2=true_pars$sig2)
    #}
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
  return(par0$pars)
}

get_abc_densities <- function(quant_df, quant, log.periog=T, periog20=F, true_vals=NULL){
  quant <- as.character(quant)
  type <- 's.logperiog'
  if (log.periog==F){
    type <- 's.periog'
  }
  if (periog20 == T){
    type <- paste0(type,'20')
  }
  else{
    type <- paste0(type,'full')
  }
  par(mfrow=c(2,2))
  plot(density(quant_df[[quant]][[type]]$d), main='d (red=true value, grey=mean)', xlim=c(0,0.5))
  if (!is.null(true_vals))abline(v=true_vals['d'], col='red')
  abline(v=mean(quant_df[[quant]][[type]]$d), col='grey')
  
  plot(density(quant_df[[quant]][[type]]$sig2), main='sig2', xlim=c(0,2))
  if (!is.null(true_vals))abline(v=true_vals['sig2'], col='red')
  abline(v=mean(quant_df[[quant]][[type]]$sig2), col='grey')
  
  plot(density(quant_df[[quant]][[type]]$ar), main='phi',xlim=c(-1,1))
  if (!is.null(true_vals))abline(v=true_vals['phi'], col='red')
  abline(v=mean(quant_df[[quant]][[type]]$ar), col='grey')
  
  plot(density(quant_df[[quant]][[type]]$ma), main='theta',xlim=c(-1,1))
  if (!is.null(true_vals))abline(v=true_vals['theta'], col='red')
  abline(v=mean(quant_df[[quant]][[type]]$ma), col='grey')
  par(mfrow=c(1,1))
}

#set.seed(d.raw)

q_pos <- as.list(round(quantile(seq(1,m),keep.quantile)))
names(q_pos) <- keep.quantile
#write.table(q_pos,file=paste('q_pos_',args[1],'.csv',sep=''), row.names=FALSE)

d_it <- 1
idx <- 1

cat('Allocating measure matrix...\n')
stt <- paste0('s.logperiog.j',J,'.l',l)

ddf <- data.frame(d.sim=numeric(m),
                  phi.sim=numeric(m),
                  theta.sim=numeric(m),
                  sig2.sim=numeric(m),
                  s.sig2=numeric(m),
                  s.arma=numeric(m),
                  s.periog20=numeric(m),
                  s.periogfull=numeric(m),
                  pll=numeric(m)
)
names(ddf)[names(ddf) == 'pll'] <- stt
ptypes <- names(ddf)[-c(1,2,3,4,5,6)]
Y.sim.measure.mat <- as.matrix(ddf)

rm(ddf)
cat('Done allocating\n')
gc()

#+2 are the following two arrays: is_zero [] and d.hat
#Y.sim.measure.is_zero <- rep(0,m)
#Y.sim.measure.d.hat <- rep(0,m)

#set.seed(d.raw*(n_rep_start + n_rep_end))
#bin_in <- file('/project/mcmc_arfima_sim1/ABC_datasets/sim_reps10m_n1000_combined.bin','rb')
#seek(bin_in, 100*m*(n_rep_start-1))

#d.sim <- rep(0,m)

accept_quant_df <- list()
accept_quant_df <- lapply(keep.quantile, function(q)q=c())
names(accept_quant_df) <- keep.quantile

for (rep in n_rep_start:n_rep_end){
  #seek(bin_in, 0)
  
  set.seed(rep)
  Y.real <- nsarfima::arfima.sim(n, ar=phi.real, ma=theta.real, d=d)
  
  Y.fft <- (1/n * Mod(fft(Y.real))^2)[1:(n/2)]
  Y.var <- var(Y.real)
  
  if (calc_logperiog){
    lp.j.l <- get_logperiog_regressor(Y.fft, J=J, l=l)
    lp.j.l <- lp.j.l[!is.na(lp.j.l)]
  }
  
  Y.mle.full <- get_arfima_mle(Y.real,p,q)
  
  if (is.null(Y.mle.full)){
    rep <- rep - 1
    next
  }
  if (p == 1){
    Y.mle <- Y.mle.full['ar.1'] 
  }
  else{
    Y.mle <- Y.mle.full['ma.1']
  }
  
  i <- 1

  while(i <= m){
    if (i%%1000==0)cat('M: ',i,'\r')
    phi.sim <- 0
    theta.sim <- 0
    
    if (p ==1){
      phi.sim <- runif(1,-1,1)
      arma.sim <- phi.sim
    }
    if (q == 1){
      theta.sim <- runif(1,-1,1)
      arma.sim <- phi.sim
    }
    
    d.sim <- runif(1,-0.5,0.5)
    sig2.sim <- rinvgamma(1,28,30)
    
    Y.sim <- nsarfima::arfima.sim(n=n,d=d.sim,ar=phi.sim,ma=theta.sim, sig2=sig2.sim)
    Y.sim.var <- var(Y.sim)
    U <- diffseries(Y.sim, d=d.sim) #diffseries just includes the first point
    
    Y.sim.fft <- (1/n * Mod(fft(Y.sim))^2)[1:(n/2)]
    Y.sim.periodogram <- Y.sim.fft
    
    if (calc_logperiog){
      Y.sim.lp.j.l <- get_logperiog_regressor(Y.sim.fft, J=2, l=0)
      Y.sim.lp.j.l <- Y.sim.lp.j.l[!is.na(Y.sim.lp.j.l)]
      #trailing nas if a non divisible n is used
    }
    
    Y.sim.mle_arma <- get_arima_mle(U, p=p, q=q)
    
    if (is.null(Y.sim.mle_arma)){
      next
    }
    
    if (p==1){
      Y.sim.mle <- Y.sim.mle_arma['ar1']
    }
    if (q == 1){
      Y.sim.mle <- Y.sim.mle_arma['ma1']
    }
    
    Y.sim.measure.mat[i,'d.sim'] <- d.sim
    if (p == 1){
      Y.sim.measure.mat[i,'phi.sim'] <- phi.sim #arma.sim
    }
    if (q == 1){
      Y.sim.measure.mat[i,'theta.sim'] <- theta.sim #arma.sim
    }
    Y.sim.measure.mat[i,'sig2.sim'] <- sig2.sim
    
    #sim.lp.j.l <- get_logperiog_regressor(Y.sim.fft, J=J, l=l)
    if (calc_logperiog){
      Y.sim.measure.mat[i,stt] <- sqrt(sum((lp.j.l - Y.sim.lp.j.l)^2))
    }
    
    Y.sim.measure.mat[i,'s.periog20'] <- sqrt(sum((Y.fft[1:num_periog] - Y.sim.fft[1:num_periog])^2))
    Y.sim.measure.mat[i,'s.periogfull'] <- sqrt(sum((Y.fft - Y.sim.fft)^2))
    
    if (p == 1 & q == 1){
      Y.sim.measure.mat[i,'s.arma'] <- sqrt((Y.sim.mle_arma['ar1'] - Y.mle.full['ar.1'])^2 + 
                                              (Y.sim.mle['ma1'] - Y.mle.full['ma.1'])^2)
    }
    else{
      Y.sim.measure.mat[i,'s.arma'] <- sqrt((Y.sim.mle - Y.mle)^2)
    }
    
    Y.sim.measure.mat[i,'s.sig2'] <- abs(Y.var - Y.sim.var)
    
    i <- i + 1
  }
  cat('\rdone m loop\n')
  #loop the 3 quantiles
  
  for (quant in keep.quantile){
    quant <- as.character(quant)
    #loop the number of periogram thresholds we are keeping
    accept_quant_df[[quant]] <- list()
    for (p_type in ptypes){
      cat('Quant: ',quant, ', Ptype: ', p_type,'             \n')
      
      d_threshold <- (Y.sim.measure.mat[order(Y.sim.measure.mat[,p_type]),p_type])[q_pos[[quant]]]
      arma_threshold <- (Y.sim.measure.mat[order(Y.sim.measure.mat[,'s.arma']),'s.arma'])[q_pos[[quant]]]
      
      sig2_threshold <- (Y.sim.measure.mat[order(Y.sim.measure.mat[,'s.sig2']),'s.sig2'])[floor(m*sig2.quant)]
      
      accept_mask <- Y.sim.measure.mat[,p_type] < d_threshold & Y.sim.measure.mat[,'s.arma'] < arma_threshold &
        Y.sim.measure.mat[,'s.sig2'] < sig2_threshold
      
      
      if (sum(accept_mask) == 0){
        cat("ERROR NOT ENOUGH VALUES TO ACCEPT")
        next
      }      
      
      Y.sim.accept <- Y.sim.measure.mat[accept_mask,c('d.sim','phi.sim','theta.sim','sig2.sim'), drop=F]
      names(Y.sim.accept) <- c('d','phi','theta','sig2')
      
      #Y.sim.accept.arma <- Y.sim.measure.mat[accept_mask,c('phi.sim','theta.sim'), drop=F]
      #Y.sim.accept.sig2 <- Y.sim.measure.mat[accept_mask, 'sig2.sim', drop=F]
      #Y.sim.phi <- Y.sim.accept.arma[,'phi.sim']
      #Y.sim.theta <- Y.sim.accept.arma[,'theta.sim']
      
      #df <- data.frame(d=Y.sim.accept[,'d.sim'],
      #                 ar=Y.sim.phi,
      #                 ma=Y.sim.theta,
      #                 sig2=Y.sim.accept.sig2[,1])
      
      file_loc <- paste(location, '/d',
                        d,'_phi',phi.real,'_theta',theta.real,
                        '_quant',quant,
                        '_sum',p_type,
                        '_rep',rep,'.csv',sep='')
      #have to write out in chunks due to memory limitations in hpc
      chunk_size <- 100
      nrows <- nrow(Y.sim.accept)
      if (nrows == 0){
        next
      }
      chunk_arr <- seq(0,nrows, by=chunk_size)
      if (nrows%%chunk_size != 0){
        chunk_arr <- c(chunk_arr,nrows)
      }
      
      for (i in 1:(length(chunk_arr)-1)){
        start_idx <- chunk_arr[i]+1
        end_idx <- chunk_arr[i+1]
        if (i == 1){
          write.table(Y.sim.accept[start_idx:end_idx,, drop=F], 
                      file_loc, sep=',', col.names = T, row.names = F)  
        }
        else{
          write.table(Y.sim.accept[start_idx:end_idx,, drop=F], 
                      file_loc, append=T, sep=',', col.names = F, row.names=F)  
        }
      }
      
      idx <- idx + 1
      #cat('done of periog loop\n')
    }
  }
  cat('Done writing\n')
}

end_time <- Sys.time()
end_time-start_time


#warnings()
