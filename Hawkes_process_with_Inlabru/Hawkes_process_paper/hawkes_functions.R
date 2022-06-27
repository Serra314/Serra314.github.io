library(ggplot2)
library(inlabru)
library(INLA)
library(foreach)
library(bayesianETAS)
library(dplyr)
library(ggstar)
library(foreach)

# time triggering function - used by bayesianETAS
gt <- function(th, t, ti, mi, M0){
  # set 0 
  output <- rep(0,length(ti))
  before <- ti < t
  # update only ti < t
  if(sum(before) > 0){
    log.out <- log(th[2]) + th[3]*(mi[before] - M0)  + log(th[5] - 1) + (th[5] - 1)*log(th[4]) - th[5]*log(t - ti[before] + th[4])
    output[before] <- exp(log.out)  
  }
  else{
    output
  }
  output
}

# time triggering function - used by Inlabru
gt.2 <- function(th, t, ti, mi, M0){
  output <- rep(0,length(ti))
  t.diff <- t - ti
  neg <- t.diff <= 0
  if(sum(!neg) > 0){
    log.out <- log(th[2]) + th[3]*(mi[!neg] - M0)  - th[5]*log(1 + t.diff[!neg]/th[4])
    output[!neg] <- exp(log.out)  
  }
  else{
    output
  }
  output
}

# conditional intensity - used by bayesianETAS
lambda_ <- function(th, t, ti.v, mi.v, M0){
  if(is.null(ti.v) | all(ti.v > t)){
    th[1]
  }
  th[1] + sum(gt(th, t, ti.v, mi.v, M0))
}

# conditional intensity - used by Inlabru
lambda_2 <- function(th, t, ti.v, mi.v, M0){
  if(is.null(ti.v) | all(ti.v > t)){
    th[1]
  }
  th[1] + sum(gt.2(th, t, ti.v, mi.v, M0))
}

# integrated triggering function - used by bayesianETAS
log.Lambda_h <- function(th, ti, mi, M0, T1, T2){
  th <- as.numeric(th)
  T.low <- sapply(ti, \(x) max(T1, x))
  log(th[2]) + th[3]*(mi - M0) + log((th[4]/(T.low - ti + th[4]))^(th[5] - 1) - (th[4]/(T2 - ti + th[4]))^(th[5] - 1))
}

# integrated triggering function - used by Inlabru
log.Lambda_h2 <- function(th, ti, mi, M0, T1, T2){
  th <- as.numeric(th)
  T.low <- sapply(ti, \(x) max(T1, x))
  
  gamma.l <- (T.low - ti)/th[4]
  gamma.u <- (T2 - ti)/th[4]
  w.l <- (1 + gamma.l)^(1-th[5])  
  w.u <- (1 + gamma.u)^(1-th[5])
  # output
  log(th[2]) + th[3]*(mi - M0) + log(th[4]) - log(th[5] - 1) + log1p(w.l - 1) + log1p(-w.u/w.l) 
  
}

# find breaks point for grid
breaks_exp <- function(tt_, T2_, coef_ = 2, delta_, N_exp_ = 10){
  
  tt_breaks <- tt_ + delta_*((1 + coef_)^(0:N_exp_))
  tt_breaks <- tt_breaks[tt_breaks < T2]
  if(T2_ - tt_ < delta_){
    return(c(tt_, T2_))
  }
  if(T2 - tt_breaks[length(tt_breaks)] < delta_){
    tt_breaks[length(tt_breaks)] = T2_
  }
  if(tt_breaks[length(tt_breaks)] < T2_){
    tt_breaks <- c(tt_breaks, T2_)
  }
  return(c(tt_,tt_breaks))
} 


time.grid <- function(data.point, coef.t, delta.t,
                            T2., displaygrid = FALSE, N.exp.){
  tt. <- data.point$ts
  idx.p <- data.point$idx.p
  # spatial bins
  if(displaygrid){
    Plot_grid(xx = xx., yy = yy., delta_ = delta., n.layer = n.layer., 
              bdy_ =  bdy., min.edge = min.edge.)
  }

  # time bins
  # find bins break points
  t_b <- breaks_exp(tt., T2., coef_ = coef.t, delta_ = delta.t, N_exp_ = N.exp.)
  time.bins <- data.frame(t.start = t_b[-length(t_b)], 
                          t.end = t_b[-1]) %>%
    mutate(t.bin.name = paste0(round(t.start,3),'-',round(t.end,3)))
  if(nrow(time.bins) - 1 == 0){
    time.bins$t.ref_layer = paste0('last-',idx.p)  
  }
  else{
    time.bins$t.ref_layer = c(1:(nrow(time.bins) - 1), paste0('last-',idx.p))
  } 
  cbind(time.bins, data.point)
}

It_df <- function(param_, time.df){
  tth <- as.numeric(time.df$ts)
  T1b <- as.numeric(time.df$t.start)
  T2b <- as.numeric(time.df$t.end)
  param_c <- param_[4]
  param_p <- param_[5]
  T.l <- sapply(1:length(tth), \(x) max(tth[x], T1b[x]))
  fun.l <- (1 + (T.l - tth)/param_c)^(1-param_p)
  fun.u <- (1 + (T2b - tth)/param_c)^(1-param_p)
  ( param_c/ (param_p - 1) )* ( fun.l - fun.u )
}


compute.grid <- function(param., gridd){
  
  t.names <- unique(gridd$t.ref_layer)
  time.sel <- gridd[sapply(t.names, \(bname) which(gridd$t.ref_layer == bname)[1]),]

  It.vec <- It_df(param_ = param., time.df = time.sel)
  
  list.ret <- list(It = It.vec[match(gridd$t.ref_layer, t.names)],
                   time.df = time.sel %>%
                     mutate(It = It.vec))
  return(list.ret)
}



# function to fit Hawkes process model
Hawkes.bru <- function(sample.s, M0, T1, T2, link.functions = NULL, 
                       coef.t., delta.t., N.max., bru.opt){
  # Expected number of background events
  df.0 <- data.frame(counts = 0, exposures = 1)
  # this is the expression of log(Lambda0)
  form.0 <- counts ~ log(link.functions$mu(th.mu)) + log(T2 - T1) 
  
  # first likelihood
  lik.0 <- inlabru::like(formula = form.0,
                         data = df.0,
                         family = 'poisson',
                         options = list(E = df.0$exposures))
  
  cat('Start creating spatial grid...', '\n')
  time.g.st <- Sys.time()
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    time.grid(data.point = sample.s[idx,], 
              coef.t = coef.t., 
              delta.t = delta.t., 
              T2. = T2, N.exp. = N.max.
              )
  } %>%
    mutate(counts = 0,
           exposures = 1)
  cat('Finished creating spatial grid, time ', Sys.time() - time.g.st, '\n')   
  
  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.p, df_grid, 
                               M0, ncore_ = ncore){
    theta_ <- c(0, 
                link.functions$K(th.K[1]), 
                link.functions$alpha(th.alpha[1]), 
                link.functions$c_(th.c[1]), 
                link.functions$p(th.p[1]))
    
    #cat('theta - LogL', theta_, '\n')
    comp.list <- compute.grid(param. = theta_, gridd = df_grid)
    #print(sum(is.na(comp.list$It)))
    #print(sum(is.infinite(comp.list$It)))
    out <- theta_[3]*(df_grid$magnitudes - M0) + log(theta_[2] + 1e-100) + log(comp.list$It + 1e-100) 
  }
  
  # creating formula for past events contributions to integrated lambda
  form.j.part <- counts ~ logLambda.h.inla(th.K = th.K, th.alpha = th.alpha,
                                           th.c = th.c, th.p = th.p,
                                           df_grid = df.j, M0 = M0)
  # second for triggered part of the integral
  lik.j.part <- inlabru::like(formula = form.j.part,
                              data = df.j,
                              family = 'poisson',
                              options = list(E = df.j$exposures)) 
  
  # third is for the sum of the log intensities
  df.s <- data.frame(ts = sample.s$ts, counts = 1, exposures = 0)
  
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, tt, th, mh, M0){
    
    if(is.null(link.functions)){
      th.p <- c(th.mu[1], th.K[1], th.alpha[1], th.c[1], th.p[1])  
    }
    else{
      th.p <- c(link.functions$mu(th.mu[1]),
                link.functions$K(th.K[1]),
                link.functions$alpha(th.alpha[1]),
                link.functions$c_(th.c[1]),
                link.functions$p(th.p[1]))
    }
    
    unlist(mclapply(1:length(tt), \(idx) 
                    log(lambda_2(th = th.p, t = tt[idx], ti.v = th[th < tt[idx]], 
                                mi.v = mh[th < tt[idx]], M0 = M0)),
                    mc.cores = 5))
  }
  
  # creating formula for summation part
  form.s.part <- counts ~ loglambda.inla(th.mu = th.mu, th.K = th.K, th.alpha = th.alpha,
                                         th.c = th.c, th.p = th.p, tt = ts, th = sample.s$ts, 
                                         mh = sample.s$magnitudes, M0 = M0)
  
  lik.s.part <- inlabru::like(formula = form.s.part,
                              data = df.s,
                              family = 'poisson',
                              options = list(E = df.s$exposures)) 
  
  cmp.part <- counts ~ -1 + 
    th.mu(1, model = 'linear', mean.linear = 0 , prec.linear = 1) + 
    th.K(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.alpha(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.c(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.p(1, model = 'linear', mean.linear = 0, prec.linear = 1) 
  
  bru(lik.0, lik.j.part, lik.s.part, components = cmp.part, 
      options = bru.opt)
  
  
}






Hawkes.bru2 <- function(sample.s, M0, T1, T2, link.functions = NULL, 
                        N.breaks.min = 3, max.length_ = 0.1, 
                        coef_exp = 1, N_exp = 5, t.lim = 10,
                        N.tail = 4, bin.strat = 'default', bru.opt){
  # Expected number of background events
  df.0 <- data.frame(counts = 0, exposures = 1)
  # this is the expression of log(Lambda0)
  form.0 <- counts ~ log(link.functions$mu(th.mu)) + log(T2 - T1) 
  
  # first likelihood
  lik.0 <- inlabru::like(formula = form.0,
                         data = df.0,
                         family = 'poisson',
                         options = list(E = df.0$exposures))
  # Expected number of triggered events
  df.j <- foreach(idx = 1:length(sample.s$ts), .combine = rbind) %do% {
    tt <- sample.s$ts[idx]
    if(bin.strat == 'default'){
      t.max <- min(tt + t.lim, T2)
      N.breaks <- max(N.breaks.min, ceiling((t.max - tt)/max.length_))
      kk <- seq_len(N.breaks) - 1
      Time.bins.matrix <- tt + cbind(kk, kk + 1) / N.breaks * (t.max - tt)
      if(t.max != T2){
        leng <- (T2 - t.max)/N.tail
        tail.b <- seq(t.max, T2, by = max(leng, max.length_))
        Time.bins.matrix <- rbind(Time.bins.matrix, 
                                  cbind(tail.b[-length(tail.b)], tail.b[-1]))  
        }
      } 
    else if(bin.strat == 'exponential'){
      tt_break <- breaks_exp(tt_ = tt, T2_ = T2, coef_ = coef_exp, delta_ = max_length_,
                             N_exp_ = N_exp)
      Time.bins.matrix <- cbind(tt_break[-length(tt_break)], tt_break[-1])
    } else {
        stop('Unknown binning strategy')
      }
    
    data.frame(ts = rep(sample.s$ts[idx], each = nrow(Time.bins.matrix)),
               mags = rep(sample.s$magnitudes[idx], each = nrow(Time.bins.matrix)),
               counts = 0,
               exposures = 1,
               bin.start = Time.bins.matrix[,1],
               bin.end = Time.bins.matrix[,2])
  }
  
  logLambda.h.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, ts, mags, T1.v, T2.v, M0){
    
    if(is.null(link.functions)){
      th.p <- c(th.mu[1], th.K[1], th.alpha[1], th.c[1], th.p[1])  
    }
    else{
      th.p <- c(link.functions$mu(th.mu[1]),
                link.functions$K(th.K[1]),
                link.functions$alpha(th.alpha[1]),
                link.functions$c_(th.c[1]),
                link.functions$p(th.p[1]))
    }
    
    out <- unlist(mclapply(1:length(ts), \(idx)
                    log.Lambda_h2(th = th.p, ti = ts[idx], mi = mags[idx],
                                  M0 = M0, T1 = T1.v[idx], T2 = T2.v[idx]),
                    mc.cores = 5))
    
  }
  
  
  # creating formula for past events contributions to integrated lambda
  form.j.part <- counts ~ logLambda.h.inla(th.mu = th.mu, th.K = th.K, th.alpha = th.alpha,
                                           th.c = th.c,th.p = th.p,
                                           ts = ts, mags = mags,
                                           T1.v = bin.start, T2.v = bin.end, M0 = M0)
  # second for triggered part of the integral
  lik.j.part <- inlabru::like(formula = form.j.part,
                              data = df.j,
                              family = 'poisson',
                              options = list(E = df.j$exposures)) 
  
  # third is for the sum of the log intensities
  df.s <- data.frame(ts = sample.s$ts, counts = 1, exposures = 0)
  
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, tt, th, mh, M0){
    
    if(is.null(link.functions)){
      th.p <- c(th.mu[1], th.K[1], th.alpha[1], th.c[1], th.p[1])  
    }
    else{
      th.p <- c(link.functions$mu(th.mu[1]),
                link.functions$K(th.K[1]),
                link.functions$alpha(th.alpha[1]),
                link.functions$c_(th.c[1]),
                link.functions$p(th.p[1]))
    }
    
    unlist(mclapply(1:length(tt), \(idx) 
                    log(lambda_2(th = th.p, t = tt[idx], ti.v = th[th < tt[idx]], 
                                 mi.v = mh[th < tt[idx]], M0 = M0)),
                    mc.cores = 5))
  }
  
  # creating formula for summation part
  form.s.part <- counts ~ loglambda.inla(th.mu = th.mu, th.K = th.K, th.alpha = th.alpha,
                                         th.c = th.c,th.p = th.p, tt = ts, th = sample.s$ts, 
                                         mh = sample.s$magnitudes, M0 = M0)
  
  lik.s.part <- inlabru::like(formula = form.s.part,
                              data = df.s,
                              family = 'poisson',
                              options = list(E = df.s$exposures)) 
  
  cmp.part <- counts ~ -1 + 
    th.mu(1, model = 'linear', mean.linear = 0 , prec.linear = 1) + 
    th.K(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.alpha(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.c(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.p(1, model = 'linear', mean.linear = 0, prec.linear = 1) 
  
  bru(lik.0, lik.j.part, lik.s.part, components = cmp.part, 
      options = bru.opt)
  
  
}


## functions to get the density from a sample and build a copula transformation

q.K.mcmc <- function(pr, lower.tail = TRUE, log.p = FALSE){
  if(lower.tail){
    if(log.p){
      return(quantile(k.bru, exp(pr)))
    }
    else{
      return(quantile(k.bru, pr))
    }
  } 
  if(!lower.tail){
    if(log.p){
      return(quantile(k.bru, 1 - exp(pr)))
    }
    else{
      return(quantile(k.bru, 1 - pr))
    }
  } 
}


mcmc.k.t <- function(x){
  bru_forward_transformation(q.K.mcmc, x)
}






