brier <- function(x, p){
  -2 * (p - x)^2
}

log.sc <- function(x,p){
  x * log(p) + (1 - x) * log(1 - p)
}


brier.score.second.der <- function(pst){-4*pst -4*(1-pst)}
log.score.second.der <- function(p, pst){-pst/(p^2) -(1-pst)/((1-p)^2)}

# parimutuel gambling score k = 2
gamb <- function(x, p1, p2){
  p0 <- (p1 + p2)/2
  x * (p1 / p0 - 1) + (1 - x) * ((1 - p1) / (1 - p0) - 1)
}

# expected parimutuel gambling score for k forecasts, 
# takes as input the mean of the other k-1 forecasts
gambling.score <- function(p1, pst, pm1, k){
  p0 = (p1 + (k-1)*pm1)/k
  score = (p1 - p0)*(pst - p0)/(p0*(1-p0))
  return(score)
}

# expected score difference as above but pm1 now is the mean of k-2 forecasts
gambling.score.diff <- function(p1, p2, pst, pm1, k){
  p0 = (p1 + p2 + (k-2)*pm1)/k
  score = (p1 - p2)*(pst - p0)/(p0*(1-p0))
  return(score)
}

# function for plotting (used for fig2)
toplot.gambling <- function(p1.v, pst, pm1.v, k){
  final <- c()
  id.df <- c() 
  p1.df <- c() 
  for(p in pm1.v){
    gamb <- sapply(p1.v, 
                   function(x) gambling.score(x, pst, p, k))
    final <- c(final, gamb)
    id.df <- c(id.df, rep(paste('p2 = ', round(p,4)), length(p1.v)))
    p1.df <- c(p1.df, p1.v)
  }
  df <- data.frame(score = -final/min(final), 
                   id = id.df,
                   p1 = p1.df)
  
  ggplot(df, aes(x = p1, y = score, linetype = id, color = id)) +
    geom_line(size = 0.8) +
    ylab("Expected Score") +
    geom_hline(yintercept = 0)
}


# function to plot the score differences between p1 and p2 in a k=3 comparison
# (used for fig3)
toplot.gambling.diff <- function(p1.v, p2, pst, p3.v, k){
  scores <- c()
  p.id <- c()
  p1 <- c()
  for(p3 in p3.v) {
    gamb <- sapply(p1.v, 
                   function(x) gambling.score.diff(x, p2, pst, p3, k))
    
    scores <- c(scores, gamb)
    p.id <- c(p.id, rep(paste('p3 = ', round(p3, 4)), length(p1.v)))
    p1 <- c(p1, p1.v)
  }  
  
  df <- data.frame(score = scores, 
                   id = p.id,
                   p1 = p1)
  
  ggplot(df, aes(x = p1, y = -score/min(score), linetype = id, color = id)) + geom_line(size = 0.8) +
    geom_hline(yintercept = 0)
}


# same thing of function above but for different ks instead of different 
# values of the average k-2 forecasts (fig 4)

toplot.gambling.k <- function(p1.v, p2, pst, p3, k.v){
  scores <- c()
  p.id <- c()
  p1 <- c()
  for(k in k.v) {
    gamb <- sapply(p1.v, 
                   function(x) gambling.score.diff(x, p2, pst, p3, k))
    
    scores <- c(scores, gamb)
    p.id <- c(p.id, rep(paste('k =', k), length(p1.v)))
    p1 <- c(p1, p1.v)
  }  
  
  p.id <- factor(p.id, levels = unique(p.id))
  df <- data.frame(score = scores, 
                   id = p.id,
                   p1 = p1)
  
  ggplot(df, aes(x = p1, y = score,
                 linetype = id, color = id)) + geom_line(size = 0.8) +
    geom_vline(xintercept = pst) + geom_hline(yintercept = 0) 
}

# function to calculate the score differences between p1, p2 
# p0 reference model for the pairwise (fig 5)

deltas <- function(p1,p2,p0,score.FUN){
  
  if(is.function(score.FUN)){
    Delta0 = score.FUN(0,p1) - score.FUN(0,p2)
    Delta1 = score.FUN(1,p1) - score.FUN(1,p2)
  }
  else if(score.FUN == 'pair.gamb'){
    Delta0 = gamb(0,p1,p0) - gamb(0,p2,p0)
    Delta1 = gamb(1,p1,p0) - gamb(1,p2,p0)
  }
  else if(score.FUN == 'full.gamb'){
    Delta0 = gamb(0,p1,p2) - gamb(0,p2,p1)
    Delta1 = gamb(1,p1,p2) - gamb(1,p2,p1)
  }
  else{
    stop('Wrong score function')
  }
  cbind(Delta0, Delta1)
}


# function to find the interval of the observed number of active bins 
# such that we are not able to take a decision. (fig 5)

Xs.extremes <- function(delta01, N.b, xs){
  
  d.lower <- delta01[1] + 
    qbeta(0.025, xs, N.b - xs + 1)*(delta01[2]- delta01[1])
  
  d.obs <- delta01[1] + 
    (xs/N.b)*(delta01[2]- delta01[1])
  
  d.upper <- delta01[1] + 
    qbeta(0.975, xs + 1, N.b - xs)*(delta01[2]- delta01[1])
  
  # check upper and lower bound
  if(d.lower[1] < d.upper[1]){
    dd <- data.frame(lower = d.lower,
                     obs = d.obs,
                     upper = d.upper, 
                     Xs = xs)  
  }
  if(d.lower[1] > d.upper[1]){
    dd <- data.frame(lower = d.upper,
                     obs = d.obs,
                     upper = d.lower, 
                     Xs = xs)  
  }
  
  # check if the score diff goes up or down (it depends if p1 >< p2)
  
  if(dd$obs[1] - dd$obs[2] < 0){
    
    X.min <- min(dd$Xs[dd$upper > 0])
    X.max <- max(dd$Xs[dd$lower < 0])
  }
  
  if(dd$obs[1] - dd$obs[2] > 0){
    X.min <- min(dd$Xs[dd$lower < 0])
    X.max <- max(dd$Xs[dd$upper > 0])
  }
  
  if(is.infinite(X.min)){X.min = NULL}
  if(is.infinite(X.max)){X.max = NULL}
  
  return(list(extremes = c(X.min, X.max),
              df = dd))
}

# Single Probability Multiple bins case
# function to calculate the preference probabilities given, number of bins,
# p1, p2, p0, pstar, score function. (fig 6)
Find.probs <- function(N.b,p1,p2,p0,pstar,score.FUN){
  
  delta01 <- deltas(p1,p2,p0, score.FUN)
  
  xs.sup <- qbinom(0.999, size = N.b, prob = pstar)
  
  xs = 1:max(xs.sup,5)
  
  res <- Xs.extremes(delta01, N.b, xs)
  
  Xs.max <-  res$extremes[2]
  Xs.min <-  res$extremes[1]
  
  if(is.na(Xs.min)){Xs.min = 0}
  if(is.na(Xs.max)){Xs.max = max(xs)}
  prob.p1 <- 1 - pbinom(Xs.max, size = N.b, prob = pstar) 
  prob.p2 <- pbinom((Xs.min - 1), size = N.b, prob = pstar)  
  prob.not <- sum(dbinom(Xs.min:Xs.max, size = N.b, prob = pstar)) 
  return(c('no pref' = prob.not, p1 = prob.p1, p2 = prob.p2))
}

# function to plot beta = prob of expressing a preference as function of pstar

toplot.beta <- function(p1, p2, p0, pstar.v, N.vec, score.FUN, leg.pos){
  
  dd <- sapply(pstar.v, 
               function(x) Find.probs(N.vec[1], p1, p2, p0, x, score.FUN))
  
  df <- data.frame(probs = dd[1,],
                   N = N.vec[1],  
                   pstar = pstar.v)
  
  if(length(N.vec) >= 2){
    for(i in 2:length(N.vec)){
      dd <- sapply(pstar.v, 
                   function(x) Find.probs(N.vec[i], p1, p2, p0, x, score.FUN))
      
      df2 <- data.frame(probs = dd[1,],
                        N = N.vec[i],  
                        pstar = pstar.v)
      df <- rbind(df, df2)
    }  
  }
  
  
  
  df <- df[order(df$N),]
  df$N <- as.factor(df$N)
  ggplot(df, aes( x= pstar, y = 1 - probs, color = N,
                  linetype = N)) + 
    
    geom_line() + 
    geom_vline(xintercept = p1, linetype = 3) +
    geom_vline(xintercept = p2, linetype = 3) +
    annotate('text', x = p2 + 0.0001/2, y = 0.9, label = TeX('p_2')) + 
    annotate('text', x = p1 + 0.0001/2, y = 0.9, label = TeX('p_1')) +
    labs(color = 'N') + 
    xlab(TeX('p^*')) + 
    ylab(expression(beta)) + 
    theme_classic() + 
    theme(legend.background = element_rect(#fill= 'lightgrey',
      linetype="solid", 
      colour ="darkgrey"),
      legend.key.size = unit(1, 'line'),
      axis.title.y = element_text(size = 14),
      legend.position = leg.pos)  
  
}


# function to calculate expected score difference considering p1 = pstar,
# p2 = alph*pstar, p0 = alph0*pstar (tab3)
exp.delta <- function(probs, alph, alph0, score.FUN){
  
  p1 = probs
  p2 = probs*alph
  p0 = probs*alph0
  
  delta01 <- deltas(p1,p2,p0,score.FUN)  
  mean(delta01[,1]) + mean(probs*(delta01[,2] - delta01[,1]))
}

# function to calculate full gambling score considering k = 3 (fig11)
fullgamb <- function(x, p1, p2, p3){
  pm <- (p1 + p2 + p3)/3
  x * (p1 / pm - 1) + (1 - x) * ((1 - p1) / (1 - pm) - 1)
}


# function to calculate the confidence intervals using normalapprox Rhoades2011 (fig12)
runn.normapprox <- function(pstar, p1, p2, p0){
  N = length(pstar)
  
  x0 = rbinom(N, 1, pstar)
  
  delta.b <- brier(x0, p1) - brier(x0, p2)
  delta.l <- log.sc(x0, p1) - log.sc(x0, p2)
  delta.pg <- gamb(x0, p1, p0) - gamb(x0, p2, p0)
  delta.fg <- gamb(x0, p1, p2) - gamb(x0, p2, p1)
  
  edelta.b <- mean(delta.b); sd.b <- sd(delta.b)
  edelta.l <- mean(delta.l); sd.l <- sd(delta.l)
  edelta.pg <- mean(delta.pg); sd.pg <- sd(delta.pg)
  edelta.fg <- mean(delta.fg); sd.fg <- sd(delta.fg)
  
  t.alpha <- qt(0.025, N - 1)
  CI.b <- sort(c(edelta.b - t.alpha*(sd.b/sqrt(N)),
                 edelta.b,
                 edelta.b + t.alpha*(sd.b/sqrt(N))))
  
  CI.l <- sort(c(edelta.l - t.alpha*(sd.l/sqrt(N)),
                 edelta.l,
                 edelta.l + t.alpha*(sd.l/sqrt(N))))
  
  CI.pg <- sort(c(edelta.pg - t.alpha*(sd.pg/sqrt(N)),
                  edelta.pg,
                  edelta.pg + t.alpha*(sd.pg/sqrt(N))))
  
  CI.fg <- sort(c(edelta.fg - t.alpha*(sd.fg/sqrt(N)),
                  edelta.fg,
                  edelta.fg + t.alpha*(sd.fg/sqrt(N))))
  
  c(CI.b, CI.l, CI.pg, CI.fg, sum(x0))
  
}

# function to store the results produced by runn.normapprox
storing <- function(sim.n){
  CI.b <- data.frame(Lower = sim.n[,1],
                     Obs = sim.n[,2],
                     Upper = sim.n[,3])
  
  CI.l <- data.frame(Lower = sim.n[,4],
                     Obs = sim.n[,5],
                     Upper = sim.n[,6])
  
  CI.pg <- data.frame(Lower = sim.n[,7],
                      Obs = sim.n[,8],
                      Upper = sim.n[,9])
  
  CI.fg <- data.frame(Lower = sim.n[,10],
                      Obs = sim.n[,11],
                      Upper = sim.n[,12])
  
  list(Brier = CI.b, Log = CI.l, PG = CI.pg, FG = CI.fg, xS = sim.n[,13])
}

# function to plot average score difference density (fig13) 
toplot.ddensity <- function(Obs, truth, x.lim){
  dd <- data.frame(Obs = Obs)
  dd.m <- mean(dd$Obs)
  dd.sd <- sd(dd$Obs)
  xx <- seq(x.lim[1], x.lim[2], length.out = 1000)
  norm.df <- data.frame(Obs = xx,
                        pdf = dnorm(xx, dd.m, dd.sd)/max(dnorm(xx, dd.m, dd.sd)))
  
  ggplot(dd, aes(x = Obs, y=..scaled..)) + geom_density() + 
    geom_line(data = norm.df, aes(y = pdf), linetype = 2, color = 2) +
    geom_vline(xintercept = truth, linetype = 3) + 
    xlab(~Delta) + 
    ylab('density') +
    xlim(x.lim)
}

## FUNCTION TO SIMULATE FROM ETAS IN PARALLEL (super effiient)
# (fig15)
new.sim.etas.par <- function(mu.par, k.par, a.par, c.par, p.par, beta.par, M0, Tlim){
  
  # sample from background
  N0 <- rpois(1, mu.par*Tlim)
  if(N0 == 0){ return('No observations') }
  prop.t <- seq(0, Tlim, length.out = Tlim*1000)
  gen0 <- sample(prop.t, size = N0)
  mag0 <- rexp(N0, beta.par) + M0
  
  Gens <- list(data.frame(ts = gen0, magnitudes = mag0, gen = 0))
  i = 1
  flag = TRUE
  while(flag){
    Prev.Gen <- Gens[[i]]
    print(c(Gen = i, N.ev = nrow(Prev.Gen)))
    
    # for each element of the previous generation
    sim.from.obs <- function(past.o){
      past.ts <- Prev.Gen$ts[past.o]
      past.ms <- Prev.Gen$magnitudes[past.o]
      N.past.o <- rpois(1, k.par*exp(a.par*(past.ms - M0)))
      
      # set max possible simulated time
      Tprop.lim <- past.ts - c.par + ( (1e-6)/((p.par - 1)*(c.par^(p.par - 1))) )^(-1/p.par)
      # set possible values
      prop.gen.t <- seq(past.ts + 1e-8, Tprop.lim, length.out = 1000*(Tprop.lim - past.ts))
      # set probabilities
      probs.t <- (c.par^(p.par - 1))*(p.par - 1)/((prop.gen.t - past.ts + c.par)^p.par)
      probs.t <- probs.t*diff(prop.gen.t)[1]
      # sample times
      Sample.past.ts <- sample(prop.gen.t, size = N.past.o, prob = probs.t)
      
      # if there is a sampled time smaller than the limit
      if(sum(Sample.past.ts < Tlim) > 0){
        # restrict to Tlim
        Sample.past.ts <- Sample.past.ts[Sample.past.ts < Tlim]
        
        # sample magnitudes
        Sample.past.ms <- rexp(length(Sample.past.ts), beta.par) + M0
        # stores in current Generation
        
        data.frame(ts = Sample.past.ts, 
                   magnitudes = Sample.past.ms,
                   gen = i)  
      }
      # else goes to another observation
      else{ NA }
    }
    
    Current.Gen.list <- mclapply(1:nrow(Prev.Gen), function(x) sim.from.obs(x),
                                 mc.cores = 5)
    Current.Gen.list <- Current.Gen.list[!is.na(Current.Gen.list)]
    Current.Gen <- c()
    # now we should have filled up the Current.Gen
    # if there are no events we can stop to simulate
    if(length(Current.Gen.list) == 0){
      flag = FALSE
    }
    # else we have to store it in Gen and update the index
    else{
      for(j in 1:length(Current.Gen.list)){
        Current.Gen <- rbind(Current.Gen, Current.Gen.list[[j]])
      }
      Gens[[i + 1]] <- Current.Gen
      i = i + 1
    }
  }
  data.frame(ts = unlist(lapply(Gens, function(x) x$ts)),
             magnitudes = unlist(lapply(Gens, function(x) x$magnitudes)),
             gen = unlist(lapply(Gens, function(x) x$gen)))
}


# function to calculate the confidence intervals for 
# the expected score diff comparing two models using catalog-based fore
# (fig16)
runner.probs <- function(n.bins, SIM.data1, SIM.data2, opt = 1, Tlim = 100){
  
  
  sum.s <- foreach(i = 1:length(SIM.data1), .combine = '+') %do% {
    
    if(is.null(SIM.data1[[i]])){
      output1 <- rep(0, n.bins)
    }
    else{
      ss.i1 <- SIM.data1[[i]]$ts#[SIM.data1[[i]]$magnitudes > 4.95]
      ss.i1 <- ss.i1[ss.i1 < Tlim]
      hh.i1 <- hist(ss.i1, breaks = seq(0,Tlim, length.out = (n.bins + 1)), plot = FALSE)
      output1 <- (hh.i1$counts > 0)
    }
    if(is.null(SIM.data2[[i]])){
      output2 <- rep(0, n.bins)
    }
    else{
      ss.i2 <- SIM.data2[[i]]$ts#[SIM.data2[[i]]$magnitudes > 4.95]
      ss.i2 <- ss.i2[ss.i2 < Tlim]
      hh.i2 <- hist(ss.i2, breaks = seq(0,Tlim,length.out = (n.bins + 1)), plot = FALSE)
      output2 <- (hh.i2$counts > 0)
    }
    
    c(output1, output2)
  }
  
  uncond.probs1 <- sum.s[1:n.bins]/length(SIM.data1)
  uncond.probs2 <- sum.s[(n.bins+1):length(sum.s)]/length(SIM.data2)
  uncond.probs1[uncond.probs1 < 1e-7] <- 1e-7; uncond.probs1[uncond.probs1 > 1 - 1e-7] <- 1 - 1e-7
  uncond.probs2[uncond.probs2 < 1e-7] <- 1e-7; uncond.probs2[uncond.probs2 > 1 - 1e-7] <- 1 - 1e-7
  
  log.sc <- .GlobalEnv$log.sc
  gamb <- .GlobalEnv$gamb

  CI.res <- foreach(i = 1:length(SIM.data1), .export = c('log.sc', 'gamb'),
                    .combine = 'rbind') %do% {
                      if(opt == 1){
                        ss.i <- SIM.data1[[i]]$ts#[SIM.data1[[i]]$magnitudes > 4.95]  
                        #print(max(ss.i))
                        ss.i <- ss.i[ss.i < Tlim]
                        #print(max(ss.i))
                      }
                      else{
                        ss.i <- SIM.data2[[i]]$ts#[SIM.data2[[i]]$magnitudes > 4.95]
                        ss.i <- ss.i[ss.i < Tlim]
                      }
                      
                      
                      hh.i <- hist(ss.i, breaks = seq(0,Tlim,length.out = (n.bins + 1)), plot = FALSE)
                      xx <- (hh.i$counts > 0)
                      #print(sum(xx))
                      t.alpha <- qt(0.025, n.bins  - 1)
                      
                      # log.score
                      diff.l <- log.sc(xx, uncond.probs1) - log.sc(xx, uncond.probs2)
                      delta.l <- mean(diff.l)
                      deltaL.l <- mean(diff.l) + t.alpha*sd(diff.l)/sqrt(n.bins)
                      deltaU.l <- mean(diff.l) - t.alpha*sd(diff.l)/sqrt(n.bins)
                      
                      # pg score
                      diff.pg <- gamb(xx, uncond.probs1, uncond.probs1*5) - gamb(xx, uncond.probs2,
                                                                                 uncond.probs1*5)
                      delta.pg <- mean(diff.pg)
                      deltaL.pg <- mean(diff.pg) + t.alpha*sd(diff.pg)/sqrt(n.bins)
                      deltaU.pg <- mean(diff.pg) - t.alpha*sd(diff.pg)/sqrt(n.bins)
                      
                      # fg score
                      diff.fg <- gamb(xx, uncond.probs1, uncond.probs2) - gamb(xx, uncond.probs2,
                                                                               uncond.probs1)
                      delta.fg <- mean(diff.fg)
                      deltaL.fg <- mean(diff.fg) + t.alpha*sd(diff.fg)/sqrt(n.bins)
                      deltaU.fg <- mean(diff.fg) - t.alpha*sd(diff.fg)/sqrt(n.bins)
                      
                      # output row
                      c(deltaL.l, delta.l, deltaU.l,
                        deltaL.pg, delta.pg, deltaU.pg,
                        deltaL.fg, delta.fg, deltaU.fg, 
                        length(ss.i), max(hh.i$counts), mean(hh.i$counts))
                      
                    }
  CI.res <- data.frame(CI.res)
  colnames(CI.res) <- c('lower.l', 'obs.l', 'upper.l',
                        'lower.pg', 'obs.pg', 'upper.pg',
                        'lower.fg', 'obs.fg', 'upper.fg',
                        'N.oss', 'Max.count', 'Av.count')
  #print("a")
  CI.res
}

# extract the preferene probabilities (fig16)
extract.probs <- function(res.run){
  prob.l <- c(p1 = mean(res.run$lower.l > 0),
              p2 = mean(res.run$upper.l < 0))
  prob.pg <- c(p1 = mean(res.run$lower.pg > 0),
               p2 = mean(res.run$upper.pg < 0))
  prob.fg <- c(p1 = mean(res.run$lower.fg > 0),
               p2 = mean(res.run$upper.fg < 0))
  list(Log = prob.l, PG = prob.pg, FG = prob.fg)
}





sim.p <- function(p1,p2,p0,pst){
  n <- length(pst)
  x <- rbinom(n,1,pst)
  d.b <- data.frame(p1 = p1,
                    p2 = p2,
                    pst = pst,
                    p0 = p0,
                    obs = x,
                    p1.brier = brier(x, p1),
                    p2.brier = brier(x, p2),
                    p1.log = log.sc(x, p1),
                    p2.log = log.sc(x, p2),
                    p1.pg = gamb(x, p1, p0),
                    p2.pg = gamb(x, p2, p0),
                    p1.fg = gamb(x, p1, p2),
                    p2.fg = gamb(x, p2, p1))
  d.b
}


# 
# pst <- agg.prob#runif(100, 1e-3, 0.2)
# p1 <- pst
# p2 <- pst/2
# p0 <- min(pst*5, 0.98)
# 
# 
# resc. <- function(x){ 
#   x01 <- (x - min(x))/(max(x) - min(x))
#   2*x01 - 1
#   }
# 
# sim.ex <- sim.p(p1,p2,p0,pst) 
# mean(sim.ex$p1.brier) - mean(sim.ex$p2.brier)
# 
# df.p <- rbind(data.frame(p = sim.ex$p1,
#                          brier = sim.ex$p1.brier,
#                          log = sim.ex$p1.log,
#                          pg = sim.ex$p1.pg,
#                          fg = sim.ex$p1.fg,
#                          obs = sim.ex$obs,
#                          fore = 'p1'),
#               data.frame(p = sim.ex$p2,
#                          brier = sim.ex$p2.brier,
#                          obs = sim.ex$obs,
#                          pg = sim.ex$p2.pg,
#                          fg = sim.ex$p1.fg,
#                          log = sim.ex$p2.log,
#                          fore = 'p2'))
# 
# df.p$int <- interaction(df.p$obs, df.p$fore)
# df.p$brier <- resc.(df.p$brier)
# df.p$log <- resc.(df.p$log)
# df.p$pg <- resc.(df.p$pg)
# df.p$fg <- resc.(df.p$fg)
# 
# 
# 
# pl.b <- ggplot(df.p, aes(x = int, y = brier, fill = fore)) + 
#   geom_boxplot() + 
#   geom_point(aes(x = 3/2,
#                  y = 
#                mean(brier[fore == 'p1'])),
#              color = 'darkred',
#              size = 3) + 
#   geom_point(aes(x = 2 + 3/2,
#                  y = 
#                    mean(brier[fore == 'p2'])),
#              color = 'darkblue',
#              size = 3)+ 
#   ylim(-1,1)
# 
# 
# pl.l <- ggplot(df.p, aes(x = int, y = log, fill = fore)) + 
#   geom_boxplot() + 
#   geom_point(aes(x = 3/2,
#                  y = 
#                    mean(log[fore == 'p1'])),
#              color = 'darkred',
#              size = 3) + 
#   geom_point(aes(x = 2 + 3/2,
#                  y = 
#                    mean(log[fore == 'p2'])),
#              color = 'darkblue',
#              size = 3)+ 
#   ylim(-1,1)
# 
# 
# pl.pg <- ggplot(df.p, aes(x = int, y = pg, fill = fore)) + 
#   geom_boxplot() + 
#   geom_point(aes(x = 3/2,
#                  y = 
#                    mean(pg[fore == 'p1'])),
#              color = 'darkred',
#              size = 3) + 
#   geom_point(aes(x = 2 + 3/2,
#                  y = 
#                    mean(pg[fore == 'p2'])),
#              color = 'darkblue',
#              size = 3)+ 
#   ylim(-1,1)
# 
# pl.fg <- ggplot(df.p, aes(x = int, y = fg, fill = fore)) + 
#   geom_boxplot() + 
#   geom_point(aes(x = 3/2,
#                  y = 
#                    mean(fg[fore == 'p1'])),
#              color = 'darkred',
#              size = 3) + 
#   geom_point(aes(x = 2 + 3/2,
#                  y = 
#                    mean(fg[fore == 'p2'])),
#              color = 'darkblue',
#              size = 3) + 
#   ylim(-1,1)
# 
# d.b <- mean(sim.ex$p1.brier) - mean(sim.ex$p2.brier)
# mean(sim.ex$p1.log) - mean(sim.ex$p2.log)
# mean(sim.ex$p1.pg) - mean(sim.ex$p2.pg)
# de.fg <- mean(sim.ex$p1.fg) - mean(sim.ex$p2.fg)
# 
# multiplot(pl.b, pl.l, pl.pg, pl.fg, 
#           layout = matrix( c(1,4,2,3), byrow = T, ncol = 2))
# 
# 
# ggplot(sim.ex, aes(x = log)) + geom_density()
# 
# boxplot(brier ~ obs,sim.ex)
# boxplot(log ~ obs,sim.ex)
# 
# ## next step try test Diaboldi-Mariano per BRIER intervals.
# 
# table(sim.15$xS)
# 
# ggplot(data.frame(x = sim.15$xS,
#                   y = sim.15$Brier$Obs), 
#        aes(x = as.factor(x), y= y)) + geom_boxplot()
# 
# ggplot(data.frame(x = sim.4$xS,
#                   y = sim.4$Log$Obs), 
#        aes(x = as.factor(x), y= y)) + geom_boxplot()
# 
# ggplot(data.frame(x = sim.4$xS,
#                   y = sim.4$PG$Obs), 
#        aes(x = as.factor(x), y= y)) + geom_boxplot()
# 
# 
# 
