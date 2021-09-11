library(ggplot2)
library(viridis)
library(inlabru)
library(INLA)
library(dplyr)
library(data.table)
library(metR)
library(matrixStats)
library(parallel)
library(mvtnorm)
library(MASS)
library(raster)

# stable magnitude triggering function.
# th.n and th.p are thresholds such that 
# if theta3 + log(m - M0) < th.n we approximate linearly around th.n
# if theta3 + log(m - M0) > th.p we approximate linearly around th.p
g.m <- function(m, theta3, M0, th.n = -29, th.p = 4){
  if(any(m - M0 < 0)){
    stop('m smaller than completeness threshold')
  }
  # if too negative
  if(any(theta3 + log(m - M0) < th.n)){
    theta3.0 <- th.n - log(m - M0)
    log.app <- theta3.0 + log(m - M0)
    exp(exp(log.app)) + exp(exp(log.app))*exp(theta3.0)*(m - M0)*(theta3 - theta3.0)
  }
  # if too positive
  else if(any(theta3 + log(m - M0) > th.p)){
    theta3.0 <- th.p - log(m - M0)
    log.app <- theta3.0 + log(m - M0)
    exp(exp(log.app)) + exp(exp(log.app))*exp(theta3.0)*(m - M0)*(theta3 - theta3.0)
  }
  else{
    log.app <- theta3 + log(m - M0)
    exp(exp(log.app))
  }
  exp(exp(theta3)*(m - M0))
}


g.t <- function(tt, theta4, theta5, Ht.times, th.p = 5, th.n = -20){
  time.diff <- tt - Ht.times
  output <- rep(0, length(time.diff))
  idx.p <- time.diff > 0
  logoutput.pos <- (-exp(theta5) - 1)*log(time.diff[idx.p] + exp(theta4))

  if(any(logoutput.pos > th.p)){
    idx.toolarge <- logoutput.pos > th.p
    output.app <- exp(logoutput.pos)
    output.app[idx.toolarge] <- exp(th.p) + (logoutput.pos[idx.toolarge] - th.p)*exp(th.p)
    output[idx.p] <- output.app
    output
  }
  else if(any(logoutput.pos < th.n)){
    idx.tooneg <- logoutput.pos < th.n
    output.app <- exp(logoutput.pos)
    output.app[idx.tooneg] <- exp(th.n) + (logoutput.pos[idx.tooneg] - th.n)*exp(th.n)
    output[idx.p] <- output.app
    output
  }
  else{
    output[idx.p] <- exp(logoutput.pos)
    output
  }
  exp(logoutput.pos)
}


# spatial triggering
g.s <- function(loc, Ht.locs, Sigma.p){
  if(!is.matrix(Sigma.p)){
    stop('Sigma not a matrix')
  }
  if(ncol(Sigma.p) != nrow(Sigma.p)){
    stop('Sigma wrong dimensions')
  }
  det.S <- det(Sigma.p)
  if(det.S <= 0){
    stop('Sigma not positive semi-definite')
  }
  out <- sapply(1:nrow(Ht.locs), function(x) 
    dmvnorm(loc, mean = as.numeric(Ht.locs[x,]), sigma = Sigma.p))
  out
}


# QUESTION : for the sampling with the thinning technique we do not account for the magnitude
#            if we would account for the magnitude of each event also, it would be very difficult to retain
#            any high magnitude event. However, I feel it is the right way to implement this thing.


# intensity function without considering magnitude of event
# this for now assumes a homogeneous background rate
logLambda <- function(ppoints, theta.vec, Ht, M0){
  tt <- as.numeric(ppoints[1])
  loc <- as.numeric(ppoints[2:3])
  
  theta1 <- theta.vec[1]
  theta2 <- theta.vec[2]
  theta3 <- theta.vec[3]
  theta4 <- theta.vec[4]
  theta5 <- theta.vec[5]
  sigmax <- exp(theta.vec[6])
  sigmay <- exp(theta.vec[7])
  sigmaxy <- exp(theta.vec[8])
  
  Sigma.p <- matrix(c(sigmax, sigmaxy, sigmaxy, sigmay), ncol = 2, byrow = TRUE)
  past.locs <- cbind(Ht$long, Ht$lat)
  past.t <- Ht$ts
  past.m <- Ht$mags
  # magnitude triggering
  
  gm.v <- as.vector(g.m(past.m, theta3, M0))
  gt.v <- as.vector(g.t(tt, theta4, theta5, past.t))
  gs.v <- as.vector(g.s(loc, past.locs, Sigma.p))
  
  # calculate summation
  log.trigs <- theta2 + log(sum(gm.v*gt.v*gs.v))
  # return intensity
  
  log(exp(theta1) + exp(log.trigs))
  #logSumExp(c(theta1, log.trigs))
  
}


## function to sample 

point_sampler <- function(loglambda, bdy, mesh, num_events, crsobj = NULL){
  ## Number of events for single catalogue from a poisson distribution with lambda = num_events
  #cat_n <- rpois(1, num_events)
  loglambda_max <- max(loglambda)
  ## number of points to sample at a time - might want to adjust depending on how many points you want to actually retain.
  n.points = 10000
  
  ## Set up a spatialpoints dataframe for our results
  samp.points <- SpatialPoints(data.frame(x = 0, y = 0))
  samp.points$mags <- 0
  samp.points <- samp.points[-1,]
  if(!is.null(crsobj)){
    proj4string(samp.points) <- crsobj}
  num <- 0
  
  ## To sample the correct number of points, keep going until the num >= cat_n
  while (num < num_events){
    points <- spsample(bdy, n.points, "random")
    ## transform to wgs84 long/lat
    if(!is.null(crsobj)){pts <- spTransform(points, crsobj)}
    else{pts <- points}
    ## next few lines modified directly from sample.lgcp
    proj <- INLA::inla.mesh.project(mesh, pts)
    lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - loglambda_max)
    keep <- proj$ok & (runif(n.points) <= lambda_ratio)
    kept <- pts[keep]
    kept$mags <- rep(0, length(kept))
    if(length(kept) > 0){
      samp.points <- rbind(samp.points, kept)
      num <- length(samp.points)
    }
    else{
      #print('no retained')
    }
  }
  
  ## Keep exactly cat_n points, choose these randomly from all of the points we've kept so far
  kp <- sample(seq(1, length(samp.points), by=1), num_events, replace=FALSE)
  samp.points <- samp.points[kp,]
  ## Get magnitudes for this catalogue
  #samp.points$mags <- TapGRM(cat_n, b_val, 8, m_min)
  
  return(samp.points)
}




# function to calculate the integral of the time-triggering function
# between Ht and Tlim. It can be seen as the points generated by the time effect between Ht.ts and Tlim

# input : Tlim --> Extreme of the interval of time on which calculate the integral
#         theta4,5 --> parameters of the Omori's law (theta5 = log(p-1))
#         Ht --> observed time points 

# output : numeric - The integral of the UNNORMALIZED Omori's law 
# notes : reliable only if AT MOST one input is a vector

I.t <- function(Tlim, theta4, theta5, Ht.ts){
  exp(-theta5)*(exp(-theta4*exp(theta5)) - (Tlim - Ht.ts + exp(theta4))^(-exp(theta5)))
}


## DOC: integral of the spatial triggering function.
# we are using a Gaussian density for which it is easy to retrieve

I.s <- function(x.lim, y.lim, Sigma, past.loc){
  pmvnorm(lower = c(x.lim[1], y.lim[1]), upper = c(x.lim[2], y.lim[2]),
          mean = as.numeric(past.loc), sigma = Sigma)[1]
}

# function that calculate I.s when Ht contains multiple points
# the observed points DO NOT need to be in the domain indicated by xlim/ylim

I.s.multi <- function(x.lim, y.lim, Sigma, past.loc.multi){
  sapply(1:nrow(past.loc.multi), function(x) 
    I.s(x.lim, y.lim, Sigma, past.loc.multi[x,]))
}

## to sample from Omoris

Inv.I.t <- function(Lambda, theta4, theta5, Ht.ts){
  (exp(-theta4*exp(theta5)) + Lambda*(-exp(theta5)))^(-1/exp(theta5)) + Ht.ts - exp(theta4) 
}


sample.omori <- function(n.ev, theta4, theta5, Ht.ts, Tlim){
  if(n.ev == 0){return('no samples')}
  bound <- I.t(Tlim, theta4, theta5, Ht.ts)
  ss <- runif(n.ev, min = 0, max = bound)
  t.sample <- Inv.I.t(ss, theta4, theta5, Ht.ts)
  if(all(t.sample > Tlim)){
    return('no samples')}
  sort(t.sample[t.sample < Tlim])
}


## two functions to sample the a GIVEN number of location according to the Gaussian 
# triggering function. 
# They dispose the number of observations generated by a past location

# the last one can be IMPROVED considering that it samples 1 point for each cell, 
# however if I already know that I need a xx samples in a cell I can sample 
# those points all together.

sample.from.idx <- function(idx, xy.grid, cell.size){
  cell <- as.numeric(xy.grid[idx, 1:2])
  c(x = runif(1, min = cell[1] - cell.size/2, max = cell[1] + cell.size/2),
    y = runif(1, min = cell[2] - cell.size/2, max = cell[2] + cell.size/2))
}


sample.spatial <- function(n.ev, x.lim, y.lim, nx, ny, Sigma, past.loc){
  # create sequences
  x.seq <- seq(x.lim[1], x.lim[2], length.out = nx)
  y.seq <- seq(y.lim[1], y.lim[2], length.out = ny)
  
  # retriece centroids
  x.cent <- x.seq[-nx] + diff(x.seq)/2
  y.cent <- y.seq[-ny] + diff(y.seq)/2
  
  # initialize grid
  xy.grid <- expand.grid(x.cent, y.cent)
  colnames(xy.grid) <- c('x', 'y')
  
  # calculate prob of being in each bin
  xy.grid$probs <- 
    sapply(1:nrow(xy.grid), function(i) 
      pmvnorm(lower = as.numeric(xy.grid[i,1:2] - diff(x.seq)/2),
              upper = as.numeric(xy.grid[i,1:2] + diff(x.seq)/2),
              mean = past.loc, sigma = Sigma)[1])
  
  xy.grid$probs[xy.grid$probs < 0] <- 0
  # sample the cells
  cell.idx <- sample(1:nrow(xy.grid), n.ev, prob = xy.grid$probs, replace = TRUE)
  
  # sample uniformly from each bin
  ss <- matrix(unlist(lapply(cell.idx, function(id) 
    sample.from.idx(id, xy.grid, diff(x.seq)))), byrow = T, ncol = 2)
  
  # return data.frame of locations 
  data.frame(x = ss[,1], y = ss[,2])
}


#crs=crs_wgs84


# main issues:

# - find a way to calculate mahalanobis distance using CRS
# - the sample points function should sample the points only where the loglambda is over a certain threshold
#   or points in a certain distance

# The function to generate an ETAS sample should : 
# 1) generate observations triggered by past seismicity
# 2) sample background
# 3) sample triggered by background and so on

# 
# g.s <- function(locs, mean.p, Sigma.p){
#   if(!is.matrix(Sigma.p)){
#     stop('Sigma not a matrix')
#   }
#   if(ncol(Sigma.p) != nrow(Sigma.p)){
#     stop('Sigma wrong dimensions')
#   }
#   det.S <- det(Sigma.p)
#   if(det.S <= 0){
#     stop('Sigma not positive semi-definite')
#   }
#   inv.S <- ginv(Sigma)
#   distances <-  pointDistance(locs, mean.p, lonlat = TRUE)/1000 #distance in km
#   
#   (1/(2*pi*sqrt(det.S)))*exp(-(1/2)*t(distances)%*%inv.S%*%distances)
# }
# 
# library(MASS)
# library(raster)
# 
# mm <- SpatialPoints(data.frame(x = 0, y = 0))
# proj4string(mm) <- crswgs
# mm <- spTransform(mm, crswgs)
# 
# pointDistance(a, mm, lonlat = TRUE)
# pointDistance(a, mm, lonlat = FALSE)
# 
# plot(bdy)
# points(a)
# points(mm)
# 
# g.s(a, mm, Sigma)
# 
# 
# #Sigma = matrix(c(0.2,0.1,0.1,0.2), byrow = T, ncol = 2)
# 
# Sigma = matrix(c(0.3, 0.29, 0.29, 0.3), byrow = T, ncol = 2)
# 
# ss <- Sys.time()
# log.l <- log(dmvnorm(mesh.n$loc[,1:2], mean = c(0.5, 0.5), sigma = Sigma))
# aa <- point_sampler(log.l, bdy, mesh.n, num_events = 1000)
# print(Sys.time() - ss)
# 
# mm <- inla.mesh.projector(mesh.n, dims = c(20, 20))
# mm2 <- inla.mesh.project(mm, exp(log.l))  
# mm3 <- pixels(mesh.n, nx = 20, ny = 20)
# mm3$canepone <- as.vector(mm2)
# 
# 
# ss <- Sys.time()
# aa2 <- sample.spatial(n.ev = 1000, x.lim = c(0,1), y.lim = c(0,1), nx = 30, ny = 30, 
#                       Sigma = Sigma, past.loc= c(0.5, 0.5))
# print(Sys.time() - ss)
# 
# 
# pl1 <- ggplot() + gg(mm3) + gg(aa) + theme(legend.position = 'bottom')
# pl2 <- ggplot() + gg(mm3) + geom_point(data = aa2, mapping = aes(x = x, y = y)) + 
#   theme(legend.position = 'bottom')
# 
# multiplot(pl1, pl2, cols = 2)
# 
# pl1 <- ggplot(data = as.data.frame(aa), aes(x = x, y = y)) + geom_density_2d()
# pl2 <- ggplot(data = aa2, aes(x = x, y = y)) + geom_density_2d()
# 
# multiplot(pl1, pl2, cols = 2)
# 
  






# now theta.par has 8 components: mu, K, alpha, c, p
# and sigmax2 sigmay2 sigmaxy

# Ht has to be such that Ht[,1:2] locations,  Ht[,3] time, Ht[,4] magnitudes
expected.n.triggered <- function(theta.par, Tlim, M0, 
                                 x.lim, y.lim, Ht.multi){
  
  theta2 = theta.par[2]
  theta3 = theta.par[3]
  theta4 = theta.par[4]
  theta5 = theta.par[5]
  
  sigmax = exp(theta.par[6])
  sigmay = exp(theta.par[7])
  sigmaxy = theta.par[8]
  Sigma <- matrix(c(sigmax, sigmaxy, sigmaxy, sigmay), byrow = TRUE, ncol = 2)
  
  past.locs <- Ht.multi[,1:2]
  past.ts <- Ht.multi[,3]
  past.mags <- Ht.multi[,4]
  
  log.n <- theta2 + exp(theta3)*(past.mags - M0) + 
    log(I.t(Tlim, theta4, theta5, past.ts)) + 
    log(I.s.multi(x.lim, y.lim, Sigma, past.locs))
 
  exp(as.numeric(log.n))
}

# function to sample points generated by a past observation

sample.triggered <- function(nn, theta.par, beta.par, Sigma, Tlim, M0,
                             x.lim, y.lim, nx, ny, Ht.single, gen.idx){
  if(nrow(Ht.single) > 1){
    stop('Error - Multiple past events')
  }
  past.ts <- Ht.single$ts
  past.loc <- c(Ht.single$x, Ht.single$y)  
  ss <- sample.spatial(nn, x.lim, y.lim, nx, ny, Sigma, as.numeric(past.loc))
  
  dd <- data.frame(x = ss[,1], y = ss[,2],
                   ts = sample.omori(nn, theta.vec[4], theta.vec[5], 
                                     past.ts, Tlim = Tlim),
                   mags = rexp(nn, beta.par) + M0,
                   gen = gen.idx)
  dd[order(dd$ts), ]
}





sample.ETAS <- function(theta.par, beta.par, M0, Tlim, 
                        x.lim, y.lim, nx, ny, n.core = 5){
  theta1 <- theta.par[1]
  theta2 <- theta.par[2]
  theta3 <- theta.par[3]
  theta4 <- theta.par[4]
  theta5 <- theta.par[5]
  
  sigmax = exp(theta.par[6])
  sigmay = exp(theta.par[7])
  sigmaxy = theta.par[8]
  Sigma <- matrix(c(sigmax, sigmaxy, sigmaxy, sigmay), byrow = TRUE, ncol = 2)
  
  # extract number of backgroud events
  Lambda.back <- exp(theta1)*Tlim*(diff(x.lim)*diff(y.lim))
  n.backg <- rpois(1, lambda = Lambda.back)
  print(c(n.backg, 0))
  # extract times and locations homogeneously
  backg.t <- runif(n.backg, min = 0, max = Tlim)
  backg.mags <- rexp(n.backg, beta.par) + M0
  backg.s <- cbind(runif(n.backg, min = x.lim[1], max = x.lim[2]),
                   runif(n.backg, min = y.lim[1], max = y.lim[2]))
  
  # store them
  backg.df <- data.frame(x = backg.s[,1],
                         y = backg.s[,2],
                         ts = backg.t,
                         mags = backg.mags,
                         gen = 0)
  
  # Initialize past.gens list with events sorted by time.
  Past.Gens <- list(backg.df[order(backg.df$ts),])
  
  flag = T
  gen.idx = 1
  while(flag){
    # take previous generation
    past.events <- Past.Gens[[gen.idx]]
    past.events <- data.frame(past.events)
    n.past.ts <- nrow(past.events)
    
   
    # for each of them sample the number of offsprings
    Lambdas <- expected.n.triggered(theta.par, Tlim, M0, x.lim, y.lim,
                                    past.events)
    
    n.off <- rpois(n.past.ts, lambda = Lambdas)
    print(n.off)
    # exit condition if none of the event in the previous generation has
    # offsprings
    if(all(n.off == 0)){
      Gens.df <- bind_rows(Past.Gens)
      Gens.df <- Gens.df[order(Gens.df$ts),]
      return(Gens.df)
    }
    # select only times which has generated events
    past.events <- past.events[n.off > 0,]
    n.off <- n.off[n.off > 0]
    
    # generate offsprings for each past event
    
    offsprings.list <- mclapply(1:length(n.off), function(x)
      sample.triggered(n.off[x], theta.par, beta.par, Sigma, Tlim, M0,
                       x.lim, y.lim, nx, ny, Past.Gens[[gen.idx]][x,] ,
                       gen.idx), mc.cores = n.core)

    # offsprings.list <- list()
    # 
    # for(i in 1:nrow(past.events)){
    #   print('###')
    #   print(n.off[i])
    #   print(past.events[i,])
    #   offsprings.list[[i]] <- sample.triggered(n.off[i], theta.par, beta.par, 
    #                                            Sigma, Tlim, M0,
    #                                            x.lim, y.lim, nx, ny, 
    #                                            past.events[i,] ,
    #                                            gen.idx)
    #   print('-------')
    # }
    # 
    
    
    offsprings.df <- rbindlist(offsprings.list)
    
    offsprings.df <-
      offsprings.df %>%
      filter(ts != 'no samples') %>%
      mutate(ts = as.numeric(ts))
    
    print(c(nrow(offsprings.df), gen.idx))
    # update gen.idx
    gen.idx <- gen.idx + 1
    Past.Gens[[gen.idx]] <- offsprings.df
  }
}


# theta.vec <- c(3, log(0.02), log(1.5), log(0.01), log(0.3), 
#                log(0.001), log(0.001), 0)
# 
# 
# Sigma <- matrix(c(0.001, 0, 0, 0.001), byrow = T, ncol = 2)
# toplot.gs(c(0.5, 0.5), Sigma, x.seq, y.seq)
# 
# 
# beta.par <- 2.3
# M0 = 2.5
# Tlim = 10
# x.lim = c(0,1)
# y.lim = c(0,1)
# nx = 100
# ny = 100
# 
# expected.n.triggered(theta.vec, Tlim, M0, x.lim, y.lim, 
#                      data.frame(x = 0.2, y = 0.2, ts = 1, mags = 3))
# 
# dd <- sample.ETAS(theta.vec, beta.par, M0, Tlim, x.lim, y.lim, nx, ny)
# 
# nrow(dd)
# 
# hist(dd$ts)
# ggplot(dd, aes(x = ts, y = mags)) + geom_point()
# 
# pl1.s <- ggplot(dd, aes(x = x, y = y, color = as.factor(gen))) + geom_point() + 
#   xlim(0,1) + ylim(0,1) + theme(legend.position = 'none') + 
#   labs(title = paste0('N = ', nrow(dd)))
# 
# 
# theta.vec2 <- c(3, log(0.02), log(1.5), log(0.01), log(0.3), 
#                log(0.01), log(0.01), 0.009)
# 
# Sigma <- matrix(c(0.01, 0.009, 0.009, 0.01), byrow = T, ncol = 2)
# toplot.gs(c(0.5, 0.5), Sigma, x.seq, y.seq)
# 
# expected.n.triggered(theta.vec2, Tlim, M0, x.lim, y.lim, 
#                      data.frame(x = 0.2, y = 0.2, ts = 1, mags = 3))
# 
# dd2 <- sample.ETAS(theta.vec2, beta.par, M0, Tlim, x.lim, y.lim, nx, ny)
# 
# nrow(dd2)
# 
# hist(dd2$ts)
# ggplot(dd2, aes(x = ts, y = mags)) + geom_point()
# 
# pl2.s <- ggplot(dd2, aes(x = x, y = y, color = as.factor(gen))) + geom_point() + 
#   xlim(0,1) + ylim(0,1) + theme(legend.position = 'none') +
#   labs(title = paste0('N = ', nrow(dd2)))
# 
# 
# multiplot(pl1.s, pl2.s, cols = 2)
# 
# 
# 
# sample.lgcp
# 
# 
# 
