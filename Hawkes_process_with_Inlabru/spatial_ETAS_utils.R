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
# mandatory input : --> m : magnitude value (has to be m > M0)
#                 : --> theta3 : value of the parameter (theta3 = log(alpha))
#                 : --> M0 : completeness magnitude
# optional input  : --> th.n : threshold s.t. if theta3 + log(m - M0) < th.n we approx linearly around th.n
#                 : --> th.p : threshold s.t. if theta3 + log(m - M0) > th.p we approx linearly around th.p

# ouput : numeric -value of g.m(m, theta3) = exp(exp(theta3)*(m - M0))

# notes : if m is scalar and theta3 is a vector the output is a vector with components g.m(m, theta3_i)
#         if theta3 is scalar and m is a vector the output is a vector with components g.m(m_i, theta3)
#         if both theta3 and m are vector (same length required) the output is a vector g.m(m_i, theta3_i)
#         same rationale applies to M0
#         parameters has to be numeric

g.m <- function(m, theta3, M0, th.n = -29, th.p = 4){
  # stop if magnitude lower than magnitude of completeness
  if(any(m - M0 < 0)){
    stop('m smaller than completeness threshold')
  }
  # if too negative
  if(any(theta3 + log(m - M0) < th.n)){
    theta3.0 <- th.n - log(m - M0)
    log.app <- theta3.0 + log(m - M0)
    # linear approximation around theta3.0
    exp(exp(log.app)) + exp(exp(log.app))*exp(theta3.0)*(m - M0)*(theta3 - theta3.0)
  }
  # if too positive
  else if(any(theta3 + log(m - M0) > th.p)){
    theta3.0 <- th.p - log(m - M0)
    log.app <- theta3.0 + log(m - M0)
    # linear approximation around theta3.0
    exp(exp(log.app)) + exp(exp(log.app))*exp(theta3.0)*(m - M0)*(theta3 - theta3.0)
  }
  else{
    # compute the log of the exponent
    log.app <- theta3 + log(m - M0)
    # return double exponential (a bit dangerous, might be changed)
    exp(exp(log.app))
  }
  exp(exp(theta3)*(m - M0))
}


# function to plot g.m as function of the magnitude for different values of theta3
toplot.gm <- function(theta3.seq, mm, M0){
  gm.list <- lapply(theta3.seq, function(x) data.frame(mags = mm,
                                                       gm = g.m(mm, x, M0),
                                                       theta3 = x))
  gm.df <- bind_rows(gm.list)
  
  ggplot(gm.df, aes(x = mags, y = gm, color = as.factor(theta3), 
                    linetype = as.factor(theta3))) + 
    geom_line() + 
    labs(color = "theta3", linetype = 'theta3')
  
}



# stable time triggering function.
# mandatory input : --> tt : time value (t)
#                 : --> theta4 : value of the parameter (theta4 = log(c))
#                 : --> theta5 : value of the parameter (theta5 = log(p - 1))
#                 : --> past obervation (t_i)
# optional input  : --> th.p : threshold s.t. if (-exp(theta5) - 1)*log(t - t_i + exp(theta4)) > th.p 
#                              we approx linearly around th.p
#                 : --> th.n : threshold s.t. if (-exp(theta5) - 1)*log(t - t_i + exp(theta4)) < th.n 
#                              we approx linearly around th.n

# ouput : numeric -value of g.t(t - t_i, theta4, theta5) = I(t-t_i > 0)*(t - t_i + exp(theta4))^(-(exp(theta5) + 1))
#         where I(t - t_i > 0) indicator funcion, it is 1 when the condition is true
# notes : if t is scalar and t_i is a vector the output is a vector with components g.m(t - t_i, theta4, theta5)
#         if t_i is scalar and t is a vector the output is a vector with components g.m(t_j - t_i, theta4, theta5)
#         if both t and t_i are vector (same length required) the output is a vector g.m(t_j - t_ij, theta4, theta5)
#         same rationale applies to theta4, theta5
#         parameters has to be numeric

g.t <- function(tt, theta4, theta5, Ht.times, th.p = 5, th.n = -20){
  # compute the difference in time (t - t_i)
  time.diff <- tt - Ht.times
  # prepare output (its 0 if t - t_i < 0)
  output <- rep(0, length(time.diff))
  # identify positive time differences
  idx.p <- time.diff > 0
  # compute log_output only for positive time differences
  logoutput.pos <- (-exp(theta5) - 1)*log(time.diff[idx.p] + exp(theta4))
  # approximate if too high
  if(any(logoutput.pos > th.p)){
    # update only positions above threshold
    idx.toolarge <- logoutput.pos > th.p
    output.app <- exp(logoutput.pos)
    output.app[idx.toolarge] <- exp(th.p) + (logoutput.pos[idx.toolarge] - th.p)*exp(th.p)
    # update only positive time differences positions
    output[idx.p] <- output.app
    output
  }
  # approximate if too low
  else if(any(logoutput.pos < th.n)){
    # update only positions below threshold
    idx.tooneg <- logoutput.pos < th.n
    output.app <- exp(logoutput.pos)
    output.app[idx.tooneg] <- exp(th.n) + (logoutput.pos[idx.tooneg] - th.n)*exp(th.n)
    # update only positive time differences positions
    output[idx.p] <- output.app
    output
  }
  else{
    # update only positive time differences positions
    output[idx.p] <- exp(logoutput.pos)
    output
  }
}


# spatial triggering function (Bidimensional Gaussian density)
# mandatory input : --> loc : 2D-location (s)
#                 : --> Ht.locs : observed 2D-locations (s_i), used as mean of the Gaussian kernel
#                 : --> Sigma.p : 2x2 covariance matrix composed by sigma^2_x, sigma^2_y, sigma_xy
#                                

# output : numeric -value of the function g.s(s - s_i, Sigma)

# notes: Sigma.p has to e a valid covariance matrix (positive-semidefinite)
#      : Ht.locs has to be a matrix or data.frame with only two columns representing the locations
#        those locations will be used as means of the Gaussian kernel          

g.s <- function(loc, Ht.locs, Sigma.p){
  # stop if Sigma not a matrix
  if(!is.matrix(Sigma.p)){
    stop('Sigma not a matrix')
  }
  # stop if Sigma not square matrix
  if(ncol(Sigma.p) != nrow(Sigma.p)){
    stop('Sigma wrong dimensions')
  }
  det.S <- det(Sigma.p)
  # stop if Sigma not a suitable covariance matrix
  if(det.S <= 0){
    stop('Sigma not positive semi-definite')
  }
  # compute the Gaussian Kernel for row of Ht.locs
  out <- sapply(1:nrow(Ht.locs), function(x) 
    dmvnorm(loc, mean = as.numeric(Ht.locs[x,]), sigma = Sigma.p))
  out
}


# function to plot the spatial triggering function g.s in space given two sequences of values x.seq and y.seq
toplot.gs <- function(Mu, Sigma, x.seq, y.seq){
  
  df <- expand.grid(x.seq, y.seq)
  colnames(df) <- c('x', 'y')
  df$density <- sapply(1:nrow(df), function(x) dmvnorm(df[x,1:2], mean = Mu,
                                                       sigma = Sigma))
  
  ggplot(df, aes(x = x, y = y, fill = density, z = density)) + 
    geom_tile() + #geom_contour() + 
    scale_fill_viridis() + 
    geom_point(aes(x = Mu[1], y = Mu[2])) + 
    labs(title = paste0('Sigma.xy = ', Sigma[2,1])) + 
    theme(legend.position = 'bottom')
}


# function to calculate the conditional log-intensity
# mandatory input : --> ppoints : points at which evaluate the conditional log-intensity (see notes)
#                 : --> theta.vec : vector of parameters (see notes)
#                 : --> Ht : history of the process (see notes)                                              
#                 : --> M0 : completeness magnitude

# output : numeric -value of the conditional log-intensity log(lambda(ppoints | Ht))

# notes : theta.vec = (theta1 = log(mu), theta2 = log(K), theta3 = log(alpha), theta4 = log(c), theta5 = log(p - 1), 
#                      theta6 = sigma^2_x, theta7 = sigma^2_y, theta8 = sigma^2_xy
#       : ppoints has to be a data.frame with time column called "ts" and location "x", "y"
#       : Ht have the same format of ppoints with the addition of the magnitude column called "mags"

logLambda <- function(ppoints, theta.vec, Ht, M0){
  # extract time and space location
  tt <- as.numeric(ppoints$ts)
  loc <- c(ppoints$x, ppoints$y)
  
  # extract parameters
  theta1 <- theta.vec[1]
  theta2 <- theta.vec[2]
  theta3 <- theta.vec[3]
  theta4 <- theta.vec[4]
  theta5 <- theta.vec[5]
  sigmax <- exp(theta.vec[6])
  sigmay <- exp(theta.vec[7])
  sigmaxy <- exp(theta.vec[8])
  
  # construct Sigma 
  Sigma.p <- matrix(c(sigmax, sigmaxy, sigmaxy, sigmay), ncol = 2, byrow = TRUE)
  
  # extract information on past observations
  past.locs <- cbind(Ht$long, Ht$lat)
  past.t <- Ht$ts
  past.m <- Ht$mags
  
  # magnitude triggering
  gm.v <- as.vector(g.m(past.m, theta3, M0))
  # time triggering
  gt.v <- as.vector(g.t(tt, theta4, theta5, past.t))
  # space triggering
  gs.v <- as.vector(g.s(loc, past.locs, Sigma.p))
  
  # calculate log of the summation
  log.trigs <- theta2 + log(sum(gm.v*gt.v*gs.v))
  # return logintensity intensity
  log(exp(theta1) + exp(log.trigs))
  # for a more stable version this should be right.
  # logSumExp(c(theta1, log.trigs))
  
}


# function to calculate the integral of the time-triggering function
# mandatory input : --> T1 : start of the time interval
#                 : -->  T2 : end of the time interval
#                 : --> theta4,5 --> parameters of the Omori's law (theta5 = log(p-1))
#                 : --> Ht.ts --> observed time points 

# output : numeric - The integral of the unnormalized Omori's law given by
# int_[T1,T2] (t - t_i + exp(theta4))^(-(exp(theta5) + 1))dt

# notes : if one of the inputs is a vector and the others are scalars the function returns a vector
#       : if more than one input is a vector they have to be of the same dimension, the ouput is a vector

I.t <- function(T1, T2, theta4, theta5, Ht.ts){
  # if the observation is after the time interval returns zero
  if(Ht.ts > T2){return(0)}
  # lower extreme is the minum between T1 and observation, before the latter the integrand is zero
  T.low <- max(Ht.ts, T1)
  # evaluate the integral
  -exp(-theta5)*((T2 - Ht.ts + exp(theta4))^(-exp(theta5)) - (T.low - Ht.ts + exp(theta4))^(-exp(theta5)))
}

## Inverse function of I.t (used to perform inversion sampling for sampling times)
# mandatory inputs : --> : It : value to be transformed
#                  : --> : theta4, theta5 : value of the parameters of gt
#                  : --> Ht.ts : past time observation (t_i)

# ouput : numeric - considering It(t1) = int_[0,t1] gt(t-t_i)dt 
#                   the outpur is InvIt(It) such that InvIt(It(t1)) = t1 

Inv.I.t <- function(It, theta4, theta5, Ht.ts){
  (exp(-theta4*exp(theta5)) + It*(-exp(theta5)))^(-1/exp(theta5)) + Ht.ts - exp(theta4) 
}


# function to place aftershocks in time
# mandatory inputs : --> n.ev : number of events to be placed
#                  : --> theta4,5 : parameters of the time-triggering function (gt(t - t_i))
#                  : --> Ht.ts : past time location (t_i)
#                  : --> T1, T2 : extremes of the time interval (see notes)

# output : numeric - sample of n.ev times from a point process with intensity lambda(t) = gt(t- t_i)

# notes : the sampled times are between T1 and T2
sample.omori <- function(n.ev, theta4, theta5, Ht.ts, T1, T2){
  if(n.ev == 0){return(NA)}
  bound.l <- I.t(Ht.ts, T1, theta4, theta5, Ht.ts)
  bound.u <- I.t(Ht.ts, T2, theta4, theta5, Ht.ts)
  unif.s <- runif(n.ev, min = bound.l, max = bound.u)
  t.sample <- Inv.I.t(unif.s, theta4, theta5, Ht.ts)
  t.sample
}



# function to calculate the integral of the space-triggering function
# mandatory input : --> x.lim : two values representing range of values of x
#                 : --> y.lim : two values representing range of values of y
#                 : --> Sigma : 2x2 covariance matrix
#                 : --> past.loc : observed locations (s_i) 

# output : numeric - The integral of the unnormalized Omori's law given by
# int_[A] gs(s - s_i)ds where A is a square with locations in the range x.lim, y.lim

# notes : the inputs have to be numeric, no matrix allowed
#       : past.loc can be outside the boundaries

I.s <- function(x.lim, y.lim, Sigma, past.loc){
  pmvnorm(lower = c(x.lim[1], y.lim[1]), upper = c(x.lim[2], y.lim[2]),
          mean = as.numeric(past.loc), sigma = Sigma)[1]
}

# function that calculate I.s when Ht contains multiple points
# past.loc.multi has to be a matrix or data.frame with only two columns representing the locations

I.s.multi <- function(x.lim, y.lim, Sigma, past.loc.multi){
  sapply(1:nrow(past.loc.multi), function(x) 
    I.s(x.lim, y.lim, Sigma, past.loc.multi[x,]))
}


# function to calculate the number of aftershock triggered by a past observation
# mandatory inputs : --> theta.v : parameters vector
#                  : --> Ht.single : a single observation (see notes)
#                  : --> T1 : start of the time interval
#                  : --> T2 : end of the time interval
#                  : --> M0 : magnitude of completeness
#                  : --> bdy : polygon representing the area (A) (see notes) 

# output : expected number of points triggered by an event (Ht.single = (t_i, s_i, m_i)) 
# given by the integral int_[T1,T2] int_[A] K*gm(m_i)*gt(t - t_i)*gs(s - s_i)dtds

# notes : Ht.single has to be a row of a data.frame with columns 
#         "ts" (time), "x" , "y" (location) and "mags" (magnitude)
#       : bdy has to be a SpatialPolygon -- for now it calculate the integral with respect to space over a square
#         region circoscribing the input polygon.


triggered.n.ev <- function(theta.v, Ht.single, T1, T2, M0, bdy){
  
  # extract parameters
  theta2 <- theta.v[2]
  theta3 <- theta.v[3]
  theta4 <- theta.v[4]
  theta5 <- theta.v[5]
  sigmax <- exp(theta.v[6])
  sigmay <- exp(theta.v[7])
  sigmaxy <- exp(theta.v[8])
  # construct Sigma
  Sigma <- matrix(c(sigmax, sigmaxy, sigmaxy, sigmay), ncol = 2, byrow = TRUE)
  
  # Try to calculate the integral of the space-triggering - sometimes FAIL unexpectedly - see weird behavior
  totry <- I.s(bdy@bbox[1,],bdy@bbox[2,],Sigma, cbind(Ht.single$x, Ht.single$y))
  # if it is NA we assume that the integral of the space-triggering is 1 (which is its maximum)
  if(is.na(totry)){
    return(exp(theta2)*exp(exp(theta3)*(as.numeric(Ht.single$mags) - M0))*
             I.t(T1, T2, theta4,theta5, as.numeric(Ht.single$ts)))
  }
  else{
    return(exp(theta2)*exp(exp(theta3)*(as.numeric(Ht.single$mags) - M0))*
             I.t(T1, T2, theta4,theta5, as.numeric(Ht.single$ts))*totry)
  }
}



## function to place in space the aftershocks of a past event
# mandatory input : --> n.ev : number of aftershocks to be generated
#                 : --> bdy : area in which they have to be generated
#                 : --> Ht.loc: past observation generating the aftershocks
#                 : --> Sigma: 2x2 covariance matrix for space-triggering function
#                 : --> crsobj : CRS object to project the locations

# output : SpatialPoints- set of n.ev points extracted from a point process with intensity 
#                         given by gs(s - s_i)


sample.given.loc <- function(n.ev, bdy, Ht.loc, Sigma, crsobj){
  
  # initialize empty SpatialPoints
  samp.points <- SpatialPoints(data.frame(x = 0, y = 0))
  samp.points <- samp.points[-1,]
  proj4string(samp.points) <- crswgs
  # initialize number of placed events
  num <- 0
  
  # until we placed all the events
  while (num < n.ev){
    # sample from the 2D Gaussian kernel without boundaries
    pts.matrix <- rmvnorm(n.ev, mean = as.numeric(Ht.loc), sigma = Sigma)
    # transform the sample in SpatialPoints and project
    pts <- data.frame(x = pts.matrix[,1], y = pts.matrix[,2])
    coordinates(pts) <- ~ x + y
    proj4string(pts) <- crsobj
    pts <- spTransform(pts, crsobj)
    # discard the ones outside bdy
    pts <- crop(pts, bdy)
    # merge sampled points
    if(length(pts) > 0){
      samp.points <- rbind(samp.points, pts)
      num <- length(samp.points)
    }
    else{
      #print('no retained')
    }
  }
  
  ## if we retained a number of points > n.ev, we select n.ev events at random
  kp <- sample(seq(1, length(samp.points), by=1), n.ev, replace=FALSE)
  samp.points <- samp.points[kp,]
  samp.points
}



# function to place in space and time aftershocks of one event
# mandatory inputs : --> theta.v --> parameters vector
#                  : --> beta.p --> parameter of the GR law for the magnitude
#                  : --> n.ev : number of aftershocks to be generated
#                  : --> Ht.single: past observation generating the aftershocks
#                  : --> Sigma: 2x2 covariance matrix for space-triggering function
#                  : --> crsobj : time interval in which events have to be placed
#                  : --> bdy : area in which events have to be placed
#                  : --> crsobj : CRS object to project the locations

# output : data.frame - columns are time "ts", magnitude "mags" and location "x", "y", they are a sample of 
#                       n.ev events from a point process with intensity 
#                       lamdba(t,s,m) = beta.p*exp(-beta.p*(m - M0))*gs(s - s_i)*g(t-t_i)

# notes : Ht.single has to be a data.frame formatted as the output ("ts" = time, "mags" = magnitude.
#                                                                   "x", "y" = location)

sample.triggered <- function(theta.v, beta.p, n.ev, Ht.single, T1, T2, M0, bdy, crsobj){
  # if the number of events to be placed is zero returns an empty data.frame
  if(n.ev == 0){
    samp.points <- data.frame(x = 1, y = 1, ts = 1, mags = 1)
    samp.points <- samp.points[-1,]
    return(samp.points)
  }
  else{
    # initialize parameters
    theta4 <- theta.v[4]
    theta5 <- theta.v[5]
    sigmax <- exp(theta.v[6])
    sigmay <- exp(theta.v[7])
    sigmaxy <- exp(theta.v[8])
    # construct Sigma
    Sigma <- matrix(c(sigmax, sigmaxy, sigmaxy, sigmay), ncol = 2, byrow = TRUE)
    # extract information about parent
    Ht.loc <- c(Ht.single$x, Ht.single$y)
    Ht.ts <- as.numeric(Ht.single$ts)
    Ht.mags <- as.numeric(Ht.single$mags)
    
    # sample times
    samp.ts <- sample.omori(n.ev, theta4, theta5, Ht.ts, T1, T2)
    # sample magnitudes
    samp.mags <- rexp(n.ev, rate = beta.p) + M0
    # sample locations
    samp.locs <- sample.given.loc(n.ev, bdy, Ht.loc, Sigma, crsobj)
    
    # build output dataset
    samp.points <- data.frame(ts = samp.ts, mags = samp.mags, x = samp.locs@coords[,1],
                              y = samp.locs@coords[,2])
    # return only the ones with time different from NA (the one with NA are outside the interval T1, T2)
    # even though it should not happen given how we built sample.omori
    return(samp.points[!is.na(samp.points$ts),])
  }
  
}

# sample a generation of aftershocks from a set of parents event Ht
# mandatory inputs : --> theta.v --> parameters vector
#                  : --> beta.p --> parameter of the GR law for the magnitude
#                  : --> n.ev : number of aftershocks to be generated
#                  : --> Ht: past observations generating the aftershocks
#                  : --> Sigma: 2x2 covariance matrix for space-triggering function
#                  : --> M0 : completeness magnitude
#                  : --> bdy : area in which events have to be placed
#                  : --> crsobj : CRS object to project the locations

# optional parameters : --> ncore : number of cores used for parallelization

# output : data.frame - columns are time "ts", magnitude "mags" and location "x", "y", they are a sample of 
#                       n.ev events from a point process with intensity 
#                       lamdba(t,s,m) = beta.p*exp(-beta.p*(m - M0))*gs(s - s_i)*g(t-t_i)

# notes : Ht.single has to be a data.frame formatted as the output ("ts" = time, "mags" = magnitude.
#                                                                   "x", "y" = location)


sample.generation <- function(theta.v, beta.p, Ht, T1, T2, M0, bdy, crsobj, ncore = 2){
  # number of parents
  n.parent <- nrow(Ht)
  
  # calculate the aftershock rate for each parent in history
  trig.rate.v <- sapply(1:n.parent, function(x) 
    triggered.n.ev(theta.v, Ht[x,], T1, T2, M0, bdy))
  # extract number of aftershock for each parent
  n.ev.v <- sapply(1:n.parent, function(x) rpois(1, trig.rate.v[x]))
  
  # if no aftershock has to be generated returns empty data.frame
  if(sum(n.ev.v) == 0){
    app <- data.frame(x = 1, y = 1, ts = 1, mags = 1)
    app <- app[-1,]
    return(app)
  }
  
  # identify parent with number of aftershocks > 0 
  idx.p <- which(n.ev.v > 0)
  # sample (in parallel) the aftershocks for each parent 
  sample.list <- mclapply(idx.p, function(x) 
    sample.triggered(theta.v, beta.p, n.ev.v[x], Ht[x,], T1, T2, M0, bdy, crsobj),
    mc.cores = ncore)
  # bind the data.frame in the list and return
  sample.pts <- bind_rows(sample.list)    
  sample.pts
}



## function to place n.ev events in space according to a given intensity
# mandatory input : --> loglambda : value of the log-intensity calculated at the mesh nodes
#                 : --> bdy : area in which the points has to be extracted
#                 : --> mesh : mesh of the area 
#                 : --> n.ev : number of events to be placed
#                 : --> crsobj : CRS object

# output : SpatialPoints - set of n.ev points from a point process with generic intensity lambda(s) 

# notes : its a thinning algorithm (we generate many points and retain them with certain probability). 
#         so the algorithm may be slow or provide undesired results if :
#                                 a) wide region with small regions of high intensity
#                                 b) mesh do not represents well changes in the intensity


point_sampler <- function(loglambda, bdy, mesh, n.ev, crsobj = NULL){
  # maximum value of loglambda calculated at the mesh nodes
  loglambda_max <- max(loglambda)
  
  ## Initialize output
  samp.points <- SpatialPoints(data.frame(x = 0, y = 0))
  samp.points <- samp.points[-1,]
  if(!is.null(crsobj)){
    proj4string(samp.points) <- crsobj}
  
  # Initialize number of events
  num <- 0
  ## Until we have retained at least n.ev events
  counter.iter <- 0
  while (num < n.ev){
    ## number of points to sample at a time - might want to adjust depending on retain rate.
    n.points = min(10000*exp(counter.iter), 500000)
    
    # place events at random in bdy
    points <- spsample(bdy, n.points, "random")
    ## transform to wgs84 long/lat
    if(!is.null(crsobj)){pts <- spTransform(points, crsobj)}
    else{pts <- points}
    
    ## calculate the value of loglambda at the extracted points projecting the values of loglambda at the mesh nodes
    proj <- INLA::inla.mesh.project(mesh, pts)
    # ratio between lambda at the extracted point and the max at the mesh nodes 
    lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - loglambda_max)
    # keep with probability lambda_ratio
    keep <- proj$ok & (runif(n.points) <= lambda_ratio)
    kept <- pts[keep]
    # if we retain at least a point
    if(length(kept) > 0){
      samp.points <- rbind(samp.points, kept)
      num <- length(samp.points)
      counter.iter = counter.iter + 1
    }
    else{
      #print('no retained')
    }
  }
  
  #if we retain too many points, we select n.ev events at random
  kp <- sample(seq(1, length(samp.points), by=1), n.ev, replace=FALSE)
  samp.points <- samp.points[kp,]
  
  return(samp.points)
}



# function to sample from an ETAS process
# mandatory inputs : --> theta.v --> parameters vector
#                  : --> beta.p --> parameter of the GR law for the magnitude
#                  : --> n.ev : number of aftershocks to be generated
#                  : --> Ht: known observations (see notes) 
#                  : --> Sigma: 2x2 covariance matrix for space-triggering function
#                  : --> M0 : completeness magnitude
#                  : --> bdy : area in which events have to be placed
#                  : --> crsobj : CRS object to project the locations
#                  : --> bk.field.list : information for spatially varying background (see notes)

# output : list of data.frames : 
#                       each data.frame is a generation of observations from an Hawkes process with intensity
#                       lambda(t,x,m) = pi(m)*(exp(theta1) + exp(theta2)*sum_{ti < t} gm(mi)*gt(t-ti)*gs(s-si))
#                       where pi(m) = beta.p*exp(-beta.p*(m - M0))
#                       columns are: time ("ts"), location ("x","y"), magnitude ("mags"), generation ("gen")
#                                    gen is such that 0 : generated from known events
#                                                     1 : background events
#                                                     2,..j.. : j-th generation, aftershocks of (j-1)-th gen 

# notes: Ht may be observations before T1 as well as observations in (T1, T2)
#      : bk.field.list has to be a list with an element called "loglambda" with the logairthm of the 
#        normalized background intensity field and an element called "mesh" containing the mesh.
#        Notice that loglambda is used only to place events in space, while the number of events is determined
#        by theta.v[1] solely.

sample.ETAS <- function(theta.v, beta.p, T1, T2, M0, bdy, crsobj, Ht = NULL, ncore = 2, 
                        bk.field.list = NULL){
  # if the upper extreme greater than lower
  if(T2 < T1){
    stop('Error - right-end of time interval greater than left-end')
  }
  # background number of events (area calculated in squared km)
  n.bkg <- rpois(1, exp(theta.v[1])*(T2 - T1)*(area(bdy)/1000000))
  #print(n.bkg)
  # if no background events are generated initialize an empty data.frame
  if(n.bkg == 0){
    bkg.df <- data.frame(x = 1, y = 1, ts = 1, mags = 1)
    bkg.df <- bkg.df[-1,]
  }
  else{
    # sample bkg events
    # if no bk.field.list element is passed it assumes uniform background rate
    if(is.null(bk.field.list)){
      bkg.locs <- spsample(bdy, n.bkg, 'random', iter = 10)
      proj4string(bkg.locs) <- crsobj
      bkg.locs <- spTransform(bkg.locs, crsobj)
      bkg.df <- data.frame(x = bkg.locs@coords[,1], 
                           y = bkg.locs@coords[,2], 
                           ts = runif(n.bkg, T1, T2), 
                           mags = rexp(n.bkg, beta.p) + M0, 
                           gen = 1)
    }
    # otherwise it samples using the information provided
    else{
      bkg.locs <- point_sampler(bk.field.list$loglambda, 
                                bdy, bk.field.list$mesh, n.bkg, crsobj)
      
      bkg.df <- data.frame(x = bkg.locs@coords[,1], 
                           y = bkg.locs@coords[,2], 
                           ts = runif(n.bkg, T1, T2), 
                           mags = rexp(n.bkg, beta.p) + M0, 
                           gen = 1)
    }
    }
    
  # if known events are provided
  if(!is.null(Ht)){
    # sample a generation from the known events
    gen.from.past <- sample.generation(theta.v, beta.p, Ht, T1, T2, M0, bdy, crsobj, ncore)
    # if at least an aftershock is produced
    if(nrow(gen.from.past) > 0){
      # set generation
      gen.from.past$gen = 0
      # Merge first generation and background events
      Gen.list <- list(rbind(gen.from.past, bkg.df))  
    }
    else{
      Gen.list <- list(bkg.df)
    }
    
  }
  else{
    Gen.list <- list(bkg.df)
  }
  # stop if we have no background events and no events generated from known observations
  if(nrow(Gen.list[[1]]) == 0){
    #print(exp(theta.v[1])*(T2 - T1)*(area(bdy)/1000000))
    #stop('No events generated - increase theta1')
    return(Gen.list)
  }
  
  # initialize flag and gen counter
  flag = TRUE
  gen = 1
  # this goes until the condition inside the loop is met
  while(flag){
    # set parents
    parents <- Gen.list[[gen]]
    # generate aftershocks
    triggered <- sample.generation(theta.v, beta.p, parents, T1, T2, M0, bdy, crsobj, ncore)
    #print(nrow(triggered))
    # stop the loop if there are no more aftershocks
    if(nrow(triggered) == 0){
      flag = FALSE}
    else{
      # set generations
      triggered$gen = gen + 1
      # store new generation
      Gen.list[[gen + 1]] = triggered
      # update generation counter
      gen = gen + 1
    }
  }
  Gen.list
}


# function to calculate the value of the background field in a given location
# mandatory inputs : --> loc : location at which evaluate the field
#                  : --> Ht : data.frame of observations to determine the field
#                  : --> Sigma.p : 2x2 covariance function for the Gaussian kernel
# optional inputs : --> ncore: number of cores to be used

# output : numeric - value of the background field. For each row in Ht we take a Gaussian kernel with mean
#                    given by the location of the event and covrariance matrix Sigma.p*m where m is the 
#                    magnitude of the event

# notes : Ht has to be a data.frame with at least one column for the magnitudes (called "mags") and 
#         two columns for the locations (called "x", "y").

background.field <- function(loc, Ht, Sigma.p, ncore = 5){
  
  Ht.locs <- cbind(Ht$x, Ht$y)
  Ht.mags <- Ht$mags
  
  n.history <- nrow(Ht)
  aa <- unlist(mclapply(1:nrow(Ht.locs), function(idx) 
    mahalanobis(loc, Ht.locs[idx,], Sigma.p), mc.cores = 5))
  
  idx.xs <- which(aa < 20)
  
  sum(unlist(mclapply(idx.xs, function(idx) 
    dmvnorm(loc, mean = Ht.locs[idx,], sigma = Sigma.p*Ht.mags[idx]), mc.cores = ncore)))
}






# function to generate a sample and extract the number of observations
# mandatory inputs : --> theta.v --> parameters vector
#                  : --> beta.p --> parameter of the GR law for the magnitude
#                  : --> n.ev : number of aftershocks to be generated
#                  : --> Ht: known observations (see notes) 
#                  : --> Sigma: 2x2 covariance matrix for space-triggering function
#                  : --> M0 : completeness magnitude
#                  : --> bdy : area in which events have to be placed
#                  : --> crsobj : CRS object to project the locations

# output : numeric - number of observations in a generated sample.

fore.number <- function(theta.v, beta.p, T1, T2, M0, bdy, crsobj, Ht, ncore = 5, bk.field.list){
  ss <- sample.ETAS(theta.v, beta.p, T1, T2, M0, bdy, crsobj, Ht, ncore, bk.field.list)
  ss1 <- bind_rows(ss)
  sum(ss1$ts >= T1 & ss1$ts < T2)
}




run.fore.experiment <- function(theta.v, beta.p, M0, bdy, crsobj, t.breaks, testing.data,
                                bk.field.list, ncore = 5, nsim = 100){
  
  store.res <- matrix(NA, ncol = 4, nrow = length(t.breaks) - 1)
  print('Completed percentage : ')
  s = Sys.time()
  for(i in 1:(length(t.breaks) - 1)){
    if(sum(ss1$ts < t.breaks[i]) == 0){
      Ht = NULL
    }
    else{
      Ht <- testing.data[testing.data$ts < t.breaks[i],]
    }
    n.sim <- sapply(1:nsim, function(x) 
      fore.number(theta.v, beta.p, T1 = t.breaks[i], T2 = t.breaks[i+1], 
                  M0, bdy, crsobj, Ht = Ht, ncore, bk.field.list))
    
    past.data <- testing.data[testing.data$ts < t.breaks[i+1],]
    
    store.res[i, ] <- c(quantile(n.sim, 0.025), # lower quantile simulated
                        quantile(n.sim, 0.5), # median simulated
                        quantile(n.sim, 0.975), # higher quantile simulated
                        sum(testing.data$ts >= t.breaks[i] & testing.data$ts < t.breaks[i + 1])) # observed 
    
    print(i/(length(t.breaks) - 1))
  } 
  print(Sys.time() - s)
  store.res  
}








