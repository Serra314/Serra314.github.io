# time triggering function
# input: tt = single time at which evaluate the triggering function
#        c.p = c parameter of Omori's law
#        p.p = p parameter of Omori's law
#        H.t = history of the process (vector of times)

trig.T <- function(tt, c.p, p.p, H.t){
  trigg <- (p.p - 1)*(c.p^(p.p - 1))/((tt - H.t + c.p)^(p.p)) 
  trigg[tt <= H.t] = 0
  trigg
}

# okay this is much more stable than the one before
# however it requires as argument p - 1
trig.T2 <- function(tt, c.p, p.pm1, H.t){
  logtrig <- log(p.pm1) + (p.pm1)*log(c.p) - (p.pm1 + 1)*log(tt - H.t + c.p)
  logtrig[tt <= H.t] = 0
  exp(logtrig)
}

# H.t has to be in order
trig.T3 <- function(tt, c.p, p.pm1, H.t, ll0 = -600){
  # compute only for positive time differeces
  diffs <- tt - H.t
  idx.p <- tt - H.t > 0
  logtrig <- log(p.pm1) + (p.pm1)*log(c.p) - (p.pm1 + 1)*log(diffs[idx.p] + c.p)
  

  # approximating for large negative logtrig
  trig = logtrig
  idx <- logtrig < ll0 
  trig[idx] <- exp(ll0)*(1 - (-logtrig[idx] + ll0))
  trig[!idx] <- exp(logtrig[!idx])
  
  toreturn <- rep(0, length(diffs))
  toreturn[idx.p] <- trig
  return(toreturn)
}


# first example of numerical instability considering varying p close to 1
agg <- seq(1e-20, 1e-15, length.out = 100)

aa <- data.frame(pp = agg, 
                 trig1 = trig.T(1.1, 0.001, 1 + agg, 1),
                 trig2 = trig.T2(1.1, 0.001, agg, 1),
                 trig3 = trig.T3(1.1, 0.001, agg, 1))

aa
ggplot(aa, aes(x = pp, y = trig1)) + 
  geom_line() + 
  geom_line(mapping = aes(y = trig2), color = 2, linetype = 2) +
  geom_line(mapping = aes(y = trig3), color = 3, linetype = 3) 

# for large pp
agg <- seq(162, 164, length.out = 100)

aa <- data.frame(pp = agg, 
                 trig1 = trig.T(1.1, 0.001, 1 + agg, 1),
                 trig2 = trig.T2(1.1, 0.001, agg, 1),
                 trig3 = trig.T3(1.1, 0.001, agg, 1))

# print aa to check stability of method 3
ggplot(aa, aes(x = pp, y = trig1)) + 
  geom_line() + 
  geom_line(mapping = aes(y = trig2), color = 2, linetype = 2) +
  geom_line(mapping = aes(y = trig3), color = 3, linetype = 3) 

## checking cp 

cp <- seq(1e40, 1e50, length.out = 100)
p.p = 1.5

aa <- data.frame(cp = cp, 
                 trig1 = trig.T(1.1, cp, p.p, 1),
                 trig2 = trig.T2(1.1, cp, p.p-1, 1),
                 trig3 = trig.T3(1.1, cp, p.p-1, 1))


ggplot(aa, aes(x = cp, y = trig1)) + 
  geom_line() + 
  geom_line(mapping = aes(y = trig2), color = 2, linetype = 2) +
  geom_line(mapping = aes(y = trig3), color = 3, linetype = 3) 

#### okay so we use trig.T3 as official triggering

## now about the integrated intensity

I.h <- function(Tt, c.p, p.p, obs.t){
  1 - (c.p^(p.p-1))*((Tt - obs.t + c.p)^(1 - p.p))
}

I.h2 <- function(Tt, c.p, p.p, obs.t, vv0 = -1e-13){
  vv <- (p.p-1)*log(c.p) + (1-p.p)*log(Tt - obs.t + c.p)
  v2 <- vv
  idx <- vv > vv0
  # we use the approx exp(x) = 1 + x for small x
  v2[idx] <- -v2[idx]
  v2[!idx] <- 1 - exp(v2[!idx])
  v2
}

# first proof of stability
ht = 1
tt = 1 + 1e-6
pp <- seq(3e-10, 3.5e-10, length.out = 100)
aa <- data.frame(p = pp,
                 trig1 = I.h(tt, c.p = 0.01, p.p = 1 + pp, ht),
                 trig2 = I.h2(tt, c.p = 0.01, p.p = 1 + pp, ht))

ggplot(aa, aes(x = p, y = trig1)) + 
  geom_line() + 
  geom_line(aes(y = trig2), color = 2, linetype = 2)


# second proof of stability
ht = 1
tt = 1 + 1e-6
pp <- seq(3e-6, 3.5e-6, length.out = 100)
aa <- data.frame(p = pp,
                 trig1 = I.h(tt, c.p = 1e10, p.p = 1 + pp, ht),
                 trig2 = I.h2(tt, c.p = 1e10, p.p = 1 + pp, ht))

ggplot(aa, aes(x = p, y = trig1)) + 
  geom_line() + 
  geom_line(aes(y = trig2), color = 2, linetype = 2)

## okay we are going to use I.h2 as stable version












# now, a function to sample from the single triggered process

# the function uses inverse sampling, 


sample.omori <- function(n.ev, K, c.p, p.p, obs.t){
    u <- runif(n.ev)
} 


