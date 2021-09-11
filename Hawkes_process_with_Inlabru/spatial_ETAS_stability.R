g.m.unstable <- function(m, theta3, M0){
  if(m - M0 < 0){
    stop('m smaller than completeness threshold')
  }
  log.app <- theta3 + log(m - M0)
  exp(exp(log.app))
}

# stability for large negative values of theta3 
th3 <- seq(4, 6, by = 0.05)
th0 <- 4
df <- rbind(data.frame(x = th3, y = g.m.unstable(9, th3, 2.5), m = 'true'),
            data.frame(x = th3, y = sapply(th3,
                                           function(x) 
                                             g.m(9, x, 2.5, th.p = th0)), m = 'approx'))
# check the differences
cbind(diff(df$y[df$m == 'true']),
      diff(df$y[df$m == 'approx']),
      th3[-1])



# stability for large values of theta3 
th3 <- seq(-32, -29, by = 0.05)
th0 <- -29
df <- rbind(data.frame(x = th3, y = g.m.unstable(2.6, th3, 2.5), m = 'true'),
            data.frame(x = th3, y = sapply(th3,
                                           function(x) 
                                             g.m(2.6, x, 2.5, th.n = th0)), m = 'approx'))

# check the differences
cbind(diff(df$y[df$m == 'true']),
      diff(df$y[df$m == 'approx']),
      th3[-1])



# stability of Omori's law 

gt.unstable <- function(tt, theta4, theta5, Ht){
  time.diff <- tt - Ht
  output <- rep(0, length(time.diff))
  
  idx.p <- time.diff > 0
  
  output[idx.p] <- (-exp(theta5) - 1)*log(time.diff[idx.p] + exp(theta4))
  exp(output)
}


# both parameters highly positive
th4 <- seq(10,15, length.out = 10)
tt <- 1 + 1e-10
ht = 1
th5 = 10

df <- rbind(data.frame(x = th4,
                       y = sapply(th4, function(x)
                         gt.unstable(tt, x, th5, ht)),
                       m = 'true'),
            data.frame(x = th4,
                       y = sapply(th4, function(x)
                         gt(tt, x, th5, ht)),
                       m = 'approx'))

cbind(th4[-1],diff(df$y[df$m == 'true']),
      diff(df$y[df$m == 'approx']))

# both parameters highly negative
th4 <- seq(-65,-55, length.out = 10)
tt <- 1 + 1e-10
ht = 1
th5 = -20

df <- rbind(data.frame(x = th4,
                       y = sapply(th4, function(x)
                         gt.unstable(tt, x, th5, ht)),
                       m = 'true'),
            data.frame(x = th4,
                       y = sapply(th4, function(x)
                         gt(tt, x, th5, ht, th.p = 2)),
                       m = 'approx'))

ggplot(df, aes(x = x, y = y, color = m)) + geom_line()

cbind(th4[-1],diff(df$y[df$m == 'true']),
      diff(df$y[df$m == 'approx']))




