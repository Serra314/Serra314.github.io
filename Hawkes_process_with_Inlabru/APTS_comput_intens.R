## ex 2 bootstrap

ss <- rnorm(1000, 2, sd = 3)

bootss <- lapply(1:1000, function(x) sample(ss, 1000, replace = TRUE))
boots.m <- unlist(lapply(bootss, mean))
summary(boots.m)

# 90 % bootstrap intervals
bootCI = quantile(boots.m, c(0.05, 0.95))

# comparison with standard confidence interval
CI = c(mean(ss) - 1.96*sqrt(var(ss)/1000),
       mean(ss) + 1.96*sqrt(var(ss)/1000))

# bootCI is slightly shorter
rbind(bootCI, CI)

## ex 4

three.steps <- function(){
  dice.s <- sample(1:6, 6, replace = T)
  return(sum(dice.s))
}


ss <- sapply(1:1000, function(x) three.steps())
hist(ss, breaks = 20)

thousand.steps <- function(){
  pos.v <- c(0)
  pos <- 0
  for(i in 1:1000){
    pos <- pos + sum(sample(1:6, 2, replace = T))
    if(pos > 40){
      pos <- pos - 40
    }
    pos.v[i + 1] <- pos
  }
  pos.v
}


ss1000 <- thousand.steps()
aa <- unlist(lapply(lapply(1:1000, function(X) thousand.steps()), mean))
hist(aa)



#############################
# practical #################
#############################

# bootstrap

library(boot)

data(Nile)

# 1. and 2. explore Nile data

hist(Nile, breaks = 20)
mean(Nile)
median(Nile)

# 3. Gaussian confidence intervals 
n = length(Nile)
zt <- qt(0.975, 99)
CI.gaus <- c(mean(Nile) - zt*(sd(Nile)/sqrt(n)),
             mean(Nile) + zt*(sd(Nile)/sqrt(n)))

CI.gaus


# 4. bootstrap
R = 10000 # replicates
boot.sample <- lapply(1:R, function(x) sample(Nile, n, replace = TRUE))
boot.means <- unlist(lapply(boot.sample, mean))
boot.medians <- unlist(lapply(boot.sample, median))

# bootstrap CI

# mean

CI.bootperc.mean <- quantile(boot.means, c(0.025, 0.975)) 
CI.bootgaus.mean <- c(mean(boot.means) - zt*(sd(boot.means)),
                      mean(boot.means) + zt*(sd(boot.means)))

rbind(CI.gaus, CI.bootperc.mean,
      CI.bootgaus.mean)
# all very close, bootstrap are slightly smaller


# median
CI.bootperc.med <- quantile(boot.medians, c(0.025, 0.975)) 
CI.bootgaus.med <- c(mean(boot.medians) - zt*(sd(boot.medians)),
                      mean(boot.medians) + zt*(sd(boot.medians)))

rbind(CI.bootperc.med,
      CI.bootgaus.med)

# here we can appreciate the difference, the percentile intervals are smaller
# that is because the bootstrap median distribution is far from normal

rr <- range(boot.medians)
mm <- seq(rr[1], rr[2], length.out = 100)
hist(boot.medians, breaks = 30, freq = FALSE)
lines(mm, dnorm(mm, mean = mean(boot.medians), sd = sd(boot.medians)))
abline(v = median(Nile), col = 'red')
abline(v = CI.bootgaus.med, col = 'blue')
abline(v = CI.bootperc.med, col = 'green')


# instead the mean distribution
rr <- range(boot.means)
mm <- seq(rr[1], rr[2], length.out = 100)
hist(boot.means, breaks = 30, freq = FALSE)
lines(mm, dnorm(mm, mean = mean(boot.means), sd = sd(boot.means)))
abline(v = mean(Nile), col = 'red')
abline(v = CI.bootgaus.mean, col = 'blue')


boot.meadian.CI <- function(R = 10000, method){
  boot.sample <- lapply(1:R, function(x) sample(Nile, n, replace = TRUE))
  
  boot.medians <- unlist(lapply(boot.sample, median))
  
  if(method == 'perc'){
    return(quantile(boot.medians, c(0.025, 0.975)))  
  }
  if(method == 'gauss'){
    return(c(mean(boot.medians) - zt*(sd(boot.medians)),
             mean(boot.medians) + zt*(sd(boot.medians))))
  }
}

# the results are stable with 10000 bootstrap replicates
matrix.median.CI <- matrix(unlist(lapply(1:100, function(x) 
  boot.meadian.CI(R = 10000, method = 'perc'))), byrow = T, ncol = 2)

matrix.median.CI

# not so stable with 10 bootstrap replicates
matrix.median.CI <- matrix(unlist(lapply(1:100, function(x) 
  boot.meadian.CI(R = 10, method = 'perc'))), byrow = T, ncol = 2)

matrix.median.CI


#### CONVERGENCE
# 1. uniform distribution

ecdf.values.unif <- function(n.v){
  ss <- rcauchy(n.v)#runif(n.v)
  F.emp <- ecdf(ss)
  dd <- data.frame(ss, F.emp(ss), n.v)
  dd[order(dd[,1]),]
}

n.v <- c(seq(10, 10000, by = 100), 10000)

ecdf.n.list <- lapply(c(10,100,1000,10000), function(x) ecdf.values.unif(x))
library(dplyr)

uu <- seq(-100,100,length.out = 1000)
ecdf.df <- bind_rows(ecdf.n.list)
df <- rbind(ecdf.df, 
            data.frame(ss = uu,
                       F.emp.ss. = pcauchy(uu),
                       n.v = 'true'))

library(ggplot2)
ggplot(df, aes(x = ss, y = F.emp.ss., color = n.v)) + 
  geom_line() + 
  xlim(-20, 20)


df$diff <- abs(df$F.emp.ss. - pcauchy(df$ss))

mm <- matrix(NA, ncol = 2, nrow = length(unique(df$n.v)))
i = 1
for(xx in unique(df$n.v)){
  df.r <- df %>%
    filter(n.v == xx)
  mm[i,] = c(xx, max(df.r$diff))
  i = i + 1
}

mm

#### 3. transformation method

kk.dens <- function(x, a, b){
  a*b*(x^(a-1))*(1 - x^a)^(b-1)
}

xx <- seq(0,1,length.out = 100)
plot(xx, kk.dens(xx, 2, 2))

kk.cdf <- function(x, a, b){
  1 - (1- x^a)^b
}

plot(xx, kk.cdf(xx, 2, 2))

kk.inv.cdf <- function(u, a, b){
  (1 - (1 - u)^(1/b))^(1/a)
}

kk.sampling <- function(size, a, b){
  u <- runif(size)
  kk.inv.cdf(u, a, b)
}

hist(kk.sampling(1000, 2, 2), breaks = 20, freq = F)
lines(xx, kk.dens(xx, 2, 2), col = 'red')

plot(ecdf(kk.sampling(1000, 2, 2)))
lines(xx, kk.cdf(xx, 2, 2), col = 'red')

timez <- rep(NA, 100)

for(i in 1:100){
  start.t <- Sys.time()
  k.s <- kk.sampling(100000, 2, 2)
  end.t <- Sys.time()
  timez[i] <- difftime(end.t, start.t, units = 'sec')
  }

summary(timez)

## rejection sampling

M = 1.55

plot(xx, M*dunif(xx), type = 'l', xlim = c(0,1), ylim = c(0,M))
lines(xx, kk.dens(xx, 2, 2), col = 'red')

kk.sampling.rej <- function(size, a, b){
  ss <- c()
  while(length(ss) < size){
    x <- runif(1)
    U <- runif(1, min = 0, max = M)
    if(U <= kk.dens(x, a, b)){
      ss <- c(ss, x)
    }
  }
  return(ss)
}


aa <- kk.sampling.rej(10000, 2, 2)
hist(aa, breaks= 20, freq = FALSE)
lines(xx, kk.dens(xx, 2, 2), col = 'red')


timez.rej <- rep(NA, 10)

for(i in 1:10){
  print(i)
  start.t <- Sys.time()
  k.s <- kk.sampling.rej(100000, 2, 2)
  end.t <- Sys.time()
  timez.rej[i] <- difftime(end.t, start.t, units = 'sec')
}

rbind(summary(timez),
      summary(timez.rej))







