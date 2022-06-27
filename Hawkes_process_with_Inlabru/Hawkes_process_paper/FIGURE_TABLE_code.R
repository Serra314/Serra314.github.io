#################################
## IMPORT FUNCTIONS & PACKAGES ##
#################################
library(readr)
library(viridis)
library(R.filesets)
source('hawkes_functions.R')
source('hawkes_utilities.R')
#################
## IMPORT DATA ##
#################

# set magnitude of completeness and time domain
M0 = 2.99
T1 = 0; T2 = 357

# import data
dd.ama <- read.csv2(file = 'data_M3.0.csv', header = TRUE, sep = ',') %>%
  mutate(time_date = as.POSIXct(paste0(year,'-',month,'-',day,' ',hr,':',min,':',sec)),
         time.diff = as.numeric(difftime(time_date, min(time_date) - 1, units = 'days')),
         Mw = as.numeric(Mw),
         Lon = as.numeric(Lon),
         Lat = as.numeric(Lat),
         Mw.class = cut(Mw, breaks = c(M0, 5, 7))) %>%
  arrange(time_date)


# histogram plot
pl.hist <- ggplot(dd.ama, aes(x = time_date)) + 
  geom_histogram(binwidth = 60*60*24*7, color = 'black', fill = 'white') + 
  geom_star(data = dd.ama[dd.ama$Mw >= 5,], mapping = aes(x = time_date, y = 0), size = 2, fill = 'red') + 
  xlab('time') + 
  ylab('counts') + 
  annotate('text', x = dd.ama$time_date[nrow(dd.ama) - 5], y = 150, label = '(a)') + 
  theme_classic() 

# time vs magnitude scatter plot
pl.mgs <- ggplot(dd.ama[dd.ama$Mw < 5,], aes(x = time_date, y = Mw)) + 
  geom_point(size = 0.2) + 
  geom_star(data = dd.ama[dd.ama$Mw > 5,], mapping = aes(x = time_date, y = Mw), size = 2, fill = 'red') + 
  xlab('time') + 
  ylab('Mw') + 
  annotate('text', x = dd.ama$time_date[nrow(dd.ama) - 5], y = 5, label = '(b)') + 
  theme_classic() 


# calculate cumulative frequencies
CumSum <- sapply(T1:T2, \(x) sum(dd.ama$time.diff < x))
CumSum5 <- sapply(T1:T2, \(x) sum(dd.ama[dd.ama$Mw > 5, ]$time.diff < x))


# plot cumulative
pl.cum <- ggplot() + 
  geom_vline(data = dd.ama[dd.ama$Mw >= 5,], mapping = aes(xintercept = time.diff), 
             linetype = 2, color = 'darkgrey', size = 0.2) + 
  geom_line(aes(x = T1:T2, y = CumSum5*50), color = 'red', linetype = 2) + 
  geom_line(aes(x = T1:T2, y = CumSum)) + 
  geom_star(data = dd.ama[dd.ama$Mw >= 5,], mapping = aes(x = time.diff, y = 0), size = 2, fill = 'red') + 
  xlab('days') + 
  annotate('text', x = dd.ama$time.diff[nrow(dd.ama) - 5], y = 250, label = '(c)') +
  annotate('text', x = dd.ama$time.diff[nrow(dd.ama) - 5], y = 900, label = 'Mw > 3') +
  annotate('text', x = dd.ama$time.diff[nrow(dd.ama) - 5], y = 550, label = 'Mw > 5', color = 'red') +
  theme_classic() + 
  scale_y_continuous("cumulative count", 
                     sec.axis = sec_axis(~ . /50, name = 'cumulative count (Mw > 5)'))


# FIGURE 1
pdf(file = 'figure1.pdf')
multiplot(pl.hist, pl.mgs, pl.cum, layout = matrix(c(1,2,3,3), byrow = TRUE, ncol = 2))
dev.off()

########################
## MCMC MODEL FITTING ##
########################

# 24 mins for 500000
# 16 mins for 300000
#  8 mins for 150000
time.st <- Sys.time()
# generate MCMC posterior samples
MCMC.ama3 <- sampleETASposterior(ts = dd.ama$time.diff, 
                                 magnitudes = dd.ama$Mw, 
                                 M0 = 2.99, 
                                 T=T2, sims = 15000, burnin = 5000, approx = TRUE)
Sys.time() - time.st
# save 
saveRDS(MCMC.ama, file = 'mcmc/MCMC.ama.samples.Rds')
MCMC.ama <- loadRDS('mcmc/MCMC.ama.samples.Rds')

##################################
## MCMC POSTERIOR DISTRIBUTIONS ##
##################################

# parameters name
par.names <- c('mu', 'K', 'alpha', 'c', 'p')

# extrat information for plotting purposes
post.MCMC <- data.frame(value = c(MCMC.ama[,1], MCMC.ama[,2], MCMC.ama[,3], MCMC.ama[,4], MCMC.ama[,5]),
                        param = rep(par.names, each = nrow(MCMC.ama)),
                        type = 'MCMC - posterior')

# data.frame containing priors
prior.MCMC <- rbind(data.frame(x = seq(min(MCMC.ama[,1]), max(MCMC.ama[,1]), 
                                       length.out = 100)) %>%
                      mutate(pdf = dgamma(x, 0.1, 0.1),
                             param = 'mu'),
                    data.frame(x = seq(min(MCMC.ama[,2]), max(MCMC.ama[,2]), 
                                       length.out = 100)) %>%
                      mutate(pdf = dunif(x, 0, 10),
                             param = 'K'),
                    data.frame(x = seq(min(MCMC.ama[,3]), max(MCMC.ama[,3]), 
                                       length.out = 100)) %>%
                      mutate(pdf = dunif(x, 0, 10),
                             param = 'alpha'),
                    data.frame(x = seq(min(MCMC.ama[,4]), max(MCMC.ama[,4]),
                                       length.out = 100)) %>%
                      mutate(pdf = dunif(x, 0, 10),
                             param = 'c'),
                    data.frame(x = seq(min(MCMC.ama[,5]), max(MCMC.ama[,5]), 
                                       length.out = 100)) %>%
                      mutate(pdf = dunif(x, 1, 10),
                             param = 'p')
) %>%
  mutate(type = 'MCMC - prior')



# FIGURE 3
pdf(file = 'figure3.pdf')
ggplot(post.MCMC, aes(x = value)) + geom_density(aes(color = type,
                                                     linetype = type)) + 
  geom_line(data = prior.MCMC, mapping = aes(x, pdf, color = type,
                                             linetype = type)) + 
  facet_wrap(facets = vars(param), scales = 'free',
             labeller = label_parsed) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom')
dev.off()


###########################
## INLABRU MODEL FITTING ##
###########################

# data for Inlabru
data.bru <- data.frame(ts = dd.ama$time.diff, 
                       magnitudes = dd.ama$Mw) %>%
  mutate(idx.p = 1:nrow(dd.ama))

# gamma copula transformation
gamma.t <- function(x, a, b){
  bru_forward_transformation(qgamma, x, a, b)
}
# uniform copula transformation
unif.t <- function(x, a, b){
  bru_forward_transformation(qunif, x, min = a, max = b)
}
# log-gaussian copula transformation
loggaus.t <- function(x, m, s){
  bru_forward_transformation(qlnorm, x, meanlog = m, sdlog = s)
}

##################################
## BINNING STRATEGY SENSITIVITY ##
##################################
# This is done considering the priors of the Inlabru replicate case

# code to find log-gaussian parameters to replicate the bayesianETAS prior for K

## find distribution of Kb as assumed by bayesianETAS priors
k.be = runif(1000000, min = 0, max = 10)
c.be = runif(1000000, min = 0, max = 10)
p.be = runif(1000000, min = 1, max = 10)

# respective
k.bru <- k.be/(c.be*(p.be-1))

# values
m.ln <- -1
sd.ln <- 2.03
# check if quantiles are close to each other
round(rbind(quantile(k.bru, c(0.01, 0.25, 0.5, 0.75, 0.99)),
            qlnorm(c(0.01, 0.25, 0.5, 0.75, 0.99), m.ln, sd.ln)), 4)
# mean and variance of the log-normal distribution
mean.ln(m.ln, sd.ln)
sd.ln(m.ln, sd.ln)

# set link functions for Inlabru replicate 
link.f.be <- list(mu = \(x) gamma.t(x, 0.1, 0.1), 
                  K = \(x) loggaus.t(x, -1, 2.03), 
                  alpha = \(x) unif.t(x, 0, 10), 
                  c_ = \(x) unif.t(x, 0, 10), 
                  p = \(x) unif.t(x, 1, 10))

# set data.frame with binning's parameters values
param.bin_exp <- expand.grid(coef_exp = c(1, 2, 3, 5, 7, 10),
                             N.exp = c(3, 10),
                             delta = c(0.1, 0.2, 0.5))


# set initial values for the parameters in the internal scale to get reasonable initial values (default is 0)
th.init <- list(th.mu = 0.5,
                th.K = 0.5,
                th.alpha = -2,
                th.c = -2,
                th.p = -2) 

# options for inlabru 
bru.opt.list <- list(bru_verbose = 4, # type of visual output 
                     bru_max_iter = 100, # maximum number of iterations
                     #bru_method = list(max_step = 0.5),
                     inla.mode = 'experimental', # type of inla algorithm
                     bru_initial = th.init) # parameters initial values

# fit the model for each binning's parameters combination
for(i in 1:nrow(param.bin_exp)){
  cat('coef.t. = ', param.bin_exp$coef_exp[i],
      'delta.t. = ', param.bin_exp$delta[i], 
      'N.max. = ', param.bin_exp$N.exp[i], '\n')
  # fit the model 
  fit_e <- Hawkes.bru(sample.s = data.bru, # data 
                      M0 = M0, # magnitude of completeness
                      T1 = 0, T2 = T2, # time domain
                      link.functions = link.f.be, # link functions
                      coef.t. = param.bin_exp$coef_exp[i], # binning parameter (delta)
                      delta.t. = param.bin_exp$delta[i], # binning parameter (Delta)
                      N.max. = param.bin_exp$N.exp[i], # binning parameter (n.max)
                      bru.opt = bru.opt.list) # bru options
  # save
  saveRDS(fit_e, file = paste0('fit_binning/',
                               'fit_c',param.bin_exp$coef_exp[i],
                               '_d', param.bin_exp$delta[i],
                               '_N', param.bin_exp$N.exp[i], '.Rds') )
  cat('Completed: ', i/nrow(param.bin_exp), '\n')
}

# load the models
fit_bin_list <- lapply(list.files(path = 'fit_binning/'), \(fore.path)
                       loadRDS(paste0('fit_binning/', fore.path)))

# check convergence - Max deviation from previous: xx% of SD has to be smaller than 1
sapply(fit_bin_list, \(fit) tail(fit$bru_iinla$log, 3) )
# extract information about models (maximum number of iterations, time, if converged)
param.bin_exp$max.iter <- sapply(fit_bin_list, \(fit) max(fit$bru_iinla$track$iteration) )
param.bin_exp$time <- round(sapply(fit_bin_list, \(fit) 
                                   as.numeric(sum(fit$bru_iinla$timings$Time), units = 'mins')), 2)
param.bin_exp$converged <- param.bin_exp$max.iter < 100
# Table 5
param.bin_exp[order(param.bin_exp$time),]


# extract posterior distributions only for delta = 2
post.list <- lapply(which(param.bin_exp$coef_exp == 2), \(idx) 
                    extract.post.df(fit_bin_list[[idx]], link.f.be) %>%
                      mutate(c.exp = param.bin_exp$coef_exp[idx],
                             N.exp = param.bin_exp$N.exp[idx],
                             delta = param.bin_exp$delta[idx]))

# FIGURE 6
pdf('figure6.pdf')
ggplot(bind_rows(post.list) %>% 
         filter(param %in% c('K', 'c', 'p')), 
       aes(x,y,color = factor(delta), linetype = factor(N.exp))) + 
  geom_line() + 
  facet_wrap(facets = vars(param), scales = 'free') + 
  scale_y_log10() + 
  labs(color = ~ Delta, linetype = ~ n[max]) +
  xlab('value') + 
  ylab('log10(pdf)') +
  theme_classic()
dev.off()


# this is the Inlabru replicate model
fit_be <- fit_bin_list[[which(param.bin_exp$delta == 0.1 &
                                param.bin_exp$N.exp == 10 &
                                param.bin_exp$coef_exp == 2
)]]


# extract posterior distribution and set priors
bru.be.post <- extract.post.df(fit_be, link.f.be) %>%
  mutate(prior.type = 'be',
         binning = 'def',
         prior = case_when(param == 'mu' ~ dgamma(x,0.1,0.1),
                           param == 'K' ~ dlnorm(x, meanlog = -1, sdlog = 2.03),
                           param == 'alpha' ~ dunif(x, 0, 10),
                           param == 'c' ~ dunif(x, 0, 10),
                           param == 'p' ~ dunif(x, 1, 10)))



#######################
## PRIOR SENSITIVITY ##
#######################

# set gamma parameters values
v.par <- c(0.1, 0.5, 1, 2, 3)

# set link function list 
link.f.gamma.list <- lapply(1:length(v.par), \(i)  
                            list(mu = \(x) gamma.t(x, force(v.par[i]), 
                                                   10*force(v.par[i])), 
                                 K = \(x) gamma.t(x, 2*force(v.par[i]), 
                                                  force(v.par[i])), 
                                 alpha = \(x) gamma.t(x, 2*force(v.par[i]), 
                                                      force(v.par[i])), 
                                 c_ = \(x) gamma.t(x, force(v.par[i]), 
                                                   10*force(v.par[i])), 
                                 p = \(x) 1 + gamma.t(x, force(v.par[i]), 
                                                      5*force(v.par[i])))
)

# set Inlabru options list
bru.opt.list.gamma <- list(bru_verbose = 3,
                           bru_max_iter = 110,
                           #bru_method = list(max_step = 0.5),
                           inla.mode = 'experimental',
                           bru_initial = list(th.mu = 1,
                                              th.K = 1,
                                              th.alpha = 1,
                                              th.c = 1,
                                              th.p = 1))
# fit the model for each parameter value
for(i in 1:length(v.par)){
  # extract link functions
  link_f <- link.f.gamma.list[[i]]
  # extract gamma value
  gamma.par <- v.par[i]
  prior.name <- paste0('gamma par : ', gamma.par)
  cat(prior.name, '\n')
  # fit the model 
  fit_ <- Hawkes.bru(sample.s = data.bru, 
                     M0 = M0,
                     T1 = 0, T2 = T2, 
                     coef.t. = 2, 
                     delta.t. = 0.1, 
                     N.max. = 10, 
                     link.functions = link_f,
                     bru.opt = bru.opt.list.gamma)
  # save
  saveRDS(fit_, file = paste0('fit_priors/gamma/', 'fit_par', gamma.par, '_N.max10.Rds') )
  cat('Perc compl', i/length(v.par), '\n')
}

# takes file names
l.files <- list.files(path = 'fit_priors/gamma/')

# load model and take parameter values from file name
fit.prior.list <- foreach(i = 1:length(l.files)) %do% {
  fit_ <- loadRDS(paste0('fit_priors/gamma/', l.files[i]))
  g.p <- parse_number(l.files[i])
  list(fit = fit_,
       gamma.par = g.p)
}

# create list of posterior distributions
post.pr.list <- lapply(1:length(fit.prior.list), \(idx) 
                       extract.post.df(fit.prior.list[[idx]]$fit, 
                                       link.f.gamma.list[[idx]]) %>%
                         mutate(gamma.p = fit.prior.list[[idx]]$gamma.par))

# bind rows and set priors
post.pr.bind <- bind_rows(post.pr.list) 
post.pr.bind <- post.pr.bind %>% 
  mutate(prior = case_when(param == 'mu' ~ dgamma(x, gamma.p, 10*gamma.p),
                           param == 'K' ~ dgamma(x, 2*gamma.p, gamma.p),
                           param == 'alpha' ~ dgamma(x, 2*gamma.p, gamma.p),
                           param == 'c' ~ dgamma(x, gamma.p, 10*gamma.p),
                           param == 'p' ~ dgamma(x - 1, gamma.p, 5*gamma.p)))

# FIGURE 7
pdf('figure7.pdf')
ggplot(post.pr.bind, 
       aes(x,y, color = factor(gamma.p))) + 
  geom_line(aes(linetype = 'Posterior')) + 
  geom_line(aes(y = prior, linetype = 'Prior')) + 
  labs(color = ~ gamma, linetype = '') + 
  scale_y_log10() +
  ylab('log10(pdf)') + 
  xlab('value') + 
  facet_wrap(facets = vars(param), scales = 'free',
             labeller = label_parsed) +
  theme_classic() + 
  theme(legend.position = 'bottom') + 
  scale_color_viridis(discrete = TRUE)
dev.off()


########################
## INLABRU GAMMA CASE ##
########################

# set link functions
link.gamma <-   list(mu = \(x) gamma.t(x, 0.1, 1), 
                     K = \(x) gamma.t(x, 1, 0.5), 
                     alpha = \(x) gamma.t(x, 1, 0.5), 
                     c_ = \(x) gamma.t(x, 0.1, 1), 
                     p = \(x) 1 + gamma.t(x, 0.1, 0.5)) 
# fit the model
fit_gamma <- Hawkes.bru(sample.s = data.bru, M0 = M0,
                        T1 = 0, T2 = T2,
                        link.functions = link.gamma,
                        coef.t. = 2, delta.t. = 0.1, N.max. = 10,
                        bru.opt = bru.opt.list.gamma)
saveRDS(fit_gamma, file = 'fit_gamma.Rds')
fit_gamma <- loadRDS('fit_gamma.Rds')

##  Prior and posterior data.frame.
bru.gamma.post <- extract.post.df(fit_gamma, link.gamma) %>%
  mutate(prior.type = 'Inlabru - gamma',
         prior = case_when(param == 'mu' ~ dgamma(x, 0.1, 1),
                           param == 'K' ~ dgamma(x, 1, 0.5),
                           param == 'alpha' ~ dgamma(x, 1, 0.5),
                           param == 'c' ~ dgamma(x, 0.1, 1),
                           param == 'p' ~ dgamma(x - 1, 0.1, 0.5)))

## posterior comparison
bru.be.post <- bru.be.post %>%
  mutate(model = 'Inlabru rep',
         prior.n = 'Prior',
         post.n = 'Posterior')

bru.gamma.post <- bru.gamma.post %>%
  mutate(model = 'Inlabru gamma',
         prior.n = 'Prior',
         post.n = 'Posterior')


bru.post.bind <- bind_rows(bru.be.post, bru.gamma.post) 

# FIGURE 4 & 5
pdf('figure5.pdf')
ggplot(bru.post.bind, aes(x)) + 
  geom_line(aes(y = y, color = model, linetype = post.n)) + 
  geom_line(aes(y = prior, color = model, linetype = prior.n)) +
  scale_y_log10() +
  facet_wrap(facets = vars(param), scales = 'free',
             labeller =  label_parsed) +
  xlab('value') + 
  ylab('log10(pdf)') + 
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom')
dev.off()

######################
## GOOODNESS-OF-FIT ##
######################

# extract 10000 values from the posterior of the parameters
# this is done in 10 batches of 1000 
# Inlabru replicate case
bru.sample.be <- foreach(i = 1:10, .combine = rbind) %do% {
  sample.p <- generate(fit_be, data.frame(1),  ~ data.frame(mu = link.f.be$mu(th.mu),
                                                            K = link.f.be$K(th.K),
                                                            alpha = link.f.be$alpha(th.alpha),
                                                            c = link.f.be$c_(th.c),
                                                            p = link.f.be$p(th.p)), 
                       n.samples = 1000)
  
  bind_rows(sample.p)
}
saveRDS(bru.be.sample, file = 'bru.sample.be.Rds')
bru.be.sample <- loadRDS('bru.sample.be.Rds')
# Inlabru gamma case
bru.sample.gamma <- foreach(i = 1:10, .combine = rbind) %do% {
  sample.p <- generate(fit_gamma, data.frame(1),  ~ data.frame(mu = link.gamma$mu(th.mu),
                                                               K = link.gamma$K(th.K),
                                                               alpha = link.gamma$alpha(th.alpha),
                                                               c = link.gamma$c_(th.c),
                                                               p = link.gamma$p(th.p)), 
                       n.samples = 1000)
  
  bind_rows(sample.p)
}
saveRDS(bru.sample.gamma, file = 'bru.sample.gamma.Rds')
bru.sample.gamma <- loadRDS('bru.sample.gamma.Rds')

# take last 10000 posterior samples
mcmc.sample.p <- tail(MCMC.ama, 10000)

# initialize empty matrices - elements would be Lambda(th)  
# different rows correspond to different posterior samples
# different columns correspond to different observations th
expect.mcmc <- matrix(NA, nrow = 10000, ncol = length(dd.ama$time.diff))
expect.bru <- matrix(NA, nrow = 10000, ncol = length(dd.ama$time.diff))
expect.bru.gamma <- matrix(NA, nrow = 10000, ncol = length(dd.ama$time.diff))

# for each posterior sample
for(i in 1:nrow(bru.sample.rep)){
  print(i/nrow(bru.sample.rep))
  # find value of Lambda(th) for each time
  expect.mcmc[i,] <- sapply(dd.ama$time.diff, \(x) mcmc.sample.p[i,1]*(x) + sum(exp(
    log.Lambda_h(th = mcmc.sample.p[i,], ti = dd.ama$time.diff[dd.ama$time.diff <= x],
                 mi = dd.ama$Mw[dd.ama$time.diff <= x], M0 = M0,T1 = 0, T2 = x)))
  )
  
  bru.p <- as.numeric(bru.sample.rep[i,])
  expect.bru[i,] <- sapply(dd.ama$time.diff, \(x) bru.p[1]*(x) + sum(exp(
    log.Lambda_h2(th = bru.p, ti = dd.ama$time.diff[dd.ama$time.diff <= x],
                  mi = dd.ama$Mw[dd.ama$time.diff <= x], M0 = M0,T1 = 0, T2 = x)))
  )  
  bru.p.gamma <- as.numeric(bru.sample.gamma[i,])
  expect.bru.gamma[i,] <- sapply(dd.ama$time.diff, \(x) bru.p.gamma[1]*(x) + sum(exp(
    log.Lambda_h2(th = bru.p.gamma, ti = dd.ama$time.diff[dd.ama$time.diff <= x],
                  mi = dd.ama$Mw[dd.ama$time.diff <= x], M0 = M0,T1 = 0, T2 = x)))
  )
}
saveRDS(expect.mcmc, file = 'expect.mcmc.Rds')
saveRDS(expect.bru, file = 'expect.bru.Rds')
saveRDS(expect.bru.gamma, file = 'expect.bru.gamma.Rds')

expect.mcmc <- loadRDS('expect.mcmc.Rds')
expect.bru <- loadRDS('expect.bru.Rds')
expect.bru.gamma <- loadRDS('expect.bru.gamma.Rds')

# extract median and quantiles value of Lambda(th) 
df.res <- data.frame(days = rep(dd.ama$time.diff,3),
                     Lambda.med = c(apply(expect.mcmc, 2, median), 
                                    apply(expect.bru, 2, median),
                                    apply(expect.bru.gamma, 2, median)),
                     Lambda.low = c( apply(expect.mcmc, 2, \(x) quantile(x, 0.025)), 
                                     apply(expect.bru, 2, \(x) quantile(x, 0.025)),
                                     apply(expect.bru.gamma, 2, \(x) quantile(x, 0.025))),
                     Lambda.up = c(apply(expect.mcmc, 2, \(x) quantile(x, 0.975)), 
                                   apply(expect.bru, 2, \(x) quantile(x, 0.975)),
                                   apply(expect.bru.gamma, 2, \(x) quantile(x, 0.975))),
                     model = rep(c('BayesianETAS', 'Inlabru - rep', 'Inlabru - gamma'), 
                                 each = nrow(dd.ama)),
                     cumfreq = rep(sapply(dd.ama$time.diff, \(x) sum(dd.ama$time.diff < x)), 3))


# find cumulative frequencies for Lambda(th)
xx <- seq(0,1000,length.out = 500)
cs.mcmc <- matrix(NA, ncol = length(xx), nrow = nrow(expect.mcmc))
cs.bru <-  matrix(NA, ncol = length(xx), nrow = nrow(expect.bru))
cs.bru.gamma <-  matrix(NA, ncol = length(xx), nrow = nrow(expect.bru.gamma))

for(i in 1:nrow(cs.mcmc)){
  cs.mcmc[i,] <- sapply(xx, \(x) sum(expect.mcmc[i,] <= x))
  cs.bru[i,] <- sapply(xx, \(x) sum(expect.bru[i,] <= x))
  cs.bru.gamma[i,] <- sapply(xx, \(x) sum(expect.bru.gamma[i,] <= x))
}


df.cs <- data.frame(Lambda = xx, 
                    cs.med = c( apply(cs.mcmc, 2, median), 
                                apply(cs.bru, 2, median),
                                apply(cs.bru.gamma, 2, median)),
                    cs.low = c( apply(cs.mcmc, 2, \(x) quantile(x, 0.025)), 
                                apply(cs.bru, 2, \(x) quantile(x, 0.025)),
                                apply(cs.bru.gamma, 2, \(x) quantile(x, 0.025))),
                    cs.up = c( apply(cs.mcmc, 2, \(x) quantile(x, 0.975)), 
                               apply(cs.bru, 2, \(x) quantile(x, 0.975)),
                               apply(cs.bru.gamma, 2, \(x) quantile(x, 0.975))),
                    model = rep(c('BayesianETAS', 'Inlabru - rep', 'Inlabru - gamma'), 
                                each = ncol(cs.mcmc))
)

# put together
df.total <- data.frame(x = c( df.res$days, df.cs$Lambda ),
                       cumfreq = c( df.res$cumfreq, rep(NA, nrow(df.cs)) ),
                       xx = c( rep(NA, nrow(df.res)), df.cs$Lambda ),
                       med = c( df.res$Lambda.med, df.cs$cs.med ),
                       low = c( df.res$Lambda.low, df.cs$cs.low ),
                       up = c( df.res$Lambda.up, df.cs$cs.up ),
                       model = c( df.res$model, df.cs$model ),
                       plot = c( rep('days', nrow(df.res)), rep('Lambda', nrow(df.cs)) )
)

# FIGURE 2 - TOP ROW
plot.gof.bru <- ggplot(df.total[df.total$model != 'BayesianETAS',], aes(x = x)) + 
  geom_point(aes(y = cumfreq), size = 1) +
  geom_line(aes(y = xx), linetype = 2) + 
  geom_line(aes(y = med, linetype = model, color = model)) +
  geom_ribbon(aes(ymin = low, ymax = up, 
                  linetype = model, color = model, fill = model), alpha = 0.2) +
  ylab('Cumulative frequencies') + 
  facet_wrap(facets = vars(plot), scales = 'free',
             strip.position = "bottom", 
             labeller = as_labeller(c(days = "days", Lambda = "Lambda") )) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  xlab(NULL) +
  geom_text(data = data.frame(x = c(250,800), y = 250, plot = c('days', 'Lambda'),
                              lab = c('(a)', '(b)')),
            mapping = aes(x, y, label = lab)) + 
  scale_color_manual(values = gg_color_hue(3)[c(3, 1)]) +
  scale_fill_manual(values = gg_color_hue(3)[c(3, 1)]) + 
  #scale_linetype_manual(values = c(1,2)) + 
  theme(legend.title = element_blank())

# FIGURE 2 - BOTTOM ROW
plot.gof.mcmc <- ggplot(df.total[df.total$model != 'Inlabru - gamma',], aes(x = x)) + 
  geom_point(aes(y = cumfreq), size = 1) +
  geom_line(aes(y = xx), linetype = 2) + 
  geom_line(aes(y = med, linetype = model, color = model)) +
  geom_ribbon(aes(ymin = low, ymax = up, 
                  linetype = model, color = model, fill = model), alpha = 0.2) +
  ylab('Cumulative frequencies') + 
  facet_wrap(facets = vars(plot), scales = 'free',
             strip.position = "bottom", 
             labeller = as_labeller(c(days = "days", Lambda = "Lambda") )) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  xlab(NULL) +
  geom_text(data = data.frame(x = c(250,800), y = 250, plot = c('days', 'Lambda'),
                              lab = c('(c)', '(d)')),
            mapping = aes(x, y, label = lab)) + 
  scale_color_manual(values = gg_color_hue(3)[c(2, 1)]) +
  scale_fill_manual(values = gg_color_hue(3)[c(2 ,1)]) + 
  scale_linetype_manual(values = c(1,3)) + 
  theme(legend.title = element_blank())

# FIGURE 2
pdf(file = 'figure2.pdf')
multiplot(plot.gof.bru, plot.gof.mcmc)
dev.off()




