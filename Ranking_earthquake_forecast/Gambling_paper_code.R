library(ggplot2)
library(scales)
library(inlabru)
library(ggpubr)
library(knitr)
library(latex2exp)
library(viridis)
library(kableExtra)
library(raster)
library(rgeos)
library(rgdal)
library(R.utils)
library(bayesianETAS)
library(foreach)
library(doParallel)
library(dplyr)
library(ggmap)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

source('score_utils.R')

####################################################
################ FIGURE 1 ##########################
####################################################


# take a value of pstar
pst = 0.001
# create sequence for plotting
p.seq = seq(0, 1,length.out = 1000)

# brier score differences
brier.s <- brier(pst, p.seq) - brier(pst, pst)
log.s <- log.sc(pst, p.seq) - log.sc(pst, pst)

# normalize brier score with second derivatives ratio
brier.s <- brier.s*(log.score.second.der(pst, pst)/brier.score.second.der(pst))

# create df for plotting 
df = data.frame(p = rep(p.seq, 2),
                Scores = c(brier.s, log.s),
                Score = rep(c('Brier', 'Log'), each = length(p.seq)))

# right plot
p.right <- ggplot(df, aes(x = p, y = Scores, 
                          color = Score, linetype = Score)) + 
  ylab(TeX("$E \\[ \\Delta \\]$"))  +
  geom_hline(yintercept = 0) +
  geom_line() +
  geom_vline(xintercept = pst, linetype = 3) + 
  annotate("text", x = 0.85, y = -2.4, 
           label="(b)", color = "black") +
  theme_classic() + 
  theme(legend.position = c(0.77, 0.8),
        legend.background = element_rect(
          colour ="darkgrey"),
        legend.key.width = unit(0.35,"cm"),
        legend.key.height = unit(0.35,"cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7)) +
  ylim(-3, 0) +
  xlim(0,1)  

# seq for left plot
p.seq = seq(0, 0.006,length.out = 100)
# compute scores
brier.s <- brier(pst, p.seq) - brier(pst, pst)
log.s <- log.sc(pst, p.seq) - log.sc(pst, pst)
# scale brier
brier.s <- brier.s*(log.score.second.der(pst, pst)/brier.score.second.der(pst))

# create df
df = data.frame(p = rep(p.seq, 2),
                Scores = c(brier.s, log.s),
                Score = rep(c('Brier', 'Log'), each = length(p.seq)))

# left plot
p.left <- ggplot(df, aes(x = p, y = Scores, 
                         color = Score, linetype = Score)) + 
  ylab(TeX("$E \\[ \\Delta \\]$")) +
  geom_hline(yintercept = 0) +
  geom_line() +
  geom_vline(xintercept = pst, linetype = 3) +
  annotate("text", x = 0.005, y = -0.004, 
           label="(a)", color = "black") +
  annotate("text", x = pst +0.0006, y = -0.005, 
           label=TeX('p^*'), color = "black", parse = TRUE) +
  theme_classic() + 
  theme(legend.position = c(0.77, 0.8),
        legend.key.width = unit(0.35,"cm"),
        legend.key.height = unit(0.35,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))+
  ylim(-0.005,0) +  
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))
  


pdf('images/figure1.pdf', width = 4, height = 3)
multiplot(p.left, p.right, cols = 2)
dev.off()

####################################################
################ FIGURE 2 ##########################
####################################################

# set pstar
pst <- 0.001
# sequence of values of p1 forecasts
p1.v <- seq(0, 4*pst, length.out = 100)
# values of p2
p2.v <- c(pst, pst*2, pst*4)

p.left <- toplot.gambling(p1.v, pst, p2.v, k = 2) + 
  geom_vline(xintercept = pst) + 
  ylim(-1,1.2) + 
  annotate("text", x = 0.0035, y = -0.3, 
           label="(a)", color = "black") +
  theme_classic() +
  theme(legend.position = c(0.77, 0.8),
        legend.key.width = unit(0.35,"cm"),
        legend.key.height = unit(0.35,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7)) +  
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))

  
p2.v <- c(pst, pst/4, pst/2)
p.right <- toplot.gambling(p1.v, pst, p2.v, k = 2) + 
  geom_vline(xintercept = pst)+
  ylim(-1,1.2) + 
  annotate("text", x = 0.0035, y = -0.3, 
           label="(b)", color = "black") +
  theme_classic() +
  theme(legend.position = c(0.77, 0.8),
        legend.key.width = unit(0.35,"cm"),
        legend.key.height = unit(0.35,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))+  
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))




pdf('images/figure2.pdf', width = 4, height = 3)
multiplot(p.left, p.right, cols = 2)
dev.off()
####################################################
################ FIGURE 3 ##########################
####################################################

# set p3 values 
p3 <- c(pst, pst/3, pst/2)

p.left <- toplot.gambling.diff(p1.v, pst, pst, p3, k = 3) +   
  # geom_vline(xintercept = pmax(2*pst - p3, 0), linetype = 2, color = hue_pal()(3) ) + 
  geom_vline(xintercept = pst)+
  ylab(TeX("$E \\[ \\Delta \\]$"))  +
  ylim(-1,0.2) + 
  annotate("text", x = 0.0005, y = -0.75, 
           label="(a)", color = "black") +
  theme_classic() +
  theme(legend.position = c(0.77, 0.8),
        legend.key.width = unit(0.35,"cm"),
        legend.key.height = unit(0.35,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))+  
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))


# set p3 values
p3 <- c(pst, pst*3, pst*2)
p.right <- toplot.gambling.diff(p1.v, pst, pst, p3, k = 3) +
  # geom_vline(xintercept = pmax(2*pst - p3,0), linetype = 2 , color = hue_pal()(3), 
  #            alpha = 0.75) +
  geom_vline(xintercept = pst)+
  ylab(TeX("$E \\[ \\Delta \\]$"))  +
  ylim(-1,0.2) + 
  annotate("text", x = 0.0005, y = -0.75, 
           label="(b)", color = "black") +
  theme_classic() +
  theme(legend.position = c(0.77, 0.8),
        legend.key.width = unit(0.35,"cm"),
        legend.key.height = unit(0.35,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))+  
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))



pdf('images/figure3.pdf', width = 4, height = 3)
multiplot(p.left, p.right, cols = 2)
dev.off()

####################################################
################ FIGURE 4 ##########################
####################################################

# set value of mean without first two fore
p3 <- pst/2
# different ks
k.v <- c(3,5,10,20)
# set left-boundary of the interval of explored probabilities
p.lim <- (k.v)*pst - ((k.v-2)*p3 + pst)
# set sequence
p1.v <- seq(0, max(p.lim),length.out = 1000)

# 
pl.left <- toplot.gambling.k(p1.v, pst, pst, p3, k.v = k.v) + 
  ylim(-0.002, 0.0015) +
  ylab(TeX("$E \\[ \\Delta \\]$")) +
  theme_classic() +
  theme(legend.position = c(0.77, 0.8),
        legend.key.width = unit(0.35,"cm"),
        legend.key.height = unit(0.35,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))+  
  annotate("text", x = 0.009, y = -0.001, 
           label="(a)", color = "black") +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))

# same 
p3 <- pst*2
k.v <- c(3,5,10,20)
p.lim <- sapply(1:length(k.v), function(x) max(0,(k.v[x])*pst - ((k.v[x]-2)*p3 + pst)))
p1.v <- seq(0, 0.005,length.out = 1000)

pl.right <- toplot.gambling.k(p1.v, pst, pst, p3, k.v = k.v) + 
  ylim(-0.002, 0.0015) +
  theme_classic() + 
  ylab(TeX("$E \\[ \\Delta \\]$")) +
  theme(legend.position = c(0.77, 0.8),
        legend.key.width = unit(0.35,"cm"),
        legend.key.height = unit(0.35,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))+  
  annotate("text", x = 0.004, y = -0.001, 
           label="(b)", color = "black") +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))


pdf('images/figure4.pdf', width = 4, height = 3)
multiplot(pl.left, pl.right, cols = 2)
dev.off()



####################################################
################ FIGURE 5 ##########################
####################################################

# set number of bins
N.bins <- 10000
# set pstar
pstar = 0.001
# expected number of active bins
Xs.exp <- pstar*N.bins
# set forecast p1
p1 = pstar
# set coefficient for forecast p2
alpha = 1/3
p2 = alpha*pstar
# set reference model for pairwise gambling score
p0 = pstar*5

# calculate score differences for different scores
delta.b <- deltas(p1,p2,p0,brier)
delta.l <- deltas(p1,p2,p0,log.sc)
delta.pg <- deltas(p1,p2,p0,'pair.gamb')
delta.fg <- deltas(p1,p2,p0,'full.gamb')

# set values of observed number of active bins
xs = 1:26
# calculate intervals for different scores
res.b <- Xs.extremes(delta.b, N.bins, xs)
res.l <- Xs.extremes(delta.l, N.bins, xs)
res.pg <- Xs.extremes(delta.pg, N.bins, xs)
res.fg <- Xs.extremes(delta.fg, N.bins, xs)

# extract confindence intervals
CI.b <- res.b$df
CI.l <- res.l$df
CI.pg <- res.pg$df
CI.fg <- res.fg$df

# extract extremes
# max
Xs.max.b <-  res.b$extremes[2]
Xs.max.l <-  res.l$extremes[2]
Xs.max.pg <-  res.pg$extremes[2]
if(is.infinite(Xs.max.pg)){Xs.max.pg = NULL}
Xs.max.fg <-  res.fg$extremes[2]
# min
Xs.min.b <-  res.b$extremes[1] 
Xs.min.l <-  res.l$extremes[1]
Xs.min.pg <-  res.pg$extremes[1]
if(is.infinite(Xs.min.pg)){Xs.min.pg = NULL}
Xs.min.fg <-  res.fg$extremes[1]

CIplot.b <- ggplot(CI.b, aes(x = Xs, y = -obs/min(lower))) +
  geom_step(size = 0.3) +
  geom_rect(aes(ymin = -lower/min(lower), 
                ymax = -upper/min(lower),
                xmin = Xs, xmax = lead(Xs)), 
            fill = 'orange', 
            alpha = 0.3, size = 0.1) +
  geom_step(aes(x = Xs, y = -lower/min(lower)), color = 'red',
            size = 0.3) +
  geom_step(aes(x = Xs, y = -upper/min(lower)), color = 'red',
            size = 0.3) +
  geom_hline(yintercept = 0) + 
  xlab('Xs') +
  ylab(TeX('$E \\[ \\Delta \\]$')) + 
  geom_vline(xintercept = Xs.max.b, linetype = 2,
             size = 0.2) +
  geom_vline(xintercept = Xs.min.b, linetype = 2,
             size = 0.2) + 
  annotate("text", x = 20, y = 4.6, label = '(a)', size = 3) + 
  annotate("text", x = 1, y = 4, label = bquote(p[2]),size = 2) +
  annotate("text", x = 7, y = 4, label = 'no pref', size = 2) +
  annotate("text", x = 15, y = 4, label = bquote(p[1]), size = 2) +
  ylim(-1,5.4) + 
  theme_classic() +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))


CIplot.l <- ggplot(CI.l, aes(x = Xs, y = -obs/min(lower))) +
  geom_step(size = 0.3) +
  geom_rect(aes(ymin = -lower/min(lower), 
                ymax = -upper/min(lower),
                xmin = Xs, xmax = lead(Xs)), 
            fill = 'orange', 
            alpha = 0.3) +
  geom_step(aes(x = Xs, y = -lower/min(lower)), color = 'red',
            size = 0.3) +
  geom_step(aes(x = Xs, y = -upper/min(lower)), color = 'red',
            size = 0.3) +
  geom_hline(yintercept = 0) + 
  xlab('Xs') + 
  ylab(TeX('$E \\[ \\Delta \\]$')) + 
  geom_vline(xintercept = Xs.max.l, linetype = 2,
             size = 0.2) +
  geom_vline(xintercept = Xs.min.l, linetype = 2,
             size = 0.2) + 
  annotate("text", x = 20, y = 4.6, label = '(c)', size = 3)+
  annotate("text", x = 0.5, y = 4, label = bquote(p[2]), size = 2) +
  annotate("text", x = 6, y = 4, label = 'no pref', size = 2) +
  annotate("text", x = 15, y = 4, label = bquote(p[1]), size = 2) +
  ylim(-1,5.4) +
  theme_classic() +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))


CIplot.pg <- ggplot(CI.pg, aes(x = Xs, y = -obs/min(lower))) +
  geom_step(size = 0.3) +
  geom_rect(aes(ymin = -lower/min(lower), 
                ymax = -upper/min(lower),
                xmin = Xs, xmax = lead(Xs)), 
            fill = 'orange', 
            alpha = 0.3) +
  geom_step(aes(x = Xs, y = -lower/min(lower)), color = 'red',
            size = 0.3) +
  geom_step(aes(x = Xs, y = -upper/min(lower)), color = 'red',
            size = 0.3) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = Xs.max.pg, linetype = 2,
             size = 0.2) +
  geom_vline(xintercept = Xs.min.pg, linetype = 2,
             size = 0.2) +
  ylab(TeX('$E \\[ \\Delta \\]$')) + 
  xlab('Xs') + 
  annotate("text", x = 2, y = 1.5, label = '(b)', size = 3) +
  annotate("text", x = 3, y = 4, label = bquote(p[2]), size = 2) +
  annotate("text", x = 14, y = 4, label = 'no pref', size = 2) +
  annotate("text", x = 25, y = 4, label = bquote(p[1]), size = 2) +
  ylim(-1,5.4) +
  theme_classic() +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))



CIplot.fg <- ggplot(CI.fg, aes(x = Xs, y = -obs/min(lower))) +
  geom_step(size = 0.3) +
  geom_rect(aes(ymin = -lower/min(lower), 
                ymax = -upper/min(lower),
                xmin = Xs, xmax = lead(Xs)), 
            fill = 'orange', 
            alpha = 0.3) +
  geom_step(aes(x = Xs, y = -lower/min(lower)), color = 'red',
            size = 0.3) +
  geom_step(aes(x = Xs, y = -upper/min(lower)), color = 'red',
            size = 0.3) +
  geom_hline(yintercept = 0) + 
  labs(x = TeX('X_S'),
       y = TeX('$E \\[ \\Delta \\]$')) + 
  geom_vline(xintercept = Xs.max.fg, linetype = 2, size = 0.2) +
  geom_vline(xintercept = Xs.min.fg, linetype = 2, size = 0.2) +
  annotate("text", x = 20, y = 4.6, label = '(d)', size = 3)+
  annotate("text", x = 1, y = 4, label = bquote(p[2]), size = 2) +
  annotate("text", x = 7, y = 4, label = 'no pref', size = 2) +
  annotate("text", x = 15, y = 4, label = bquote(p[1]), size = 2) +
  ylim(-1,5.4) +
  theme_classic() +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))



pdf('images/figure5.pdf', width = 4, height = 3)
multiplot(CIplot.b, CIplot.l, CIplot.pg, CIplot.fg,
          cols = 2)
dev.off()






####################################################
################ TABLE 1 ##########################
####################################################

tt <- data.frame(t(matrix(c(Xs.min.b, Xs.min.l, Xs.min.pg, Xs.min.fg,
                            Xs.max.b, Xs.max.l, Xs.max.pg, Xs.max.fg), byrow = T, ncol = 4)))
colnames(tt) <- c('x~min~', 'x~max~')
rownames(tt) <- c('Brier','Log','PG','FG')
tt

####################################################
################ TABLE 2 ##########################
####################################################

if(is.null(Xs.max.pg)){Xs.max.pg = N.bins}

# calculate probs considering p1 true
prob.p1.b <- 1 - pbinom(Xs.max.b, size = N.bins, prob = pstar) 
prob.p1.l <- 1 - pbinom(Xs.max.l, size = N.bins, prob = pstar) 
prob.p1.pg <- 1 - pbinom(Xs.max.pg, size = N.bins, prob = pstar)
prob.p1.fg <- 1 - pbinom(Xs.max.fg, size = N.bins, prob = pstar) 

prob.p2.b <- pbinom((Xs.min.b-1), size = N.bins, prob = pstar)  
prob.p2.l <- pbinom((Xs.min.l-1), size = N.bins, prob = pstar)  
prob.p2.pg <- pbinom((Xs.min.pg-1), size = N.bins, prob = pstar)  
prob.p2.fg <- pbinom((Xs.min.fg-1), size = N.bins, prob = pstar)  

prob.not.b <- sum(dbinom(Xs.min.b:Xs.max.b, size = N.bins, prob = pstar)) 
prob.not.l <- sum(dbinom(Xs.min.l:Xs.max.l, size = N.bins, prob = pstar))
prob.not.pg <- sum(dbinom(Xs.min.pg:Xs.max.pg, size = N.bins, prob = pstar))
prob.not.fg <- sum(dbinom(Xs.min.fg:Xs.max.fg, size = N.bins, prob = pstar)) 


tt <- data.frame(round(rbind(c(prob.not.b, prob.p1.b,  prob.p2.b),
                             c(prob.not.l, prob.p1.l,  prob.p2.l),
                             c(prob.not.pg, prob.p1.pg,  prob.p2.pg),
                             c(prob.not.fg, prob.p1.fg,  prob.p2.fg)),4))

# same thing but considering p2 true
pstar = p2

prob.p1.b <- 1 - pbinom(Xs.max.b, size = N.bins, prob = pstar) 
prob.p1.l <- 1 - pbinom(Xs.max.l, size = N.bins, prob = pstar) 
prob.p1.pg <- 1 - pbinom(Xs.max.pg, size = N.bins, prob = pstar)
prob.p1.fg <- 1 - pbinom(Xs.max.fg, size = N.bins, prob = pstar) 

prob.p2.b <- pbinom((Xs.min.b-1), size = N.bins, prob = pstar)  
prob.p2.l <- pbinom((Xs.min.l-1), size = N.bins, prob = pstar)  
prob.p2.pg <- pbinom((Xs.min.pg-1), size = N.bins, prob = pstar)  
prob.p2.fg <- pbinom((Xs.min.fg-1), size = N.bins, prob = pstar)  

prob.not.b <- sum(dbinom(Xs.min.b:Xs.max.b, size = N.bins, prob = pstar)) 
prob.not.l <- sum(dbinom(Xs.min.l:Xs.max.l, size = N.bins, prob = pstar))
prob.not.pg <- sum(dbinom(Xs.min.pg:Xs.max.pg, size = N.bins, prob = pstar))
prob.not.fg <- sum(dbinom(Xs.min.fg:Xs.max.fg, size = N.bins, prob = pstar)) 


# merge
tt2 <- data.frame(round(rbind(c(prob.not.b, prob.p1.b,  prob.p2.b),
                              c(prob.not.l, prob.p1.l,  prob.p2.l),
                              c(prob.not.pg, prob.p1.pg,  prob.p2.pg),
                              c(prob.not.fg, prob.p1.fg,  prob.p2.fg)),4))

tt$case = rep('$p^* = p_1$', nrow(tt))

tt2$case = rep('$p^* = p_2$', nrow(tt2))
tt.f <- rbind(tt, tt2)

tt.f <- as.matrix(tt.f)
colnames(tt.f) <- c('No pref', 'p1', 'p2', 'case')
rownames(tt.f) <- rep(c('Brier', 'Log', 'PG', 'FG'),2)
tt.f


####################################################
################ FIGURE 6 ##########################
####################################################

# set p1
p1 = 0.001
# set p2
p2 = p1/3
# set p0 for pairwise gambling
p0 = p1*5
# set sequence of values of pstar
pstar.v = seq(1e-6, 2e-3, length.out = 100)

# calculate the probs for each pstar for different brier and pairwise gambling score
dd.b <- sapply(pstar.v, function(x) Find.probs(N.bins, p1, p2, p0, x, brier))
dd.pg <- sapply(pstar.v, function(x) Find.probs(N.bins, p1, p2, p0, x, 'pair.gamb'))

df.b <- data.frame(probs = c(dd.b[1,], dd.b[2,], dd.b[3,]),
                   Preferences = rep(rownames(dd.b), each = ncol(dd.b)),
                   pstar = rep(pstar.v, 3))

df.pg <- data.frame(probs = c(dd.pg[1,], dd.pg[2,], dd.pg[3,]),
                    Preferences = rep(rownames(dd.pg), each = ncol(dd.pg)),
                    pstar = rep(pstar.v, 3))


prob.plot.b <- ggplot(df.b, aes(x = pstar, y = probs, 
                                linetype = Preferences, color = Preferences)) + 
  geom_line() +   
  geom_vline(xintercept = p2, linetype = 4) +
  geom_vline(xintercept = p1, linetype = 4) +
  #geom_vline(xintercept = (p1 + p2)/2, linetype = 2) +
  annotate('text', x = 0.0018, y = 0.5, label = "(a)") + 
  annotate('text', x = p2 + 0.0001, y = 0.5, label = TeX('p_2'), size = 3) + 
  annotate('text', x = p1 + 0.0001, y = 0.5, label = TeX('p_1'), size = 3) +
  xlab(TeX('p^*')) + 
  ylab('Pr') +
  theme_classic() + 
  theme(legend.background = element_rect(#fill= 'lightgrey',
    linetype="solid", 
    colour ="darkgrey"),
    legend.key.size = unit(1, 'line'),
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7)) 





prob.plot.pg <- ggplot(df.pg, aes(x = pstar, y = probs, 
                                  linetype = Preferences, color = Preferences)) + 
  geom_line() +   
  geom_vline(xintercept = p2, linetype = 4) +
  geom_vline(xintercept = p1, linetype = 4) +
  #geom_vline(xintercept = (p1 + p2)/2, linetype = 2) +
  annotate('text', x = 0.0018, y = 0.5, label = "(b)") + 
  annotate('text', x = p2 + 0.0001, y = 0.5, label = TeX('p_2'), size = 3) + 
  annotate('text', x = p1 + 0.0001, y = 0.5, label = TeX('p_1'), size = 3) +
  xlab(TeX('p^*')) + 
  ylab('Pr') +
  theme_classic() + 
  theme(legend.background = element_rect(#fill= 'lightgrey',
    linetype="solid", 
    colour ="darkgrey"),
    legend.key.size = unit(1, 'line'),
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7))  


pdf('images/figure6.pdf', width = 4, height = 3)
multiplot(prob.plot.b, prob.plot.pg)
dev.off()




####################################################
################ FIGURE 7 ##########################
####################################################

# set pstar vector and number of bins
pstar.v = seq(1e-6, 2e-3, length.out = 100)
N.bins.vec <- c(2000, 5000, 10000, 20000)

library(RColorBrewer)
pal.green <- colorRampPalette(brewer.pal(9, "Greens"))(5)[2:5]

# create plot for brier (ignore warnings)
beta.b <- toplot.beta(p1, p2, p0, pstar.v, N.bins.vec, brier, 'right') +
  annotate('text', x = 0, y = 0.15, label = "(a)") +
  scale_color_manual(values = pal.green) 
  

# finds prob for the scores
db.b <- sapply(pstar.v, 
               function(x) Find.probs(5000, p1, p2, p0, x, brier))
db.l <- sapply(pstar.v, 
               function(x) Find.probs(5000, p1, p2, p0, x, log.sc))
db.pg <- sapply(pstar.v, 
                function(x) Find.probs(5000, p1, p2, p0, x, 'pair.gamb'))
db.fg <- sapply(pstar.v, 
                function(x) Find.probs(5000, p1, p2, p0, x, 'full.gamb'))
beta.v <- c(1 - db.b[1,],
            1 - db.l[1,],
            1 - db.pg[1,],
            1 - db.fg[1,])


# set names
idx.v <- rep(c('Brier', 'Log', 'PG', 'FG'), each = length(pstar.v))

# top plot
pl1 <- ggplot(data.frame(beta = beta.v, 
                         score = idx.v,
                         pstar = rep(pstar.v, 4)), 
              aes(x = pstar, y = beta, linetype = score,, color = score)) + geom_line() + 
  scale_color_brewer(palette="Dark2") + 
  annotate('text', x = 0.00005, y = 0.15, label = "(b)", 
           size = 3) + 
  geom_vline(xintercept = p1, linetype = 3) +
  geom_vline(xintercept = p2, linetype = 3) +
  xlab(TeX('p^*')) + 
  ylab(expression(beta)) + 
  theme_classic() + 
  theme(legend.position = c(0.82, 0.8),
        legend.key.width = unit(0.22,"cm"),
        legend.key.height = unit(0.22,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 4),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))+  
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))



# find probs for different number of bins
db.b <- sapply(pstar.v, 
               function(x) Find.probs(20000, p1, p2, p0, x, brier))
db.l <- sapply(pstar.v, 
               function(x) Find.probs(20000, p1, p2, p0, x, log.sc))
db.pg <- sapply(pstar.v, 
                function(x) Find.probs(20000, p1, p2, p0, x, 'pair.gamb'))
db.fg <- sapply(pstar.v, 
                function(x) Find.probs(20000, p1, p2, p0, x, 'full.gamb'))
beta.v <- c(1 - db.b[1,],
            1 - db.l[1,],
            1 - db.pg[1,],
            1 - db.fg[1,])

pl2 <- ggplot(data.frame(beta = beta.v, 
                         score = idx.v,
                         pstar = rep(pstar.v, 4)), 
              aes(x = pstar, y = beta, 
                  linetype = score, color = score)) + geom_line() + 
  #scale_color_viridis(discrete = T)+
  scale_color_brewer(palette="Dark2") +
  geom_vline(xintercept = p1, linetype = 3) +
  geom_vline(xintercept = p2, linetype = 3) +
  xlab(TeX('p^*')) + 
  annotate('text', x = 0.00005, y = 0.15, label = "(c)",
           size = 3) +
  ylab(expression(beta)) + 
  theme_classic() + 
  theme(legend.position = c(0.82, 0.8),
        legend.key.width = unit(0.22,"cm"),
        legend.key.height = unit(0.22,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 4),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7)) +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))

pdf('images/figure7.pdf', width = 4, height = 3)
multiplot(beta.b, pl1, pl2, layout = matrix(c(1,1,2,3), byrow = T, 
                                                 ncol = 2))

dev.off()




####################################################
################ FIGURE 8 ##########################
####################################################

# set pstar and forecasts p1, p0 (reference model for pairwise gambling score)
p1 = 0.001
pstar <- p1
p0 = 5*p1

# set values of the coefficient determining p2
alphas <- c(seq(1/10, 0.99, length.out = 50),seq(1.01, 4, length.out = 50)) 
p2.v <- alphas*p1

# set number of bins
N.bins.vec <- c(2000, 5000, 10000, 20000)


# find probabilities for the brier score
dd.b <- sapply(p2.v, 
               function(x) Find.probs(N.bins.vec[1], p1, x, p0, pstar, brier))

# store in dataframe
df.b <- data.frame(probs = dd.b[1,],
                   N = N.bins.vec[1],  
                   alpha = alphas)

# repeat for different N.bins
for(i in 2:length(N.bins.vec)){
  dd.b <- sapply(p2.v, 
                 function(x) Find.probs(N.bins.vec[i], p1, x, p0, pstar, brier))
  
  df <- data.frame(probs = dd.b[1,],
                   N = N.bins.vec[i],  
                   alpha = alphas)
  df.b <- rbind(df.b, df)
}

# store
df.b <- df.b[!is.na(df.b$N), ]
df.b <- df.b[order(df.b$N),]
df.b$N <- as.factor(df.b$N)

# create top plot
alpha.beta.pl <- ggplot(df.b, 
                        aes(x = alpha, y = 1 - probs, 
                            linetype = N, color = N)) + 
  geom_line() + 
  geom_vline(xintercept = 1, linetype = 2) +
  labs(color = 'N') + 
  xlab(expression(omega)) + 
  ylab(expression(beta)) +
  theme_classic() + 
  theme(legend.position = c(0.82, 0.68),
        legend.key.width = unit(0.29,"cm"),
        legend.key.height = unit(0.29,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7)) + 
  annotate('text', x = 0, y = 0.85, label = "(a)") +
  scale_color_manual(values = pal.green) 




# repeat using other scoring rules
db.b2 <- sapply(p2.v, 
                function(x) Find.probs(5000, p1, x, p0, pstar, brier))
db.l2 <- sapply(p2.v, 
                function(x) Find.probs(5000, p1, x, p0, pstar, log.sc))
db.pg2 <- sapply(p2.v, 
                 function(x) Find.probs(5000, p1, x, p0, pstar, 'pair.gamb'))
db.fg2 <- sapply(p2.v, 
                 function(x) Find.probs(5000, p1, x, p0, pstar, 'full.gamb'))

# store
beta.v2 <- c(1 - db.b2[1,],
             1 - db.l2[1,],
             1 - db.pg2[1,],
             1 - db.fg2[1,])


# set score names
idx.v2 <- rep(c('Brier', 'Log', 'PG', 'FG'), each = length(p2.v))

# left plot
pl12 <- ggplot(data.frame(beta = beta.v2, 
                          score = idx.v2,
                          p2v = rep(p2.v, 4)), 
               aes(x = p2v/p1, y = beta, color = score, linetype = score)) + geom_line() + 
  scale_color_brewer(palette="Dark2") + 
  geom_vline(xintercept = 1, linetype = 3) +
  xlab(expression(omega)) +
  ylab(expression(beta)) + 
  theme_classic() + 
  ylim(0,1) +
  annotate('text', x = 3.5, y = 0.15, label = "(b)",
           size = 3) +
  theme(legend.position = c(0.82, 0.8),
        legend.key.width = unit(0.22,"cm"),
        legend.key.height = unit(0.22,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 4),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7)) +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))



# repeat everything with N = 20000

db.b2 <- sapply(p2.v, 
                function(x) Find.probs(20000, p1, x, p0, pstar, brier))
db.l2 <- sapply(p2.v, 
                function(x) Find.probs(20000, p1, x, p0, pstar, log.sc))
db.pg2 <- sapply(p2.v, 
                 function(x) Find.probs(20000, p1, x, p0, pstar, 'pair.gamb'))
db.fg2 <- sapply(p2.v, 
                 function(x) Find.probs(20000, p1, x, p0, pstar, 'full.gamb'))

beta.v2 <- c(1 - db.b2[1,],
             1 - db.l2[1,],
             1 - db.pg2[1,],
             1 - db.fg2[1,])



idx.v2 <- rep(c('Brier', 'Log', 'PG', 'FG'), each = length(p2.v))

pl22 <- ggplot(data.frame(beta = beta.v2, 
                          score = idx.v2,
                          p2v = rep(p2.v, 4)), 
               aes(x = p2v/p1, y = beta, color = score, linetype = score)) + geom_line() + 
  scale_color_brewer(palette="Dark2") + 
  geom_vline(xintercept = 1, linetype = 3) +
  xlab(expression(omega)) + 
  ylab(expression(beta)) + 
  theme_classic() +
  ylim(0,1) +
  annotate('text', x = 3.5, y = 0.15, label = "(c)",
           size = 3) +
  theme(legend.position = c(0.82, 0.8),
        legend.key.width = unit(0.22,"cm"),
        legend.key.height = unit(0.22,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 4),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7)) +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))


pdf('images/figure8.pdf', width = 4, height = 3)
multiplot(alpha.beta.pl, pl12, pl22, 
          layout = matrix(c(1,1,2,3), 
                          byrow = T, ncol = 2))

dev.off()




####################################################
################ FIGURE 9 ##########################
####################################################

# import forecast
italy.fore <- read.delim(file = "werner.HiResSmoSeis-m2.italy.5yr.dat",
                         header = FALSE, sep = ' ')

# transform it in numeric and reshape data
italy.f <- as.numeric(italy.fore[1,!is.na(italy.fore[1,])])
for(i in 2:nrow(italy.fore)){
  italy.f <- rbind(italy.f, 
                   as.numeric(italy.fore[i,!is.na(italy.fore[i,])]))
}

# set min magnitude
start = 4.95
col.names <- c('longitude', 'latitude')
# trasform magnitudes in mag class
for(i in 3:ncol(italy.f)){
  end <- round(start + 0.1, 2)
  mag.class <- paste0('[',start,',',end,')')
  col.names <- c(col.names, mag.class)
  start = end
}

# create dataframe storing data
italy.f <- as.data.frame(italy.f)
colnames(italy.f) <- col.names

# initialize probabilities
italy.prob <- italy.f

# calculate probabilities using Poisson distribution
for(i in 3:ncol(italy.f)){
  italy.prob[,i] <-  1 - dpois(0, italy.f[,i])
}


# aggregate 
agg.prob <- as.numeric(apply(italy.prob[,3:ncol(italy.prob)], 1, sum))

# add aggregate probabilities to dataframe and convert it in SpatialPointsDataFrame
df.sp <- cbind(italy.f[,1:2], agg.prob)
coordinates(df.sp) <- ~ longitude + latitude
proj4string(df.sp) <- CRS("+proj=longlat +datum=WGS84")
df.sp <- spTransform(df.sp, CRS("+proj=longlat +datum=WGS84"))
# transform it in SpatialGridDataFrame
gridded(df.sp) <- TRUE

# create boundary box for the costline map
fig.box <- as(extent(df.sp), 'SpatialPolygons')
proj4string(fig.box) <- CRS("+proj=longlat +datum=WGS84")
fig.box <- spTransform(fig.box, CRS("+proj=longlat +datum=WGS84"))
fig.box <- as(fig.box, 'SpatialPolygonsDataFrame')


# take the natural log of the probabilities
df.sp$agg.prob <- log(df.sp$agg.prob)

# take European map
world <- ne_countries(continent = 'Europe', returnclass = "sf", scale = 'medium')

# plot
pdf('images/figure9.pdf', width = 3, height = 3)
ggplot() +  
  gg(df.sp) +# gg(w.shape2) + 
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="lightgrey") + 
  ylim(fig.box@bbox[2,]) + 
  scale_fill_viridis(option = 'magma') + 
  labs(fill = "ln(p)") + 
  theme_classic() +
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7)) + 
  scale_x_continuous(breaks = round(seq(fig.box@bbox[1,1],
                                  fig.box@bbox[1,2],
                                  length.out = 4),1), 
                     limits = fig.box@bbox[1,])
dev.off()



####################################################
################# TABLE 3 ##########################
####################################################

exp.delta.b <- exp.delta(agg.prob, 1/3, 5, brier)
exp.delta.l <- exp.delta(agg.prob, 1/3, 5, log.sc)
exp.delta.pg <- exp.delta(agg.prob, 1/3, 5, 'pair.gamb')
exp.delta.fg <- exp.delta(agg.prob, 1/3, 5, 'full.gamb')


tt3 <- as.matrix(c(Brier = exp.delta.b, 
                   Log = exp.delta.l, 
                   PG = exp.delta.pg, 
                   FG = exp.delta.fg), nrow  = 4)


####################################################
################ FIGURE 10 ##########################
####################################################

# set pstar, forecast p1 and p0 reference model for pairwise gambling score
pstar = agg.prob
p1 = pstar
p0 = pstar*5

# set values of coefficient to calculate p2
alphas <- c(seq(1e-3, 1, length.out = 50),seq(1.01, 4,length.out = 50))

# initialize expected score vector
scored.b <- rep(NA, length(alphas))
scored.l <- rep(NA, length(alphas))
scored.pg <- rep(NA, length(alphas))
scored.fg <- rep(NA, length(alphas))

# iterate over coefficient vector and store
for(i in 1:length(alphas)){
  p2 = alphas[i]*pstar
  
  delta.b <- deltas(p1,p2,p0,brier)
  delta.l <- deltas(p1,p2,p0,log.sc)
  delta.pg <- deltas(p1,p2,p0,'pair.gamb')
  delta.fg <- deltas(p1,p2,p0,'full.gamb')
  
  scored.b[i] <- mean(delta.b[,1]) + mean(pstar*(delta.b[,2] - delta.b[,1]))
  scored.l[i] <- mean(delta.l[,1]) + mean(pstar*(delta.l[,2] - delta.l[,1]))
  scored.pg[i] <- mean(delta.pg[,1]) + mean(pstar*(delta.pg[,2] - delta.pg[,1]))
  scored.fg[i] <- mean(delta.fg[,1]) + mean(pstar*(delta.fg[,2] - delta.fg[,1]))
}

# create a dataframe for plotting
df.p <- rbind(data.frame(alpha = alphas,
                         scored = -scored.b/max(scored.b), score = 'Brier'),
              data.frame(alpha = alphas,
                         scored = -scored.l/max(scored.l), score = 'Log'),
              data.frame(alpha = alphas,
                         scored = -scored.pg/max(scored.pg), score = 'PG'),
              data.frame(alpha = alphas,
                         scored = -scored.fg/max(scored.fg), score = 'FG'))

df.p$score = factor(df.p$score)
# plot


pdf('images/figure10.pdf', width = 4, height = 3)
ggplot(df.p, aes(x = alpha, y = scored, 
                 color = score, linetype = score)) + 
  geom_line() + 
  geom_hline(yintercept = 0) + 
  ylim(-1, 0.3) + 
  geom_vline(xintercept = 1, linetype = 2)+
  scale_color_brewer(palette="Dark2") + 
  ylab(bquote("E["*Delta*"]")) +
  xlab(~omega)  + 
  theme_classic() + 
  theme(legend.position = c(0.82, 0.78),
        legend.key.width = unit(0.35,"cm"),
        legend.key.height = unit(0.35,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7)) 
  
dev.off()



####################################################
################ FIGURE 11 ##########################
####################################################

# set p3
p0 = 5*pstar
# set coefficient for p2
alphas <- c(alphas, seq(4,7,length.out = 100))

# initialize matrix for each score
Brier = matrix(NA, ncol = 3, nrow = length(alphas))
Log = matrix(NA, ncol = 3, nrow = length(alphas))
PG = matrix(NA, ncol = 3, nrow = length(alphas))
FG = matrix(NA, ncol = 3, nrow = length(alphas))

# iterate over coefficient values
for(i in 1:length(alphas)){
  p2 = alphas[i]*pstar
  Brier[i,] = c(p1 = mean(brier(pstar, p1)),
                p2 = mean(brier(pstar, p2)),
                p0 = mean(brier(pstar, p0)))
  
  Log[i,] = c(p1 = mean(log.sc(pstar, p1)),
              p2 = mean(log.sc(pstar, p2)),
              p0 = mean(log.sc(pstar, p0)))
  
  
  PG[i, ] = c(p1 = mean(gamb(pstar, p1, p0)),
              p2 = mean(gamb(pstar, p2, p0)),
              p0 = mean(gamb(pstar, p0, p0)))
  
  FG[i,] = c(p1 = mean(fullgamb(pstar, p1, p2, p0)),
             p2 = mean(fullgamb(pstar, p2, p1, p0)),
             p0 = mean(fullgamb(pstar, p0, p1, p2)))
  
}

# store
Brier.df <- rbind( 
  data.frame(score = Brier[,2] - Brier[,1], 
             alpha = alphas, forecast = 'p2 - p1'),
  data.frame(score = Brier[,2] - Brier[,3], 
             alpha = alphas, forecast = 'p2 - p0'))

Log.df <- rbind( 
  data.frame(score = Log[,2] - Log[,1], 
             alpha = alphas, forecast = 'p2 - p1'),
  data.frame(score = Log[,2] - Log[,3], 
             alpha = alphas, forecast = 'p2 - p0'))

PG.df <- rbind( data.frame(score = PG[,2] - PG[,1], 
                           alpha = alphas, forecast = 'p2 - p1'),
                data.frame(score = PG[,2] - PG[,3], 
                           alpha = alphas, forecast = 'p2 - p0'))

FG.df <- rbind( data.frame(score = FG[,2] - FG[,1], 
                           alpha = alphas, forecast = 'p2 - p1'),
                data.frame(score = FG[,2] - FG[,3], 
                           alpha = alphas, forecast = 'p2 - p0'))


# labels for legend
lab.g <- list(bquote('E['*Delta*'('*bold(p)[2]*','*bold(p)[0]*')]'),
              bquote('E['*Delta*'('*bold(p)[2]*','*bold(p)[1]*')]'))

# legend plot
pl.leg3 <- ggplot(Brier.df, aes(x = alpha, y = -score/min(score), 
                                linetype = forecast, color = forecast)) +
  geom_line() + 
  scale_colour_manual(values = c("#f8766d", "#00b0f6"),
                      labels = lab.g)+
  scale_linetype_manual(values = 1:2,
                        labels =lab.g) + 
  theme(legend.text=element_text(size=18),
        legend.title = element_blank(),
        legend.position = 'bottom',
        legend.background = element_rect(#fill= 'lightgrey',
          linetype = "solid",
          colour = "darkgrey"),
        legend.key.size = unit(3, 'line'))+ 
  guides(linetype = guide_legend(override.aes = list(size = 2)))

# extract legend
leg3 <- as_ggplot(cowplot::get_legend(pl.leg3))


# plot brier
pl.b3 <- ggplot(Brier.df, aes(x = alpha, y = -score/min(score), 
                              linetype = forecast, color = forecast)) +
  geom_line() + 
  geom_hline(yintercept = 0) + 
  ylim(-1, 0.5) + 
  ylab(bquote("E["*Delta*"]")) +
  xlab(~omega) + 
  theme_classic() +
  scale_color_discrete(labels = lab.g) +
  scale_linetype_discrete(labels = lab.g) + 
  annotate("text", x = 5.9, y = 0.4, label = '(a)', size = 3) +
  theme(legend.title = element_text(size = 0),
        legend.position = c(0.78, 0.22),
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.3,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7)) + 
  geom_vline(xintercept = 5, linetype = 3) + 
  geom_vline(xintercept = 1, linetype = 3)


pl.l3 <- ggplot(Log.df, aes(x = alpha, y = -score/min(score), 
                            linetype = forecast, color = forecast)) +
  geom_line()+ 
  geom_hline(yintercept = 0) + 
  ylim(-1, 0.5) + 
  ylab(bquote("E["*Delta*"]")) +
  xlab(~omega) +
  theme_classic() +
  scale_color_discrete(labels = lab.g) +
  scale_linetype_discrete(labels = lab.g) + 
  annotate("text", x = 5.9, y = 0.4, label = '(c)', size = 3) +
  theme(legend.title = element_text(size = 0),
        legend.position = c(0.78, 0.22),
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.3,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7)) + 
  geom_vline(xintercept = 5, linetype = 3)+ 
  geom_vline(xintercept = 1, linetype = 3)



pl.pg3 <- ggplot(PG.df, aes(x = alpha, y = -score/min(score),
                            linetype = forecast, color = forecast)) +
  geom_line() + 
  geom_hline(yintercept = 0) + 
  ylim(-1, 0.5) +
  ylab(bquote("E["*Delta*"]")) +
  xlab(~omega) + 
  theme_classic() +
  scale_color_discrete(labels = lab.g) +
  scale_linetype_discrete(labels = lab.g) + 
  annotate("text", x = 5.9, y = 0.4, label = '(b)', size = 3) +
  theme(legend.title = element_text(size = 0),
        legend.position = c(0.78, 0.22),
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.3,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7)) + 
  geom_vline(xintercept = 5, linetype = 3) + 
  geom_vline(xintercept = 1, linetype = 3)




pl.fg3 <- ggplot(FG.df, aes(x = alpha, y = -score/min(score), 
                            linetype = forecast, color = forecast)) +
  geom_line() + 
  geom_hline(yintercept = 0) + 
  ylim(-1, 0.5) + 
  ylab(bquote("E["*Delta*"]")) +
  xlab(~omega) + 
  theme_classic() +
  scale_color_discrete(labels = lab.g) +
  scale_linetype_discrete(labels = lab.g) + 
  annotate("text", x = 5.9, y = 0.4, label = '(d)', size = 3) +
  theme(legend.title = element_text(size = 0),
        legend.position = c(0.78, 0.22),
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.3,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7)) + 
  geom_vline(xintercept = 5, linetype = 3) + 
  geom_vline(xintercept = 1, linetype = 3)



pdf('images/figure11.pdf', width = 4, height = 3)
multiplot(pl.b3, pl.pg3, pl.l3, pl.fg3, 
          layout = matrix(c(1,2,
                            3,4), byrow = T, ncol = 2))
dev.off()


####################################################
################ FIGURE 12 #########################
####################################################

# set pstar, forecasts p1, p0, number of bins and coefficients vector for p2
pstar = agg.prob
p1 = pstar
p0 = 5*pstar
N = length(pstar)
alphas <- seq(1e-5, 7, length.out = 500)

# set number of simulations
max_loop = 10000
# initialize matrices for the results
Coverage = matrix(NA, ncol = 4, nrow = length(alphas))
Probs.p1 = matrix(NA, ncol = 4, nrow = length(alphas))
Probs.p2 = matrix(NA, ncol = 4, nrow = length(alphas))
Probs.np = matrix(NA, ncol = 4, nrow = length(alphas))

# initialize list storing simulated data
sim.list <- list()

###################
# WARNING NOT RUN # it takes 16 mins results can be loaded later on
###################


# create cluster for parallel computing
cl <- makeCluster(5); registerDoParallel(cl)
# iterate (non in parallel) over coefficients
st.time = Sys.time()
for(i in 1:length(alphas)){
  print(alphas[i])
  # set second forecast
  p2 = alphas[i]*pstar
  
  # find true expected score differences
  truth.b <- mean(brier(pstar, p1) - brier(pstar, p2))
  truth.l <- mean(log.sc(pstar, p1) - log.sc(pstar, p2))
  truth.pg <- mean(gamb(pstar, p1, p0) - gamb(pstar, p2, p0))
  truth.fg <- mean(gamb(pstar, p1, p2) - gamb(pstar, p2, p1))
  
  # run the simulations (in parallel)
  sim.normapprox <- foreach(loop = 1:max_loop, .combine = rbind) %dopar%
    runn.normapprox(pstar,p1,p2,p0)
  
  # storing
  save(sim.normapprox, file =
         paste0('sim.normap.',round(alphas[i],3),'.RData'))
  sim.list[[i]] <- storing(sim.normapprox)
  sim.l <- sim.list[[i]]
  
  ## coverage
  truth.v <- c(truth.b, truth.l, truth.pg, truth.fg)
  Coverage[i,] = sapply(1:4, function(x) mean(truth.v[x] >= sim.l[[x]]$Lower & 
                                                truth.v[x] <= sim.l[[x]]$Upper ))

  # probability of taking a decision
  Probs.p1[i,] <- sapply(1:4, function(x) mean(sim.l[[x]]$Lower > 0))
  Probs.p2[i,] <- sapply(1:4, function(x) mean(sim.l[[x]]$Upper < 0))
  Probs.np[i,] <- 1 - Probs.p1[i,] - Probs.p2[i,]
  
}
print(Sys.time() - st.time)
# stop cluster
stopCluster(cl)
names(sim.list) <- alphas
# save results
save(sim.list, file  = 'simulation.list.RData')
save(Coverage, file = 'coverage.csv')
save(Probs.p1, file = 'Probs.p1.csv')
save(Probs.p2, file = 'Probs.p2.csv')
save(Probs.np, file = 'Probs.np.csv')



# load results
load(file = 'coverage.csv')
load(file = 'Iscore.csv')
load(file = 'Probs.p1.csv')
load(file = 'Probs.p2.csv')
load(file = 'Probs.np.csv')


df.pl <- rbind(data.frame(alpha = alphas,
                          coverage = Coverage[,1],
                          score = 'Brier'),
               data.frame(alpha = alphas,
                          coverage = Coverage[,2],
                          score = 'Log'),
               data.frame(alpha = alphas,
                          coverage = Coverage[,3],
                          score = 'PG'),
               data.frame(alpha = alphas,
                          coverage = Coverage[,4],
                          score = 'FG'))

df.pl$score <- factor(df.pl$score, levels = c('Brier', 'Log',
                                              'FG', 'PG'))

pdf('images/figure12.pdf', width = 4, height = 3)
ggplot(data = df.pl, 
       aes(x = alpha, y = coverage, linetype = score, 
           color = score)) + 
  geom_line() +
  geom_hline(yintercept = 0.95, color = 'darkgrey') +
  scale_color_brewer(palette="Dark2") +
  xlab(~omega) +
  ylab('coverage probability') +
  theme_classic() +
  theme(legend.title = element_text(size = 7),
        legend.position = c(0.78, 0.22),
        legend.key.width = unit(0.3,"cm"),
        legend.key.height = unit(0.3,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 5))
dev.off()


####################################################
################ FIGURE 13 #########################
####################################################

# set p2
p2 = 1e-3*pstar

# find true expected score difference
truth.b01 <- mean(brier(pstar, p1) - brier(pstar, p2))
truth.l01 <- mean(log.sc(pstar, p1) - log.sc(pstar, p2))

#############################################################
# NEED TO RUN ONLY THE LINES ABOUT TRUE EXPECTED SCORE DIFF #
#############################################################

# run the simulation and store the results for gamma = 0.001
sim.na0.001 <- foreach(loop = 1:max_loop, .combine = rbind) %dopar%
   runn.normapprox(pstar,p1,p2,p0)
sim.0.001 <- storing(sim.na0.001)
save(sim.0.001, file = 'sim01.csv')



# find true expected score difference
p2 = 1.5*pstar
truth.b15 <- mean(brier(pstar, p1) - brier(pstar, p2))
truth.l15 <- mean(log.sc(pstar, p1) - log.sc(pstar, p2))

sim.na15 <- foreach(loop = 1:max_loop, .combine = rbind) %dopar%
   runn.normapprox(pstar,p1,p2,p0)

sim.15 <- storing(sim.na15)
save(sim.15, file = 'sim15.csv')

# find true expected score difference
p2 = 4*pstar
truth.b4 <- mean(brier(pstar, p1) - brier(pstar, p2))
truth.l4 <- mean(log.sc(pstar, p1) - log.sc(pstar, p2))

#run the simulation
sim.na4 <- foreach(loop = 1:max_loop, .combine = rbind) %dopar%
  runn.normapprox(pstar,p1,p2,p0)

sim.4 <- storing(sim.na4)
save(sim.4, file = 'sim.4.csv')

# LOAD THE RESULTS
load('sim15.csv')
load('sim01.csv')
load('sim4.csv')

size.let <- 3
pl01 <- toplot.ddensity(sim.0.001$Brier$Obs,
                        truth.b01, c(-3e-5,7e-5)) + theme_classic() +
  annotate("text", x = -2.2e-5, y = 0.9, label = '(a)', 
           size = size.let) +
  annotate("text", x = 5e-5, y = 0.9, label = TeX('$\\omega = 0.001$'), 
           size = size.let - 1, color = 'blue') +
  theme(axis.title = element_text(size = 5),
        axis.text = element_text(size = 4))
  


pl15 <- toplot.ddensity(sim.15$Brier$Obs,
                        truth.b15, c(-4e-5,5e-5)) + theme_classic() + 
  annotate("text", x = -2.7e-5, y = 0.9, label = '(b)', 
           size = size.let) +
  annotate("text", x = 4e-5, y = 0.9, label = TeX('$\\omega = 1.5$'), 
           size = size.let - 1, color = 'blue') +
  theme(axis.title = element_text(size = 5),
        axis.text = element_text(size = 4))

pl4 <- toplot.ddensity(sim.4$Brier$Obs,
                       truth.b4, c(-1.2e-4,1.5e-4)) + theme_classic()+ 
  annotate("text", x = -1e-4, y = 0.9, label = '(c)', 
           size = size.let) + 
  annotate("text", x = 1.3e-4, y = 0.9, label = TeX('$\\omega = 4$'), 
                                    size = size.let - 1, color = 'blue') +
  theme(axis.title = element_text(size = 5),
        axis.text = element_text(size = 4))


pl01.l <- toplot.ddensity(sim.0.001$Log$Obs,
                          truth.l01, c(-1e-3,1e-2)) + theme_classic() + 
  annotate("text", x = 0, y = 0.9, label = '(d)', 
           size = size.let) + 
  annotate("text", x = 0.008, y = 0.9, label = TeX('$\\omega = 0.001$'), 
           size = size.let - 1, color = 'blue') +
  theme(axis.title = element_text(size = 5),
        axis.text = element_text(size = 4))


pl15.l <- toplot.ddensity(sim.15$Log$Obs,
                          truth.l15, c(-1e-3,1e-3)) + theme_classic() + 
  annotate("text", x = -0.7e-3, y = 0.9, label = '(e)',
           size = size.let)+ 
  annotate("text", x = 0.0008, y = 0.9, label = TeX('$\\omega = 1.5$'), 
           size = size.let - 1, color = 'blue') +
  theme(axis.title = element_text(size = 5),
        axis.text = element_text(size = 4))


pl4.l <- toplot.ddensity(sim.4$Log$Obs,
                         truth.l4, c(-1e-3,0.3e-2)) + theme_classic() + 
  annotate("text", x = -0.6e-3, y = 0.9, label = '(f)', 
           size = size.let) + 
  annotate("text", x = 0.0025, y = 0.9, label = TeX('$\\omega = 4$'), 
           size = size.let - 1, color = 'blue') +
  theme(axis.title = element_text(size = 5),
        axis.text = element_text(size = 4))

pdf('images/figure13.pdf', width = 4, height = 3)
multiplot(pl01, pl15, pl4, pl01.l, pl15.l, pl4.l, 
          layout = matrix(1:6, byrow = T, ncol = 3))
dev.off()


####################################################
################ FIGURE 14 #########################
####################################################

# extract the results about probabilities
df.log.pref <- rbind(data.frame(prob = Probs.np[,2], pref = 'no pref',
                                alpha = alphas),
                     data.frame(prob = Probs.p1[,2], pref = 'p1',
                                alpha = alphas),
                     data.frame(prob = Probs.p2[,2], pref = 'p2',
                                alpha = alphas))

df.pg.pref <- rbind(data.frame(prob = Probs.np[,3], pref = 'no pref',
                               alpha = alphas),
                    data.frame(prob = Probs.p1[,3], pref = 'p1',
                               alpha = alphas),
                    data.frame(prob = Probs.p2[,3], pref = 'p2',
                               alpha = alphas))

df.fg.pref <- rbind(data.frame(prob = Probs.np[,4], pref = 'no pref',
                               alpha = alphas),
                    data.frame(prob = Probs.p1[,4], pref = 'p1',
                               alpha = alphas),
                    data.frame(prob = Probs.p2[,4], pref = 'p2',
                               alpha = alphas))


# plot them
pl.prob.l <- ggplot(df.log.pref, aes(x = alpha, y = prob, 
                                     linetype = pref, color = pref)) + 
  annotate("text", x = 5.9, y = 0.5, label = '(a)', size = 3) +
  geom_line() + 
  xlab(~omega)+
  ylab('Pr') +
  labs(color = 'Preferences',
       linetype = 'Preferences') + 
  theme_classic() + 
  theme(#legend.position = c(0.7, 0.7), 
        legend.title = element_text(size = 6),
        legend.key.width = unit(0.3,"cm"),
        legend.key.height = unit(0.3,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))


pl.prob.pg <- ggplot(df.pg.pref, aes(x = alpha, y = prob,
                                     linetype = pref, color = pref)) + 
  geom_line() + 
  annotate("text", x = 5.9, y = 0.5, label = '(b)', size = 3) +
  xlab(~omega)+
  ylab('Pr') +
  labs(color = 'Preferences',
       linetype = 'Preferences') + 
  theme_classic() + 
  theme(#legend.position = c(0.7, 0.7), 
    legend.title = element_text(size = 6),
    legend.key.width = unit(0.3,"cm"),
    legend.key.height = unit(0.3,"cm"),
    legend.background = element_rect(colour ="darkgrey"),
    legend.text = element_text(size = 5),
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7))


pl.prob.fg <- ggplot(df.fg.pref, aes(x = alpha, y = prob,
                                     linetype = pref, color = pref)) + 
  geom_line() + 
  annotate("text", x = 5.9, y = 0.5, label = '(c)', size = 3) +
  xlab(~omega)+
  ylab('Pr') + 
  labs(color = 'Preferences',
       linetype = 'Preferences') +
  theme_classic() + 
  theme(#legend.position = c(0.7, 0.7), 
    legend.title = element_text(size = 6),
    legend.key.width = unit(0.3,"cm"),
    legend.key.height = unit(0.3,"cm"),
    legend.background = element_rect(colour ="darkgrey"),
    legend.text = element_text(size = 5),
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7))


pdf('images/figure14.pdf', width = 4, height = 3)
multiplot(pl.prob.l, pl.prob.pg, pl.prob.fg)
dev.off()

png('images/figure14_red.png', width = 480*1.5, height = 480)
multiplot(pl.prob.l, pl.prob.pg)
dev.off()


####################################################
################ FIGURE 15 #########################
####################################################

# import hauksson catalog
cal.cat <- read.table('~/SEISMICITY/confirmation/houksson_catalog.txt', fill = T)

# select columns to be used
cal.cat <- cal.cat[, 1:11]

# name them
colnames(cal.cat) <- c('yyear', 'mmonth', 'dday', "hhour" ,"mminute", "ssecond", "ID", 
                       "latitude", "longitude" , "depth", "magnitude")

# create a vector of times
time.format <- paste(paste(cal.cat$yyear, cal.cat$mmonth, cal.cat$dday, sep = '-'), paste(cal.cat$hhour, cal.cat$mminute, cal.cat$ssecond, sep = ':'))

# convert it in Date obj
cal.cat$time <- as.POSIXct(time.format, format = "%Y-%m-%d %H:%M:%OS")
cal.cat <- cal.cat[!is.na(cal.cat$time),]

# select only the interval of time needed
cal.10y <- cal.cat[cal.cat$yyear >= 2000 & cal.cat$yyear < 2010,] 
cal.10y <- cal.10y[cal.10y$magnitude > 3.95,]

# time will be expressed in weeks
ss.time <- as.numeric(difftime(cal.10y$time, cal.10y$time[1], unit = 'weeks'))

# find ML estimates
ml.est <- maxLikelihoodETAS(ss.time, cal.10y$magnitude, M0 = 3.95, T = (max(ss.time) + 1),
                            displayOutput = FALSE)

# store them
m.p = ml.est$params[1] 
k.p = ml.est$params[2] 
a.p = ml.est$params[3]
cp = ml.est$params[4]
pp = ml.est$params[5]
beta.p = ml.est$params[6]
Tlim = (max(ss.time) + 1)

##################
# NO NEED TO RUN #
##################

# simulate using the ML estimates as params
# simulate from forecast p1
sim1 <- list()
for(i in 1:1000){
  sim1[[i]] <- new.sim.etas.par(m.p, k.p, a.p, cp, pp, beta.p, M0, Tlim)
}
save(sim1, file = 'sim.final1.RData')

# simulate from forecast p2
sim2 <- list()
for(i in 1:1000){
  sim2[[i]] <- new.sim.etas.par(m.p/3, k.p, a.p, cp, pp, beta.p, M0, Tlim)
}
save(sim2, file = 'sim.final2.RData')

# simulate from forecast p3
sim3 <- list()
for(i in 1:1000){
  sim3[[i]] <- new.sim.etas.par(m.p, k.p/3, a.p, cp, pp, beta.p, M0, Tlim)
}
save(sim3, file = 'sim.final3.RData')

# load the simulated res
load(file = 'sim.final1.RData')
load(file = 'sim.final2.RData')
load(file = 'sim.final3.RData')

# select number of events T < 104 
NN1.104 <- sapply(sim1, function(x) sum(x$ts < 104))
NN2.104 <- sapply(sim2, function(x) sum(x$ts < 104))
NN3.104 <- sapply(sim3, function(x) sum(x$ts < 104))

# select number of events T < 260 
NN1.260 <- sapply(sim1, function(x) sum(x$ts < 260))
NN2.260 <- sapply(sim2, function(x) sum(x$ts < 260))
NN3.260 <- sapply(sim3, function(x) sum(x$ts < 260))

# store them
df.num.104 <- rbind(data.frame(x = NN1.104, model = 'model 1'),
                    data.frame(x = NN2.104, model = 'model 2'),
                    data.frame(x = NN3.104, model = 'model 3'))

df.num.260 <- rbind(data.frame(x = NN1.260, model = 'model 1'),
                    data.frame(x = NN2.260, model = 'model 2'),
                    data.frame(x = NN3.260, model = 'model 3'))

# legend plot
plot.leg.d <- ggplot(df.num.104, aes(x = x, linetype = model, 
                                     color = model)) +
  stat_density(geom = 'line') + 
  labs(color = 'Model', linetype = 'Model') + 
  theme(legend.text=element_text(size=18),
        legend.title = element_blank(),
        legend.position = 'bottom',
        legend.background = element_rect(#fill= 'lightgrey',
          linetype = "solid",
          colour = "darkgrey"),
        legend.key.size = unit(2, 'line'))+ 
  guides(linetype = guide_legend(override.aes = list(size = 1.5)))

# extract legend
leg.d <- as_ggplot(cowplot::get_legend(plot.leg.d))


hist(df.num.104$x[df.num.104$model == 'model 3'], 
     breaks = 300)
df.num.104$model <- as.factor(df.num.104$model)

plot.104 <- ggplot(df.num.104, aes(x = x,
                                   y = ..scaled.. ,
                                   colour = model, 
                                   linetype = model)) +
  stat_density(geom = 'line', position = 'identity') +
  xlim(0,350) + 
  ylim(0,1.05) +
  xlab('Number of events') + 
  ylab('Empirical Density') +
  theme_classic() + 
  annotate("text", x = 300, y = 0.2, label = '(a)', size = 4) +
  theme(legend.position = c(0.8, 0.8), 
    legend.title = element_text(size = 0),
    legend.key.width = unit(0.4,"cm"),
    legend.key.height = unit(0.4,"cm"),
    legend.background = element_rect(colour ="darkgrey"),
    legend.text = element_text(size = 5)) 

plot.260 <- ggplot(df.num.260, aes(x = x, y = ..scaled.., 
                                   color = model, linetype = model)) +
  stat_density(geom = 'line', position = 'identity') +
  xlim(0,350) + 
  ylim(0,1.05) +
  xlab('Number of events') + 
  ylab('Empirical Density') + 
  theme_classic() + 
  annotate("text", x = 300, y = 0.2, label = '(b)', size = 4) +
  theme(legend.position = c(0.8, 0.8), 
        legend.title = element_text(size = 0),
        legend.key.width = unit(0.4,"cm"),
        legend.key.height = unit(0.4,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.text = element_text(size = 5)) 

pdf('images/figure15.pdf', width = 4, height = 3)
multiplot(plot.104, plot.260, cols = 2)
dev.off()

#######################################################
################ FIGURE 15 slides #####################
#######################################################

df.num.104.slide <- df.num.104[df.num.104$model != "model 3",]
df.num.260.slide <- df.num.260[df.num.260$model != "model 3",]

# legend plot
plot.leg.d.s <- ggplot(df.num.104.slide, aes(x = x, linetype = model, color = model)) +
  stat_density(geom = 'line') + 
  labs(color = 'Model', linetype = 'Model') + 
  theme(legend.text=element_text(size=18),
        legend.title = element_blank(),
        legend.position = 'bottom',
        legend.background = element_rect(#fill= 'lightgrey',
          linetype = "solid",
          colour = "darkgrey"),
        legend.key.size = unit(2, 'line'))+ 
  guides(linetype = guide_legend(override.aes = list(size = 1.5)))

# extract legend
leg.d.s <- as_ggplot(cowplot::get_legend(plot.leg.d.s))


plot.104.s <- ggplot(df.num.104.slide, aes(x = x, color = model, linetype = model)) +
  geom_density(size = 1.5)  +
  xlim(0,350) + 
  xlab('Number of events') + 
  ylab('Empirical Density') +
  theme_classic() + 
  theme(legend.position = 'none',
        text = element_text(size = 12))

plot.260.s <- ggplot(df.num.260.slide, aes(x = x, color = model, linetype = model)) +
  geom_density(size = 1.5)  +
  xlim(0,350) + 
  xlab('Number of events') + 
  ylab('Empirical Density') + 
  theme_classic() + 
  theme(legend.position = 'none',
        text = element_text(size = 12))

png('images/figure15_s.png', width = 480*2, height = 480)
multiplot(plot.104.s, plot.260.s, leg.d.s, layout = matrix(c(1,1,2,2,
                                                       1,1,2,2,
                                                       1,1,2,2,
                                                       3,3,3,3), byrow = T, ncol = 4))

dev.off()



#######################################################
################ FIGURE 16-17 #########################
#######################################################

# the code is exactly the same, it changes the TL argument

# set number of bins vector
N.vec <- c(seq(10,200, by = 5), as.integer(seq(219, 500, by = 10)))

TL = 52*2

res12.100 <- mclapply(N.vec, function(x) extract.probs(
  runner.probs(x, sim1, sim2, opt = 1, Tlim = TL)), mc.cores = 5)

res13.100 <- mclapply(N.vec, function(x) extract.probs(
  runner.probs(x, sim1, sim3, opt = 1, Tlim = TL)), mc.cores = 5)

res23.100 <- mclapply(N.vec, function(x) extract.probs(
  runner.probs(x, sim2, sim3, opt = 1, Tlim = TL)), mc.cores = 5)

res21.100 <- mclapply(N.vec, function(x) extract.probs(
  runner.probs(x, sim1, sim2, opt = 2, Tlim = TL)), mc.cores = 5)

res31.100 <- mclapply(N.vec, function(x) extract.probs(
  runner.probs(x, sim1, sim3, opt = 2, Tlim = TL)), mc.cores = 5)

res32.100 <- mclapply(N.vec, function(x) extract.probs(
  runner.probs(x, sim2, sim3, opt = 2, Tlim = TL)), mc.cores = 5)


pp12.100.l <- sapply(res12.100, function(x) x$Log[1])
pp21.100.l <- sapply(res21.100, function(x) x$Log[2])

pp13.100.l <- sapply(res13.100, function(x) x$Log[1])
pp31.100.l <- sapply(res31.100, function(x) x$Log[2])

pp23.100.l <- sapply(res23.100, function(x) x$Log[1])
pp32.100.l <- sapply(res32.100, function(x) x$Log[2])

pp12.100.pg <- sapply(res12.100, function(x) x$PG[1])
pp21.100.pg <- sapply(res21.100, function(x) x$PG[2])

pp13.100.pg <- sapply(res13.100, function(x) x$PG[1])
pp31.100.pg <- sapply(res31.100, function(x) x$PG[2])

pp23.100.pg <- sapply(res23.100, function(x) x$PG[1])
pp32.100.pg <- sapply(res32.100, function(x) x$PG[2])

pp12.100.fg <- sapply(res12.100, function(x) x$FG[1])
pp21.100.fg <- sapply(res21.100, function(x) x$FG[2])

pp13.100.fg <- sapply(res13.100, function(x) x$FG[1])
pp31.100.fg <- sapply(res31.100, function(x) x$FG[2])

pp23.100.fg <- sapply(res23.100, function(x) x$FG[1])
pp32.100.fg <- sapply(res32.100, function(x) x$FG[2])



## MODEL 1 AGAINST MODEL 2

df.prob.1vs2.100 <- rbind(data.frame(N = N.vec, prob.p1 = pp12.100.l, prob.p2 = pp21.100.l,
                                     score = 'Log'),
                          data.frame(N = N.vec, prob.p1 = pp12.100.pg, prob.p2 = pp21.100.pg, 
                                     score = 'PG'),
                          data.frame(N = N.vec, prob.p1 = pp12.100.fg, prob.p2 = pp21.100.fg, 
                                     score = 'FG'))

# model 1 against model 2
df.1vs2.100.l <- rbind(data.frame(x = N.vec, y = pp12.100.l, truth = 'p1'),
                       data.frame(x = N.vec, y = pp21.100.l, truth = 'p2'))

df.1vs2.100.pg <- rbind(data.frame(x = N.vec, y = pp12.100.pg, truth = 'p1'),
                        data.frame(x = N.vec, y = pp21.100.pg, truth = 'p2'))

df.1vs2.100.fg <- rbind(data.frame(x = N.vec, y = pp12.100.fg, truth = 'p1'),
                        data.frame(x = N.vec, y = pp21.100.fg, truth = 'p2'))

# model 1 against model 3
df.1vs3.100.l <- rbind(data.frame(x = N.vec, y = pp13.100.l, truth = 'p1'),
                       data.frame(x = N.vec, y = pp31.100.l, truth = 'p3'))

df.1vs3.100.pg <- rbind(data.frame(x = N.vec, y = pp13.100.pg, truth = 'p1'),
                        data.frame(x = N.vec, y = pp31.100.pg, truth = 'p3'))

df.1vs3.100.fg <- rbind(data.frame(x = N.vec, y = pp13.100.fg, truth = 'p1'),
                        data.frame(x = N.vec, y = pp31.100.fg, truth = 'p3'))

# model 2 against model 3

df.2vs3.100.l <- rbind(data.frame(x = N.vec, y = pp23.100.l, truth = 'p2'),
                       data.frame(x = N.vec, y = pp32.100.l, truth = 'p3'))

df.2vs3.100.pg <- rbind(data.frame(x = N.vec, y = pp23.100.pg, truth = 'p2'),
                        data.frame(x = N.vec, y = pp32.100.pg, truth = 'p3'))

df.2vs3.100.fg <- rbind(data.frame(x = N.vec, y = pp23.100.fg, truth = 'p2'),
                        data.frame(x = N.vec, y = pp32.100.fg, truth = 'p3'))


## producing legends
pl.leg.1vs2 <- ggplot(df.1vs2.100.l, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line() +
  labs(color = 'Data generating model', linetype = 'Data generating model') + 
  theme_classic() +
  theme(#legend.position = c(0.8, 0.8), 
        legend.title = element_text(size = 0),
        legend.key.width = unit(0.4,"cm"),
        legend.key.height = unit(0.4,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.text = element_text(size = 5)) 

leg.1vs2 <- as_ggplot(cowplot::get_legend(pl.leg.1vs2))

pl.leg.1vs3 <- ggplot(df.1vs3.100.l, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line() +
  labs(color = 'Data generating model', linetype = 'Data generating model') + 
  theme_classic() +
  theme(#legend.position = c(0.8, 0.8), 
    legend.title = element_text(size = 0),
    legend.key.width = unit(0.4,"cm"),
    legend.key.height = unit(0.4,"cm"),
    legend.background = element_rect(colour ="darkgrey"),
    legend.text = element_text(size = 5)) 

leg.1vs3 <- as_ggplot(cowplot::get_legend(pl.leg.1vs3))


pl.leg.2vs3 <- ggplot(df.2vs3.100.l, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line() +
  labs(color = 'Data generating model', linetype = 'Data generating model') + 
  theme_classic() +
  theme(#legend.position = c(0.8, 0.8), 
    legend.title = element_text(size = 0),
    legend.key.width = unit(0.4,"cm"),
    legend.key.height = unit(0.4,"cm"),
    legend.background = element_rect(colour ="darkgrey"),
    legend.text = element_text(size = 5)) 

leg.2vs3 <- as_ggplot(cowplot::get_legend(pl.leg.2vs3))


# model 1 against model 2

plot.1vs2.100.l.prob <- 
  ggplot(df.1vs2.100.l, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line() +
  ylim(0,1) +
  annotate("text", x = 450, y = 0.5, label = '(a)', size = 3) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text=element_text(size = 4),
        axis.title = element_text(size = 5)) + 
  xlab('N') + 
  ylab(~beta) + 
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))


df.1vs2.100.pg <- df.1vs2.100.pg[!is.na(df.1vs2.100.pg$y),]
plot.1vs2.100.pg.prob <- 
  ggplot(df.1vs2.100.pg, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line() +
  ylim(0,1) + 
  annotate("text", x = 450, y = 0.5, label = '(b)', size = 3) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text=element_text(size = 4),
        axis.title = element_text(size = 5))+ 
  xlab('N') + 
  ylab(~beta)+ 
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))


plot.1vs2.100.fg.prob <- 
  ggplot(df.1vs2.100.fg, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line() + 
  ylim(0,1) + 
  annotate("text", x = 450, y = 0.5, label = '(c)', size = 3) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text=element_text(size = 4),
        axis.title = element_text(size = 5))+ 
  xlab('N') + 
  ylab(~beta)+ 
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))


# model 1 against model 3

plot.1vs3.100.l.prob <- 
  ggplot(df.1vs3.100.l, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line() + 
  ylim(0,1) + 
  annotate("text", x = 450, y = 0.5, label = '(d)', size = 3) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text=element_text(size = 4),
        axis.title = element_text(size = 5)) + 
  xlab('N') + 
  ylab(~beta) + 
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))


df.1vs3.100.pg <- df.1vs3.100.pg[!is.na(df.1vs3.100.pg$y),]
plot.1vs3.100.pg.prob <- 
  ggplot(df.1vs3.100.pg, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line() +
  ylim(0,1) + 
  annotate("text", x = 450, y = 0.5, label = '(e)', size = 3) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text=element_text(size = 4),
        axis.title = element_text(size = 5))+ 
  xlab('N') + 
  ylab(~beta)+ 
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))


plot.1vs3.100.fg.prob <- 
  ggplot(df.1vs3.100.fg, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line() + 
  ylim(0,1) + 
  annotate("text", x = 450, y = 0.5, label = '(f)', size = 3) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text=element_text(size = 4),
        axis.title = element_text(size = 5))+ 
  xlab('N') + 
  ylab(~beta) + 
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))



# model 2 against model 3

plot.2vs3.100.l.prob <- 
  ggplot(df.2vs3.100.l, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line() + 
  ylim(0,1) + 
  annotate("text", x = 450, y = 0.5, label = '(g)', size = 3) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text=element_text(size = 4),
        axis.title = element_text(size = 5)) + 
  xlab('N') + 
  ylab(~beta)+ 
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))


df.2vs3.100.pg <- df.2vs3.100.pg[!is.na(df.2vs3.100.pg$y),]
plot.2vs3.100.pg.prob <- 
  ggplot(df.2vs3.100.pg, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line() + 
  ylim(0,1) + 
  annotate("text", x = 450, y = 0.5, label = '(h)', size = 3) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text=element_text(size = 4),
        axis.title = element_text(size = 5))+ 
  xlab('N') + 
  ylab(~beta)+ 
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))


plot.2vs3.100.fg.prob <- 
  ggplot(df.2vs3.100.fg, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line() + 
  ylim(0,1) + 
  annotate("text", x = 450, y = 0.5, label = '(i)', size = 3) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text=element_text(size = 4),
        axis.title = element_text(size = 5)) + 
  xlab('N') + 
  ylab(~beta) + 
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) 
  


pdf('images/figure16.pdf', width = 4, height = 4)
multiplot(plot.1vs2.100.l.prob, plot.1vs2.100.pg.prob, plot.1vs2.100.fg.prob, leg.1vs2,
          plot.1vs3.100.l.prob, plot.1vs3.100.pg.prob, plot.1vs3.100.fg.prob, leg.1vs3,
          plot.2vs3.100.l.prob, plot.2vs3.100.pg.prob, plot.2vs3.100.fg.prob, leg.2vs3,
          layout = matrix(c(1,1,2,2,3,3,4,
                            5,5,6,6,7,7,8,
                            9,9,10,10,11,11,12), byrow = T, nrow = 3))
dev.off()

###################################################
############ FIGURE 17 slides #####################
###################################################


# model 1 against model 2

plot.1vs2.200.l.prob <- 
  ggplot(df.1vs2.100.l, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line(size = 1.5) +
  ylim(0,1) +
  theme_classic() + 
  theme(legend.position = 'none') + 
  xlab('N') + 
  ylab('Prob') 


plot.1vs2.200.pg.prob <- 
  ggplot(df.1vs2.100.pg[
    !is.na(df.1vs2.100.pg$y),], aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line(size = 1.5) +
  ylim(0,1) + 
  theme_classic() + 
  theme(legend.position = 'none')+ 
  xlab('N') + 
  ylab('Prob')

plot.1vs2.200.fg.prob <- 
  ggplot(df.1vs2.100.fg, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line(size = 1.5) + 
  ylim(0,1) + 
  theme_classic() + 
  theme(legend.position = 'none')+ 
  xlab('N') + 
  ylab('Prob')

# model 1 against model 3

plot.1vs3.200.l.prob2 <- 
  ggplot(df.1vs3.100.l, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line(size = 1.5) + 
  ylim(0,1) + 
  theme_classic() + 
  theme(legend.position = 'none') + 
  xlab('N') + 
  ylab('Prob') 


plot.1vs3.200.pg.prob2 <- 
  ggplot(df.1vs3.100.pg, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line(size = 1.5) +
  ylim(0,1) + 
  theme_classic() + 
  theme(legend.position = 'none')+ 
  xlab('N') + 
  ylab('Prob')

plot.1vs3.200.fg.prob2 <- 
  ggplot(df.1vs3.100.fg, aes(x = x, y = y, color = truth, linetype = truth)) +
  geom_line(size = 1.5) + 
  ylim(0,1) + 
  theme_classic() + 
  theme(legend.position = 'none')+ 
  xlab('N') + 
  ylab('Prob')



png('images/figure17_slid.png', width = 480*2, height = 480)
multiplot(plot.1vs2.100.l.prob2, plot.1vs2.100.pg.prob2, plot.1vs2.100.fg.prob2, leg.1vs2,
          plot.1vs2.200.l.prob2, plot.1vs2.200.pg.prob, plot.1vs2.200.fg.prob2, leg.1vs2,
          layout = matrix(c(1,1,2,2,3,3,4,
                            5,5,6,6,7,7,8), byrow = T, nrow = 2))
dev.off()




####################################################
################ FIGURE 18 #########################
####################################################

# set Time limit
TL = 52*5
# find probabilities
res13.100 <- mclapply(N.vec, function(x) extract.probs(
  runner.probs(x, sim1, sim3, opt = 1, Tlim = TL)), mc.cores = 5)

# log score
pref.l <- lapply(res13.100, function(x) c(p1 = x$Log[1], p2 = x$Log[2]))  
pref1.l <- sapply(pref.l, function(x) x[1])
pref2.l <- sapply(pref.l, function(x) x[2])

df.l <- rbind(data.frame(N = N.vec, p = pref1.l, pref = 'p1'),
              data.frame(N = N.vec, p = pref2.l, pref = 'p3'),
              data.frame(N = N.vec, p = 1 - pref1.l - pref2.l, pref = 'no pref'))

# full gambling score
pref.fg <- lapply(res13.100, function(x) c(p1 = x$FG[1], p2 = x$FG[2]))  
pref1.fg <- sapply(pref.fg, function(x) x[1])
pref2.fg <- sapply(pref.fg, function(x) x[2])

df.fg <- rbind(data.frame(N = N.vec, p = pref1.fg, pref = 'p1'),
               data.frame(N = N.vec, p = pref2.fg, pref = 'p3'),
               data.frame(N = N.vec, p = 1 - pref1.fg - pref2.fg, pref = 'no pref'))

pl.leg.pref <- ggplot(df.l, aes(x = N, y  = p, color = pref, linetype = pref)) +
  geom_line() +
  labs(color = 'Preferences', linetype = 'Preferences') + 
  theme(legend.text=element_text(size=18),
        legend.title = element_blank(),
        legend.position = 'bottom',
        legend.background = element_rect(#fill= 'lightgrey',
          linetype = "solid",
          colour = "darkgrey"),
        legend.key.size = unit(2, 'line'))+ 
  guides(linetype = guide_legend(override.aes = list(size = 1.5)))

leg.pref <- as_ggplot(cowplot::get_legend(pl.leg.pref))


plot.pref.l <- ggplot(df.l, aes(x = N, y  = p, color = pref, linetype = pref)) + 
  geom_line() + 
  theme_classic() +
  theme(legend.position = c(0.8, 0.6), 
    legend.title = element_text(size = 0),
    legend.key.width = unit(0.4,"cm"),
    legend.key.height = unit(0.4,"cm"),
    legend.background = element_rect(colour ="darkgrey"),
    legend.text = element_text(size = 6)) +
  ylab('Probs')


plot.pref.fg <- ggplot(df.fg, aes(x = N, y  = p, color = pref, linetype = pref)) + 
  geom_line() + 
  theme_classic() +
  theme(legend.position = c(0.8, 0.6), 
        legend.title = element_text(size = 0),
        legend.key.width = unit(0.4,"cm"),
        legend.key.height = unit(0.4,"cm"),
        legend.background = element_rect(colour ="darkgrey"),
        legend.text = element_text(size = 6)) +
  ylab('Probs')

pdf('images/figure18.pdf', width = 4, height = 3)
multiplot(plot.pref.l, plot.pref.fg, cols = 2)

dev.off()




