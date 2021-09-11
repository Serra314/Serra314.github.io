source('omori_utils.R')

# first derivative
# mu.v

mu.v <- seq(-6,10,length.out = 1000)

th = c(1e-15, 1e-10)

n.der.list <- lapply(th, function(tr)
  sapply(mu.v, function(x)
    (log.lambda(1, c(exp(x + tr), 1, 0.1, 0.5), 0.5) - 
     log.lambda(1, c(exp(x), 1, 0.1, 0.5), 0.5))/tr))

n.der.2 <- sapply(mu.v, function(x)
  exp(x)/exp(log.lambda(1, c(exp(x), 1, 0.1, 0.5), 0.5)))



# so if INLA uses numeric derivatives computed as I have just done
# then it may troubling. 
df.list <- lapply(n.der.list,
                  function(x) 
                    rbind(data.frame(mu = mu.v, 
                                     der = x,
                                     type = 1),
                          data.frame(mu = mu.v, 
                                     der = n.der.2,
                                     type = 2)))
gg.list <- lapply(1:length(df.list), function(idx) 
  ggplot(df.list[[idx]], 
         aes(x = mu, y = der, color = factor(type),
             linetype = factor(type))) + 
    geom_line() +
    labs(color = 'method',
         linetype = 'method',
         title = paste0('th = ', th[idx])) +
    scale_color_discrete(labels = c('numerical', 
                                    'analytical'))+
    scale_linetype_discrete(labels = c('numerical', 
                                       'analytical')))

multiplot(plotlist = gg.list, cols = 2)

###############################
### with respect to k #########
###############################

k.v <- seq(-5, 5, length.out = 100)

th = c(1e-15, 1e-10)

n.der.list <- lapply(th, function(tr)
  sapply(k.v, function(x)
    (log.lambda(1, c(1, exp(x + tr), 0.1, 0.5), 0.5) - 
       log.lambda(1, c(1, exp(x), 0.1, 0.5), 0.5))/tr))

n.der.2 <- sapply(k.v, function(x)
  theta2.der(1, x, 0.1, 0.5, 0.5)/exp(log.lambda(1, c(1, exp(x), 0.1, 0.5), 0.5)))



# so if INLA uses numeric derivatives computed as I have just done
# then it may troubling. 
df.list <- lapply(n.der.list,
                  function(x) 
                    rbind(data.frame(k = k.v, 
                                     der = x,
                                     type = 1),
                          data.frame(k = k.v, 
                                     der = n.der.2,
                                     type = 2)))

gg.list <- lapply(1:length(df.list), function(idx) 
  ggplot(df.list[[idx]], 
         aes(x = k, y = der, color = factor(type),
             linetype = factor(type))) + 
    geom_line() +
    labs(color = 'method',
         linetype = 'method',
         title = paste0('th = ', th[idx])) +
    scale_color_discrete(labels = c('numerical', 
                                    'analytical'))+
    scale_linetype_discrete(labels = c('numerical', 
                                    'analytical')))

multiplot(plotlist = gg.list, cols = 2)


###############################
### with respect to c #########
###############################

cp.v <- seq(-10, 2, length.out = 100)

th = c(1e-15, 1e-10)

n.der.list <- lapply(th, function(tr)
  sapply(cp.v, function(x)
    (log.lambda(1, c(1, 1, exp(x + tr), 0.5), 0.5) - 
       log.lambda(1, c(1, 1, exp(x), 0.5), 0.5))/tr))

n.der.2 <- sapply(cp.v, function(x)
  theta3.der(1, 1, x, 0.5, 0.5)/exp(log.lambda(1, c(1, 1, exp(x), 0.5), 0.5)))


# so if INLA uses numeric derivatives computed as I have just done
# then it may troubling. 
df.list <- lapply(n.der.list,
                  function(x) 
                    rbind(data.frame(cp = cp.v, 
                                     der = x,
                                     type = 1),
                          data.frame(cp = cp.v, 
                                     der = n.der.2,
                                     type = 2)))

gg.list <- lapply(1:length(df.list), function(idx) 
  ggplot(df.list[[idx]], 
         aes(x = cp, y = der, color = factor(type),
             linetype = factor(type))) + 
    geom_line() +
    labs(color = 'method',
         linetype = 'method',
         title = paste0('th = ', th[idx])) +
    scale_color_discrete(labels = c('numerical', 
                                    'analytical'))+
    scale_linetype_discrete(labels = c('numerical', 
                                       'analytical')))

multiplot(plotlist = gg.list, cols = 2)


###############################
### with respect to p #########
###############################

p.v <- seq(-10, 10, length.out = 100)

th = c(1e-15, 1e-10)

n.der.list <- lapply(th, function(tr)
  sapply(p.v, function(x)
    (log.lambda(1, c(1, 1, 0.1, exp(x + tr)), 0.5) - 
       log.lambda(1, c(1, 1, 0.1, exp(x)), 0.5))/tr))

n.der.2 <- sapply(p.v, function(x)
  theta4.der(1, 1, 0.1, x, 0.5)/exp(log.lambda(1, c(1, 1, 0.1, exp(x)), 0.5)))


# so if INLA uses numeric derivatives computed as I have just done
# then it may troubling. 
df.list <- lapply(n.der.list,
                  function(x) 
                    rbind(data.frame(p = p.v, 
                                     der = x,
                                     type = 1),
                          data.frame(p = p.v, 
                                     der = n.der.2,
                                     type = 2)))

gg.list <- lapply(1:length(df.list), function(idx) 
  ggplot(df.list[[idx]], 
         aes(x = p, y = der, color = factor(type),
             linetype = factor(type))) + 
    geom_line() +
    labs(color = 'method',
         linetype = 'method',
         title = paste0('th = ', th[idx])) +
    scale_color_discrete(labels = c('numerical', 
                                    'analytical'))+
    scale_linetype_discrete(labels = c('numerical', 
                                       'analytical')))

multiplot(plotlist = gg.list, cols = 2)

