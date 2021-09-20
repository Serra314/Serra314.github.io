library(inlabru)
library(INLA)
library(ggplot2)
## Import data

# page 56 - 70
data(SPDEtoy)

## create a polygon to define the study region

pl.dom <- cbind(c(0, 1, 1, 0.7, 0), c(0, 0, 0.7, 1, 1))
mesh5 <- inla.mesh.2d(loc.domain = pl.dom, max.e = c(0.092, 0.2))

# function gg() plots automatically various kind of data.
ggplot() + gg(mesh5)

# create an spde object with pc.prior

spde5 <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter
  mesh = mesh5, alpha = 2,
  # P(practic.range < 0.3) = 0.5
  prior.range = c(0.3, 0.5),
  # P(sigma > 1) = 0.01
  prior.sigma = c(10, 0.01)) 


##### we are going to fit a 'simple' model
## where the mean of y is (pag56) is given by beta0 + Au
# with beta0 intercept
# u GF(0, Sigma)
# A projection matrix

#' to do this in INLA we have to explictly build the projection matrix A
#' that is easy (just one line) but then, we have to build the stack
#' that is tricky because we have to REMOVE the intercept
#' and then add it as a fixed effect, this also means that we have to add
#' a covariate (rep(1)) to the dataset.
#' and specify in the function inla.stack that it has not to be projected 
#' imposing A = list(A,1) and the effects in the right order.
#' 
#' With inlabru we don't have to think about it at all!
#' we can just run the model as it is written
#' BUT we have to transform the data into a SpatialPointDataFrame

coordinates(SPDEtoy) = c('s1', 's2')
class(SPDEtoy)

cmp <- y ~ Intercept + gmrf(main = coordinates, model = spde5) 

fit.bru <- bru(components = cmp,
               family = 'gaussian',
               data = SPDEtoy)

fit.bru$summary.fixed
fit.bru$summary.hyperpar


## and very easy to extract posterior distros
p1 <- plot(fit.bru, 'Intercept') + xlab(expression(alpha))
p2 <- plot(fit.bru, 'Range for gmrf') + xlab(expression(rho))
p3 <- plot(fit.bru, 'Stdev for gmrf') + xlab(expression(sigma[epsilon]))
p4 <- plot(fit.bru, 'Precision for the Gaussian observations') + xlab(expression(tau))
# inlabru function to merge different ggplot obj. Quite useful
multiplot(p1, p2, p3, p4, layout = matrix(c(1,2,3,4), byrow = TRUE, ncol = 2))


## no way to get the posterior?
## if we want the posterior distros of GF hyperpar, this also is way much easier
## than INLA that uses functions like emarginal, tmarginal...ecc
## once we have the posteriors we can still use them if we want something
## not straightforward

post.range <- spde.posterior(fit.bru, 'gmrf', 'range')
post.logrange <- spde.posterior(fit.bru, 'gmrf', 'log.range')
post.var <- spde.posterior(fit.bru, 'gmrf', 'variance')
post.logvar <- spde.posterior(fit.bru, 'gmrf', 'log.variance')
post.corr <- spde.posterior(fit.bru, 'gmrf', 'matern.correlation')
post.cov <- spde.posterior(fit.bru, 'gmrf', 'matern.covariance')

# and we can make nice plots in an easy way
multiplot(plot(post.range), plot(post.logrange), plot(post.var),
          plot(post.logvar), plot(post.corr), plot(post.cov),
          layout = matrix(c(1,2,3,4,5,6), byrow = TRUE, ncol = 2))




##----------------- PREDICTION --------------------- ##

#'The point is to predict the expected value in a point that is not been obverved
#' in INLA to do that we have to build another stack relative to the points
#' to be predicted, so we need another A matrix (l111), and to define the 
#' intercept only on those values (l130) and then merge togheter the two
#' stacks (l137). Then, we have to run INLA again.
#' The role played by the option control.mode is not clear and not explained!
#'  
#' After that we have to retrieve those values paying great attention to select
#' the right indexes (l145)
#' 
#' in inlabru instead it is way much easier
#' we need just to define where we want to predict

pts3 <- rbind(c(0.1, 0.1), c(0.5, 0.55), c(0.7, 0.9))

## it doesn't work with a matrix
pts3 <- as.data.frame(pts3); coordinates(pts3) <- c('V1', 'V2')

par(mfrow = c(1,1))
plot(coordinates(SPDEtoy), cex = 0.5)
points(pts3, col = 'red', cex = 0.5, lwd = 2)

pred3 <- predict(fit.bru, pts3, ~ gmrf + Intercept)

## and we don't have to care about indexes or call other function to get the 
## posterior quantity of interest.

pred3 

## ex of function of the field

pred.square <- predict(fit.bru, pts3, ~ gmrf^2)
pred.prod <- predict(fit.bru, pts3, ~gmrf*Intercept)
pred.square
pred.prod

## We can also generate posterior samples in the same way

field.sample <- generate(fit.bru, pts3,
                         ~ gmrf, n.samples = 100)


#### the situation becomes worste if we want to predict the expected value on 
#### regular grid! 
#### first of all, in INLA we have to build a inla.mesh.projector object 
#### then, again, build another stack, (paying attention to the Intercept) and 
#### to the tag used to be able to retrieve what we want, merge together the 
#### two stacks and run again inla.

## in inlabru we just need to build a SpatialPixels and call predict
## To build the SpatialPixels from a mesh there is a specific function

# or we can build it by hand if we are familiar with them
pix <- pixels(mesh5, nx = 101, ny = 101)
class(pix)

pred.grid <- predict(fit.bru, pix, ~ gmrf + Intercept)
pred.grid.gmrf <- predict(fit.bru, pix, ~ gmrf )
# in the same way is possible also to generate samples from the posterior
# also, the output is a SpatialPixelsDataFrame that makes our life way much
# easier.

class(pred.grid)

# also obtaining these plots is not straightforward in INLA

pl_posterior_mean <- ggplot() + 
  gg(pred.grid) + 
  ggtitle("Posterior mean") + 
  coord_fixed()

pl_posterior_sd <- ggplot() + 
  gg(pred.grid, mapping=aes_string(x="x", y="y", fill="sd")) + 
  ggtitle("Posterior sd") + 
  coord_fixed()

gmrf_posterior_mean <- ggplot() + 
  gg(pred.grid.gmrf) + 
  ggtitle("Posterior field mean") + 
  coord_fixed()

gmrf_posterior_sd <- ggplot() + 
  gg(pred.grid.gmrf, mapping=aes_string(x="x", y="y", fill="sd")) + 
  ggtitle("Posterior field sd") + 
  coord_fixed()

multiplot(pl_posterior_mean, gmrf_posterior_mean, 
          pl_posterior_sd, gmrf_posterior_sd, 
          layout = matrix(c(1,2,3,4), byrow = TRUE, ncol = 2))

























