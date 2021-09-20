#setwd("~/Desktop/INLABRU/SPDE_BOOK/inlabru_rep")
library(lattice) 
library(INLA)
library(inlabru)

# page 4 - 36
# function to adjust color scale
colsc <- function(...) {
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                       limits = range(..., na.rm=TRUE))
}
# Import Data
data(SPDEtoy)

SPDEtoy.sp <- SPDEtoy
coordinates(SPDEtoy.sp) <- ~ s1 + s2

# figure 1.1
bubble(SPDEtoy.sp, "y", key.entries = c(5, 7.5, 10, 12.5, 15), 
       maxsize = 2, xlab = "s1", ylab = "s2")

## ------------------------------------------------------------------------
fm = y ~ s1 + s2

## first difference : without specifing a family it gives an Error

fit0 <- bru(fm, family = 'Gaussian', data = SPDEtoy)

# fiure 1.2 easier with Inlabru
p1 <- plot(fit0, 'Intercept') + xlab(expression(alpha))
p2 <- plot(fit0, 's1') + xlab(expression(beta[1]))
p3 <- plot(fit0, 's2') + xlab(expression(beta[2]))
p4 <- plot(fit0, 'Precision for the Gaussian observations') + xlab(expression(tau))
multiplot(p1, p2, p3, p4, layout = matrix(c(1,2,3,4), byrow = TRUE, ncol = 2))

# Summary :
# Much more clean than INLA in which you have to set "by hand" all the graphical 
# parameters, but, we have to specify GAUSSIAN in the bru arguments while INLA
# does it automatically.


## NON-LINEAR EFFECTS OF COVARIATES

# This works pretty similarly and produce similar results

fm1 = y ~ Intercept + 
  smooth1(main = s1, values = sort(unique(SPDEtoy$s1)),  
          model = 'rw1', scale.model = T) +
  smooth2(main = s2, values = sort(unique(SPDEtoy$s2)),  
          model = 'rw1', scale.model = T) 


fit1 <- bru(fm1, family = 'Gaussian', data = SPDEtoy)
summary(fit1)

# easier to extract information, figure 1.3
plot(fit1, 'Intercept')
plot(fit1, 'Precision for the Gaussian observations')
plot(fit1, 'Precision for smooth1')
plot(fit1, 'Precision for smooth2')
plot(fit1, 'smooth1')
plot(fit1, 'smooth2')


########### PREDICTIONS ##########################################################

##' Here the differences are more significant, 'cause inlabru has a specific function
##' to predict, this, in my view, makes things a easier as we don't have to pay
##' attention on setting the right NAs

## INLA instead works providing a dataset with NA where we want to have predictions,
## so in inlabru we can't just add a line on the dataset we have to build up a new one

##' MAIN DIFFERENCES:
##' 
##' - INLA needs to run the entire model with the additional rows to be predicted,
##' 
##' - INLABRU doesn't give detailed information but just summaries


head(SPDEtoy)
# prediction at location (0.5 0.5) figure 1.4
topred0 <- data.frame(s1 = 0.5, s2 = 0.5)
pred0 <- predict(fit0, topred0, y ~ s1 + s2 + Intercept)


## ADDITIONAL ARGUMENTS ##

# They work the same but we need to specify them in a separate input.
# example 

# to use control.compute = list(dic = TRUE)

# we need to specify options = list(control.compute = list(dic = TRUE))


## PRIORS ##

#' They work exactly in the same way

## SEVERAL LIKELIHOODS

#' Here also there are substantial differences between INLA and Inlabru.
#' In INLA we can use two likelihood specifying only one data.frame with NAs in the correct positions
I
# In inlabru instead, we use two different data.frames 1 per likelihood

df1 <- data.frame(Y = SPDEtoy$y,
                  s1 = SPDEtoy$s1,
                  s2 = SPDEtoy$s2)

# create a likelihood

lik1 <- like(formula = Y ~ s1 + s2 + Intercept,
             data = df1,
             family = 'Gaussian')

# create second data.frame and likelihood
df2 <- data.frame(Y = SPDEtoy$y + rnorm(nrow(SPDEtoy), sd = 2),
                  s1 = SPDEtoy$s1,
                  s2 = SPDEtoy$s2)

# create a likelihood
lik2 <- like(formula = Y ~ s1 + s2 + Intercept,
             data = df2,
             family = 'Gaussian')


# fit the model using the likelihoods

fit2 <- bru(lik1, lik2, components = Y ~ s1 + s2)
summary(fit2)



# with didn't explore the other features reported in the book.













