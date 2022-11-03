rm(list = ls())
library(deSolve)

## This does a simple monod growth model, in 2d space. 

# this "inner" model would look something like this:
monod = function(time, state, pars){
  # this is the other way to get variables (as opposed to using named lists)
  # i do it this way here because when we move to 2d, thats how we'll get the data
  ecoli = state[1]
  glucose = state[2]
  
  with(as.list(pars),{ # still want to just use named parameters
    decoli = ecoli * vmax * glucose / (glucose + km)
    dglucose = -ecoli * vmax * glucose / (glucose + km)
    
    return(list(c(decoli,dglucose)))
  })
}

# now we put that into a 2d space. We'll make a square, 5cm x 5cm, with 50 boxes per side
fullwidth = 5 # cm
N = 50   # boxes
dx = fullwidth/N# boxwidth (cm)
NN = N * N # total # of boxes

# we'll make the resource diffuse at 0.018 cm2/h (=5e-6cm2/s), and the bacteria
# at 1/50th of that
D_glucose = 0.018   # cm2/h, disffusion const (=5e-6 cm2/s)
D_ecoli = D_glucose / 50

# we supply the spatial variables after the parameters
monod2D <- function (time, state, pars, N, D_ecoli, D_glucose, dx) {
  # this function needs to return a vector of data, not a matrix. That means,
  # each time, the data must be put back into matrix form to do reaction (growth)
  # and diffusion, then put back into vector form. Here we grab the state variables
  # from "state" and separate them into the two different matrices
  NN <- N*N
  ecoli <- matrix(nrow = N, ncol = N,state[1:NN])
  glucose <- matrix(nrow = N, ncol = N,state[(NN+1):(2*NN)])
  # a check, in case of numerical issues. this is good practice for anything which can approach zero, but shouldn't go negative
  glucose[glucose < 0] = 0
  
  with (as.list(pars), {
    ## the "reaction" still happens the same way, since the math is vectorized
    decoli_dt = ecoli * vmax * glucose / (glucose + km)
    dglucose_dt = -ecoli * vmax * glucose / (glucose + km)
    
    # "zero" is a helper deSolve suggested to make calculating diffusion fluxes easier
    zero = rep(0, N)
    
    # next set of equations does the diffusion calculations
    
    # if you imagine data like this:   c(1,3,2,4). Focus on the 3. It changes due to 
    # diffusion in one axis by adding the data on both sides, then subtracting 2x of the focal.
    # so that "3" changes like:  +1 + 2 -3 -3.  Then, this value is made spatially accurate
    # by multiplying by D/(dx^2)
    #
    # in practice its slightly trickier, because the edges should reflect. We deal
    # with this by adding a copied row above and below before performing that calculation.
    ecoli_rows_bigger = rbind(ecoli[1,], ecoli, ecoli[N,])
    ecoli_rows_diff = ecoli_rows_bigger[1:N,] + ecoli_rows_bigger[3:(N+2),] - 2 * ecoli
    ecoli_rows_diff = D_ecoli * (1 / dx^2) * ecoli_rows_diff
    
    glucose_rows_bigger = rbind(glucose[1,], glucose, glucose[N,])
    glucose_rows_diff = glucose_rows_bigger[1:N,] + glucose_rows_bigger[3:(N+2),] - 2 * glucose
    glucose_rows_diff = D_glucose * (1 / dx^2) * glucose_rows_diff
    
    ## Add flux gradient across rows to rate of change, aka add the diffusion to the reaction
    decoli_dt    = decoli_dt + ecoli_rows_diff
    dglucose_dt  = dglucose_dt + glucose_rows_diff
    
    ## 2. Now repeat, but across columns
    ecoli_cols_bigger = cbind(ecoli[,1], ecoli, ecoli[,N])
    ecoli_cols_diff = ecoli_cols_bigger[,1:N] + ecoli_cols_bigger[,3:(N+2)] - 2 * ecoli
    ecoli_cols_diff = D_ecoli * (1 / dx^2) * ecoli_cols_diff
    
    glucose_cols_bigger = cbind(glucose[,1], glucose, glucose[,N])
    glucose_cols_diff = glucose_cols_bigger[,1:N] + glucose_cols_bigger[,3:(N+2)] - 2 * glucose
    glucose_cols_diff = D_glucose * (1 / dx^2) * glucose_cols_diff
    ## Add flux gradient across columns to rate of change
    decoli_dt   = decoli_dt + ecoli_cols_diff
    dglucose_dt = dglucose_dt + glucose_cols_diff
    
    # reconvert the state variables back into vectors, bind them together, and return as list
    return(list(c(as.vector(decoli_dt), as.vector(dglucose_dt))))
  })
}

pars    <- c(vmax   = 1,    # /hr
             km = 0.01)     # half-saturation constant, [resource]



## initial conditions--need NxN for ecoli AND glucose
# get random locs for ecoli, mostly zeros
ecoli_ini = sample(c(rep(0, 250), 1), size = N*N, replace = TRUE)
# uniform resource concentration
glucose_ini = rep(10, N*N)
# all data
yini = c(ecoli_ini, glucose_ini)

# 18 hours
times   <- seq(0, 18, length.out = 50)

## solve model
# note how we can just supply the named spatial parms
# we also use Cash-Karp Runge-Kutta method for integration cuz its faster than lsoda
# 
# this takes about a minute on a 2022 mid-range desktop
out <- ode.2D(y = yini, times = times, func = monod2D, parms = pars,
              dimens = c(N, N), names = c("ecoli", "glucose"),
              N = N, D_ecoli = D_ecoli, D_glucose = D_glucose, dx = dx, 
              method = rkMethod("rk45ck"))

diagnostics(out)
summary(out)

# use "select" to grab specific state variables. this nicely makes arrays for you
bacteria <- subset(out, select = "ecoli", arr = TRUE)
resource <- subset(out, select = "glucose", arr = TRUE)

dim(bacteria)
par(mfrow = c(1,2))
# plot totals through time. margin 3 is time (1 = x, 2 = y)
mean_bacteria <- apply(bacteria, MARGIN = 3, FUN = sum)
plot(times, mean_bacteria)
plot(times, apply(resource, MARGIN = 3, FUN = sum))


# here is an image plot through time, arbitrarily choosing 6 slides
par(mfrow = c(4,3))
for (i in seq(from = 1, to = 50, length.out = 6)){
  image(bacteria[,,i])
}
# and the resource below
for (i in seq(from = 1, to = 50, length.out = 6)){
  image(resource[,,i])
}

# as a note, back when I did this for the Isme paper, I found that having the 
# diffusion calculations inside the solver sometimes resulted in out-of-memory
# issues when I was using (large) worlds of 100x100, with three different state variables.
#
# If you run into that, it is pretty straightforward to separate the steps.
# Then, you still use deSolve for the "reaction" (growth), but just manually
# do the diffusion calculations. The only real difference here is you multiply
# by a pre-specified dt, and you run all of this in a loop. Hopefully this isn't 
# necessary anymore, as it requires a very small dt. 
