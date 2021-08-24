library(ggplot2)
library(gridExtra)
library(mvtnorm)

# Function for P(||X|| < r), where X is a d-dimensional t variable with v d.f.
pT = function(r, d, v){
  pf(r^2/d, d, v)
}

# Function for the mass of the Gaussian approximation to f_X over ||x|| < r
# Note that this converges to the LA as r -> Inf
pGauss = function(r, d, v){
  r0 = sqrt((v+d)/v)*r
  exp(lgamma((v+d)/2) - lgamma(v/2) + d*log(2/(v+d))/2)*pchisq(r0^2, d)
}

par(mfrow = c(2,1))
plot(function(x) pT(x, 72, 25921), xlim = c(0, 12), lwd = 2)
plot(function(x) pGauss(x, 72, 25921), col = 'red', add = TRUE, xlim = c(0,12))

tau = function(x, d, v){
  dmvt(cbind(matrix(x), matrix(0, length(x), d - 1)), df = v, log = FALSE)
}

gauss_approx = function(x, d, v){
  exp(lgamma((v+d)/2) - lgamma(v/2) - d*log(v*pi)/2 - x^2*(v+d)/(2*v))
}

plot(function(x) (tau(x, 72, 25921) - gauss_approx(x, 72, 25921))/tau(0,72,25921), xlim = c(0,12))
