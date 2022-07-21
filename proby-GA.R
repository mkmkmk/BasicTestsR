# ------------------------------------
# https://cran.r-project.org/web/packages/GA/vignettes/GA.html
# ------------------------------------

if (!("GA" %in% installed.packages()))
    install.packages("GA")
library(GA)


# ------------------------------------

f <- function(x)  (x^2 + x) * cos(x)
lbound <- -10; ubound <- 10
curve(f, from = lbound, to = ubound, n = 1000)


GA <- ga(type = "real-valued", fitness = f, lower = c(th = lbound), upper = ubound)

summary(GA)

curve(f, from = lbound, to = ubound, n = 1000)
points(GA@population, GA@fitness, pch = 20, col = 4)
points(GA@solution, GA@fitnessValue, col = 2, pch = 19)
rug(GA@population, col = 2)


# ------------------------------------
Rastrigin <- function(x1, x2)
{
    20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
}

x1 <- x2 <- seq(-5.12, 5.12, by = 0.1)
f <- outer(x1, x2, Rastrigin)
persp3D(x1, x2, f, theta = 50, phi = 20, col.palette = bl2gr.colors)

filled.contour(x1, x2, f, color.palette = bl2gr.colors)

GA <- ga(type = "real-valued", 
         fitness =  function(x) -Rastrigin(x[1], x[2]),
         lower = c(-5.12, -5.12), upper = c(5.12, 5.12), 
         popSize = 50, maxiter = 1000, run = 100)


summary(GA)

plot(GA)
filled.contour(x1, x2, f, color.palette = bl2gr.colors, 
               plot.axes = { axis(1); axis(2); 
                   points(GA@solution[,1], GA@solution[,2], 
                          pch = 3, cex = 2, col = "white", lwd = 2) }
                )

monitor <- function(obj) 
{ 
    contour(x1, x2, f, drawlabels = FALSE, col = grey(0.5))
    title(paste("iteration =", obj@iter), font.main = 1)
    points(obj@population, pch = 20, col = 2)
    Sys.sleep(0.2)
}

if (F)
GA <- ga(type = "real-valued", 
         fitness =  function(x) -Rastrigin(x[1], x[2]),
         lower = c(-5.12, -5.12), upper = c(5.12, 5.12), 
         popSize = 50, maxiter = 100, 
         monitor = monitor)


suggestedSol <- matrix(c(0.2,1.5,-1.5,0.5), nrow = 2, ncol = 2, byrow = TRUE)
GA1 <- ga(type = "real-valued", 
          fitness =  function(x) -Rastrigin(x[1], x[2]),
          lower = c(-5.12, -5.12), upper = c(5.12, 5.12), 
          suggestions = suggestedSol,
          popSize = 50, maxiter = 1)
head(GA1@population)


# ------------------------------------
f <- function(x)
{ 100 * (x[1]^2 - x[2])^2 + (1 - x[1])^2 }

c1 <- function(x) 
{ x[1]*x[2] + x[1] - x[2] + 1.5 }

c2 <- function(x) 
{ 10 - x[1]*x[2] }

ngrid <- 250
x1 <- seq(0, 1, length = ngrid)
x2 <- seq(0, 13, length = ngrid)
x12 <- expand.grid(x1, x2)
col <- adjustcolor(bl2gr.colors(4)[2:3], alpha = 0.2)
plot(x1, x2, type = "n", xaxs = "i", yaxs = "i")
image(x1, x2, matrix(ifelse(apply(x12, 1, c1) <= 0, 0, NA), ngrid, ngrid), 
      col = col[1], add = TRUE)
image(x1, x2, matrix(ifelse(apply(x12, 1, c2) <= 0, 0, NA), ngrid, ngrid), 
      col = col[2], add = TRUE)
contour(x1, x2, matrix(apply(x12, 1, f), ngrid, ngrid), 
        nlevels = 21, add = TRUE)

fitness <- function(x) 
{ 
    f <- -f(x)                         # we need to maximise -f(x)
    pen <- sqrt(.Machine$double.xmax)  # penalty term
    penalty1 <- max(c1(x),0)*pen       # penalisation for 1st inequality constraint
    penalty2 <- max(c2(x),0)*pen       # penalisation for 2nd inequality constraint
    f - penalty1 - penalty2            # fitness function value
}


GA <- ga("real-valued", fitness = fitness, 
         lower = c(0,0), upper = c(1,13), 
         # selection = GA:::gareal_lsSelection_R,
         maxiter = 1000, run = 200, seed = 123)


summary(GA)


points(GA@solution[1], GA@solution[2], col = "dodgerblue3", pch = 3)  # GA solution


# ------------------------------------


AQL   <- 0.01; alpha <- 0.05
LTPD  <- 0.06; beta  <- 0.10
plot(0, 0, type="n", xlim=c(0,0.2), ylim=c(0,1), bty="l", xaxs="i", yaxs="i", 
     ylab="Prob. of acceptance", xlab=expression(p))
lines(c(0,AQL), rep(1-alpha,2), lty=2, col="grey")
lines(rep(AQL,2), c(1-alpha,0), lty=2, col="grey")
lines(c(0,LTPD), rep(beta,2), lty=2, col="grey")
lines(rep(LTPD,2), c(beta,0), lty=2, col="grey")
points(c(AQL, LTPD), c(1-alpha, beta), pch=16)
text(AQL, 1-alpha, labels=expression(paste("(", AQL, ", ", 1-alpha, ")")), pos=4)
text(LTPD, beta, labels=expression(paste("(", LTPD, ", ", beta, ")")), pos=4)



decode1 <- function(x)
{ 
    x <- gray2binary(x)
    n <- binary2decimal(x[1:l1])
    c <- min(n, binary2decimal(x[(l1+1):(l1+l2)]))
    out <- structure(c(n,c), names = c("n", "c"))
    return(out)
}

fitness1 <- function(x) 
{ 
    par <- decode1(x)
    n <- par[1]  # sample size
    c <- par[2]  # acceptance number
    Pa1 <- pbinom(c, n, AQL)
    Pa2 <- pbinom(c, n, LTPD)
    Loss <- (Pa1-(1-alpha))^2 + (Pa2-beta)^2
    -Loss
}

if (F)
{
    sapply(0:16, function(x) decimal2binary(x, 5))
    sapply(0:16, function(x) binary2gray(decimal2binary(x, 5)))
    str(structure(1:2, names = c("n", "c")))
    l1 = 5
    l2 = 5
    decimal2binary(666, l1 + l2)
    decimal2binary(3, l1 + l2)
    decode1(binary2gray(decimal2binary(666, l1 + l2)))
    decode1(binary2gray(decimal2binary(6 + 16 * 2^5, l1 + l2)))
}


n  <- 2:200                  # range of values to search
b1 <- decimal2binary(max(n)) # max number of bits requires
l1 <- length(b1)             # length of bits needed for encoding
c  <- 0:20                   # range of values to search
b2 <- decimal2binary(max(c)) # max number of bits requires
l2 <- length(b2)             # length of bits needed for encoding

GA1 <- ga(type = "binary", fitness = fitness1, 
          nBits = l1 + l2, 
          popSize = 100, maxiter = 1000, run = 100)
summary(GA1)

sol = decode1(GA1@solution)
sol
n <- sol[1]; c <- sol[2]
p <- seq(0, 0.2, by = 0.001)
Pa <- pbinom(c, n, p)
lines(p, Pa, col = 2)









