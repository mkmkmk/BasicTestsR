
len = 20000

# cand = exp(1i * runif(len, 0, pi / 2))
# curMax = runif(len, 0, 2) * exp(1i * runif(len, 0, pi / 2))

cand = runif(len, 0, 2)
curMax = runif(len, 0, 2)

doPlot = T
doPlot = F 

if(doPlot)
{
    plot(cand, col = "blue", type = 'p')
    lines(curMax, col = "blue", type = 'p')
}

anyErr = F

for(i in 1:len)
{
     reject = abs(Re(cand[i])) < abs(Re(curMax[i])) & abs(Im(cand[i])) < abs(Im(curMax[i]))     
    
     error = reject & Mod(cand[i]) > Mod(curMax[i])
    
     anyErr = anyErr | error
     
     if(doPlot & error)
     {
         lines(c(cand[i], curMax[i]), type = 'l')
     }
}

stopifnot(!anyErr)