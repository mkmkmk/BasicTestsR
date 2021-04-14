rm(list = ls())
library(signal)

# y[n] = (x[n] - y[n-1]) * K + y[n-1]
# 
# y[n] = x[n] * K - y[n-1] * K + y[n-1]
# 
# y[n] = K * x[n] + (1 - K) * y[n-1]


# y[n] = (x[n]/K - y[n-1]) * K + y[n-1]
# y[n] = (x[n] - y[n-1] * K) + y[n-1]

refFlt = function(x, K)
{
    K_INV = 1/K
    y = rep(0, length(x))
    for(n in 2:length(sig))
       y[n] = y[n-1] + x[n] - trunc(y[n-1] / K_INV)
       
       #y[n] =(x[n] - y[n-1]) / K_INV + y[n-1]
       #y[n] = (x[n] - y[n-1]) * K + y[n-1]
       #y[n] = K * x[n] + (1 - K) * y[n-1]   
    # y 
    y / K_INV
}


mul = 1
K = 0.0247 / mul
K = 1 / 32
K = 0.0247

mul = 20
K = 0.0247 / mul

mul = 20
K = 1 / 128
K = 1 / 2^12
{

    fsamp = 50 * mul
    
    fz = freqz(K, c(1, -(1 - K)), Fs = fsamp, n = 2^16)
    fz
    fz$db = 10 * log10 (Mod(fz$h))
    plot(fz$f, fz$db, type = 'l', ylim = c(-3, 0))
    
    c1 = fz$f[which.min(abs(fz$db - (-3)))]
    cat(sprintf("c1 = %.3g Hz\n", c1))
    
    sig = rep(1, 20 * fsamp)*10
    fsig = filter(K, c(1,-(1 - K)), sig)
    fsigRef = refFlt(sig, K)
    tm = 1:length(sig) / fsamp
    plot(tm, fsig, type = 'l', col='blue3', xlab = "t [s]", xlim = c(0, 20))
    lines(tm, fsigRef, type = 'l', col='red4', xlab = "t [s]")
}





