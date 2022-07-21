rm(list = ls())

eps = 1e-10
N = 18
alpha = -N:N

x = cos(alpha * pi / N)
y = sin(alpha * pi / N)

plot(x, y, col = "blue")

myatan = function(z)
    z - z^3/3 + z^5/5 - z^7/7 + z^9/9 - z^10/10

myatan2 = function(yt,xt)
{
    # atan(1e-3)
    # if(abs(x)<1e-5)return
    # res = myatan(y / x)
    # res[x < 0 & y < 0] = res[x < 0 & y < 0] - pi
    # res[x < 0 & y > 0] = res[x < 0 & y > 0] + pi
    # res
    n = length(x)
    stopifnot(n == length(y))
    res = xt
    for(i in 1:n)
    {
        x = xt[i]
        y = yt[i]
        if (abs(y) < abs(x) * 1e-3)
            res[i] = y / x
        else if(abs(y) > abs(x) * 1e3 )
            res[i] = if (y > 0) pi / 2 else -pi / 2
        else
            res[i] = atan(y / x)

        if(x < 0)
            res[i] = res[i] + (if (y>0) pi else -pi)
    }
    res
}

if (F)
{
    atan2(y,x)*N/pi
    plot(round(atan(y/x)*N/pi))
    plot(round(myatan2(y,x)*N/pi))
}

plot(alpha, col="blue")
lines(atan2(y,x)*N/pi, col="red", type = "p")
lines(myatan2(y,x)*N/pi, col="black", type = "p")


stopifnot(atan2(y,x) - myatan2(y,x) < eps)


atan(116)/ pi * 180
atan(1) / pi * 180

at = atan(seq(0, 120, 0.05)) / pi * 180

a = seq(0, 120, 0.05)
a[diff(round(atan(a) / pi * 180)) != 0]

at[abs(at-90)<0.5]

at[abs(at-0)<0.5]
at[abs(at-1)<0.5]
at[abs(at-2)<0.5]
at[abs(at-3)<0.5]



length(at)

abs(at-90)

plot(a, atan(a)/pi*180, type = 'l', col = "red")


