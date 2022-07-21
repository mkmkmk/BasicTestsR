rm(list = ls())

x = 0:1000*2*pi/1000
if (F)
    plot(x, atan(x), type = 'l', col='blue')

x = 0:1000/1000
dt = .25
for(i in 0:3)
{
    st = i*.25
    xr = x[x>=st & x<st+dt]
    yr = atan(xr)
    alm = lm(yr~xr)
    print(c(st, st+dt, as.numeric(alm$coefficients) * 180 / pi))
}

if (F)
{
    ap0 = function(x)
    {
        sg = sign(x)
        x = abs(x)
        inv = x > 1
        x[inv] = 1 / x
    
        y = x
        y[x<.25] = x[x<.25] * 56
        y[x>.25 & x<.5] = x[x>.25 & x<.5] * 50 + 1.5
        y[x>.5 & x<.75] = x[x>.5 & x<.75] * 40 + 6.5
        y[x>.75] = x[x>.75] * 32 + 13
    
        y[inv] = 90 - y[inv]
        y * pi / 180
    }
}

x = 0.7763656
ap0_one = function(x)
{
    sg = sign(x)
    x = abs(x)
    inv = x > 1
    if (inv)
        x = 1 / x
    if (x <= .25)
        y = x * 56
    else if (x <= .5)
        y = x * 50 + 1.5
    else if (x <= .75)
        y = x * 40 + 6.5
    else 
        y = x * 32 + 13
    if (inv)    
        y = 90 - y
    sg * y * pi / 180
}

ap0 = function(x)
    sapply(x, ap0_one)


x = 0:1000 * 0.3 * pi / 1000
stopifnot(max(abs(ap0(x) - atan(x))) < 0.010)
max(abs(ap0(x) - atan(x)))

x = 0:1000 * 2 * pi / 1000
stopifnot(max(abs(ap0(x) - atan(x))) < 0.010)
max(abs(ap0(x) - atan(x)))

if (F)
{
    plot(x, atan(x), type = 'l', col='blue')
    lines(x, ap0(x), type = 'l', col='red')
}





