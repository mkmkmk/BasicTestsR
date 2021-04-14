rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(getwd())

library(signal)

source('doublePSD.R')


fsample = 1000
fullTm = 2

carrFreq = 210
fr = 20
band = 100
# band = 40


if(F)
{
    pd = fsample * (1 / fr  - 1 / band)
    
    stopifnot(pd > 0)
    
    modSig = 0
    for(i in 1:(fullTm * fr))
        modSig = c(modSig, rep(0, pd), rep(1, fsample / band))
    modSig = modSig[1:(fullTm * fsample)]

}else
{

    modSig = 0
    for(i in -250:250)
    {
        if (i == 0) 
            modSig = modSig + fr / band
        else
            modSig = modSig + 1 / i / pi * sin(pi * fr / band * i) * cos(2 * pi * 1:(fullTm * fsample) * (fr * i) / fsample)
        
    }
    
    if(F)
    {
        N = length(modSig)
        ft = fft(modSig)
        ft = c(ft[(N / 2 + 1):N],  ft[1:(N / 2)])
        plot(Re(ft), type = 'l', col = 'blue')
        plot(Im(ft), type = 'l', col = 'red')
    }
    
    #plot(modSig[1:(fsample * 4 / fr)], type='l')
    
    combSig = modSig * cos(2 * pi * 1:(fullTm * fsample) * carrFreq / fsample)

}

plot(modSig[1:(fsample * 4 / fr)], type='l')

plot(modSig[1:(fsample * 1 / fr)], type='o')

if(F)
{
    #ps = welchPSD(ts((modSig), frequency = fsample), seglength = 1024 / fsample)
    ps = doublePSD(modSig, fsample, 1024)
    plot(ps$frequency, ps$power, type = 'l')
}


combSig = modSig * cos(2 * pi * 1:(fullTm * fsample) * carrFreq / fsample)

plot(combSig[1:(fsample * 4 / fr)], type='l')
plot(combSig, type='l')



ps = welchPSD(ts((combSig), frequency = fsample), seglength = 1024 / fsample)
plot(ps$frequency, ps$power, type = 'l')

rn = which(ps$frequency >carrFreq - band * 1.1 & ps$frequency < carrFreq + band * 1.1)
plot(ps$frequency[rn], ps$power[rn], type = 'l')


if(F)
    plot(ps$frequency, 10 * log10(ps$power), type = 'l')

if(F)
{
    f = fft(combSig)
    f = f[1:(length(f)/2)]
    a = (Arg(f)) %% (2*pi) - pi
    plot(1:(length(f)) / 2, a, type = 'l', col='blue')
    lines(1:(length(f)) / 2, Mod(f)/100, type = 'l', col='red')
}  
    

