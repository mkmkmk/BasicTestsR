# test pokazuje generowanie ciągu impulsów ze składowych sinusoidalnych
#
# patrz >> ad DSP
#
#' dla przebiegu pulse-train w funkcji sinx/x w argumencie sinusa 
#' jest współczynnik wypełnienia, czyli jego odwrotność to okres sinusa, 
#' czyli jedno z pierwszych zer sinx/x, więc zgrubnie ta liczba razy odstęp 
#' prążków określa pasmo, w każdym razie współczynnik wypełnienia jest powiązany 
#' z pasmem, gdy współczynnik wypełnienia zmierza do zera, to impulsy (w granicy) 
#' stają się impulsami Diracka, wtedy też okres funkcji sinx/x staje się nieskończony 
#' amplitudy wszystkich prążków (dziedzina częstotliwości) stają się podobne 
#' i ciągną się do nieskończoności
#'


rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(getwd())

library(signal)

source('doublePSD.R')

# częstotliwość próbkowania (Hz)
fsample = 2000

# czas (s)
fullTm = 2

# liczba próbek
nsamp = fullTm * fsample


# częstotliwość sygnału - impulsów (Hz)
fr = 20

# współczynnik wypełnienia
duty = 0.05
duty = 0.2
duty = 0.1

# liczba prążków
nbar = 1/duty
nbar = 100
nbar = 200
nbar = 20

# wypełnienie zerami w próbkach
pd = fsample * (1 / fr  - duty / fr)

stopifnot(pd > 0)

#klasycznie generuję ciąg impulsów
modSig0 = 0
for(i in 1:(fullTm * fr))
    modSig0 = c(modSig0, rep(1, fsample * duty / fr), rep(0, pd))
modSig0 = modSig0[1:(fullTm * fsample)]

plot(modSig0[1:(fsample * 4 / fr)], type='l')

if(F)
{
    ps = doublePSD(modSig0, fsample, 1024)
    plot(ps$frequency, ps$power, type = 'l')
}


# https://en.wikipedia.org/wiki/Fourier_series#Table_of_common_Fourier_series
modSig = 0
for(i in 0:nbar)
{
    if (i == 0) 
    {
        ai = duty 
        bi = 0
    }   
    else 
    {
        ai = 1/i/pi * sin(2*pi * i * duty) 
        bi = 2/i/pi * sin(pi * i * duty)^2 
    }
    carg = 2*pi * 1:nsamp * (fr * i) / fsample
    modSig = modSig + ai * cos(carg) + bi * sin(carg)
}

ylim = range(modSig) * 1.1
plot(modSig0[1:(fsample * 1 / fr)], type='l', col = 'blue', ylim = ylim)
lines(modSig[1:(fsample * 1 / fr)], type='l', col = 'red')

plot(modSig0[1:(fsample * 4 / fr)], type='l', col = 'blue', ylim = ylim)
lines(modSig[1:(fsample * 4 / fr)], type='l', col = 'red')

# + prostsza wersja, symetryczna, bi = 0, wypełnienie jest 2x większe, bo są 2 połówki
# zobacz w google "fourier series for common signals"
modSig2 = 0
for(i in 0:nbar)
{
    if (i == 0) 
        ai = duty 
    else 
        ai = 1/i/pi * sin(2*pi * i * duty) 
    carg = 2*pi * 1:nsamp * (fr * i) / fsample
    modSig2 = modSig2 + 2 * ai * cos(carg)
}

ylim = range(modSig2) * 1.1
plot(modSig0[1:(fsample * 1 / fr)], type='l', col = 'blue', ylim = ylim)
lines(modSig2[1:(fsample * 1 / fr)], type='l', col = 'red')

plot(modSig0[1:(fsample * 4 / fr)], type='l', col = 'blue', ylim = ylim)
lines(modSig2[1:(fsample * 4 / fr)], type='l', col = 'red')


# -------- -------- analiza widma -------- -------- 

# x prążków
xf = (-nsamp/2+1):(nsamp/2)

f0 = fft(modSig0)
# zamiana połówek
f0 = c(f0[(nsamp/2+1):nsamp], f0[1 : (nsamp/2)] )

f = fft(modSig)
# zamiana połówek
f = c(f[(nsamp/2+1):nsamp], f[1 : (nsamp/2)] )

plot(xf, Re(f0), type = 'l', col = 'red', ylim = max(Mod(f0))*c(-1,1))
lines(xf, Re(f), type = 'l', col = 'blue')

plot(xf, Im(f0), type = 'l', col = 'red', ylim = max(Mod(f0))*c(-1,1))
lines(xf, Im(f), type = 'l', col = 'blue')



