#
# skąd ja to wziąłem ?????/!!!!!
# pierwsza wersja chyba z 2019.01.14
#
# nie napisałem sobie jak to działa jak to skleciłem
# komentarze dopisuję po fakcie próbując to zrozumieć
#
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(getwd())

library(signal)

source('doublePSD.R')


fsample = 1000
fullTm = 2

# liczba próbek
nsamp = fullTm * fsample

carrFreq = 210
fr = 20

# współczynnik wypełnienia
duty = 0.2
duty = 0.1

# w funkcji sinx/x przy generowaniu modSig, w argumencie sinusa jest 1/duty 
# czyli to okres sinusa czyli jedno z pierwszych zer sinx/x
# więc zgrubnie ta liczba razy odstęp prążków to pasmo
# band = fr/duty
# w każdym razie współczynnik wypełnienia jest powiązany z pasmem

# liczba prążków
nbar = 1 / duty
nbar = 20

# patrz pulse-gen.R
# prostsza wersja pulse train, symetryczna, bi = 0, 
# wypełnienie jest 2x większe, bo są 2 połówki
# zobacz w google "fourier series for common signals"
# w https://en.wikipedia.org/wiki/Fourier_series#Table_of_common_Fourier_series
# jest tylko ta bardziej skomplikowana
modSig = 0
for(i in 0:nbar)
{
    if (i == 0) 
        ai = duty 
    else 
        ai = 1/i/pi * sin(2*pi * i * duty) 
    carg = 2*pi * 1:(fullTm * fsample) * (fr * i) / fsample
    modSig = modSig + 2 * ai * cos(carg)
}
    
if(F)
{
    ft = fft(modSig)
    ft = c(ft[(nsamp / 2 + 1):nsamp],  ft[1:(nsamp / 2)])
    plot(Re(ft), type = 'l', col = 'blue')
    plot(Im(ft), type = 'l', col = 'red')
}

#plot(modSig[1:(fsample * 4 / fr)], type='l')

#combSig = modSig * cos(2 * pi * 1:(fullTm * fsample) * carrFreq / fsample)

plot(modSig[1:(fsample * 4 / fr)], type='l')
plot(modSig[1:(fsample * 1 / fr)], type='o')

if(F)
{
    #ps = welchPSD(ts((modSig), frequency = fsample), seglength = 1024 / fsample)
    ps = doublePSD(modSig, fsample, 1024)
    plot(ps$frequency, ps$power, type = 'l')
}

# modulacja --> przeniesienie na wyższą częstotliwość
# czyli obwiednia to sygnał impulsu, 
# im krótszy impuls tym większe pasmo / więcej prążków
combSig = modSig * cos(2 * pi * 1:nsamp * carrFreq / fsample)

plot(combSig[1:(fsample * 4 / fr)], type='l')
plot(combSig, type='l')



ps = welchPSD(ts((combSig), frequency = fsample), seglength = 1024 / fsample)
plot(ps$frequency, ps$power, type = 'l')

rn = which(ps$frequency >carrFreq - fr/duty * 1.1 & ps$frequency < carrFreq + fr/duty * 1.1)
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
    

