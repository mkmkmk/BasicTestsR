
# to ma przypominać sytuację jaka jest z próbkowaniem sygnału bez i z filtrem antyaliasingowym, 
# ciąg z b. wysoką cz. próbkowania ma przypominać sygnał analogowy, który jest próbkowany;
# szum biały z definicji ma równą moc w całym paśmie i jeżeli nie ma filtra anyaliasingowego 
# to szum z całego dostępnego pasma sygnału oryginalnego przez aliasing wchodzi 
# w pasmo sygnału spróbkowanego (od 0 do połowy cz. próbkowania) i dodatkowo zakłóca 
# sygnał spróbkowany - przebieg czerwony; na niebiesko ten sam sygnał, 
# który został spróbkowany po zastosowaniu filtra antyaliasingowego

rm(list = ls())

library(signal)
library(bspec)

dur = 5
# duża cz. próbkowania sygnału udającego sygnał analogowy
fsamp = 100e3
fsig = 3.3333333
len = dur * fsamp

t = 1:len / fsamp

range(t)

sig = 1/8 * sin(t * 2 * pi * fsig)
sig = 1/8 * sign(sin(t * 2 * pi * fsig))
sig = 1.1 * sign(sin(t * 2 * pi * fsig))

noise = runif(len) - .5
noise  = rnorm(len)

noise = noise / sd(noise)
sig1 = noise + sig

sd(sig) / sd(sig1)

q = 1000

# nowa cz. próbkowania
fsamp / q
fsamp / 2 * .5 / q

sampSeq = seq(1, length(sig1), q)

# próbkowanie bez filtra
sig2 = sig1[sampSeq]

length(sig1) / length(sig2)

if(F)
    plot(sig1, type = 'l', col = 'blue3')

sd(sig1)
sd(sig2)


decFilter = butter(5, .8 / q)

# próbkowanie z filtrem
sig2ft = filter(decFilter, sig1)[sampSeq]

sd(sig2ft)

gn = sd(sig2) / sd(sig2ft)


plot(sig2, type = 'l', col = 'red3')

plot(sig2, type = 'l', col = 'red3')
lines(gn * sig2ft, type = 'l', col = 'blue3')


# plot(sig2ft, type = 'l', col = 'blue3')

if(F)
{
    sigIn = sig2ft
    sigIn = sig2
    psd = welchPSD(ts((sigIn), frequency = fsamp / q ),
             seglength = 2^8 / (fsamp / q),
             two.sided = T)
    
    diff(psd$frequency)[1]
    
    frg = psd$frequency < 50
    plot(psd$frequency[frg], psd$power[frg], col='red3', type = 'l')
}







