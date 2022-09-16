# “Upsampling” is the process of inserting zero-valued samples 
# between original samples to increase the sampling rate. (This is called “zero-stuffing”.) 
# Upsampling adds to the original signal undesired spectral images 
# which are centered on multiples of the original sampling rate.
# “Interpolation”, in the DSP sense, is the process of upsampling followed by filtering. 
# (The filtering removes the undesired spectral images.) 
# As a linear process, the DSP sense of interpolation is somewhat different 
# from the “math” sense of interpolation, but the result is conceptually similar: 
# to create “in-between” samples from the original samples. 
# The result is as if you had just originally sampled your signal at the higher rate.
# http://www.dspguru.com/dsp/faqs/multirate/interpolation
#
# >> ad DSP
#

rm(list = ls())

library(signal)
library(bspec)

dur = 5
fsamp = 100
fsig = 3.3333333
len = dur * fsamp
tm = 1:len / fsamp

up = 6

range(tm)

sig = 1.1 * sign(sin(tm * 2 * pi * fsig))
sig = 1/8 * sin(tm * 2 * pi * fsig)
sig = 1/8 * (sin(tm * 2 * pi * fsig) + 0.5 * sin(pi/4 + tm * 2 * pi * 2 * fsig))
sig = 1/8 * sign(sin(tm * 2 * pi * fsig))


#if (F)
# sygnał oryginalny
plot(tm, sig, type = 'l')


uplen = up * len
upfsamp = up * fsamp

upsig = rep(0, uplen)
uptm = 1:uplen / upfsamp

upsig[seq(1, uplen, up)] = sig

# sygnał z upchanymi zerami
plot(uptm, upsig, type = 'l')

# widmo syg. oryginalnego
ps = welchPSD(ts(sig, frequency = fsamp), seglength = dur / 2)
plot(ps$frequency, ps$power, type = 'l')

# widmo syg. z upchanymi zerami
ps = welchPSD(ts(upsig, frequency = upfsamp), seglength = dur / 2)
plot(ps$frequency, ps$power, type = 'l')

decFilter = butter(4, .6 / up)

upsigf = filter(decFilter, upsig)

# sygnał upsamplowany
plot(uptm, upsigf, type = 'l')
