library(signal)

sig_freq = 12.8

fsample = 100

duration = 0.4
1/.4
N = duration * fsample

print(N)

tm = 1:N / fsample


sig = exp(2i * pi * sig_freq * tm)

if(F)
    sig = Re(sig)


if(F)
{
    # sig = hanning(N) * sig
    sig = bartlett(N) * sig
}

if(F)
    plot(bartlett(N))


length(sig)

plot(tm, Re(sig), type = 'l')

spec = Mod(fft(sig))
maxSpec = max(spec)

if(F)
    spec = 10 * log10(spec) - 10 * log10(maxSpec)


plot(spec, type = 'o')


print(length(which(spec / max(spec) > 0.1)))




