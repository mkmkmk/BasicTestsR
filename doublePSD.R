library(bspec)

doublePSD = function(sig, fsample, N = 128)
{
  specR <-
    welchPSD(ts((sig), frequency = fsample),
             seglength = N / fsample,
             two.sided = T)
  specL <-
    welchPSD(
      ts(Conj(sig), frequency = fsample),
      seglength = N / fsample,
      two.sided = T
    )
  freqT = c(rev(-specL$frequency), specR$frequency)
  powerT = c(rev(specL$power), specR$power)
  list(frequency = freqT,  power = powerT)
}
