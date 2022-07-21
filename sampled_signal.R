
fsig = 10
fsamp = 1000
idx = 0:500 - 1
tm = idx / fsamp
sig = 3 * sin(2 * pi * fsig / fsamp * idx)
plot(tm, sig, type = 'l')



