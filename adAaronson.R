bufSize = 128
resampPeriod = 1/25


bpmMax = 60 / 2 / resampPeriod
bpmMax

bpm1 = 2 * bpmMax / bufSize
bpm1

60 / resampPeriod / bufSize


bpm1*bufSize/2




bufSize = 128
resampPeriod = 1/12


bufSize = 256
resampPeriod = 1/25

n = 20
n2 = n^2

n/log(n)
n2/log(n2)

2^n/log(2^n)

2^n2/log(2^n2,2)
2^n2/n2

# The famous Prime Number Theorem liczba liczb pierwszych w x : x/log(x)
# n^2 cyfr to - dla sys. binarnego - max liczba jest 2^(n^2) i jak tę liczbę wstawimy do wzoru x/log(x)
# + zastosujemy podstawę binarną logarytmu (bo to bardzo już nie zmienia wyniku) to mamy
2^(n^2)/log2(2^(n^2))
# +
log2(2^(n^2)) == n^2
# czyli 
2^(n^2)/n^2

# --------------------
    
(1/3)^2 + 2*1/3*2/3 + (2/3)^2
1/9 + 4/9 + 4/9

# ----------------------

matrix(c(1,-1,1,1), nrow=2) 


# ---------------------- ad quantum 

A = t(t(c(1,1))) / sqrt(2)
B = t(t(c(1,-1))) / sqrt(2)
1/2 * A %*% t(A) + 1/2 * B %*% t(B) 

C = t(t(c(1,0))) / sqrt(2)
D = t(t(c(0,1))) / sqrt(2)
1/2 * A %*% t(A) + 1/2 * B %*% t(B) 


2^1000
10^300


2^100
10^30

2^100
10^(100 / log2(10))

log(10,2)
301/1000

n=100
2^(-2*n)
2^(2*n)

----------------
# w obliczeniach rzędu złożoności w ogóle nie ma znaczenia jaka jest podstawa logarytmu
# bo logarytmy różnia się tylko stałą np. dla log10() vs log2() jest to log2(10) == log10(2)^-1 
 
log10(4)
log2(4)

log2(10) - log10(2)^-1
log2(exp(1)) - log(2)^-1

l = runif(1, 100, 10e10)
log10(l)
log2(l)/log2(10)

log10    
    
2^1000

log10(2^1000) / log10(2)
log2(2^1000) / log2(10)

# --------------------------
# ad lecture 15
C = 1e30
eps = 1e-3
dt = 0.1
m = 1/eps * log(C/dt)
# ==
C*(1 - eps)^m    
dt
# ==
(1 - eps)^m    
dt/C
# ==
m * log(1 - eps)
log(dt / C)
# ==
m
log(C / dt) / -log(1 - eps)
# ==
m
log(C / dt) / log(1 / (1 - eps))
#| log(1 / (1 - eps)) ==> eps, bo 1 / (1 - eps) ==> 1 i wtedy
#| log() ==> 0, czyli to też jest jakiś eps
# ==
m
log(C / dt) / eps

-log(1 - eps)

# ----------------------
# reszta była o 
# https://www.scottaaronson.com/qclec.pdf
# poszła do adAaronsonQlec.R
# ----------------------

