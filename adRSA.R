# ad szyfrowanie RSA
#
# https://www.scottaaronson.com/democritus/lec8.html
#
#

library(Rcpp)

epsilon = 1e-12

c = 11
c = 7
c = 5
c = 3

mx = 10000
mx = 1000
mx = 100 
mx = 50

# https://en.wikipedia.org/wiki/Modular_exponentiation
cfun = 
    "
        // [[Rcpp::export]]
        double modpows(int base, int exponent, int modulus)
        {
            if (modulus <= 1)
               return 0;
            long long c = 1LL;
            for (int i = 0; i < exponent; ++i)
               c = (c * base) % modulus;
            return c;
        }
    "
sourceCpp(code = cfun)

if (FALSE)
{
    modpows = function(base, exponent, modulus)
    {
        if (modulus == 1) 
            return(0)
        c = 1
        eprime = 0
        while (eprime < exponent) 
        {
            c = (c * base) %% modulus
            eprime = eprime + 1
        }    
        return(c)
    }
}

modpow = function(base, exponent, modulus) 
    sapply(exponent, function(it) modpows(base, it, modulus))

stopifnot(c(
    modpow(5, 3, 13) == 5^3 %% 13,
    modpow(2, 100, 7) == 2,
    modpow(3, 100, 7) == 4,
    modpow(10, 333, 5) == 0,
    modpow(10, 333, 10) == 0,
    modpow(10, 333, 2) == 0,
    modpow(666, 0, 666) == 1
    ))


if (FALSE)
{
    if (c == 3)
    {
        p = 11
        q = 5
        x = 41
        
        p = 11
        q = 17
        x = 128
        
        p = 17
        q = 23
        x = 333
    } else if (c == 5)
    {
        p = 17 
        q = 13 
        x = 111
    } else
    {
        stop("not implemented")
    }
}else
{
    while(1)
    {
        p = trunc(runif(1, c, mx))
        q = trunc(runif(1, c, mx))
        if (p==q)
            next
        if (!all(abs(q / 2:(q-1) - trunc(q / 2:(q-1))) > epsilon))
            next
        if (!all(abs(p / 2:(p-1) - trunc(p / 2:(p-1))) > epsilon))
            next
        Npq = (p-1) * (q-1)
        if (Npq/c - trunc(Npq/c) != 0)
        {
            x = trunc(runif(1, p*q/2, p*q))
            cat(sprintf("p=%g q=%g x=%g\n", p, q, x))
            break
        }
    }    
}

N = p * q
N   

Npq = (p-1)*(q-1)
Npq
Npq / c

stopifnot(x < N)
stopifnot(Npq / c - trunc(Npq / c) != 0)

if (FALSE)
{
    modpow(x, 0:20, N)
    modpow(x, 0:20 + Npq, N)
    modpow(x, 1, N)
    modpow(x, c, N)
}
retval = modpow(x, c, N)
retval

if(FALSE)
{
    modpow(x, c, N)
    modpow(retval, 1, N)
}

# podnosząc retval do kolejnych potęg przeskakuję o c pozycje w szeregu
# ale k jest takie że (c*k) stoi następne po końcu okresu (p-1)(q-1), właśnie tam gdzie siedzi x
# (po końcu okresu siedzi 1, bo to jest kopia x^0, następne jest x^1)
k = which((c*(1:Npq)) %% Npq == 1)[1]
stopifnot(!is.na(k))
k
(c*k) %% Npq
(c*(1:Npq)) %% Npq

#to już kończy sprawę !! to jest x !!
res = modpow(retval, k, N)
stopifnot(res == x)
cat("OK!\n")

# -------
# + wygląda że okres jest kilka mniejszy niż Npq - to ciekawe!
if(F)
{
    which(modpow(retval, 1:Npq, N) == x)
    diff(which(modpow(retval, 1:Npq, N) == x))
    f = diff(which(modpow(retval, 1:Npq, N) == x))[1]
    print(Npq/f)
    modpow(x, 0:k, N)
    modpow(retval, 0:(k+10), N)
}
