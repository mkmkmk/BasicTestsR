
library(Rcpp)

sourceCpp(code=
    "
    // [[Rcpp::export]]
    int rshift(int val, int n){ return val >> n; }")


sourceCpp(code=
    "
    #include <stdlib.h>
    #define _SIGN(value) (value >= 0 ? 1 : -1)
    #define _SCALE_1357(value, mul_Q10) (1 + (((abs(value) * (mul_Q10 >> 1)) >> 10) << 1)) * _SIGN(value)

    // [[Rcpp::export]]
    int q_mul(int samp, int sig_amp_q10)
    { 
        //return (1 + (((abs(samp) * (sig_amp_q10 >> 1)) >> 10) << 1)) * (samp >= 0 ? 1 : -1); 
        return _SCALE_1357(samp, sig_amp_q10);
    }")





fsamp = 100
fsig = 4
idx = (-5):5

#sig = trunc(128 * sin(2 * pi * idx * fsig/ fsamp))
#sig = 128 * idx

inp = (-128):128
# plot(idx, sig, type='o', col='blue4')
# sig2 = trunc(sig / 32)
# sig2i = sapply(sig, function(s) rshift(s, 5))
# plot(idx, sig2, type='o', col='red4')
# plot(idx, sig2i, type='o', col='red4')

# out = inp
# out[inp>=0] = trunc((inp[inp>=0]+31) / 32)
# out[inp<0] = trunc((inp[inp<0]-31) / 32)

# out2 = inp
# out2[inp>=0] = trunc((inp[inp>=0]+31) / 32)
# out2[inp<0] = -trunc((abs(inp[inp<0])+31) / 32)
# all(out2==out)

inp = 0:128
inp = (-1027):1027
inp = (-127):127
inp = (-15):15
inp = (-31):31

out = trunc(inp / 32)
plot(inp, out, type='l', col='red4')

quant_div = function(inp, div) 
    sapply(inp, function(ii) sign(ii) * (1 + 2 * trunc(abs(ii) / div / 2)))

quant_div_one = function(inp, div) 
    sign(inp) * (1 + 2 * trunc(abs(inp) / div / 2))

quant_div = function(inp, div) 
    sapply(inp, function(ii) quant_div_one(ii, div))

quant_mul_one = function(inp, mul) 
    sign(inp) * (1 + 2 * trunc(abs(inp) * mul / 2))

quant_mul = function(inp, mul) 
    sapply(inp, function(ii) quant_mul_one(ii, mul))


div = 4
out = quant_div(inp, div)
outx = quant_mul(inp, 1/div)
stopifnot(outx==out)
out2 = sapply(inp, function(s) q_mul(s, 1024 / div))

plot(inp, out, type='l', col='red4') # , xlim=c(900,1040), ylim=c(50, 70))
lines(inp, inp/div, col = 'blue4')
plot(inp, out2, col = 'blue4', type = 'l')
data.frame(inp, out, out2)[which(out != out2), ]

#--
inp = (-127):127
q10_mul = 50
level = 2048/q10_mul
level
out = sapply(inp, function(s) quant_mul(s, q10_mul / 1024))
out2 = sapply(inp, function(s) q_mul(s, q10_mul))
data.frame(inp, out, out2)






