#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double floormod(double a, double b) 
{
  return a - b * floor(a / b);
}

// [[Rcpp::export]]
NumericVector floormodtb(NumericVector a, double b) 
{
    NumericVector outp = NumericVector(a.size());
    for(int i = 0; i < a.size(); i++)
    {
        outp[i] = floormod(a[i], b);
    }
    return outp;
}


// [[Rcpp::export]]
double floormod2(double a, double b) 
{
    return a < 0 ? fmod(a, b) + b : fmod(a, b);
}

// [[Rcpp::export]]
NumericVector floormodtb2(NumericVector a, double b) 
{
    NumericVector outp = NumericVector(a.size());
    for(int i = 0; i < a.size(); i++)
    {
        outp[i] = floormod2(a[i], b);
    }
    return outp;
}


/*** R

x = seq(-7, 7, 0.01)

# floormodtb = function(a,b)sapply(a,function(aa)floormod(aa, b))
# floormodtb2 = function(a,b)sapply(a,function(aa)floormod2(aa, b))

plot(x, floormodtb(x, 3), type = 'l', col = 'red', ylim = c(-4, 4))
lines(x, (x %% 3), type = 'l', col = 'blue')

library(microbenchmark)
microbenchmark(floormodtb(x,3), times =  10000)
microbenchmark(floormodtb2(x,3), times = 10000)

# czyli błędnie działa i zysku na czasie nie ma
all(floormodtb(x, 3) == (x %% 3))
all(floormodtb2(x, 3) == (x %% 3))

*/

