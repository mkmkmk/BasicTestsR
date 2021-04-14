#include <Rcpp.h>
using namespace Rcpp;

#define Q15 15

// [[Rcpp::export]]
double mulQ15(double a, double b)
{
  int qa = a * (1 << Q15);
  int qb = b * (1 << Q15);
  
  //int res = ((int64_t)qa * qb) >> Q15;
  int res = (qa * qb) >> Q15;
  
  return (double)res / (1 << Q15);
  
}

// [[Rcpp::export]]
double divQ15(double a, double b)
{
  int qa = a * (1 << Q15);
  int qb = b * (1 << Q15);
  
  //int res = ((int64_t)qa << Q15) / qb;
  int res = (qa << Q15) / qb;
  
  
  return (double)res / (1 << Q15);
  
}




/*** R

2^15
mulQ15(.25, 4)
mulQ15(.25, 8)

mulQ15(8, 8)

(2^15 - 1) * 2 == mulQ15(2^15 - 1, 2)

divQ15(1, 8)
divQ15(100, 8)

*/