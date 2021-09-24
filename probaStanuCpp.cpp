#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

int state = 0;

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) 
{
    state++;
    printf("state = %d \n", state);
    return x * 2;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
if (F)
{
    timesTwo(42)
    timesTwo(42)
    timesTwo(42)
    timesTwo(42)
}
*/
