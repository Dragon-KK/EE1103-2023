#include <stdio.h>
#include <stdlib.h>
// LU Decomposition
// Calculate a0 and a1
int *f(){
    int c[1];
    *c=5;
    return c;
}
int main(){
    int *ptr = f();
    printf("%d",*ptr);
}

// Interpolation vs Curve fitting
// (cubic) Spline interpolation moment

// The following is pretty much just taylor series written in a different way
// Play with it a bit to get completely satisfied
// f(x) = f(x0) + (x-x0)f[x1, x0] + (x - x0)(x - x1)f[x2, x1, x0]
// f[x2, x1, x0] = (f[x2, x1] - f[x1, x0]) / (x2 - x0)
// Now we may construct the entire polynomial

// Lagrange polynomial
// sum<0, n>  [ Li(x)*f(xi) ]
// Where Li(x) = product<0, n> (x - xj) / (xi - xj)

// Quadratic splines

// Unkowns: 4*n - 4
// n-2 + 2*(n-2) + 2
// 3n - 4


// Take a truncated lorentzian with amplitude noise,
// Fit it to f(x) = 1/(1 + 25x^2)
// Use the gnuplot `fit` function to get R^2
// Then we can plot R^2 vs noise
// Fit it to a gaussian then find the sigma (in gaussian formula) e^-x^2/2sig

// Write a c program to call gnuplot and extract the fitted values

// Me when I have to do EE1103 :<


// ode
