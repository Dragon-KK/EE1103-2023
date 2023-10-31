/**
 * EE23B135 Kaushik G Iyer
 * 18/08/2023
 * 
 * Calculates sin(x) using a very naive method
 * 
 * Expected method of calling the program:
 * ./sin {precision} {x}
 * Where precision is the number of terms used from the taylor expansion of sin(x) and x is the value of which sin is found
 * 
 * Outputs {sin(x)},{sin(x) - trueSin(x)}
 * 
 * NOTE: Run gcc with -lm (to link with the math.h library)
*/ 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/** 
 * Calculates the factorial vro :)
 * 
 * @param n The number whose factorial is returned
*/
long long int factorial(long long int n){ return (n == 1)? 1 : n * factorial(n - 1); }
long double mySin(long double x, int precision);

int main(int argc, char** argv){
    if (argc < 3) {
        printf("Please provide {precision} and {x} vro\n");
        return 1;
    }
    int precision = atoi(argv[1]);
    float x = atof(argv[2]);

    long double mySinOfX = mySin(x, precision);
    
    printf("%.3Lf,%.3Lf\n", mySinOfX, sinl(x) - mySinOfX);
}

/**
 * Calculates the value of sin(x) using the taylor series
 * @note sin(x) = x - (x^3/3!) + (x^5/5!) - (x^7/7!) ...
 * 
 * @param x The number whose sin is calculated
 * @param precision The number of terms used from the taylor expansion
*/
long double mySin(long double x, int precision){
    long double ans = 0;
    for (int i=0; i<precision; i++){
        if (i%2 == 0) ans += (powl(x, 2*i + 1)/factorial(2*i + 1));
        else ans -= (powl(x, 2*i + 1)/factorial(2*i + 1));
    }
    return ans;
}