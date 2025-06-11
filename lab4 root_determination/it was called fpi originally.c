/**
 * EE23B135 Kaushik G Iyer
 * 01/09/2023
 * 
 * Figures out the root of an already predetermined function (between predetermined xl and xu)
 * 
 * Expected input...
 * One argument either `1` or `2`
 * If the argument is `1` then the root is calculatd using the Bisection method
 * If the argument is `2` then the root is calculated using the False Position method
 *
 * Outputs:
 *  The calculated value of the root
 * 
 * NOTE: Link with -lm if required
 * 
*/ 
#include <string.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define BISECTION_X_START 0.0f
#define BISECTION_X_END 1.0f
#define BISECTION_ERROR 0.001f // i.e 10%

#define FALSE_POSITION_X_START 0.0f
#define FALSE_POSITION_X_END 1.0f
#define FALSE_POSITION_ERROR 0.002f // i.e 0.2%

double findRootByBisection(double xl, double xu);
double findRootsByFalsePosition(double xl, double xu);

int main(int argc, char const *argv[])
{
    if (argc < 2){
        printf("Please provide method to use vro\n");
        return 1;
    }
    
    if (!strcmp(argv[1], "1")){
        printf("%f\n", findRootByBisection(BISECTION_X_START, BISECTION_X_END));
    }
    else if (!strcmp(argv[1], "2")){
        printf("%f\n", findRootsByFalsePosition(FALSE_POSITION_X_START, FALSE_POSITION_X_END));        
    }
    else{
        printf("Unrecognized method! Please put either `1` for bisection or `2` for false position");
        return 1;
    }
    return 0;
}

double function(double x){
    
    return powl(M_E, -x) - x;
}

double findRootByBisection(double xl, double xu){
    double xr = (xl + xu) / 2; // NOTE: I am using the same notation as the textbook :)
    double currError = fabs(1 - (xl/xr));
    double nxr; // This will store the new xr value (Calculated in the loop)
    double mult; // Will store the sign of f(xr)*f(xl)

    if (function(xl) == 0){return xl;} // If we dont check this and f(xl) == 0, our logic breaks (since we return xr)
    
    while (currError > BISECTION_ERROR){ 
        mult = function(xr) * function(xl);
        if (mult > 0){ xl = xr; } // The root is between xr and xu
        else if (mult < 0){ xu = xr; } // The root is between xl and xr
        else {return xr;} // In case we somehow found the exact root

        nxr = (xl + xu) / 2; // Set the new xr to the mean
        currError = fabs(1 - (xr/nxr)); // Calculate the error
        xr = nxr; // Update the xr
    }

    return xr;
}

double findRootsByFalsePosition(double xl, double xu){
    double fxl = function(xl); // Stores f(xr)
    double fxu = function(xu); // Stores f(xu)
    double fxr; // Will store f(xr)

    if (fxl == 0){return xl;} // If we dont check this and f(xl) == 0, our logic breaks (since we return xr)

    double xr = xu - (fxu*(xu - xl)/(fxu - fxl)); // Refer to textbook for details on this formula
    double currError = fabs(1 - (xl/xr));
    double nxr; // Will store the newly calculated xr
    double mult; // Will store the sign of f(xr) * f(xu)

    while (currError > FALSE_POSITION_ERROR){
        fxr = function(xr);
        mult = fxr*fxl;
        if (mult > 0){ // Roots are between xr and xu
            xl = xr;
            fxl = fxr;
        }
        else if (mult < 0){ // Roots are between xl and xu
            xu = xr;
            fxu = fxr;
        }
        else{ // The root is xr
            return xr;
        }
        nxr = xu - (fxu*(xu - xl)/(fxu - fxl)); // Calculate the new xr value (Refer to book for formula)
        currError = fabs(1 - (xr/nxr));
        xr = nxr; // Store the xr value
    }
    return xr;
}