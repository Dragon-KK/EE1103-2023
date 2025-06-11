/**
 * EE23B135 Kaushik G Iyer
 * 20/10/2023
 * 
 * An exercise on interpolation (Newton's method and Lagrange's method)
 * We were asked to first sample the true function to create L (the polynomial created using the lagrange interpolation method)
 * Then sample this function with newton's method and find the final value at xp (provided in input)
 * 
 * Expected input... (2 Arguments)
 *  - order<int> Degree of the polynomial ( = number of points - 1)
 *  - xstart<double> Defines the range in which we wish to interpolate
 *  - xend<double> Defines the range in which we wish to interpolate
 *  - xp<double> The x value at which we wish to find the y for by interpolation
 * 
 * NOTE: The file is expected to have N+1 Columns (Corresponding to the N+1 coefficients) and N rows (Corresponding to the N linear equations)
 *
 * Outputs...
 *  - Final interpolated value (Basically  = N(L(f)))
 * 
 * NOTE: Both Newton and Lagrange create the same polynomiald
 * NOTE: Polynomials created again based on a polynomial created by interpolation just provides the same polynomial :)
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

struct Options{
    int order; // Degree of the polynomial ( = number of points - 1)
    double xstart; // Defines the range in which we wish to interpolate
    double xend; // Defines the range in which we wish to interpolate
    double xp; // The x value at which we wish to find the y for
};

struct Options getOptions(int argc, char** argv);

double lagrange(double * x, double* y, int n, double xp);
double newton(double * x, double* y, int n, double xp);

double function(double x);

int main(int argc, char** argv){
    struct Options options = getOptions(argc, argv);

    double* x = malloc((options.order + 1) * sizeof(double)); // Stores the x values of the samples of the true function
    double* y = malloc((options.order + 1) * sizeof(double)); // Stores the y values of the samples of the true function
    double* xl = malloc((options.order + 1) * sizeof(double)); // Stores the x values of `L` (the polynomial formed from the samples above)
    double* yl = malloc((options.order + 1) * sizeof(double)); // Stores the y values of `L` (the polynomial formed from the samples above)
    if (x == NULL || y == NULL || xl == NULL || yl == NULL){
        fprintf(stderr, "ERROR! Could not allocate enough memory for storing the points\n");
        abort();
    }

    for (int i = 0; i <= options.order; i++){ // Populate the x and y values from the true function
        x[i] = options.xstart + ((options.xend - options.xstart) * i / options.order);
        y[i] = function(x[i]);
    }

    // I find the value at (x<i> + x<i+1>) / 2 by interpolating (Thus creating my function `L`)
    for (int i = 0; i < options.order; i++){ // Populate the x and y values from the lagrange samples
        xl[i] = x[i] + ((options.xend - options.xstart) * 0.5 / options.order);
        yl[i] = lagrange(x, y, options.order + 1, xl[i]);        
    }
    // For the last point I just the same x<n> (Since x<n+1> is not defined)
    xl[options.order] = x[options.order];
    yl[options.order] = lagrange(x, y, options.order + 1, xl[options.order]);

    // Print value at xp for the function `L`
    printf("%f\n", newton(xl, yl, options.order + 1, options.xp));

    free(x);
    free(xl);
    free(y);
    free(yl);
}


/**
 * Interpolates using Newton's  method
 * Refer to chapra for the pseudo code
*/
double newton(double* x, double* y, int n, double xp){
    // They called this fdd in the pseudo code (It's basically stores the derivative type of thing that's used in Newton's method)
    // It basically allows us to calculate f[xi, xj...] relatively easily
    double** matrix = malloc(n * sizeof(double*));
    if (matrix == NULL){
        fprintf(stderr, "ERROR! Could not allocate enough memory to create matrix in newton's method\n");
        abort();
    }
    for (int i = 0; i < n; i++){ 
        matrix[i] = malloc(n * sizeof(double));
        if (matrix[i] == NULL) {
            fprintf(stderr, "ERROR! Could not allocate enough memory to create matrix in newton method\n");
            abort();
        }
        matrix[i][0] = y[i];
    }

    for (int j = 1; j < n; j++){
        for (int i = 0; i < n - j; i++){
            matrix[i][j] = (matrix[i + 1][j - 1] - matrix[i][j - 1]) / (x[i + j] - x[i]);
            // This is a slightly non trivial thing to understand, but the math checks out :)
        }
    }

    double xterm = 1; // Stores the (xp - xi) x ... part
    double yp = matrix[0][0]; // Stores the final answer
    for (int order = 1; order < n; order++){
        xterm *= (xp - x[order - 1]);
        yp += matrix[0][order] * xterm;
    }

    for (int i = 0; i < n; i++){ free(matrix[i]); }
    free(matrix);

    return yp;
}

/**
 * Interpolates using Lagrange's  method
 * Refer to chapra for the pseudo code
*/
double lagrange(double* x, double* y, int n, double xp){
    double sum = 0;
    for (int i = 0; i < n; i++){
        double product = y[i];
        for (int j = 0; j < n; j++){
            if (i == j) continue;
            product = product * (xp - x[j]) / (x[i] - x[j]);
            // It is very easy to see how this produces the formula given in Lagrange's method
        }
        sum += product;
    }
    return sum;
}

/**
 * The true function we use for this whole exercise
*/
double function(double x){
    return 1 / (1 + (25*pow(x, 2)));
}

/**
 * Parses the command line arguments to read the options
*/
struct Options getOptions(int argc, char** argv){
    if (argc < 5){
        fprintf(stderr, "ERROR! Expected atleast 4 arguments (order, xstart, xend, xp)\n");
        abort();
    }
    struct Options options;
    options.order = atoi(argv[1]);
    options.xstart = atof(argv[2]);
    options.xend = atof(argv[3]);
    options.xp = atof(argv[4]);

    // Remember atoi returns 0 if the argument was not an actual integer
    if (options.order <= 0){
        fprintf(stderr, "ERROR! Expected order (argv[1]) to be a natural number (that can be stored in an int)\n");
        abort();
    }

    if (options.xend <= options.xstart){
        fprintf(stderr, "ERROR! Expected xend (argv[3]) > xstart (argv[2])\n");
        abort();
    }

    // Ig we can allow xp not in between xstart and xend?
    return options;
}
