/**
 * EE23B135 Kaushik G Iyer
 * 24/10/2023
 * 
 * Uses cubic spline interpolation to interpolate values :)
 *  
 * Expected input... (2 Arguments)
 *  - order<int> Degree of the polynomial ( = number of points - 1)
 *  - xp<double> The x value at which we wish to find the y for by interpolation
 *  - xstart? <double> Defines the range in which we wish to interpolate
 *  - xend? <double> Defines the range in which we wish to interpolate
 * 
 * Outputs...
 *  - The interpolated value at xp
 * 
 * NOTE: When the x values are very high (Thus giving low values of coefficients), the gaussian elimination does not provide roots (Due to the tolerance set)
 * So either decrease the tolerance or change the function to give higher y values at bigger x to avoid this problem
 * 
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

// Used to define how close to 0 is 0 (Used in gaussian elimination)
#define TOLERANCE 0.0001

struct Options{
    int order; // The number of the polynomials ( = number of points - 1)
    double xstart; // Defines the range in which we wish to interpolate
    double xend; // Defines the range in which we wish to interpolate
    double xp; // The x value at which we wish to find the y for
};

struct Cubic{
    double a0; // Coefficient of constant
    double a1; // Coefficient of x term
    double a2; // Coefficient of x^2 term
    double a3; // Coefficient of x^3 term
};

bool findCubicEquations(double* x, double* y, int n, struct Cubic* cubics);

double function(double x);
struct Options getOptions(int argc, char** argv);
bool solve(double** coeffs, double* roots, int n);

int main(int argc, char** argv){
    struct Options options = getOptions(argc, argv);

    double* x = malloc((options.order + 1) * sizeof(double)); // Stores the x values of the samples of the true function
    double* y = malloc((options.order + 1) * sizeof(double)); // Stores the y values of the samples of the true function
    struct Cubic* cubics = malloc(options.order * sizeof(struct Cubic));

    if (x == NULL || y == NULL){
        fprintf(stderr, "ERROR! Could not allocate enough memory for storing the points\n");
        abort();
    }
    if (cubics == NULL){
        fprintf(stderr, "ERROR! Could not allocate enough memory for storing the coefficients of the cubic equations\n");
        abort();
    }

    for (int i = 0; i <= options.order; i++){ // Populate the x and y values from the true function
        x[i] = options.xstart + ((options.xend - options.xstart) * i / options.order);
        y[i] = function(x[i]);
    }

    if (findCubicEquations(x, y, options.order + 1, cubics)){
        int cubicIndex = (options.xp - options.xstart) * options.order / (options.xend - options.xstart);

        // Clamp it just in case xp is not in between start and end
        if (cubicIndex > options.order - 1) cubicIndex = options.order - 1;
        if (cubicIndex < 0) cubicIndex = 0;
        
        double val = cubics[cubicIndex].a0 + cubics[cubicIndex].a1*options.xp + cubics[cubicIndex].a2*pow(options.xp, 2) + cubics[cubicIndex].a3*pow(options.xp, 3);
        
        printf("%f\n", val);        
    }
    else{
        fprintf(stderr, "ERROR! Could not solve the system of equations to find the coefficients of the cubic polynomials\n");
    }

    free(cubics);
    free(x);
    free(y);
}


/**
 * Uses mafs to find the coefficients of the n-1 cubic equations
 * There are 4*(n - 1) unkowns
 * By saying f(xi) = yi we get 2*(n - 2) + 2 equations
 * By equation the slopes and concativity of adjacent polynomials we get another 2*(n-2) equations
 * Then for the last 2 equations, we assume the concativity of the edge polynomials to be 0 at the edge points
 * 
 * NOTE: This function returns false if the cubics could not be found (Due to the system of equations being unsolvable)
*/
bool findCubicEquations(double* x, double* y, int n, struct Cubic* cubics){
    double** matrix = malloc(sizeof(double*) * 4 * (n - 1));
    if (matrix == NULL){
        fprintf(stderr, "ERROR! Could not allocate enough memory to create matrix (To solve for coefficients)\n");
        abort();
    }
    for (int i = 0; i < 4*(n-1); i++){
        matrix[i] = malloc(sizeof(double) * (4*(n - 1) + 1));
        if (matrix[i] == NULL){
            fprintf(stderr, "ERROR! Could not allocate enough memory to create matrix (To solve for coefficients)\n");
            abort();
        }

        for (int j = 0; j < 4*(n-1) + 1; j++){
            matrix[i][j] = 0; // Initially set every coefficient to 0
        }
    }
    double *roots = malloc(sizeof(double) * 4 * (n - 1));
    if (roots == NULL){
        fprintf(stderr, "ERROR! Could not allocate enough memory to store the roots of the equation (To solve for coefficients)\n");
        abort();
    }

    int eqn = 0; // Keeps track of the equatino that we are currently writing
    
    // Use f(xi) = yi
    for (int i = 0; i < 4; i++){ matrix[eqn][i] = pow(x[0], i); }
    matrix[eqn][4 * (n - 1)] = y[0];
    eqn++;

    
    for (int p = 1; p < n-1; p++){
        for (int i = 0; i < 4; i++){ matrix[eqn][4*p - 4 + i] = pow(x[p], i); }
        matrix[eqn][4 * (n - 1)] = y[p];
        eqn++;
        for (int i = 0; i < 4; i++){ matrix[eqn][4*p + i] = pow(x[p], i); }
        matrix[eqn][4 * (n - 1)] = y[p];
        eqn++;
    }

    for (int i = 0; i < 4; i++){ matrix[eqn][4*(n - 1) - 4 + i] = pow(x[n-1], i); }
    matrix[eqn][4 * (n - 1)] = y[n-1];
    eqn++;

    // Equate slopes
    for (int p = 1; p < n-1; p++){
        // a1 + 2a2*x + 3a3*x^2 = a1 + 2a2*x + 3a3*x^2
        for (int i = 1; i < 4; i++){ matrix[eqn][4*p - 4 + i] = i * pow(x[p], i - 1); }
        for (int i = 1; i < 4; i++){ matrix[eqn][4*p + i] =  - i * pow(x[p], i - 1); }
        eqn++;
    }

    // Equate concativity
    for (int p = 1; p < n-1; p++){
        // 2a2 + 6a3*x = 2a2 + 6a3x
        for (int i = 2; i < 4; i++){ matrix[eqn][4*p - 4 + i] = i * (i - 1) * pow(x[p], i - 2); }
        for (int i = 2; i < 4; i++){ matrix[eqn][4*p + i] = - i * (i - 1) * pow(x[p], i - 2); }
        eqn++;
    }

    // Assume the concativity of the edge polynomials to be 0 at the edge points
    // 2a2 + 6a3*x = 0
    for (int i = 2; i < 4; i++){ matrix[eqn][i] = i * (i - 1) * pow(x[0], i - 2); }
    eqn++;
    for (int i = 2; i < 4; i++){ matrix[eqn][4*(n-1) - 4 + i] = i * (i - 1) * pow(x[n-1], i - 2); }
    eqn++;

    // Just for debugging (prints the system of equations)
    // for (int i = 0; i < 4*(n-1); i++){
    //     for (int j = 0; j < 4*(n - 1); j++){
    //         printf("(%f) + ", matrix[i][j]);
    //     }
    //     printf(" = %f\n", matrix[i][4*(n-1)]);
    // }

    bool succesful = solve(matrix, roots, 4*(n-1));
    for (int i = 0; i < 4*(n-1); i++){ free(matrix[i]); }
    free(matrix);

    for (int i = 0; i < n-1; i++){
        cubics[i].a0 = roots[4*i + 0];
        cubics[i].a1 = roots[4*i + 1];
        cubics[i].a2 = roots[4*i + 2];
        cubics[i].a3 = roots[4*i + 3];
    }

    free(roots);

    return succesful;
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
    if (argc < 3){
        fprintf(stderr, "ERROR! Expected atleast 2 arguments (order, xp, xstart?, xend?)\n");
        abort();
    }
    struct Options options;
    options.order = atoi(argv[1]);
    options.xp = atof(argv[2]);
    options.xstart = (argc>3)?atof(argv[3]):-2;
    options.xend = (argc>4)?atof(argv[4]):2;

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


// Taken from GaussianElimination.c (A previous assignment)

/**
 * Pls refer to chapra vro
*/
void pivot(double** coeffs, double* maxValueArray, int n, int k){
    int p = k;
    double big = fabs(coeffs[k][k] / maxValueArray[k]);
    double tmp;
    for (int i = k + 1; i < n; i++){
        tmp = fabs(coeffs[i][k] / maxValueArray[i]);
        if (tmp > big){
            big = tmp;
            p = i;
        }
    }

    if (p == k){ return; }

    for (int j = k; j < n + 1; j++){
        tmp = coeffs[p][j];
        coeffs[p][j] = coeffs[k][j];
        coeffs[k][j] = tmp;
    }
    tmp = maxValueArray[p];
    maxValueArray[p] = maxValueArray[k];
    maxValueArray[k] = tmp;    
}

/**
 * Refer to chapra vro
*/
bool eliminate(double** coeffs, double* maxValueArray, int n){
    for (int k = 0; k < n-1; k++){
        pivot(coeffs, maxValueArray, n, k);
        if (fabs(coeffs[k][k]/maxValueArray[k]) < TOLERANCE){
            return false;
        }

        for (int i = k + 1; i < n; i++){
            double factor = coeffs[i][k] / coeffs[k][k];
            for (int j = k; j < n + 1; j++){
                coeffs[i][j] -= factor*coeffs[k][j];
            }
        }
    }

    if (fabs(coeffs[n-1][n-1]/maxValueArray[n-1]) < TOLERANCE){
        return false;
    }
    return true;
}

/**
 * Refer to chapra lmao
*/
void substitute(double** coeffs, double* roots, int n){
    roots[n-1] = coeffs[n-1][n] / coeffs[n-1][n-1];
    for (int i = n-2; i >-1; i--){
        double sum = 0;
        for (int j = i + 1; j < n; j++){
            sum += coeffs[i][j] * roots[j];
        }
        roots[i] = (coeffs[i][n] - sum) / coeffs[i][i];
    }
}

/**
 * Finds the roots of the system of equations expressed by the coefficients using the Gauss Elimination Method
 * (Refer to Chapre pg 245 for pseudo code)
*/
bool solve(double** coeffs, double* roots, int n){
    double* maxValueArray = malloc(sizeof(double) * n);

    for (int i = 0; i < n; i++){
        maxValueArray[i] = fabs(coeffs[i][0]);
        for (int j = 1; j < n; j++){
            maxValueArray[i] = fmax(maxValueArray[i], fabs(coeffs[i][j]));
        }
    }

    bool succesful = eliminate(coeffs, maxValueArray, n);
    free(maxValueArray);

    if (succesful){
        substitute(coeffs, roots, n);
        return true;
    }

    fprintf(stderr, "ERROR! No unique solution for system of equations\n");
    return false;
}
