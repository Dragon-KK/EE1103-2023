/**
 * EE23B135 Kaushik G Iyer
 * 13/10/2023
 * 
 * Figures out the roots of n dimensional series of linear equations
 * More of just an exercise in converting pseudo code that somebody else wrote to code
 * Refer to chapra Guass Elimination for the pseudo code
 * 
 * Expected input... (2 Arguments)
 *  - file_name <string> (The name of the file with the inputs)
 *  - N <int> (The number of variables)
 * 
 * NOTE: The file is expected to have N+1 Columns (Corresponding to the N+1 coefficients) and N rows (Corresponding to the N linear equations)
 *
 * Outputs...
 *  - Space separated values of xi <float>
 * 
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define TOLERANCE 0.0001

struct Options{
    int n; // Number of variables
    char* inputFileName; // The name of the file with the matrix
};

bool read(double** coeffs, char* inputFile, int n);
bool solve(double** coeffs, double* roots, int n);

struct Options getOptions(int argc, char** argv);

int main(int argc, char** argv){
    struct Options options = getOptions(argc, argv);

    double** matrix = malloc(sizeof(double*) * options.n);
    for (int i = 0; i < options.n; i++){
        matrix[i] = malloc(sizeof(double) * (options.n + 1));
    }
    double *roots = malloc(sizeof(double) * options.n);

    bool succesful = read(matrix, options.inputFileName, options.n);
    succesful = succesful && solve(matrix, roots, options.n);

    for (int i = 0; i < options.n; i++){ free(matrix[i]); }
    free(matrix);

    if (succesful){
        for (int i = 0; i < options.n - 1; i++){
            printf("%f ", roots[i]);
        }
        printf("%f\n", roots[options.n - 1]);
    }

    free(roots);

    if (!succesful){
        fprintf(stderr, "ERROR! Could not find roots\n");
        abort();
    }
}

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

/**
 * Reads the matrix of coefficients from the input file
*/
bool read(double** coeffs, char* inputFile, int n){
    FILE* fd = fopen(inputFile, "r");
    if (fd == NULL){
        fprintf(stderr, "ERROR! Could not open file to read coefficients\n");
        return false;
    }

    float tmp;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n + 1; j++){
            if (fscanf(fd, "%f ", &tmp) != 1){
                fprintf(stderr, "ERROR! Could not read Coefficients[%d][%d]\n", i, j);
                return false;
            };
            coeffs[i][j] = tmp;
        }
    }
    return true;
}


/**
 * Parses the command line arguments to read the options
*/
struct Options getOptions(int argc, char** argv){
    if (argc < 3){
        fprintf(stderr, "ERROR! Expected atleast 2 arguments (fileName, N)\n");
        abort();
    }
    struct Options options;
    options.inputFileName = argv[1];
    options.n = atoi(argv[2]);

    // Remember atoi returns 0 if the argument was not an actual integer
    if (options.n <= 0){
        fprintf(stderr, "ERROR! Expected N (argv[2]) to be a natural number\n");
        abort();
    }

    return options;
}
