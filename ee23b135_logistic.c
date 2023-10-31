/**
 * EE23B135 Kaushik G Iyer
 * 08/09/2023
 * 
 * Plots the logistic map (xt+1 = r * xt(1- xt)) for different values of r
 * 
 * Input:
 *  ./prog {iterations?} {step?} {error %?}
 * 
 * Where iterations is the maximum number of iterations made for a given `r`
 * Step is the difference between subsequent `r`s
 * Error % is the range within which numbers are considered to be the `same`
 * 
 * Outputs:
 *  A text file `points.txt` containing all the points present on the graph
 * 
 * 
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <unistd.h>

struct Options{
    double start; // The start value of r
    double end; // The end value of r
    double step; // The step value of r


    int points; // Max number of iterations per points
    int convergenceIterations; // Number of iterations to let the points settle

    double x0; // Initial value of x (x0)

    double tolerance; // If two consecutive terms are within this tolerance range, we will stop iterations
};

struct Options options = {
    0.0f, 4.0f, 0.005f,
    500, 500,
    0.5f,
    0.1
};

void getOptions(int argc, char** argv);

double next(double xt, double r){
    return r*xt*(1-xt);
}

bool errorIsAcceptable(double xt, double xold){
    return (fabs(xt - xold) / xold) < options.tolerance;
}

void storeValues(FILE* fd, double r){
    double xt = options.x0;
    int iterCount = 0;
    for (int i=0; i<options.convergenceIterations; i++){ // The values usually take a while to become little `clean` (this prevents some random offshoots that arise in the beginning)
        xt = next(xt, r);
        if (xt == 0){ // Once it reaches 0, it will always stay at 0
            fprintf(fd, "%.8f %.8f\n", r, 0);
            return;
        }
    }
    double xold = xt;
    xt = next(xt, r);

    fprintf(fd, "%.8f %.8f\n", r, xt);

    // Keep checking for more points as long as the new points we see are satisfactory (not too close to old value) and we have not yet reached max iteration count (i.e options.points)
    while ((!errorIsAcceptable(xt, xold)) && ++iterCount < options.points){
        xold = xt;
        xt = next(xt, r);
        fprintf(fd, "%.8f %.8f\n", r, xt);
    }
}


int main(int argc, char** argv){
    getOptions(argc, argv);
    float r = options.start;

    FILE* file = fopen("points.txt", "w");
    if (file == NULL){
        printf("Could not open the file vro ;-;\n");
        return 1;
    }

    while (r <= options.end){
        storeValues(file, r);
        r += options.step;        
    }

    fclose(file);
}


/**
 * Gets options by parsing the command line arguments :)
*/
void getOptions(int argc, char** argv){
    if (argc>1)options.points = atoi(argv[1]);
    if (argc>2)options.step = atof(argv[2]);
    if (argc>3)options.tolerance = atof(argv[3]) / 100;   

    // int c;
    
    // while ((c = getopt(argc, argv, "f:l:s:p:c:x:t:")) != -1){
    //     switch (c)
    //     {
    //     case 'f': // Start of r
    //         options.start = atof(optarg);            
    //         break;

    //     case 'l': // End of r
    //         options.end = atof(optarg);
    //         break;

    //     case 's': // Step size of r
    //         options.step = atof(optarg);
    //         break;

    //     case 'p': // Number of iterations
    //         options.points = atoi(optarg);         
    //         break;

    //     case 'c': // Number of iterations to allow the values to stabilize
    //         options.convergenceIterations = atoi(optarg);         
    //         break;
        
    //     case 'x': // Starting value (x0)
    //         options.x0 = atof(optarg);
    //         break;

    //     case 't':
    //         options.tolerance = atof(optarg) / 100;
    //         break;

    //     default:
    //         abort();
    //         break;
    //     }
    // }
}
