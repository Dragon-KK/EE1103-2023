/**
 * EE23B135 Kaushik G Iyer
 * 03/11/2023 + 19/11/2023
 * Version 2 ig?
 * 
 * Differential equations are based on the rotation of the Earth's magnetic field
 * Compares different ODE solvers (RK45, Euler, Huen)
 * Then calculates the R^2 difference (in the theta vs t graph)
 * 
 * References:
 *  - RK45 Code yeeted from https://web.ma.utexas.edu/CNA/cheney-kincaid/Ccode/CHP10/rk45.c
 *  - Physics copy pasted from https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=875251
 * 
 * Input:
 *  ./prog thetaStart thetaEnd alpha delT
 *  - thetaStart <double> The starting angle of the magnetic field
 *  - thetaEnd <double> The ending angle of the magnetic field
 *  - alpha <double> Basically a constant that defines how quickly the rotation occurs
 *  - delT <double> The shortest time interval in our 'simulation'
 * 
 * Output:
 *  alpha delT Euler_R^2 Huen_R^2 (If 'PRINT_TAU' is set to false)
 *  OR
 *  tau_rk_golden tau_euler tau_huen (If 'PRINT_TAU' is set to true)
 * 
 * NOTE: Certain values of deltaT cause the theta to get stuck below thetaMax, please use small deltaT
*/

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

#define SAVE_PLOT 0
#define PRINT_TAU 0
#define RK45_PLOT_FILE_NAME "rk45.dat"
#define EULER_PLOT_FILE_NAME "euler.dat"
#define HUEN_PLOT_FILE_NAME "huen.dat"

#define RK45_DELT 0.01f

#define GAMMA 1.76E11f
#define H 1E-10f
#define Hk 0

struct Options{
    double thetaStart;
    double thetaStop;
    double delT;
    double alpha;
};

struct Vector{
    double* arr;
    long long int MAX_SIZE;
    long long int size;
};

void setOptions(int argc, char** argv, struct Options* options);

double getGoldenRK45(struct Vector* vector, double t);

void populatePoints(struct Vector* vector, double thetaStart, double thetaStop, double delT, double alpha, void odeSolver(double f(double, double, double), double h, double *x, double y, double alpha));
void savePlot(struct Options options, char* fileName, void odeSolver(double f(double, double, double), double h, double *x, double y, double alpha));

double thetaDifferential(double theta, double _phi, double alpha);
double phiDifferential(double _phi, double theta, double alpha);

struct Vector* VEC__create(long long int max_size);
void VEC__delete(struct Vector* vector);
void VEC__append(struct Vector* vector, double value);
double VEC__sum(struct Vector* vector);

void rk45(double f(double, double, double), double h, double *x, double y, double alpha);
void euler(double f(double, double, double), double h, double *x, double y, double alpha);
void huen(double f(double, double, double), double h, double *x, double y, double alpha);


void main(int argc, char** argv)
{
    struct Options options;
    setOptions(argc, argv, &options);

    struct Vector* rk45Points = VEC__create(1000); // 1000 feels like a good starting size :)
    populatePoints(rk45Points, options.thetaStart, options.thetaStop, RK45_DELT, options.alpha, rk45);
    struct Vector* eulerPoints = VEC__create(1000); // 1000 feels like a good starting size :)
    populatePoints(eulerPoints, options.thetaStart, options.thetaStop, options.delT, options.alpha, euler);
    struct Vector* huenPoints = VEC__create(1000); // 1000 feels like a good starting size :)
    populatePoints(huenPoints, options.thetaStart, options.thetaStop, options.delT, options.alpha, huen);

    if (SAVE_PLOT){
        // The save plot function also finds the value of phi (Which is not done by the populatePoints function)
        // Additional note, the rk45 saved is of the 'golden' delT
        double tmp = options.delT;
        options.delT = RK45_DELT;
        savePlot(options, RK45_PLOT_FILE_NAME, rk45);
        options.delT = tmp;
        savePlot(options, EULER_PLOT_FILE_NAME, euler);
        savePlot(options, HUEN_PLOT_FILE_NAME, huen);
    }

    double eulerMean = VEC__sum(eulerPoints) / eulerPoints->size;
    double huenMean = VEC__sum(huenPoints) / huenPoints->size;

    double eulerTSS = 0;
    double eulerRSS = 0;
    for (long long int i = 0; i < eulerPoints->size; i++){
        eulerTSS += pow(eulerPoints->arr[i] - eulerMean, 2); // TSS = sum[ (yi - ymean)^2 ]
        eulerRSS += pow(eulerPoints->arr[i] - getGoldenRK45(rk45Points, i * options.delT), 2); // RSS = sum[ (yi - f(ti))^2 ]
    }
    double huenTSS = 0;
    double huenRSS = 0;
    for (long long int i = 0; i < huenPoints->size; i++){
        huenTSS += pow(huenPoints->arr[i] - huenMean, 2); // TSS = sum[ (yi - ymean)^2 ]
        huenRSS += pow(huenPoints->arr[i] - getGoldenRK45(rk45Points, i * options.delT), 2); // RSS = sum[ (yi - f(ti))^2 ]
    }

    double eulerR2 = 1 - (eulerRSS / eulerTSS);
    double huenR2 = 1 - (huenRSS / huenTSS);

    if (PRINT_TAU){
        printf("%f %f %f %f\n", options.alpha, (rk45Points->size - 1) * RK45_DELT, (eulerPoints->size - 1) * options.delT, (huenPoints->size - 1) * options.delT);
    }
    else{
        // The default output specified in the assignment
        printf("%f %f %f %f\n", options.alpha, options.delT, eulerR2, huenR2);
    }

    VEC__delete(rk45Points);
    VEC__delete(eulerPoints);
    VEC__delete(huenPoints);
}

/**
 * Gets the value of the 'golden' rk45 function at a given t
 * 
 * NOTE: We linearly interpolate to find the value of the function if we can't find the value
*/
double getGoldenRK45(struct Vector* points, double t){
    int i = (t / RK45_DELT);
    if (i + 1 >= points->size){
        // I assume that after theta reaches theta stop, it remains constant
        return points->arr[points->size - 1];
    }
    double lvalue = points->arr[i];   
    double rvalue = points->arr[i + 1];
    // // Just linear interpolation
    return lvalue + ((t - (RK45_DELT*i)) * (rvalue - lvalue)/RK45_DELT);
}

/**
 * Stores the values of theta using the ode
 * NOTE: We consider phi differential also just in case the user wants to use non zero Hk :)
*/
void populatePoints(struct Vector* vector, double thetaStart, double thetaStop, double delT, double alpha, void odeSolver(double f(double, double, double), double h, double *x, double y, double alpha)){
    double theta = thetaStart;
    VEC__append(vector, theta);
    while(theta < thetaStop){
        // NOTE: theta is independant of phi, therefor we don't need to update phi here :)
        odeSolver(thetaDifferential, delT, &theta, 0, alpha);
        VEC__append(vector, theta);
    }
}

/**
 * Saves the plot as theta goes from thetaStart to thetaEnd for a given odeSolver
 * NOTE: The plot saved might have different number of points saved for different odeSolvers (Since the loop here stops only when theta reaches thetaEnd)
 * NOTE: This doesn't cause problems as the plots saved are expected to only be used for graphing purposes (and not to find R^2, for that you will need to cutoff all graphs at the same time)
*/
void savePlot(struct Options options, char* fileName, void odeSolver(double f(double, double, double), double h, double *x, double y, double alpha)){
    double theta = options.thetaStart;
    double phi = 0;
    double time = 0;

    FILE* fd = fopen(fileName, "w");
    if (fd == NULL){
        fprintf(stderr, "ERROR! Could not open file `%s` to write plot\n", fileName);
        abort();
    }
 
    fprintf(fd, "%f %f %f\n", time, theta, phi);
    while(theta < options.thetaStop)
    {
        odeSolver(thetaDifferential, options.delT, &theta, phi, options.alpha);
        odeSolver(phiDifferential, options.delT, &phi, theta, options.alpha);
        time += options.delT;
        fprintf(fd, "%f %f %f\n", time, theta, phi);
    }
    fclose(fd);
}

/**
 * Allocates memory and sets the initial values for your vector :)
 * NOTE: You must free your vector (call VEC__delete)
*/
struct Vector* VEC__create(long long int max_size){
    struct Vector* vec = malloc(sizeof(struct Vector));
    vec->arr = malloc(max_size * sizeof(double));
    if (vec->arr == NULL){
        fprintf(stderr, "ERROR! Could not allocate memory for vector vro\n");
        abort();
    }
    vec->MAX_SIZE = max_size;
    vec->size = 0;
    return vec;
}

/**
 * Frees your vector :)
*/
void VEC__delete(struct Vector* vector){
    free(vector->arr);
    free(vector);
}

/**
 * Appends an item to our vector (Grows it if required)
*/
void VEC__append(struct Vector* vector, double value){
    if (vector->size == vector->MAX_SIZE){
        // Grow it by a factor of 2
        double* grown_arr = malloc(vector->MAX_SIZE * 2 * sizeof(double));
        if (grown_arr == NULL){
            fprintf(stderr, "ERROR! Could not allocate memory for vector vro\n");
            abort();
        }
        memcpy(grown_arr, vector->arr, vector->MAX_SIZE * sizeof(double));
        free(vector->arr);
        vector->arr = grown_arr;
        vector->MAX_SIZE = vector->MAX_SIZE * 2;
    }

    vector->arr[vector->size] = value;
    vector->size++;
}

/**
 * Finds the sum of all elements in the vector
*/
double VEC__sum(struct Vector* vector){
    double sum = 0;
    for (long long int i = 0; i < vector->size; i++){
        sum += vector->arr[i];
    }
    return sum;
}

/**
 * This method was taught in class
 * x(i + 1) = x(i) +x '(i) * h
*/
void euler(double f(double, double, double), double h, double *x, double y, double alpha){
    *x += f(*x, y, alpha) * h;
}

/**
 * This method was taught in class
 * x(guessed) = x(i) + x'(i) * h
 * x(i + 1) = x(i) + ((x'(i) + x'(guessed)) / 2) * h
*/
void huen(double f(double, double, double), double h, double *x, double y, double alpha){
    double xGuess = *x + f(*x, y, alpha) * h;
    *x += (f(*x, y, alpha) +  f(xGuess, y, alpha)) * h / 2;
}

/**
 * Refer to page 5 of https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=875251
*/
double thetaDifferential(double theta, double _phi, double alpha){
    // return GAMMA * alpha * ((H * sin(theta)) - (Hk * sin(theta) * cos(theta))) / (1 + (alpha*alpha));
    return (GAMMA * H) * alpha * sin(theta) / (1 + (alpha*alpha));
}

/**
 * Refer to page 5 of https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=875251
*/
double phiDifferential(double _phi, double theta, double alpha){
    // return -GAMMA * (H - (Hk * cos(theta))) / (1 + (alpha*alpha));
    return -(GAMMA*H)/ (1 + (alpha*alpha));
}

void setOptions(int argc, char** argv, struct Options* options){
    if (argc < 5){
        fprintf(stderr, "ERROR! Expected atleast 4 arguments (theta_start, theta_stop, alpha, del_t)\n");
        abort();
    }
    options->thetaStart = atof(argv[1]);
    options->thetaStop = atof(argv[2]);
    options->alpha = atof(argv[3]);
    options->delT = atof(argv[4]);

    if (options->thetaStart >= options->thetaStop){
        fprintf(stderr, "ERROR! Expected thetaStart (argv[1]) < thetaEnd (argv[2])\n");
        abort();
    }
    if (options->thetaStart <= 0){
        fprintf(stderr, "ERROR! Expected thetaStart (argv[1]) to be >0\n");
        abort();
    }
    if (options->thetaStop >= M_PI){
        fprintf(stderr, "ERROR! Expected thetaStop (argv[2]) to be < PI\n");
        abort();
    }
    if (options->delT <= 0){
        fprintf(stderr, "ERROR! Expected del_t (argv[4]) to be >0\n");
        abort();
    }
}

/**
 * Taken from https://web.ma.utexas.edu/CNA/cheney-kincaid/Ccode/CHP10/rk45.c
*/
void rk45(double f(double, double, double), double h, double* x, double y, double alpha)
{
    double c20 = 0.25, c21 = 0.25;
    double c30 = 0.375, c31 = 0.09375, c32 = 0.28125;
    double c40,c41, c42,c43;
    double c51, c52 = -8.0, c53, c54;
    double c60 = 0.5, c61, c62 = 2, c63, c64;
    double c65 = -0.275;
    // double a1, a2 = 0, a3, a4, a5 = -0.2;
    double b1, b2 = 0, b3, b4, b5= -0.18, b6;
    double F1, F2, F3, F4, F5, F6;
    // double x4;

    c40 = (double) 12/ (double) 13;
    c41 = (double) 1932/(double) 2197;
    c42 = (double) -7200/(double) 2197;
    c43 = (double) 7296/(double) 2197;
    c51 = c53 = (double) 439/ (double) 216;
    c54 = (double) -845/(double) 4104;
    c61 = (double) -8/(double) 27;
    c63 = (double) -3544/(double) 2565;
    c64 = (double) 1859/(double) 4104;
    // a1 = (double) 25/(double) 216;
    // a3 = (double) 1408/(double) 2565;
    // a4 = (double) 2197/(double) 4104;
    b1 = (double) 16/(double) 135;
    b3 = (double) 6656/(double) 12825;
    b4 = (double) 28561/(double) 56430;
    b6 = (double) 2/(double) 55;


    F1 = h * f(*x, y, alpha);
    F2 = h * f(*x + c21 * F1, y, alpha);
    F3 = h * f(*x + c31 * F1 + c32 * F2, y, alpha);
    F4 = h * f(*x + c41 * F1 + c42 * F2 + c43 * F3, y, alpha);
    F5 = h * f(*x + c51 * F1 + c52 * F2 + c53 * F3 + c54 * F4 , y, alpha);
    F6 = h * f(*x + c61 * F1 + c62 * F2 + c63 * F3 + c64 * F4 + c65 * F5, y, alpha);

    // x4 = *x + a1 * F1 + a3 * F3 + a4 * F4 + a5 * F5;
    *x += b1 * F1 + b3 * F3 + b4 * F4 + b5 * F5 + b6 * F6;
    // *epsilon = fabs(*x - x4);
}