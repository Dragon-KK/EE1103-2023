/**
 * EE23B135 Kaushik G Iyer
 * 03/11/2023
 * 
 * Differential equations are based on the rotation of the Earth's magnetic field
 * Compares different ODE solvers (RK45, Euler, Huen)
 * Then calculates the R^2 difference  (in the theta vs phi graph)
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
 *  alpha delT Euler_R^2 Huen_R^2
 * 
 * NOTE: Certain values of deltaT cause the theta to get stuck below thetaMax, please use small deltaT
*/

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define SAVE_PLOT 0
#define RK45_PLOT_FILE_NAME "rk45.dat"
#define EULER_PLOT_FILE_NAME "euler.dat"
#define HUEN_PLOT_FILE_NAME "huen.dat"

#define GAMMA 1.76E11f
#define H 1E-10f
#define Hk 0

struct Options{
    double thetaStart;
    double thetaStop;
    double delT;
    double alpha;
};

void setOptions(int argc, char** argv, struct Options* options);

void savePlot(struct Options options, char* fileName, void odeSolver(double f(double, double, double), double h, double *x, double y, double alpha));

double thetaDifferential(double theta, double _phi, double alpha);
double phiDifferential(double _phi, double theta, double alpha);

// void rk45(double f(double, double), double h, double *x, double alpha, double *epsilon);
void rk45(double f(double, double, double), double h, double *x, double y, double alpha);
void euler(double f(double, double, double), double h, double *x, double y, double alpha);
void huen(double f(double, double, double), double h, double *x, double y, double alpha);


void main(int argc, char** argv)
{
    struct Options options;
    setOptions(argc, argv, &options);   
    
    if (SAVE_PLOT){
        savePlot(options, RK45_PLOT_FILE_NAME, rk45);
        savePlot(options, EULER_PLOT_FILE_NAME, euler);
        savePlot(options, HUEN_PLOT_FILE_NAME, huen);
    }

    double thetaRk, thetaEuler, thetaHuen;

    thetaRk = thetaEuler = thetaHuen = options.thetaStart;
    
    double RSS_Euler = 0; // The residual^2 for initial point is 0 only
    double RSS_Huen = 0; // The residual^2 for initial point is 0 only
    
    long long int count = 1; // Counts the number of points we are considering (Set to 1 as the initial point is also counted)

    double eulerMean = thetaEuler; // Initially we find the sum of all thetas then we divide by count later (Used in TSS)
    double huenMean = thetaHuen; // Initially we find the sum of all thetas then we divide by count later (Used in TSS)
    
    // NOTE: Since phi does not affect theta v t relation, we may just ignore it :)
    // Keep going until all the thetas reach thetaStop
    while(thetaRk < options.thetaStop || thetaEuler < options.thetaStop || thetaHuen < options.thetaStop)
    {
        rk45(thetaDifferential, options.delT, &thetaRk, 0, options.alpha);
        euler(thetaDifferential, options.delT, &thetaEuler, 0, options.alpha);
        huen(thetaDifferential, options.delT, &thetaHuen, 0, options.alpha);
        
        RSS_Euler += pow(thetaEuler - thetaRk, 2);
        RSS_Huen += pow(thetaHuen - thetaRk, 2);
        
        eulerMean += thetaEuler;
        huenMean += thetaHuen;

        count++;
    }

    eulerMean = eulerMean / count;
    huenMean = huenMean / count;

    // Reset the loop again
    thetaRk = thetaEuler = thetaHuen = options.thetaStart;
    
    double TSS_Euler = pow(thetaEuler - eulerMean, 2);
    double TSS_Huen = pow(thetaHuen - huenMean, 2);

    while(thetaRk < options.thetaStop || thetaEuler < options.thetaStop || thetaHuen < options.thetaStop)
    {
        rk45(thetaDifferential, options.delT, &thetaRk, 0, options.alpha);
        euler(thetaDifferential, options.delT, &thetaEuler, 0, options.alpha);
        huen(thetaDifferential, options.delT, &thetaHuen, 0, options.alpha);
        
        TSS_Euler += pow(thetaEuler - eulerMean, 2);
        TSS_Huen += pow(thetaHuen - huenMean, 2);
    }

    // R^2 = 1 - (RSS/TSS)
    double eulerR2 = 1 - (RSS_Euler / TSS_Euler);
    double huenR2 = 1 - (RSS_Huen / TSS_Huen);

    printf("%f %f %f %f\n", options.alpha, options.delT, eulerR2, huenR2);
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