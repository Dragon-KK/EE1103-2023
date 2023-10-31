/**
 * EE23B135 Kaushik G Iyer
 * 27/10/2023
 * 
 * Uses gnuplot to figure out goodness of fit (Fitting data created from a lorentzian to a gaussian)
 * 
 * Input:
 *  ./prog {N} {sigma_N} {x_domain?}
 * 
 * Where:
 *  - N is the number of points
 *  - sigma_N is the standard deviation of the noise (the noise that is added on the lorentzian)
 *  - x_domain defines the domain of the graph (x values will range from -x_domain to +x_domain)
 * 
 * Outputs:
 *  {sigma_N} {amplitude (from fit data)} {sigma_G (from fit data)} {R^2 (Goodness of fit)}
 *  
*/ 
#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define GRAPH_TEMP_FILE_NAME "graph.dat"
#define OUT_TEMP_FILE_NAME "output.dat"
#define MAX_FLOAT_LEN 1000


struct Options{
    int N; // The number of points
    double sigma_N; // The standard deviation of the base noise
    double gap; // Scale is determined based on N and graph domain (It is the x gap between two points in the graph)
    double domain; // The graph is platted between -domain and domain
};

double randn(double sigma);
double function(double x);
void setOptions(int argc, char** argv, struct Options* options);

// Scale the noise to be bw -3sig and 3sig (sig taken from the input
int main(int argc, char** argv){
    struct Options options;
    setOptions(argc, argv, &options);
    
    FILE* fd = fopen(GRAPH_TEMP_FILE_NAME, "w");
    if (fd == NULL){
        fprintf(stderr, "ERROR! Could not open file to write\n");
        abort();
    }

    for (int xi = 0; xi < options.N; xi++){
        double x = -options.domain + (xi * options.gap);
        double y = randn(options.sigma_N);
        y += function(x);
        fprintf(fd, "%f %f\n", x, y);
    }
    fclose(fd);

    FILE* pipe = popen("gnuplot -persistent", "w");
    if (pipe == NULL){
        fprintf(stderr, "ERROR! Could not open file to write\n");
        abort();
    }

    // This makes the output from any `print` statement in gnuplot save itself to `OUT_TEMP_FILE_NAME`
    fprintf(pipe, "set print \"%s\"\n", OUT_TEMP_FILE_NAME);

    // This prevents the gnuplot fit function from printing stuff to stdout
    fprintf(pipe, "set fit quiet\n");

    // R^2 of our fit is basically 1 - (RSS/TSS)
    // Here TSS is the total sum of squares = Sum of (y - mean)^2
    // RSS is the sum of the squares of the residuals = Sum of (y - f(xi))^2

    // Now gnuplot provides a variable FIT_STDFIT, which is equal to sqrt((sum of (y - f(xi))^2) / N)
    // Note that when we fit the points to a constant function, it will make the constant function the mean
    // Therefore sqrt((sum of (y - f(xi))^2) / N) = sqrt((sum of (y - mean)^2) / N) = sqrt(TSS/N)

    // Also notice the FIT_STDFIT on the gaussian is equal to sqrt(RSS/N)

    // Now also notice that the ratio required RSS/TSS is then just the FIT_STDFIT(g)**2 / FIT_STDFIT(constant)**2
    // (Since the N parts get cancelled out)

    fprintf(pipe, "m(x) = m\n");
    fprintf(pipe, "fit m(x) '%s' via m\n", GRAPH_TEMP_FILE_NAME);
    fprintf(pipe, "print FIT_STDFIT**2\n"); // This is = the TSS/N

    fprintf(pipe, "g(x) = A*exp(-x**2 / (2*s*s))\n");
    fprintf(pipe, "fit g(x) '%s' via A,s\n", GRAPH_TEMP_FILE_NAME);
    fprintf(pipe, "print FIT_STDFIT**2\n"); // This is =  RSS/N

    fprintf(pipe, "print A\n"); // The A of the fitted gaussian
    fprintf(pipe, "print s\n"); // The sigma_g of the fitted gaussian
    pclose(pipe);
    
    fd = fopen(OUT_TEMP_FILE_NAME, "r");
    if (fd == NULL){
        fprintf(stderr, "ERROR! Could not open file to read output of gnuplot\n");
        abort();
    }

    // For some reason fscanf refused to read the floats, so I am using fgets along with atof to read the values
    char buff[MAX_FLOAT_LEN];

    fgets(buff, MAX_FLOAT_LEN, fd);
    double tss = atof(buff); // Actually TSS/N but the N gets cancelled out anyways
    fgets(buff, MAX_FLOAT_LEN, fd);
    double rss = atof(buff); // Actually RSS/N but the N gets cancelled out anyways
    fgets(buff, MAX_FLOAT_LEN, fd);
    double A = atof(buff);
    fgets(buff, MAX_FLOAT_LEN, fd);
    double sigma_g = atof(buff);

    double R2 = 1 - (rss/tss); // R^2 = 1 - (RSS/TSS)
    printf("%f %f %f %f\n", options.sigma_N, A, sigma_g, R2);    
}

/**
 * The functions used to create the data
 * NOTE: Doing the same exercise for the rayleigh function is as simple as changing this function :)
*/
double function(double x){
    return 1/(1+(25*x*x));
}

/**
 * Returns gaussian noise with standard deviation sigma
*/
double randn(double sigma){
    // https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform

    double u1 = (double)rand()/(double)RAND_MAX;
    double u2 = (double)rand()/(double)RAND_MAX;;

    double res = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2) * sigma;

    if (fabs(res) < 6*sigma){
        return res;
    }

    return (res>0)? 6*sigma : -6*sigma;
}

/**
 * Parses the command line arguments to read the options
*/
void setOptions(int argc, char** argv, struct Options* options){
    if (argc < 3){
        fprintf(stderr, "ERROR! Expected atleast 2 arguments (N, sigma_N)\n");
        abort();
    }
    options->N = atoi(argv[1]);
    if (options->N <= 1){
        fprintf(stderr, "ERROR! Expected N (argv[1]) to be atleast 2\n");
        abort();
    }
    options->sigma_N = atof(argv[2]);
    if (options->sigma_N < 0){
        fprintf(stderr, "ERROR! Expected sigma_N (argv[2]) to be non negative\n");
        abort();
    }
    
    options->domain = fabs((argc > 3)? atof(argv[3]) : 2.0f);
    options->gap = 2*options->domain / (options->N - 1);
}