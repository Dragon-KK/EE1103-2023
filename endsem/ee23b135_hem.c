/**
 * EE23B135 Kaushik G Iyer
 * 14/11/2023
 * 
 * Outputs the Mth value of the Hemachandra sequence
 * Implements a recursive function the find the 'nth' number of the Hemachandra (Fibonacci) sequence
 * 
 * Input:
 *  M <int> (Specifies the term required)
 * 
 * Outputs:
 *  hem(M) <int> (The Mth integer in the hemachandra sequence)
 * 
 * NOTE: The program expected a 1-indexed index
 * i.e. ./hem 1 -> 0
 * ./hem 2 -> 1 
 * And so on
*/ 

#include <stdio.h>
#include <stdlib.h>

// Do I reallly need a struct for just one input?
// No... But I like consistency :)

struct Options{
    int M;
};

void setOptions(struct Options* options, int argc, char** argv);
int hem(int n);


int main(int argc, char** argv){
    struct Options options;
    setOptions(&options, argc, argv);

    printf("%d\n", hem(options.M));
}

/**
 * Returns numbers from the hemachandra sequence (basically the fibonacci sequence)
 * NOTE: Care for overflows
 * NOTE: If you are calling hem(n) multiple times, we must try and cache our results in perhaps an arrray :)
*/
int hem(int n){
    if (n == 1) {
        return 0;
    }
    if (n == 2){
        return 1;
    }
    return hem(n-1) + hem(n-2);
}

/**
 * Parses the command line arguments to read the options
*/
void setOptions(struct Options* options, int argc, char** argv){
    if (argc < 2){
        fprintf(stderr, "ERROR! Expected 1 argument M (int)\n");
        abort();
    }
    options->M = atoi(argv[1]);

    if (options->M <= 0){
        fprintf(stderr, "ERROR! Expected M (argv[1]) to be a positive integer\n");
        abort();
    }
}