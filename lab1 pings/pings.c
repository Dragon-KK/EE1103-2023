/**
 * EE23B135 Kaushik G Iyer
 * 
 * Reads from a file (name provided as argument) that follows a specific format
 * The program then calculates the mean and standard deviation of the data
 * 
 * NOTE: Run gcc with -lm (to link with the math.h library)
 * NOTE: Format expected:
 *  - 8th column (space separated) must contain the relevant "ping info"
 *  - "ping info" must consist of 5 non space characters followed by a float
 *      eg: `time=123.3` is valid (and is expected) 
*/ 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#pragma message "WARNING! This code only works if the data is cleaned (all the lines in the file contain data in the format descripted in the header)"

int main(int argc, char** argv){
    if (argc < 2){
        printf("Please give file name vro\n");
        return 1;
    }
    char* fileName = argv[1];    

    FILE* file = fopen(fileName, "r");

    if (file == NULL){
        printf("Please give valid file name vro\n");
        return 1;
    }

    float sum = 0;
    float sqsum = 0;
    int lines = 0;

    // Potential location of error (if time stamp is longer than 20 characters)
    char pingDat[20]; // Buffer to store data from scanf

    while (fscanf(file, "%*s %*s %*s %*s %*s %*s %*s %s %*s ", pingDat) == 1){ // The second last column contains the "ping info"
        float ping = atof(pingDat + 5); // We have to skip the first 5 characters of the timestamp (basically right after the `time=` part)
        
        sum += ping;
        sqsum += ping*ping;
        lines++;
    }

    // The formulae were provided by my seatmate :)
    printf("mean is %f | stdev is %f\n", sum/lines, sqrtf((sqsum/lines) - (sum/lines)*(sum/lines)));    
}