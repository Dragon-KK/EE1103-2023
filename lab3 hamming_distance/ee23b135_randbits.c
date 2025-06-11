/**
 * EE23B135 Kaushik G Iyer
 * 25/08/2023
 * 
 * Creates random (ASCII) bit sequences and stores it in "randbits.txt"
 * AND?? calculates hamming distance between two input files
 * 
 * Expected input...
 * NOTE: When -h is used, the program does not create random bit sequence
 * 
 * For creating random bit
 *  -n {number of bits} // Provide length of random bit sequence required
 *  -s {seed} // Provide seed to set to srand
 *  -t // Use time as seed for srand NOTE: This is given preference over -s
 * 
 * For calculating hamming distance
 *  -h {file1} {file2} // The program will read from both file1 and file2 and calculate the hamming distance
 *      NOTE: If the data in the files are not of same size, stops calculating distance at the smaller file
 * 
 * Outputs:
 *  A file called "randbits.txt" in the same directory when -h is not set
 *  Prints `Hamming distance is {..}` when -h is set
 * 
 * NOTE: I am using an int to store hamming distance so beware of overflow for biiiiiiig files
 * 
*/ 
#include <stdbool.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define OUT_FILE_NAME "randbits.txt"

struct Options{
    int numberOfBits;
    
    bool useTimeAsSeed;
    int seed;
    
    bool calculateHammingDistance;
    char* hammingInput1;
    char* hammingInput2;
};

struct Options getOptions(int argc, char** argv);
void outputHammingDistance(char* path1, char* path2);

int main(int argc, char** argv){
    struct Options options = getOptions(argc, argv); // Parse command line options

    if (options.calculateHammingDistance){ // NOTE: If -h was set we don't generate the bit sequence
        outputHammingDistance(options.hammingInput1, options.hammingInput2);
        return 0;
    }

    // If we have to use the time flag, use the time otherwise use the seed set with -s
    int seed = options.useTimeAsSeed? time(NULL) : options.seed;
    if (seed != -1) srand(seed); // note that ff seed was not set (i.e. neither -s nor -t were used) we just don't set seed

    // We create the random sequence of bits only if we recieve N > 0
    if (options.numberOfBits > 0){
        char bits[options.numberOfBits + 1];

        for (int i = 0; i < options.numberOfBits; i++){
            bits[i] = '0' + rand() % 2; // bits[i] will be either '0' or '1'
        }
        bits[options.numberOfBits] = '\0'; // NULL temrinated string moment

        FILE* fp = fopen(OUT_FILE_NAME, "w");
        if (fp == NULL){
            printf("Could not open file to write ;-;\n");
            return 1;
        }
        else{
            fprintf(fp, "%s", bits);
            fclose(fp);
        }
    }
    else{
        printf("Please give an n>0 to create random bit sequence\n");
    }
}

/**
 * Reads data from two files and prints out the hamming distance between them
 * 
 * @param path1 Path to the first file
 * @param path2 Path to the second
*/
void outputHammingDistance(char* path1, char* path2){
    FILE* fp1 = fopen(path1, "r");
    FILE* fp2 = fopen(path2, "r");

    if (fp1 == NULL){
        printf("`%s` is not a valid file vro :(\n", path1);
        return;
    }
    if (fp2 == NULL){
        printf("`%s` is not a valid file vro :(\n", path2);
        return;
    }
    int hammingDistance = 0;
    
    char c1 = fgetc(fp1);
    char c2 = fgetc(fp2);

    // If strings are of unequal length we just stop at the shorter string :)
    while (c1 != EOF && c2 != EOF){
        if (c1 != c2){ // Hamming distance is basically the number of different `bits`
            hammingDistance++;
        }
        
        c1 = fgetc(fp1);
        c2 = fgetc(fp2); // Read character by character
    }

    printf("Hamming distance is %d\n", hammingDistance);
}


/**
 * Gets options by parsing the command line arguments :)
*/
struct Options getOptions(int argc, char** argv){
    struct Options options = {
        0, // numberOfBits
        false, -1, // useTimeAsSeed, seed 
        false, NULL, NULL // calculateHammingDistance, hammingInput1, hammingInput2
    };
    int c;
    
    while ((c = getopt(argc, argv, "n:s:th:")) != -1){
        switch (c)
        {
        case 'n': // Flag for setting number of bits
            options.numberOfBits = atoi(optarg);            
            break;

        case 's': // Flag for setting seed
            options.seed = atoi(optarg);
            break;

        case 't': // Flag for setting seed based on time
            options.useTimeAsSeed = true;
            break;

        case 'h': // Flag for calculating hamming distance
            options.calculateHammingDistance = true;
            options.hammingInput1 = optarg;
            
            // This part is a little messy, basically we try to read the next argument
            // It works bro full trust :)
            if (optind >= argc || argv[optind][0] == '-'){ 
                // If we reach the last argument or the next argument starts with a `-` (therefore its possibly a different flag)
                printf("We expect 2 file names after `-h` bro ;_;\n");
                abort();
            }
            options.hammingInput2 = argv[optind];            
            break;
        
        default:
            abort();
            break;
        }
    }
    return options;
}
