/**
 * EE23B135 Kaushik G Iyer
 * 15/09/2023
 * Quiz 1 (Question 1)
 * 
 * Creates to random bit sequences of size N and M (N > M) and then figures out the minimum hamming distance we can get and the location (1 indexed) of where it can be found
 * 
 * Expected input... (4 Arguments)
 *  - N (An integer that is greater than 0), the size of the larger bit sequence
 *  - NSeed (An integer), the seed to use while generating the larger bit sequence
 *  - M (An integer that is greater than 0), the size of the smaller bit sequence
 *  - MSeed (An integer), the seed to use while generating the smaller bit sequence
 * 
 * NOTE: N must be strictly greater than M (As specified in the assignment)
 *
 * Outputs...
 * {location of first minimum (1 indexed)} {minimum hamming distance}
 * 
*/

#include <stdio.h>
#include <stdlib.h>

#define N_BITS_FILE_NAME "Nfile.dat"
#define M_BITS_FILE_NAME "Mfile.dat"

struct Options{
    int n; // Size of the larger buffer
    int nSeed; // Seed to use while creating the larger buffer
    int m; // Size of the smaller buffer
    int mSeed; // Seed to use while creating the smaller buffer
};

struct Options getOptions(int argc, char** argv);
void storeRandomBitsIntoBuffer(int* buffer, int size, int seed);
void writeBufferToFile(int* buffer, int size, const char* fileName);
void outputMinimumHammingDistance(int* nbits, int n, int* mbits, int m);
int getHammingDistance(int* buff1, int* buff2, int size);

int main(int argc, char** argv){
    struct Options options = getOptions(argc, argv);

    int nbits[options.n]; // Create a buffer of N integers
    int mbits[options.m]; // Create a buffer of M integers
    
    // Populate our buffers
    storeRandomBitsIntoBuffer(nbits, options.n, options.nSeed);
    storeRandomBitsIntoBuffer(mbits, options.m, options.mSeed);

    // Write the value of the buffers to their respective files
    writeBufferToFile(nbits, options.n, N_BITS_FILE_NAME);
    writeBufferToFile(mbits, options.m, M_BITS_FILE_NAME);

    // This one is kind of self descriptive ngl
    outputMinimumHammingDistance(nbits, options.n, mbits, options.m);
}

/**
 * Basically moves the smaller array like a slider over the larger array, and finds where the hamming distance is the minimum
*/
void outputMinimumHammingDistance(int* nbits, int n, int* mbits, int m){
    // First assume that the smallest hamming distance was when our slider was initially at index 0
    int minHammingDistance = getHammingDistance(nbits, mbits, m);
    int minHammingDistanceLocation = 0;
    nbits++;

    // i is basically the 0 indexed starting location of our slider
    for (int i=1; i <= n-m; i++){
        int hammingDistance = getHammingDistance(nbits, mbits, m); // Get the hamming distance (Imagine the slider to be moved)
        
        // If we find a smaller hammingDistance, update the location
        if (hammingDistance < minHammingDistance){ // Since we require the first instance of minimum, we use `<`
            minHammingDistance = hammingDistance;
            minHammingDistanceLocation = i;
        }
        nbits++; // Move the slider by 1
    }

    // NOTE: We do minHammingDistanceLocation + 1 since the location that we output is expected to be 1 indexed
    printf("%d %d\n", minHammingDistanceLocation + 1, minHammingDistance);
}

/**
 * Calculates hamming distance between two equally sized buffers
 * NOTE: Expects buff1 and buff2 to be of atleast `size` length (i.e. buff1[size-1] and buff2[size-1] are well defined)
*/
int getHammingDistance(int* buff1, int* buff2, int size){
    int hammingDistance = 0;
    for (int i=0; i<size; i++){
        // Hamming distance increases by one for every two elements that aren't the same
        hammingDistance += (*buff1)^(*buff2);

        buff1++; // Go to the next element
        buff2++; // Go to the next element
    }

    return hammingDistance;
}

/**
 * Writes a given integer buffer to a file (space separated)
*/
void writeBufferToFile(int* buffer, int size, const char* fileName){
    FILE* fd = fopen(fileName, "w");

    if (fd == NULL){
        printf("ERROR! Could not open file to write\n");
        abort();
    }

    for (int i = 0; i < size - 1; i++){
        fprintf(fd, "%d ", buffer[i]);
    }
    fprintf(fd, "%d", buffer[size - 1]); // Just to ensure we dont put a space after printing the last element
    
    fclose(fd);
}

/**
 * Populates the given buffer of given size with random bits using the given seed
*/
void storeRandomBitsIntoBuffer(int* buffer, int size, int seed){
    srand(seed);
    for (int i = 0; i < size; i++){
        buffer[i] = rand() % 2;
    }
}

/**
 * Parses the command line arguments to read the options
*/
struct Options getOptions(int argc, char** argv){
    if (argc < 5){
        printf("ERROR! Expected atleast 4 arguments (N, N Seed, M, M Seed)\n");
        abort();
    }
    struct Options options;
    options.n = atoi(argv[1]);
    options.nSeed = atoi(argv[2]);
    options.m = atoi(argv[3]);
    options.mSeed = atoi(argv[4]);

    // Remember atoi returns 0 if the argument was not an actual integer
    if (options.n <= 0){
        printf("ERROR! Expected N (argv[1]) to be a natural number\n");
        abort();
    }
    if (options.m <= 0){
        printf("ERROR! Expected M (argv[3]) to be a natural number\n");
        abort();
    }
    if (options.m >= options.n){
        printf("ERROR! Expected M (argv[3]) < N (argv[1])\n");
        abort();
    }

    return options;
}