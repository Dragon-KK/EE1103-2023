/**
 * EE23B135 Kaushik G Iyer
 * 22/09/2023
 * Quiz 1 (Question 2, Part 1)
 * 
 * Solves the question described in the assignment using a naive method:
 *  - We store the `group` each planet is in (initially all are in separate groups)
 *  - For each query we check if the planets are in the same group, in that case they are connected (print 0)
 *  - If not we merge the two groups together, therefore we create a new connection (print 1)
 * 
 * Expected input... (1 Argument)
 *  - input_file (A text file that conforms to the format described in the problem statement)
 * 
 * Outputs a file called ee23b135_quiz2_q1_output.txt that contains...
 * {(1 if road_must_be_built else 0) for connection in connections}\n
 * {time_taken_since_start_of_program} ms
 * --------------------------------------------
 * Format of input file:
 *  first line: N(int) M(int)
 *  next M lines: planet1(int) planet2(int)
 * 
 * Here N is the number of `planets` in our system
 * M is the number of queries
 * Each query asks wheter a connection must be made between planet1 and planet2 to make them connected
 * We must output `1` if a connection must be made or `0` if the planets are already connected
 * 
 * Format of output file:
 *  The first line of output contains an M length long sequence of 1s and 0s
 *  The second line contains the time taken by the progam to run
 * 
 * NOTE: I am using the clock() function due to the reasoning provided here (https://www.gnu.org/software/libc/manual/html_node/Processor-And-CPU-Time.html)
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define OUT_FILE_NAME "ee23b135_quiz2_q1_output.txt"

int main(int argc, char** argv){
    if (argc < 2){
        printf("Please provide the input file name as an argument :)\n");
        return 1;
    }

    FILE* in_fd = fopen(argv[1], "r");
    if (in_fd == NULL){
        printf("Could not open the given file!\n");
        return 1;
    }

    FILE* out_fd = fopen(OUT_FILE_NAME, "w");
    if (out_fd == NULL){
        printf("Could not open file to write!\n");
        return 1;
    }

    clock_t start = clock();

    int n, m; // n contains the number of planets, and m contains the number of queries asked
    fscanf(in_fd, "%d %d", &n, &m);

    /** Stores the `group` of the i'th planet */
    int* group = malloc(sizeof(int) * n);
    if (group == NULL){
        printf("Not enough space to allocated memory!\n");
        return 1;
    }

    for (int i = 0; i < n; i++){
        group[i] = i; // Initially each node is in its own group
    }

    int planet1, planet2, group1, group2; // These are set in the loop below

    for (int _=0; _<m; _++){ // Go through all the queries        
        fscanf(in_fd, "%d %d", &planet1, &planet2);

        // Get which group each planet belongs to
        group1 = group[planet1];
        group2 = group[planet2];
        if (group1 == group2){
            // The planets are of the same group, there is no need to build a new road here
            fprintf(out_fd, "0");
        }
        else{            
            // Now we basically build a road between these two planets, thus connecting the two disconnected groups
            // I arbitrarily decide that the group of the second planet of the query (planet2) shall be merged into the group of the first planet (planet1)
            // i.e. All the planets in group2 are now converted to group1
            for (int i = 0; i < n; i++){
                if (group[i] == group2){
                    group[i] = group1;
                }
            }
            fprintf(out_fd, "1");
        }
    }

    free(group);
    fprintf(out_fd, "\n%d ms", ((clock() - start) * 1000) / CLOCKS_PER_SEC); // Since CLOCKS_PER_SEC is an integer (in windows and posix systems atleast) I don't really need to cast them to doubles before dividing
    
    fclose(in_fd);
    fclose(out_fd);
}