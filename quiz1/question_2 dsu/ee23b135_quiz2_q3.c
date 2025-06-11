/**
 * EE23B135 Kaushik G Iyer
 * 22/09/2023
 * Quiz 1 (Question 2, Part 3)
 *
 * Solves the question described in the assignment using a decent method:
 *  - We store the `group` each planet is in (initially all are in separate groups)
 *  - For each query we check if the planets are in the same group, in that case they are connected (print 0)
 *  - If not we merge the two groups together, therefore we create a new connection (print 1)
 *
 * Expected input... (1 Argument)
 *  - input_file (A text file that conforms to the format described in the problem statement)
 *
 * Outputs a file called ee23b135_quiz2_q3_output.txt that contains...
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

#define OUT_FILE_NAME "ee23b135_quiz2_q3_output.txt"

struct Node{
    int weight;
    struct Node* parent;
};

struct Node* getRoot(struct Node* node);

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

    int n, m;
    fscanf(in_fd, "%d %d", &n, &m);

    struct Node* planets = malloc(n * sizeof(struct Node));

    if (planets == NULL){
        printf("Not enough space to allocated memory!\n");
        return 1;
    }

    for (int i = 0; i < n; i++){
        // Initially all nodes are separated (They are their own tree and therefore their weight is 1)
        planets[i].weight = 1;
        planets[i].parent = &planets[i];
    }

    int planet1, planet2; // These are set in the loop below
    struct Node *root1, *root2;

    for (int _=0; _<m; _++){ // Go through all the queries
        fscanf(in_fd, "%d %d", &planet1, &planet2);

        root1 = getRoot(&planets[planet1]);
        root2 = getRoot(&planets[planet2]);

        if (root1 == root2){
            // The planets are of the same tree, there is no need to build a new road here
            fprintf(out_fd, "0");
        }
        else{
            // Now we basically build a road between these two planets, thus connecting the two disconnected groups
            // We connect the smaller tree as a child of the bigger tree
            if (root1->weight > root2->weight){
                root2->parent = root1;
                root1->weight += root2->weight;
            }
            else{
                root1->parent = root2;
                root2->weight += root1->weight;
            }
            fprintf(out_fd, "1");
        }
    }

    free(planets);
    fprintf(out_fd, "\n%d ms", ((clock() - start) * 1000) / CLOCKS_PER_SEC); // Since CLOCKS_PER_SEC is an integer (in windows and posix systems atleast) I don't really need to cast them to doubles before dividing

    fclose(in_fd);
    fclose(out_fd);
}

struct Node* getRoot(struct Node* node){
    // NOTE: We are never going to hit any recursion limits since the maximum depth of our tree is atmost log(n)
    // Since n <= 10^8 (as provided in the addendum) the maximum depth is approximately 27

    if (node->parent == node){return node;}

    struct Node* parent = getRoot(node->parent);
    node->parent = parent;

    return parent;
}