/**
 * EE23B135 Kaushik G Iyer
 * 14/11/2023
 * 
 * Plots a nice graph `ee23b135_hemplot.jpg` using gnuplot on the points produced by the algorithm given in the question
 * Input:
 *  No inputs expected
 * 
 * Outputs:
 *  No output to stdout
 *  Spits out `ee23b135_hemplot.jpg`
 * 
 * NOTE: The plot seems very chaotic when viewed in a scale that shows all points uptill a certain iteration.
 * When more iterations are done the seemingly chaotic points form a closed loop (although by doing this new 'chaotic' points are introduced)
 * 
 * That is why I set the scale of the axis in such a way so as to 'hide' the chaotic points
*/ 

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define ITER_COUNT 50
#define THETA_STEP 0.005f

void hemplot(double x0, double y0, double r0);

int main(int argc, char** argv){
    hemplot(0, 0, 1);
}

/**
 * Implements the algorithm specified in the question
 * NOTE: The plot obtained seems very random but with more and more iterations we get to see some distinct 'loops' being formed
*/
void hemplot(double x0, double y0, double r0){
    long long int a = 0;
    long long int b = 0;
    long long int c = 1;

    double r = c * r0;
    double x = x0, y = y0, theta, thetaEnd;

    FILE* pipe = popen("gnuplot -persistent", "w");
    if (pipe == NULL){
        fprintf(stderr, "ERROR! Could not open pipe to gnuplot\n");
        abort();
    }

    fprintf(pipe, "set term jpeg size 1440, 1440 font \"Arial,20\";\n");
    fprintf(pipe, "set output \"ee23b135_hemplot.jpg\";\n");

    fprintf(pipe, "set title \"Hemplot (50 Iterations, thetaStep=0.005)\";\n");   
    fprintf(pipe, "plot '-' w p ls 7 ps 0.2 pt 7 title \"\";\n");   
    
    // Basically for the portion below I used simple logic that I formed after viewing the graph that needs to be formed
    // The graph in https://www.quora.com/How-Fibonacci-spiral-can-be-drawn-with-a-C-algorithm#:~:text=cout%3C%3C%22enter%20the%20number,%2F%2Fget%20the%20initial%20center
    // Was very helpful

    // We basically have to move are centres as required and then basically plot intervals of PI/2
    for (long long int i = 0; i < ITER_COUNT; i++){

        if (i%4 == 0){
            x0 = x0 + a;
            thetaEnd = 3 * M_PI_2;
        }
        else if (i%4 == 1){
            y0 = y0 + a;
            thetaEnd = 4 * M_PI_2;
        }
        else if (i%4 == 2){
            x0 = x0 - a;
            thetaEnd = 1* M_PI_2;
        }
        else if (i%4 == 3){
            y0 = y0 - a;
            thetaEnd = 2 * M_PI_2;
        }
        
        r = r0 * c;
        for (theta = thetaEnd - M_PI_2; theta < thetaEnd; theta += THETA_STEP){
            x = x0 + (r * cos(theta));
            y = y0 + (r * sin(theta));  
            fprintf(pipe, "%f %f\n", x, y);
        }

        a = b;
        b = c;
        c = a + b;
    }
    fprintf(pipe, "e;\n");   

    pclose(pipe);
}