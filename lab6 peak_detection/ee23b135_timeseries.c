/**
 * EE23B135 Kaushik G Iyer
 * 04/10/2023
 * 
 * An exercise on faking experimental data and estimating values based on the faked data lmao
 * 
 * Input:
 *  ./prog {M} {T} {a} {graphType?}
 * 
 * Where:
 *  - M is the number of peaks
 *  - T is the average distance between the peaks
 *  - a is half of the width at half height
 *  - graphType is either "Gaussian" | "Lorentzian"
 * 
 * Outputs:
 *  {estimated_mean_T} {estimated_mean_a} {stdev_T} {stdev_a}
 * 
 * NOTE: Setting too low values of `a` lead to there not being enough points to properly define a curve, please refrain from doing so :(
 * NOTE: Setting a small T/a ratio leads to the curves overlapping thus causing errors, pls don't do this boss :)
 *  
*/ 

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdbool.h>

// Set to truthy value to save the graph into `graph.dat`
#define SAVE_GRAPH 0

// Set to truthy value to filter the graph and also save the filtered graph into `clean_graph.dat`
#define FILTER_GRAPH 0

#define GAUSSIAN_FILTER_HALF_WIDTH 3.0f
#define LORENTZIAN_FILTER_HALF_WIDTH 3.0f

// Defines the clamp size for the gaussian random numbers produced
#define GAUSSIAN_CLAMP 6

// The x scale basically
#define SCALE 0.1f

// The relative error in A
#define REL_DEL_A 0.05f
// The relative error in the position of the peaks
#define REL_DEL_T 0.05f
// The relative error in the amplitude of the graph
#define REL_DEL_AMP 0.008f
// Defines what pretty much zero means (At what height can we cut off our peaks. This just leads to a cleaner graph in general)
#define PEAK_CUTOFF_THRESHOLD (REL_DEL_AMP / 2)

// Defines what error is acceptable error when searching for points near our threshold (Used in the estimation of peaks)
#define DISTANCE_THRESHOLD 0.05f
// This defines the minimum distance the mean of the in between points should be (Used to detec wheter two points of interest are on different sides of a peak)
#define MEAN_ABOVE_DISTANCE_THRESHOLD 0.2f


enum GraphType{
    LORENTZIAN, GAUSSIAN
};

struct Options{
    int M; // Number of peaks
    double T; // Distance between each lorentzian
    double a; // Half width of a lorentzian
    enum GraphType graphType; // "Lorentzian" | "Gaussian"
};

struct Peak{ // Stores information about the curves present in the graph
    double location; // mT
    double half_width; // a
};

struct Stats{ // Stores mean and standard deviation
    double mean;
    double standardDeviation;
};


void filter(double* graph, int graphSize, enum GraphType graphType);
struct Peak* estimate(const double* graph, int graphSize, int* detectedPeakCount, int);

void populatePeaks(struct Peak*, struct Options);
double getCurveValue(double x, struct Peak peak, enum GraphType graphType);
void populateGraph(double* graph, const struct Peak*, int, int, enum GraphType);

double randn();
struct Stats generateStatistics(const double* array, int size);
struct Stats generateLocationStatisticsForPeaks(const struct Peak* peaks, int size);
struct Stats generateHalfWidthStatisticsForPeaks(const struct Peak* peaks, int size);

struct Options getOptions(int argc, char** argv);
void save(double* graph, int graphSize, char* filename);

int main(int argc, char** argv){
    struct Options options = getOptions(argc, argv);

    struct Peak* peaks = malloc(sizeof(struct Peak) * options.M);
    if (peaks == NULL){
        printf("ERROR! Could not allocate memory for storing peaks\n");
        abort();
    }
    populatePeaks(peaks, options);

    // M*T is the location of the last peak, and since we expect T to be bigger than the curve width anyways, (M+1)*T should be enough widht to contain all the points
    int graphSize = (options.M + 1) * options.T / SCALE;    
    
    double* graph = malloc(sizeof(double) * graphSize);   
    if (graph == NULL){
        printf("ERROR! Could not allocate memory for storing graph\n");
        abort();
    } 
    populateGraph(graph, peaks, options.M, graphSize, options.graphType);    

    if (SAVE_GRAPH){ save(graph, graphSize, "graph.dat"); }

    if (FILTER_GRAPH){
        filter(graph, graphSize, options.graphType);
        save(graph, graphSize, "clean_graph.dat");
    }

    int detectedPeakCount;
    struct Peak* detectedPeaks = estimate(graph, graphSize, &detectedPeakCount, ((options.T/4) / SCALE));
    // To estimate the peaks, I am first telling my function what all x values to consider to find the `ambient noise`
    // Now, isn't it cheating to tell how far to look in terms of T (the parameter that is supposed to be found?)
    // Well, it doesn't really have to be exact T/4, any rough approximation of T/4 will hold
    // (It is expected that the rough approximation can be found just by viewing the graph)

    if (detectedPeakCount != options.M){
        printf("Warning! detected %d peaks while actual number of peaks is %d\nConsider playing with the T/a ratio\n", detectedPeakCount, options.M);
    }

    struct Stats aStats = generateHalfWidthStatisticsForPeaks(detectedPeaks, detectedPeakCount);
    struct Stats tStats = generateLocationStatisticsForPeaks(detectedPeaks, detectedPeakCount);

    printf("%f %f %f %f\n", tStats.mean, aStats.mean, tStats.standardDeviation, aStats.standardDeviation);

    free(graph);
    free(peaks);
    free(detectedPeaks);
}

#pragma region Estimation

/**
 * Filters the graph with a suitably shaped window filter
 * NOTE: The amplitude of the graphs decreases, this might cause problems for some specific input
 * ALso the estimated values with or without filtering come out to be the same
*/
void filter(double* graph, int graphSize, enum GraphType graphType){
    struct Peak peak;
    double* window;
    double* filtered;

    if (graphType == GAUSSIAN){
        peak.half_width = GAUSSIAN_FILTER_HALF_WIDTH; // Some arbitrary value (Too high values lead to the amplitude decreasing like craaaazy, If you want to do this, decrease MEAN_ABOVE_DISTANCE_THRESHOLD to a suitable value)
        peak.location = peak.half_width; // I want the peak to be in the middle of the window      
    }
    else if (graphType == LORENTZIAN) {
        peak.half_width = LORENTZIAN_FILTER_HALF_WIDTH; // Just some arbitrary value
        peak.location = peak.half_width;
    }
    else{
        printf("ERROR! Unsupported graph type `%d`\n", graphType);
        abort();
    }
    int halfWindowWidth =  peak.half_width / SCALE; // I define the window size based on my half width
    window = malloc(sizeof(double) * halfWindowWidth * 2);
    for (int i = 0; i < halfWindowWidth * 2; i++){
        window[i] = getCurveValue(i * SCALE, peak, graphType);
    }
    filtered = malloc(sizeof(double) * graphSize);
    for (int xi = 0; xi < graphSize; xi++){
        double sum = 0;
        for (int i = -halfWindowWidth; i < halfWindowWidth; i++){
            if (xi + i < 0 || xi + i > graphSize){ continue; }
            sum += window[i + halfWindowWidth] * graph[xi + i];
        }
        filtered[xi] = sum/(2*halfWindowWidth);
    }
    for (int xi = 0; xi < graphSize; xi++){
        graph[xi] = filtered[xi];
    }
    free(filtered);
    free(window);  
}

/**
 * Estimates the peaks present in the graph
 * 
 * We do this by taking advantage of the fact that our peaks are symmetric
 * Therfore we simply find two points on different sides of the peak that have the same y value
 * The location of the peak is then given by the average of the x value of the two points taken in the step before :)
 * 
 * NOTE: The returned array must be freed!
*/
struct Peak* estimate(const double* graph, int graphSize, int* detectedPeakCount, int noiseSampleSize){
    
    // First let us a find a suitable y value to check for intersections
    // Since our graph still might contain some noise, ideally we would want a value that is well above the noisy 0 level
    struct Stats stats = generateStatistics(graph, noiseSampleSize);
    
    // Therefore a good level to check would be the (mean+6*sigma)
    // (mean and sigma are of the initial 0 level part in the graph)
    double y = stats.mean + (6*stats.standardDeviation);

    // NOTE: I call it prefixSum, but its really 'prefixSum of points that are of interest (close to the y level)'
    double* prefixSum = malloc(sizeof(double) * graphSize);
    if (prefixSum == NULL){
        printf("ERROR! Could not allocate memory to estimate the peaks\n");
        abort();
    }
    // Stores the prefix sums (to make calculating the mean fast)
    // This takes advantage of the fact that the maximum sum of all points in my array is atmost:
    // GLOBAL_MAX * NUMBER_OF_POINTS
    // 1 * (M + 1) * T / SCALE
    // = (M + 1) * T * 10
    // Which is expected to be less that 10^9
    // TL: Our prefix sum never overflows, so we lose nothing by doing this

    double currSum = 0;

    for (int xi = 0; xi < graphSize; xi++){
        currSum += graph[xi];
        prefixSum[xi] = currSum;

        if (fabs(graph[xi] - y) > DISTANCE_THRESHOLD){
            prefixSum[xi] = 0; // Basically only the prefix sum of points close to the y level are stored
            continue;
        }
    }
    // Right now, all the non zero indices in the prefixSum array correspond to points of interest
    // Now we go through the points of interest and only choose 2 points per peak (one on the left and one on the right)

    int firstPointOfInterest = -1;
    int lastPointOfInterest = -1;
    // Just find the first point of interest
    for (int xi = 0; xi < graphSize; xi++){
        if (prefixSum[xi] == 0) {continue;} // It is not a point of interest
        lastPointOfInterest = xi;
        firstPointOfInterest = xi;
        break;
    }

    if (firstPointOfInterest == -1){
        printf("ERROR! Could not find any point of interest!\n");
        abort();
    }

    // We are forced to iterate through our graph to find the number of peaks (in order to store the peaks in an array we need to know the sizes)
    *detectedPeakCount = 0;
    for (int xi = lastPointOfInterest + 1; xi < graphSize; xi++){
        if (prefixSum[xi] == 0) {continue;} // It is not a point of interest

        // The mean of two points of interest (will be relatively higher than the y level (Since on average all points in between are greater than y))
        if ((prefixSum[xi] - prefixSum[lastPointOfInterest]) / (xi - lastPointOfInterest) > MEAN_ABOVE_DISTANCE_THRESHOLD + y){
            // We detect a peak!
            (*detectedPeakCount)++;
        }
        // If this is not true, then both points of interest are on the same side of the peak
        lastPointOfInterest = xi;
    }

    struct Peak* peaks = malloc(sizeof(struct Peak) * (*detectedPeakCount));
    if (peaks == NULL){
        printf("ERROR! Could not allocate enough memory to store estimated peaks\n");
        abort();
    }

    int currPeak = 0;
    lastPointOfInterest = firstPointOfInterest;
    for (int xi = lastPointOfInterest + 1; xi < graphSize; xi++){
        if (prefixSum[xi] == 0) {continue;} // It is not a point of interest

        // The mean of two points of interest (will be relatively higher than the y level (Since on average all points in between are greater than y))
        if ((prefixSum[xi] - prefixSum[lastPointOfInterest]) / (xi - lastPointOfInterest) > MEAN_ABOVE_DISTANCE_THRESHOLD + y){
            // We detect a peak!
            peaks[currPeak].location = ((double)(xi + lastPointOfInterest)*SCALE)/2;
            
            double halfAmplitude = graph[(xi + lastPointOfInterest) / 2] / 2;
            // Now we must find the x values that are close to our halfAmplitudeValue
            // We expect the xp corresponding to the halfAmplitudeValue to be between our lastPointOfInterst and xi
            
            int firstValidXp = -1; // The first xp that is close to halfAmplitude (It will be on the left of the peak)
            int lastValidXp = -1; // The last xp that is close to halfAmplitude (It will be on the right of the peak)
            for (int xp = lastPointOfInterest+1; xp < xi; xp++){
                if (fabs(graph[xp] - halfAmplitude) > DISTANCE_THRESHOLD){ continue; }
                firstValidXp = xp;
                break;
            }
            if (firstValidXp == -1){
                printf("ERROR! Could not find any points corresponding to half amplitude of peak number %d\nConsider increasing the `DISTANCE_THRESHOLD`, using a value of `a` that is bigger than 1 or playing with the T/a ratio\n", currPeak);
                abort();
            }

            for (int xp = firstValidXp + 1; xp < xi; xp++){
                if (fabs(graph[xp] - halfAmplitude) > DISTANCE_THRESHOLD){ continue; }

                lastValidXp = xp;
            }

            if (lastValidXp == -1){
                printf("ERROR! Could not find two points corresponding to half amplitude of peak number %d\nConsider increasing the `DISTANCE_THRESHOLD`, using a value of `a` that is bigger than 1 or playing with the T/a ratio\n", currPeak);
                abort();
            }

            // The width is then given by (lastValidXp - firstValidXp)*SCALE
            peaks[currPeak].half_width =  ((double)(lastValidXp - firstValidXp)*SCALE)/2;

            currPeak++;
        }
        // If this is not true, then both points of interest are on the same side of the peak
        lastPointOfInterest = xi;
    }

    return peaks;
}
#pragma endregion


#pragma region Data Faking

/**
 * Gets the y value of the function (function is defined based on the values stored in the Peak struct)
*/
double getCurveValue(double x, struct Peak peak, enum GraphType graphType){
    if (graphType == LORENTZIAN){
        double a2 = (peak.half_width * peak.half_width);
        double locationFactor = (x - peak.location);
        return a2 / ((locationFactor*locationFactor) + a2);
    }
    else if (graphType == GAUSSIAN) {
        double locationFactor = (x - peak.location);
        double power = locationFactor*locationFactor*log(2)/(peak.half_width * peak.half_width);
        return powl(M_E, -power);
    }
    else{
        printf("ERROR! Unsupported graph type `%d`\n", graphType);
        abort();
    }
}


/**
 * Basically creates the graph (Using the peaks)
 * NOTE: I don't just superimpose all the curves, inorder to better simulate the experiment we are modelling after, I individually add the curves
 * (At any point only the y value of a single curve is considered)
 * NOTE: This adds the noise to the amplitude too
*/
void populateGraph(double* graph, const struct Peak* peaks, int peakCount, int graphSize, enum GraphType graphType){    
    int currPeak = 0;
    for (int xi = 0; xi < graphSize; xi++){
        // graph[xi] is the ith point in our graph it is NOT the graph at x=xi
        // graph[xi] is the graph at x=xi * SCALE

        double currPeakValue = (currPeak<peakCount)? getCurveValue(xi*SCALE, peaks[currPeak], graphType) : 0;
        double nextPeakValue = (currPeak+1<peakCount)? getCurveValue(xi*SCALE, peaks[currPeak+1], graphType) : 0;

        if (currPeakValue < PEAK_CUTOFF_THRESHOLD){
            currPeakValue = 0; // This just makes the function look a lot cleaner
            // The idea being if the value of the function itself is very small (so small that the noise added only is bigger),
            // There is really no difference in not adding that small value (Except for the fact that the graph looks sort of like the experiment we are modellings)
        }
        if (nextPeakValue > currPeakValue + PEAK_CUTOFF_THRESHOLD){
            // The next peak is now more important than the current peak
            currPeak = (currPeak>=peakCount)? peakCount : (currPeak+1);
            currPeakValue = nextPeakValue;
        }
        
        graph[xi] = currPeakValue;

        graph[xi] += randn() * REL_DEL_AMP;
    }
}


/**
 * Sets the location and width of the peaks
 * NOTE: This also adds the errors of widths and location
*/
void populatePeaks(struct Peak* peaks, struct Options options){
    for (int m = 1; m <= options.M; m++){
        peaks[m - 1].half_width = (options.a) + (randn() * REL_DEL_A * options.a);
        peaks[m - 1].location = (m * options.T) + (randn() * REL_DEL_T * options.T);
    }
}
#pragma endregion


#pragma region Maths

/**
 * Outputs random numbers that are part of a gaussian distribution (Uses box muller transformation)
 * NOTE: The output is clamped between -GAUSSIAN_CLAMP and +GAUSSIAN_CLAMP
*/
double randn(){
    // https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform

    double u1 = (double)rand()/(double)RAND_MAX;
    double u2 = (double)rand()/(double)RAND_MAX;;

    double res = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);

    if (fabs(res) < GAUSSIAN_CLAMP){
        return res;
    }

    return (res>0)? GAUSSIAN_CLAMP : -GAUSSIAN_CLAMP;
}

/**
 * Calculates the standard deviation and mean of the data
*/
struct Stats generateStatistics(const double* array, int size){
    double sum = 0;
    double sum2 = 0; // Sum of the squares
    for (int i = 0; i < size; i++){
        sum+=array[i];
        sum2+=array[i]*array[i];
    }
    double mean = sum/size;

    struct Stats stats;
    stats.mean = mean;
    stats.standardDeviation = sqrt((sum2/size) - (mean*mean));
    return stats;
}

struct Stats generateLocationStatisticsForPeaks(const struct Peak* peaks, int size){
    // You may ask, is there really a need for this to be a seperate function if I'm only using it once?
    // No, I just thought I'd use it more than I actually had tos
    double sum = 0;
    double sum2 = 0; // Sum of the squares
    double lastPeak = 0;
    for (int i = 0; i < size; i++){
        double T = peaks[i].location - lastPeak;
        sum+=T;
        sum2+=T*T;

        lastPeak = peaks[i].location;
    }
    double mean = sum/size;

    struct Stats stats;
    stats.mean = mean;
    stats.standardDeviation = sqrt((sum2/size) - (mean*mean));
    return stats;
}

struct Stats generateHalfWidthStatisticsForPeaks(const struct Peak* peaks, int size){
    double sum = 0;
    double sum2 = 0; // Sum of the squares
    for (int i = 0; i < size; i++){
        sum+=peaks[i].half_width;
        sum2+=peaks[i].half_width*peaks[i].half_width;
    }
    double mean = sum/size;

    struct Stats stats;
    stats.mean = mean;
    stats.standardDeviation = sqrt((sum2/size) - (mean*mean));
    return stats;
}
#pragma endregion


/**
 * Parses the command line arguments to read the options
*/
struct Options getOptions(int argc, char** argv){
    if (argc < 4){
        printf("ERROR! Expected atleast 3 arguments (M, T, a)\n");
        abort();
    }
    struct Options options;
    options.M = atoi(argv[1]);
    options.T = atof(argv[2]);
    options.a = atof(argv[3]);
    options.graphType = LORENTZIAN;

    if (argc > 4){
        if (strcmp(argv[4], "Lorentzian") == 0){
            options.graphType = LORENTZIAN;
        }
        else if (strcmp(argv[4], "Gaussian") == 0){
            options.graphType = GAUSSIAN;
        }
        else{
            printf("Warning! Invalid graph type `%s`, expected either `Lorentzian` or `Gaussian`\n", argv[4]);
        }
    }

    if (options.M <= 0){
        printf("ERROR! Expected M (argv[1]) to be a natural number\n");
        abort();
    }
    if (options.T <= 0){
        printf("ERROR! Expected T (argv[2]) to be positive\n");
        abort();
    }
    if (options.a <= 0){
        printf("ERROR! Expected a (argv[3]) to be positive\n");
        abort();
    }   

    return options;
}

/**
 * Saves the graph into a file
 * Each line of output file contains `{xi} {yi}`
*/
void save(double* graph, int graphSize, char* filename){
    FILE* fd = fopen(filename, "w");
    if (fd == NULL){
        printf("ERROR! Could not open file to write\n");
        abort();
    }
    
    for (int xi = 0; xi < graphSize; xi++){
        fprintf(fd, "%f %f\n", xi*SCALE, graph[xi]);
    }
    fclose(fd);
}
