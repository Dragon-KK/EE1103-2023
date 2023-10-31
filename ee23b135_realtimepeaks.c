/**
 * EE23B135 Kaushik G Iyer
 * 12/10/2023
 * 
 * An exercise on faking experimental data and estimating values (simulating a realtime environment) based on the faked data lmao
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
 * NOTE: When I run this code on my linux machine, it runs a lot slower. Specifically, saving the graph to the file takes forever.
 * This is probably due to the lack of buffering in file io (I don't know how to force it to buffer)
 * 
 * NOTE: When running with very less peaks, there is a chance the standard deviation may come out to be imaginary (Due to some error in precision)
*/ 

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdbool.h>

#pragma region Constants
// The name of the file the graph is saved into
#define INPUT_FILE_NAME "graph.dat"

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
// This defines the minimum distance the mean of the in between points should be (Used to detect wheter two points of interest are on different sides of a peak)
#define MEAN_ABOVE_DISTANCE_THRESHOLD 0.2f

// The size by which the lists grow
#define LIST_CHUNK_SIZE 1000
#pragma endregion


#pragma region Structs
enum GraphType{
    LORENTZIAN, GAUSSIAN
};

struct Options{
    int M; // Number of peaks
    double T; // Distance between each lorentzian
    double a; // Half width of a lorentzian
    enum GraphType graphType; // "Lorentzian" | "Gaussian"
    double smootheningConstant; // Used in the exponential filter
};

struct Peak{ // Stores information about the curves present in the graph
    double location; // mT
    double half_width; // a
    double amplitude;
};

struct Stats{ // Stores mean and standard deviation
    double mean;
    double standardDeviation;
};

struct Point{
    double x;
    double y;
};

struct PointList{
    struct Point* chunk;
    int size;
    struct PointList* next;
};

struct PeakList{
    struct Peak* chunk;
    int size;
    struct PeakList* next;
};

struct PointListIterator{
    struct PointList* workingList;
    int index;
};
struct PeakListIterator{
    struct PeakList* workingList;
    int index;
};
#pragma endregion


#pragma region Definitions
struct PeakList* estimate(char* input_file, int* detectedPeakCount, int noiseSampleSize, double smootheningConstant);

void populatePeaks(struct Peak*, struct Options);
double getCurveValue(double x, struct Peak peak, enum GraphType graphType);
void populateGraph(double* graph, const struct Peak*, int, int, enum GraphType);

double randn();
struct Stats generateStatistics(const double* array, int size);
struct Stats generateLocationStatisticsForPeaks(struct PeakList* peaks, int size);
struct Stats generateHalfWidthStatisticsForPeaks(struct PeakList* peaks, int size);

struct PointList* createPointList();
void appendPoint(double x, double y, struct PointList* workingList);
void clearPointList(struct PointList*);
void freePointList(struct PointList*);

struct PeakList* createPeakList();
void appendPeak(double location, double half_width, double amplitude, struct PeakList* workingList);
void freePeakList(struct PeakList*);

bool pointIteratorNotFinished(struct PointListIterator* iter);
struct Point* nextPoint(struct PointListIterator* iter);
bool peakIteratorNotFinished(struct PeakListIterator* iter);
struct Peak* nextPeak(struct PeakListIterator* iter);

struct Options getOptions(int argc, char** argv);
void save(double* graph, int graphSize, char* filename);
#pragma endregion


int main(int argc, char** argv){
    struct Options options = getOptions(argc, argv);

    struct Peak* peaks = malloc(sizeof(struct Peak) * options.M);
    if (peaks == NULL){
        fprintf(stderr, "ERROR! Could not allocate memory for storing peaks\n");
        abort();
    }
    populatePeaks(peaks, options);

    // M*T is the location of the last peak, and since we expect T to be bigger than the curve width anyways, (M+1)*T should be enough widht to contain all the points
    int graphSize = (options.M + 1) * options.T / SCALE;    
    
    double* graph = malloc(sizeof(double) * graphSize);   
    if (graph == NULL){
        fprintf(stderr, "ERROR! Could not allocate memory for storing graph\n");
        abort();
    } 
    populateGraph(graph, peaks, options.M, graphSize, options.graphType);    
    free(peaks);

    save(graph, graphSize, INPUT_FILE_NAME);
    free(graph);

    int detectedPeakCount;
    struct PeakList* detectedPeaks = estimate(INPUT_FILE_NAME, &detectedPeakCount, ((options.T/4) / SCALE), options.smootheningConstant);
    // To estimate the peaks, I am first telling my function what all x values to consider to find the `ambient noise`
    // Now, isn't it cheating to tell how far to look in terms of T (the parameter that is supposed to be found?)
    // Well, it doesn't really have to be exact T/4, any rough approximation of T/4 will hold
    // (It is expected that the rough approximation can be found just by viewing the graph)

    if (detectedPeakCount != options.M){
        fprintf(stderr, "Warning! detected %d peaks while actual number of peaks is %d\nConsider playing with the T/a ratio\n", detectedPeakCount, options.M);
    }

    struct Stats aStats = generateHalfWidthStatisticsForPeaks(detectedPeaks, detectedPeakCount);
    struct Stats tStats = generateLocationStatisticsForPeaks(detectedPeaks, detectedPeakCount);

    freePeakList(detectedPeaks);

    printf("%f %f %f %f\n", tStats.mean, aStats.mean, tStats.standardDeviation, aStats.standardDeviation);
}

#pragma region Estimation

/**
 * Estimates the peaks present in the graph
 * 
 * We do this by taking advantage of the fact that our peaks are symmetric
 * Therfore we simply find two points on different sides of the peak that have the same y value
 * The location of the peak is then given by the average of the x value of the two points taken in the step before :)
 * 
 * NOTE: The returned list must be freed!
*/
struct PeakList* estimate(char* input_file, int* detectedPeakCount, int noiseSampleSize, double smootheningConstant){
    FILE* fd = fopen(input_file, "r");
    if (fd == NULL){
        fprintf(stderr, "ERROR! Could not open input file\n");
        abort();
    }
    double x;
    double y;

    int amountOfNoiseSampled = 0;
    double *noiseSample = malloc(sizeof(double) * noiseSampleSize);
    if (noiseSample == NULL){
        fprintf(stderr, "ERROR! Could not allocate memory for storing noise sample\n");
        abort();
    }
    
    double y0 = 0; // The y level we check for intersestions (Preferably above all the noise)
    
    bool firstPointOfInterestIsSet = false; // Initially we don't know where our first point of interest is, this flag just helps ensure we don't get a false peak
    double lastPointOfInterest; // Stores the x value of the last point that was close to y0 level

    double pointSinceLastPointOfInterest = 0; // Used to calculate the mean along with cumSum
    double cumSum = 0; // The cumulative sum starting from the last point of interest (Used to check whether two points of interest are on different sides of a peak)

    *detectedPeakCount = 0;
    struct PointList* points = createPointList();
    struct PeakList* peaks = createPeakList();

    double lastY = 0;

    while (fscanf(fd, "%lf %lf ", &x, &y) == 2)
    {
        y = (lastY * smootheningConstant) + (y * (1 - smootheningConstant));
        lastY = y;

        // We assume that the initial input to us does not immediately contain a peak
        // It is a huge assumption yes, but it makes stuff a lot easier
        // An alternative (the most general solution) would be to keep recalculating the mean and std at different intervals,
        // Then we can safely assume the lowest mean + 6*std to correspond to the base level of noise
        if (amountOfNoiseSampled < noiseSampleSize){
            noiseSample[amountOfNoiseSampled] = y;
            amountOfNoiseSampled++;

            if (amountOfNoiseSampled == noiseSampleSize){
                struct Stats stats = generateStatistics(noiseSample, noiseSampleSize);
                y0 = stats.mean + (6*stats.standardDeviation); // mean + 6sigma will ensure that most of the noise (if any) is disregarded
                free(noiseSample);
            }
            continue;
        }

        cumSum += y;
        appendPoint(x, y, points);
        pointSinceLastPointOfInterest++;

        if (fabs(y - y0) > DISTANCE_THRESHOLD){
            continue; // This point is not close to our y0 level so it is not a point of interest
        }

        if (firstPointOfInterestIsSet && (cumSum / pointSinceLastPointOfInterest) > MEAN_ABOVE_DISTANCE_THRESHOLD + y0){
            // i.e There is a peak bw the last point of interest and x
            (*detectedPeakCount)++;

            // Now we need to check the points that we stored to find the amplitude of the peak
            double xpeak = (lastPointOfInterest + x) / 2; // Location of the amplitudes
            
            double amplitude = 0;
            double closestX = lastPointOfInterest; // Let it store the minimum x distance from xpeak
            
            // Iterate through to find point closest to xpeak
            struct PointListIterator iter = {points, 0};
            struct Point* point;
            while(pointIteratorNotFinished(&iter)){
                point = nextPoint(&iter);
                if (fabs(point->x - xpeak) > fabs(closestX - xpeak)){
                    // We are now going further away from the peak
                    break;
                }
                amplitude = point->y;
                closestX = point->x;                
            }

            // Reset the iterator
            iter.workingList = points;
            iter.index = 0;

            double half_amplitude = amplitude/2;

            double firstXp = -1;
            double lastXp = -1;
            // Iterate through the list again to find the x values of points close to the half amplitude
            while(pointIteratorNotFinished(&iter)){
                point = nextPoint(&iter);
                if (fabs(point->y - half_amplitude) > DISTANCE_THRESHOLD){
                    continue;
                }
                firstXp = point->x;       
                break;
            }
            if (firstXp == -1){
                fprintf(stderr, "ERROR! Could not find any points corresponding to half amplitude of peak number %d\nConsider increasing the `DISTANCE_THRESHOLD`, using a value of `a` that is bigger than 1 or playing with the T/a ratio\n", *detectedPeakCount);
                abort();
            }
            while(pointIteratorNotFinished(&iter)){
                point = nextPoint(&iter);
                if (fabs(point->y - half_amplitude) > DISTANCE_THRESHOLD){
                    continue;
                }
                lastXp = point->x;
            }
            if (lastXp == -1){
                fprintf(stderr, "ERROR! Could not find two points corresponding to half amplitude of peak number %d\nConsider increasing the `DISTANCE_THRESHOLD`, using a value of `a` that is bigger than 1 or playing with the T/a ratio\n", *detectedPeakCount);
                abort();
            }

            appendPeak(xpeak, (lastXp - firstXp)/2, amplitude, peaks);
        }

        firstPointOfInterestIsSet = true;
        lastPointOfInterest = x;
        cumSum = 0;
        pointSinceLastPointOfInterest = 0;    
        clearPointList(points);
    }

    if (amountOfNoiseSampled < noiseSampleSize){
        free(noiseSample); // In case our file ends before we stop reading the sample
    }
    freePointList(points);
    
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
        fprintf(stderr, "ERROR! Unsupported graph type `%d`\n", graphType);
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
        peaks[m - 1].amplitude = 1;
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

struct Stats generateLocationStatisticsForPeaks(struct PeakList* peaks, int size){
    // You may ask, is there really a need for this to be a seperate function if I'm only using it once?
    // No, I just thought I'd use it more than I actually had tos
    double sum = 0;
    double sum2 = 0; // Sum of the squares
    double lastPeak = 0;

    struct PeakListIterator iter = {peaks, 0};
    struct Peak* peak;
    while (peakIteratorNotFinished(&iter)){
        peak = nextPeak(&iter);
        double T = peak->location - lastPeak;
        sum+=T;
        sum2+=T*T;

        lastPeak = peak->location;
    }
    double mean = sum/size;

    struct Stats stats;
    stats.mean = mean;
    stats.standardDeviation = sqrt((sum2/size) - (mean*mean));
    return stats;
}

struct Stats generateHalfWidthStatisticsForPeaks(struct PeakList* peaks, int size){
    double sum = 0;
    double sum2 = 0; // Sum of the squares
    struct PeakListIterator iter = {peaks, 0};
    struct Peak* peak;
    while (peakIteratorNotFinished(&iter)){
        peak = nextPeak(&iter);
        sum+=peak->half_width;
        sum2+=peak->half_width*peak->half_width;
    }
    double mean = sum/size;

    struct Stats stats;
    stats.mean = mean;
    stats.standardDeviation = sqrt((sum2/size) - (mean*mean));
    return stats;
}
#pragma endregion


#pragma region List
/**
 * Creates a new point list
 * NOTE: Both the point list and the PointList.chunk must be freed
*/
struct PointList* createPointList(){
    struct PointList* pointList = malloc(sizeof(struct PointList));
    pointList->next = NULL;
    pointList->size = 0;
    pointList->chunk = malloc(sizeof(struct Point) * LIST_CHUNK_SIZE);

    if (pointList->chunk == NULL){
        fprintf(stderr, "ERROR! Could not allocate memory for point list\n");
        abort();
    }

    return pointList;
}

/**
 * Appends a point to the list
*/
void appendPoint(double x, double y, struct PointList* workingList){
    if (workingList->size == LIST_CHUNK_SIZE){
        if (workingList->next == NULL){
            workingList->next = createPointList();
        }
        return appendPoint(x, y, workingList->next);        
    }
    workingList->chunk[workingList->size].x = x;
    workingList->chunk[workingList->size].y = y;
    workingList->size++;
}

/**
 * Clears all elements in the PointList
 * NOTE: This does not free anything
*/
void clearPointList(struct PointList* head){
    head->size = 0;
    if (head->next != NULL){
        clearPointList(head->next);
    }
}

/**
 * Frees all memory allocated in the point list
*/
void freePointList(struct PointList* head){
    if (head->next != NULL){
        freePointList(head->next);
    }
    free(head->chunk);
    free(head);
}

/**
 * Creates a new peak list
 * NOTE: Both the peak list and the PeakList.chunk must be freed
*/
struct PeakList* createPeakList(){
    struct PeakList* peakList = malloc(sizeof(struct PeakList));
    peakList->next = NULL;
    peakList->size = 0;
    peakList->chunk = malloc(sizeof(struct Peak) * LIST_CHUNK_SIZE);

    if (peakList->chunk == NULL){
        fprintf(stderr, "ERROR! Could not allocate memory for peak list\n");
        abort();
    }

    return peakList;
}

/**
 * Clears all elements in the PeakList
 * NOTE: This does not free anything
*/
void appendPeak(double location, double half_width, double amplitude, struct PeakList* workingList){
    if (workingList->size == LIST_CHUNK_SIZE){
        if (workingList->next == NULL){
            workingList->next = createPeakList();
        }
        return appendPeak(location, half_width, amplitude, workingList->next);        
    }
    workingList->chunk[workingList->size].half_width = half_width;
    workingList->chunk[workingList->size].location = location;
    workingList->chunk[workingList->size].amplitude = amplitude;
    workingList->size++;
}

/**
 * Frees all memory allocated in the peak list
*/
void freePeakList(struct PeakList* head){
    if (head->next != NULL){
        freePeakList(head->next);
    }
    free(head->chunk);
    free(head);
}

bool pointIteratorNotFinished(struct PointListIterator* iter){
    return (iter->workingList != NULL) && (iter->index < iter->workingList->size);
}
struct Point* nextPoint(struct PointListIterator* iter){
    struct Point* ret =  &iter->workingList->chunk[iter->index];
    iter->index++;
    if (iter->index >= iter->workingList->size){
        iter->workingList = iter->workingList->next;
        iter->index = 0;
    }
    return ret;    
}

bool peakIteratorNotFinished(struct PeakListIterator* iter){
    return (iter->workingList != NULL) && (iter->index < iter->workingList->size);
}
struct Peak* nextPeak(struct PeakListIterator* iter){
    struct Peak* ret =  &iter->workingList->chunk[iter->index];
    iter->index++;
    if (iter->index >= iter->workingList->size){
        iter->workingList = iter->workingList->next;
        iter->index = 0;
    }
    return ret;    
}

#pragma endregion


/**
 * Parses the command line arguments to read the options
*/
struct Options getOptions(int argc, char** argv){
    if (argc < 4){
        fprintf(stderr, "ERROR! Expected atleast 3 arguments (M, T, a)\n");
        abort();
    }
    struct Options options;
    options.M = atoi(argv[1]);
    options.T = atof(argv[2]);
    options.a = atof(argv[3]);
    options.graphType = LORENTZIAN;
    options.smootheningConstant = 0.01f; // Just a nice default value that seemed to work well (Just by observing the mean values for different smootheningConstants)

    if (argc > 4){
        if (strcmp(argv[4], "Lorentzian") == 0){
            options.graphType = LORENTZIAN;
        }
        else if (strcmp(argv[4], "Gaussian") == 0){
            options.graphType = GAUSSIAN;
        }
        else{
            fprintf(stderr, "Warning! Invalid graph type `%s`, expected either `Lorentzian` or `Gaussian`\n", argv[4]);
        }
    }

    if (argc > 5){ options.smootheningConstant = atof(argv[5]); }

    if (options.M <= 0){
        fprintf(stderr, "ERROR! Expected M (argv[1]) to be a natural number\n");
        abort();
    }
    if (options.T <= 0){
        fprintf(stderr, "ERROR! Expected T (argv[2]) to be positive\n");
        abort();
    }
    if (options.a <= 0){
        fprintf(stderr, "ERROR! Expected a (argv[3]) to be positive\n");
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
        fprintf(stderr, "ERROR! Could not open file to write\n");
        abort();
    }
    
    for (int xi = 0; xi < graphSize; xi++){
        fprintf(fd, "%f %f\n", xi*SCALE, graph[xi]);
    }
    fclose(fd);
}
