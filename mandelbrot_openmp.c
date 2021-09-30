#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include <unistd.h>

#include "matrix.h"
#include "util.h"

#include <omp.h>

// Compile take 4 arguements
// iterations
// width
// height
// Set Num (0 = Julia, 1 = Mandlebrot, 2 = Mltibrot)
// Number of threads
// gcc -Wall -O3 -fopenmp -march=native mandelbrot_openmp.c matrix.c util.c -o mandel -lm
// srun --cpus-per-task=8 ./mandel 1000 1000 1000 1 32
// Plot after you compile and run
// python3 plot.py

void julia_set(Matrix* pixels, int max_iteration, int num_threads) {

    int R = 20;  // choose R > 0 such that R**2 - R >= sqrt(cx**2 + cy**2)
    float cx = -0.7269;
    float cy = 0.1889;
    size_t width = pixels->cols;
    size_t height = pixels->rows;

    #pragma omp parallel num_threads(num_threads)
    for (int x = 0; x < width; x++) {
        #pragma omp for schedule(dynamic, 32)
        for (int y = 0; y < height; y++) {

            float scaled_x = (1.f - (-2.5f)) / width * x - 2.5f;
            float scaled_y = (1.f - (-1.f)) / height * y - 1.f;
            int iteration = 0;
            while (scaled_x*scaled_x + scaled_y*scaled_y <= R*R && iteration < max_iteration) {
                float xtemp = scaled_x*scaled_x - scaled_y*scaled_y;
                scaled_y = 2*scaled_x*scaled_y + cy;
                scaled_x = xtemp + cx;
                iteration += 1;
            }
            pixels->data[y*pixels->rows + x] = iteration;
        }
        
    }

}

void mandelbrot(Matrix* pixels, int max_iteration, int num_threads) {

    size_t width = pixels->cols;
    size_t height = pixels->rows;

    #pragma omp parallel num_threads(num_threads)
    for (int x = 0; x < width; x++) {
        #pragma omp for schedule(dynamic, 32)
        for (int y = 0; y < height; y++) {

            float scaled_x = (1.f - (-2.5f)) / width * x - 2.5f;
            float scaled_y = (1.f - (-1.f)) / height * y - 1.f;
            float real = 0.f;
            float imaginary = 0.f;
            int iteration = 0;
            while (real*real + imaginary*imaginary <= 2*2 && iteration <= max_iteration) {
                float xtemp = real*real - imaginary*imaginary + scaled_x;
                imaginary = 2*real*imaginary + scaled_y;
                real = xtemp;
                iteration += 1;
            }
            pixels->data[y*pixels->rows + x] = iteration;
        }
        
    }
}

void multibrot(Matrix* pixels, int max_iteration) {

    size_t width = pixels->cols;
    size_t height = pixels->rows;

    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {

            // printf("%d, %d\n", x, y);
            float scaled_x = (1.f - (-2.5f)) / width * x - 2.5f;
            float scaled_y = (1.f - (-1.f)) / height * y - 1.f;
        
            int iteration = 0;
            while (scaled_x*scaled_x + scaled_y*scaled_y <= (2*2) && iteration < max_iteration) {

                int a = scaled_x;
                int b = scaled_y;

                float xtmp= (scaled_x * scaled_x * scaled_x * scaled_x * scaled_x) -
                            10 * (scaled_x * scaled_x * scaled_x) * (scaled_y * scaled_y) + 5 * scaled_x * 
                            (scaled_y * scaled_y * scaled_y * scaled_y) + a;

                scaled_y = 5 * (scaled_x * scaled_x * scaled_x * scaled_x)*scaled_y-
                10 * (scaled_x * scaled_x) * (scaled_y * scaled_y * scaled_y) + (scaled_y * scaled_y * scaled_y * scaled_y * scaled_y) + b;
                scaled_x = xtmp;

                iteration = iteration + 1;
            }

            pixels->data[y*pixels->rows + x] = iteration;

            // printf("%d\n", iteration);
            // printf("%d, %d\n", x, y);
        
            // if (iteration = max_iteration)
            //     colour = black;
            // else
            //     colour = iteration;
        
            // plot(scaled_x, scaled_y, colour)
        }
    }

}

int main(int argn, const char* argv[])
{
    // default command-line options
    int iterations = 100;   
    int width = 100; 
    int height = 100; 
    // Specify set to plot (0 Julia, 1 Mandelbrot, 2 Multibrot)
    int sets = 2;
    size_t num_threads = 2;

    if (argn >= 2) {        
        int n = atoi(argv[1]);        
        if (n < 1) {
            printf("Iterations must be positive integer\n"); 
            return 1; 
        }        
        iterations = n;
    }
    if (argn >= 3) {        
        int n = atoi(argv[2]);        
        if (n < 1) {
            printf("Width must be positive integer\n"); 
            return 1; 
        }        
        width = n;
    }
    if (argn >= 4) {        
        int n = atoi(argv[3]);        
        if (n < 1) {
            printf("Height must be positive integer\n"); 
            return 1; 
        }
        height = n;
    }
    if (argn >= 5) {        
        int n = atoi(argv[4]);        
        if (n < 0 || n > 2) {
            printf("Specify 0 for Julia set, 1 for Mandlebrot set, 2 for Multibrot \n"); 
            return 1; 
        }
        sets = n;
    }
    if (argn >= 6) {        
        int n = atoi(argv[5]);        
        if (n < 0) {
            printf("Number of threads must be a positive integer \n"); 
            return 1; 
        }
        num_threads = n;
    }
    // seed random number generator
    srand(time(NULL));

    // start the timer
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    Matrix* pixels = matrix_zeros(width, height);
    printf("SETS: %d\n", sets);
    if(sets == 0) {
        julia_set(pixels, iterations, num_threads);
    } else if(sets == 1) {
        mandelbrot(pixels, iterations, num_threads);
    } else if(sets == 2) {
        multibrot(pixels, iterations);
    }
    matrix_to_npy_path("pixels.npy", pixels);
    // get the end and computation time
    clock_gettime(CLOCK_MONOTONIC, &end);
    double time = get_time_diff(&start, &end);
    printf("%f secs\n", time);

    matrix_free(pixels);

	return 0;
}
