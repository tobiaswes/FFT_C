#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define PI 3.14159265358979323846

void iterative_fft(double complex *x, int N) {
    // 1. BIT-REVERSAL PERMUTATION
    int j = 0;
    for (int i = 1; i < N; i++) {
        int bit = N >> 1;
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j) {
            double complex temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
    }

    // 2. BUTTERFLY BERÄKNINGAR
    for (int length = 2; length <= N; length <<= 1) {
        double ang = -2.0 * PI / length;
        double complex w_m = cexp(I * ang); // I är den imaginära enheten i C
        
        for (int i = 0; i < N; i += length) {
            double complex w = 1.0 + 0.0 * I;
            for (int k = i; k < i + length / 2; k++) {
                // Butterfly-matematik
                double complex u = x[k];
                double complex t = w * x[k + length / 2];
                
                x[k] = u + t;
                x[k + length / 2] = u - t;
                
                w *= w_m;
            }
        }
    }
}

int main() {
    int N = 1024;
    double complex *test_signal = malloc(N * sizeof(double complex));

    // Skapa testsignal (sinusvåg)
    for (int i = 0; i < N; i++) {
        test_signal[i] = sin(2 * PI * i / 32.0) + 0.0 * I;
    }

    // Tidmätning
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    iterative_fft(test_signal, N);

    clock_gettime(CLOCK_MONOTONIC, &end);

    // Beräkna tid i millisekunder
    double time_taken = (end.tv_sec - start.tv_sec) * 1e3 + (end.tv_nsec - start.tv_nsec) / 1e6;
    printf("FFT klar pa %.3f ms\n", time_taken);

    free(test_signal);
    return 0;
}