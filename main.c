#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>

#define PI 3.14159265358979323846

long long current_iteration_memory = 0;

typedef struct {
    char riff[4]; int size; char wave[4]; char fmt[4]; int fmt_len;
    short format_tag; short channels; int sample_rate; int byte_rate;
    short block_align; short bits_per_sample; char data[4]; int data_len;
} WavHeader;

// --- HJÄLPFUNKTIONER ---

double complex* read_wav(const char* filename, int N, int *sample_rate) {
    FILE *f = fopen(filename, "rb");
    if (!f) return NULL;
    WavHeader header;
    fread(&header, sizeof(WavHeader), 1, f);
    *sample_rate = header.sample_rate;
    double complex *signal = malloc(N * sizeof(double complex));
    short temp;
    for (int i = 0; i < N; i++) {
        if (fread(&temp, sizeof(short), 1, f) == 1) {
            signal[i] = (double)temp + 0.0 * I;
            if (header.channels == 2) fseek(f, sizeof(short), SEEK_CUR);
        } else {
            signal[i] = 0.0 + 0.0 * I;
        }
    }
    fclose(f);
    return signal;
}

void write_csv(const char* filename, int N, int fs, double complex *result) {
    FILE *f = fopen(filename, "w");
    fprintf(f, "Frequency,Magnitude\n");
    for (int i = 0; i < N / 2; i++) {
        fprintf(f, "%.2f,%.6f\n", (double)i * fs / N, cabs(result[i]));
    }
    fclose(f);
}

// --- ALGORITMEN (ITERATIV) ---

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
        double complex w_m = cexp(I * ang);
        
        for (int i = 0; i < N; i += length) {
            double complex w = 1.0 + 0.0 * I;
            for (int k = i; k < i + length / 2; k++) {
                double complex u = x[k];
                double complex t = w * x[k + length / 2];
                x[k] = u + t;
                x[k + length / 2] = u - t;
                w *= w_m;
            }
        }
    }
}

// --- MAIN MED STATISTIK ---

int main() {
    int N = 4096; 
    int fs;
    int iterations = 100;
    
    double *times = malloc(iterations * sizeof(double));
    long long *mem_usages = malloc(iterations * sizeof(long long));

    // 1. Läs in originalet en gång
    double complex *original_signal = read_wav("A5_test.wav", N, &fs);
    if (!original_signal) {
        printf("Kunde inte lasa filen!\n");
        return 1;
    }

    printf("Startar %d iterationer av Iterativ FFT i C (N=%d)...\n", iterations, N);

    // 2. TEST-LOOP
    for (int i = 0; i < iterations; i++) {
        // Skapa en kopia för varje körning (In-place ändrar arrayen)
        double complex *working_copy = malloc(N * sizeof(double complex));
        memcpy(working_copy, original_signal, N * sizeof(double complex));

        struct timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);

        // KÖR ALGORITMEN
        iterative_fft(working_copy, N);

        clock_gettime(CLOCK_MONOTONIC, &end);
        
        double time_ms = (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_nsec - start.tv_nsec) / 1000000.0;
        times[i] = time_ms;
        mem_usages[i] = current_iteration_memory; // Spara mätningen

        free(working_copy);
    }

    // 3. BERÄKNA STATISTIK
    double sum_time = 0, min_time = times[0], max_time = times[0];
    long long sum_mem = 0, min_mem = mem_usages[0], max_mem = mem_usages[0];

    for (int i = 0; i < iterations; i++) {
        sum_time += times[i];
        if (times[i] < min_time) min_time = times[i];
        if (times[i] > max_time) max_time = times[i];

        sum_mem += mem_usages[i];
        if (mem_usages[i] < min_mem) min_mem = mem_usages[i];
        if (mem_usages[i] > max_mem) max_mem = mem_usages[i];
    }

    // 4. PRESENTERA RESULTAT
    printf("\n====================================================\n");
    printf("   RESULTAT: C ITERATIV FFT (N=%d)\n", N);
    printf("====================================================\n");
    printf("ANTAL KORNINGAR: %d\n\n", iterations);
    
    printf("TID (Millisekunder)\n");
    printf("  Medel: %.4f ms\n", sum_time / iterations);
    printf("  Min:   %.4f ms\n", min_time);
    printf("  Max:   %.4f ms\n\n", max_time);

    printf("MINNE (Kilobytes allokerat under körning)\n");
    printf("  Medel: %.2f KB\n", (double)sum_mem / iterations / 1024.0);
    printf("  Min:   %.2f KB\n", (double)min_mem / 1024.0);
    printf("  Max:   %.2f KB\n", (double)max_mem / 1024.0);
    printf("====================================================\n");

    // Spara en CSV för validering
    iterative_fft(original_signal, N);
    write_csv("fft_results_c.csv", N, fs, original_signal);

    free(original_signal);
    free(times);
    free(mem_usages);
    return 0;
}