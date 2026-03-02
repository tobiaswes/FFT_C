#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define PI 3.14159265358979323846

// Struktur för att hålla koll på WAV-headern (förenklad)
typedef struct {
    char riff[4];
    int size;
    char wave[4];
    char fmt[4];
    int fmt_len;
    short format_tag;
    short channels;
    int sample_rate;
    int byte_rate;
    short block_align;
    short bits_per_sample;
    char data[4];
    int data_len;
} WavHeader;

// Läs WAV-fil (Hanterar 16-bit PCM, Mono/Stereo)
double complex* read_wav(const char* filename, int N, int *sample_rate) {
    FILE *f = fopen(filename, "rb");
    if (!f) return NULL;

    WavHeader header;
    fread(&header, sizeof(WavHeader), 1, f);
    *sample_rate = header.sample_rate;

    double complex *signal = malloc(N * sizeof(double complex));
    short temp_sample;

    for (int i = 0; i < N; i++) {
        // Läs vänster kanal (eller mono)
        if (fread(&temp_sample, sizeof(short), 1, f) == 1) {
            signal[i] = (double)temp_sample + 0.0 * I;
            
            // Om Stereo: hoppa över höger kanal för att hålla rätt frekvens
            if (header.channels == 2) {
                fseek(f, sizeof(short), SEEK_CUR);
            }
        } else {
            signal[i] = 0.0 + 0.0 * I; // Padding
        }
    }

    fclose(f);
    return signal;
}

// Skriv resultat till CSV
void write_csv(const char* filename, int N, int fs, double complex *result) {
    FILE *f = fopen(filename, "w");
    fprintf(f, "Frequency,Magnitude\n");

    for (int i = 0; i < N / 2; i++) {
        double freq = (double)i * fs / N;
        double mag = cabs(result[i]); // cabs beräknar magnituden för komplexa tal
        fprintf(f, "%.2f,%.6f\n", freq, mag);
    }

    fclose(f);
}

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
    int N = 4096; 
    int fs;
    
    // 1. Läs in ljudfilen
    double complex *test_signal = read_wav("A5_test.wav", N, &fs);
    if (!test_signal) {
        printf("Kunde inte lasa filen!\n");
        return 1;
    }

    // 2. Tidmätning START
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    iterative_fft(test_signal, N);

    clock_gettime(CLOCK_MONOTONIC, &end);
    // 3. Tidmätning SLUT

    double time_taken = (end.tv_sec - start.tv_sec) * 1e3 + (end.tv_nsec - start.tv_nsec) / 1e6;
    printf("C FFT (N=%d) klar pa %.4f ms\n", N, time_taken);

    // 4. Spara till CSV för analys
    write_csv("fft_results_c.csv", N, fs, test_signal);

    free(test_signal);
    return 0;
}