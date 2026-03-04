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

void recursive_fft(double complex *x, int N) {
    // Basfall
    if (N <= 1) return;

    // 1. Dela upp i jämna och udda index
    // Här skapar vi nya minnesblock, precis som "new" i Java
    double complex *even = malloc((N / 2) * sizeof(double complex));
    double complex *odd = malloc((N / 2) * sizeof(double complex));

    for (int i = 0; i < N / 2; i++) {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    // 2. Rekursiva anrop
    recursive_fft(even, N / 2);
    recursive_fft(odd, N / 2);

    // 3. Slå ihop resultaten (Butterfly)
    for (int k = 0; k < N / 2; k++) {
        double ang = -2.0 * PI * k / N;
        double complex w = cexp(I * ang); // Twiddle factor
        
        double complex t = w * odd[k];
        
        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }

    // Mycket viktigt i C: Frigör minnet vi skapade!
    // Java gör detta automatiskt via Garbage Collector, men här gör vi det manuellt.
    free(even);
    free(odd);
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

    recursive_fft(test_signal, N);

    clock_gettime(CLOCK_MONOTONIC, &end);
    // 3. Tidmätning SLUT

    double time_taken = (end.tv_sec - start.tv_sec) * 1e3 + (end.tv_nsec - start.tv_nsec) / 1e6;
    printf("C FFT (N=%d) klar pa %.4f ms\n", N, time_taken);

    // 4. Spara till CSV för analys
    write_csv("recursive_fft_results_c.csv", N, fs, test_signal);

    free(test_signal);
    return 0;
}