#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>

#define PI 3.14159265358979323846

// Global variabel för att spåra allokerat minne i den rekursiva funktionen
long long current_iteration_memory = 0;

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

    // --- HÄR ÄR FIXEN ---
    // Vi räknar ut hur mycket minne dessa två malloc kommer att kräva
    size_t size_per_array = (N / 2) * sizeof(double complex);
    current_iteration_memory += (size_per_array * 2);

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
    int iterations = 100;
    
    double *times = malloc(iterations * sizeof(double));
    long long *mem_usages = malloc(iterations * sizeof(long long));

    double complex *original_signal = read_wav("A5_test.wav", N, &fs);
    if (!original_signal) return 1;

    printf("Startar %d iterationer i C...\n", iterations);

    for (int i = 0; i < iterations; i++) {
        // Skapa en kopia inför varje körning (motsvarar clone() i Java)
        double complex *working_copy = malloc(N * sizeof(double complex));
        memcpy(working_copy, original_signal, N * sizeof(double complex));

        current_iteration_memory = 0; // Nollställ minnesmätaren för denna runda

        struct timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);

        recursive_fft(working_copy, N);

        clock_gettime(CLOCK_MONOTONIC, &end);
        
        double time_ms = (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_nsec - start.tv_nsec) / 1000000.0;
        
        times[i] = time_ms;
        mem_usages[i] = current_iteration_memory;

        free(working_copy);
    }

    // Beräkna statistik
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

    // Presentation (Samma format som Java)
    printf("\n====================================================\n");
    printf("   RESULTAT: C REKURSIV FFT (N=%d)\n", N);
    printf("====================================================\n");
    printf("ANTAL KORNINGAR: %d\n\n", iterations);
    
    printf("TID (Millisekunder)\n");
    printf("  Medel: %.4f ms\n", sum_time / iterations);
    printf("  Min:   %.4f ms\n", min_time);
    printf("  Max:   %.4f ms\n\n", max_time);

    printf("MINNE (Kilobytes allokerat under korning)\n");
    printf("  Medel: %.2f KB\n", (sum_mem / iterations) / 1024.0);
    printf("  Min:   %.2f KB\n", min_mem / 1024.0);
    printf("  Max:   %.2f KB\n", max_mem / 1024.0);
    printf("====================================================\n");

    // Spara sista körningen till CSV för validering
    recursive_fft(original_signal, N);
    write_csv("recursive_fft_results_c.csv", N, fs, original_signal);

    free(original_signal);
    free(times);
    free(mem_usages);
    return 0;
}