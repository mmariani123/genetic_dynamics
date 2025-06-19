#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Optional OpenMP support
#ifdef _OPENMP
#include <omp.h>
#endif

// Structure to hold SIR model parameters
typedef struct {
    double beta;   // infection rate
    double gamma;  // recovery rate
    int S;        // number of susceptible individuals
    int I;        // number of infected individuals
    int R;        // number of recovered individuals
    double t;     // current time
} SIRState;

// Structure to hold simulation results
typedef struct {
    double *times;
    int *S_counts;
    int *I_counts;
    int *R_counts;
    int length;
} SimulationResult;

// Function to generate random number between 0 and 1
double random_uniform() {
    return (double)rand() / RAND_MAX;
}

// Function to run a single trajectory of the Gillespie algorithm
SimulationResult run_gillespie_trajectory(SIRState initial_state, double t_max) {
    // Initialize state
    SIRState state = initial_state;

    // Allocate memory for results (with initial capacity)
    int capacity = 1000;
    double *times = malloc(capacity * sizeof(double));
    int *S_counts = malloc(capacity * sizeof(int));
    int *I_counts = malloc(capacity * sizeof(int));
    int *R_counts = malloc(capacity * sizeof(int));

    // Store initial state
    int idx = 0;
    times[idx] = state.t;
    S_counts[idx] = state.S;
    I_counts[idx] = state.I;
    R_counts[idx] = state.R;
    idx++;

    // Main simulation loop
    while (state.t < t_max && state.I > 0) {
        // Calculate propensities
        double a1 = state.beta * state.S * state.I;  // infection
        double a2 = state.gamma * state.I;          // recovery
        double a_total = a1 + a2;

        // Break if no reactions possible
        if (a_total == 0) break;

        // Generate random numbers
        double r1 = random_uniform();
        double r2 = random_uniform();

        // Update time
        state.t += -log(r1) / a_total;

        // Choose and execute reaction
        if (r2 < a1/a_total) {
            // Infection
            state.S--;
            state.I++;
        } else {
            // Recovery
            state.I--;
            state.R++;
        }

        // Store state
        if (idx >= capacity) {
            capacity *= 2;
            times = realloc(times, capacity * sizeof(double));
            S_counts = realloc(S_counts, capacity * sizeof(int));
            I_counts = realloc(I_counts, capacity * sizeof(int));
            R_counts = realloc(R_counts, capacity * sizeof(int));
        }

        times[idx] = state.t;
        S_counts[idx] = state.S;
        I_counts[idx] = state.I;
        R_counts[idx] = state.R;
        idx++;
    }

    // Create result structure
    SimulationResult result = {
        .times = times,
        .S_counts = S_counts,
        .I_counts = I_counts,
        .R_counts = R_counts,
        .length = idx
    };

    return result;
}

// Function to run multiple trajectories in serial
void run_serial_simulations(int n_trajectories, SIRState initial_state,
                          double t_max, const char* output_file) {
    // Open output file
    FILE *fp = fopen(output_file, "w");
    if (fp == NULL) {
        printf("Error opening file!\n");
        return;
    }

    fprintf(fp, "Trajectory,Time,S,I,R\n");

    printf("Running %d trajectories in serial\n", n_trajectories);

    // Set random seed
    srand(time(NULL));

    // Run trajectories sequentially
    for (int i = 0; i < n_trajectories; i++) {
        // Run single trajectory
        SimulationResult result = run_gillespie_trajectory(initial_state, t_max);

        // Write results to file
        for (int j = 0; j < result.length; j++) {
            fprintf(fp, "%d,%f,%d,%d,%d\n",
                    i, result.times[j],
                    result.S_counts[j],
                    result.I_counts[j],
                    result.R_counts[j]);
        }

        // Free memory
        free(result.times);
        free(result.S_counts);
        free(result.I_counts);
        free(result.R_counts);

        // Print progress
        if ((i + 1) % 10 == 0) {
            printf("Completed %d trajectories\n", i + 1);
        }
    }

    fclose(fp);
}

#ifdef _OPENMP
// Function to run multiple trajectories in parallel (only compiled if OpenMP is available)
void run_parallel_simulations(int n_trajectories, SIRState initial_state,
                            double t_max, const char* output_file) {
    // Open output file
    FILE *fp = fopen(output_file, "w");
    if (fp == NULL) {
        printf("Error opening file!\n");
        return;
    }

    fprintf(fp, "Trajectory,Time,S,I,R\n");

    // Set number of threads
    int n_threads = omp_get_max_threads();
    printf("Running %d trajectories using %d threads\n", n_trajectories, n_threads);

    // Parallel region
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_trajectories; i++) {
        // Set different random seed for each thread
        unsigned int seed = time(NULL) + omp_get_thread_num();
        srand(seed);

        // Run single trajectory
        SimulationResult result = run_gillespie_trajectory(initial_state, t_max);

        // Write results to file (use critical section to avoid conflicts)
        #pragma omp critical
        {
            for (int j = 0; j < result.length; j++) {
                fprintf(fp, "%d,%f,%d,%d,%d\n",
                        i, result.times[j],
                        result.S_counts[j],
                        result.I_counts[j],
                        result.R_counts[j]);
            }

            if ((i + 1) % 10 == 0) {
                printf("Completed %d trajectories\n", i + 1);
            }
        }

        // Free memory
        free(result.times);
        free(result.S_counts);
        free(result.I_counts);
        free(result.R_counts);
    }

    fclose(fp);
}
#endif

int main() {
    // Set simulation parameters
    SIRState initial_state = {
        .beta = 0.3,    // infection rate
        .gamma = 0.1,   // recovery rate
        .S = 999,       // initial susceptible population
        .I = 1,         // initial infected population
        .R = 0,         // initial recovered population
        .t = 0.0        // initial time
    };

    double t_max = 100.0;        // maximum simulation time
    int n_trajectories = 100;    // number of trajectories to simulate

    // Run simulations
    #ifdef _OPENMP
        run_parallel_simulations(n_trajectories, initial_state, t_max, "sir_results.csv");
    #else
        run_serial_simulations(n_trajectories, initial_state, t_max, "sir_results.csv");
    #endif

    printf("Simulation complete. Results written to sir_results.csv\n");

    return 0;
}
