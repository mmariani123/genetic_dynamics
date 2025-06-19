#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

// Structure to hold reaction information
typedef struct {
    int* reactants;  // Stoichiometric coefficients for reactants
    int* products;   // Stoichiometric coefficients for products
    double rate;     // Reaction rate constant
    int* affected_species; // Change in molecule numbers when reaction occurs
} Reaction;

// Structure to hold system state
typedef struct {
    double* species;     // Current number of molecules for each species
    int n_species;       // Number of species
    Reaction* reactions; // Array of reactions
    int n_reactions;     // Number of reactions
} System;

// Calculate propensity function for a reaction
double calculate_propensity(Reaction* reaction, double* species, int n_species) {
    double prop = reaction->rate;
    for (int i = 0; i < n_species; i++) {
        if (reaction->reactants[i] > 0) {
            for (int j = 0; j < reaction->reactants[i]; j++) {
                prop *= (species[i] - j);
            }
            prop /= pow(species[i], reaction->reactants[i]);
        }
    }
    return prop;
}

// Generate Poisson random number
int poisson_random(double lambda) {
    if (lambda > 700) {
        // For large lambda, use normal approximation
        double normal = sqrt(-2.0 * log(rand() / (double)RAND_MAX)) *
                       cos(2.0 * M_PI * rand() / (double)RAND_MAX);
        return (int)(lambda + sqrt(lambda) * normal + 0.5);
    }

    double L = exp(-lambda);
    double p = 1.0;
    int k = 0;

    do {
        k++;
        p *= rand() / (double)RAND_MAX;
    } while (p > L);

    return k - 1;
}

// Calculate tau leap step size
double calculate_tau(System* sys, double epsilon, double* propensities) {
    double* mean_change = calloc(sys->n_species, sizeof(double));
    double* variance_change = calloc(sys->n_species, sizeof(double));

    // Calculate mean and variance of change for each species
    for (int i = 0; i < sys->n_reactions; i++) {
        for (int j = 0; j < sys->n_species; j++) {
            mean_change[j] += sys->reactions[i].affected_species[j] * propensities[i];
            variance_change[j] += pow(sys->reactions[i].affected_species[j], 2) * propensities[i];
        }
    }

    // Find minimum tau according to leap condition
    double tau = DBL_MAX;
    for (int i = 0; i < sys->n_species; i++) {
        if (sys->species[i] > 0) {
            double tau1 = epsilon * sys->species[i] / fabs(mean_change[i]);
            double tau2 = pow(epsilon * sys->species[i], 2) / variance_change[i];
            tau = fmin(tau, fmin(tau1, tau2));
        }
    }

    free(mean_change);
    free(variance_change);
    return tau;
}

// Perform one tau leap step
void tau_leap_step(System* sys, double epsilon, double* time) {
    // Calculate propensities
    double* propensities = calloc(sys->n_reactions, sizeof(double));
    for (int i = 0; i < sys->n_reactions; i++) {
        propensities[i] = calculate_propensity(&sys->reactions[i], sys->species, sys->n_species);
    }

    // Calculate tau
    double tau = calculate_tau(sys, epsilon, propensities);

    // Generate number of times each reaction occurs
    int* reaction_counts = calloc(sys->n_reactions, sizeof(int));
    for (int i = 0; i < sys->n_reactions; i++) {
        reaction_counts[i] = poisson_random(propensities[i] * tau);
    }

    // Update species counts
    for (int i = 0; i < sys->n_reactions; i++) {
        if (reaction_counts[i] > 0) {
            for (int j = 0; j < sys->n_species; j++) {
                sys->species[j] += sys->reactions[i].affected_species[j] * reaction_counts[i];
                // Ensure non-negative populations
                if (sys->species[j] < 0) sys->species[j] = 0;
            }
        }
    }

    *time += tau;

    free(propensities);
    free(reaction_counts);
}

// Example usage
void simulate_system(System* sys, double t_final, double epsilon) {
    double t = 0.0;

    while (t < t_final) {
        tau_leap_step(sys, epsilon, &t);

        // Print current state (for demonstration)
        printf("Time: %.6f | Species:", t);
        for (int i = 0; i < sys->n_species; i++) {
            printf(" %.0f", sys->species[i]);
        }
        printf("\n");
    }
}

// Helper function to initialize a reaction
Reaction create_reaction(int n_species, double rate, int* reactants, int* products) {
    Reaction r;
    r.rate = rate;
    r.reactants = malloc(n_species * sizeof(int));
    r.products = malloc(n_species * sizeof(int));
    r.affected_species = malloc(n_species * sizeof(int));

    for (int i = 0; i < n_species; i++) {
        r.reactants[i] = reactants[i];
        r.products[i] = products[i];
        r.affected_species[i] = products[i] - reactants[i];
    }
    return r;
}

const float alpha = 4.65e-1;
const float beta  = 2.22e-3;
const int S0      = 763;
const int I0      =   3;
const int R0      =   0;
const int tmin    =   0;
const int tmax    =  15;

const int react    = 1;
const int prod     = 1;
const int affSpec  = 1;
const double curNum = 1;
const Reaction infection    = { .reactants=&react, .products=&prod, .rate=alpha, .affected_species=&affSpec };
const Reaction removal      = { .reactants=&react, .products=&prod, .rate=beta, .affected_species=&affSpec };
const Reaction reactionsArray[2] = { infection, removal };
const System sysIn    = { .species = &curNum, .n_species = 3, .reactions = &reactionsArray, .n_reactions=2 };

const double finalTime = 15;
const double epsilon = 0.1;

int main(void)
{

simulate_system(&sysIn, finalTime, epsilon);

return 0;

}

// Structure to hold reaction information
//typedef struct {
//    int* reactants;  // Stoichiometric coefficients for reactants
//    int* products;   // Stoichiometric coefficients for products
//    double rate;     // Reaction rate constant
//    int* affected_species; // Change in molecule numbers when reaction occurs
//} Reaction;

// Structure to hold system state
//typedef struct {
//   double* species;     // Current number of molecules for each species
//    int n_species;       // Number of species
//    Reaction* reactions; // Array of reactions
//    int n_reactions;     // Number of reactions
//} System;
