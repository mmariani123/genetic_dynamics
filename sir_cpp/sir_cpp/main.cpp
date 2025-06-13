#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <cmath>
#include <iomanip>
#include <type_traits>

struct SIRState {
    int S, I, R;
    double time;

    SIRState(int s, int i, int r, double t) : S(s), I(i), R(r), time(t) {}

};

class GillespieSSA {
private:
    double beta;      // infection rate
    double gamma;     // recovery rate
    int N;           // total population
    int num_threads;

public:
    GillespieSSA(double b, double g, int population, int threads = 4)
        : beta(b), gamma(g), N(population), num_threads(threads) {
        omp_set_num_threads(num_threads);
    }

    // Single simulation run
    std::vector<SIRState> simulate_single(int S0, int I0, int R0, double t_max,
                                         std::mt19937& rng) {
        std::vector<SIRState> trajectory;

        int S = S0, I = I0, R = R0;
        double t = 0.0;

        std::uniform_real_distribution<double> uniform(0.0, 1.0);

        trajectory.push_back(SIRState(S, I, R, t));

        while (t < t_max && I > 0) {
            // Calculate propensities
            double a1 = beta * S * I / static_cast<double>(N);  // infection
            double a2 = gamma * I;                              // recovery
            double a_total = a1 + a2;

            if (a_total == 0.0) break;

            // Generate random numbers
            double r1 = uniform(rng);
            double r2 = uniform(rng);

            // Calculate time to next event
            double tau = -log(r1) / a_total;
            t += tau;

            if (t > t_max) break;

            // Determine which reaction occurs
            if (r2 * a_total < a1) {
                // Infection event
                S--;
                I++;
            } else {
                // Recovery event
                I--;
                R++;
            }

            trajectory.push_back(SIRState(S, I, R, t));
        }

        return trajectory;
    }

    // Parallel ensemble simulation

    //omp_set_num_threads(1);
    //std::cout << "before parallel section: " << std::endl;
    //std::cout << "Num threads = " << omp_get_num_threads() << std::endl;
    //std::cout << "Max threads = " << omp_get_max_threads() << std::endl;

    //https://stackoverflow.com/questions/1448318/omp-parallel-vs-omp-parallel-for

    std::vector<std::vector<SIRState>> simulate_ensemble(int S0, int I0, int R0,
                                                        double t_max, int num_runs) {
        std::vector<std::vector<SIRState>> all_trajectories(num_runs);

        std::cout << "Running " << num_runs << " simulations with "
                  << num_threads << " threads..." << std::endl;

        auto start_time = std::chrono::high_resolution_clock::now();

        #pragma omp parallel
        {
            // Each thread gets its own random number generator
            std::random_device rd;
            std::mt19937 rng(rd() + omp_get_thread_num());

            #pragma omp for schedule(dynamic)
            for (int run = 0; run < num_runs; run++) {
                all_trajectories[run] = simulate_single(S0, I0, R0, t_max, rng);

                if (run % 100 == 0 && omp_get_thread_num() == 0) {
                    std::cout << "Completed " << run << "/" << num_runs
                              << " simulations" << std::endl;
                }
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        //auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
        double duration = std::chrono::duration<double>(end_time - start_time).count();

        //std::cout << std::fixed << iomanip::setprecision(20);
        //double duration = start_time.time_since_epoch();


        //using namespace std::chrono;

        //std::cout << std::boolalpha;
        //std::cout << std::experimental::is_same_v<system_clock, high_resolution_clock> << std::endl;

        std::cout << "\nGillespie SSA Algorithm Timing: "
                  << static_cast<double>(duration) << " milliseconds" << std::endl;
                  //   std::cout << duration << "\n";
        std::cout << "Average time per simulation: "
                  << static_cast<double>(static_cast<double>(duration) / num_runs)
                  //<< static_cast<double>(duration) / num_runs
                  << " milliseconds" << std::endl;

        return all_trajectories;
    }

    // Calculate ensemble statistics at regular time intervals
    void calculate_ensemble_stats(const std::vector<std::vector<SIRState>>& trajectories,
                                 double t_max, double dt, const std::string& filename) {
        std::vector<double> time_points;
        std::vector<double> S_mean, I_mean, R_mean;
        std::vector<double> S_std, I_std, R_std;

        int num_points = static_cast<int>(t_max / dt) + 1;

        for (int i = 0; i < num_points; i++) {
            double t = i * dt;
            time_points.push_back(t);

            std::vector<int> S_values, I_values, R_values;

            // Interpolate each trajectory at time t
            for (const auto& traj : trajectories) {
                if (traj.empty()) continue;

                // Find the state at time t
                auto it = std::lower_bound(traj.begin(), traj.end(), t,
                    [](const SIRState& state, double time) {
                        return state.time < time;
                    });

                if (it == traj.end()) {
                    // Use last state
                    S_values.push_back(traj.back().S);
                    I_values.push_back(traj.back().I);
                    R_values.push_back(traj.back().R);
                } else if (it == traj.begin()) {
                    // Use first state
                    S_values.push_back(traj.front().S);
                    I_values.push_back(traj.front().I);
                    R_values.push_back(traj.front().R);
                } else {
                    // Use the previous state (step function)
                    --it;
                    S_values.push_back(it->S);
                    I_values.push_back(it->I);
                    R_values.push_back(it->R);
                }
            }

            // Calculate means
            double S_sum = 0, I_sum = 0, R_sum = 0;
            for (size_t j = 0; j < S_values.size(); j++) {
                S_sum += S_values[j];
                I_sum += I_values[j];
                R_sum += R_values[j];
            }

            double n = static_cast<double>(S_values.size());
            S_mean.push_back(S_sum / n);
            I_mean.push_back(I_sum / n);
            R_mean.push_back(R_sum / n);

            // Calculate standard deviations
            double S_var = 0, I_var = 0, R_var = 0;
            for (size_t j = 0; j < S_values.size(); j++) {
                S_var += (S_values[j] - S_mean.back()) * (S_values[j] - S_mean.back());
                I_var += (I_values[j] - I_mean.back()) * (I_values[j] - I_mean.back());
                R_var += (R_values[j] - R_mean.back()) * (R_values[j] - R_mean.back());
            }

            S_std.push_back(sqrt(S_var / n));
            I_std.push_back(sqrt(I_var / n));
            R_std.push_back(sqrt(R_var / n));
        }

        // Write results to file
        std::ofstream file(filename);
        file << "Time,S_mean,I_mean,R_mean,S_std,I_std,R_std\n";
        for (size_t i = 0; i < time_points.size(); i++) {
            file << time_points[i] << "," << S_mean[i] << "," << I_mean[i] << ","
                 << R_mean[i] << "," << S_std[i] << "," << I_std[i] << ","
                 << R_std[i] << "\n";
        }
        file.close();

        std::cout << "Ensemble statistics saved to " << filename << std::endl;
    }

    // Save individual trajectories for visualization
    void save_trajectories(const std::vector<std::vector<SIRState>>& trajectories,
                          const std::string& filename, int max_trajectories = 50) {
        std::ofstream file(filename);
        file << "Time,S,I,R,Trajectory\n";

        int traj_count = 0;
        for (const auto& traj : trajectories) {
            if (traj_count >= max_trajectories) break;

            for (const auto& state : traj) {
                file << state.time << "," << state.S << "," << state.I << ","
                     << state.R << "," << traj_count << "\n";
            }
            traj_count++;
        }
        file.close();

        std::cout << "Individual trajectories saved to " << filename << std::endl;
    }

    // Generate Python visualization script
    void generate_visualization_script(const std::string& stats_file,
                                      const std::string& traj_file) {
        std::ofstream script("visualize_sir.py");
        script << R"(import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read ensemble statistics
stats_df = pd.read_csv(')" << stats_file << R"(')

# Read individual trajectories
traj_df = pd.read_csv(')" << traj_file << R"(')

# Create figure with subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Plot 1: Ensemble statistics with error bars
ax1.errorbar(stats_df['Time'], stats_df['S_mean'], yerr=stats_df['S_std'],
             label='Susceptible', alpha=0.7, capsize=2)
ax1.errorbar(stats_df['Time'], stats_df['I_mean'], yerr=stats_df['I_std'],
             label='Infected', alpha=0.7, capsize=2)
ax1.errorbar(stats_df['Time'], stats_df['R_mean'], yerr=stats_df['R_std'],
             label='Recovered', alpha=0.7, capsize=2)

ax1.set_xlabel('Time')
ax1.set_ylabel('Population')
ax1.set_title('SIR Model - Ensemble Statistics (Mean ± Std)')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: Individual trajectories
for traj_id in traj_df['Trajectory'].unique()[:20]:  # Show first 20 trajectories
    traj_data = traj_df[traj_df['Trajectory'] == traj_id]
    ax2.plot(traj_data['Time'], traj_data['S'], 'b-', alpha=0.3, linewidth=0.5)
    ax2.plot(traj_data['Time'], traj_data['I'], 'r-', alpha=0.3, linewidth=0.5)
    ax2.plot(traj_data['Time'], traj_data['R'], 'g-', alpha=0.3, linewidth=0.5)

# Add legend for individual trajectories
from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], color='b', label='Susceptible'),
                  Line2D([0], [0], color='r', label='Infected'),
                  Line2D([0], [0], color='g', label='Recovered')]
ax2.legend(handles=legend_elements)

ax2.set_xlabel('Time')
ax2.set_ylabel('Population')
ax2.set_title('SIR Model - Individual Trajectories')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('sir_gillespie_results.png', dpi=300, bbox_inches='tight')
plt.show()

print("Visualization saved as 'sir_gillespie_results.png'")
)";
        script.close();

        std::cout << "Python visualization script generated: visualize_sir.py" << std::endl;
        std::cout << "Run 'python visualize_sir.py' to generate plots" << std::endl;
    }
};

int main() {
    // SIR model parameters
    double beta = 2.22e-3;    // infection rate
    double gamma = 4.65e-1;   // recovery rate
    int N = 766;              // total population
    int S0 = 763;             // initial susceptible
    int I0 = 3;               // initial infected
    int R0 = 0;               // initial recovered
    double t_max = 15.0;      // simulation time
    int num_runs = 20000;     // number of ensemble runs
    int num_threads = omp_get_max_threads();

    std::cout << "=== Gillespie SSA for SIR Epidemiological Model ===" << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "  Population (N): " << N << std::endl;
    std::cout << "  Infection rate (β): " << beta << std::endl;
    std::cout << "  Recovery rate (γ): " << gamma << std::endl;
    std::cout << "  Basic reproduction number (R₀): " << beta/gamma << std::endl;
    std::cout << "  Initial conditions: S=" << S0 << ", I=" << I0 << ", R=" << R0 << std::endl;
    std::cout << "  Simulation time: " << t_max << std::endl;
    std::cout << "  Number of runs: " << num_runs << std::endl;
    std::cout << "  Available threads: " << num_threads << std::endl;
    std::cout << std::endl;

    // Create Gillespie SSA simulator
    GillespieSSA simulator(beta, gamma, N, num_threads);

    // Run ensemble simulations
    auto trajectories = simulator.simulate_ensemble(S0, I0, R0, t_max, num_runs);

    std::cout << "\nPost-processing results..." << std::endl;

    // Calculate and save ensemble statistics
    simulator.calculate_ensemble_stats(trajectories, t_max, 1.0, "sir_ensemble_stats.csv");

    // Save individual trajectories for visualization
    simulator.save_trajectories(trajectories, "sir_trajectories.csv", 50);

    // Generate visualization script
    simulator.generate_visualization_script("sir_ensemble_stats.csv", "sir_trajectories.csv");

    std::cout << "\n=== Simulation Complete ===" << std::endl;
    std::cout << "Files generated:" << std::endl;
    std::cout << "  - sir_ensemble_stats.csv: Ensemble statistics" << std::endl;
    std::cout << "  - sir_trajectories.csv: Individual trajectories" << std::endl;
    std::cout << "  - visualize_sir.py: Python visualization script" << std::endl;

    return 0;
}
