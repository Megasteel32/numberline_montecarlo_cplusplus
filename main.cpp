#include <iostream>
#include <random>
#include <vector>
#include <chrono>
#include <thread>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>

struct SimulationResult {
    double min_sum;
    double max_sum;
};

SimulationResult simulate_points(int num_simulations) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    double min_sum = 0.0;
    double max_sum = 0.0;

    for (int i = 0; i < num_simulations; ++i) {
        auto point1 = dis(gen);
        auto point2 = dis(gen);
        min_sum += std::min(point1, point2);
        max_sum += std::max(point1, point2);
    }

    return {min_sum, max_sum};
}

std::pair<double, double> parallel_simulate(int total_simulations, int num_threads) {
    if (num_threads == 1) {
        auto result = simulate_points(total_simulations);
        return {result.min_sum / total_simulations, result.max_sum / total_simulations};
    }

    const int chunk_size = total_simulations / num_threads;
    std::vector<int> chunks(num_threads, chunk_size);
    chunks.back() += total_simulations % num_threads;

    std::vector<SimulationResult> results(num_threads);
    std::vector<std::thread> threads;

    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back([&results, i, &chunks]() {
            results[i] = simulate_points(chunks[i]);
        });
    }

    for (auto& thread : threads) {
        thread.join();
    }

    auto total_result = std::accumulate(results.begin(), results.end(), SimulationResult{0.0, 0.0},
        [](SimulationResult acc, const SimulationResult& res) {
            return SimulationResult{acc.min_sum + res.min_sum, acc.max_sum + res.max_sum};
        });

    return {total_result.min_sum / total_simulations, total_result.max_sum / total_simulations};
}

std::string format_double(double value, int precision = 8) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

void print_usage() {
    std::cout << "Usage: ./program [options]\n"
              << "Options:\n"
              << "  -s <num>    Set the number of simulations (default: 100000000)\n"
              << "  -t <num>    Set the number of threads (default: 1)\n";
}

int main(int argc, char* argv[]) {
    int total_simulations = 100'000'000;
    int num_threads = 1;

    std::unordered_map<std::string, std::string> args;
    for (int i = 1; i < argc; i += 2) {
        if (i + 1 < argc) {
            args[argv[i]] = argv[i + 1];
        } else {
            std::cerr << "Error: Missing value for " << argv[i] << "\n";
            print_usage();
            return 1;
        }
    }

    if (args.count("-s")) {
        try {
            total_simulations = std::stoi(args["-s"]);
            if (total_simulations <= 0) throw std::runtime_error("Number of simulations must be positive");
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid argument for number of simulations. Using default value of 100,000,000.\n";
            total_simulations = 100'000'000;
        }
    }

    if (args.count("-t")) {
        try {
            num_threads = std::stoi(args["-t"]);
            if (num_threads < 1) throw std::runtime_error("Number of threads must be at least 1");
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid argument for number of threads. Using 1 thread.\n";
            num_threads = 1;
        }
    }

    std::cout << "Running " << total_simulations << " simulations using " << num_threads << " thread(s)...\n";

    auto start_time = std::chrono::high_resolution_clock::now();

    auto [expected_min, expected_max] = parallel_simulate(total_simulations, num_threads);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration<double>(end_time - start_time).count();

    std::cout << "\n\nSimulation completed in " << format_double(elapsed_time, 2) << " seconds\n"
              << "Number of simulations: " << total_simulations << "\n"
              << "Number of threads used: " << num_threads << "\n"
              << "Expected value of minimum point: " << format_double(expected_min) << "\n"
              << "Expected value of maximum point: " << format_double(expected_max) << "\n\n"
              << "Theoretical expected value of minimum: " << format_double(1.0/3.0) << "\n"
              << "Theoretical expected value of maximum: " << format_double(2.0/3.0) << "\n"
              << "Difference from theoretical (min): " << format_double(std::abs(expected_min - 1.0/3.0)) << "\n"
              << "Difference from theoretical (max): " << format_double(std::abs(expected_max - 2.0/3.0)) << "\n";

    return 0;
}