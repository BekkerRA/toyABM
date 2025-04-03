//
// Created by Rebecca Bekker on 3/26/25.
//
#include "Environment.h"
#include "Cell.h"
#include "RNG.h"
#include <array>
#include <iostream>
#include <set>
#include <vector>
#include <thread>


// Environment::Environment(int seed,int xD, int yD, double prob_div, double prob_die, double prob_mig): rng(seed) {
//     divProb = prob_div;
//     deathProb = prob_die;
//     moveProb = prob_mig;
//     xDim = xD;
//     yDim = yD;
//
//     rng = RNG(seed);
//     // Generate a sequence of unique seeds for each thread
//     std::vector<unsigned int> thread_seeds(omp_get_max_threads());
//     std::generate(thread_seeds.begin(), thread_seeds.end(), std::mt19937(seed));
//
//     // Create a thread-local variable to store the thread's seed
//     thread_local unsigned int thread_seed;
//
// #pragma omp parallel
//     {
//         thread_seed = thread_seeds[omp_get_thread_num()]; // Assign a unique seed to each thread
//     }
//
// };

Environment::Environment(int seed, int xD, int yD, double prob_div, double prob_die, double prob_mig)
    : rng(seed), xDim(xD), yDim(yD), divProb(prob_div), deathProb(prob_die), moveProb(prob_mig) {
    // No need to manually assign thread-local seeds anymore
}

void Environment::initializeCells(int numCells) {
    for (int i = 0; i < numCells; i++) {
        Cell newCell = Cell({rng.uniform(0,xDim),rng.uniform(0,yDim)}, 1);
        cell_list.push_back(newCell);
    }

    std::cout<<"Model initialized." <<std::endl;
    std::cout<<"Time: 0 Cells: "<<cell_list.size()<<std::endl;
}

void Environment::update_parallel() {
    double totalProb = divProb + deathProb + moveProb;

    int master_seed = 100;
    std::vector<unsigned int> thread_seeds(omp_get_max_threads());
    std::mt19937 master_gen(master_seed);  // Master RNG to generate per-thread seeds
    std::generate(thread_seeds.begin(), thread_seeds.end(), [&master_gen]() {
        return master_gen();
    });

#pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        std::mt19937 local_gen(thread_seeds[thread_id]); //create a local generator
        std::uniform_real_distribution<double> dist(0.0, 1.0); //create a local distribution

#pragma omp for schedule (static)
        for (size_t i = 0; i < cell_list.size(); ++i){
            Cell& cell = cell_list[i];
            double prob = dist(local_gen); //use local generator and distribution.

            if (prob < divProb / totalProb) {
                divide(cell);
            } else if (prob < (divProb + moveProb) / totalProb) {
                move(cell);
            } else {
                die(cell);
            }
        }
    }
}


void Environment::update_sequential() {

    double totalProb = divProb + deathProb + moveProb;
    for (Cell& cell : cell_list) {
        double prob = rng.uniform(0,1);
        if (prob < divProb / totalProb) {
            divide(cell);
        } else if (prob < (divProb + moveProb) / totalProb) {
            move(cell);
        } else {
            die(cell);
        }
    }

}


void Environment::simulate(int finalTimeStep) {
    for (int i = 1; i < finalTimeStep;++i) {
        shuffleCells();
        update_parallel();
        //update_sequential();
        cleanCells();
        std::cout<<"Time: "<< i << " Cells: "<<cell_list.size()<<std::endl;
    }
}


void Environment::cleanCells() {
    cell_list.erase(std::remove_if(cell_list.begin(), cell_list.end(),
        [](const Cell& cell) {return cell.state == -1;}), cell_list.end());
}

void Environment::shuffleCells() {
    std::shuffle(cell_list.begin(), cell_list.end(), rng.getGenerator());
}

void Environment::divide(Cell& cell) {
    if (cell.state !=-1) {
        Cell daughterCell = Cell({rng.uniform(0,xDim),rng.uniform(0,yDim)}, cell.state);
        cell_list.push_back(daughterCell);
    }
}

void Environment::move(Cell& cell) {
    if (cell.state !=-1) {
        cell.x[0] = cell.x[0] + rng.uniform(-1,1);
        cell.x[1] = cell.x[1] + rng.uniform(-1,1);
    }
}

void Environment::die(Cell& cell) {
    cell.state = -1;
}

void Environment::testRNGThreadSafety(int numThreads, int numSamples, int seed) {
    std::vector<std::vector<double>> results(numThreads);
    std::vector<std::thread> threads;

    auto rngTest = [&](int threadIndex) {
        RNG rng(seed);  // Each thread gets its own instance
        for (int i = 0; i < numSamples; ++i) {
            results[threadIndex].push_back(rng.uniform(0.0, 1.0));
        }
    };

    for (int i = 0; i < numThreads; ++i) {
        threads.emplace_back(rngTest, i);
    }

    for (auto& thread : threads) {
        thread.join();
    }

    // Verify results
    std::set<double> uniqueValues;
    bool duplicates = false;
    for (const auto& threadResults : results) {
        for (double value : threadResults) {
            if (!uniqueValues.insert(value).second) {
                duplicates = true;
            }
        }
    }

    if (duplicates) {
        std::cout << "❌ RNG thread safety test failed! Some values are duplicated between threads.\n";
    } else {
        std::cout << "✅ RNG thread safety test passed. Each thread has an independent RNG.\n";
    }
}

void Environment::testRNGSeedingAndThreadSafety(int numIterations) {
    int numThreads = omp_get_max_threads();
    std::vector<std::set<double>> threadResults(numThreads);  // Store unique values per thread
    std::mutex mtx;  // Ensure safe output printing

#pragma omp parallel
    {
        int threadId = omp_get_thread_num();
        std::set<double> localValues;

#pragma omp for
        for (int i = 0; i < numIterations; ++i) {
            double val = rng.uniform(0.0, 1.0);
            localValues.insert(val);  // Unique values for this thread
        }

        // Merge results safely
#pragma omp critical
        threadResults[threadId] = std::move(localValues);
    }

    // Check for per-thread uniqueness
    bool passed = true;
    std::set<double> globalValues;

    for (int i = 0; i < numThreads; ++i) {
        for (double value : threadResults[i]) {
            if (!globalValues.insert(value).second) {  // Duplicate found
                passed = false;
            }
        }
    }

    if (passed) {
        std::cout << "✅ RNG seeding and thread safety test passed." << std::endl;
    } else {
        std::cout << "❌ RNG seeding and thread safety test failed! Some values are duplicated." << std::endl;
    }

    // Print results safely
    std::lock_guard<std::mutex> lock(mtx);
    for (int i = 0; i < numThreads; ++i) {
        std::cout << "Thread " << i << ": ";
        for (double value : threadResults[i]) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }
}