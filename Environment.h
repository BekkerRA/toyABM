//
// Created by Rebecca Bekker on 3/26/25.
//

#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H
#include <random>
#include <vector>
#include <algorithm>
#include <string>
#include <filesystem>
#include "Cell.h"
#include "RNG.h"
#include "/opt/homebrew/opt/libomp/include/omp.h"

class Environment {
public:
    Environment(int seed, int xDim, int yDim, double divProb, double deathProb, double migProb);
    void initializeCells(int numCells);
    void simulate(int finalTimeStep);
    void cleanCells();
    void shuffleCells();
    void update_parallel();
    void update_sequential();
    void divide(Cell& cell);
    void move(Cell& cell);
    void die(Cell& cell);
    void testRNGThreadSafety(int numThreads, int numSamples, int seed);
    void testRNGSeedingAndThreadSafety(int numIterations);

private:
    std::vector<Cell> cell_list;
    int xDim;
    int yDim;
    RNG rng;
    double divProb;
    double deathProb;
    double moveProb;

};
#endif //ENVIRONMENT_H
