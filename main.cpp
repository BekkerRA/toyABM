#include <iostream>
#include "Environment.h"
#include "Cell.h"
#include "RNG.h"

int main() {
    int numCells = 10;
    int xDim = 50;
    int yDim = 50;
    int totalTimeSteps = 10;
    int seed = 9999999;

    double divProb = 0.5;
    double deathProb = 0.4;
    double migProb = 0.2;

    //RNG rng(seed);
    Environment model(seed,xDim,yDim,divProb,deathProb,migProb);
    model.initializeCells(numCells);
    model.simulate(totalTimeSteps);
 // model.testRNGThreadSafety(1000,42, 42);
 //    model.testRNGSeedingAndThreadSafety(1000);

    return 0;
}

