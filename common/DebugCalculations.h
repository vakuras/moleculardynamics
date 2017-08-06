///
/// Output Header
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef DEBUGCALCULATIONS_HEADER
#define DEBUGCALCULATIONS_HEADER

#include "results.h"
#include "Global.h"
#include "Constants.h"
#include "arvdwfc.h"
#include "getters.h"

void PerformCalculations(real4 * posArray, real3 * velocityArray, int NumberOfParticles, double & centerOfMassx, double & centerOfMassy, double & centerOfMassz, double & kineticEnergy, double & momentumx, double & momentumy, double & momentumz, double & temperature);
void performOutput(vector<results> * resultVec, int NumberOfParticles, real4 * posArray, real3 * velocityArray, int timesteps, real dt, bool print);
void calculatePotentional(real4 * posArray, double * vibLJ, double * vibMB, int NumberOfParticles);

#endif