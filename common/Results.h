///
/// Results Header
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef RESULTS_HEADER
#define RESULTS_HEADER

#include "Global.h"

//result structure
struct results
{
	double time;
	double temperature;
	double e; //total energy
	double ek; //kinetic energy
	double eu; //potentional energy
	double centerOfMassx; //center of mass
	double centerOfMassy; //center of mass
	double centerOfMassz; //center of mass
	double momentumx; //momentum x
	double momentumy; //momentum y
	double momentumz; //momentum z
};

void pushAnim(vector<float*> & posVec, vector<float*> & velVec, real4* posArray, real3* velocityArray, int nop);
void flushVectors(vector<results> & resultVec, vector<float*> & posVec, vector<float*> & velVec);
void writeAnimationBinaryData(vector<float*> & posVec, vector<float*> & velVec, configuration * config, float vxcm);
void writeResults(vector<results> & resultVec, configuration * config, float elapsed, string tail, float vxcm);
string tailResults(real4 * posArray, int ManyBodyParticles);

#endif