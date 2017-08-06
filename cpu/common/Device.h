///
/// Device Header
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef DEVICE_HEADER
#define DEVICE_HEADER

#include "..\..\common\Global.h"
#include "..\..\common\Constants.h"
#include "..\..\common\arvdwfc.h"
#include "..\..\common\getters.h"

#define __fmul_rn(a,b) (a*b)			//cuda redefine

void predict(real4 *r,real3 *v,real3 *a,real3 *b,real3 *c, real dt, int nop);
void correct(real4 *r,real3 *v,real3 *f,real3 *a, real3 *b,real3 *c, real dt, int nop, real nPos, bool & flag, float energyLoss);
void lennardJonesForces(real4 * posArray, real3 * forceArray, int NumberOfParticles, int nlargestsize, int nlistsize, int * list, real rcutsq);
void manyBodyForces(real4 * posArray, real3 * forceArray, int NumberOfParticles, int nlargestsize, int nlistsize, int * list, real rcutsq, bool ljexist);
void calculateAccelerations(real4 * posArray, real3 * accelerationArray, real3 * forceArray, int NumberOfParticles);

#endif