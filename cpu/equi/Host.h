///
/// Host Header
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef HOST_HEADER
#define HOST_HEADER

#include "..\..\common\Global.h"
#include "..\..\common\NeighborList.h"
#include "..\common\Device.h"
#include "..\..\common\ConfigParser.h"
#include "..\..\common\Constants.h"
#include "..\..\common\DebugCalculations.h"
#include "Statistics.h"
#include "..\..\common\Loader.h"
#include "..\..\common\results.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

#define MD_RANDMAX (RAND_MAX + 2)

int hostMain(configuration * config);
void readyList(real4 * posArray, real3 * velocityArray, int NumberOfParticles, real & highestVelocity, real & boxSize);
void buildList(real4 * posArray, listSettings * listsettings, real boxSize, int buckets, real highestVelocity, int NumberOfParticles, int & nextBuild, int currentTimestep, int ** list, real dt);
bool choosePositionsNoble(real4 * posArray, real boxl, configuration * config);
bool choosePositionsManyBody(real4 * posArray, real boxl, configuration * config);
void chooseVelocities(real4 * posArray, real3 * velocityArray, real temp, configuration * config);

#endif