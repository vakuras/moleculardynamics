///
/// Neighbor List Header
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef NEIGHBORLIST_HEADER
#define NEIGHBORLIST_HEADER

#define H1				2971215073
#define H2				433494437
#define H3				1073807359

int * buildNeighborList(float4 * posArray, listSettings * listsettings, float boxSize, int buckets, int NumberOfParticles);
int computeBucket(int cx, int cy, int cz, int3 * bucketlist, int buckets, bool pop);

#endif