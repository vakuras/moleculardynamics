///
/// Neighbor List Header
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef NEIGHBORLIST_HEADER
#define NEIGHBORLIST_HEADER

#include "Global.h"

#define H1				2971215073
#define H2				433494437
#define H3				1073807359

//list settings structure
struct listSettings
{
	float maxnlmove; //lennard jones maximum neighbour list move
	float maxnlmovesq; //lennard jones maximum neighbour list move squared
	float rcut;
	float rs;
	float rcutsq; //lennard jones rcut squared
	int nlistsize; //neightbour list size per particle
	int nlargestsize; //neightbour largest list size
};

int * buildNeighborList(real4 * posArray, listSettings * listsettings, real boxSize, int buckets, int NumberOfParticles);
int computeBucket(int cx, int cy, int cz, int3 * bucketlist, int buckets, bool pop);

#endif