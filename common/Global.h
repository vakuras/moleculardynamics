///
/// Global Header
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef GLOBAL_HEADER
#define GLOBAL_HEADER

#include <cmath>
#include <iostream>
#include <string>
#include <windows.h>
#include <fstream>
#include <vector>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <sstream>

using namespace std;

#ifdef EQUI //euqilibrate
typedef double real;
#elif CPUMD //cpu simulation
typedef float real;
#elif CUDAMD //cuda simulation
typedef float real;
#elif UI //ui
typedef double real;
#endif

#ifdef CUDAMD
#include <cutil.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

//add libraries from cuda
#pragma comment(lib, "cuda.lib")
#pragma comment(lib, "cudart.lib")

#ifdef _DEBUG
#pragma comment(lib, "cutil64d.lib")
#else
#pragma comment(lib, "cutil64.lib")
#endif
#else
struct int3
{
	int x;
	int y;
	int z;
};
#endif

struct real4
{
	real x;
	real y;
	real z;
	real w;
};


struct real3
{
	real x;
	real y;
	real z;
};

//configuration structure
struct configuration
{
	int LennardJonesParticles; //number of lennard jones particles in the simulation
	int ManyBodyParticles;	//number of many body particles
	float DT; //delta-t
	int Timesteps; //timesteps to make
	float Temperature; //temperature
	float LennardJonesRS; //LennardJones rs
	float LennardJonesRCUT; //LennardJones rcut
	float ManyBodyRS; //many body rs
	float ManyBodyRCUT; //many body rcut
	int CudaBlocks; //amount of blocks for cuda
	bool UseCuda; // use cuda - ui only
	bool useLennardJones; //to use the lennard-jones force
	bool useManyBody; //to use the many-body force
	bool lennardJonesBPP; //use the block per particle version of lennard jones forces
	bool manyBodyBPP; //use the block per particle version of many body forces
	int OutputTimesteps; //when to perform output
	bool Debug; //to output during simulation
	bool Fallback; //to change back to tpp version if bpp is inperformable
	float vxcmFrom; //velocity to add - x direction (from)
	float vxcmTo; //velocity to add - x direction (to)
	float vxcmStep; //velocity to add - x direction (step)
	float energyLoss; //energyLoss for box hits
	int animts; //each ts to save animation data
	string Input; //input file
	string Output; //output file
	string Filename; //filename template for results
	string AnimFile; //filename for animation data
};

//safe calls
#ifdef UI
#include <msclr/marshal.h>

#ifndef CFGPARSER_H
#define MSGERR(msg) MessageBox::Show(msg, "Error", MessageBoxButtons::OK, MessageBoxIcon::Error)
#define TC(newcode) try {newcode;}catch(::bad_alloc){ MSGERR("Unable to allocate memory"); return (false);}
#define CPPSAFE_CALL(code,err) if(code) { MSGERR(err); return (false);}
#endif
#else
#define TC(newcode) try {newcode;}catch(bad_alloc){cerr << "Error: Unable to allocate memory.\nFile: '" << __FILE__ << "'\n:Line: " << __LINE__ << endl; exit(EXIT_FAILURE);}
#define CPPSAFE_CALL(code,err) if(code) {cerr << "Error: " << err << "\nFile: '" << __FILE__ << "'\nLine: " << __LINE__ << endl; exit(EXIT_FAILURE);}
#endif

#endif