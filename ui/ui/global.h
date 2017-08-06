///
/// Global definitions.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef GLOBALH
#define GLOBALH

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <msclr/marshal.h>

using namespace std;

//config file reader
#define NOP "Lennard Jones Particles"
#define MBP "Many Body Particles"
#define TIMESTEPS "Timesteps"
#define DTQ "DT"
#define TEMP "Temperature"
#define LJRSQ "LennardJonesRS"
#define LJRCUTQ "LennardJonesRCUT"
#define MBRSQ "ManyBodyRS"
#define MBRCUTQ "ManyBodyRCUT"
#define INPUT "Input"
#define OUTPUT "Output"
#define FILENAME "Filename"
#define CUDABLOCKS "CUDA Blocks"
#define USECUDA "Use Cuda"
#define USELJ "UseLennardJones"
#define USEMB "UseManyBody"
#define LJBPP "LennardJonesBPP"
#define MBBPP "ManyBodyBPP"
#define OUTPUTTS "OutputTimesteps"
#define QDEBUG "Debug"
#define FALLBACK "Fallback"
#define VXCMFROM "VXCMFROM"
#define VXCMTO "VXCMTO"
#define VXCMSTEP "VXCMSTEP"
#define ENERGYLOSS "EnergyLoss"
#define ANIMTS "AnimationTimesteps"
#define ANIMDATA "AnimationData"

//safe calls
#define MSGERR(msg) MessageBox::Show(msg, "Error", MessageBoxButtons::OK, MessageBoxIcon::Error)
#define TC(newcode) try {newcode;}catch(::bad_alloc){ MSGERR("Unable to allocate memory"); return (false);}
#define CPPSAFE_CALL(code,err) if(code) { MSGERR(err); return (false);}

struct configuration
{
	int LennardJonesParticles; //number of particles in the simulation
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
	string Filename; //result filename template
	string AnimFile; //filename for animation data
};

namespace ui
{
	//
	// Native functions (NOT .NET)
	//
	string MakeCSTR(System::String ^ ); //.NET String to C++ String
	bool ReadConfiguration(string filename, configuration * pconfig); //Read configuration from ini file
	bool WriteConfiguration(string filename, configuration * pconfig); //Write configuration to ini file
	bool ReadInput(const char * filename, int nop, float * posArray, float * velocityArray);	//Read input file
};

#endif