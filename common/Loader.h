///
/// Loader Header
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef LOADER_HEADER
#define LOADER_HEADER

#include "Global.h"
#include "ConfigParser.h"

#ifdef UI
using namespace System::Windows::Forms;
#endif

void readConfiguration(string filename, configuration * pconfig);
void writeOutput(configuration * config, real * posArray, real * velocityArray);
void readInput(configuration * config, real * posArray, real * velocityArray);
void Loader(int argc, char ** argv, char * appexe, configuration & config);

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

#endif