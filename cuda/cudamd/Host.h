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
#include "..\..\common\ConfigParser.h"
#include "..\..\common\Constants.h"
#include "..\..\common\Loader.h"
#include "..\..\common\NeighborList.h"
#include "..\..\common\results.h"

//#define USECUDADEBUG

#ifndef USECUDADEBUG
#include "..\..\common\DebugCalculations.h"
#endif

int hostMain(CUdevice device, char * module_path, configuration * config);
void readyList(real4 * posArray, real3 * velocityArray, int NumberOfParticles, float & highestVelocity, float & boxSize);
void buildList(real4 * posArray, listSettings * listsettings, float boxSize, int buckets, float highestVelocity, int NumberOfParticles, int & nextBuild, int currentTimestep, int ** list, float dt, char * lpTexRef, CUarray * cu_array, CUtexref * cu_texref, CUmodule * module);
#ifdef USECUDADEBUG
void performOutput(vector<results> & resultVec, int BlocksPerGrid, int informationSize, float * informationMemory, CUdeviceptr & devInformation, CUfunction & performCalculations, CUfunction & calculatePotentional, int timesteps, float dt, bool print);
#else
void performOutput(vector<results> * resultVec, int NumberOfParticles, real4 * posArray, real3 * velocityArray, CUdeviceptr & devPosArray, CUdeviceptr & devVelocityArray, int timesteps, real dt, bool print);
#endif

//cude func-build macros
#define CUDA_DEF_SET(v)				offset = (offset + __alignof(v) - 1) & ~(__alignof(v) - 1)
#define OFFSET_ADD(v)				offset += sizeof(v)

#define CUDA_RESET_OFFSET			offset=0
#define CUDA_SET_OFFSET(x)			offset = x
#define CUDA_GET_OFFSET(x)			x = offset
#define CUDA_POINTER_ALLOC(f,p)		ptr = (void*)(size_t)p; offset = (offset + __alignof(ptr) - 1) & ~(__alignof(ptr) - 1); CU_SAFE_CALL(cuParamSetv(f, offset, &ptr, sizeof(ptr));	offset+=sizeof(ptr))
#define CUDA_FLOAT_ALLOC(f,v)		CUDA_DEF_SET(v); CU_SAFE_CALL(cuParamSetf( f, offset, v )); OFFSET_ADD(v)
#define CUDA_UINT_ALLOC(f,v)		CUDA_DEF_SET(v); CU_SAFE_CALL(cuParamSeti( f, offset, v )); OFFSET_ADD(v)

#define THREADSPERBLOCK				512

#endif