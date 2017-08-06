///
/// CUDA Device Functions
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "..\..\common\Constants.h"

texture<int, 2, cudaReadModeElementType> ljTexRef; //lennrad jones texture reference
texture<int, 2, cudaReadModeElementType> mbTexRef; //many body texture reference

//more functions -> merged together to one cuda module
#include "Getters.cuh"
#include "arvdwfc.cuh"
#include "BondOrderFuncs.cuh"
#include "LennardJones.cuh"
#include "Information.cuh"

///
/// Acceleration Calculation Kernel
///
extern "C" __global__ void calculateAccelerations(float4 * posArray, float3 * forceArray, float3 * aaccArray, int NumberOfParticles)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	if (id>=NumberOfParticles) //more thread than particles
		return;

	//read global memory
	float3 f = forceArray[id];
	float3 a;
	float4 r = posArray[id];

	float mass = getMass((int)r.w); //get mass for particle

	//calculate accelerations
	a.x = f.x/mass;
	a.y = f.y/mass;
	a.z = f.z/mass;

	//write memory
	aaccArray[id] = a;
}

///
/// Kernel to handle the predict function.
///
extern "C" __global__ void predict(float4 * posArray, float3 * velocityArray, float3 * aaccArray, float3 * baccArray, float3 * caccArray, float dt, int nop)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	if (id>=nop) //more thread than particles
		return;

	//read memory
	float4 r = posArray[id];
	float3 v = velocityArray[id];
	float3 a = aaccArray[id];
	float3 b = baccArray[id];
	float3 c = caccArray[id];

	//predict...
	float c1 = dt;
    float c2 = __fmul_rn(c1,dt)/2.0f;
    float c3 = __fmul_rn(c2,dt)/3.0f;
    float c4 = __fmul_rn(c3,dt)/4.0f;

    r.x += __fmul_rn(c1,v.x) + __fmul_rn(c2,a.x) + __fmul_rn(c3,b.x) + __fmul_rn(c4,c.x);
	r.y += __fmul_rn(c1,v.y) + __fmul_rn(c2,a.y) + __fmul_rn(c3,b.y) + __fmul_rn(c4,c.y);
	r.z += __fmul_rn(c1,v.z) + __fmul_rn(c2,a.z) + __fmul_rn(c3,b.z) + __fmul_rn(c4,c.z);
	v.x += __fmul_rn(c1,a.x) + __fmul_rn(c2,b.x) + __fmul_rn(c3,c.x);
	v.y += __fmul_rn(c1,a.y) + __fmul_rn(c2,b.y) + __fmul_rn(c3,c.y);    
	v.z += __fmul_rn(c1,a.z) + __fmul_rn(c2,b.z) + __fmul_rn(c3,c.z);    
	a.x += __fmul_rn(c1,b.x) + __fmul_rn(c2,c.x);
	a.y += __fmul_rn(c1,b.y) + __fmul_rn(c2,c.y); 
	a.z += __fmul_rn(c1,b.z) + __fmul_rn(c2,c.z);
	b.x += __fmul_rn(c1,c.x);
	b.y += __fmul_rn(c1,c.y);
	b.z += __fmul_rn(c1,c.z);

	//write memory
	posArray[id] = r;
	velocityArray[id] = v;
	aaccArray[id] = a;
	baccArray[id] = b;
	caccArray[id] = c;
}

///
/// Kernel to handle the correct function.
///
extern "C" __global__ void correct(float4 * posArray, float3 * velocityArray, float3 * forceArray, float3 * aaccArray, float3 * baccArray, float3 * caccArray, float dt, int nop, float nPos, bool * flag, float energyLoss)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	if (id>=nop) //more thread than particles
		return;

	//read memory
	float4 r = posArray[id];
	float3 v = velocityArray[id];
	float3 f = forceArray[id];
	float3 a = aaccArray[id];
	float3 b = baccArray[id];
	float3 c = caccArray[id];

	float mass = getMass(r.w); //get mass

	float c1 = dt;
    float c2 = __fmul_rn(c1,dt)/2.0f;
    float c3 = __fmul_rn(c2,dt)/3.0f;
    float c4 = __fmul_rn(c3,dt)/4.0f;

	float cr = __fmul_rn(GEAR1,c2);
	float cv = __fmul_rn(GEAR2,c2)/c1;
	float cb = __fmul_rn(GEAR3,c2)/c3;
	float cc = __fmul_rn(GEAR4,c2)/c4;

	float axi = f.x/mass;
  	float ayi = f.y/mass;
  	float azi = f.z/mass;

	float corrx = axi - a.x;
	float corry = ayi - a.y;
	float corrz = azi - a.z;

	r.x += __fmul_rn(cr,corrx);
	r.y += __fmul_rn(cr,corry);
	r.z += __fmul_rn(cr,corrz);
	v.x += __fmul_rn(cv,corrx);
	v.y += __fmul_rn(cv,corry);
	v.z += __fmul_rn(cv,corrz);
	a.x = axi;
	a.y = ayi;
	a.z = azi;
	b.x += __fmul_rn(cb,corrx);
	b.y += __fmul_rn(cb,corry);
	b.z += __fmul_rn(cb,corrz);
	c.x += __fmul_rn(cc,corrx);
	c.y += __fmul_rn(cc,corry);
	c.z += __fmul_rn(cc,corrz);

	//out of the box check
	if (r.x > nPos && r.y > 0)
	{
		*flag = true;
		v.x = -v.x*energyLoss;
		r.x = nPos;
	}

	if (r.x > nPos && r.y < 0)
	{
		*flag = true;
		r.x = -nPos;
	}

	if (r.x < -nPos)
	{
		*flag = true;
		v.x = -v.x;
		r.x = -nPos;
	}

	if (fabs(r.y) > nPos)
	{
		*flag = true;
		v.y = -v.y;
		r.y = __fmul_rn(nPos, (r.y/fabs(r.y)));
	}

	if (fabs(r.z) > nPos)
	{
		*flag = true;
		v.z = -v.z;
		r.z = __fmul_rn(nPos, (r.z/fabs(r.z)));
	}

	//write memory
	posArray[id] = r;
	velocityArray[id] = v;
	aaccArray[id] = a;
	baccArray[id] = b;
	caccArray[id] = c;
}