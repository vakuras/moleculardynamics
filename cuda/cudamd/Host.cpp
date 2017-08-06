///
/// Host Implementation
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "Host.h"

///
/// Simulation main host function
///
int hostMain(CUdevice device, char * module_path, configuration * config)
{
	real4 * posArray; //positions (host)
	real3 * velocityArray; //velocity (host)
	bool * flag; //box hit flag (host)
	int * ljList = NULL; //lennard jones neigbor list pointer
	int * mbList = NULL; //many body neigbor list pointer
	listSettings lj; //lennard jones neigbor list settings
	listSettings mb; //many body neigbor list settings
	float highestVelocity; //heighest velocity
	int timesteps; //timestep counter
	int ljNextBuild; //on which timestep to call the build of the lennard jones list
	int mbNextBuild; //on which timestep to call the build of the many body list
	float boxSize=pow((float)(config->LennardJonesParticles + config->ManyBodyParticles),1.0f/3.0f)*52.8f; //box size
	float boxSizeList; //box size for the neighbor lists
	int buckets; //bucket count for the neighbor list
	bool fallbackmb = false; //fallback from tpp to bpp for many body
	bool fallbacklj = false; //fallback from tpp to bpp for lennard jones

	vector<results> resultVec; //results vector
	vector<float*> posVec; //position vector for animation
	vector<float*> velVec; //velocity vector for animation

#ifdef USECUDADEBUG
	float * informationMemory; //information memory pointer
	CUdeviceptr devInformation; //information memory pointer for cuda

	//vibLJ, vibMB, kinetic, temperature, centerOfMass, momentum
	const int informationSize = sizeof(float)+sizeof(float)+sizeof(float)+sizeof(float)+sizeof(float3)+sizeof(real3);
#endif

	//device variables
	CUdeviceptr devPosArray; //positions cuda pointer
	CUdeviceptr devVelocityArray; //velocity cuda pointer
	CUdeviceptr devAAccArray; //acceleration cuda pointer
	CUdeviceptr devForceArray; //force cuda pointer
	CUdeviceptr devBAccArray; //b-acc cuda pointer
	CUdeviceptr devCAccArray; //c-acc cuda pointer
	CUdeviceptr devMemAlloc; //device block-shared memory pointer
	CUdeviceptr devFlag; //device box-hit-flag
	CUcontext context; //device context
	CUmodule module; //cude module (.ptx file)
	//function pointers from module:
#ifdef USECUDADEBUG
	CUfunction performCalculations;
	CUfunction calculatePotentional;
#endif
	CUfunction correct;
	CUfunction predict;
	CUfunction lennardJonesForces;
	CUfunction lennardJonesForcesBPP;
	CUfunction manyBodyForces1;
	CUfunction manyBodyForces2;
	CUfunction manyBodyForcesBPP1;
	CUfunction manyBodyForcesBPP2;
	CUfunction calculateAccelerations;

	//textures for neighbor lists
	CUarray devLjList = NULL;
	CUarray devMbList = NULL;
	CUtexref devLjTexRef = NULL;
	CUtexref devMbTexRef = NULL;

	//function build helpers - for use by consts defined in Host.h (CUDA_DEF_SET, CUDA_RESET_OFFSET... and so)
	void * ptr;
	int offset;
	unsigned int val;
	float fval;

	int nlsoffset; //neighbor list 'largest list size' byte offset in the byte array of the force calculation functions

	//block configuration
	int BlocksPerGrid = config->CudaBlocks;
	int ThreadsPerBlocks = (config->LennardJonesParticles + config->ManyBodyParticles)/BlocksPerGrid + ((config->LennardJonesParticles + config->ManyBodyParticles)%BlocksPerGrid == 0 ? 0:1);
	int mbThreadsPerBlocks = config->ManyBodyParticles/BlocksPerGrid + (config->ManyBodyParticles%BlocksPerGrid == 0 ? 0:1);

    //time-measure events
	CUevent start;
	CUevent stop;
	float elapsedTime;

	config->energyLoss = sqrt(config->energyLoss); //fix energy loss

	//create context
	CU_SAFE_CALL(cuCtxCreate(&context, 0, device));

	//load module
	CU_SAFE_CALL(cuModuleLoad(&module, module_path));

	//get functions
#ifdef USECUDADEBUG
	CU_SAFE_CALL(cuModuleGetFunction(&performCalculations, module, "performCalculations"));
	CU_SAFE_CALL(cuModuleGetFunction(&calculatePotentional, module, "calculatePotentional"));
#endif
	CU_SAFE_CALL(cuModuleGetFunction(&predict, module, "predict"));
	CU_SAFE_CALL(cuModuleGetFunction(&correct, module, "correct"));
	CU_SAFE_CALL(cuModuleGetFunction(&calculateAccelerations, module, "calculateAccelerations"));

	if (config->useLennardJones) //if using lennard jones
	{
		CU_SAFE_CALL(cuModuleGetFunction(&lennardJonesForcesBPP, module, "lennardJonesForcesBPP"));
		CU_SAFE_CALL(cuModuleGetFunction(&lennardJonesForces, module, "lennardJonesForces"));
	}

	if (config->useManyBody) //if using many body
	{
		CU_SAFE_CALL(cuModuleGetFunction(&manyBodyForcesBPP1, module, "manyBodyForcesBPP1"));
		CU_SAFE_CALL(cuModuleGetFunction(&manyBodyForcesBPP2, module, "manyBodyForcesBPP2"));
		CU_SAFE_CALL(cuModuleGetFunction(&manyBodyForces1, module, "manyBodyForces1"));
		CU_SAFE_CALL(cuModuleGetFunction(&manyBodyForces2, module, "manyBodyForces2"));
		CU_SAFE_CALL(cuMemAlloc(&devMemAlloc,sizeof(float2) * config->ManyBodyParticles * config->ManyBodyParticles)); //allocate memory used by many body
	}
	
	//set events
	CU_SAFE_CALL(cuEventCreate(&start, CU_EVENT_DEFAULT));
	CU_SAFE_CALL(cuEventCreate(&stop, CU_EVENT_DEFAULT));

	//allocate memory for data on host (for future device mapping)
	CU_SAFE_CALL(cuMemAllocHost((void**)&posArray, (config->LennardJonesParticles + config->ManyBodyParticles) * sizeof(real4)));
	CU_SAFE_CALL(cuMemAllocHost((void**)&velocityArray, (config->LennardJonesParticles + config->ManyBodyParticles) * sizeof(real3)));
#ifdef USECUDADEBUG
	CU_SAFE_CALL(cuMemAllocHost((void**)&informationMemory, informationSize));
#endif
	CU_SAFE_CALL(cuMemAllocHost((void**)&flag, sizeof(bool)));
	
	//allocate memory for data on device
	CU_SAFE_CALL(cuMemAlloc(&devPosArray,sizeof(real4) * (config->LennardJonesParticles + config->ManyBodyParticles)));
	CU_SAFE_CALL(cuMemAlloc(&devVelocityArray,sizeof(real3) * (config->LennardJonesParticles + config->ManyBodyParticles)));
	CU_SAFE_CALL(cuMemAlloc(&devAAccArray,sizeof(real3) * (config->LennardJonesParticles + config->ManyBodyParticles)));
	CU_SAFE_CALL(cuMemAlloc(&devForceArray,sizeof(real3) * (config->LennardJonesParticles + config->ManyBodyParticles)));
	CU_SAFE_CALL(cuMemAlloc(&devBAccArray,sizeof(real3) * (config->LennardJonesParticles + config->ManyBodyParticles)));
	CU_SAFE_CALL(cuMemAlloc(&devCAccArray,sizeof(real3) * (config->LennardJonesParticles + config->ManyBodyParticles)));
#ifdef USECUDADEBUG
	CU_SAFE_CALL(cuMemAlloc(&devInformation,informationSize));
#endif
	CU_SAFE_CALL(cuMemAlloc(&devFlag,sizeof(bool)));

	//number of buckets for list hashs
	buckets = (config->LennardJonesParticles + config->ManyBodyParticles) * 3;

	mb.nlargestsize = lj.nlargestsize = 0;

	if (config->useLennardJones) //set lennard jones list settings
	{
		lj.maxnlmove = ((config->LennardJonesRS) - (config->LennardJonesRCUT)) * 0.5f;
		lj.maxnlmovesq = pow(lj.maxnlmove, 2);
		lj.rcutsq = pow(config->LennardJonesRCUT,2);
		lj.rcut = config->LennardJonesRCUT;
		lj.rs = config->LennardJonesRS;
	}

	if (config->useManyBody) //set many body list settings
	{
		mb.maxnlmove = ((config->ManyBodyRS) - (config->ManyBodyRCUT)) * 0.5f;
		mb.maxnlmovesq = pow(mb.maxnlmove, 2);
		mb.rcutsq = pow(config->ManyBodyRCUT,2);
		mb.rcut = config->ManyBodyRCUT;
		mb.rs = config->ManyBodyRS;
	}

#ifdef USECUDADEBUG
	//parameter setup for performCalculations
	CUDA_RESET_OFFSET; //reset offset to zero
	CUDA_POINTER_ALLOC(performCalculations, devPosArray); //allocate pointer in the byte array of the function
	CUDA_POINTER_ALLOC(performCalculations, devVelocityArray);
	val = (config->LennardJonesParticles + config->ManyBodyParticles);
	CUDA_UINT_ALLOC(performCalculations, val); //allocate unsigned int in the byte array of the function
	CUDA_POINTER_ALLOC(performCalculations, devInformation);
	CU_SAFE_CALL(cuParamSetSize( performCalculations, offset )); //set byte array size
    CU_SAFE_CALL(cuFuncSetBlockShape( performCalculations, 1, 1, 1 )); //set block grid

	//parameter setup for calculatePotentional
	CUDA_RESET_OFFSET;
	CUDA_POINTER_ALLOC(calculatePotentional, devPosArray);
	val = (config->LennardJonesParticles + config->ManyBodyParticles);
	CUDA_UINT_ALLOC(calculatePotentional, val);
	CUDA_POINTER_ALLOC(calculatePotentional, devInformation);
	CU_SAFE_CALL(cuParamSetSize( calculatePotentional, offset ));
    CU_SAFE_CALL(cuFuncSetBlockShape( calculatePotentional, ThreadsPerBlocks, 1, 1 ));
#endif

	//parameter setup for predict
	CUDA_RESET_OFFSET;
	CUDA_POINTER_ALLOC(predict, devPosArray);
	CUDA_POINTER_ALLOC(predict, devVelocityArray);
	CUDA_POINTER_ALLOC(predict, devAAccArray);
	CUDA_POINTER_ALLOC(predict, devBAccArray);
	CUDA_POINTER_ALLOC(predict, devCAccArray);
	fval = config->DT;
	CUDA_FLOAT_ALLOC(predict, fval); //allocate float in the byte array of the function
	val = (config->LennardJonesParticles + config->ManyBodyParticles);
	CUDA_UINT_ALLOC(predict, val);
	CU_SAFE_CALL(cuParamSetSize( predict, offset ));
    CU_SAFE_CALL(cuFuncSetBlockShape( predict, ThreadsPerBlocks, 1, 1 ));

	//parameter setup for correct
	CUDA_RESET_OFFSET;
	CUDA_POINTER_ALLOC(correct, devPosArray);
	CUDA_POINTER_ALLOC(correct, devVelocityArray);
	CUDA_POINTER_ALLOC(correct, devForceArray);
	CUDA_POINTER_ALLOC(correct, devAAccArray);
	CUDA_POINTER_ALLOC(correct, devBAccArray);
	CUDA_POINTER_ALLOC(correct, devCAccArray);
	fval = config->DT;
	CUDA_FLOAT_ALLOC(correct, fval);
	val = (config->LennardJonesParticles + config->ManyBodyParticles);
	CUDA_UINT_ALLOC(correct, val);
	fval = boxSize;
	CUDA_FLOAT_ALLOC(correct, fval);
	CUDA_POINTER_ALLOC(correct, devFlag);
	fval = config->energyLoss;
	CUDA_FLOAT_ALLOC(correct, fval);
	CU_SAFE_CALL(cuParamSetSize( correct, offset ));
    CU_SAFE_CALL(cuFuncSetBlockShape( correct, ThreadsPerBlocks, 1, 1 ));

	//parameter setup for calculateAccelerations
	CUDA_RESET_OFFSET;
	CUDA_POINTER_ALLOC(calculateAccelerations, devPosArray);
	CUDA_POINTER_ALLOC(calculateAccelerations, devForceArray);
	CUDA_POINTER_ALLOC(calculateAccelerations, devAAccArray);
	val = (config->LennardJonesParticles + config->ManyBodyParticles);
	CUDA_UINT_ALLOC(calculateAccelerations, val);
	CU_SAFE_CALL(cuParamSetSize( calculateAccelerations, offset ));
    CU_SAFE_CALL(cuFuncSetBlockShape( calculateAccelerations, ThreadsPerBlocks, 1, 1 ));

	if (config->useLennardJones)
	{
		//parameter setup for lennardJonesForcesBPP
		CUDA_RESET_OFFSET;
		CUDA_POINTER_ALLOC(lennardJonesForcesBPP, devPosArray);
		CUDA_POINTER_ALLOC(lennardJonesForcesBPP, devForceArray);
		fval = lj.rcutsq;
		CUDA_FLOAT_ALLOC(lennardJonesForcesBPP, fval);
		CU_SAFE_CALL(cuParamSetSize( lennardJonesForcesBPP, offset ));
		
		//parameter setup for lennardJonesForces
		CUDA_RESET_OFFSET;
		CUDA_POINTER_ALLOC(lennardJonesForces, devPosArray);
		CUDA_POINTER_ALLOC(lennardJonesForces, devForceArray);
		val = 0; //temporarily
		CUDA_GET_OFFSET(nlsoffset);
		CUDA_UINT_ALLOC(lennardJonesForces, val);
		val = (config->LennardJonesParticles + config->ManyBodyParticles);
		CUDA_UINT_ALLOC(lennardJonesForces, val);
		fval = lj.rcutsq;
		CUDA_FLOAT_ALLOC(lennardJonesForces, fval);
		CU_SAFE_CALL(cuParamSetSize( lennardJonesForces, offset ));
		CU_SAFE_CALL(cuFuncSetBlockShape( lennardJonesForces, ThreadsPerBlocks, 1, 1 ));
	}

	if (config->useManyBody)
	{
		//parameter setup for manyBodyForcesBPP1
		CUDA_RESET_OFFSET;
		CUDA_POINTER_ALLOC(manyBodyForcesBPP1, devPosArray);
		CUDA_POINTER_ALLOC(manyBodyForcesBPP1, devForceArray);
		fval = mb.rcutsq;
		CUDA_FLOAT_ALLOC(manyBodyForcesBPP1, fval);
		val = config->ManyBodyParticles;
		CUDA_UINT_ALLOC(manyBodyForcesBPP1, val);
		CUDA_POINTER_ALLOC(manyBodyForcesBPP1, devMemAlloc);
		val = config->useLennardJones;
		CUDA_UINT_ALLOC(manyBodyForcesBPP1, val);
		CU_SAFE_CALL(cuParamSetSize( manyBodyForcesBPP1, offset ));

		//parameter setup for manyBodyForcesBPP2
		CUDA_RESET_OFFSET;
		CUDA_POINTER_ALLOC(manyBodyForcesBPP2, devPosArray);
		CUDA_POINTER_ALLOC(manyBodyForcesBPP2, devForceArray);
		fval = mb.rcutsq;
		CUDA_FLOAT_ALLOC(manyBodyForcesBPP2, fval);
		val = config->ManyBodyParticles;
		CUDA_UINT_ALLOC(manyBodyForcesBPP2, val);
		CUDA_POINTER_ALLOC(manyBodyForcesBPP2, devMemAlloc);
		CU_SAFE_CALL(cuParamSetSize( manyBodyForcesBPP2, offset ));

		//parameter setup for manyBodyForces1
		CUDA_RESET_OFFSET;
		CUDA_POINTER_ALLOC(manyBodyForces1, devPosArray);
		CUDA_POINTER_ALLOC(manyBodyForces1, devForceArray);
		val = 0; //temporarily
		CUDA_GET_OFFSET(nlsoffset);
		CUDA_UINT_ALLOC(manyBodyForces1, val);
		val = config->ManyBodyParticles;
		CUDA_UINT_ALLOC(manyBodyForces1, val);
		fval = mb.rcutsq;
		CUDA_FLOAT_ALLOC(manyBodyForces1, fval);
		CUDA_POINTER_ALLOC(manyBodyForces1, devMemAlloc);
		val = config->useLennardJones;
		CUDA_UINT_ALLOC(manyBodyForces1, val);
		CU_SAFE_CALL(cuParamSetSize( manyBodyForces1, offset ));
		CU_SAFE_CALL(cuFuncSetBlockShape( manyBodyForces1, mbThreadsPerBlocks, 1, 1 ));

		//parameter setup for manyBodyForces2
		CUDA_RESET_OFFSET;
		CUDA_POINTER_ALLOC(manyBodyForces2, devPosArray);
		CUDA_POINTER_ALLOC(manyBodyForces2, devForceArray);
		val = 0; //temporarily
		CUDA_UINT_ALLOC(manyBodyForces2, val);
		val = config->ManyBodyParticles;
		CUDA_UINT_ALLOC(manyBodyForces2, val);
		fval = mb.rcutsq;
		CUDA_FLOAT_ALLOC(manyBodyForces2, fval);
		CUDA_POINTER_ALLOC(manyBodyForces2, devMemAlloc);
		CU_SAFE_CALL(cuParamSetSize( manyBodyForces2, offset ));
		CU_SAFE_CALL(cuFuncSetBlockShape( manyBodyForces2, mbThreadsPerBlocks, 1, 1 ));
	}

	for(float vxcm=config->vxcmFrom; vxcm<=config->vxcmTo; vxcm+=config->vxcmStep)
	{
		//init
		timesteps = 0;
		ljNextBuild = mbNextBuild = 0;

		//information
		int printout = config->OutputTimesteps; //next debug output time step
		int animts = config->animts; //next animation save time step

		//read input file
		readInput(config, (float*) posArray, (float*) velocityArray);

		//push the particles outside the box into the box & vxcm
		for (int id=0; id< (config->LennardJonesParticles + config->ManyBodyParticles); id++)
		{
			if (fabs(posArray[id].x) > (boxSize/2))
				posArray[id].x = (boxSize/2) * (posArray[id].x/fabs(posArray[id].x));

			if (fabs(posArray[id].y) > (boxSize/2))
				posArray[id].y = (boxSize/2) * (posArray[id].y/fabs(posArray[id].y));

			if (fabs(posArray[id].z) > (boxSize/2))
				posArray[id].z = (boxSize/2) * (posArray[id].z/fabs(posArray[id].z));

			velocityArray[id].x += vxcm;
		}

		//copy position & velocity to device memory 
		CU_SAFE_CALL(cuMemcpyHtoD(devPosArray, posArray, sizeof(real4) * (config->LennardJonesParticles + config->ManyBodyParticles)));
		CU_SAFE_CALL(cuMemcpyHtoD(devVelocityArray, velocityArray, sizeof(real3) * (config->LennardJonesParticles + config->ManyBodyParticles)));

		//zero some device memory
		CU_SAFE_CALL(cuMemsetD32(devAAccArray, 0, 3 * (config->LennardJonesParticles + config->ManyBodyParticles)));
		CU_SAFE_CALL(cuMemsetD32(devBAccArray, 0, 3 * (config->LennardJonesParticles + config->ManyBodyParticles)));
		CU_SAFE_CALL(cuMemsetD32(devCAccArray, 0, 3 * (config->LennardJonesParticles + config->ManyBodyParticles)));
		CU_SAFE_CALL(cuMemsetD8(devFlag, 0, sizeof(bool)));

		//before sim output
		cout << fixed << setprecision(6); //set maximal precision for output
		cout.setf(ios::fixed,ios::floatfield); //zero padding
		cout << "Before Simulation [VXCM = " << vxcm << "]:" << endl;

		//calculate highest velocity and boxSizeList
		readyList(posArray, velocityArray, (config->LennardJonesParticles + config->ManyBodyParticles), highestVelocity, boxSizeList);

		if (config->useLennardJones)
		{
			//build lennard jones list
			buildList(posArray, &lj, boxSizeList, buckets, highestVelocity, (config->LennardJonesParticles + config->ManyBodyParticles), ljNextBuild, 0, &ljList, config->DT, "ljTexRef", &devLjList, &devLjTexRef, &module);
			//check if a fallback is to be used
			fallbacklj = (lj.nlargestsize > THREADSPERBLOCK) && config->Fallback;

			//set the 'largest list size' according to the new data
			if (config->lennardJonesBPP && !fallbacklj)
			{
				if (lj.nlargestsize > 0)
				{
					CU_SAFE_CALL(cuFuncSetBlockShape( lennardJonesForcesBPP, lj.nlargestsize, 1, 1 ));
					CU_SAFE_CALL(cuFuncSetSharedSize( lennardJonesForcesBPP, lj.nlargestsize*sizeof(real4)));
				}
			}
			else
			{
				val = lj.nlargestsize;
				CUDA_SET_OFFSET(nlsoffset);
				CUDA_UINT_ALLOC(lennardJonesForces, val);
			}
		}
		if (config->useManyBody) 
		{
			//build many body list
			buildList(posArray, &mb, boxSizeList, buckets, highestVelocity, config->ManyBodyParticles, mbNextBuild, 0, &mbList, config->DT, "mbTexRef", &devMbList, &devMbTexRef, &module);
			//check if a fallback is to be used
			fallbackmb = (mb.nlargestsize > THREADSPERBLOCK) && config->Fallback;

			//set the 'largest list size' according to the new data
			if (config->manyBodyBPP && !fallbackmb)
			{
				if (mb.nlargestsize>0)
				{
					CU_SAFE_CALL(cuFuncSetBlockShape( manyBodyForcesBPP1, mb.nlargestsize, 1, 1 ));
					CU_SAFE_CALL(cuFuncSetSharedSize( manyBodyForcesBPP1, mb.nlargestsize*sizeof(real4)));
					CU_SAFE_CALL(cuFuncSetBlockShape( manyBodyForcesBPP2, mb.nlargestsize, 1, 1 ));
					CU_SAFE_CALL(cuFuncSetSharedSize( manyBodyForcesBPP2, mb.nlargestsize*sizeof(real4)));
				}
			}
			else
			{
				val = mb.nlargestsize;
				CUDA_SET_OFFSET(nlsoffset);
				CUDA_UINT_ALLOC(manyBodyForces1, val);
				CUDA_SET_OFFSET(nlsoffset);
				CUDA_UINT_ALLOC(manyBodyForces2, val);
			}
		}

		//perform the output function
#ifdef USECUDADEBUG
		performOutput(resultVec, BlocksPerGrid, informationSize, informationMemory, devInformation, performCalculations, calculatePotentional, timesteps, config->DT, true);
#else
		performOutput(&resultVec, (config->LennardJonesParticles + config->ManyBodyParticles), posArray, velocityArray, devPosArray, devVelocityArray, timesteps, config->DT, true);
#endif
		cout << endl;

		if (config->animts>-1)
			pushAnim(posVec, velVec, posArray, velocityArray, (config->LennardJonesParticles + config->ManyBodyParticles));

		//measure & perform
		CU_SAFE_CALL(cuEventRecord(start,0));

		if (config->useLennardJones && lj.nlargestsize>0) //calculate force for lennard jones
		{
			if (config->lennardJonesBPP && !fallbacklj) //bpp and no fallback
			{
				CU_SAFE_CALL(cuLaunchGrid(lennardJonesForcesBPP, (config->LennardJonesParticles + config->ManyBodyParticles), 1));
				CU_SAFE_CALL(cuCtxSynchronize()); //synchronize context
			}
			else //tpp
			{
				CU_SAFE_CALL(cuLaunchGrid(lennardJonesForces, BlocksPerGrid, 1));
				CU_SAFE_CALL(cuCtxSynchronize()); //synchronize context
			}
		}

		if (config->useManyBody && mb.nlargestsize>0) //calculate force for many body
		{
			if (config->manyBodyBPP && !fallbackmb) //bpp and no fallback
			{
				CU_SAFE_CALL(cuLaunchGrid(manyBodyForcesBPP1, config->ManyBodyParticles, 1)); //b-o 1
				CU_SAFE_CALL(cuCtxSynchronize());
				CU_SAFE_CALL(cuLaunchGrid(manyBodyForcesBPP2, config->ManyBodyParticles, 1)); //b-o 2
				CU_SAFE_CALL(cuCtxSynchronize());
			}
			else //tpp
			{
				CU_SAFE_CALL(cuLaunchGrid(manyBodyForces1, BlocksPerGrid, 1)); //b-o 1
				CU_SAFE_CALL(cuCtxSynchronize());
				CU_SAFE_CALL(cuLaunchGrid(manyBodyForces2, BlocksPerGrid, 1)); //b-o 2
				CU_SAFE_CALL(cuCtxSynchronize());
			}
		}

		//calculate the accelerations after the force
		CU_SAFE_CALL(cuLaunchGrid(calculateAccelerations, BlocksPerGrid, 1));
		CU_SAFE_CALL(cuCtxSynchronize());

		//while not reached the timestep goal
		while(timesteps<config->Timesteps)
		{
			if (ljNextBuild == timesteps || mbNextBuild == timesteps) //is is time to build any of the lists?
			{
				//copy memory from device
				CU_SAFE_CALL(cuMemcpyDtoH(posArray, devPosArray, sizeof(real4) * (config->LennardJonesParticles + config->ManyBodyParticles)));
				CU_SAFE_CALL(cuMemcpyDtoH(velocityArray, devVelocityArray, sizeof(real3) * (config->LennardJonesParticles + config->ManyBodyParticles)));

				readyList(posArray, velocityArray, (config->LennardJonesParticles + config->ManyBodyParticles), highestVelocity, boxSizeList);

				if (config->useLennardJones && ljNextBuild == timesteps) //is it time to build the lj list?
				{
					buildList(posArray, &lj, boxSizeList, buckets, highestVelocity, (config->LennardJonesParticles + config->ManyBodyParticles), ljNextBuild, timesteps, &ljList, config->DT, "ljTexRef", &devLjList, &devLjTexRef, &module);

					fallbacklj = (lj.nlargestsize > THREADSPERBLOCK) && config->Fallback; //check for fallback

					//set the 'largest list size' according to the new data
					if (config->lennardJonesBPP && !fallbacklj)
					{
						if (lj.nlargestsize > 0)
						{
							CU_SAFE_CALL(cuFuncSetBlockShape( lennardJonesForcesBPP, lj.nlargestsize, 1, 1 ));
							CU_SAFE_CALL(cuFuncSetSharedSize( lennardJonesForcesBPP, lj.nlargestsize*sizeof(real4)));
						}
					}
					else
					{
						val = lj.nlargestsize;
						CUDA_SET_OFFSET(nlsoffset);
						CUDA_UINT_ALLOC(lennardJonesForces, val);
					}
				}

				if (config->useManyBody && mbNextBuild == timesteps) // is it time to build the mb list?
				{
					buildList(posArray, &mb, boxSizeList, buckets, highestVelocity, config->ManyBodyParticles, mbNextBuild, timesteps, &mbList, config->DT, "mbTexRef", &devMbList, &devMbTexRef, &module);		

					fallbackmb = (mb.nlargestsize > THREADSPERBLOCK) && config->Fallback; //check for fallback

					//set the 'largest list size' according to the new data
					if (config->manyBodyBPP  && !fallbackmb)
					{
						if (mb.nlargestsize>0)
						{
							CU_SAFE_CALL(cuFuncSetBlockShape( manyBodyForcesBPP1, mb.nlargestsize, 1, 1 ));
							CU_SAFE_CALL(cuFuncSetSharedSize( manyBodyForcesBPP1, mb.nlargestsize*sizeof(real4)));
							CU_SAFE_CALL(cuFuncSetBlockShape( manyBodyForcesBPP2, mb.nlargestsize, 1, 1 ));
							CU_SAFE_CALL(cuFuncSetSharedSize( manyBodyForcesBPP2, mb.nlargestsize*sizeof(real4)));
						}
					}
					else
					{
						val = mb.nlargestsize;
						CUDA_SET_OFFSET(nlsoffset);
						CUDA_UINT_ALLOC(manyBodyForces1, val);
						CUDA_SET_OFFSET(nlsoffset);
						CUDA_UINT_ALLOC(manyBodyForces2, val);
					}
				}
			}

			//predict
			CU_SAFE_CALL(cuLaunchGrid(predict, BlocksPerGrid,1));
			CU_SAFE_CALL(cuCtxSynchronize());

			//force lj
			if (config->useLennardJones && lj.nlargestsize>0)
			{
				if (config->lennardJonesBPP && !fallbacklj)
				{
					CU_SAFE_CALL(cuLaunchGrid(lennardJonesForcesBPP, (config->LennardJonesParticles + config->ManyBodyParticles),1));
					CU_SAFE_CALL(cuCtxSynchronize());
				}
				else
				{
					CU_SAFE_CALL(cuLaunchGrid(lennardJonesForces, BlocksPerGrid,1));
					CU_SAFE_CALL(cuCtxSynchronize());
				}
			}

			//force mb
			if (config->useManyBody && mb.nlargestsize>0)
			{
				if (config->manyBodyBPP && !fallbackmb) //bpp and no fallback
				{
					CU_SAFE_CALL(cuLaunchGrid(manyBodyForcesBPP1, config->ManyBodyParticles,1));
					CU_SAFE_CALL(cuCtxSynchronize());
					CU_SAFE_CALL(cuLaunchGrid(manyBodyForcesBPP2, config->ManyBodyParticles,1));
					CU_SAFE_CALL(cuCtxSynchronize());
				}
				else
				{
					CU_SAFE_CALL(cuLaunchGrid(manyBodyForces1, BlocksPerGrid,1));
					CU_SAFE_CALL(cuCtxSynchronize());
					CU_SAFE_CALL(cuLaunchGrid(manyBodyForces2, BlocksPerGrid,1));
					CU_SAFE_CALL(cuCtxSynchronize());
				}
			}

			//correct
			CU_SAFE_CALL(cuLaunchGrid(correct, BlocksPerGrid,1));
			CU_SAFE_CALL(cuCtxSynchronize());

			//copy flag mem
			CU_SAFE_CALL(cuMemcpyDtoH(flag, devFlag, sizeof(bool)));

			if (*flag)
			{
				if (config->useLennardJones && lj.nlargestsize>0) //calculate force for lennard jones
				{
					if (config->lennardJonesBPP && !fallbacklj) //bpp and no fallback
					{
						CU_SAFE_CALL(cuLaunchGrid(lennardJonesForcesBPP, (config->LennardJonesParticles + config->ManyBodyParticles), 1));
						CU_SAFE_CALL(cuCtxSynchronize()); //synchronize context
					}
					else //tpp
					{
						CU_SAFE_CALL(cuLaunchGrid(lennardJonesForces, BlocksPerGrid, 1));
						CU_SAFE_CALL(cuCtxSynchronize()); //synchronize context
					}
				}

				if (config->useManyBody && mb.nlargestsize>0) //calculate force for many body
				{
					if (config->manyBodyBPP && !fallbackmb) //bpp and no fallback
					{
						CU_SAFE_CALL(cuLaunchGrid(manyBodyForcesBPP1, config->ManyBodyParticles, 1)); //b-o 1
						CU_SAFE_CALL(cuCtxSynchronize());
						CU_SAFE_CALL(cuLaunchGrid(manyBodyForcesBPP2, config->ManyBodyParticles, 1)); //b-o 2
						CU_SAFE_CALL(cuCtxSynchronize());
					}
					else //tpp
					{
						CU_SAFE_CALL(cuLaunchGrid(manyBodyForces1, BlocksPerGrid, 1)); //b-o 1
						CU_SAFE_CALL(cuCtxSynchronize());
						CU_SAFE_CALL(cuLaunchGrid(manyBodyForces2, BlocksPerGrid, 1)); //b-o 2
						CU_SAFE_CALL(cuCtxSynchronize());
					}
				}

				//calculate the accelerations after the force
				CU_SAFE_CALL(cuLaunchGrid(calculateAccelerations, BlocksPerGrid, 1));
				CU_SAFE_CALL(cuCtxSynchronize());

				//clear flag
				CU_SAFE_CALL(cuMemsetD8(devFlag, 0, sizeof(bool)));
			}

			//results
			if (printout==timesteps)
			{
	#ifdef USECUDADEBUG
				performOutput(resultVec, BlocksPerGrid, informationSize, informationMemory, devInformation, performCalculations, calculatePotentional, timesteps, config->DT, config->Debug);
	#else
				performOutput(&resultVec, (config->LennardJonesParticles + config->ManyBodyParticles), posArray, velocityArray, devPosArray, devVelocityArray, timesteps, config->DT, config->Debug);
	#endif
				printout += config->OutputTimesteps;
			}

			//animation
			if (animts==timesteps)
			{
				//copy memory from device
				CU_SAFE_CALL(cuMemcpyDtoH(posArray, devPosArray, sizeof(real4) * (config->LennardJonesParticles + config->ManyBodyParticles)));
				CU_SAFE_CALL(cuMemcpyDtoH(velocityArray, devVelocityArray, sizeof(real3) * (config->LennardJonesParticles + config->ManyBodyParticles)));

				pushAnim(posVec, velVec, posArray, velocityArray, (config->LennardJonesParticles + config->ManyBodyParticles));
				animts += config->animts;
			}
			
			timesteps++;
		}

		//stop measuring
		CU_SAFE_CALL(cuEventRecord(stop,0));
		CU_SAFE_CALL(cuEventSynchronize(stop));

		//calculate time
		CU_SAFE_CALL(cuEventElapsedTime(&elapsedTime,start,stop));

		//after sim output
		cout << "\nAfter Simulation [VXCM = " << vxcm << "]:" << endl;
#ifdef USECUDADEBUG
		performOutput(resultVec, BlocksPerGrid, informationSize, informationMemory, devInformation, performCalculations, calculatePotentional, timesteps, config->DT, true);
#else
		performOutput(&resultVec, (config->LennardJonesParticles + config->ManyBodyParticles), posArray, velocityArray, devPosArray, devVelocityArray, timesteps, config->DT, true);
#endif
		cout << "\nTime took: " << elapsedTime << "ms" << endl;

		//copy memory from device
		CU_SAFE_CALL(cuMemcpyDtoH(posArray, devPosArray, sizeof(real4) * (config->LennardJonesParticles + config->ManyBodyParticles)));
		CU_SAFE_CALL(cuMemcpyDtoH(velocityArray, devVelocityArray, sizeof(real3) * (config->LennardJonesParticles + config->ManyBodyParticles)));
	    
		//write the output
		writeOutput(config, (float*) posArray, (float*) velocityArray);

		//output info to log & results
		writeResults(resultVec, config, elapsedTime, tailResults(posArray, config->ManyBodyParticles), vxcm);

		if (config->animts>-1)
			writeAnimationBinaryData(posVec, velVec, config, vxcm);

		flushVectors(resultVec, posVec, velVec);
	}

	//release lists memory
	if (ljList)
	{
		delete [] ljList;
		ljList = NULL;
	}

	if (mbList)
	{
		delete [] mbList;
		mbList = NULL;
	}
    
	//free events
	CU_SAFE_CALL(cuEventDestroy(start));
	CU_SAFE_CALL(cuEventDestroy(stop));

	//free device memory
	CU_SAFE_CALL(cuMemFree(devPosArray));
	CU_SAFE_CALL(cuMemFree(devVelocityArray));
	CU_SAFE_CALL(cuMemFree(devAAccArray));
	CU_SAFE_CALL(cuMemFree(devForceArray));
	CU_SAFE_CALL(cuMemFree(devBAccArray));
	CU_SAFE_CALL(cuMemFree(devCAccArray));
#ifdef USECUDADEBUG
	CU_SAFE_CALL(cuMemFree(devInformation));
#endif
	CU_SAFE_CALL(cuMemFree(devFlag));

	if (config->useManyBody)
		CU_SAFE_CALL(cuMemFree(devMemAlloc));

	//free memory
	CU_SAFE_CALL(cuMemFreeHost(posArray));
	CU_SAFE_CALL(cuMemFreeHost(velocityArray));
	CU_SAFE_CALL(cuMemFreeHost(flag));
#ifdef USECUDADEBUG
	CU_SAFE_CALL(cuMemFreeHost(informationMemory));
#endif

	//drop context
	CU_SAFE_CALL(cuCtxDetach(context));

	//free module path
	cutFree(module_path);

	return EXIT_SUCCESS;
}

#ifdef USECUDADEBUG
void performOutput(vector<results> & resultVec, int BlocksPerGrid, int informationSize, float * informationMemory, CUdeviceptr & devInformation, CUfunction & performCalculations, CUfunction & calculatePotentional, int timesteps, float dt, bool print)
{
	real3 cmassTmp;
	real3 momentumTmp;

	//reset information memory on device
	CU_SAFE_CALL(cuMemsetD32(devInformation, 0, informationSize/sizeof(float)));

	//kinetic, memontum an so
	CU_SAFE_CALL(cuLaunchGrid(performCalculations, 1, 1));
	CU_SAFE_CALL(cuCtxSynchronize());

	//potentional energy calculation
	CU_SAFE_CALL(cuLaunchGrid(calculatePotentional, BlocksPerGrid, 1));
	CU_SAFE_CALL(cuCtxSynchronize());

	//copy from device to host
	CU_SAFE_CALL(cuMemcpyDtoH(informationMemory, devInformation, informationSize));
	
	//fill results structure and push into vector
	results res;
	res.time = timesteps * dt;
	res.ek = (double)*(informationMemory+2);
	res.eu = (double)*informationMemory + (double)*(informationMemory+1);
	res.e = (double)res.ek + res.eu;
	res.temperature = (double)*(informationMemory+3);
	cmassTmp = *(real3*)(informationMemory+4);
	momentumTmp = *(real3*)(informationMemory+7);
	res.centerOfMassx = (double)cmassTmp.x;
	res.centerOfMassy = (double)cmassTmp.y;
	res.centerOfMassz = (double)cmassTmp.z;
	res.momentumx = (double)momentumTmp.x;
	res.momentumy = (double)momentumTmp.y;
	res.momentumz = (double)momentumTmp.z;
	
	if (print) //debug-output
		cout << "t=" << res.time 
			<< "\tT=" << res.temperature 
			<< "\tE=" << res.e
			<< "\tEK=" << res.ek
			<< "\tEU=" << res.eu
			<< endl;

	resultVec.push_back(res);
}
#else
void performOutput(vector<results> * resultVec, int NumberOfParticles, real4 * posArray, real3 * velocityArray, CUdeviceptr & devPosArray, CUdeviceptr & devVelocityArray, int timesteps, real dt, bool print)
{
	CU_SAFE_CALL(cuMemcpyDtoH(posArray, devPosArray, sizeof(real4) * NumberOfParticles));
	CU_SAFE_CALL(cuMemcpyDtoH(velocityArray, devVelocityArray, sizeof(real3) * NumberOfParticles));
	performOutput(resultVec, NumberOfParticles, posArray, velocityArray, timesteps, dt, print);
}
#endif

///
/// Get Highest Velocity & Box Size
///
void readyList(real4 * posArray, real3 * velocityArray, int NumberOfParticles, float & highestVelocity, float & boxSize)
{
	highestVelocity = 0;
	boxSize = 0;

	for(int i=0; i<NumberOfParticles; i++)
	{
		//calculate highest velocity
		float velocity = sqrt(velocityArray[i].x * velocityArray[i].x + velocityArray[i].y * velocityArray[i].y + velocityArray[i].z * velocityArray[i].z); 

		if (velocity>highestVelocity)
			highestVelocity = velocity;

		//boxsize = (the particle that is far the most from 0,0,0) * 2

		if (boxSize<abs(posArray[i].x))
			boxSize = abs(posArray[i].x);

		if (boxSize<abs(posArray[i].y))
			boxSize = abs(posArray[i].y);

		if (boxSize<abs(posArray[i].z))
			boxSize = abs(posArray[i].z);
	}

	boxSize*=2; //final box size
}

///
/// Build Neighbor List
///
void buildList(real4 * posArray, listSettings * listsettings, float boxSize, int buckets, float highestVelocity, int NumberOfParticles, int & nextBuild, int currentTimestep, int ** list, float dt, char * lpTexRef, CUarray * cu_array, CUtexref * cu_texref, CUmodule * module)
{
	CUDA_ARRAY_DESCRIPTOR desc;
	CUDA_MEMCPY2D copyParam;

	//free memory if needed
	if (*list)
	{
		delete [] *list;
		*list = NULL;
	}

	//build list
	*list = buildNeighborList(posArray, listsettings, boxSize, buckets, NumberOfParticles);

	//calculte next build time step
	nextBuild = (int) ((listsettings->maxnlmove / highestVelocity) / fabs(dt));
	nextBuild += currentTimestep;

	//create texture
	if (*cu_array)
		CU_SAFE_CALL( cuArrayDestroy( *cu_array));
	desc.Format = CU_AD_FORMAT_SIGNED_INT32;
	desc.NumChannels = 1;
	desc.Width = listsettings->nlistsize;
	desc.Height = NumberOfParticles;
	CU_SAFE_CALL( cuArrayCreate( cu_array, &desc ));
	
	//set memory copy params host-to-cuda-array
	memset(&copyParam, 0, sizeof(copyParam));
	copyParam.dstMemoryType = CU_MEMORYTYPE_ARRAY;
	copyParam.dstArray = *cu_array;
	copyParam.srcMemoryType = CU_MEMORYTYPE_HOST;
	copyParam.srcHost = *list;
	copyParam.srcPitch = listsettings->nlistsize * sizeof(int);
	copyParam.WidthInBytes = copyParam.srcPitch;
	copyParam.Height = NumberOfParticles;
	CU_SAFE_CALL(cuMemcpy2D(&copyParam));

	//build texture
	CU_SAFE_CALL(cuModuleGetTexRef(cu_texref, *module, lpTexRef));
	CU_SAFE_CALL(cuTexRefSetArray(*cu_texref, *cu_array, CU_TRSA_OVERRIDE_FORMAT));
	CU_SAFE_CALL(cuTexRefSetAddressMode(*cu_texref, 0, CU_TR_ADDRESS_MODE_CLAMP));
	CU_SAFE_CALL(cuTexRefSetAddressMode(*cu_texref, 1, CU_TR_ADDRESS_MODE_CLAMP));
	CU_SAFE_CALL(cuTexRefSetFilterMode(*cu_texref, CU_TR_FILTER_MODE_POINT));
	CU_SAFE_CALL(cuTexRefSetFlags(*cu_texref, CU_TRSF_READ_AS_INTEGER));
	CU_SAFE_CALL(cuTexRefSetFormat(*cu_texref, CU_AD_FORMAT_SIGNED_INT32, 1));

	//free memory
	if (*list)
	{
		delete [] *list;
		*list = NULL;
	}
}