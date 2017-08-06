///
/// CPU Host
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "Host.h"

///
/// Simulation main host function
///
int hostMain(configuration * config)
{
	real4 * posArray; //positions
	real3 * velocityArray; //velocity
	real3 * forceArray; //force
	real3 * aAccArray; //acceleration
	real3 * bAccArray; //b-acceleration
	real3 * cAccArray; //c-acceleration

	int * ljList = NULL; //lennard jones list pointer
	int * mbList = NULL; //many body list pointer
	listSettings lj; //lennard jones list settings
	listSettings mb; //many body list settings
	real highestVelocity; //heighest velocity
	int timesteps; //timestep counter
	int ljNextBuild; //on which timestep to call the build of the lennard jones list
	int mbNextBuild; //on which timestep to call the build of the many body list
	real boxSize=pow((config->LennardJonesParticles + config->ManyBodyParticles),1.0/3.0)*52.8; //box size
	real boxSizeList; //box size for the neighbor lists
	int buckets; //bucket count for the neighbor list
	
	vector<results> resultVec; //results vector
	vector<float*> posVec; //position vector for animation
	vector<float*> velVec; //velocity vector for animation

    //time-measure events
    DWORD start, stop, elapsed;

    //allocate memory for data on host
	TC(posArray = new real4[(config->LennardJonesParticles + config->ManyBodyParticles)]);
    TC(velocityArray = new real3[(config->LennardJonesParticles + config->ManyBodyParticles)]);
    TC(aAccArray = new real3[(config->LennardJonesParticles + config->ManyBodyParticles)]);
    TC(forceArray = new real3[(config->LennardJonesParticles + config->ManyBodyParticles)]);
	TC(bAccArray = new real3[(config->LennardJonesParticles + config->ManyBodyParticles)]);
	TC(cAccArray = new real3[(config->LennardJonesParticles + config->ManyBodyParticles)]);

	//number of buckets for list hashs
	buckets = (config->LennardJonesParticles + config->ManyBodyParticles) * 3;

	lj.nlargestsize = mb.nlargestsize = 0;

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

	config->energyLoss = sqrt(config->energyLoss); //fix energy loss

	for(float vxcm=config->vxcmFrom; vxcm<=config->vxcmTo; vxcm+=config->vxcmStep)
	{
		//init
		timesteps = 0;
		ljNextBuild = mbNextBuild = 0;

		//information
		int printout = config->OutputTimesteps; //next debug output time step
		int animts = config->animts; //next animation save time step

		//set the allocated memory to 0
		memset(aAccArray,0, sizeof(real3)*(config->LennardJonesParticles + config->ManyBodyParticles));
		memset(bAccArray,0, sizeof(real3)*(config->LennardJonesParticles + config->ManyBodyParticles));
		memset(cAccArray,0, sizeof(real3)*(config->LennardJonesParticles + config->ManyBodyParticles));

		//read input file
		readInput(config, (real*) posArray, (real*) velocityArray);

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

		//before sim output
		cout << fixed << setprecision(6); //set maximal precision for output
		cout.setf(ios::fixed,ios::floatfield); //zero padding
		cout << "Before Simulation [VXCM = " << vxcm << "]:" << endl;

		//calculate highest velocity and boxSizeList
		readyList(posArray, velocityArray, (config->LennardJonesParticles + config->ManyBodyParticles), highestVelocity, boxSizeList);

		if (config->useLennardJones) //build lennard jones list
			buildList(posArray, &lj, boxSizeList, buckets, highestVelocity, (config->LennardJonesParticles + config->ManyBodyParticles), ljNextBuild, 0, &ljList, config->DT);

		if (config->useManyBody) //build many body list
			buildList(posArray, &mb, boxSizeList, buckets, highestVelocity, config->ManyBodyParticles, mbNextBuild, 0, &mbList, config->DT);

		//perform the output function
		performOutput(&resultVec, (config->LennardJonesParticles + config->ManyBodyParticles), posArray, velocityArray, timesteps, config->DT, true);
		cout << endl;

		if (config->animts>-1)
			pushAnim(posVec, velVec, posArray, velocityArray, (config->LennardJonesParticles + config->ManyBodyParticles));

		//measure & perform
		start = GetTickCount();

		if (config->useLennardJones) //calculate force for lennard jones
			lennardJonesForces(posArray, forceArray, (config->LennardJonesParticles + config->ManyBodyParticles), lj.nlargestsize, lj.nlistsize, ljList, lj.rcutsq);
		if (config->useManyBody) //calculate force for many body
			manyBodyForces(posArray, forceArray, config->ManyBodyParticles, mb.nlargestsize, mb.nlistsize, mbList, mb.rcutsq, config->useLennardJones);

		//calculate the accelerations after the force
		calculateAccelerations(posArray, aAccArray, forceArray, (config->LennardJonesParticles + config->ManyBodyParticles));

		//while not reached the timestep goal
		while(timesteps<config->Timesteps)
		{
			if (ljNextBuild == timesteps || mbNextBuild == timesteps) //is is time to build any of the lists?
			{
				readyList(posArray, velocityArray, (config->LennardJonesParticles + config->ManyBodyParticles), highestVelocity, boxSizeList);

				if (config->useLennardJones && ljNextBuild == timesteps)
					buildList(posArray, &lj, boxSizeList, buckets, highestVelocity, (config->LennardJonesParticles + config->ManyBodyParticles), ljNextBuild, timesteps, &ljList, config->DT);

				if (config->useManyBody && mbNextBuild == timesteps)
					buildList(posArray, &mb, boxSizeList, buckets, highestVelocity, config->ManyBodyParticles, mbNextBuild, timesteps, &mbList, config->DT);
			}

			//predict
			predict(posArray, velocityArray, aAccArray, bAccArray, cAccArray, config->DT, (config->LennardJonesParticles + config->ManyBodyParticles));

			if (config->useLennardJones) //force lj
				lennardJonesForces(posArray, forceArray, (config->LennardJonesParticles + config->ManyBodyParticles), lj.nlargestsize, lj.nlistsize, ljList, lj.rcutsq);

			if (config->useManyBody) //force mb
				manyBodyForces(posArray, forceArray, config->ManyBodyParticles, mb.nlargestsize, mb.nlistsize, mbList, mb.rcutsq, config->useLennardJones);

			bool flag = false;

			//correct
			correct(posArray, velocityArray, forceArray, aAccArray, bAccArray, cAccArray, config->DT, (config->LennardJonesParticles + config->ManyBodyParticles), boxSize, flag, config->energyLoss);

			if (flag)
			{
				if (config->useLennardJones) //force lj
				lennardJonesForces(posArray, forceArray, (config->LennardJonesParticles + config->ManyBodyParticles), lj.nlargestsize, lj.nlistsize, ljList, lj.rcutsq);

				if (config->useManyBody) //force mb
					manyBodyForces(posArray, forceArray, config->ManyBodyParticles, mb.nlargestsize, mb.nlistsize, mbList, mb.rcutsq, config->useLennardJones);

				//calculate the accelerations after the force
				calculateAccelerations(posArray, aAccArray, forceArray, (config->LennardJonesParticles + config->ManyBodyParticles));
			}

			//results
			if (printout==timesteps)
			{
				performOutput(&resultVec, (config->LennardJonesParticles + config->ManyBodyParticles), posArray, velocityArray, timesteps, config->DT, config->Debug);
				printout += config->OutputTimesteps;
			}

			//animation
			if (animts==timesteps)
			{
				pushAnim(posVec, velVec, posArray, velocityArray, (config->LennardJonesParticles + config->ManyBodyParticles));
				animts += config->animts;
			}

			timesteps++;
		}

		//stop measuring
		stop = GetTickCount();

		//calculate the time took to simulate
		elapsed = stop - start;   

		//after sim output
		cout << "\nAfter Simulation [VXCM = " << vxcm << "]:" << endl;
		performOutput(&resultVec, (config->LennardJonesParticles + config->ManyBodyParticles), posArray, velocityArray, timesteps, config->DT, true);
		cout << "\nTime took: " << elapsed << "ms" << endl << endl;
	    
		//write the output
		writeOutput(config, (real*) posArray, (real*) velocityArray);
		writeResults(resultVec, config, (float)elapsed, tailResults(posArray, config->ManyBodyParticles), vxcm);

		if (config->animts > -1)
			writeAnimationBinaryData(posVec, velVec, config, vxcm);

		flushVectors(resultVec, posVec, velVec);
	}

    //free memory
    delete [] posArray;
    delete [] velocityArray;
    delete [] aAccArray;
    delete [] forceArray;
	delete [] bAccArray;
	delete [] cAccArray;

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

	return EXIT_SUCCESS;
}

///
/// Get Highest Velocity & Box Size
///
void readyList(real4 * posArray, real3 * velocityArray, int NumberOfParticles, real & highestVelocity, real & boxSize)
{
	highestVelocity = 0;
	boxSize = 0;

	for(int i=0; i<NumberOfParticles; i++)
	{
		//calculate highest velocity
		real velocity = sqrt(velocityArray[i].x * velocityArray[i].x + velocityArray[i].y * velocityArray[i].y + velocityArray[i].z * velocityArray[i].z); 

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
void buildList(real4 * posArray, listSettings * listsettings, real boxSize, int buckets, real highestVelocity, int NumberOfParticles, int & nextBuild, int currentTimestep, int ** list, real dt)
{
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
}