///
/// CPU Host
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "Host.h"

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

	//init
	timesteps = 1;
	ljNextBuild = mbNextBuild = 0;

	//information
	int printout = config->OutputTimesteps; //next debug output time step

	//set the allocated memory to 0
	memset(posArray, 0, sizeof(real4)*(config->LennardJonesParticles + config->ManyBodyParticles));
	memset(velocityArray, 0, sizeof(real3)*(config->LennardJonesParticles + config->ManyBodyParticles));
	memset(aAccArray, 0, sizeof(real3)*(config->LennardJonesParticles + config->ManyBodyParticles));
	memset(forceArray, 0, sizeof(real3)*(config->LennardJonesParticles + config->ManyBodyParticles));
	memset(bAccArray, 0, sizeof(real3)*(config->LennardJonesParticles + config->ManyBodyParticles));
	memset(cAccArray, 0, sizeof(real3)*(config->LennardJonesParticles + config->ManyBodyParticles));

	srand((unsigned int)time(NULL));

	real rb = 8.314e-7;
	real temp = config->Temperature*rb;  
	real boxl = boxSize;

	cout << boxl << endl;

	if (config->useManyBody)
		if (!choosePositionsManyBody(posArray, boxl, config))
			return EXIT_FAILURE;

	if (config->useLennardJones)
		if (!choosePositionsNoble(posArray, boxl, config))
			return EXIT_FAILURE;

	chooseVelocities(posArray, velocityArray, temp, config);

	//before sim output
	cout << fixed << setprecision(6); //set maximal precision for output
	cout.setf(ios::fixed,ios::floatfield); //zero padding
	cout << "Before Simulation" << endl;

	//calculate highest velocity and boxSizeList
	readyList(posArray, velocityArray, (config->LennardJonesParticles + config->ManyBodyParticles), highestVelocity, boxSizeList);

	if (config->useLennardJones) //build lennard jones list
		buildList(posArray, &lj, boxSizeList, buckets, highestVelocity, (config->LennardJonesParticles + config->ManyBodyParticles), ljNextBuild, 0, &ljList, config->DT);

	if (config->useManyBody) //build many body list
		buildList(posArray, &mb, boxSizeList, buckets, highestVelocity, config->ManyBodyParticles, mbNextBuild, 0, &mbList, config->DT);

	//perform the output function
	performOutput(&resultVec, (config->LennardJonesParticles + config->ManyBodyParticles), posArray, velocityArray, timesteps, config->DT, true);
	cout << endl;

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

		if(timesteps%50==0)
			chooseVelocities(posArray, velocityArray, temp, config);

		timesteps++;
	}

	//stop measuring
	stop = GetTickCount();

	//calculate the time took to simulate
	elapsed = stop - start;   

	//after sim output
	cout << "\nAfter Equi" << endl;
	performOutput(&resultVec, (config->LennardJonesParticles + config->ManyBodyParticles), posArray, velocityArray, timesteps, config->DT, true);
	cout << "\nTime took: " << elapsed << "ms" << endl << endl;
    
	//write the output
	writeOutput(config, (real*) posArray, (real*) velocityArray);
	writeResults(resultVec, config, (float)elapsed, tailResults(posArray, config->ManyBodyParticles), 0);

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

bool choosePositionsManyBody(real4 * posArray, real boxl, configuration * config)
{
	real nl=0;
	real re = reN;
	real mass = N;

	real mindist = 10;

	real rsqminNN = mindist*mindist*sigNvdw*sigNvdw/(boxl*boxl);
	real rsqminOO = mindist*mindist*sigOvdw*sigOvdw/(boxl*boxl);
	real rsqminNO = mindist*mindist*0.25*(sigNvdw+sigOvdw)*(sigNvdw+sigOvdw)/(boxl*boxl);

	for(int i=0;i<config->ManyBodyParticles;)
	{
		if (i==config->ManyBodyParticles/2)
		{
			re = reO;
			mass = O;
		}

		posArray[i].x = (((real)rand() + 1) / (real)MD_RANDMAX) - 0.5;
		posArray[i].y = (((real)rand() + 1) / (real)MD_RANDMAX) - 0.5;
		posArray[i].z = (((real)rand() + 1) / (real)MD_RANDMAX) - 0.5;	
		posArray[i].w = mass;

		real the = (((real)rand() + 1) / (real)MD_RANDMAX)*M_PI;
		real phi = (((real)rand() + 1) / (real)MD_RANDMAX)*2.0*M_PI;
		real xsi = (((real)rand() + 1) / (real)MD_RANDMAX)*2.0*M_PI;

		posArray[i+1].x = posArray[i].x + re/boxl*sin(the)*cos(phi);
		posArray[i+1].y = posArray[i].y + re/boxl*sin(the)*sin(phi);
		posArray[i+1].z = posArray[i].z + re/boxl*cos(the);
		posArray[i+1].w = mass;

		int j;
		for(j=0;j<i;)
		{
			real rxij = posArray[i].x - posArray[j].x;
			real ryij = posArray[i].y - posArray[j].y;
			real rzij = posArray[i].z - posArray[j].z;
			real dist1 = rxij*rxij + ryij*ryij + rzij*rzij;

			rxij = posArray[i+1].x - posArray[j].x;
			ryij = posArray[i+1].y - posArray[j].y;
			rzij = posArray[i+1].z - posArray[j].z;
			real dist2 = rxij*rxij + ryij*ryij + rzij*rzij;

			if (i<config->ManyBodyParticles/2)
			{
				if((dist1>rsqminNN) && (dist2>rsqminNN))
					j++;
				else
					break;
			}
			else
			{
				if(j<config->ManyBodyParticles/2)
				{
					if((dist1>rsqminNO) && (dist2>rsqminNO))
					  j++;
					else
					  break;
				}
				else
				{
					if((dist1>rsqminOO) && (dist2>rsqminOO))
					  j++;
					else
					  break;
				}   
			}
		}

		if(j==i)
			i = i + 2;

		nl++;

		if(nl>10e11)
		{
			cerr << "Can't put the particles in the box...\n";
			return false;
		}
	}
	
	for(int i=0;i<config->ManyBodyParticles;i++)
	{
		posArray[i].x = posArray[i].x*boxl;
		posArray[i].y = posArray[i].y*boxl;
		posArray[i].z = posArray[i].z*boxl;
	}
	
	return true;
}

bool choosePositionsNoble(real4 * posArray, real boxl, configuration * config)
{
	real rsqminNGNG = 10.0*10.0*sigNG*sigNG/(boxl*boxl);
	int nl=0;
	int i,j;

	for(i=config->ManyBodyParticles;i<(config->LennardJonesParticles + config->ManyBodyParticles);i++)
	{
		posArray[i].x = (((real)rand() + 1) / (real)MD_RANDMAX) - 0.5;
		posArray[i].y = (((real)rand() + 1) / (real)MD_RANDMAX) - 0.5;
		posArray[i].z = (((real)rand() + 1) / (real)MD_RANDMAX) - 0.5;	
		posArray[i].w = NG;

		for(j=0;j<i;)
		{
			real rxij = posArray[i].x - posArray[j].x;
			real ryij = posArray[i].y - posArray[j].y;
			real rzij = posArray[i].z - posArray[j].z;
			real dist1 = rxij*rxij + ryij*ryij + rzij*rzij;

			if(dist1>rsqminNGNG)
				j++;
			else
				break;   
		}

		if(j==i)
			i = i++;

		nl++;

		if(nl>10e11)
		{
			cerr << "Can't put the particles in the box...\n";
			return false;
		}
	}   

	for(i=config->ManyBodyParticles;i<(config->LennardJonesParticles + config->ManyBodyParticles);i++)
	{
		posArray[i].x = posArray[i].x*boxl;
		posArray[i].y = posArray[i].y*boxl;
		posArray[i].z = posArray[i].z*boxl;
	}

	return true;
}

void chooseVelocities(real4 * posArray, real3 * velocityArray, real temp, configuration * config)
{
	real cmvx=0;
	real cmvy=0;
	real cmvz=0;
	real totalm=0;
	real t=0;

	for(int i=0; i<(config->LennardJonesParticles + config->ManyBodyParticles); i++)
	{
		real mass = 0;

		switch((int)posArray[i].w)
		{
		case N:
			mass = mN;
			break;
		case O:
			mass = mO;
			break;
		case NG:
			mass = mNG;
			break;
		}

		real y = ((real)rand() + 1) / (real)MD_RANDMAX;

		y = 2.0*y - 1.0;

		y = inverf(y);

		velocityArray[i].x = y*sqrt(2.0*temp/mass);

		y = ((real)rand() + 1) / (real)MD_RANDMAX;

		y = 2.0*y - 1.0;

		y = inverf(y);

		velocityArray[i].y = y*sqrt(2.0*temp/mass);

		y = ((real)rand() + 1) / (real)MD_RANDMAX;
		
		y = inverf(y);

		y = 2.0*y - 1.0;

		velocityArray[i].z = y*sqrt(2.0*temp/mass);

		if (i<config->ManyBodyParticles)
		{
			velocityArray[i+1].x = velocityArray[i].x;
			velocityArray[i+1].y = velocityArray[i].y;
			velocityArray[i+1].z = velocityArray[i].z;

			i++;
		}

		cmvx = cmvx + mass*velocityArray[i].x;
		cmvy = cmvy + mass*velocityArray[i].y;

		cmvz = cmvz + mass*velocityArray[i].z;
		totalm+=mass;
	}

	cmvx = cmvx/totalm;
	cmvy = cmvy/totalm;
	cmvz = cmvz/totalm;

	while((fabs(cmvx)>1.0e-15)||(fabs(cmvy)>1.0e-15)||(fabs(cmvz)>1.0e-15))
	{
		for(int i=0;i<(config->LennardJonesParticles + config->ManyBodyParticles);i++)
		{
			velocityArray[i].x = velocityArray[i].x - cmvx;
			velocityArray[i].y = velocityArray[i].y - cmvy;
			velocityArray[i].z = velocityArray[i].z - cmvz;
		}

		cmvx = 0.0;
		cmvy = 0.0;
		cmvz = 0.0;

		for(int i=0; i< (config->LennardJonesParticles + config->ManyBodyParticles); i++)
		{

			real mass = 0;
			
			switch((int)posArray[i].w)
			{
			case N:
				mass = mN;
				break;
			case O:
				mass = mO;
				break;
			case NG:
				mass = mNG;
				break;
			}

			cmvx = cmvx + velocityArray[i].x*mass;
			cmvy = cmvy + velocityArray[i].y*mass;
			cmvz = cmvz + velocityArray[i].z*mass;
		}
    
		cmvx = cmvx/totalm;
		cmvy = cmvy/totalm;
		cmvz = cmvz/totalm;
	}

	for(int i=0;i<(config->LennardJonesParticles + config->ManyBodyParticles);i++)
	{
		real mass = 0;
			
			switch((int)posArray[i].w)
			{
			case N:
				mass = mN;
				break;
			case O:
				mass = mO;
				break;
			case NG:
				mass = mNG;
				break;
			}
	
			t+=mass*(velocityArray[i].x*velocityArray[i].x 
			+ velocityArray[i].y*velocityArray[i].y
			+ velocityArray[i].z*velocityArray[i].z);
	}

	t = t/(3.0*(real)((config->LennardJonesParticles + config->ManyBodyParticles)));
	real scalev = sqrt(temp/t);

	for(int i=0;i<(config->LennardJonesParticles + config->ManyBodyParticles);i++)
	{
		velocityArray[i].x = velocityArray[i].x*scalev;
		velocityArray[i].y = velocityArray[i].y*scalev;
		velocityArray[i].z = velocityArray[i].z*scalev;
	}

	/*static int index = 0;

	stringstream ss;

	ss << "vel" << index++ << ".txt";

	ofstream f(ss.str().c_str());

	for(int i=0;i<config->ManyBodyParticles;i++)
	{
		real vv = sqrt(velocityArray[i].x*velocityArray[i].x + velocityArray[i].y*velocityArray[i].y + velocityArray[i].z*velocityArray[i].z);
		f << velocityArray[i].x << " " << velocityArray[i].y << " " << velocityArray[i].z << " " << vv << endl;
	}

	f.close();*/
}