///
/// CUDA Device Information Function
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#define VIBLJ memory
#define VIBMB memory+1
#define KINETIC memory[2]
#define TEMPERATURE memory[3]
#define CENTEROFMASS *(float3*)(memory+4)
#define MOMENTUM *(float3*)(memory+7)

///
/// atomicAdd for floats
///
extern "C" __device__ inline void atomicAdd(float* address, float value)
{
	while ((value = atomicExch(address, atomicExch(address, 0.0f)+value))!=0.0f);
}

///
/// Perform Calculations
///
extern "C" __global__ void performCalculations(float4 * posArray, float3 * velocityArray, int NumberOfParticles, float * memory)
{
	float mass;
	float3 cmv = {0.0f, 0.0f, 0.0f};
	float3 centerOfMass = {0.0f, 0.0f, 0.0f};
	float total = 0;
	float kinetic=0;
	float temp=0;

	float4 pos;
	float3 vel;

    for(int i=0;i<NumberOfParticles;i++)
	{
		pos = posArray[i];
		vel = velocityArray[i];

		mass = getMass(pos.w); //get mass for particle

		//kinetic energy
		kinetic+=mass*(vel.x*vel.x+vel.y*vel.y+vel.z*vel.z);

		//momentum
		cmv.x += mass*vel.x;
		cmv.y += mass*vel.y;
		cmv.z += mass*vel.z;

		//total for temperature & center of mass
		total += mass;

		centerOfMass.x += mass*pos.x;
		centerOfMass.y += mass*pos.y;
		centerOfMass.z += mass*pos.z;
	}

	MOMENTUM = cmv; //save momentum to memory

	centerOfMass.x /= total;
	centerOfMass.y /= total;
	centerOfMass.z /= total;

	CENTEROFMASS = centerOfMass; //save center of mass to memory

	cmv.x /= total;
	cmv.y /= total;
	cmv.z /= total;

	for(int i=0; i<NumberOfParticles; i++)
	{
		pos = posArray[i];
		vel = velocityArray[i];

		mass = getMass(pos.w);
		temp += mass*((vel.x-cmv.x)*(vel.x-cmv.x) + (vel.y-cmv.y)*(vel.y-cmv.y) + (vel.z-cmv.z)*(vel.z-cmv.z));
	}

	//save to global
	TEMPERATURE = temp / (3.0f*((float)NumberOfParticles)*8.314e-7f);
	KINETIC = kinetic / 2;
}

///
/// Calculate Potentional Energy (we calculate it again, cause we need to have all the particles)
///
extern "C" __global__ void calculatePotentional(float4 * posArray, int NumberOfParticles, float * memory)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	if (id>=NumberOfParticles) //more thread than particles
		return;

	float vibLJ = 0;
	float vibMB = 0;
	float3 ij;
	float rijsq;

	float4 ipos = posArray[id];

	int itype = (int)ipos.w;

	for(int j=0; j<NumberOfParticles; j++)
	{
		if (id==j) //same particle
			continue;

		float4 jpos = posArray[j];

		int jtype = (int)jpos.w;

		ij.x = ipos.x - jpos.x;
		ij.y = ipos.y - jpos.y;
		ij.z = ipos.z - jpos.z;

		rijsq = ij.x*ij.x + ij.y*ij.y + ij.z*ij.z;

		if (jtype == NG || itype == NG)
		{
			int select = (itype == NG) ? jtype : itype;  //select type

			float sig;
			float eps;
			
			if (select==N) //select parameters
			{
				eps = epsNGN;
				sig = sigNGN;
			}
			else if (select==O)
			{
				eps = epsNGO;
				sig = sigNGO;
			}
			else
			{
				eps = epsNG;
				sig = sigNG;
			}

			float sr2 = sig*sig/rijsq;
			float sr6 = sr2*sr2*sr2;
			vibLJ += 2.0f*eps*(sr6*sr6 - sr6);
		}
		else
		{
			float va=0, dva, vr=0, dvr, vdw=0, dvdw;
			float rij = sqrt(rijsq);
			float bo=0;
			float fcut;
			float dfcut;

			if (rij > CUTOFF_RIJ) //distance is too large... ignoring...
				continue;

			cutoff(rij/10, fcut, dfcut);

			if (fcut>0)
			{
				float de = getDE(itype, jtype);
				float re = getRE(itype, jtype);
				float beta = getbeta(itype, jtype);

				attractive(va, dva, rij, de, beta , re);
				repulsive(vr, dvr, rij, de, beta, re);

				float3 im;
				float3 jm;
				float rim;
				float rjm;
				float2 xs = {0.0f, 0.0f};

				for(int m=0; m<NumberOfParticles; m++)
				{
					float4 mpos = posArray[m];

					int mtype = (int)mpos.w;

					if (mtype == NG)
						continue;

					if (m==j || m==id)
						continue;

					im.x = ipos.x - mpos.x;
					im.y = ipos.y - mpos.y;
					im.z = ipos.z - mpos.z;

					rim = sqrt(im.x*im.x + im.y*im.y + im.z*im.z);

					jm.x = jpos.x - mpos.x;
					jm.y = jpos.y - mpos.y;
					jm.z = jpos.z - mpos.z;

					rjm = sqrt(jm.x*jm.x + jm.y*jm.y + jm.z*jm.z);

					float fci;
					float dfci;
					float fcj;
					float dfcj;
		
					cutoff(rim, fci, dfci);
					cutoff(rjm, fcj, dfcj);

					if (fci>0)
					{
						float icosthe = (rij*rij+rim*rim-rjm*rjm)/ (2.0f*rij*rim);
						float igthe = 1.0f +C/D -C/(D+(H-icosthe)*(H-icosthe));
						xsi += fci*igthe*exp(LAMDA3*(rij - re));
					}

					if (fcj>0)
					{
						float jcosthe = (rij*rij+rjm*rjm-rim*rim)/ (2.0f*rij*rjm);
						float jgthe = 1.0f +C/D -C/(D+(H-jcosthe)*(H-jcosthe));	
						xsj += fcj*jgthe*exp(LAMDA3*(rij - re));
					}
				}	

				float bij = 1.0f/(1.0f + xsi);
				float bji = 1.0f/(1.0f + xsj);

				bo = 0.5f*(bij + bji);
			}

			float r0v = getR0V(itype, jtype);
			vanderWaals(vdw, dvdw, rij, r0v);

			vibMB += 0.5f*((vr - bo*va) + vdw);
		}
	}

	atomicAdd(VIBMB, vibMB);
	atomicAdd(VIBLJ, vibLJ);
}