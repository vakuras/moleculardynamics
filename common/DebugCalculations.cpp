///
/// Debug Calculations Functions
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "DebugCalculations.h"

///
/// Output information
///
void performOutput(vector<results> * resultVec, int NumberOfParticles, real4 * posArray, real3 * velocityArray, int timesteps, real dt, bool print)
{
	results res;
	double vibLJ;
	double vibMB;

	//kinetic, memontum an so
	PerformCalculations(posArray, velocityArray, NumberOfParticles, res.centerOfMassx, res.centerOfMassy, res.centerOfMassz, res.ek, res.momentumx, res.momentumy, res.momentumz, res.temperature);

	//potentional energy calculation
	calculatePotentional(posArray, &vibLJ, &vibMB, NumberOfParticles);

	//fill results structure and push into vector
	res.eu = (double)vibMB + (double)vibLJ;
	res.e = (double)res.ek + res.eu;
	res.time = timesteps * dt;;

	if (resultVec)
		resultVec->push_back(res);
	
	if (print)
		cout << "t=" << res.time 
			<< "\tT=" << res.temperature 
			<< "\tE=" << res.e
			<< "\tEK=" << res.ek
			<< "\tEU=" << res.eu
			<< endl;
}

///
/// Perform Calculations
///
void PerformCalculations(real4 * posArray, real3 * velocityArray, int NumberOfParticles, double & centerOfMassx, double & centerOfMassy, double & centerOfMassz, double & kineticEnergy, double & momentumx, double & momentumy, double & momentumz, double & temperature)
{
	double mass;
	double cmvx = 0; 
	double cmvy = 0; 
	double cmvz = 0;
	double total = 0;

	kineticEnergy = 0;
	temperature = 0;
	centerOfMassx = 0;
	centerOfMassy = 0;
	centerOfMassz = 0;

    for(int i=0;i<NumberOfParticles;i++)
	{
		switch((int)posArray[i].w) //select mass based on particle index
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

		//kinetic energy
		kineticEnergy+=mass*(velocityArray[i].x*velocityArray[i].x+velocityArray[i].y*velocityArray[i].y+velocityArray[i].z*velocityArray[i].z);

		//momentum
		cmvx += mass*velocityArray[i].x;
		cmvy += mass*velocityArray[i].y;
		cmvz += mass*velocityArray[i].z;

		//total for temperature & center of mass
		total += mass;

		centerOfMassx += mass*posArray[i].x;
		centerOfMassy += mass*posArray[i].y;
		centerOfMassz += mass*posArray[i].z;
	}

	centerOfMassx /= total;
	centerOfMassy /= total;
	centerOfMassz /= total;

	momentumx = cmvx;
	momentumy = cmvy;
	momentumz = cmvz;

	cmvx /= total;
	cmvy /= total;
	cmvz /= total;

	for(int i=0; i<NumberOfParticles; i++)
	{
		switch((int)posArray[i].w)
		{
		case N:
			mass = mN;
			break;
		case O:
			mass = mO;
			break;
		default:
			mass = mNG;
		}

		temperature += mass*((velocityArray[i].x-cmvx)*(velocityArray[i].x-cmvx) +
                      		(velocityArray[i].y-cmvy)*(velocityArray[i].y-cmvy) +
              		  		(velocityArray[i].z-cmvz)*(velocityArray[i].z-cmvz));
	}

	temperature /= (3.0f*((real)NumberOfParticles)*8.314e-7f);
	kineticEnergy /=2;
}

///
/// Calculate Potentional Energy
///
void calculatePotentional(real4 * posArray, double * vibLJ, double * vibMB, int NumberOfParticles)
{
	*vibLJ = 0;
	*vibMB = 0;
	double ijx;
	double ijy;
	double ijz;
	double rijsq;

	for(int i=0; i<NumberOfParticles; i++)
	{
		int itype = (int)posArray[i].w;

		for(int j=0; j<NumberOfParticles; j++)
		{
			if (i==j)
				continue;

			int jtype = (int)posArray[j].w;

			ijx = posArray[i].x - posArray[j].x;
			ijy = posArray[i].y - posArray[j].y;
			ijz = posArray[i].z - posArray[j].z;

			rijsq = ijx*ijx + ijy*ijy + ijz*ijz;

			if (jtype == NG || itype == NG)
			{
				int select = (itype == NG) ? jtype : itype; 

				double sig;
				double eps;
				
				if (select==N)
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

				double sr2 = sig*sig/rijsq;
				double sr6 = sr2*sr2*sr2;
				*vibLJ += 2.0f*eps*(sr6*sr6 - sr6);
			}
			else
			{
				double va=0, dva, vr=0, dvr, vdw, dvdw;
				double rij = sqrt(rijsq);
				double bo=0;

				if (rij > CUTOFF_RIJ)
					continue;

				double fcut, dfcut;
				cutoff(rij/10, fcut, dfcut);

				if (fcut > 0)
				{
					double de = getDE<double>(itype, jtype);
					double re = getRE<double>(itype, jtype);
					double beta = getbeta<double>(itype, jtype);

					attractive(va, dva, rij, de, beta , re);
					repulsive(vr, dvr, rij, de, beta, re);

					double imx;
					double imy;
					double imz;
					double jmx;
					double jmy;
					double jmz;
					double rim;
					double rjm;
					double xsi = 0, xsj = 0;

					for(int m=0; m<NumberOfParticles; m++)
					{
						int mtype = (int)posArray[m].w;

						if (mtype == NG)
							continue;

						if (m==j || m==i)
							continue;

						imx = posArray[i].x - posArray[m].x;
						imy = posArray[i].y - posArray[m].y;
						imz = posArray[i].z - posArray[m].z;

						rim = sqrt(imx*imx + imy*imy + imz*imz);

						jmx = posArray[j].x - posArray[m].x;
						jmy = posArray[j].y - posArray[m].y;
						jmz = posArray[j].z - posArray[m].z;

						rjm = sqrt(jmx*jmx + jmy*jmy + jmz*jmz);

						double fci;
						double dfci;
						double fcj;
						double dfcj;
			
						cutoff(rim, fci, dfci);
						cutoff(rjm, fcj, dfcj);

						if (fci>0)
						{
							double icosthe = (rij*rij+rim*rim-rjm*rjm)/ (2.0f*rij*rim);
							double igthe = 1.0f +C/D -C/(D+(H-icosthe)*(H-icosthe));
							xsi += fci*igthe*exp(LAMDA3*(rij - re));
						}

						if (fcj>0)
						{
							double jcosthe = (rij*rij+rjm*rjm-rim*rim)/ (2.0f*rij*rjm);
							double jgthe = 1.0f +C/D -C/(D+(H-jcosthe)*(H-jcosthe));	
							xsj += fcj*jgthe*exp(LAMDA3*(rij - re));
						}
					}	

					double bij = 1.0f/(1.0f + xsi);
					double bji = 1.0f/(1.0f + xsj);

					bo = 0.5f*(bij + bji);
				}

				vanderWaals(vdw, dvdw, rij, getR0V<double>(itype, jtype));

				*vibMB += 0.5f*((vr - bo*va) + vdw);
			}
		}
	}
}