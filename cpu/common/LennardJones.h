///
/// CPU Device Lennard Jones Function
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

///
/// Lennard-Jones potentional
///
void lennardJones(real3 & force, real dist, real3 ij, real sig, real eps)
{
	real sigsq = sig*sig;
	real con = 24.0f*eps/sigsq;
    real dist2 = dist * dist;
	dist2 /= sigsq;

    real dist4 = dist2*dist2;
    real dist8 = dist4*dist4;
    real dist14 = dist2*dist4*dist8;
    real invdist8= 1.0f/dist8;
    real invdist14= 1.0f/dist14;
    real s = __fmul_rn(2,invdist14)-invdist8;
	real fij = s * con;

	force.x += ij.x * fij;
    force.y += ij.y * fij;
    force.z += ij.z * fij;
}

///
/// Lennard Jones Forces
///
void lennardJonesForces(real4 * posArray, real3 * forceArray, int NumberOfParticles, int nlargestsize, int nlistsize, int * list, real rcutsq)
{
	real3 force;
	real3 ij;

	for(int i=0; i<NumberOfParticles; i++)
	{
		force.x = force.y = force.z = 0.0;

		for(int j=0; j<nlargestsize && list[nlistsize*i+j]!=-1; j++)
		{
			int jp = list[nlistsize*i+j];
				
			ij.x = posArray[i].x - posArray[jp].x;
			ij.y = posArray[i].y - posArray[jp].y;
			ij.z = posArray[i].z - posArray[jp].z;

			real rij = ij.x*ij.x + ij.y*ij.y + ij.z*ij.z;

			if (rij>=rcutsq)
				continue;

			rij = sqrt(rij);

			ij.x = posArray[i].x - posArray[jp].x;
			ij.y = posArray[i].y - posArray[jp].y;
			ij.z = posArray[i].z - posArray[jp].z;
			
			int jtype = (int)posArray[jp].w;

			if (jtype == NG || (int)posArray[i].w == NG)
			{
				int select = ((int)posArray[i].w == NG) ? jtype : (int)posArray[i].w; 

				real sig;
				real eps;
				
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

				lennardJones(force, rij, ij, sig, eps); 
			}
		}

		forceArray[i] = force;
	}
}