///
/// CPU Device functions
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "Device.h"

#include "BondOrderFuncs.h"
#include "LennardJones.h"

///
/// Predict function - predicts where the particle will move
///
void predict(real4 *r,real3 *v,real3 *a,real3 *b,real3 *c, real dt, int nop)
{
	real c1 = dt;
    real c2 = c1*dt/2.0f;
    real c3 = c2*dt/3.0f;
    real c4 = c3*dt/4.0f;

	for(int i = 0; i<nop; i++)
	{
    	r[i].x = r[i].x + c1*v[i].x + c2*a[i].x + c3*b[i].x + c4*c[i].x;
    	r[i].y = r[i].y + c1*v[i].y + c2*a[i].y + c3*b[i].y + c4*c[i].y;
    	r[i].z = r[i].z + c1*v[i].z + c2*a[i].z + c3*b[i].z + c4*c[i].z;
		v[i].x = v[i].x + c1*a[i].x + c2*b[i].x + c3*c[i].x;
		v[i].y = v[i].y + c1*a[i].y + c2*b[i].y + c3*c[i].y;    
		v[i].z = v[i].z + c1*a[i].z + c2*b[i].z + c3*c[i].z;    
		a[i].x = a[i].x + c1*b[i].x + c2*c[i].x;
		a[i].y = a[i].y + c1*b[i].y + c2*c[i].y; 
		a[i].z = a[i].z + c1*b[i].z + c2*c[i].z;
		b[i].x = b[i].x + c1*c[i].x;
		b[i].y = b[i].y + c1*c[i].y;
		b[i].z = b[i].z + c1*c[i].z;
	}
}

///
/// Correct function - based on the calculated force it corrects the particle projection
///
void correct(real4 *r,real3 *v,real3 *f,real3 *a, real3 *b,real3 *c, real dt, int nop, real nPos, bool & flag, float energyLoss)
{
	flag = false;

	real c1 = dt ;
	real c2 = c1*dt/2.0f;
	real c3 = c2*dt/3.0f; 
	real c4 = c3*dt/4.0f;

	real cr = GEAR1*c2;
	real cv = GEAR2*c2/c1;
	real cb = GEAR3*c2/c3;
	real cc = GEAR4*c2/c4;

	for(int i = 0; i<nop; i++)
	{
		real mass = getMass<real>((int)r[i].w);

		real axi = f[i].x/mass;
      	real ayi = f[i].y/mass;
      	real azi = f[i].z/mass;

		real corrx = axi - a[i].x;
    	real corry = ayi - a[i].y;
    	real corrz = azi - a[i].z;

		r[i].x = r[i].x + cr*corrx;
		r[i].y = r[i].y + cr*corry;
		r[i].z = r[i].z + cr*corrz;
		v[i].x = v[i].x + cv*corrx;
		v[i].y = v[i].y + cv*corry;
		v[i].z = v[i].z + cv*corrz;
		a[i].x = axi;
		a[i].y = ayi;
		a[i].z = azi;
		b[i].x = b[i].x + cb*corrx;
		b[i].y = b[i].y + cb*corry;
		b[i].z = b[i].z + cb*corrz;
		c[i].x = c[i].x + cc*corrx;
		c[i].y = c[i].y + cc*corry;
    	c[i].z = c[i].z + cc*corrz;

		if (r[i].x > nPos && r[i].y > 0)
		{
			flag = true;
			v[i].x = -v[i].x*energyLoss;
			r[i].x = nPos;
		}

		if (r[i].x > nPos && r[i].y < 0)
		{
			flag = true;
			r[i].x = -nPos;
		}

		if (r[i].x < -nPos)
		{
			flag = true;
			v[i].x = -v[i].x;
			r[i].x = -nPos;
		}

		if (fabs(r[i].y) > nPos)
		{
			flag = true;
			v[i].y = -v[i].y;
			r[i].y = __fmul_rn(nPos, (r[i].y/fabs(r[i].y)));
		}

		if (fabs(r[i].z) > nPos)
		{
			flag = true;
			v[i].z = -v[i].z;
			r[i].z = __fmul_rn(nPos, (r[i].z/fabs(r[i].z)));
		}
	}
}

///
/// Calculate Accelerations
///
void calculateAccelerations(real4 * posArray, real3 * accelerationArray, real3 * forceArray, int NumberOfParticles)
{
	for(int i=0; i<NumberOfParticles; i++) 
	{
		real mass = getMass<real>((int)posArray[i].w);
		accelerationArray[i].x = forceArray[i].x / mass;
		accelerationArray[i].y = forceArray[i].y / mass;
		accelerationArray[i].z = forceArray[i].z / mass;
	}
}