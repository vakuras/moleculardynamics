///
/// CPU Device Bond Order Functions
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

///
/// Bond Order 1
///
void bondOrder1(real4 * posArray, real re, real rij, int NumberOfParticles, int nlargestsize, int nlistsize, int * list, int i, int j, real3 ij, real & bo, real3 & dbo, real rcutsq, real * xsi, real * xsj)
{
	real drijx = (ij.x)/rij;
	real drijy = (ij.y)/rij;
	real drijz = (ij.z)/rij;

	xsi[i*NumberOfParticles+j] = 0;
	real dxsix = 0.0f;
	real dxsiy = 0.0f;
	real dxsiz = 0.0f;

	xsj[i*NumberOfParticles+j] = 0;
	real dxsjx = 0.0f;
	real dxsjy = 0.0f;
	real dxsjz = 0.0f;

	real3 im;
	real3 jm;
	real rim;
	real rjm;

	for(int m=0; m<nlargestsize && list[nlistsize*i+m]!=-1; m++)
	{
		int mp = list[nlistsize*i+m];

		int mtype = (int)posArray[mp].w;

		if (mtype == NG)
			continue;

		if (mp==j)
			continue;

		im.x = posArray[i].x - posArray[mp].x;
		im.y = posArray[i].y - posArray[mp].y;
		im.z = posArray[i].z - posArray[mp].z;

		rim = im.x*im.x + im.y*im.y + im.z*im.z;

		if (rim>=rcutsq) //distance is too large... ignoring...
			continue;

		rim = sqrt(rim);

		jm.x = posArray[j].x - posArray[mp].x;
		jm.y = posArray[j].y - posArray[mp].y;
		jm.z = posArray[j].z - posArray[mp].z;

		rjm = sqrt(jm.x*jm.x + jm.y*jm.y + jm.z*jm.z);

		real fci;
		real dfci;
		real fcj;
		real dfcj;
		
		cutoff(rim, fci, dfci);
		cutoff(rjm, fcj, dfcj);

		real drimx = (posArray[i].x - posArray[mp].x)/rim;
		real drimy = (posArray[i].y - posArray[mp].y)/rim;
		real drimz = (posArray[i].z - posArray[mp].z)/rim;     

		if (fci>0)
		{
			real icosthe = (rij*rij+rim*rim-rjm*rjm)/ (2.0f*rij*rim);
			real igthe = 1.0f +C/D -C/(D+(H-icosthe)*(H-icosthe));
			xsi[i*NumberOfParticles+j] += fci*igthe*exp(LAMDA3*(rij - re));
			real dicosthex = ((2.0f*rij*drijx+2.0f*rim*drimx)*
					2.0f*rij*rim - (2.0f*rim*drijx + 2.0f*rij*drimx)*
					(rij*rij+rim*rim-rjm*rjm))/(4.0f*rij*rij*rim*rim);

			real dicosthey = ((2.0f*rij*drijy+2.0f*rim*drimy)*
					2.0f*rij*rim - (2.0f*rim*drijy + 2.0f*rij*drimy)*
					(rij*rij+rim*rim-rjm*rjm))/(4.0f*rij*rij*rim*rim);

			real dicosthez = ((2.0f*rij*drijz+2.0f*rim*drimz)*
					2.0f*rij*rim - (2.0f*rim*drijz + 2.0f*rij*drimz)*
					(rij*rij+rim*rim-rjm*rjm))/(4.0f*rij*rij*rim*rim);

			real digthex = -C/((D+(H-icosthe)*(H-icosthe))*(D+(H-icosthe)*(H-icosthe)))* 2.0f*(H-icosthe)*dicosthex;
			real digthey = -C/((D+(H-icosthe)*(H-icosthe))*(D+(H-icosthe)*(H-icosthe)))* 2.0f*(H-icosthe)*dicosthey;
			real digthez = -C/((D+(H-icosthe)*(H-icosthe))*(D+(H-icosthe)*(H-icosthe)))* 2.0f*(H-icosthe)*dicosthez;

			dxsix += (dfci*drimx*igthe*exp(LAMDA3*(rij - re)) + fci*digthex*exp(LAMDA3*(rij - re)) + fci*igthe*LAMDA3*exp(LAMDA3*(rij - re))*drijx);
			dxsiy += (dfci*drimy*igthe*exp(LAMDA3*(rij - re)) + fci*digthey*exp(LAMDA3*(rij - re)) + fci*igthe*LAMDA3*exp(LAMDA3*(rij - re))*drijy);
			dxsiz += (dfci*drimz*igthe*exp(LAMDA3*(rij - re)) + fci*digthez*exp(LAMDA3*(rij - re)) + fci*igthe*LAMDA3*exp(LAMDA3*(rij - re))*drijz);
		}

		if (fcj>0)
		{
			real jcosthe = (rij*rij+rjm*rjm-rim*rim)/ (2.0f*rij*rjm);
			real jgthe = 1.0f +C/D -C/(D+(H-jcosthe)*(H-jcosthe));
			xsj[i*NumberOfParticles+j] += fcj*jgthe*exp(LAMDA3*(rij - re));

			real djcosthex = ((2.0f*rij*drijx-2.0f*rim*drimx)*
						  2.0f*rij*rjm -2.0f*rjm*drijx*
						  (rij*rij+rjm*rjm-rim*rim))/(4.0f*rij*rij*rjm*rjm);

			real djcosthey = ((2.0f*rij*drijy-2.0f*rim*drimy)*
						  2.0f*rij*rjm -2.0f*rjm*drijy*
						  (rij*rij+rjm*rjm-rim*rim))/(4.0f*rij*rij*rjm*rjm);

			real djcosthez = ((2.0f*rij*drijz-2.0f*rim*drimz)*
						  2.0f*rij*rjm -2.0f*rjm*drijz*
						  (rij*rij+rjm*rjm-rim*rim))/(4.0f*rij*rij*rjm*rjm);

			real djgthex = -C/((D+(H-jcosthe)*(H-jcosthe))*(D+(H-jcosthe)*(H-jcosthe)))* 2.0f*(H-jcosthe)*djcosthex;
			real djgthey = -C/((D+(H-jcosthe)*(H-jcosthe))*(D+(H-jcosthe)*(H-jcosthe)))* 2.0f*(H-jcosthe)*djcosthey;
			real djgthez = -C/((D+(H-jcosthe)*(H-jcosthe))*(D+(H-jcosthe)*(H-jcosthe)))* 2.0f*(H-jcosthe)*djcosthez;
			dxsjx += (fcj*(djgthex*exp(LAMDA3*(rij-re)) + LAMDA3*exp(LAMDA3*(rij-re))*drijx*jgthe));
			dxsjy += (fcj*(djgthey*exp(LAMDA3*(rij-re)) + LAMDA3*exp(LAMDA3*(rij-re))*drijy*jgthe));
			dxsjz += (fcj*(djgthez*exp(LAMDA3*(rij-re)) + LAMDA3*exp(LAMDA3*(rij-re))*drijz*jgthe));
		}
	}	

	real bij = 1.0f/(1.0f + xsi[i*NumberOfParticles+j]);

	real dbijx = (1.0f/((1.0f + xsi[i*NumberOfParticles+j])*(1.0f + xsi[i*NumberOfParticles+j])))*dxsix;
	real dbijy = (1.0f/((1.0f + xsi[i*NumberOfParticles+j])*(1.0f + xsi[i*NumberOfParticles+j])))*dxsiy;
	real dbijz = (1.0f/((1.0f + xsi[i*NumberOfParticles+j])*(1.0f + xsi[i*NumberOfParticles+j])))*dxsiz;

	real bji = 1.0f/(1.0f + xsj[i*NumberOfParticles+j]);

	real dbjix = (1.0f/((1.0f + xsj[i*NumberOfParticles+j])*(1.0f + xsj[i*NumberOfParticles+j])))*dxsjx;
	real dbjiy = (1.0f/((1.0f + xsj[i*NumberOfParticles+j])*(1.0f + xsj[i*NumberOfParticles+j])))*dxsjy;
	real dbjiz = (1.0f/((1.0f + xsj[i*NumberOfParticles+j])*(1.0f + xsj[i*NumberOfParticles+j])))*dxsjz;

	bo = 0.5f*(bij + bji);
	dbo.x = 0.5f*(dbijx + dbjix);
	dbo.y = 0.5f*(dbijy + dbjiy);
	dbo.z = 0.5f*(dbijz + dbjiz);
}

///
/// Bond Order 2
///
void bondOrder2(real4 * posArray, real re, real rij, real rik, int NumberOfParticles, int nlargestsize, int nlistsize, int * list, int k, int i, int j, real3 & dbo, real rcutsq, real xsi, real xsj)
{
	dbo.x = dbo.y = dbo.z = 0;
	real3 jk;
	real rjk;
	
	jk.x = posArray[j].x - posArray[k].x;
	jk.y = posArray[j].y - posArray[k].y;
	jk.z = posArray[j].z - posArray[k].z;

	rjk = sqrt(jk.x*jk.x + jk.y*jk.y + jk.z*jk.z);

	real drjkx = (posArray[k].x - posArray[j].x)/rjk;
  	real drjky = (posArray[k].y - posArray[j].y)/rjk;
  	real drjkz = (posArray[k].z - posArray[j].z)/rjk;
  
  	real drikx = (posArray[k].x - posArray[i].x)/rik;
  	real driky = (posArray[k].y - posArray[i].y)/rik;
  	real drikz = (posArray[k].z - posArray[i].z)/rik; 	

	real fci;
	real dfci;
	real fcj;
	real dfcj;
	
	cutoff(rik, fci, dfci);
	cutoff(rjk, fcj, dfcj);

	if (fci>0)
	{
		real icosthe = (rij*rij+rik*rik-rjk*rjk)/(2.0f*rij*rik);
		real igthe = 1.0f +C/D -C/(D+(H-icosthe)*(H-icosthe));

		real idcosthex = ((2.0f*rik*drikx - 2.0f*rjk*drjkx)*2.0f*rij*rik - (2.0f*rij*drikx)*
	      (rij*rij+rik*rik-rjk*rjk))/(4.0f*rij*rij*rik*rik);

		real idcosthey = ((2.0f*rik*driky - 2.0f*rjk*drjky)*2.0f*rij*rik - (2.0f*rij*driky)*
	      (rij*rij+rik*rik-rjk*rjk))/(4.0f*rij*rij*rik*rik);

		real idcosthez = ((2.0f*rik*drikz - 2.0f*rjk*drjkz)*2.0f*rij*rik - (2.0f*rij*drikz)*
	      (rij*rij+rik*rik-rjk*rjk))/(4.0f*rij*rij*rik*rik);

		real idgthex = -C/((D+(H-icosthe)*(H-icosthe))*(D+(H-icosthe)*(H-icosthe)))*2.0f*(H-icosthe)*idcosthex;
		real idgthey = -C/((D+(H-icosthe)*(H-icosthe))*(D+(H-icosthe)*(H-icosthe)))*2.0f*(H-icosthe)*idcosthey;
		real idgthez = -C/((D+(H-icosthe)*(H-icosthe))*(D+(H-icosthe)*(H-icosthe)))*2.0f*(H-icosthe)*idcosthez;

		real idxsix = dfci*drikx*igthe*exp(LAMDA3*(rij - re)) + fci*idgthex*exp(LAMDA3*(rij - re));
		real idxsiy = dfci*driky*igthe*exp(LAMDA3*(rij - re)) + fci*idgthey*exp(LAMDA3*(rij - re));
  		real idxsiz = dfci*drikz*igthe*exp(LAMDA3*(rij - re)) + fci*idgthez*exp(LAMDA3*(rij - re));

		real dbijx = (1.0f/((1.0f + xsi)*(1.0f + xsi)))*idxsix;
		real dbijy = (1.0f/((1.0f + xsi)*(1.0f + xsi)))*idxsiy;
  		real dbijz = (1.0f/((1.0f + xsi)*(1.0f + xsi)))*idxsiz;

		dbo.x += __fmul_rn(0.5f, dbijx);
		dbo.y += __fmul_rn(0.5f, dbijy);
		dbo.z += __fmul_rn(0.5f, dbijz);
	}

	if(fcj>0)
	{
		real jcosthe = (rij*rij+rjk*rjk-rik*rik)/(2.0f*rij*rjk);
		real jgthe = 1.0f +C/D -C/(D+(H-jcosthe)*(H-jcosthe));

		real jdcosthex = ((2.0f*rjk*drjkx - 2.0f*rik*drikx)*2.0f*rij*rjk - (2.0f*rij*drjkx)*
		  (rij*rij+rjk*rjk-rik*rik))/(4.0f*rij*rij*rjk*rjk);

		real jdcosthey = ((2.0f*rjk*drjky - 2.0f*rik*driky)*2.0f*rij*rjk - (2.0f*rij*drjky)*
		  (rij*rij+rjk*rjk-rik*rik))/(4.0f*rij*rij*rjk*rjk);

		real jdcosthez = ((2.0f*rjk*drjkz - 2.0f*rik*drikz)*2.0f*rij*rjk - (2.0f*rij*drjkz)*
		  (rij*rij+rjk*rjk-rik*rik))/(4.0f*rij*rij*rjk*rjk);

		real jdgthex = -C/((D+(H-jcosthe)*(H-jcosthe))*(D+(H-jcosthe)*(H-jcosthe)))*2.0f*(H-jcosthe)*jdcosthex;
		real jdgthey = -C/((D+(H-jcosthe)*(H-jcosthe))*(D+(H-jcosthe)*(H-jcosthe)))*2.0f*(H-jcosthe)*jdcosthey;
		real jdgthez = -C/((D+(H-jcosthe)*(H-jcosthe))*(D+(H-jcosthe)*(H-jcosthe)))*2.0f*(H-jcosthe)*jdcosthez;

		real jdxsix = dfcj*drjkx*jgthe*exp(LAMDA3*(rij - re)) + fcj*jdgthex*exp(LAMDA3*(rij - re));
		real jdxsiy = dfcj*drjky*jgthe*exp(LAMDA3*(rij - re)) + fcj*jdgthey*exp(LAMDA3*(rij - re));
  		real jdxsiz = dfcj*drjkz*jgthe*exp(LAMDA3*(rij - re)) + fcj*jdgthez*exp(LAMDA3*(rij - re));

		real dbjix = (1.0f/((1.0f + xsj)*(1.0f + xsj)))*jdxsix;
		real dbjiy = (1.0f/((1.0f + xsj)*(1.0f + xsj)))*jdxsiy;
		real dbjiz = (1.0f/((1.0f + xsj)*(1.0f + xsj)))*jdxsiz;
	
		dbo.x += __fmul_rn(0.5f, dbjix);
		dbo.y += __fmul_rn(0.5f, dbjiy);
		dbo.z += __fmul_rn(0.5f, dbjiz);
	}
}

///
/// Bond Order Forces
///
void bondForces(real3 & force, real4 * posArray, real rij, int NumberOfParticles, int i, int j, int itype, int jtype, real3 ij, int nlargestsize, int nlistsize, int * list, real rcutsq, real * xsi, real * xsj)
{
	real va=0;
	real dva=0;
	real vr=0;
	real dvr=0;
	real vdw;
	real dvdw;

	real bo=0;
	real3 dbo={0.0f, 0.0f, 0.0f};

	real fcut, dfcut;
	cutoff(rij/10, fcut, dfcut);

	if (fcut > 0)
	{
		real de = getDE<real>(itype, jtype);
		real re = getRE<real>(itype, jtype);
		real beta = getbeta<real>(itype, jtype);

		attractive(va, dva, rij, de, beta , re);
		repulsive(vr, dvr, rij, de, beta, re);
		bondOrder1(posArray, re, rij, NumberOfParticles, nlargestsize, nlistsize, list, i, j, ij, bo, dbo, rcutsq, xsi, xsj);
	}

	vanderWaals(vdw, dvdw, rij, getR0V<real>(itype, jtype));
	real fij = dvr - bo*dva + dvdw;

	real fxij = fij*(ij.x)/rij - dbo.x*va;
	real fyij = fij*(ij.y)/rij - dbo.y*va;
	real fzij = fij*(ij.z)/rij - dbo.z*va;

	force.x += fxij;
	force.y += fyij;
	force.z += fzij;
}

///
/// Many Body Forces
///
void manyBodyForces(real4 * posArray, real3 * forceArray, int NumberOfParticles, int nlargestsize, int nlistsize, int * list, real rcutsq, bool ljexist)
{
	real3 force;
	real3 ij;
	real3 ik;
	real rij;
	
	real * xsi;
	real * xsj;

	TC(xsi = new real[NumberOfParticles*NumberOfParticles]);
	TC(xsj = new real[NumberOfParticles*NumberOfParticles]);

	for(int i=0; i<NumberOfParticles; i++)
	{
		int type = (int)posArray[i].w;

		force.x = force.y = force.z = 0.0;

		for(int j=0; j<nlargestsize && list[nlistsize*i+j]!=-1; j++)
		{
			int jp = list[nlistsize*i+j];
				
			ij.x = posArray[i].x - posArray[jp].x;
			ij.y = posArray[i].y - posArray[jp].y;
			ij.z = posArray[i].z - posArray[jp].z;

			rij = ij.x*ij.x + ij.y*ij.y + ij.z*ij.z;

			if (rij>=rcutsq || rij>CUTOFF_RIJSQ) //distance is too large... ignoring...
				continue;

			rij = sqrt(rij);
			
			int jtype = (int)posArray[jp].w;

			if (jtype != NG)
			{
				bondForces(force, posArray, rij, NumberOfParticles, i, jp, type, jtype, ij, nlargestsize, nlistsize, list, rcutsq, xsi, xsj);
			}
		}

		if (ljexist)
		{
			forceArray[i].x += force.x;
			forceArray[i].y += force.y;
			forceArray[i].z += force.z;
		}
		else
		{
			forceArray[i].x = force.x;
			forceArray[i].y = force.y;
			forceArray[i].z = force.z;
		}
	}

	for(int k=0; k<NumberOfParticles; k++)
	{
		for(int i=0; i<nlargestsize && list[nlistsize*k+i]!=-1; i++)
		{
			int ip = list[nlistsize*k+i];

			int itype = (int)posArray[ip].w;

			ik.x = posArray[k].x - posArray[ip].x;
			ik.y = posArray[k].y - posArray[ip].y;
			ik.z = posArray[k].z - posArray[ip].z;

			real rik = ik.x*ik.x + ik.y*ik.y + ik.z*ik.z;

			if (rik>=rcutsq)
				continue;

			rik = sqrt(rik);

			for(int j=0; j<nlargestsize && list[nlistsize*k+j]!=-1; j++)
			{
				int jp = list[nlistsize*k+j];

				int jtype = (int)posArray[jp].w;

				if (ip == jp)
					continue;

				ij.x = posArray[ip].x - posArray[jp].x;
				ij.y = posArray[ip].y - posArray[jp].y;
				ij.z = posArray[ip].z - posArray[jp].z;

				rij = ij.x*ij.x + ij.y*ij.y + ij.z*ij.z;

				if (rij>=rcutsq || rij>CUTOFF_RIJSQ) //distance is too large... ignoring...
					continue;

				rij = sqrt(rij);

				real de = getDE<real>(itype, jtype);
				real re = getRE<real>(itype, jtype);
				real beta = getbeta<real>(itype, jtype);

				real va;
				real dva;
				real3 dbo;

				attractive(va, dva, rij, de, beta, re); 
				bondOrder2(posArray, re, rij, rik, NumberOfParticles, nlargestsize, nlistsize, list, k, ip, jp, dbo, rcutsq, xsi[ip*NumberOfParticles+jp], xsj[ip*NumberOfParticles+jp]);
				forceArray[k].x -= 0.5f*dbo.x*va;
	    		forceArray[k].y -= 0.5f*dbo.y*va;
	    		forceArray[k].z -= 0.5f*dbo.z*va;
			}
		}
	}

	delete [] xsi;
	delete [] xsj;
}