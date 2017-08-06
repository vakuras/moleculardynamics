///
/// CUDA Device Bond Order Functions
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#define xsi			xs.x
#define xsj			xs.y

extern "C" __device__ void bondOrder1(float4 * posArray, float4 ipos, float4 jpos, float re, float rij, float3 ij, int i, int j, float & bo, float3 & dbo, int nlargestsize, float rcutsq, float2 * memfloat, int NumberOfParticles)
{
	float drijx = (ij.x)/rij;
	float drijy = (ij.y)/rij;
	float drijz = (ij.z)/rij;

	float2 xs = {0.0f,0.0f};
	float dxsix = 0.0f;
	float dxsiy = 0.0f;
	float dxsiz = 0.0f;

	float dxsjx = 0.0f;
	float dxsjy = 0.0f;
	float dxsjz = 0.0f;

	float3 im;
	float3 jm;
	float rim;
	float rjm;
	float4 mpos;

	for(int m=0; m<nlargestsize; m++) //for all the particles in the list of i
	{
		int mp = tex2D(mbTexRef, m, i); //get index from texture (m in the list of i)

		if (mp<0) //stop on dead index
			break;

		if (mp==j) //ignore same particle
			continue;

		mpos = posArray[mp]; //get position from global memory

		im.x = ipos.x - mpos.x;
		im.y = ipos.y - mpos.y;
		im.z = ipos.z - mpos.z;

		rim = __fmul_rn(im.x,im.x) + __fmul_rn(im.y,im.y) + __fmul_rn(im.z,im.z); //calculate distance

		if (rim>=rcutsq) //distance is bigger than rcut... ignoring...
			continue;

		rim = __fsqrt_rn(rim); //calculate distance finally

		jm.x = jpos.x - mpos.x;
		jm.y = jpos.y - mpos.y;
		jm.z = jpos.z - mpos.z;

		rjm = __fsqrt_rn(__fmul_rn(jm.x,jm.x) + __fmul_rn(jm.y,jm.y) + __fmul_rn(jm.z,jm.z)); //same for j-m

		float fci;
		float dfci;
		float fcj;
		float dfcj;
		
		cutoff(rim, fci, dfci);
		cutoff(rjm, fcj, dfcj);

		float drimx = (ipos.x - mpos.x)/rim;
		float drimy = (ipos.y - mpos.y)/rim;
		float drimz = (ipos.z - mpos.z)/rim;     

		if (fci>0)
		{
			float icosthe = (__fmul_rn(rij,rij)+__fmul_rn(rim,rim)-__fmul_rn(rjm,rjm))/ (__fmul_rn(__fmul_rn(2.0f,rij),rim));
			float igthe = 1.0f +C/D -C/(D+__fmul_rn(H-icosthe,H-icosthe));
			xsi += __fmul_rn(__fmul_rn(fci,igthe),exp(__fmul_rn(LAMDA3,rij - re)));

			float dicosthex = (__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rij),drijx)+__fmul_rn(__fmul_rn(2.0f,rim),drimx),__fmul_rn(2.0f,__fmul_rn(rij,rim)))
					- __fmul_rn(__fmul_rn(__fmul_rn(2.0f,rim),drijx) + __fmul_rn(__fmul_rn(2.0f,rij),drimx),
					(__fmul_rn(rij,rij)+__fmul_rn(rim,rim)-__fmul_rn(rjm,rjm))))/__fmul_rn(4.0f,__fmul_rn(rij,__fmul_rn(rij,__fmul_rn(rim,rim))));
					
			float dicosthey = (__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rij),drijy)+__fmul_rn(__fmul_rn(2.0f,rim),drimy),__fmul_rn(2.0f,__fmul_rn(rij,rim)))
					- __fmul_rn(__fmul_rn(__fmul_rn(2.0f,rim),drijy) + __fmul_rn(__fmul_rn(2.0f,rij),drimy),
					(__fmul_rn(rij,rij)+__fmul_rn(rim,rim)-__fmul_rn(rjm,rjm))))/__fmul_rn(4.0f,__fmul_rn(rij,__fmul_rn(rij,__fmul_rn(rim,rim))));
					
			float dicosthez = (__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rij),drijz)+__fmul_rn(__fmul_rn(2.0f,rim),drimz),__fmul_rn(2.0f,__fmul_rn(rij,rim)))
					- __fmul_rn(__fmul_rn(__fmul_rn(2.0f,rim),drijz) + __fmul_rn(__fmul_rn(2.0f,rij),drimz),
					(__fmul_rn(rij,rij)+__fmul_rn(rim,rim)-__fmul_rn(rjm,rjm))))/__fmul_rn(4.0f,__fmul_rn(rij,__fmul_rn(rij,__fmul_rn(rim,rim))));

			float digthex = __fmul_rn(-C/(__fmul_rn((D+__fmul_rn(H-icosthe,H-icosthe)),(D+__fmul_rn(H-icosthe,H-icosthe)))),__fmul_rn(2.0f,__fmul_rn(H-icosthe,dicosthex)));
			float digthey = __fmul_rn(-C/(__fmul_rn((D+__fmul_rn(H-icosthe,H-icosthe)),(D+__fmul_rn(H-icosthe,H-icosthe)))),__fmul_rn(2.0f,__fmul_rn(H-icosthe,dicosthey)));
			float digthez = __fmul_rn(-C/(__fmul_rn((D+__fmul_rn(H-icosthe,H-icosthe)),(D+__fmul_rn(H-icosthe,H-icosthe)))),__fmul_rn(2.0f,__fmul_rn(H-icosthe,dicosthez)));

			dxsix += __fmul_rn(__fmul_rn(__fmul_rn(dfci,drimx),igthe),exp(__fmul_rn(LAMDA3,rij - re))) + __fmul_rn(fci,__fmul_rn(digthex,exp(__fmul_rn(LAMDA3,rij - re)))) + __fmul_rn(fci,__fmul_rn(igthe,__fmul_rn(LAMDA3,__fmul_rn(exp(__fmul_rn(LAMDA3,rij - re)),drijx))));
			dxsiy += __fmul_rn(__fmul_rn(__fmul_rn(dfci,drimy),igthe),exp(__fmul_rn(LAMDA3,rij - re))) + __fmul_rn(fci,__fmul_rn(digthey,exp(__fmul_rn(LAMDA3,rij - re)))) + __fmul_rn(fci,__fmul_rn(igthe,__fmul_rn(LAMDA3,__fmul_rn(exp(__fmul_rn(LAMDA3,rij - re)),drijy))));
			dxsiz += __fmul_rn(__fmul_rn(__fmul_rn(dfci,drimz),igthe),exp(__fmul_rn(LAMDA3,rij - re))) + __fmul_rn(fci,__fmul_rn(digthez,exp(__fmul_rn(LAMDA3,rij - re)))) + __fmul_rn(fci,__fmul_rn(igthe,__fmul_rn(LAMDA3,__fmul_rn(exp(__fmul_rn(LAMDA3,rij - re)),drijz))));
		}

		if (fcj>0)
		{
			float jcosthe = (__fmul_rn(rij,rij)+__fmul_rn(rjm,rjm)-__fmul_rn(rim,rim))/ (__fmul_rn(__fmul_rn(2.0f,rij),rjm));
			float jgthe = 1.0f +C/D -C/(D+__fmul_rn(H-jcosthe,H-jcosthe));
			xsj += __fmul_rn(__fmul_rn(fcj,jgthe),exp(__fmul_rn(LAMDA3,rij - re)));

			float djcosthex = (__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rij),drijx)-__fmul_rn(__fmul_rn(2.0f,rim),drimx),
						  __fmul_rn(__fmul_rn(2.0f,rij),rjm)) -__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rjm),drijx),
						  (__fmul_rn(rij,rij)+__fmul_rn(rjm,rjm)-__fmul_rn(rim,rim))))/__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(4.0f,rij),rij),rjm),rjm);

			float djcosthey = (__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rij),drijy)-__fmul_rn(__fmul_rn(2.0f,rim),drimy),
						  __fmul_rn(__fmul_rn(2.0f,rij),rjm)) -__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rjm),drijy),
						  (__fmul_rn(rij,rij)+__fmul_rn(rjm,rjm)-__fmul_rn(rim,rim))))/__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(4.0f,rij),rij),rjm),rjm);
						  
			float djcosthez = (__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rij),drijz)-__fmul_rn(__fmul_rn(2.0f,rim),drimz),
						  __fmul_rn(__fmul_rn(2.0f,rij),rjm)) -__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rjm),drijz),
						  (__fmul_rn(rij,rij)+__fmul_rn(rjm,rjm)-__fmul_rn(rim,rim))))/__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(4.0f,rij),rij),rjm),rjm);
			
			float djgthex = __fmul_rn(-C/__fmul_rn(D+__fmul_rn(H-jcosthe,H-jcosthe),D+__fmul_rn(H-jcosthe,H-jcosthe)),__fmul_rn(__fmul_rn(2.0f,H-jcosthe),djcosthex));
			float djgthey = __fmul_rn(-C/__fmul_rn(D+__fmul_rn(H-jcosthe,H-jcosthe),D+__fmul_rn(H-jcosthe,H-jcosthe)),__fmul_rn(__fmul_rn(2.0f,H-jcosthe),djcosthey));
			float djgthez = __fmul_rn(-C/__fmul_rn(D+__fmul_rn(H-jcosthe,H-jcosthe),D+__fmul_rn(H-jcosthe,H-jcosthe)),__fmul_rn(__fmul_rn(2.0f,H-jcosthe),djcosthez));

			
			dxsjx += __fmul_rn(fcj,__fmul_rn(djgthex,exp(__fmul_rn(LAMDA3,rij-re))) + __fmul_rn(__fmul_rn(__fmul_rn(LAMDA3,exp(__fmul_rn(LAMDA3,rij-re))),drijx),jgthe));
			dxsjy += __fmul_rn(fcj,__fmul_rn(djgthey,exp(__fmul_rn(LAMDA3,rij-re))) + __fmul_rn(__fmul_rn(__fmul_rn(LAMDA3,exp(__fmul_rn(LAMDA3,rij-re))),drijy),jgthe));
			dxsjz += __fmul_rn(fcj,__fmul_rn(djgthez,exp(__fmul_rn(LAMDA3,rij-re))) + __fmul_rn(__fmul_rn(__fmul_rn(LAMDA3,exp(__fmul_rn(LAMDA3,rij-re))),drijz),jgthe));
		}
	}	

	float bij = 1.0f/(1.0f + xsi);

	float dbijx = __fmul_rn(1.0f/(__fmul_rn(1.0f + xsi,1.0f + xsi)),dxsix);
	float dbijy = __fmul_rn(1.0f/(__fmul_rn(1.0f + xsi,1.0f + xsi)),dxsiy);
	float dbijz = __fmul_rn(1.0f/(__fmul_rn(1.0f + xsi,1.0f + xsi)),dxsiz);

	float bji = 1.0f/(1.0f + xsj);

	float dbjix = __fmul_rn(1.0f/(__fmul_rn(1.0f + xsj,1.0f + xsj)),dxsjx);
	float dbjiy = __fmul_rn(1.0f/(__fmul_rn(1.0f + xsj,1.0f + xsj)),dxsjy);
	float dbjiz = __fmul_rn(1.0f/(__fmul_rn(1.0f + xsj,1.0f + xsj)),dxsjz);


	memfloat[i*NumberOfParticles+j] = xs; //save to global memory for use in bond order 2

	bo = __fmul_rn(0.5f,bij + bji); //finally calculate bo-1
	dbo.x = __fmul_rn(0.5f,dbijx + dbjix);
	dbo.y = __fmul_rn(0.5f,dbijy + dbjiy);
	dbo.z = __fmul_rn(0.5f,dbijz + dbjiz);
}

///
/// Bond Order 2
///
extern "C" __device__ void bondOrder2(float4 * posArray, float4 ipos, float4 jpos, float4 kpos, float re, float rij, float rik, int k, int i, int j, float3 & dbo, float rcutsq, float2 * memfloat, int NumberOfParticles)
{
	dbo.x = dbo.y = dbo.z = 0;

	float3 jk;
	float rjk;
	float2 xs = memfloat[i*NumberOfParticles+j]; //read from global memory the calculation made in b-o-1

	jk.x = jpos.x - kpos.x;
	jk.y = jpos.y - kpos.y;
	jk.z = jpos.z - kpos.z;

	rjk = __fsqrt_rn(__fmul_rn(jk.x,jk.x) + __fmul_rn(jk.y,jk.y) + __fmul_rn(jk.z,jk.z)); //calculate distance

	float drjkx = (kpos.x - jpos.x)/rjk;
  	float drjky = (kpos.y - jpos.y)/rjk;
  	float drjkz = (kpos.z - jpos.z)/rjk;
  
  	float drikx = (kpos.x - posArray[i].x)/rik;
  	float driky = (kpos.y - posArray[i].y)/rik;
  	float drikz = (kpos.z - posArray[i].z)/rik; 

	float fci;
	float dfci;
	float fcj;
	float dfcj;
	
	cutoff(rik, fci, dfci);
	cutoff(rjk, fcj, dfcj);

	if (fci>0)
	{
		float icosthe = (__fmul_rn(rij,rij)+__fmul_rn(rik,rik)-__fmul_rn(rjk,rjk))/(__fmul_rn(__fmul_rn(2.0f,rij),rik));
		float igthe = 1.0f +C/D -C/(D+__fmul_rn(H-icosthe,H-icosthe));

		float idcosthex = (__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rik),drikx) - __fmul_rn(__fmul_rn(2.0f,rjk),drjkx),2.0f),rij),rik) - 
						__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rij),drikx),__fmul_rn(rij,rij)+__fmul_rn(rik,rik)-__fmul_rn(rjk,rjk)))/(__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(4.0f,rij),rij),rik),rik));

		float idcosthey = (__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rik),driky) - __fmul_rn(__fmul_rn(2.0f,rjk),drjky),2.0f),rij),rik) - 
						__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rij),driky),__fmul_rn(rij,rij)+__fmul_rn(rik,rik)-__fmul_rn(rjk,rjk)))/(__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(4.0f,rij),rij),rik),rik));

		float idcosthez = (__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rik),drikz) - __fmul_rn(__fmul_rn(2.0f,rjk),drjkz),2.0f),rij),rik) - 
						__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rij),drikz),__fmul_rn(rij,rij)+__fmul_rn(rik,rik)-__fmul_rn(rjk,rjk)))/(__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(4.0f,rij),rij),rik),rik));

		float idgthex = __fmul_rn(__fmul_rn(__fmul_rn(-C/(__fmul_rn(D+__fmul_rn(H-icosthe,H-icosthe),D+__fmul_rn(H-icosthe,H-icosthe))),2.0f),H-icosthe),idcosthex);
		float idgthey = __fmul_rn(__fmul_rn(__fmul_rn(-C/(__fmul_rn(D+__fmul_rn(H-icosthe,H-icosthe),D+__fmul_rn(H-icosthe,H-icosthe))),2.0f),H-icosthe),idcosthey);
		float idgthez = __fmul_rn(__fmul_rn(__fmul_rn(-C/(__fmul_rn(D+__fmul_rn(H-icosthe,H-icosthe),D+__fmul_rn(H-icosthe,H-icosthe))),2.0f),H-icosthe),idcosthez);

		float idxsix = __fmul_rn(__fmul_rn(__fmul_rn(dfci,drikx),igthe),exp(__fmul_rn(LAMDA3,rij - re))) + __fmul_rn(__fmul_rn(fci,idgthex),exp(__fmul_rn(LAMDA3,rij - re)));
		float idxsiy = __fmul_rn(__fmul_rn(__fmul_rn(dfci,driky),igthe),exp(__fmul_rn(LAMDA3,rij - re))) + __fmul_rn(__fmul_rn(fci,idgthey),exp(__fmul_rn(LAMDA3,rij - re)));
		float idxsiz = __fmul_rn(__fmul_rn(__fmul_rn(dfci,drikz),igthe),exp(__fmul_rn(LAMDA3,rij - re))) + __fmul_rn(__fmul_rn(fci,idgthez),exp(__fmul_rn(LAMDA3,rij - re)));

		float dbijx = __fmul_rn((1.0f/(__fmul_rn(1.0f + xsi,1.0f + xsi))),idxsix);
		float dbijy = __fmul_rn((1.0f/(__fmul_rn(1.0f + xsi,1.0f + xsi))),idxsiy);
		float dbijz = __fmul_rn((1.0f/(__fmul_rn(1.0f + xsi,1.0f + xsi))),idxsiz);

		dbo.x += __fmul_rn(0.5f, dbijx);
		dbo.y += __fmul_rn(0.5f, dbijy);
		dbo.z += __fmul_rn(0.5f, dbijz);
	}

	if (fcj>0)
	{
		float jcosthe = (__fmul_rn(rij,rij)+__fmul_rn(rjk,rjk)-__fmul_rn(rik,rik))/(__fmul_rn(__fmul_rn(2.0f,rij),rjk));
		float jgthe = 1.0f +C/D -C/(D+__fmul_rn(H-jcosthe,H-jcosthe));

		float jdcosthex = (__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rjk),drjkx) - __fmul_rn(__fmul_rn(2.0f,rik),drikx),2.0f),rij),rjk) - 
						__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rij),drjkx),__fmul_rn(rij,rij)+__fmul_rn(rjk,rjk)-__fmul_rn(rik,rik)))/(__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(4.0f,rij),rij),rjk),rjk));

		float jdcosthey = (__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rjk),drjky) - __fmul_rn(__fmul_rn(2.0f,rik),driky),2.0f),rij),rjk) - 
						__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rij),drjky),__fmul_rn(rij,rij)+__fmul_rn(rjk,rjk)-__fmul_rn(rik,rik)))/(__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(4.0f,rij),rij),rjk),rjk));

		float jdcosthez = (__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rjk),drjkz) - __fmul_rn(__fmul_rn(2.0f,rik),drikz),2.0f),rij),rjk) - 
						__fmul_rn(__fmul_rn(__fmul_rn(2.0f,rij),drjkz),__fmul_rn(rij,rij)+__fmul_rn(rjk,rjk)-__fmul_rn(rik,rik)))/(__fmul_rn(__fmul_rn(__fmul_rn(__fmul_rn(4.0f,rij),rij),rjk),rjk));

		
		float jdgthex = __fmul_rn(__fmul_rn(__fmul_rn(-C/(__fmul_rn(D+__fmul_rn(H-jcosthe,H-jcosthe),D+__fmul_rn(H-jcosthe,H-jcosthe))),2.0f),H-jcosthe),jdcosthex);
		float jdgthey = __fmul_rn(__fmul_rn(__fmul_rn(-C/(__fmul_rn(D+__fmul_rn(H-jcosthe,H-jcosthe),D+__fmul_rn(H-jcosthe,H-jcosthe))),2.0f),H-jcosthe),jdcosthey);
		float jdgthez = __fmul_rn(__fmul_rn(__fmul_rn(-C/(__fmul_rn(D+__fmul_rn(H-jcosthe,H-jcosthe),D+__fmul_rn(H-jcosthe,H-jcosthe))),2.0f),H-jcosthe),jdcosthez);

		
		float jdxsix = __fmul_rn(__fmul_rn(__fmul_rn(dfcj,drjkx),jgthe),exp(__fmul_rn(LAMDA3,rij - re))) + __fmul_rn(__fmul_rn(fcj,jdgthex),exp(__fmul_rn(LAMDA3,rij - re)));
		float jdxsiy = __fmul_rn(__fmul_rn(__fmul_rn(dfcj,drjky),jgthe),exp(__fmul_rn(LAMDA3,rij - re))) + __fmul_rn(__fmul_rn(fcj,jdgthey),exp(__fmul_rn(LAMDA3,rij - re)));
		float jdxsiz = __fmul_rn(__fmul_rn(__fmul_rn(dfcj,drjkz),jgthe),exp(__fmul_rn(LAMDA3,rij - re))) + __fmul_rn(__fmul_rn(fcj,jdgthez),exp(__fmul_rn(LAMDA3,rij - re)));


		float dbjix = __fmul_rn((1.0f/(__fmul_rn(1.0f + xsj,1.0f + xsj))),jdxsix);
  		float dbjiy = __fmul_rn((1.0f/(__fmul_rn(1.0f + xsj,1.0f + xsj))),jdxsiy);
  		float dbjiz = __fmul_rn((1.0f/(__fmul_rn(1.0f + xsj,1.0f + xsj))),jdxsiz);

		dbo.x += __fmul_rn(0.5f, dbjix);
		dbo.y += __fmul_rn(0.5f, dbjiy);
		dbo.z += __fmul_rn(0.5f, dbjiz);
	}
}

///
/// Bond Forces
///
extern "C" __device__ void bondForces(float3 & force, float4 * posArray, float4 ipos, float4 jpos, float rij, int i, int jp, int itype, int jtype, float3 ij, int nlargestsize, float rcutsq, float2 * memfloat, int NumberOfParticles)
{
	float va=0;
	float dva=0;
	float vr=0;
	float dvr=0;
	float vdw=0;
	float dvdw=0;

	float bo=0;
	float3 dbo={0,0,0};

	float fcut, dfcut;
	cutoff(rij/10, fcut, dfcut);

	if (fcut > 0) //only if cutoff is > 0
	{
		float de = getDE(itype, jtype);
		float re = getRE(itype, jtype);
		float beta = getbeta(itype, jtype);

		attractive(va, dva, rij, de, beta , re);
		repulsive(vr, dvr, rij, de, beta, re);
		bondOrder1(posArray, ipos, jpos, re, rij, ij, i, jp, bo, dbo, nlargestsize, rcutsq, memfloat, NumberOfParticles);
	}

	vanderWaals(vdw, dvdw, rij, getR0V(itype, jtype)); //do always vdw
	float fij = dvr - __fmul_rn(bo,dva) + dvdw;

	float fxij = __fmul_rn(fij,ij.x)/rij - __fmul_rn(dbo.x,va);
	float fyij = __fmul_rn(fij,ij.y)/rij - __fmul_rn(dbo.y,va);
	float fzij = __fmul_rn(fij,ij.z)/rij - __fmul_rn(dbo.z,va);

	force.x += fxij;
	force.y += fyij;
	force.z += fzij;
}

///
/// Many Body Forces 1 Kernel
///
extern "C" __global__ void manyBodyForces1(float4 * posArray, float3 * forceArray, int nlargestsize, int NumberOfParticles, float rcutsq, float2 * devXS, int ljexist)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	if (id>=NumberOfParticles) //more thread than particles
		return;

	float3 force;
	float4 idpos = posArray[id];

	if (ljexist==0) //zero force if not using lennard jones in the simulation
		force.x = force.y = force.z = 0.0f;
	else //if using lennard jones in the simulation, then read it's result
		force = forceArray[id];

	float3 ij;
	float4 jpos;

	for(int j=0; j<nlargestsize; j++)
	{
		int jp = tex2D(mbTexRef, j, id); //get index from texture

		if (jp<0) //stop on dead index
			break;
		
		jpos = posArray[jp];

		ij.x = idpos.x - jpos.x;
		ij.y = idpos.y - jpos.y;
		ij.z = idpos.z - jpos.z;

		float rij = __fmul_rn(ij.x,ij.x) + __fmul_rn(ij.y,ij.y) + __fmul_rn(ij.z,ij.z);

		if (rij>=rcutsq || rij>CUTOFF_RIJSQ) //distance is too large or bigger than rcut... ignoring...
			continue;

		rij = __fsqrt_rn(rij);

		bondForces(force, posArray, idpos, jpos, rij, id, jp, (int)idpos.w, (int)jpos.w, ij, nlargestsize, rcutsq, devXS, NumberOfParticles);
	}

	forceArray[id] = force;
}

///
/// Many Body Forces 2 Kernel
///
extern "C" __global__ void manyBodyForces2(float4 * posArray, float3 * forceArray, int nlargestsize, int NumberOfParticles, float rcutsq, float2 * devXS)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	if (id>=NumberOfParticles) //more thread than particles
		return;

	float3 force = forceArray[id];
	float4 idpos = posArray[id];

	float3 ij;
	float4 jpos;
	float3 ik;
	float4 ipos;

	for(int i=0; i<nlargestsize; i++)
	{
		int ip = tex2D(mbTexRef, i, id); //get index from texture

		if (ip<0) //stop on dead index
			break;

		ipos = posArray[ip];

		ik.x = ipos.x - idpos.x;
		ik.y = ipos.y - idpos.y;
		ik.z = ipos.z - idpos.z;

		float rik = __fmul_rn(ik.x,ik.x) + __fmul_rn(ik.y,ik.y) + __fmul_rn(ik.z,ik.z);

		if (rik>=rcutsq || rik>CUTOFF_RIJSQ) //distance is too large or bigger than rcut... ignoring...
			continue;

		rik = __fsqrt_rn(rik);

		for(int j=0; j<nlargestsize; j++)
		{
			int jp = tex2D(mbTexRef, j, id); //get index from texture


			if (jp<0) //stop on dead index
				break;

			if (jp == ip) //same particle -> ignore
				continue;


			jpos = posArray[jp];

			ij.x = ipos.x - jpos.x;
			ij.y = ipos.y - jpos.y;
			ij.z = ipos.z - jpos.z;

			float rij = __fmul_rn(ij.x,ij.x) + __fmul_rn(ij.y,ij.y) + __fmul_rn(ij.z,ij.z);

			if (rij>=rcutsq || rij>CUTOFF_RIJSQ) //distance is too large... ignoring...
				continue;

			rij = __fsqrt_rn(rij);

			float va;
			float dva;
			float3 dbo;

			float de = getDE(ipos.w, jpos.w);
			float re = getRE(ipos.w, jpos.w);
			float beta = getbeta(ipos.w, jpos.w);

			attractive(va, dva, rij, de, beta, re); 
			bondOrder2(posArray, ipos, jpos, idpos, re, rij, rik, id, ip, jp, dbo, rcutsq, devXS, NumberOfParticles);

			force.x -= __fmul_rn(__fmul_rn(0.5f,dbo.x),va);
			force.y -= __fmul_rn(__fmul_rn(0.5f,dbo.y),va);
			force.z -= __fmul_rn(__fmul_rn(0.5f,dbo.z),va);
		}
	}

	forceArray[id] = force;
}


///
/// Many Body Forces Kernel 1 BPP Version
///
extern "C" __global__ void manyBodyForcesBPP1(float4 * posArray, float3 * forceArray, float rcutsq, int NumberOfParticles, float2 * devXS, bool ljexist)
{
	extern __shared__ float4 shforce[];

	float3 force = {0.0f, 0.0f, 0.0f};
	float3 bt;
	float4 tpos;
	float4 bpos;

	int tp = tex2D(mbTexRef, threadIdx.x, blockIdx.x); //get index from texture

	if (tp > -1) //do only if the index is a dead one
	{				//not exiting on tp<0 cause if there is not particle we still want to set the force to 0
		tpos = posArray[tp];
		bpos = posArray[blockIdx.x];

		bt.x = bpos.x - tpos.x;
		bt.y = bpos.y - tpos.y;
		bt.z = bpos.z - tpos.z;

		float rbt = __fmul_rn(bt.x,bt.x) + __fmul_rn(bt.y,bt.y) + __fmul_rn(bt.z,bt.z);

		if (rbt<rcutsq && rbt<CUTOFF_RIJSQ) //distance is too large or bigger than rcut... ignoring...
		{
			rbt = __fsqrt_rn(rbt);
			bondForces(force, posArray, bpos, tpos, rbt, blockIdx.x, tp, (int)bpos.w, (int)tpos.w, bt, blockDim.x, rcutsq, devXS, NumberOfParticles);
		}
	}

	shforce[threadIdx.x].x = force.x;
	shforce[threadIdx.x].y = force.y;
	shforce[threadIdx.x].z = force.z;

	__syncthreads();  //wait for all thread to finish

	if (threadIdx.x==0) //first thread makes the accumulation
	{
		if (ljexist==0) //zero force if not using lennard jones in the simulation
			force.x = force.y = force.z = 0.0f;
		else //if using lennard jones in the simulation, then read it's result
			force = forceArray[blockIdx.x];

		float4 tmp;

		for(int i=0; i<blockDim.x; i++) //accumulate
		{
			tmp = shforce[i];
			force.x += tmp.x;
			force.y += tmp.y;
			force.z += tmp.z;
		}

		forceArray[blockIdx.x] = force; //save force result to memory
	}
}

///
/// Many Body Forces Kernel 2 BPP Version
///
extern "C" __global__ void manyBodyForcesBPP2(float4 * posArray, float3 * forceArray, float rcutsq, int NumberOfParticles, float2 * devXS)
{
	extern __shared__ float4 shforce[];

	float3 force = {0.0f, 0.0f, 0.0f};
	float3 bt;
	float3 tj;
	float4 tpos;
	float4 bpos;

	int tp = tex2D(mbTexRef, threadIdx.x, blockIdx.x); //get index from texture

	if (tp > -1) //do only if the index is a dead one
	{				//not exiting on tp<0 cause if there is not particle we still want to set the force to 0
		tpos = posArray[tp];
		bpos = posArray[blockIdx.x];

		bt.x = tpos.x - bpos.x;
		bt.y = tpos.y - bpos.y;
		bt.z = tpos.z - bpos.z;

		float rbt = __fmul_rn(bt.x,bt.x) + __fmul_rn(bt.y,bt.y) + __fmul_rn(bt.z,bt.z);

		if (rbt<rcutsq && rbt<CUTOFF_RIJSQ) //distance is too large or bigger than rcut... ignoring...
		{
			rbt = __fsqrt_rn(rbt);

			for(int j=0; j<blockDim.x; j++)
			{
				int jp = tex2D(mbTexRef, j, blockIdx.x); //get index from texture

				if (jp<0) //stop on dead index
					break;

				if (jp == tp) //same particle -> ignore
					continue;

				float4 jpos = posArray[jp];

				tj.x = tpos.x - jpos.x;
				tj.y = tpos.y - jpos.y;
				tj.z = tpos.z - jpos.z;

				float rtj = __fmul_rn(tj.x,tj.x) + __fmul_rn(tj.y,tj.y) + __fmul_rn(tj.z,tj.z);

				if (rtj>=rcutsq || rtj>CUTOFF_RIJSQ) //distance is too large or bigger than rcut... ignoring...
					continue;

				rtj = __fsqrt_rn(rtj);

				float va;
				float dva;
				float3 dbo;

				float de = getDE(tpos.w, jpos.w);
				float re = getRE(tpos.w, jpos.w);
				float beta = getbeta(tpos.w, jpos.w);

				attractive(va, dva, rtj, de, beta, re); 
				bondOrder2(posArray, tpos, jpos, bpos, re, rtj, rbt, blockIdx.x, tp, jp, dbo, rcutsq, devXS, NumberOfParticles);

				force.x -= __fmul_rn(__fmul_rn(0.5f,dbo.x),va);
				force.y -= __fmul_rn(__fmul_rn(0.5f,dbo.y),va);
				force.z -= __fmul_rn(__fmul_rn(0.5f,dbo.z),va);
			}
		}
	}

	shforce[threadIdx.x].x = force.x;
	shforce[threadIdx.x].y = force.y;
	shforce[threadIdx.x].z = force.z;

	__syncthreads(); //wait for all thread to finish their work

	if (threadIdx.x==0) //first thread makes the accumulation
	{
		force = forceArray[blockIdx.x];
		float4 tmp;

		for(int i=0; i<blockDim.x; i++) //accumulate
		{
			tmp = shforce[i];
			force.x += tmp.x;
			force.y += tmp.y;
			force.z += tmp.z;
		}

		forceArray[blockIdx.x] = force; //save force result to memory
	}
}

