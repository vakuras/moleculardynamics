///
/// CUDA Device Lennard Jones Function
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

///
/// Lennard-Jones potentional
///
extern "C" __device__ void lennardJones(float3 & force, float rij, float3 ij, float sig, float eps)
{
	float sigsq = __fmul_rn(sig,sig);
	float con = __fmul_rn(24.0f,eps)/sigsq;
    float dist2 = __fmul_rn(rij,rij);
	dist2 = dist2 / sigsq;

    float dist4 = __fmul_rn(dist2,dist2);
    float dist8 = __fmul_rn(dist4,dist4);
    float dist14 = __fmul_rn(__fmul_rn(dist2,dist4),dist8);
    float invdist8= 1.0f/dist8;
    float invdist14= 1.0f/dist14;
    float s = __fmul_rn(2.0f,invdist14)-invdist8;
	float fij = __fmul_rn(s,con);

	force.x += __fmul_rn(ij.x,fij);
    force.y += __fmul_rn(ij.y,fij);
    force.z += __fmul_rn(ij.z,fij);
}

///
/// Lennard Jones Forces Kernel
///
extern "C" __global__ void lennardJonesForces(float4 * posArray, float3 * forceArray, int nlargestsize, int NumberOfParticles, float rcutsq)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	if (id>=NumberOfParticles) //more thread than particles
		return;

	float3 force = {0.0f, 0.0f, 0.0f}; 
	float3 ij;
	float4 ipos = posArray[id];
	float4 jpos;

	int idtype = (int)ipos.w;

	for(int j=0; j<nlargestsize; j++)
	{
		int jp = tex2D(ljTexRef, j, id); //get index from texture

		if (jp<0) //ignore dead index
			break;
		
		jpos = posArray[jp]; //read memory

		ij.x = ipos.x - jpos.x; //calculate distance
		ij.y = ipos.y - jpos.y;
		ij.z = ipos.z - jpos.z;

		float rij = __fmul_rn(ij.x,ij.x) + __fmul_rn(ij.y,ij.y) + __fmul_rn(ij.z,ij.z);

		if (rij>=rcutsq)
			continue;

		rij = __fsqrt_rn(rij);
		
		int jtype = (int)jpos.w;

		if (jtype == NG || idtype == NG)
		{
			int select = (idtype == NG) ? jtype : idtype; //select type

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

			lennardJones(force, rij, ij, sig, eps); 
		}
	}

	forceArray[id] = force;
}

///
/// Lennard Jones Forces Kernel BPP Version
///
extern "C" __global__ void lennardJonesForcesBPP(float4 * posArray, float3 * forceArray, float rcutsq)
{
	extern __shared__ float4 shforce[];

	float4 tmp = {0.0f, 0.0f, 0.0f, 0.0f};
	float3 force = {0.0f, 0.0f, 0.0f};
	float3 bt;
	float4 tpos;
	float4 bpos;

	shforce[threadIdx.x] = tmp; //clear shared memory index for thread

	__syncthreads(); //synchornize

	int tp = tex2D(ljTexRef, threadIdx.x, blockIdx.x); //get index from texture

	if (tp > -1) //do only if the index is a dead one
	{				//not exiting on tp<0 cause if there is not particle we still want to set the force to 0
		tpos = posArray[tp];
		bpos = posArray[blockIdx.x];

		bt.x = bpos.x - tpos.x;
		bt.y = bpos.y - tpos.y;
		bt.z = bpos.z - tpos.z;

		int btype = (int)bpos.w;
		int ttype = (int)tpos.w;

		float rbt = __fmul_rn(bt.x,bt.x) + __fmul_rn(bt.y,bt.y) + __fmul_rn(bt.z,bt.z);

		if (rbt<rcutsq && (btype == NG || ttype == NG))
		{
			int select = (btype == NG) ? ttype : btype; //select type

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

			float sigsq = __fmul_rn(sig,sig);
			float con = __fmul_rn(24.0f,eps)/sigsq;
			rbt = rbt / sigsq;

			float dist4 = __fmul_rn(rbt,rbt);
			float dist8 = __fmul_rn(dist4,dist4);
			float dist14 = __fmul_rn(__fmul_rn(rbt,dist4),dist8);
			float invdist8= 1.0f/dist8;
			float invdist14= 1.0f/dist14;
			float s = __fmul_rn(2.0f,invdist14)-invdist8;
			float fij = __fmul_rn(s,con);

			shforce[threadIdx.x].x = __fmul_rn(bt.x,fij);
			shforce[threadIdx.x].y = __fmul_rn(bt.y,fij);
			shforce[threadIdx.x].z = __fmul_rn(bt.z,fij);
		}
	}

	__syncthreads(); //wait for all thread to finish their work

	if (threadIdx.x==0) //the first thread makes the additions
	{
		force.x = force.y = force.z = 0;

		for(int i=0; i<blockDim.x; i++)
		{
			tmp = shforce[i]; //get result from shared memory
			force.x += tmp.x; //add
			force.y += tmp.y;
			force.z += tmp.z;
		}

		forceArray[blockIdx.x] = force; //save accumulated result to memory
	}
}