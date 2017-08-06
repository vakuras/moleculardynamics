///
/// CUDA Device Getter Functions
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

extern "C" __device__ float getMass(int id)
{
	if (id == N)
		return mN;
	else if (id == O)
		return mO;
	else
		return mNG;
}

extern "C" __device__ float getRE(int ti, int tj)
{
	if (ti==tj)
	{
		if (ti==N)
			return reN;
		else
			return reO;
	}
	else
		return reNO;
}

extern "C" __device__ float getDE(int ti, int tj)
{
	if (ti==tj)
	{
		if (ti==N)
			return deN;
		else
			return deO;
	}
	else
		return deNO;
}

extern "C" __device__ float getbeta(int ti, int tj)
{
	if (ti==tj)
	{
		if (ti==N)
			return betaN;
		else
			return betaO;
	}
	else
		return betaNO;
}

extern "C" __device__ float getR0V(int ti, int tj)
{
	if (ti==tj)
	{
		if (ti==N)
			return reNvdw;
		else
			return reOvdw;
	}
	else
		return reNOvdw;
}