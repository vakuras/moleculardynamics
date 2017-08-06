///
/// CUDA Device Functions: Attractive, Repulsive, van der Waals & Cutoff
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

///
/// Attractive force
///
extern "C" __device__ void attractive(float & va, float & dva, float rij, float De, float beta, float re)
{
	va = __fmul_rn(__fmul_rn(S,De)/(S-1.0f),exp(__fmul_rn(-beta,__fmul_rn(sqrt(2.0f/S),rij - re))));
	dva = __fmul_rn(beta,__fmul_rn(sqrt(2.0f/S),va));
}

///
/// Repulsive force
///
extern "C" __device__ void repulsive(float & vr, float & dvr, float rij, float De, float beta, float re)
{
	vr = __fmul_rn(De/(S-1.0f),exp(__fmul_rn(-beta,__fmul_rn(sqrt(S*2.0f),rij-re)))) + __fmul_rn(7.0e-8f,exp(__fmul_rn(-10.0f,__fmul_rn(beta,rij-re))));
	dvr = __fmul_rn(beta,__fmul_rn(sqrt(__fmul_rn(2.0f,S)),vr));
}

///
/// van der Waals force
///
extern "C" __device__ void vanderWaals(float & vdw, float & dvdw, float rij, float r0v)
{
	float alfavm = __fmul_rn(ALFAV,(rij-r0v));
	float gmexp = __fmul_rn(GM,exp(-alfavm));
	float a = exp(alfavm) + gmexp;
	float b = __fmul_rn(ALFAV,(exp(alfavm) - gmexp));
	float db = __fmul_rn(ALFAVSQ,a);

	vdw = 1.0f+__fmul_rn(50.0f,GM)+GMSQ+__fmul_rn(__fmul_rn(10.0f,GM-1.0f),b);
	dvdw = __fmul_rn(-1.0e-4f,(__fmul_rn(10.0f,__fmul_rn(GM-1.0f,__fmul_rn(db,a)))-__fmul_rn(2.0f,__fmul_rn(b,vdw)))/__fmul_rn(a,__fmul_rn(a,a)));
	vdw = 1.0e-4f*vdw/__fmul_rn(a,a);
}

///
/// Cutoff
///
extern "C" __device__ void cutoff(float rij, float & fc, float & dfc)
{
	fc = tanh(__fmul_rn(ALFA,rij-R0));
	dfc = __fmul_rn(-0.5f,__fmul_rn(ALFA,1.0f - __fmul_rn(fc,fc)));
	fc = __fmul_rn(0.5f,1.0f - fc);
}