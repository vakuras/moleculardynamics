///
/// Getter Functions Implementation
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef GF_CPP
#define GF_CPP
#include "Getters.h"

template <typename T>
T getMass(int id)
{
	if (id == N)
		return mN;
	else if (id == O)
		return mO;
	else
		return mNG;
}

template <typename T>
T getRE(int ti, int tj)
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

template <typename T>
T getDE(int ti, int tj)
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

template <typename T>
T getbeta(int ti, int tj)
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

template <typename T>
T getR0V(int ti, int tj)
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

#endif