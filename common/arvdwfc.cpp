///
/// Attractive, Repulsive, van der Waals & Cutoff Implementation
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef ARVDWFC_CPP
#define ARVDWFC_CPP
#include "arvdwfc.h"

///
/// Attractive force
///
template <typename T>
void attractive(T & va, T & dva, T rij, T De, T beta, T re)
{
	va = S*De/(S-1.0f)*exp(-beta*sqrt(2.0f/S)*(rij - re));
	dva = beta*sqrt(2.0f/S)*va;
}

///
/// Repulsive force
///
template <typename T>
void repulsive(T & vr, T & dvr, T rij, T De, T beta, T re)
{
	vr = De/(S-1.0f)*exp(-beta*sqrt(S*2.0f)*(rij-re)) + 7.0e-8f*exp(-10.0f*beta*(rij-re));
	dvr = beta*sqrt(2.0f*S)*vr;
}

///
/// van der Waals force
///
template <typename T>
void vanderWaals(T & vdw, T & dvdw, T rij, T r0v)
{
	T alfavm = ALFAV*(rij-r0v);
	T gmexp = GM*exp(-alfavm);
	T a = exp(alfavm) + gmexp;
	T b = ALFAV*(exp(alfavm) - gmexp);
	T db = ALFAVSQ*a;

	vdw = 1.0f+50.0f*GM+GMSQ+10.0f*(GM-1.0f)*b;
	dvdw = -1.0e-4f*((10.0f*(GM-1.0f)*db*a-2.0f*b*vdw)/(a*a*a));
	vdw = 1.0e-4f*vdw/(a*a);
}

///
/// Cutoff
///
template <typename T>
void cutoff(T  rij, T & fc, T & dfc)
{
	fc = tanh(ALFA*(rij-R0));
	dfc = -0.5f*ALFA*(1.0f - fc*fc);
	fc = 0.5f*(1.0f - fc);
}

#endif