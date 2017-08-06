///
/// Attractive, Repulsive, van der Waals & Cutoff Definitions
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef ARVDWFC_HEADER
#define ARVDWFC_HEADER

#include "global.h"
#include "constants.h"

template <typename T>
inline void attractive(T & va, T & dva, T rij, T De, T beta, T re);
template <typename T>
void repulsive(T & vr, T & dvr, T rij, T De, T beta, T re);
template <typename T>
void vanderWaals(T & vdw, T & dvdw, T rij, T r0v);
template <typename T>
void cutoff(T  rij, T & fc, T & dfc);

#include "arvdwfc.cpp" //c++ templates can't have companion cpp files... so we just include it in the header

#endif
