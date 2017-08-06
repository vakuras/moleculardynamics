///
/// Getter Funtions Definitions
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef GF_HEADER
#define GF_HEADER

#include "global.h"
#include "constants.h"

template <typename T>
T getMass(int id);
template <typename T>
T getRE(int ti, int tj);
template <typename T>
T getDE(int ti, int tj);
template <typename T>
T getbeta(int ti, int tj);
template <typename T>
T getR0V(int ti, int tj);

#include "getters.cpp" //c++ templates can't have companion cpp files... so we just include it in the header

#endif