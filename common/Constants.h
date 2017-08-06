///
/// Constants Header
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef CONSTANTS_H
#define CONSTANTS_H

//Particles index
#define N 0
#define O 1
#define NG 2

//UI particle radii
#define radN 0.71f
#define radO 0.63f
#define radNG 0.58f

//the units are Ang.
#define reN    1.097f
#define reO    1.2075f
#define reNO   1.15f

//VDW parameters

#define sigNvdw   3.31f
#define sigOvdw   2.95f
#define reNvdw   3.72f
#define reOvdw   3.31f
#define reNOvdw  3.51f

//VDW parameters for Ne - Ne and Ne - N/O
#define sigNG    2.72f
#define sigNGN   3.015f
#define sigNGO   2.835f

//VDW parameters for Ne - Ne and Ne - N/O
#define epsNG   0.391e-4f
#define epsNGN  0.328e-4f
#define epsNGO  0.447e-4f

#define deN   941.4e-4f
#define deO   493.55309e-4f
#define deNO  626.734e-4f

// The units are 1/Ang
#define betaN  2.71009f
#define betaO  2.2478f
#define betaNO 2.24068f

// Masses: The units are gr/mol
#define mN 14.00674f
#define mO 15.9994f
#define mNG 20.1797f

// Cutoffs
#define CUTOFF_RIJ 15*sigNvdw
#define CUTOFF_RIJSQ 225*sigNvdw*sigNvdw

// Device Constants
#define GEAR1			0.15833333333333333333333333333333f
#define GEAR2			0.75f
#define GEAR3			0.5f
#define GEAR4			0.08333333333333333333333333333333f

#define S					2.0f
#define GM				0.3f
#define GMSQ			0.09f
#define ALFAV			3.0f
#define ALFAVSQ	9.0f
#define ALFA				2.0f
#define R0					1.3f
#define LAMDA3		1.33f
#define C					25.f
#define D					4.0417f
#define H					0.0f

#endif