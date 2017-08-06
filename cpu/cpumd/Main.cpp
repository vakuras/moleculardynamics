///
/// Main for CPU Simulation
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#define APPEXE "cpumd.exe"

#include "..\..\common\global.h"
#include "..\..\common\loader.h"
#include "host.h"

///
/// Main function.
///
int main(int argc, char ** argv)
{
	configuration config; //simulation configuration structure
	Loader(argc, argv, APPEXE, config); //load configuration

	//output logo...
	cout << "Molecular::Dynamics\n";

	return hostMain(&config); //begin
}