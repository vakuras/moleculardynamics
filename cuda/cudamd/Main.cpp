///
/// Main for CUDA Simulation
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#define APPEXE "cudamd.exe"

#include "..\..\common\global.h"
#include "..\..\common\loader.h"
#include "host.h"

///
/// Main function.
///
int main(int argc, char ** argv)
{
	configuration   config; //simulation configuration structure
	int				deviceCount; //cuda device counter
	CUdevice		device; //cuda device
	char* 			deviceName=(char*)malloc(200*sizeof(char)); //cuda device name
	int				deviceNameLen=200*sizeof(char); //cuda device name length
	unsigned int    deviceMem; //cuda device memory
	int				selectedDevice=-1; //cuda selected device
	char*			module_path;

	Loader(argc, argv, APPEXE, config); //load configuration

	//output logo...
	std::cout << "Molecular::Dynamics\n";

	//init cuda
	CU_SAFE_CALL(cuInit(0));

	//get device count
	CU_SAFE_CALL(cuDeviceGetCount(&deviceCount));

	//no devices
	if(deviceCount == 0)
    {
        std::cerr << "No devices available!\n";
		return EXIT_FAILURE;
    }

	//list available CUDA devices
	for(int i = 0; i < deviceCount; i++)
    {
		CU_SAFE_CALL(cuDeviceGet(&device, i));
        CU_SAFE_CALL(cuDeviceGetName(deviceName, deviceNameLen, device));
		std::cout << "Found device " << i << ": " << deviceName << std::endl;
    }
        
	//ask user to select cuda device to work with if there are more then one device
	if (deviceCount != 1)
	{
		std::cout << "Select a device from the list above: ";

		while (selectedDevice<0 || selectedDevice>deviceCount-1)
			std::cin >> selectedDevice;
	}
	else //only one device available - so select it
		selectedDevice = 0;

	CU_SAFE_CALL(cuDeviceGet(&device, selectedDevice));
	CU_SAFE_CALL(cuDeviceGetName(deviceName, deviceNameLen, device));
	CU_SAFE_CALL(cuDeviceTotalMem(&deviceMem, device));

	//module path - not working correctly on windows xp x64, not sure for xp x32
	//				works correctly on win7/vista x64,x32
	module_path = cutFindFilePath("Device.ptx", argv[0]);

	/*
	workaround for cutFindFilePath
	*/
	/*
	TCHAR dirBuf[MAX_PATH+1];
	if (!GetCurrentDirectory(MAX_PATH, dirBuf))
		std::cerr << "\nGetCurrentDirectory failed!\n";

	stringstream dirBuilder;
	dirBuilder << dirBuf << "\\data\\Device.ptx";
	module_path = (char*)(dirBuilder.str().c_str());*/
	
	//print selected device
	std::cout << "\nSelected device is:\n(" << selectedDevice << ") " <<  deviceName << " Memory: " << deviceMem/1048576 << "MB\n";

	return hostMain(device, module_path, &config); //begin
}