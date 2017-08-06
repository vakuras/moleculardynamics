///
/// Loader Implementation
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "Loader.h"

///
/// Main function.
///
void Loader(int argc, char ** argv, char * appexe, configuration & config)
{
	//default configuration file
	string configfile = "default.cfg";

	if (argc==2)
		configfile = string(argv[1]); //get the configuration file from command line
	else if (argc>2)
	{
		cerr << "usage: " << appexe << " [configuration-file].\n"; //wrong usage
		exit(EXIT_FAILURE);
	}

	//read the configuration file
	readConfiguration(configfile, &config);
}

///
/// Read configuration from file
///
void readConfiguration(string filename, configuration * pconfig)
{
	bool value = true;

	//create an instance of the configuration parser & parse configuration file
	ConfigParser cp;
	if (!cp.Parse(filename)) 
	{
		cerr << "Error in configuration file!" << endl;
		exit(EXIT_FAILURE);
	}

	memset(pconfig, 0, sizeof(configuration)); //clear configuration structure

	//read configuration data from the parser
	value &= cp.GetInt(NOP, pconfig->LennardJonesParticles);
	value &= cp.GetInt(TIMESTEPS, pconfig->Timesteps);
	value &= cp.GetFloat(DTQ, pconfig->DT);
	value &= cp.GetFloat(TEMP, pconfig->Temperature);
	value &= cp.GetFloat(LJRSQ, pconfig->LennardJonesRS);
	value &= cp.GetFloat(LJRCUTQ, pconfig->LennardJonesRCUT);
	value &= cp.GetFloat(MBRSQ, pconfig->ManyBodyRS);
	value &= cp.GetFloat(MBRCUTQ, pconfig->ManyBodyRCUT);
	value &= cp.GetString(INPUT, pconfig->Input);
	value &= cp.GetString(OUTPUT, pconfig->Output);
	value &= cp.GetString(FILENAME, pconfig->Filename);
	value &= cp.GetInt(CUDABLOCKS, pconfig->CudaBlocks);
	value &= cp.GetBoolean(USECUDA, pconfig->UseCuda);
	value &= cp.GetInt(OUTPUTTS, pconfig->OutputTimesteps);
	value &= cp.GetBoolean(QDEBUG, pconfig->Debug);
	value &= cp.GetBoolean(FALLBACK, pconfig->Fallback);
	value &= cp.GetBoolean(USELJ, pconfig->useLennardJones);
	value &= cp.GetBoolean(USEMB, pconfig->useManyBody);
	value &= cp.GetBoolean(LJBPP, pconfig->lennardJonesBPP);
	value &= cp.GetBoolean(MBBPP, pconfig->manyBodyBPP);
	value &= cp.GetFloat(VXCMFROM, pconfig->vxcmFrom);
	value &= cp.GetFloat(VXCMTO, pconfig->vxcmTo);
	value &= cp.GetFloat(VXCMSTEP, pconfig->vxcmStep);
	value &= cp.GetFloat(ENERGYLOSS, pconfig->energyLoss);
	value &= cp.GetInt(MBP, pconfig->ManyBodyParticles);
	value &= cp.GetInt(ANIMTS, pconfig->animts);
	value &= cp.GetString(ANIMDATA, pconfig->AnimFile);

	//one of the configuration attributes is missing or wrong
	if (!value)
	{
		cerr << "Unable to read configuration file!" << endl;
		exit(EXIT_FAILURE);
	}
}

///
/// Write output file
///
void writeOutput(configuration * config, real * posArray, real * velocityArray)
{
	ofstream output;

    //open input for reading
	output.open(config->Output.c_str(), ofstream::out | ofstream::trunc);
	CPPSAFE_CALL(!output.good(), "Error opening output file.");

	output << scientific; //scientific output 0.0e0

    //copy data to output file
    for (int i=0;i<config->LennardJonesParticles + config->ManyBodyParticles;i++)
    {
        //position & mass
        for (int j=0;j<4;j++)
        {
			output << *posArray << " ";
            posArray++;
        }

        //velocity
        for (int j=0;j<3;j++)
        {
			output << *velocityArray << " ";
            velocityArray++;
        }

		//check for errors after write
		CPPSAFE_CALL(!output.good(), "Error writing to output file file.");
    }

    output.close();
}

///
/// Read input file
///
void readInput(configuration * config, real * posArray, real * velocityArray)
{
	ifstream input;

    //open input for reading
	input.open(config->Input.c_str(), ifstream::in);
	CPPSAFE_CALL(!input.good(), "Error opening input file.");

    //copy data from input file
    for (int i=0;i<config->LennardJonesParticles + config->ManyBodyParticles;i++)
    {
        //position & mass
        for (int j=0;j<4;j++)
        {
			input >> *posArray;
            posArray++;
        }

        //velocity
        for (int j=0;j<3;j++)
        {
			input >> *velocityArray;
            velocityArray++;
        }

		//check for errors after read
		CPPSAFE_CALL(!input.good(), "Error reading from input file.");
    }

    input.close();
}