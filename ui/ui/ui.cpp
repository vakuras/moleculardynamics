///
/// Molecular::Dynamics main class.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "stdafx.h"
#include "LogoWindow.h"
#include "MainWindow.h"
#include "ConfigParser.h"

namespace ui
{
	using namespace msclr::interop;

	///
	/// Read input file
	///
	bool ReadInput(const char * filename, int nop, float * posArray, float * velocityArray)
	{
		ifstream input;

		//open input for reading
		input.open(filename, ifstream::in);
		CPPSAFE_CALL(!input.good(), "Error opening input file.");

		//copy data from input file
		for (int i=0;i<nop;i++)
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

		return true;
	}

	///
	/// Convert .NET String to C String
	///
	string MakeCSTR(System::String ^ from)
	{
		marshal_context ^ context = gcnew marshal_context();
		const char * str = context->marshal_as<const char*>(from);
		string sstr = string(str);
		delete context;
		return sstr;
	}

	///
	/// Read configuration from file
	///
	bool ReadConfiguration(string filename, configuration * pconfig)
	{
		bool value = true;

		ConfigParser cp;
		if (!cp.Parse(filename))
			return false;

		memset(pconfig, 0, sizeof(configuration));

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

		return value;
	}

	///
	/// Write configuration to ini file
	///
	bool WriteConfiguration(string filename, configuration * pconfig)
	{
		ofstream ofs (filename.c_str(), ifstream::out);
		string line;

		CPPSAFE_CALL(ofs.fail(), "Error writing configuration file.");

		ofs << NOP << "=" << pconfig->LennardJonesParticles << endl;
		ofs << MBP << "=" << pconfig->ManyBodyParticles << endl;
		ofs << TIMESTEPS << "=" << pconfig->Timesteps << endl;
		ofs << scientific << DTQ << "=" << pconfig->DT << endl;
		ofs << TEMP << "=" << pconfig->Temperature << endl;
		ofs << LJRSQ << "=" << pconfig->LennardJonesRS << endl;
		ofs << LJRCUTQ << "=" << pconfig->LennardJonesRCUT << endl;
		ofs << MBRSQ << "=" << pconfig->ManyBodyRS << endl;
		ofs << MBRCUTQ << "=" << pconfig->ManyBodyRCUT << endl;
		ofs << USELJ << "=" << pconfig->useLennardJones << endl;
		ofs << USEMB << "=" << pconfig->useManyBody << endl;
		ofs << QDEBUG << "=" << pconfig->Debug << endl;
		ofs << OUTPUTTS << "=" << pconfig->OutputTimesteps << endl;
		ofs << INPUT << "=" << pconfig->Input << endl;
		ofs << OUTPUT << "=" << pconfig->Output << endl;
		ofs << FILENAME << "=" << pconfig->Filename << endl;
		ofs << USECUDA << "=" << pconfig->UseCuda << endl;
		ofs << CUDABLOCKS << "=" << pconfig->CudaBlocks << endl;
		ofs << FALLBACK << "=" << pconfig->Fallback << endl;
		ofs << LJBPP << "=" << pconfig->lennardJonesBPP << endl;
		ofs << MBBPP << "=" << pconfig->manyBodyBPP << endl;
		ofs << VXCMFROM << "=" << pconfig->vxcmFrom << endl;
		ofs << VXCMTO << "=" << pconfig->vxcmTo << endl;
		ofs << VXCMSTEP << "=" << pconfig->vxcmStep << endl;
		ofs << ENERGYLOSS << "=" << pconfig->energyLoss << endl;
		ofs << ANIMTS << "=" << pconfig->animts << endl;
		ofs << ANIMDATA << "=" << pconfig->AnimFile << endl;

		CPPSAFE_CALL(ofs.fail(), "Error writing configuration file.");

		ofs.close();

		return true;
	}
}

using namespace ui;

///
/// Main .NET function
///
[STAThreadAttribute]
int main(array<System::String ^> ^)
{
	// Enabling Windows XP visual effects before any controls are created
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false);

	//create main window
	MainWindow ^ main = MainWindow::GetInstance();
	main->Show();

	//pop logo window
	LogoWindow::GetInstance(main);

	//run the application as the main window
	Application::Run(main);
	return 0;
}
