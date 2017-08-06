///
/// Results Implementation
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "Results.h"
#include "Constants.h"

void pushAnim(vector<float*> & posVec, vector<float*> & velVec, real4* posArray, real3* velocityArray, int nop)
{
	float * pAcp;
	float * vAcp;

	TC(pAcp = new float[nop*4]);
	TC(vAcp = new float[nop*3]);

	for(int i=0; i<nop; i++)
	{
		pAcp[i*4] = (float)posArray[i].x;
		pAcp[i*4+1] = (float)posArray[i].y;
		pAcp[i*4+2] = (float)posArray[i].z;
		pAcp[i*4+3] = (float)posArray[i].w;
		vAcp[i*3] = (float)velocityArray[i].x;
		vAcp[i*3+1] = (float)velocityArray[i].y;
		vAcp[i*3+2] = (float)velocityArray[i].z;
	}

	posVec.push_back(pAcp);
	velVec.push_back(vAcp);
}

void flushVectors(vector<results> & resultVec, vector<float*> & posVec, vector<float*> & velVec)
{
	resultVec.clear();

	for(vector<float*>::iterator it=posVec.begin() ; it < posVec.end(); it++)
		delete [] *it;

	for(vector<float*>::iterator it=velVec.begin() ; it < velVec.end(); it++)
		delete [] *it;

	posVec.clear();
	velVec.clear();
}

void writeAnimationBinaryData(vector<float*> & posVec, vector<float*> & velVec, configuration * config, float vxcm)
{
	ofstream animdata;

	stringstream adname;
	adname << config->AnimFile << "(vxcm-" << vxcm << ").anm";

	//open the corresponding files for output
	animdata.open(adname.str().c_str(), fstream::binary | fstream::out | fstream::trunc);
	CPPSAFE_CALL(animdata.fail(), "Error opening result file.");

	int nop = config->LennardJonesParticles + config->ManyBodyParticles;
	int count = (int)posVec.size();

	animdata.write((char*) &nop, sizeof(int));
	animdata.write((char*) &(config->animts), sizeof(int));
	animdata.write((char*) &count, sizeof(int));

	//check for errors after write
	CPPSAFE_CALL(animdata.fail(), "Error writing to result file.");

	for(vector<float*>::iterator it=posVec.begin() ; it<posVec.end(); it++)
	{
		float * posArray = *it;

		for(int i=0; i<nop*4; i++)
			animdata.write((char*) (posArray+i), sizeof(float));
	}
			

	for(vector<float*>::iterator it=velVec.begin() ; it<velVec.end(); it++)
	{
		float * velArray = *it;

		for(int i=0; i<nop*3; i++)
			animdata.write((char*) (velArray+i), sizeof(float));
	}
}

string tailResults(real4 * posArray, int ManyBodyParticles)
{
	int pairs = 0;
	stringstream ssPairs;
	stringstream res;

	res << "ManyBody Pairs:\n";

	for(int i=0;i<ManyBodyParticles; i++)
		for(int j=i+1; j<ManyBodyParticles; j++)
		{
			if ((int)posArray[i].w == (int)posArray[j].w)
				continue;

			double dx = posArray[i].x - posArray[j].x;
			double dy = posArray[i].y - posArray[j].y;
			double dz = posArray[i].z - posArray[j].z;

			double dist = sqrt(dx*dx+dy*dy+dz*dz);

			if (dist<2*reNO)
			{
				pairs++;
				ssPairs << "[" << i << ", " << j << "] ";
			}
		}

	res << "Amount: " << pairs << "\n" << ssPairs.str();

	return res.str();
}

void writeResults(vector<results> & resultVec, configuration * config, float elapsed, string tail, float vxcm)
{
	ofstream result;
	ofstream logfile;

	stringstream rfname;
	stringstream lfname;
	rfname << config->Filename << "(vxcm-" << vxcm << ").dat";
	lfname << config->Filename << "(vxcm-" << vxcm << ").txt";

	//open the corresponding files for output
	result.open(rfname.str().c_str(), fstream::binary | fstream::out | fstream::trunc);
	CPPSAFE_CALL(result.fail(), "Error opening result file.");
	logfile.open(lfname.str().c_str(), ofstream::out | ofstream::trunc);
	CPPSAFE_CALL(logfile.fail(), "Error opening log file.");

	logfile << fixed << setprecision(6);
	logfile.setf(ios::fixed,ios::floatfield);

	logfile << "Number Of Particles = " << config->LennardJonesParticles + config->ManyBodyParticles << ", Timesteps = " << config->Timesteps << "\n";
	logfile << "DT = " << config->DT;

	if (config->useLennardJones) 
		logfile << ", LennardJonesRS = " << config->LennardJonesRS << ", LennardJonesRCUT = " << config->LennardJonesRCUT;

	if (config->useManyBody)
		logfile << ", ManyBodyRS = " << config->ManyBodyRS << ", ManyBodyRCUT = " << config->ManyBodyRCUT;

	logfile << "\n\n";

	logfile << "\nVXCM [" << vxcm << "]: Time took: " << elapsed << " ms\n\n";
	logfile << "t\t\tT\t\tE\t\tEK\t\tEU\t\tR.x\t\tR.y\t\tR.z\t\tp.x\t\tp.y\t\tp.z\n";

	//check for errors after write
	CPPSAFE_CALL(logfile.fail(), "Error writing to log file.");

	unsigned char usebyte;
	int nop =config->LennardJonesParticles + config->ManyBodyParticles;

	result.write((char*) &nop, sizeof(int));
	result.write((char*) &config->Timesteps, sizeof(int));
	result.write((char*) &config->DT, sizeof(float));

	usebyte = (config->useLennardJones) ? 0xFF : 0x00;

	result.write((char*)&usebyte, sizeof(unsigned char));

	if (config->useLennardJones)
	{
		result.write((char*) &config->LennardJonesRS, sizeof(float));
		result.write((char*) &config->LennardJonesRCUT, sizeof(float));
	}

	usebyte = (config->useManyBody) ? 0xFF : 0x00;

	result.write((char*)&usebyte, sizeof(unsigned char));

	if (config->useManyBody)
	{
		result.write((char*) &config->ManyBodyRS, sizeof(float));
		result.write((char*) &config->ManyBodyRCUT, sizeof(float));
	}

	result.write((char*)&vxcm,sizeof(float));

	//check for errors after write
	CPPSAFE_CALL(result.fail(), "Error writing to result file.");

	for(vector<results>::iterator it=resultVec.begin() ; it < resultVec.end(); it++)
	{
		logfile << endl;

		result.write((char *)&(it->time), sizeof(double));
		result.write((char *)&(it->temperature), sizeof(double));
		result.write((char *)&(it->e), sizeof(double));
		result.write((char *)&(it->ek), sizeof(double));
		result.write((char *)&(it->eu), sizeof(double));
		result.write((char *)&(it->centerOfMassx), sizeof(double));
		result.write((char *)&(it->centerOfMassy), sizeof(double));
		result.write((char *)&(it->centerOfMassz), sizeof(double));
		result.write((char *)&(it->momentumx), sizeof(double));		
		result.write((char *)&(it->momentumy), sizeof(double));		
		result.write((char *)&(it->momentumz), sizeof(double));		

		//check for errors after read
		CPPSAFE_CALL(result.fail(), "Error writing to result file.");
		
		logfile << it->time;
		logfile << scientific;
		logfile << "\t" << it->temperature;
		logfile << "\t" << it->e;
		logfile << "\t" << it->ek;
		logfile << "\t" << it->eu;
		logfile << "\t" << it->centerOfMassx;
		logfile << "\t" << it->centerOfMassy;
		logfile << "\t" << it->centerOfMassz;
		logfile << "\t" << it->momentumx;
		logfile << "\t" << it->momentumy;
		logfile << "\t" << it->momentumz;
		logfile << fixed;

		//check for errors after write
		CPPSAFE_CALL(logfile.fail(), "Error writing to log file.");
	}

	logfile << "\n\n" << tail;

	result.close();
	logfile.close();
}