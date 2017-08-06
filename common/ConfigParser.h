///
/// Configuration Parser Definitions
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef CFGPARSER_H
#define CFGPARSER_H

#include <string>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;
using namespace tr1;

#ifdef UI
#define MSGERR(msg) MessageBox::Show(msg, "Error", MessageBoxButtons::OK, MessageBoxIcon::Error)
#define CPPSAFE_CALL(code,err) if(code) { MSGERR(err);}
#else
#define CPPSAFE_CALL(code,err) if(code) {cerr << "Error: " << err << "\nFile: '" << __FILE__ << "'\nLine: " << __LINE__ << endl; exit(EXIT_FAILURE);}
#endif

typedef unordered_map<string, string> datamap; 

class ConfigParser
{
private:
	datamap * data; //data holder

public:
	ConfigParser(); //constructor
	bool Parse(string filename); //parse configuration file
	bool GetString(string key, string & value); //gets a string from the data holder
	bool GetInt(string key, int & value); //gets an int from the data holder
	bool GetFloat(string key, float & value); //gets a float from the data holder
	bool GetFloats(string key, float ** value); //gets a set of floats from the data holder
	bool GetBoolean(string key, bool & value); //gets a boolean from the data holder
	~ConfigParser();
private:
	//extracts values from c++ string
	template <class T> bool from_string(T& t, const string& s, ios_base& (*f)(ios_base&));
};


#endif