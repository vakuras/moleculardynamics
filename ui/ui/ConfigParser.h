///
/// Configuration Parser definitions.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef CFGPARSER_H
#define CFGPARSER_H

#include "global.h"
#include <unordered_map>

using namespace std;
using namespace tr1;

namespace ui
{
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
}

#endif