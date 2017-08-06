///
/// Configuration Parser implementation.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///
#include "stdafx.h"
#include "ConfigParser.h"

using namespace System::Windows::Forms;

namespace ui
{
	///
	/// Constructor
	///
	ConfigParser::ConfigParser()
	{
		data = new unordered_map<string, string>(); //create data holder
	}

	///
	/// Get a boolean value from the data holder
	///
	bool ConfigParser::Parse(string filename)
	{
		size_t pos;
		int lineid=0;
		string key;
		string value;

		data->clear();	//clear data holder

		ifstream ifs (filename.c_str(), ifstream::in); //open file
		string line;

		CPPSAFE_CALL(ifs.fail(), "Error opening configuration file.");

		for(;;)
		{
			line = "";

			getline(ifs, line);	//get a line

			if (ifs.eof() && line.empty()) //check for eof
				break;

			CPPSAFE_CALL(ifs.fail(), "Unable to read line from configuration file.");
			
			if (!line.substr(0,2).compare("//") || line.empty()) //skip comments and empty lines
				continue;	

			while((pos = line.rfind("//"))!=string::npos) //remove comments from the end
				line = line.substr(0,pos);

			pos = line.find("="); //locate equality sign

			if (pos==string::npos) //no equality sign - error!
			{
				MSGERR("Error in configuration file!\nLine: " + lineid);
				return false;
			}

			string key = line.substr(0,pos); //get key
			string value = line.substr(pos+1); //get value

			if (value.empty() || key.empty()) //check if key or value is/are empty
			{
				MSGERR("Error in configuration file!\nLine: " + lineid);
				return false;
			}

			//remove spaces and tabs from the end
			pos = value.size() - 1; 
			while (value[pos]==' ' || value[pos]=='\t')
			{
				pos--;

				if (pos<0)
					break;
			}
			if (pos>0 && pos!=value.size()-1)
				value = value.substr(0,pos+1);

			//uppercase
			for (unsigned int i = 0; i < key.size(); i++){key[i]=(char)toupper(key[i]);}

			//insert value to data holder
			data->insert(datamap::value_type(key,value));

			//count lines
			lineid++;
		}

		ifs.close();

		return true;
	}

	///
	/// Gets a string value from the data holder
	///
	bool ConfigParser::GetString(string key, string & value)
	{
		//uppercase
		for (unsigned int i = 0; i < key.size(); i++){key[i]=(char)toupper(key[i]);}

		if (data->find(key)==data->end())
			return false;

		value = (*data)[key];

		return true;
	}

	///
	/// Gets an int value from the data holder
	///
	bool ConfigParser::GetInt(string key, int & value)
	{
		//uppercase
		for (unsigned int i = 0; i < key.size(); i++){key[i]=(char)toupper(key[i]);}

		if (data->find(key)==data->end())
			return false;

		return from_string<int>(value, (*data)[key], dec);
	}

	///
	/// Gets a float value from the data holder
	///
	bool ConfigParser::GetFloat(string key, float & value)
	{
		//uppercase
		for (unsigned int i = 0; i < key.size(); i++){key[i]=(char)toupper(key[i]);}
		
		if (data->find(key)==data->end())
			return false;

		return from_string<float>(value, (*data)[key], dec);
	}

	///
	/// Gets a boolean value from the data holder
	///
	bool ConfigParser::GetBoolean(string key, bool & value)
	{
		//uppercase
		for (unsigned int i = 0; i < key.size(); i++){key[i]=(char)toupper(key[i]);}
		
		if (data->find(key)==data->end())
			return false;

		return from_string<bool>(value, (*data)[key], dec);
	}

	///
	/// Gets a set of floats from the data holder
	///
	bool ConfigParser::GetFloats(string key, float ** value)
	{
		float val;
		int counter = 0;

		//uppercase
		for (unsigned int i = 0; i < key.size(); i++){key[i]=(char)toupper(key[i]);}
		
		if (data->find(key)==data->end())
			return false;

		if (*value)
		{
			free (*value);
			*value = NULL;
		}

		string line = (*data)[key];

		while(from_string<float>(val, line, dec))
		{
			*value = (float *) realloc(*value, ++counter*sizeof(float));
			(*value)[counter-1] = val;

			size_t pos = line.find(" ");

			if (pos==string::npos)
				break;

			line = line.substr(pos+1);
		}

		return (counter>0);
	}

	///
	/// ConfigParser destructor
	///
	ConfigParser::~ConfigParser()
	{
		delete data;
	}

	///
	/// Extracting from c++ string
	///
	template <class T>
	bool ConfigParser::from_string(T& t, const string& s, ios_base& (*f)(ios_base&))
	{
		istringstream iss(s);
		return !(iss >> f >> t).fail();
	}
}