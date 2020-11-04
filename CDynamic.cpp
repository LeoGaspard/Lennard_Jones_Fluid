///////////////////////////////////////////////////////
//NAME:			CDynamic.cpp
//
//PURPOSE:		Definition of the CDynamic
//			class
//
//FUNCTIONS/OBJECTS:	CDynamic
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////


#include "CDynamic.h"

//Constructor
CDynamic::CDynamic()
{

} //CDynamic

//Destructor
CDynamic::~CDynamic()
{

}//~CDynamic

// This setup the dynamic simulation using the provided input
void CDynamic::Setup(int argc, const char * argv[])
{
	this->ParseCommandLineOptions(argc, argv);
	//DEBUG
	std::cout << "Input  file       : " << m_strInputFile << std::endl;
	std::cout << "Output file       : " << m_strOutputFile << std::endl;
	std::cout << "Number of threads : " << m_iNumberThreads << std::endl << std::endl;
	//END DEBUG
	this->ParseInputFile();

}//Setup

//Parsing and checking the validity of the input file
void CDynamic::ParseInputFile()
{
	// These are the possible input variables
	double 					dA(0.0), dB(0.0), dC(0.0), dAlpha(90.0), dBeta(90.0), dGamma(90.0);
	std::vector<unsigned int> 		iNatom;
	std::vector<double>			dSigma, dEpsilon;
	
	//This regex will remove blank spaces before and after the first and last char
	std::regex trim("^\\s+|\\s+$");

	std::ifstream input(m_strInputFile);
	if(input.is_open())
	{
		// The input file has been successfully opened

		for(std::string strLine; getline(input,strLine); )
		{
			// This trims the line
			strLine = regex_replace(strLine,trim,"");
			// This removes the comments (starting with "!")
			strLine = strLine.substr(0,strLine.find("!"));

			std::istringstream 	streamLine(strLine);
			std::string 		strFirstParam;

			std::getline(streamLine,strFirstParam,' ');

			if(strcmp(strFirstParam.data(),"Box")==0)
			{

				for(std::string strBoxParam; strcmp(strBoxParam.data(),"End")!=0; getline(input, strBoxParam))
				{
					strBoxParam = regex_replace(strBoxParam,trim,""); 
					strBoxParam = strBoxParam.substr(0,strBoxParam.find("!"));

					std::istringstream streamBoxParam(strBoxParam);
					std::string strFirstParam;

					std::getline(streamBoxParam,strFirstParam,' ');

					if(strcmp(strFirstParam.data(),"a")==0)
					{
						std::string strValue;
						std::getline(streamBoxParam, strValue, ' ');

						dA = std::stod(strValue);

						// We check if there was a unit specification
						std::getline(streamBoxParam, strValue, ' ');
						if(strcmp(strValue.data(),"A")==0)
						{
							dA *= ANGSTROM_TO_BOHR;
						}
					}
					else if(strcmp(strFirstParam.data(),"b")==0)
					{
						std::string strValue;
						std::getline(streamBoxParam, strValue, ' ');

						dB = std::stod(strValue);

						// We check if there was a unit specification
						std::getline(streamBoxParam, strValue, ' ');
						if(strcmp(strValue.data(),"A")==0)
						{
							dB *= ANGSTROM_TO_BOHR;
						}
					}
					else if(strcmp(strFirstParam.data(),"c")==0)
					{
						std::string strValue;
						std::getline(streamBoxParam, strValue, ' ');

						dC = std::stod(strValue);

						// We check if there was a unit specification
						std::getline(streamBoxParam, strValue, ' ');
						if(strcmp(strValue.data(),"A")==0)
						{
							dC *= ANGSTROM_TO_BOHR;
						}
					}
					else if(strcmp(strFirstParam.data(),"alpha")==0)
					{
						std::string strValue;
						std::getline(streamBoxParam, strValue, ' ');

						dAlpha = std::stod(strValue);
					}
					else if(strcmp(strFirstParam.data(),"beta")==0)
					{
						std::string strValue;
						std::getline(streamBoxParam, strValue, ' ');

						dBeta = std::stod(strValue);
					}
					else if(strcmp(strFirstParam.data(),"gamma")==0)
					{
						std::string strValue;
						std::getline(streamBoxParam, strValue, ' ');

						dGamma = std::stod(strValue);
					}
					else if(strcmp(strBoxParam.data(),"")!=0)
					{
						std::cerr << "Fatal Error : Bad input with line " << std::endl << strBoxParam << std::endl;
						exit(1);
					}
				}

				//DEBUG
				std::cout << "a is " << dA << " bohr" << std::endl << "b is " << dB << " bohr" << std::endl << "c is " << dC  << " bohr" << std::endl << "alpha is " << dAlpha << "°" << std::endl << "beta is " << dBeta << "°" << std::endl << "gamma is " << dGamma << "°" << std::endl << std::endl;
				//END OF DEBUG
			}
		}
		input.close();
	}
	else
	{
		// The input file can not be opened
		std::cerr << "Fatal Error : The input file " << m_strInputFile << " can not be opened" << std::endl << "The program will stop." << std::endl;
		exit(1);
	}
	
} //m_ParseInputFile

