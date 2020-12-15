///////////////////////////////////////////////////////////////
//NAME				CApplication.cpp
//
//PURPOSE:			Definition of the CApplication
//				class
//
//FUNCTIONS/OBJECTS:		CApplication
//
//AUTHOR:			Léo Gaspard
///////////////////////////////////////////////////////////////


#include "CApplication.h"

//Constructor of the CApplication object
CApplication::CApplication()
{
	m_iNumberThreads = 1;
	m_strInputFile = "";
	m_strOutputFile = "";
} //CApplication

//Destructor
CApplication::~CApplication()
{

} // ~CApplication

//Parsing the command line options 
void CApplication::ParseCommandLineOptions(int argc,const char * argv[])
{
	// The help message to be printed if -h|--help is provided or if
	// no argument is given
	if((argc==2 and (strcmp(argv[1],"-h")==0 or strcmp(argv[1],"--help")==0)) or argc == 1)
	{
		std::cout << argv[0] << ":" << std::endl;
		std::cout << "	Usage : Projet -i <Input file name> -o <Output file name> [options]" << std::endl << std::endl;	
		std::cout << "	Options:" << std::endl;
		std::cout << " -h|--help		The help text" << std::endl;
		std::cout << " -i|--input		The input file name" << std::endl;
		std::cout << " -o|--output		The output file name" << std::endl;
		std::cout << " -n|--nthread <Number>	The number of threads" << std::endl;
		
		// The program will not run if the help needs to be printed
		exit(0);
	}
	// Reading all the command line options 
	for(int i=1; i<argc; i++)
	{
		if(strcmp(argv[i],"-i")==0 or strcmp(argv[i],"--input")==0)
		{
			m_strInputFile = argv[i+1];
			i++;
		}
		else if(strcmp(argv[i],"-o")==0 or strcmp(argv[i],"--output")==0)
		{
			m_strOutputFile = argv[i+1];
			i++;
		}
		else if(strcmp(argv[i],"-n")==0 or strcmp(argv[i],"--nthread")==0)
		{
			try
			{
				m_iNumberThreads = std::stoi(argv[i+1]);
				i++;
			}
			catch(...)
			{
				std::stringstream errMsg;
				errMsg << "Invalid number of threads  \"" <<  argv[i+1] << "\"";
				throw std::runtime_error(errMsg.str());
			}
		}
	}

} //m_ParseCommandLineOptions
