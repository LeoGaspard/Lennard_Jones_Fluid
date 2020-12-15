///////////////////////////////////////////////////////
//NAME:			CApplication.h
//
//PURPOSE:		Definition of the CApplication
//			class
//
//FUNCTIONS/OBJECTS:	CAppplication
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#ifndef CAPPLICATION_H_INCLUDED
#define CAPPLICATION_H_INCLUDED

#include <cstring>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <regex>
#include <vector>

#include "Constants.h"

class CApplication
{
	public:
		         		CApplication();
		virtual 		~CApplication();

	protected:
		void 			ParseCommandLineOptions(int argc,const char * argv[]);
		void			ParseInputFile();

	protected :
		std::ofstream		m_streamOutput;
		unsigned int		m_iNumberThreads;
		std::string		m_strInputFile, m_strOutputFile;

};

#endif  // CAPPLICATION_H_INCLUDED
