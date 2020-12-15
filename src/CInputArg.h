///////////////////////////////////////////////////////
//NAME:			CInputArg.h
//
//PURPOSE:		Definition of the CInputArg
//			class
//
//FUNCTIONS/OBJECTS:	CInputArg
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#ifndef CINPUTARG_H_INCLUDED
#define CINPUTARG_H_INCLUDED

#include <iostream>
#include <vector>
#include <cstring>

class CInputArg
{
	private:
		std::string 			m_sName;
		unsigned int			m_iNValue;
		std::vector<std::string>	m_vKeywords, m_vValues;	

	public:

					CInputArg(std::string s,unsigned int n);
		virtual 		~CInputArg();

		// Check if a string is the same as the name of the object
		int			Compare(std::string s) {return strcmp(m_sName.data(),s.data());}; //m_Compare


		//Getters
		int			GetNValue() {return m_iNValue;};
		std::string		GetName() {return m_sName;};

		//Setters
		void			SetKeyword(unsigned int n,std::string s){m_vKeywords[n]=s;};
		void			SetValue(unsigned int n,std::string s){m_vValues[n]=s;};

};

#endif // CINPUTARG_H_INCLUDED
