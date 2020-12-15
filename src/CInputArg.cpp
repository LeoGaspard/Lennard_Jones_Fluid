///////////////////////////////////////////////////////
//NAME:			CInputArg.cpp
//
//PURPOSE:		Definition of the CInputArg
//			class
//
//FUNCTIONS/OBJECTS:	CInputArg
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#include "CInputArg.h"

// Constructor of the CInputArg class
CInputArg::CInputArg(std::string s,unsigned int n)
{
	m_sName = s;
	m_iNValue = n;

	if(n > 1)
	{
		m_vKeywords.resize(m_iNValue);
		m_vValues.resize(m_iNValue);
	}
	else
	{
		m_vValues.resize(1);
	}
} //CInputArg

// Destructor
CInputArg::~CInputArg()
{
	m_vKeywords.clear();
	m_vValues.clear();
} // ~CInputArg

