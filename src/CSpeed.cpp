///////////////////////////////////////////////////////
//NAME:			CSpeed.cpp
//
//PURPOSE:		Definition of the CSpeed
//			class
//
//FUNCTIONS/OBJECTS:	CSpeed
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#include "CSpeed.h"

// Overloading the assignment operator
CSpeed&	CSpeed::operator=(C3Vec v)
{
	m_dX = v.GetX();
	m_dY = v.GetY();
	m_dZ = v.GetZ();

	return *this;
}//operator=

// Overloading the assignment operator
CSpeed&	CSpeed::operator=(CSpeed v)
{
	m_dX = v.GetX();
	m_dY = v.GetY();
	m_dZ = v.GetZ();

	return *this;
}//operator=

