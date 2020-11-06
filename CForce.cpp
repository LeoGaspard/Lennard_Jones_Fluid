///////////////////////////////////////////////////////
//NAME:			CForce.cpp
//
//PURPOSE:		Definition of the CForce
//			class
//
//FUNCTIONS/OBJECTS:	CForce
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#include "CForce.h"

// Overloading the assignment operator
CForce&	CForce::operator=(C3Vec v)
{
	m_dX = v.GetX();
	m_dY = v.GetY();
	m_dZ = v.GetZ();

	return *this;
}//operator=

CForce&	CForce::operator=(CForce v)
{
	m_dX = v.GetX();
	m_dY = v.GetY();
	m_dZ = v.GetZ();

	return *this;
}//operator=
