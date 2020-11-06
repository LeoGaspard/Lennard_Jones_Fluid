///////////////////////////////////////////////////////
//NAME:			CPos.cpp
//
//PURPOSE:		Definition of the CPos
//			class
//
//FUNCTIONS/OBJECTS:	CPos
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#include "CPos.h"

// Computes and returns the distance between this and p
double	CPos::Distance(CPos p)
{
	double x,y,z;

	x = m_dX - p.GetX();
	y = m_dY - p.GetY();
	z = m_dZ - p.GetZ();

	return sqrt(x*x+y*y+z*z);
}//Distance

// Computes and returns the square of the distance between this and p
double	CPos::Distance2(CPos p)
{
	double x,y,z;

	x = m_dX - p.GetX();
	y = m_dY - p.GetY();
	z = m_dZ - p.GetZ();

	return x*x+y*y+z*z;
}//Distance2

// Overloading the assignment operator
CPos&	CPos::operator=(C3Vec v)
{
	m_dX = v.GetX();
	m_dY = v.GetY();
	m_dZ = v.GetZ();

	return *this;
}//operator=

// Overloading the assignment operator
CPos&	CPos::operator=(CPos v)
{
	m_dX = v.GetX();
	m_dY = v.GetY();
	m_dZ = v.GetZ();

	return *this;
}//operator=
