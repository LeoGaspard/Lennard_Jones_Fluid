///////////////////////////////////////////////////////
//NAME:			CAtom.cpp
//
//PURPOSE:		Definition of the CAtom
//			class
//
//FUNCTIONS/OBJECTS:	CAtom
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#include "CAtom.h"

//Constructor
CAtom::CAtom(double inSigma, double inE)
{
	m_dSigma = inSigma;
	m_dEpsilon = inE;

	CPos p;
	CSpeed s;
	CForce f;

	m_Position = p;
	m_Speed = s;
	m_Forces = f;

}//CAtom

//Constructor
CAtom::CAtom(double inSigma, double inE, CPos inP)
{
	m_dSigma = inSigma;
	m_dEpsilon = inE;
	m_Position = inP;

	CSpeed s;
	CForce f;

	m_Speed = s;
	m_Forces = f;
}//CAtom

//Constructor
CAtom::CAtom(double inSigma, double inE, CPos inP, CSpeed inS)
{
	m_dSigma = inSigma;
	m_dEpsilon = inE;
	m_Position = inP;
	m_Speed = inS;

	CForce f;

	m_Forces = f;
}//CAtom

//Constructor
CAtom::CAtom(double inSigma, double inE, CPos inP, CSpeed inS, CForce inF)
{
	m_dSigma = inSigma;
	m_dEpsilon = inE;
	m_Position = inP;
	m_Speed = inS;
	m_Forces = inF;
}//CAtom

//Destructor
CAtom::~CAtom()
{

}//~CAtom

//// Prints the position vector in the stream f
//void	CAtom::OutPos(std::ofstream& f)
//{
//	f << m_Position << std::endl;
//} // OutPos
//
//// Prints the speed vector in the stream f
//void	CAtom::OutSpeed(std::ofstream& f)
//{
//	f << m_Speed << std::endl;
//} // OutSpeed
//
//// Prints the force vector in the stream f
//void	CAtom::OutForce(std::ofstream& f)
//{
//	f << m_Forces << std::endl;
//} // OutSpeed


