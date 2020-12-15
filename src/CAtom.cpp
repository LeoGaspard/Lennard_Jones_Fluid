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
CAtom::CAtom(double inMass,double inSigma, double inE)
{
	m_dSigma = inSigma;
	m_dEpsilon = inE;
	m_dMass = inMass;

	CPos p;
	CSpeed s;
	CForce f;

	m_Position = p;
	m_Speed = s;
	m_Forces = f;

}//CAtom

//Constructor
CAtom::CAtom(double inMass,double inSigma, double inE, CPos inP)
{
	m_dSigma = inSigma;
	m_dMass = inMass;
	m_dEpsilon = inE;
	m_Position = inP;

	CSpeed s;
	CForce f;

	m_Speed = s;
	m_Forces = f;
}//CAtom

//Constructor
CAtom::CAtom(double inMass,double inSigma, double inE, CPos inP, CSpeed inS)
{
	m_dSigma = inSigma;
	m_dMass = inMass;
	m_dEpsilon = inE;
	m_Position = inP;
	m_Speed = inS;

	CForce f;

	m_Forces = f;
}//CAtom

//Constructor
CAtom::CAtom(double inMass,double inSigma, double inE, CPos inP, CSpeed inS, CForce inF)
{
	m_dSigma = inSigma;
	m_dMass = inMass;
	m_dEpsilon = inE;
	m_Position = inP;
	m_Speed = inS;
	m_Forces = inF;
}//CAtom

//Destructor
CAtom::~CAtom()
{

}//~CAtom

//Computes the kinetic energy E = 1/2 mv^2
void	CAtom::ComputeKineticEnergy()
{
	// Dal.bohr.fs^-2
	m_dKineticEnergy = 0.5*m_dMass*m_Speed.Norm2();
}//ComputeKineticEnergy
