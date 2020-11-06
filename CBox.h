///////////////////////////////////////////////////////
//NAME:			CBox.h
//
//PURPOSE:		Definition of the CBox
//			class
//
//FUNCTIONS/OBJECTS:	CBox
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#ifndef CBOX_H_INCLUDED
#define CBOX_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/format.hpp>

#include <vector>
#include <math.h>
#include <random>

#include "Constants.h"
#include "CAtom.h"
#include "C3Vec.h"
#include "CPos.h"
#include "CSpeed.h"
#include "CForce.h"


class CBox
{
	private:
		double 			m_dA, m_dB, m_dC, m_dAlpha, m_dBeta, m_dGamma, m_dVolume;
		std::vector<CAtom>	m_vAtomList;
		C3Vec			m_AVec, m_BVec, m_CVec, m_UVec, m_VVec, m_WVec;

	public:
		// Constructors & destructor
					CBox(double inA=0.0, double inB=0.0, double inC=0.0, double inAlpha=90.0, double inBeta=90.0, double inGamma=90.0);
					~CBox();

		// Methods
		void			Wrap();
		void			InitPosFromRandomDistribution(unsigned int inNPoints, double inDMin);
		CAtom&			GetAtom(unsigned int n);

		// Output methods
		void			OutBoxParam(std::ofstream& f);
		void			OutAtomPos(std::ofstream& f);
		void			OutAtomSpeed(std::ofstream& f);
		void			OutAtomForces(std::ofstream& f);
};

#endif // CBOX_H_INCLUDED
