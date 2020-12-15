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
#include <omp.h>

#include "Constants.h"
#include "CAtom.h"
#include "C3Mat.h"
#include "C3Vec.h"
#include "CPos.h"
#include "CSpeed.h"
#include "CForce.h"


class CBox
{
	private:
		double 					m_dA, m_dB, m_dC, m_dAlpha, m_dBeta, m_dGamma, m_dVolume;
		double					m_dPotEnergy,m_dKinEnergy,m_dDensity;
		std::vector<CAtom>			m_vAtomList;
		C3Mat					m_H;
		std::vector<unsigned int>		m_vNeighborList;
		std::vector<CPos>			m_vPosList;
		std::vector<double>			m_vRadialDistributionFunction;

	public:
		// Constructors & destructor
					CBox(double inA=0.0, double inB=0.0, double inC=0.0, double inAlpha=90.0, double inBeta=90.0, double inGamma=90.0);
					~CBox();

		// Methods
		void			Setup();
		void			AddAtoms(unsigned int inNumber, double inSigma, double inEpsilon,double inMass);
		void			Wrap();
		void			InitPosFromRandomDistribution(double inDMin);
		void			InitSpeedRandom(double inTemperature);
		double			ComputeTemperature();
		double			ComputeForces();
		void			NeighborList(double cutoff, double neighbor);
		bool			CheckNeighborList(double neighbor);
		void			UpdatePositions(double deltaT);
		double			UpdateSpeeds(double deltaT);
		void			ComputeRadialDistributionFunction(double inStep, double inDMax,int inSteps);

		// Getters
		CAtom&			GetAtom(unsigned int n);
		double			GetTemperature(){return this->ComputeTemperature();};
		double			GetPotEnergy(){return m_dPotEnergy;};
		double			GetKinEnergy(){return m_dKinEnergy;};
		double			GetVolume(){return m_dVolume;};
		double			GetDensity(){return m_vAtomList.size()/m_dVolume;};
		unsigned int		GetNAtom(){return m_vAtomList.size();};

		// Setters
		void			SetA(double inA){m_dA = inA;};
		void			SetB(double inB){m_dB = inB;};
		void			SetC(double inC){m_dC = inC;};
		void			SetAlpha(double inAlpha){m_dAlpha = inAlpha;};
		void			SetBeta(double inBeta){m_dBeta = inBeta;};
		void			SetGamma(double inGamma){m_dGamma = inGamma;};

		// Output methods
		void			OutBoxParam(std::ofstream& f);
		void			OutAtomPos(std::ofstream& f);
		void			OutAtomSpeed(std::ofstream& f);
		void			OutAtomForces(std::ofstream& f);
		void			OutRadialDistributionFunction(double inStep, std::ofstream& f);
};

#endif // CBOX_H_INCLUDED
