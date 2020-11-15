///////////////////////////////////////////////////////
//NAME:			CAtom.h
//
//PURPOSE:		Definition of the CAtom
//			class
//
//FUNCTIONS/OBJECTS:	CAtom
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#ifndef CATOM_H_INCLUDED
#define CATOM_H_INCLUDED

#include "CPos.h"
#include "CSpeed.h"
#include "CForce.h"
#include "Constants.h"

class CAtom
{
	private:
		double		m_dSigma, m_dEpsilon,m_dMass,m_dKineticEnergy;
		CPos		m_Position;
		CSpeed		m_Speed;
		CForce		m_Forces;

	public:
				CAtom(double inMass,double inSigma=0.0, double inZ=0.0);
				CAtom(double inMass,double inSigma, double inZ,CPos inP);
				CAtom(double inMass,double inSigma, double inZ,CPos inP, CSpeed inS);
				CAtom(double inMass,double inSigma, double inZ,CPos inP, CSpeed inS, CForce inF);
				~CAtom();
		void		ComputeKineticEnergy();

		// Getters
		CPos		GetPos(){return m_Position;};
		CSpeed		GetSpeed(){return m_Speed;};
		CForce		GetForces(){return m_Forces;};
		double		GetSigma(){return m_dSigma;};
		double		GetEpsilon(){return m_dEpsilon;};
		double		GetKineticEnergy(){return m_dKineticEnergy;};
		double		GetMass(){return m_dMass;};

		// Setters
		void		SetPos(CPos p){m_Position = p;};
		void		SetSpeed(CSpeed s){m_Speed = s;};
		void		SetForces(CForce f){m_Forces = f;};
		void		Move(CPos p){m_Position += p;};
		void		AddForce(CForce f){m_Forces += f;};
		void		ChangeSpeed(CSpeed s){m_Speed += s;};
		void		SetMass(double m){m_dMass = m;};
};

#endif // CATOM_H_INCLUDED
