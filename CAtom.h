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

class CAtom
{
	private:
		double		m_dSigma, m_dEpsilon;
		CPos		m_Position;
		CSpeed		m_Speed;
		CForce		m_Forces;

	public:
				CAtom(double inSigma=0.0, double inZ=0.0);
				CAtom(double inSigma, double inZ,CPos inP);
				CAtom(double inSigma, double inZ,CPos inP, CSpeed inS);
				CAtom(double inSigma, double inZ,CPos inP, CSpeed inS, CForce inF);
				~CAtom();

		// Getters
		CPos		GetPos(){return m_Position;};
		CSpeed		GetSpeed(){return m_Speed;};
		CForce		GetForces(){return m_Forces;};
		double		GetSigma(){return m_dSigma;};
		double		GetEpsilon(){return m_dEpsilon;};

		// Setters
		void		SetPos(CPos p){m_Position = p;};
		void		SetSpeed(CSpeed s){m_Speed = s;};
		void		SetForces(CForce f){m_Forces = f;};
		void		Move(CPos p){m_Position += p;};
		void		ChangeSpeed(CSpeed s){m_Speed += s;};


		// Output methods
//		void		OutPos(std::ofstream& f);
//		void		OutSpeed(std::ofstream& f);
//		void		OutForce(std::ofstream& f);
	
};

#endif // CATOM_H_INCLUDED
