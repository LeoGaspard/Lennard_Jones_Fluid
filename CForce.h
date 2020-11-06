///////////////////////////////////////////////////////
//NAME:			CForce.h
//
//PURPOSE:		Definition of the CForce
//			class
//
//FUNCTIONS/OBJECTS:	CForce
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#ifndef CFORCE_H_INCLUDED
#define CFORCE_H_INCLUDED

#include "C3Vec.h"

class CForce : public C3Vec
{
	public:
				CForce(double inX,double inY,double inZ) : C3Vec(inX,inY,inZ){};
				CForce() : C3Vec(){};
				~CForce(){}; 

		// Operator
		
		CForce&		operator=(C3Vec v);
		CForce&		operator=(CForce v);

	private:	
};

#endif //CFORCE_H_INCLUDED
