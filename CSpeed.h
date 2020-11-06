///////////////////////////////////////////////////////
//NAME:			CSpeed.h
//
//PURPOSE:		Definition of the CSpeed
//			class
//
//FUNCTIONS/OBJECTS:	CSpeed
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#ifndef CSPEED_H_INCLUDED
#define CSPEED_H_INCLUDED

#include "C3Vec.h"

class CSpeed : public C3Vec
{
	public:
				CSpeed(double inX,double inY,double inZ) : C3Vec(inX,inY,inZ){};
				CSpeed() : C3Vec(){};
				~CSpeed(){}; 

		// Operator
		
		CSpeed&		operator=(C3Vec v);
		CSpeed&		operator=(CSpeed v);

	private:	
};

#endif //CSPEED_H_INCLUDED
