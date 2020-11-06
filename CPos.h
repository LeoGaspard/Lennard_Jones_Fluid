///////////////////////////////////////////////////////
//NAME:			CPos.h
//
//PURPOSE:		Definition of the CPos
//			class
//
//FUNCTIONS/OBJECTS:	CPos
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#ifndef CPOS_H_INCLUDED
#define CPOS_H_INCLUDED

#include "C3Vec.h"

class CPos : public C3Vec
{
	public:
				CPos(double inX,double inY,double inZ) : C3Vec(inX,inY,inZ){};
				CPos() : C3Vec(){};
				~CPos(){}; 

		double		Distance(CPos p);
		double		Distance2(CPos p);

		// Operator
		
		CPos&		operator=(C3Vec v);
		CPos&		operator=(CPos v);

	private:	
};

#endif //CPOS_H_INCLUDED
