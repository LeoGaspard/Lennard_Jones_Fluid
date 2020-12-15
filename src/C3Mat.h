///////////////////////////////////////////////////////
//NAME:			C3Mat
//
//PURPOSE:		Definition of the C3Mat
//			class
//
//FUNCTIONS/OBJECTS:	C3Mat
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#ifndef C3MAT_H_INCLUDED
#define C3MAT_H_INCLUDED

#include <iostream>
#include <boost/format.hpp>
#include <math.h>
#include "C3Vec.h"


class C3Mat
{
	protected:
		double		m_dValues[9];

	public:
				C3Mat();
				C3Mat(double inValues[9]);
				~C3Mat(){};
		C3Mat		Transpose();
		double		Determinant();
		C3Mat		Inverse();
		std::ostream&	Print(std::ostream& out) const;

		C3Vec		operator*(const C3Vec& v);
		C3Mat		operator*(const C3Mat& m);	
		C3Mat		operator*(const double& d);
		double&		operator()(unsigned int row, unsigned int col);
		double		operator()(unsigned int row, unsigned int col) const;

		friend std::ostream& operator<<(std::ostream& os, const C3Mat& m);
};


#endif //C3MAT_H_INCLUDED
