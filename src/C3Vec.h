///////////////////////////////////////////////////////
//NAME:			C3Vec.h
//
//PURPOSE:		Definition of the C3Vec
//			class
//
//FUNCTIONS/OBJECTS:	C3Vec
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#ifndef C3VEC_H_INCLUDED
#define C3VEC_H_INCLUDED

#include <iostream>
#include <boost/format.hpp>
#include <math.h>


class C3Mat;

class C3Vec
{
	friend class C3Mat;
	protected:
		double		m_dX,m_dY,m_dZ;
	public:
				C3Vec(double inX, double inY, double inZ);
				C3Vec();
				~C3Vec();
		double		Norm();
		double		Norm2();
		double		Dot(C3Vec v);
		C3Vec		Cross(C3Vec v);
		std::ostream&	Print(std::ostream& out) const;

				// Getters
		double		GetX(){return m_dX;};
		double		GetY(){return m_dY;};
		double		GetZ(){return m_dZ;};

				// Setters
		void		SetX(double inX){m_dX = inX;};
		void		SetY(double inY){m_dY = inY;};
		void		SetZ(double inZ){m_dZ = inZ;};
		
				// Operators
		C3Vec		operator+(const C3Vec& v);
		void		operator+=(const C3Vec& v);
		C3Vec 		operator-(const C3Vec& v);
		void		operator-=(const C3Vec& v);
		C3Vec		operator*(const double& d);
		void		operator*=(const double& d);
		C3Vec		operator/(const double& d);
		void		operator/=(const double& d);

		friend std::ostream& operator<<(std::ostream& os, const C3Vec& v);
		friend C3Mat operator*(const C3Vec& v);




	private:

};

#endif // C3VEC_H_INCLUDED
