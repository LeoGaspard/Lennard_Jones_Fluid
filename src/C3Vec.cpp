///////////////////////////////////////////////////////
//NAME:			C3Vec.coo
//
//PURPOSE:		Definition of the C3Vec
//			class
//
//FUNCTIONS/OBJECTS:	C3Vec
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#include "C3Vec.h"

//Constructor
C3Vec::C3Vec(double inX, double inY, double inZ)
{
	m_dX = inX;
	m_dY = inY;
	m_dZ = inZ;
}//C3Vec

//Constructor
C3Vec::C3Vec()
{
	m_dX = 0.0;
	m_dY = 0.0;
	m_dZ = 0.0;
}//C3Vec

//Destructor
C3Vec::~C3Vec()
{
}//~C3Vec

// Returns the norm of the vector
double	C3Vec::Norm()
{
	return sqrt(m_dX*m_dX+m_dY*m_dY+m_dZ*m_dZ);
} //Norm

// Returns the square of the norm of the vector
double	C3Vec::Norm2()
{
	return m_dX*m_dX+m_dY*m_dY+m_dZ*m_dZ;
} //Norm2

// Returns the dot product between this and v
double	C3Vec::Dot(C3Vec v)
{
	return m_dX*v.GetX()+m_dY*v.GetY()+m_dZ*v.GetZ();
}//Dot

// Returns the cross product between this and v
C3Vec	C3Vec::Cross(C3Vec v)
{
	double x,y,z;

	x = m_dY*v.GetZ() - m_dZ*v.GetY();
	y = m_dZ*v.GetX() - m_dX*v.GetZ();
	z = m_dX*v.GetY() - m_dY*v.GetX();

	return C3Vec(x,y,z);
} // Cross

std::ostream&	C3Vec::Print(std::ostream& out) const
{
	out << boost::format("% 8.6e     % 8.6e     % 8.6e")%m_dX%m_dY%m_dZ;
	return out;
}

// Overloads the + operator between two C3Vec
C3Vec	C3Vec::operator+(const C3Vec& v)
{
	double x,y,z;

	x = this->m_dX + v.m_dX;
	y = this->m_dY + v.m_dY;
	z = this->m_dZ + v.m_dZ;

	return C3Vec(x,y,z);
}//operator+

// Overloads the += operator between two C3Vec
void	C3Vec::operator+=(const C3Vec& v)
{
	this->m_dX += v.m_dX;
	this->m_dY += v.m_dY;
	this->m_dZ += v.m_dZ;
}//operator+=

// Overloads the - operator between two C3Vec
C3Vec	C3Vec::operator-(const C3Vec& v)
{
	double x,y,z;

	x = this->m_dX - v.m_dX;
	y = this->m_dY - v.m_dY;
	z = this->m_dZ - v.m_dZ;

	return C3Vec(x,y,z);
}//operator-

// Overloads the -= operator between two C3Vec
void	C3Vec::operator-=(const C3Vec& v)
{
	this->m_dX -= v.m_dX;
	this->m_dY -= v.m_dY;
	this->m_dZ -= v.m_dZ;
}//operator-=

// Overloads the * operator between a C3Vec and a scalar
C3Vec	C3Vec::operator*(const double& d)
{
	double x,y,z;

	x = this->m_dX*d;
	y = this->m_dY*d;
	z = this->m_dZ*d;

	return C3Vec(x,y,z);
}//operator*

// Overloads the *= between a C3Vec and a scalar
void	C3Vec::operator*=(const double& d)
{
	this->m_dX *= d;
	this->m_dY *= d;
	this->m_dZ *= d;
}//operator*=

// Overloads the / operator between a C3Vec and a scalar
C3Vec	C3Vec::operator/(const double& d)
{
	double x,y,z;

	x = this->m_dX/d;
	y = this->m_dY/d;
	z = this->m_dZ/d;

	return C3Vec(x,y,z);
}//operator*

// Overloads the /= between a C3Vec and a scalar
void	C3Vec::operator/=(const double& d)
{
	this->m_dX /= d;
	this->m_dY /= d;
	this->m_dZ /= d;
}//operator*=

// Overloads the printing operator for a C3Vec
std::ostream&	operator<<(std::ostream& os, const C3Vec& v)
{
	return v.Print(os);
}//operator<<
