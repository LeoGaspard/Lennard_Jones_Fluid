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

#include "C3Mat.h"

//Constructor
C3Mat::C3Mat()
{
	for(int i=0;i<9;i++)
	{
		m_dValues[i] = 0.0;
	}
} //C3Mat

//Constructor
C3Mat::C3Mat(double inValues[9])
{
	for(int i=0;i<9;i++)
	{
		m_dValues[i] = inValues[i];
	}
} //C3Mat

//Returns a C3Mat object containing the transpose of the matrix
C3Mat	C3Mat::Transpose()
{
	C3Mat m;

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			m(i,j) = m_dValues[3*i+j];						
		}
	}

	return m;
} //Transpose

//Returns the determinant of the matrix
double	C3Mat::Determinant()
{
	double det(0.0);

	det += m_dValues[0] * (m_dValues[4]*m_dValues[8]-m_dValues[5]*m_dValues[7]);
	det -= m_dValues[1] * (m_dValues[3]*m_dValues[8]-m_dValues[5]*m_dValues[6]);
	det += m_dValues[2] * (m_dValues[3]*m_dValues[7]-m_dValues[4]*m_dValues[6]);

	return det;
}

//Returns the inverse of the matrix using the algorithm described in
C3Mat	C3Mat::Inverse()
{
	C3Mat m;
	double det = this->Determinant();

	if(det == 0)
	{
		return m;
	}

	m(0,0) = m_dValues[4]*m_dValues[8]-m_dValues[5]*m_dValues[7]; 
	m(0,1) = m_dValues[2]*m_dValues[7]-m_dValues[1]*m_dValues[8]; 
	m(0,2) = m_dValues[1]*m_dValues[5]-m_dValues[2]*m_dValues[4]; 
	m(1,0) = m_dValues[6]*m_dValues[5]-m_dValues[3]*m_dValues[8]; 
	m(1,1) = m_dValues[0]*m_dValues[8]-m_dValues[2]*m_dValues[6]; 
	m(1,2) = m_dValues[2]*m_dValues[3]-m_dValues[0]*m_dValues[5]; 
	m(2,0) = m_dValues[3]*m_dValues[7]-m_dValues[4]*m_dValues[6]; 
	m(2,1) = m_dValues[1]*m_dValues[6]-m_dValues[0]*m_dValues[7]; 
	m(2,2) = m_dValues[0]*m_dValues[4]-m_dValues[1]*m_dValues[3]; 

	return m*(1/det);
}

//Prints the matrix to the ostream 
std::ostream&	C3Mat::Print(std::ostream& out) const
{
	out << boost::format("% 8.6f     % 8.6f     % 8.6f")%m_dValues[0]%m_dValues[1]%m_dValues[2] << std::endl;
	out << boost::format("% 8.6f     % 8.6f     % 8.6f")%m_dValues[3]%m_dValues[4]%m_dValues[5] << std::endl;
	out << boost::format("% 8.6f     % 8.6f     % 8.6f")%m_dValues[6]%m_dValues[7]%m_dValues[8];

	return out;
}//Print

//Overloads the * operator for a C3Mat and a C3Vec
C3Vec	C3Mat::operator*(const C3Vec& v)
{
	double x,y,z;

	x = m_dValues[0]*v.m_dX+m_dValues[1]*v.m_dY+m_dValues[2]*v.m_dZ;
	y = m_dValues[3]*v.m_dX+m_dValues[4]*v.m_dY+m_dValues[5]*v.m_dZ;
	z = m_dValues[6]*v.m_dX+m_dValues[7]*v.m_dY+m_dValues[8]*v.m_dZ;
	return C3Vec(x,y,z);
}//operator*

//Overloads the * operator between two C3Mat
C3Mat	C3Mat::operator*(const C3Mat& m)
{
	C3Mat res;

	res(0,0) = m_dValues[0]*m(0,0) + m_dValues[1]*m(1,0) + m_dValues[2]*m(2,0); 
	res(0,1) = m_dValues[0]*m(0,1) + m_dValues[1]*m(1,1) + m_dValues[2]*m(2,1); 
	res(0,2) = m_dValues[0]*m(0,2) + m_dValues[1]*m(1,2) + m_dValues[2]*m(2,2); 
	res(1,0) = m_dValues[3]*m(0,0) + m_dValues[4]*m(1,0) + m_dValues[5]*m(2,0); 
	res(1,1) = m_dValues[3]*m(0,1) + m_dValues[4]*m(1,1) + m_dValues[5]*m(2,1); 
	res(1,2) = m_dValues[3]*m(0,2) + m_dValues[4]*m(1,2) + m_dValues[5]*m(2,2); 
	res(2,0) = m_dValues[6]*m(0,0) + m_dValues[7]*m(1,0) + m_dValues[8]*m(2,0); 
	res(2,1) = m_dValues[6]*m(0,1) + m_dValues[7]*m(1,1) + m_dValues[8]*m(2,1); 
	res(2,2) = m_dValues[6]*m(0,2) + m_dValues[7]*m(1,2) + m_dValues[8]*m(2,2);  

	return res;
} //operator*

//Overloads the * operator between a C3Mat and a double
C3Mat	C3Mat::operator*(const double& d)
{
	C3Mat res;	
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			res(i,j) = d*m_dValues[3*i+j];
		}
	}	

	return res;
}//operator*

//Overloads the () operator to get the coefficient
double&	C3Mat::operator()(unsigned int row, unsigned int col) 
{
	if(row>2 or col>2)
	{
		throw std::out_of_range("");
	}
	return m_dValues[3*row+col];			
}//operator()

//Overloads the () operator to get the coefficient
double	C3Mat::operator()(unsigned int row, unsigned int col) const
{
	if(row>2 or col>2)
	{
		throw std::out_of_range("");
	}
	return m_dValues[3*row+col];
}

//Overloads the printing operator for a C3Mat
std::ostream&	operator<<(std::ostream& os, const C3Mat& m)
{
	return m.Print(os);
}//operator<<






