///////////////////////////////////////////////////////
//NAME:			CBox.cpp
//
//PURPOSE:		Definition of the CBox
//			class
//
//FUNCTIONS/OBJECTS:	CBox
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#include "CBox.h"

//Constructor
CBox::CBox(double inA, double inB, double inC, double inAlpha, double inBeta, double inGamma)
{
	m_dA = inA;
	m_dB = inB;
	m_dC = inC;
	m_dAlpha = inAlpha;
	m_dBeta = inBeta;
	m_dGamma = inGamma;

	double ca = cos(inAlpha*DEGREE_TO_RADIAN);
	double cb = cos(inBeta*DEGREE_TO_RADIAN);
	double cc = cos(inGamma*DEGREE_TO_RADIAN);
	double sc = sin(inGamma*DEGREE_TO_RADIAN);

	m_dVolume = inA*inB*inC*sqrt(1-ca*ca-cb*cb-cc*cc+2*ca*cb*cc);

	//Building the unit cell vectors in the cartesian basis
	
	// a = [a,0,0]
	double x = inA;
	double y = 0.0;
	double z = 0.0;

	m_AVec = C3Vec(x,y,z);

	// b = [b*cos(gamma),b*sin(gamma),0]
	x = inB*cc;
	y = inB*sc;
	z = 0.0;

	m_BVec = C3Vec(x,y,z);

	// c = [c*cos(beta),c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma),volume/(a*b*sin(gamma))]
	x = inC*cb;
	y = inC*(ca-cb*cc)/sc;
	z = m_dVolume/(inA*inB*sc);

	m_CVec = C3Vec(x,y,z);

	//Building the unit cell vectors in the fractionnal basis
	
	// u = [1/a,0,0]
	x = 1/inA;
	y = 0.0;
	z = 0.0;

	m_UVec = C3Vec(x,y,z);
	
	// v = [-cos(gamma)/(a*sin(gamma)),1/(b*sin(gamma)),0]
	x = -cc/(inA*sc);
	y = 1/(inB*sc);
	z = 0.0;

	m_VVec = C3Vec(x,y,z);
	
	// w = [b*c*(cos(alpha)*cos(gamma)-cos(beta))/(volume*sin(gamma)),a*c*(cos(beta)*cos(gamma)-cos(alpha))/(volume*sin(gamma)),a*b*sin(gamma)/volume]
	x = inB*inC*(ca*cc-cb)/(m_dVolume*sc);
	y = inA*inC*(cb*cc-ca)/(m_dVolume*sc);
	z = inA*inB*sc/m_dVolume;

	m_WVec = C3Vec(x,y,z);
}//CBox

//Destructor
CBox::~CBox()
{

}//~CBox

// Generates random positions separated by inDMin inside the cell
// The cell is placed inside of a grid of cubes of great diagonal = inDMin such
// that each cube can contain one and only one point.
// For each trial position, only the surrounding cubes are tested
void	CBox::InitPosFromRandomDistribution(unsigned int inNPoints, double inDMin)
{
	double 						r2 = inDMin*inDMin;
	double 						dCellSize = inDMin/sqrt(3); // The great diagonal of a cube of sides a is a*sqrt(3)
	std::vector<std::vector<std::vector<int>>>	vGrid;

	// In order to find the bounds of the grid, we need to find the minimum
	// and maximum values of x,y and z.
	
	double dXMin(0.0), dXMax(0.0), dYMin(0.0), dYMax(0.0), dZMin(0.0), dZMax(0.0);

	for(unsigned int i=0;i<2;i++)
	{
		for(unsigned int j=0;j<2;j++)
		{
			for(unsigned int k=0;k<2;k++)
			{
				double x = i*m_AVec.GetX() + j*m_BVec.GetX() + k*m_CVec.GetX();
				double y = i*m_AVec.GetY() + j*m_BVec.GetY() + k*m_CVec.GetY();
				double z = i*m_AVec.GetZ() + j*m_BVec.GetZ() + k*m_CVec.GetZ();

				dXMin = (x < dXMin) ? x : dXMin; 
				dXMax = (x > dXMax) ? x : dXMax; 
				dYMin = (y < dYMin) ? y : dYMin; 
				dYMax = (y > dYMax) ? y : dYMax; 
				dZMin = (z < dZMin) ? z : dZMin; 
				dZMax = (z > dZMax) ? z : dZMax; 
			}
		}
	}

	// Determining the number of division of the grid in each directions
	unsigned int iNDivX = static_cast<unsigned int>(ceil((dXMax-dXMin)/dCellSize));
	unsigned int iNDivY = static_cast<unsigned int>(ceil((dYMax-dYMin)/dCellSize));
	unsigned int iNDivZ = static_cast<unsigned int>(ceil((dZMax-dZMin)/dCellSize));
	
	vGrid.resize(iNDivX);
	for(unsigned int i=0;i<iNDivX;i++)
	{
		vGrid[i].resize(iNDivY);
		for(unsigned int j=0;j<iNDivY;j++)
		{
			vGrid[i][j].resize(iNDivZ);
			for(unsigned int k=0;k<iNDivZ;k++)
			{
				vGrid[i][j][k] = -1; // -1 means that this cube is empty
			}
		}
	}

	// Loop to generate the positions
	
	std::random_device			rd;
	std::default_random_engine		eng(rd());
	std::uniform_real_distribution<double>	seed(0,1);
	unsigned int				iNAccepted(0);
	while(iNAccepted < inNPoints)
	{
		bool	bAccepted = true;

		// Generate the trial fractional coordinates
		double	dTrialU = seed(eng);
		double	dTrialV = seed(eng);
		double	dTrialW = seed(eng);

		// Express the trial coordinates in the cartesian basis
		double dTrialX = dTrialU*m_AVec.GetX() + dTrialV*m_BVec.GetX() + dTrialW*m_CVec.GetX();
		double dTrialY = dTrialU*m_AVec.GetY() + dTrialV*m_BVec.GetY() + dTrialW*m_CVec.GetY();
		double dTrialZ = dTrialU*m_AVec.GetZ() + dTrialV*m_BVec.GetZ() + dTrialW*m_CVec.GetZ();


		CPos	p1(dTrialX,dTrialY,dTrialZ);

		// Find the grid position of this point
		int iGridX = static_cast<int>(floor(dTrialX/dCellSize));
		int iGridY = static_cast<int>(floor(dTrialY/dCellSize));
		int iGridZ = static_cast<int>(floor(dTrialZ/dCellSize));

		// Applying PBC
		iGridX = (iGridX < 0) ? iGridX+iNDivX : iGridX;
		iGridY = (iGridY < 0) ? iGridY+iNDivY : iGridY;
		iGridZ = (iGridZ < 0) ? iGridZ+iNDivZ : iGridZ;

		// Looping on the neighbouring cubes
		for(int x = iGridX-2; x<=iGridX+2 && bAccepted; x++)
		{
			for(int y = iGridY-2; y <=iGridY+2 && bAccepted; y++)
			{
				for(int z = iGridZ-2; z <=iGridZ+2 && bAccepted; z++)
				{
					// Applying PBC to x y and z
					unsigned int iGx = (x < 0) ? x+iNDivX : x;
					iGx = (iGx >= iNDivX) ? iGx-iNDivX : iGx; 
					unsigned int iGy = (y < 0) ? y+iNDivY : y;
					iGy = (iGy >= iNDivY) ? iGy-iNDivY : iGy; 
					unsigned int iGz = (z < 0) ? z+iNDivZ : z;
					iGz = (iGz >= iNDivZ) ? iGz-iNDivZ : iGz; 

					int index = vGrid[iGx][iGy][iGz];
					if(index != -1)
					{
						CPos p2 = m_vAtomList[index].GetPos();

						double d2 = p1.Distance2(p2);
						if(d2 < r2)
						{
							bAccepted = false;
						}
					}
				}
			}
		}
		if(bAccepted)
		{
			vGrid[iGridX][iGridY][iGridZ] = iNAccepted;
			m_vAtomList[iNAccepted].SetPos(p1);
			iNAccepted++;
		}
	}
}// InitPosFromRandomDistribution

CAtom&	CBox::GetAtom(unsigned int n)
{
	if(n < m_vAtomList.size())
	{
		return m_vAtomList[n];
	}
	else
	{
		std::cerr << "Fatal Error : Requested access to atom n°" << n+1 << " where there are only " << m_vAtomList.size() << std::endl;
		exit(1);
	}
}//GetAtom

// Wraps all the atoms inside the box using PBC
// This works by converting all the positions to fractional coordinates
// and then back again to cartesian coordinates
void	CBox::Wrap()
{
	for(unsigned int i=0; i<m_vAtomList.size(); i++)
	{
		CPos p = m_vAtomList[i].GetPos();

		double x = p.GetX();
		double y = p.GetY();
		double z = p.GetZ();

		// Transform x, y and z into fractional coordinates

		double u = x*m_UVec.GetX() + y*m_VVec.GetX() + z*m_WVec.GetX();
		double v = x*m_UVec.GetY() + y*m_VVec.GetY() + z*m_WVec.GetY();
		double w = x*m_UVec.GetZ() + y*m_VVec.GetZ() + z*m_WVec.GetZ();

		// Apply PBC

		u = (u<0 || u>1) ? u-floor(u) : u;
		v = (v<0 || v>1) ? v-floor(v) : v;
		w = (w<0 || w>1) ? w-floor(w) : w;

		// Converting back to cartesian coordinates

		x = u*m_AVec.GetX() + v*m_BVec.GetX() + w*m_CVec.GetX();
		y = u*m_AVec.GetY() + v*m_BVec.GetY() + w*m_CVec.GetY();
		z = u*m_AVec.GetZ() + v*m_BVec.GetZ() + w*m_CVec.GetZ();

		m_vAtomList[i].SetPos(CPos(x,y,z));
	}
}//Wrap

// Write the box parameters in the stream f
void	CBox::OutBoxParam(std::ofstream& f)
{
	// Box parameters
	f << std::string(150,'=') << std::endl;
	f << boost::format("%-150s")%"Characteristic lengths and angles of the box :" << std::endl << std::endl;
	f << boost::format("%-7s  %-10.8d bohr")%"a"%m_dA << std::endl;
	f << boost::format("%-7s  %-10.8d bohr")%"b"%m_dB << std::endl;
	f << boost::format("%-7s  %-10.8d bohr")%"c"%m_dC << std::endl;
	f << boost::format("%-7s  %-10.8d °")%"alpha"%m_dAlpha << std::endl;
	f << boost::format("%-7s  %-10.8d °")%"beta"%m_dBeta << std::endl;
	f << boost::format("%-7s  %-10.8d °")%"gamma"%m_dGamma << std::endl;
	f << boost::format("%-7s  %-10.8d bohr^3")%"volume"%m_dVolume << std::endl;

	// Box vectors
	f << std::string(150,'-') << std::endl;
	f << boost::format("%-150s")%"a b and c vectors in the cartesian basis :" << std::endl << std::endl;
	f << m_AVec << std::endl;
	f << m_BVec << std::endl;
	f << m_CVec << std::endl;
} // OutBoxParam

// Write the atomic coordinates in the stream f
void	CBox::OutAtomPos(std::ofstream& f)
{
	for(unsigned int i = 0; i<m_vAtomList.size();i++)
	{
		f << m_vAtomList[i].GetPos() << std::endl;
	}
}//OutAtomPos

// Write the atomic speeds in the stream f
void	CBox::OutAtomSpeed(std::ofstream& f)
{
	for(unsigned int i = 0; i<m_vAtomList.size();i++)
	{
		f << m_vAtomList[i].GetSpeed() << std::endl;
	}
}//OutAtomSpeed

// Write the forces vectors in the stream f
void	CBox::OutAtomForces(std::ofstream& f)
{
	for(unsigned int i = 0; i<m_vAtomList.size();i++)
	{
		f << m_vAtomList[i].GetForces() << std::endl;
	}
}//OutAtomForces
