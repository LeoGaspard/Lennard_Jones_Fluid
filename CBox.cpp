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
}//CBox

//Destructor
CBox::~CBox()
{

}//~CBox

// Computes the volume and vectors of the box
void	CBox::Setup()
{
	double ca = cos(m_dAlpha*DEGREE_TO_RADIAN);
	double cb = cos(m_dBeta*DEGREE_TO_RADIAN);
	double cc = cos(m_dGamma*DEGREE_TO_RADIAN);
	double sc = sin(m_dGamma*DEGREE_TO_RADIAN);

	//Building the unit cell vectors in the cartesian basis
	m_H(0,0) = m_dA;
	m_H(1,0) = 0.0;
	m_H(2,0) = 0.0;
	m_H(0,1) = m_dB*cc;
	m_H(1,1) = m_dB*sc;
	m_H(2,1) = 0.0;
	m_H(0,2) = m_dC*cb;
	m_H(1,2) = m_dC*(ca-cb*cc)/sc;
	m_H(2,2) = m_dC*sqrt(1-cb*cb-(ca-cb*cc)*(ca-cb*cc)/(sc*sc));

	m_dVolume = m_H.Determinant();
}//Setup

// Generates random positions separated by inDMin inside the cell
// The cell is placed inside of a grid of cubes of great diagonal = inDMin such
// that each cube can contain one and only one point.
// For each trial position, only the surrounding cubes are tested
void	CBox::InitPosFromRandomDistribution(double inDMin)
{
	unsigned int					iNPoints(m_vAtomList.size());
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
				double x = i*m_H(0,0) + j*m_H(0,1) + k*m_H(0,2);
				double y = i*m_H(1,0) + j*m_H(1,1) + k*m_H(1,2);
				double z = i*m_H(2,0) + j*m_H(2,1) + k*m_H(2,2);


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
	while(iNAccepted < iNPoints)
	{
		bool	bAccepted = true;

		// Generate the trial fractional coordinates
		CPos	p1(seed(eng),seed(eng),seed(eng));

		// Express the trial coordinates in the cartesian basis
		p1 = m_H*p1;

		// Find the grid position of this point
		int iGridX = static_cast<int>(floor(p1.GetX()/dCellSize));
		int iGridY = static_cast<int>(floor(p1.GetY()/dCellSize));
		int iGridZ = static_cast<int>(floor(p1.GetZ()/dCellSize));

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
		std::stringstream errMsg;
		errMsg << "Requested access to atom n°" << n+1 << " where there are only " << m_vAtomList.size();
		throw std::length_error(errMsg.str());
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

		p = m_H.Inverse()*p;

		double u = p.GetX();
		double v = p.GetY();
		double w = p.GetZ();

		p = m_H.Inverse()*p;

		// Apply PBC

		u = (u<0 || u>1) ? u-floor(u) : u;
		v = (v<0 || v>1) ? v-floor(v) : v;
		w = (w<0 || w>1) ? w-floor(w) : w;

		p.SetX(u);
		p.SetY(v);
		p.SetZ(w);

		// Converting back to cartesian coordinates
		p= m_H*p;

		m_vAtomList[i].SetPos(p);
	}
}//Wrap

void	CBox::AddAtoms(unsigned int inNumber, double inSigma, double inEpsilon,double inMass)
{
	for(unsigned int i=0; i<inNumber; i++)
	{
		m_vAtomList.push_back(CAtom(inMass,inSigma,inEpsilon));
	}
}

// Initialise the speeds with random values
void	CBox::InitSpeedRandom(double inTemperature)
{
	std::random_device			rd;
	std::default_random_engine		eng(rd());
	std::uniform_real_distribution<double>	seed(0,1);

	for(unsigned int i=0; i<m_vAtomList.size(); i++)
	{
		double dX = seed(eng);
		double dY = seed(eng);
		double dZ = seed(eng);

		CSpeed s(dX,dY,dZ);

		// scaling the speed using <V>²=3Kb<T>/m
		double dScaling = (s.Norm()*BOHR_TO_ANGSTROM*ANGSTROM_TO_M)/sqrt(3*KB*inTemperature/(m_vAtomList[i].GetMass()*DAL_TO_KG));
		s /= dScaling;
		m_vAtomList[i].SetSpeed(s);
	}
}//InitSpeedRandom

double	CBox::ComputeTemperature()
{
	double dT(0.0);

	for(unsigned int i=0; i<m_vAtomList.size();i++)
	{
		m_vAtomList[i].ComputeKineticEnergy();
		
		dT += 2*m_vAtomList[i].GetKineticEnergy()/(3*KB);
	}

	dT /= m_vAtomList.size();

	return dT;
}//ComputeTemperature

void	CBox::ComputeForces()
{
	double dPotential(0.0);

	#pragma omp parallel for reduction(+:dPotential)
	for(unsigned int i=0;i<(m_vAtomList.size()-1);i++)
	{
		for(unsigned int k=0;k<m_vNeighborList[i].size();k++)
		{
			unsigned int j = m_vNeighborList[i][k];
			C3Vec ri = m_vAtomList[i].GetPos();
			C3Vec rj = m_vAtomList[j].GetPos();

			// Algorithm for the minimum image convention in
			// Appendix B, Eq. B.9. 
			// M. E. Tuckerman. Statistical Mechanics : Theory and Molecular Simulation
			// Oxford University Press, Oxford, UK, 2010
			ri = m_H.Inverse()*ri;
			rj = m_H.Inverse()*rj;

			C3Vec rij = ri-rj;

			rij.SetX(rij.GetX()-round(rij.GetX()));
			rij.SetX(rij.GetY()-round(rij.GetY()));
			rij.SetX(rij.GetZ()-round(rij.GetZ()));

			rij = m_H*rij;

			// End of the algorithm

			double r2 = rij.Norm2();
			double sigma(0.0), epsilon(0.0);

			if(m_vAtomList[i].GetSigma() == m_vAtomList[j].GetSigma())
			{
				sigma = m_vAtomList[i].GetSigma();
			}
			else
			{
				sigma = (m_vAtomList[i].GetSigma()+m_vAtomList[j].GetSigma())/2.0;
			}
			if(m_vAtomList[i].GetEpsilon() == m_vAtomList[j].GetEpsilon())
			{
				epsilon = m_vAtomList[i].GetEpsilon();
			}
			else
			{
				epsilon = sqrt(m_vAtomList[i].GetEpsilon()*m_vAtomList[j].GetEpsilon());
			}

			double r6 = r2*r2*r2;
			double s6 = sigma*sigma*sigma*sigma*sigma*sigma;

			// V = 4E[s^12/r^12-s^6/r^6]
			dPotential += 4*epsilon*((s6*s6)/(r6*r6)-s6/r6);

			CForce f;

			// dV/dx1 = (x1-x2)*E*[24*s^6/r^8-48*s^12/r^14]
			f.SetX((ri.GetX()-rj.GetX())*epsilon*(24*s6/(r6*r2)-46*s6*s6/(r6*r6*r2)));
			f.SetY((ri.GetY()-rj.GetY())*epsilon*(24*s6/(r6*r2)-46*s6*s6/(r6*r6*r2)));
			f.SetZ((ri.GetZ()-rj.GetZ())*epsilon*(24*s6/(r6*r2)-46*s6*s6/(r6*r6*r2)));

			#pragma omp critical
			m_vAtomList[i].AddForce(f);

			f = CForce() - f;

			#pragma omp critical
			m_vAtomList[j].AddForce(f);
		}
	}
}

// Builds the neighbor list in the box
void	CBox::NeighborList(double cutoff, double neighbor)
{
	m_vNeighborList.clear();
	m_vPosList.clear();

	m_vNeighborList.resize(m_vAtomList.size());

	double 	r2 = (cutoff+neighbor)*(cutoff+neighbor);
	#pragma omp parallel for schedule(guided)
	for(unsigned int i=0; i<m_vAtomList.size()-1; i++)
	{
		for(unsigned int j=i+1; j<m_vAtomList.size();j++)
		{
			CPos ri = m_vAtomList[i].GetPos();
			CPos rj = m_vAtomList[j].GetPos();

			// Algorithm for the minimum image convention in
			// Appendix B, Eq. B.9. 
			// M. E. Tuckerman. Statistical Mechanics : Theory and Molecular Simulation
			// Oxford University Press, Oxford, UK, 2010
			ri = m_H.Inverse()*ri;
			rj = m_H.Inverse()*rj;

			C3Vec rij = ri-rj;

			rij.SetX(rij.GetX()-round(rij.GetX()));
			rij.SetX(rij.GetY()-round(rij.GetY()));
			rij.SetX(rij.GetZ()-round(rij.GetZ()));

			rij = m_H*rij;

			// End of the algorithm
			
			double	d2 = rij.Norm2();


			if(d2<r2)
			{
				#pragma omp critical
				{
					m_vNeighborList[i].push_back(j);
				}
			}
		}
	}
	for(unsigned int i=0;i<m_vAtomList.size();i++)
	{
		m_vPosList.push_back(m_vAtomList[i].GetPos());
	}
}//NeighborList

//Check if the maximum displacement from the m_vPosList is 
//greater than neighbor/2, return true if it is
bool	CBox::CheckNeighborList(double neighbor)
{
	bool check = true;
	#pragma omp parallel for
	for(unsigned int i=0;i<m_vAtomList.size();i++)
	{
		CPos ri = m_vAtomList[i].GetPos();
		CPos rj = m_vPosList[i];

		C3Vec rij = ri - rj;

		if(rij.Norm2() > neighbor/2)
		{
			std::cout << rij.Norm2() << std::endl;
			check = false;
		}
	}
	return check;
}


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
	f << std::endl << m_H << std::endl;
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


