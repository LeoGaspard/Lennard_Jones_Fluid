///////////////////////////////////////////////////////
//NAME:			CDynamic.h
//
//PURPOSE:		Definition of the CDynamic
//			class
//
//FUNCTIONS/OBJECTS:	CDynamic
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////

#include <boost/format.hpp>

#include "CApplication.h"
#include "CBox.h"

class CDynamic : public CApplication
{
	public:
		// Methods
				CDynamic();
				~CDynamic();
		void		Setup(int, const char **);
		void		Run();
		


			// Thermostats
		double		Berendsen();


		// Output methods
		void		OutHeader();

	protected:
		void		ParseInputFile();

	private:
		CBox		m_Box;

		//MD attributes
		double		m_dTimeStep, m_dInitTemperature, m_dNeighbor, m_dCutoff, m_dRadialStep, m_dRadialMax;
		unsigned int	m_iNStep, m_iNEquilibrationStep, m_iThermostatStep, m_iPrintEvery;

		//Thermostats attributes
		double		m_dBerendsenTau;
};
