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
		


			// Thermostats
		void		Berendsen();


		// Output methods
		void		OutHeader();

	protected:
		void		ParseInputFile();

	private:
		CBox		m_Box;

		//MD attributes
		double		m_dTimeStep, m_dInitTemperature, m_dNeighbor, m_dCutoff;
		int		m_iNStep;

		//Thermostats attributes
		double		m_dBerendsenTau;
};
