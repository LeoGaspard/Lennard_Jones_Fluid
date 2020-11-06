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
				CDynamic();
				~CDynamic();
		void		Setup(int, const char **);


		// Output methods
		void		OutHeader();

	protected:
		void		ParseInputFile();

	private:
		int		m_iNStep;
		double		m_dTimeStep;
		CBox		m_Box;
};
