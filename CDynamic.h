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


#include "CApplication.h"

class CDynamic : public CApplication
{
	public:
				CDynamic();
				~CDynamic();
		void		Setup(int, const char **);

	protected:
		void		ParseInputFile();

	private:
		int		m_iNStep;
		double		m_dTimeStep;
//		Box		m_Box
};
