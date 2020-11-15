#include <string>
#include <exception>

class CError: public std::exception
{
	public:
		CError(std::string const& inMessage="") throw() : m_sMessage(inMessage) {}
		virtual ~CError() {}

		virtual const char* what() const throw()
		{
			return m_sMessage.c_str();
		}

	private:
		std::string	m_sMessage;
};
