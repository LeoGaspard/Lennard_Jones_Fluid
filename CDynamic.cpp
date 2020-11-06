///////////////////////////////////////////////////////
//NAME:			CDynamic.cpp
//
//PURPOSE:		Definition of the CDynamic
//			class
//
//FUNCTIONS/OBJECTS:	CDynamic
//
//AUTHOR:		Léo Gaspard
///////////////////////////////////////////////////////


#include "CDynamic.h"

//Constructor
CDynamic::CDynamic()
{
} //CDynamic

//Destructor
CDynamic::~CDynamic()
{

}//~CDynamic

// This setup the dynamic simulation using the provided input
void CDynamic::Setup(int argc, const char * argv[])
{
	this->ParseCommandLineOptions(argc, argv);
	//DEBUG
	std::cout << "Input  file       : " << m_strInputFile << std::endl;
	std::cout << "Output file       : " << m_strOutputFile << std::endl;
	std::cout << "Number of threads : " << m_iNumberThreads << std::endl << std::endl;
	//END DEBUG
	m_streamOutput.open(m_strOutputFile,std::ios::out);
	if(m_streamOutput.is_open())
	{
		this->OutHeader();
		this->ParseInputFile();
	}
	else
	{
		std::cerr << "Fatal Error : Can't open the output file : " << m_strOutputFile;
	}

}//Setup

//Parsing and checking the validity of the input file
void CDynamic::ParseInputFile()
{
	// These are the possible input variables
	double 					dA(0.0), dB(0.0), dC(0.0), dAlpha(90.0), dBeta(90.0), dGamma(90.0);
	std::vector<unsigned int> 		iNatom;
	std::vector<double>			dSigma, dEpsilon;
	
	//This regex will remove blank spaces before and after the first and last char
	std::regex trim("^\\s+|\\s+$");

	std::ifstream input(m_strInputFile);
	if(input.is_open())
	{
		// The input file has been successfully opened
		
		m_streamOutput << "Reading the input file : " << m_strInputFile << std::endl << std::endl;

		for(std::string strLine; getline(input,strLine); )
		{
			// This trims the line
			strLine = regex_replace(strLine,trim,"");
			// This removes the comments (starting with "!")
			strLine = strLine.substr(0,strLine.find("!"));

			if(strLine.size()>0)
			{
				m_streamOutput << strLine << std::endl; 
			}

			std::istringstream 	streamLine(strLine);
			std::string 		strFirstParam;

			std::getline(streamLine,strFirstParam,' ');

			if(strcmp(strFirstParam.data(),"Box")==0)
			{

				for(std::string strBoxParam; strcmp(strBoxParam.data(),"End")!=0; getline(input, strBoxParam))
				{
					strBoxParam = regex_replace(strBoxParam,trim,""); 
					strBoxParam = strBoxParam.substr(0,strBoxParam.find("!"));

					std::istringstream streamBoxParam(strBoxParam);
					std::string strFirstParam;

					std::getline(streamBoxParam,strFirstParam,' ');

					if(strcmp(strFirstParam.data(),"a")==0)
					{
						std::string strValue;
						std::getline(streamBoxParam, strValue, ' ');

						dA = std::stod(strValue);

						// We check if there was a unit specification
						std::getline(streamBoxParam, strValue, ' ');
						if(strcmp(strValue.data(),"A")==0)
						{
							dA *= ANGSTROM_TO_BOHR;
						}
					}
					else if(strcmp(strFirstParam.data(),"b")==0)
					{
						std::string strValue;
						std::getline(streamBoxParam, strValue, ' ');

						dB = std::stod(strValue);

						// We check if there was a unit specification
						std::getline(streamBoxParam, strValue, ' ');
						if(strcmp(strValue.data(),"A")==0)
						{
							dB *= ANGSTROM_TO_BOHR;
						}
					}
					else if(strcmp(strFirstParam.data(),"c")==0)
					{
						std::string strValue;
						std::getline(streamBoxParam, strValue, ' ');

						dC = std::stod(strValue);

						// We check if there was a unit specification
						std::getline(streamBoxParam, strValue, ' ');
						if(strcmp(strValue.data(),"A")==0)
						{
							dC *= ANGSTROM_TO_BOHR;
						}
					}
					else if(strcmp(strFirstParam.data(),"alpha")==0)
					{
						std::string strValue;
						std::getline(streamBoxParam, strValue, ' ');

						dAlpha = std::stod(strValue);
					}
					else if(strcmp(strFirstParam.data(),"beta")==0)
					{
						std::string strValue;
						std::getline(streamBoxParam, strValue, ' ');

						dBeta = std::stod(strValue);
					}
					else if(strcmp(strFirstParam.data(),"gamma")==0)
					{
						std::string strValue;
						std::getline(streamBoxParam, strValue, ' ');

						dGamma = std::stod(strValue);
					}
					else if(strcmp(strBoxParam.data(),"")!=0)
					{
						std::cerr << "Fatal Error : Bad input with line " << std::endl << strBoxParam << std::endl;
						exit(1);
					}
				}

				if(dAlpha+dBeta+dGamma > 360)
				{
					std::cerr << "Fatal Error : Bad angle value for the box\nCheck https://doi.org/10.1107/S0108767310044296 for more informations" << std::endl;
					exit(1);
				}
				if(dAlpha-dBeta+dGamma > 360)
				{
					std::cerr << "Fatal Error : Bad angle value for the box\nCheck https://doi.org/10.1107/S0108767310044296 for more informations" << std::endl;
					exit(1);
				}
				if(dAlpha+dBeta-dGamma > 360)
				{
					std::cerr << "Fatal Error : Bad angle value for the box\nCheck https://doi.org/10.1107/S0108767310044296 for more informations" << std::endl;
					exit(1);
				}
				if(-dAlpha+dBeta+dGamma > 360)
				{
					std::cerr << "Fatal Error : Bad angle value for the box\nCheck https://doi.org/10.1107/S0108767310044296 for more informations" << std::endl;
					exit(1);
				}


				m_Box = CBox(dA,dB,dC,dAlpha,dBeta,dGamma);
			}
		}
		input.close();
		
		m_streamOutput << std::endl << "Input file is consistent" << std::endl;
		m_streamOutput << std::string(150,'=') << std::endl;
		m_streamOutput << boost::format("%=150s")%"PRINTING THE PARAMETERS AT THE BEGINNING OF THE COMPUTATION" << std::endl;
		m_Box.OutBoxParam(m_streamOutput);
		m_streamOutput << std::string(150,'-') << std::endl;
		m_streamOutput << boost::format("%-25s %i")%"Number of threads"%m_iNumberThreads << std::endl;
		m_streamOutput << boost::format("%-25s %i")%"Number of steps"%m_iNStep << std::endl;
		m_streamOutput << boost::format("%-25s %d picoseconds")%"Timestep"%m_iNStep << std::endl;

		m_streamOutput << std::string(150,'=') << std::endl;
	}
	else
	{
		// The input file can not be opened
		std::cerr << "Fatal Error : The input file " << m_strInputFile << " can not be opened" << std::endl << "The program will stop." << std::endl;
		exit(1);
	}
} //m_ParseInputFile

// This writes the header (TCCM logo + code name) 
// at the top of the output file
void	CDynamic::OutHeader()
{
	std::string strLogo("                                                                                                           .:/+/`                                    \n                                                                                                        `:yhdddds                                     \n                                                                                                       `sdddyhddd                                     \n                                                                                                       oddho+dddy                                     \n                                                                                                      :ddds/yddd:                                     \n                                                                                                     `sddh/odddo`                                     \n                                                                                                     /dddo+hddh.                                      \n                                                                                                    `hddy/yddd:                                       \n                                                                                                    /ddhosdddo`                                       \n                                                                                                   .yddyohddy`                                        \n                                                                                                   +dddohddh.                                        \n \n                                                                  ```.``                         .hddyhddh/                                   ```    \n                                                               `-/osyhhyyo:`                      +dddhddd+`                             `.-/ossss+:  \n                                                            `:sydddhhhhhddd+                     .hdddddds`                           `-+shddhy+:.`   \n                                                          .ohddhys++////sddd-                    +ddddddy.                          .+ydddhs/.        \n                                                        .ohddhs+/:::::::/hdd+                   -hdddddh-                         :shdddho-           \n                                                      `/hddho////:::::::/hddo                   odddddd/                        :sddddh+.             \n                                                  `.``sdddhsssyyyso+//::/hdd+                  .ddddddo                       -sddddh+.               \n                                                 -/`.sdddddddhhhhdddhyo/+hdd:                  oddddds.                     `ohddddy-                 \n                                                /o``sddhyo/:--..--:ohddhyddh`                 -hddddh.                     :hddddd/`                  \n                                               /h` -y+-``.-:///:-`  `sddddd+                  oddddd:                    `odddddy-                   \n                                              -ds   `-+shddddddddh/  .ydddy.                 -ddddd+                    -ydddddo`                     \n                                             `ydds+oso/-.....-/yddh-  oddd+                  sdddds`                   /hddddh/`                      \n       `..```````                             ohh+-`            odd+  oddh`                 -hdddh.                  `odddddh:                        \n      .hdddddddhhyyyss++/:-.``                 `                -hd: `sdd/                 `sdddh:                  .ydddddy.                         \n      `hdd/:////+oossyyhddddddhso+/-.`                          .hh. -hds`                 :dddd+                  :hdddhhs`                          \n       ydd.     ://-.````..:/+osyhdddhhs+/-`                    :ds  +dh:                  ydddy`                `+hdhhhh+`                           \n       odd:      `.+hyso/`      ``.-/+oyhddhyo/-`               +d- `hds                  :dddh-                `shhhhhh+                             \n       +dd/       .hd:``..:::.         ``-:+shhdhy+:`           hh` -dd.                 `yddd/                -yhhhhhy/                              \n       :dd+      `ods`  /y/ .-   `.-.        `-/oyddhs/`       -d+  odo                  /ddds                :hhhhhhh:                               \n       -ddo      :dd- `oh/     `+ys-./           `./sys.       oh: .yd:                 .hddh.              `+hhhhhyy:                                \n       .dds     `ydo  ods     -ydo`     .ss-         ``       .yy` /dh-                 +ddh:              `shhhhhyy-                                 \n       .hds     +dh. -dh-    .yds`     `sdd+  `:::            /d+ `sdd:                -hddo              -shhyyyys-                                  \n       `yds`   .yd+ `sdo`    +dd.      +dhs/ -o:+:           `sd. .ddd+               `sddh`            `/yhhyyyyo-                                   \n       `ydy`   +dh. .yh:     hds      .hh:+:/s.-o            .dy  /dddh-              :ddh:            .ohhysyyyo.                                    \n       `sdy`  .dd/  .hy.    .hd/      odo`oos.`s:            +d+  yddddh:            `hdd+           `/yhhssyyy/`                                     \n       `sdy`  :hh`  `yy` .- .hh-     -hh- os` /s`            yd:  hddhddho-`         +ddy`         `:yhhyosyys:                                       \n       `ody.  `..    -s+:o- `yh.  .  +d+  `` `y-            .hd- `hddsohddhy/-`     :hdd:       `./yhhyooyyyo-                                        \n       `ody.          `..`   -s:`-o``hd.     /o             :hd-  hddy//oyhddhyo+::ohddd:````.:+yhdhyo+syyy+`                                         \n       `ody.                  `--:` .y+     .y-             /dd:  oddh/::/+oshhhddddhhddhysyyhddhyo+/syyys-                                           \n       `odh.  -.`        `.-`               :s`             :ddo  -hdds/::::://+oooo+/osyyyyyso+//+syhys/`                                            \n       `odh.  `````   ``````                .-              .hdd.  /hdds/:::::::::::::://////:/+oyhhhs/`                                              \n       `ody.     `--.``                                      sdds.  :hddho//::::::::::::///+osyhhhys:`                                                \n       `ody. `...``.`                                        -hddy-  .ohddhyo+//::///+oosyhhddhhs/.`                                                  \n       `ody. ```    .           `...`                         :hddh/`  ./shdddhs++syhddddhhhs+:.                                                      \n       `ody`        ``         `/. ..    `                     .ydddh+-`  .:oyddddddhso+:-.                                                           \n       `sds`        `..`      `:.`.- .. .-`                     `:yddddy+:`  `+dddd+`                                                                 \n       `yds`        ```      `-. -:  --`-`- ````    `` `          `-+yhdddho  .ddd+                                                                   \n       `hds                      .`  .-:`. `.....  `:``.`   `         ./shds  oddy.                                                                   \n       .dd+                          `-      ```` `/.-` .- -.`           oh- :hdd:                                                                    \n       -dd/        ```                            --`.` :.-`..          .h+ `ydds                                                                     \n       :dd-           ```                               `-.``          `oy  +ddh-                                                                     \n       +dd`              ````                            `             :h- -ddd+                                                                      \n       ydh                  ``.`      ..                              `h/ .yddh.                                                                      \n      `hdo    `` `` `. `- `-   ``..```/-                             `os` oddd+                                                                       \n      :dd:      ` `..``..``-`-``   ` +/..``                          /h. :dddh.                                                                       \n      odh.           ..`.`..`-...`` -/...  ``````````               .h: -hddd+                                                                        \n     .hdo     ``.`.````.. -` .  `  `/:-`             `````         `s+ `ydddh.                                                                        \n     odh-            ````....```` `::``    -   `                   +s  odddd+                                                                         \n    -ddh++///:--.``            ```-+.:.        `  ``              :y. /ddddd`                                                                         \n    +ddddddddddhhhyso+/:.`       `/-  `..``       `   `   `      .y: -hdddd+                                                                          \n    `......--:/++osyhdddhhys+:-.  `       `..``                 `s/ .hddddh-                                                                          \n                  ```.-:+oyhhdhhyo/-.`        `````             oo``sdddddo`                                                                          \n                          `.-:+syhdhhyo/-.        ``           :h- :hddddh.                                                                           \n                                `.:+syhdhhs+:.`               -hy` /ddddh:                                                                            \n                                     `.:+shddhy+:.`           ody` .syy+-                                                                             \n                                         `.-+syddhyo:.`      `sdh/` ```                                                                               \n                                              `-/syhdhy:      +dddyso+`      -so.                                                                     \n                                                  `-+hdy`     `+syso:`       odd:                                                                     \n                                                    `ydy.                   `odd/                                                                     \n                                                    `oddo.                 `oddh-                                                                     \n                                                     `+hdh/     `....``   :ydds.                                                                      \n                                                       -yddo-:osyhhhhhyo/+ddh+`                                                                       \n                                                        `ohdhdhhssooooshdddy-                                                                         \n                                                          :hddy+:::::/+ydh+`                                                                          \n                                                           .sddh+/::/ohdy:                                                                            \n                                                             /hdds//yddo.                                                                             \n                                                              -sddhhdh/                                                                               \n                                                               `+hdds.                                                                                \n                                                                 -+:`                                                                               ");

	m_streamOutput << strLogo << std::endl << std::endl << std::endl;
	m_streamOutput << boost::format("%=150s")%"Lennard-Jones fluid simulation" << std::endl;
	m_streamOutput << boost::format("%=150s")%"TCCM M2 Computational Chemistry Programming Project course" << std::endl << std::endl;
	m_streamOutput << boost::format("%=150s")%"Léo Gaspard" << std::endl;
	m_streamOutput << std::endl << std::string(150,'=') << std::endl;
}//OutHeader

