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
	omp_set_num_threads(m_iNumberThreads);
	m_streamOutput.open(m_strOutputFile,std::ios::out);
	if(m_streamOutput.is_open())
	{
		this->OutHeader();
		this->ParseInputFile();
	}
	else
	{
		throw std::runtime_error("Can't open the output file");
	}

}//Setup

//Parsing and checking the validity of the input file
void CDynamic::ParseInputFile()
{
	// Boolean to check if the mandatory input variables have been provided
	bool					bBox(false),bAtoms(false);
	bool					bRandom(false); // Check if the positions have to be randomly set
	// These are the possible input variables
	double 					dA(0.0), dB(0.0), dC(0.0), dAlpha(90.0), dBeta(90.0), dGamma(90.0),dDist(0.0);
	
	std::ifstream input(m_strInputFile);
	if(input.is_open())
	{
		// The input file has been successfully opened
		
		m_streamOutput << "Reading the input file : " << m_strInputFile << std::endl << std::endl;

		for(std::string strLine; getline(input,strLine);)
		{
			m_streamOutput << strLine << std::endl;
			// This removes the comments (starting with "!")
			strLine = strLine.substr(0,strLine.find("!"));

			std::istringstream 	streamLine(strLine);
			std::string 		strFirstParam;

			streamLine >> strFirstParam;
			//std::getline(streamLine,strFirstParam,' ');
			
			std::string strParam;
			if(strcmp(strFirstParam.data(),"Box")==0)
			{
				for(std::string strBoxParam; strcmp(strBoxParam.data(),"End")!=0; getline(input, strBoxParam),m_streamOutput<<strBoxParam<<std::endl)
				{
					strBoxParam = strBoxParam.substr(0,strBoxParam.find("!"));
					std::istringstream streamBoxParam(strBoxParam);

					streamBoxParam >> strParam;

					if(strcmp(strParam.data(),"a")==0)
					{
						streamBoxParam >> dA;
						if(streamBoxParam.fail())
						{
							std::stringstream errMsg;
							errMsg << "Bad input with line : " << strBoxParam << std::endl;
							throw std::runtime_error(errMsg.str());
						}
						// We check if there was a unit specification
						streamBoxParam >> strParam;
						if(strcmp(strParam.data(),"A")==0)
						{
							dA *= ANGSTROM_TO_BOHR;
						}
					}
					else if(strcmp(strParam.data(),"b")==0)
					{
						streamBoxParam >> dB;
						if(streamBoxParam.fail())
						{
							std::stringstream errMsg;
							errMsg << "Bad input with line : " << strBoxParam << std::endl;
							throw std::runtime_error(errMsg.str());
						}
						// We check if there was a unit specification
						streamBoxParam >> strParam;
						if(strcmp(strParam.data(),"A")==0)
						{
							dB *= ANGSTROM_TO_BOHR;
						}
					}
					else if(strcmp(strParam.data(),"c")==0)
					{
						streamBoxParam >> dC;
						if(streamBoxParam.fail())
						{
							std::stringstream errMsg;
							errMsg << "Bad input with line : " << strBoxParam << std::endl;
							throw std::runtime_error(errMsg.str());
						}
						// We check if there was a unit specification
						streamBoxParam >> strParam;
						if(strcmp(strParam.data(),"A")==0)
						{
							dC *= ANGSTROM_TO_BOHR;
						}
					}
					else if(strcmp(strParam.data(),"alpha")==0)
					{
						streamBoxParam >> dAlpha;
						if(streamBoxParam.fail())
						{
							std::stringstream errMsg;
							errMsg << "Bad input with line : " << strBoxParam << std::endl;
							throw std::runtime_error(errMsg.str());
						}
					}
					else if(strcmp(strParam.data(),"beta")==0)
					{
						streamBoxParam >> dBeta;
						if(streamBoxParam.fail())
						{
							std::stringstream errMsg;
							errMsg << "Bad input with line : " << strBoxParam << std::endl;
							throw std::runtime_error(errMsg.str());
						}
					}
					else if(strcmp(strParam.data(),"gamma")==0)
					{
						streamBoxParam >> dGamma;
						if(streamBoxParam.fail())
						{
							std::stringstream errMsg;
							errMsg << "Bad input with line : " << strBoxParam << std::endl;
							throw std::runtime_error(errMsg.str());
						}
					}
					else if(strcmp(strBoxParam.data(),"")!=0)
					{
						std::stringstream errMsg;
						errMsg << "Bad input with line : " << strBoxParam << std::endl;
						throw std::runtime_error(errMsg.str());
					}
				}

				if(dAlpha+dBeta+dGamma > 360 or dAlpha+dBeta+dGamma < 0)
				{
					throw std::runtime_error("Bad angle value for the box\nCheck https://doi.org/10.1107/S0108767310044296 for more informations");
				}
				if(dAlpha-dBeta+dGamma > 360 or dAlpha-dBeta+dGamma < 0)
				{
					throw std::runtime_error("Bad angle value for the box\nCheck https://doi.org/10.1107/S0108767310044296 for more informations");
				}
				if(dAlpha+dBeta-dGamma > 360 or dAlpha+dBeta-dGamma < 0)
				{
					throw std::runtime_error("Bad angle value for the box\nCheck https://doi.org/10.1107/S0108767310044296 for more informations");
				}
				if(-dAlpha+dBeta+dGamma > 360 or -dAlpha+dBeta+dGamma < 0)
				{
					throw std::runtime_error("Bad angle value for the box\nCheck https://doi.org/10.1107/S0108767310044296 for more informations");
				}
				m_Box.SetA(dA);
				m_Box.SetB(dB);
				m_Box.SetC(dC);
				m_Box.SetAlpha(dAlpha);
				m_Box.SetBeta(dBeta);
				m_Box.SetGamma(dGamma);
				bBox = true;
			}
			else if(strcmp(strFirstParam.data(),"Atoms")==0)
			{
				streamLine >> strParam;
				if(strcmp(strParam.data(),"random")==0)
				{
					// random keyword specified, reading the number and atomic properties
					bRandom = true;

					streamLine >> dDist; 
					if(streamLine.fail())
					{
						std::stringstream errMsg;
						errMsg << "Bad input with line : " << strLine << std::endl;
						throw std::runtime_error(errMsg.str());
					}
					streamLine >> strParam;
					if(strcmp(strParam.data(),"A")==0)
					{
						// We convert from Angstrom to Bohr
						dDist *= ANGSTROM_TO_BOHR;
					}
					for(std::string strAtomParam; strcmp(strAtomParam.data(),"End"); getline(input,strAtomParam),m_streamOutput<<strAtomParam<<std::endl)
					{
						strAtomParam.substr(0,strAtomParam.find("!"));
						std::istringstream streamAtomParam(strAtomParam);

						unsigned int iN(0);
						double dEpsilon(0.0), dSigma(0.0), dMass(0.0);
						if(strcmp(strAtomParam.data(),"")!=0)
						{
							streamAtomParam >> iN >> dEpsilon; 
							std::stringstream::pos_type pos = streamAtomParam.tellg();
							streamAtomParam >> strParam;
							if(strcmp(strParam.data(),"K")==0)
							{
								// We convert from Kelvin to eV
								dEpsilon *= K_TO_EV;
								streamAtomParam >> dSigma >> strParam;
							}
							else if(strcmp(strParam.data(),"H")==0)
							{
								// We convert from Hartree to eV
								dEpsilon *= HARTREE_TO_EV;
								streamAtomParam >> dSigma;
							        pos = streamAtomParam.tellg();
								streamAtomParam	>> strParam;
							}
							else if(strcmp(strParam.data(),"J")==0)
							{
								// We convert from Joules to eV
								dEpsilon *= J_TO_EV;
								streamAtomParam >> dSigma;
							        pos = streamAtomParam.tellg();
								streamAtomParam	>> strParam;
							}
							else if(strcmp(strParam.data(),"cm")==0)
							{
								// We convert from cm-1 to eV
								dEpsilon *= CM_TO_EV;
								streamAtomParam >> dSigma;
							        pos = streamAtomParam.tellg();
								streamAtomParam	>> strParam;
							}
							else if(strcmp(strParam.data(),"kcalmol")==0)
							{
								// We convert from kcal/mol to eV
								dEpsilon *= KCALMOL_TO_EV;
								streamAtomParam >> dSigma;
							        pos = streamAtomParam.tellg();
								streamAtomParam	>> strParam;
							}
							else if(strcmp(strParam.data(),"kjmol")==0)
							{
								// We convert from kJ/mol to eV
								dEpsilon *= KJMOL_TO_EV;
								streamAtomParam >> dSigma;
							        pos = streamAtomParam.tellg();
								streamAtomParam	>> strParam;
							}
							else
							{
								streamAtomParam.clear();
								streamAtomParam.seekg(pos,streamAtomParam.beg);
								streamAtomParam >> dSigma;
							        pos = streamAtomParam.tellg();
								streamAtomParam	>> strParam;
							}
							if(streamAtomParam.fail())
							{
								std::stringstream errMsg;
								errMsg << "Bad input with line : " << streamAtomParam.str() << std::endl; 
								throw std::runtime_error(errMsg.str());
							}


							if(strcmp(strParam.data(),"A")==0)
							{
								// We convert from Angstrom to bohr
								dSigma *= ANGSTROM_TO_BOHR;
								streamAtomParam >> dMass;
							}
							else
							{
								streamAtomParam.clear();
								streamAtomParam.seekg(pos, streamAtomParam.beg);
								streamAtomParam >> dMass;
							}
							m_Box.AddAtoms(iN,dSigma,dEpsilon,dMass);
							bAtoms = true;
						}
					}
				}
			}
			else if(strcmp(strFirstParam.data(),"MD")==0)
			{
				unsigned int	iSteps(0);
				double		dTimeStep(0.0),dTemperature(0.0),dNeighbor(0.0),dCutoff(0.0);
				for(std::string strMDParam; strcmp(strMDParam.data(),"End")!=0;getline(input,strMDParam),m_streamOutput<<strMDParam<<std::endl)
				{
					strMDParam = strMDParam.substr(0,strMDParam.find("!"));
					std::istringstream streamMDParam(strMDParam);

					streamMDParam >> strParam;

					if(strcmp(strParam.data(),"Timestep")==0)
					{
						streamMDParam >> dTimeStep;
						if(streamMDParam.fail())
						{
							std::stringstream errMsg;
							errMsg << "Bad input with line : " << strMDParam << std::endl;
							throw std::runtime_error(errMsg.str());
						}
						streamMDParam >> strParam;
						// Need to make a change to the possible units 
					}
					else if(strcmp(strParam.data(),"Steps")==0)
					{
						streamMDParam >> iSteps;
						if(streamMDParam.fail())
						{
							std::stringstream errMsg;
							errMsg << "Bad input with line : " << strMDParam << std::endl;
							throw std::runtime_error(errMsg.str());
						}
					}
					else if(strcmp(strParam.data(),"Temperature")==0)
					{
						streamMDParam >> dTemperature;
						if(streamMDParam.fail())
						{
							std::stringstream errMsg;
							errMsg << "Bad input with line : " << strMDParam << std::endl;
							throw std::runtime_error(errMsg.str());
						}
					}
					else if(strcmp(strParam.data(),"Neighbor")==0)
					{
						streamMDParam >> dNeighbor;
						if(streamMDParam.fail())
						{
							std::stringstream errMsg;
							errMsg << "Bad input with line : " << strMDParam << std::endl;
							throw std::runtime_error(errMsg.str());
						}
						streamMDParam >> strParam;
						if(strcmp(strParam.data(),"A")==0)
						{
							dNeighbor *= ANGSTROM_TO_BOHR;
						}
					}
					else if(strcmp(strParam.data(),"Cutoff")==0)
					{
						streamMDParam >> dCutoff;
						if(streamMDParam.fail())
						{
							std::stringstream errMsg;
							errMsg << "Bad input with line : " << strMDParam << std::endl;
							throw std::runtime_error(errMsg.str());
						}
						streamMDParam >> strParam;
						if(strcmp(strParam.data(),"A")==0)
						{
							dCutoff *= ANGSTROM_TO_BOHR;
						}
					}

				}

				m_dTimeStep = dTimeStep;
				m_dInitTemperature = dTemperature;
				m_dNeighbor = dNeighbor;
				m_dCutoff = dCutoff;
				m_iNStep = iSteps;

			}
			else if(strcmp(strFirstParam.data(),"")!=0)
			{
				std::stringstream errMsg;
				errMsg << "Bad input with line : " << strLine ;
				throw std::runtime_error(errMsg.str());
			}
		}
		input.close();
		if(!bBox)
		{
			throw std::runtime_error("No box specified");
		}
		if(!bAtoms)
		{
			throw std::runtime_error("No atoms specified");
		}
		
		m_Box.Setup();
		if(bRandom)
		{
			m_Box.InitPosFromRandomDistribution(dDist);
			m_Box.InitSpeedRandom(m_dInitTemperature);
			std::ofstream	streamInitPos("Initial.pos");
			if(streamInitPos.is_open())
			{
				m_Box.OutAtomPos(streamInitPos);				
			}
			else
			{
				throw std::runtime_error("Can't open the file Initial.pos");
			}
		}
		m_Box.NeighborList(m_dCutoff,m_dNeighbor);

		
		m_streamOutput << std::endl << "Input file is consistent" << std::endl;
		m_streamOutput << std::string(150,'=') << std::endl;
		m_streamOutput << boost::format("%=150s")%"PRINTING THE PARAMETERS AT THE BEGINNING OF THE COMPUTATION" << std::endl;
		m_Box.OutBoxParam(m_streamOutput);
		m_streamOutput << std::string(150,'-') << std::endl;
		m_streamOutput << boost::format("%-25s %i")%"Number of threads"%m_iNumberThreads << std::endl;
		m_streamOutput << boost::format("%-25s %i")%"Number of steps"%m_iNStep << std::endl;
		m_streamOutput << boost::format("%-25s %d picoseconds")%"Timestep"%m_dTimeStep << std::endl;
		m_streamOutput << boost::format("%-25s %d K")%"Temperature"%m_dInitTemperature << std::endl;
		m_streamOutput << std::string(150,'=') << std::endl;
		m_Box.ComputeForces();
	}
	else
	{
		// The input file can not be opened
		throw std::runtime_error("The input file can not be opened");
	}
} //m_ParseInputFile

// This applies the Berendsen thermostat to the system
void	CDynamic::Berendsen()
{
	double dTargetT(0.0),dEffectiveT(m_Box.GetTemperature());

	dTargetT = dEffectiveT+(m_dInitTemperature-dEffectiveT)*m_dTimeStep/m_dBerendsenTau;

	for(unsigned int i=0;i<m_Box.GetNAtom();i++)
	{
		CAtom& a = m_Box.GetAtom(i);
		CSpeed s = a.GetSpeed();
		double dScaling = (s.Norm()*BOHR_TO_ANGSTROM*ANGSTROM_TO_M)/sqrt(3*KB*dTargetT/(a.GetMass()*DAL_TO_KG));
		s /= dScaling;	   
	    	a.SetSpeed(s);	
	}
}//Berendsen


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

