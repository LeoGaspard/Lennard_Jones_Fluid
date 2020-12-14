///////////////////////////////////////////////////////////////
//NAME			Constants.h
//PURPOSE		Define all the constants that will
//			be used by the program
//
//FUNCTIONS/OBJECTS	None
//
//AUTHOR:		Léo Gaspard
//////////////////////////////////////////////////////////////

#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED

// Constants

#define PI 			3.14159265358979324
#define C_LIGHT 		299792458 		// m.s^{-1}
#define H_PLANCK 		6.626070040E-34 	// J.s
#define AVOGADRO 		6.022140857E+23 	// mol^{-1}
#define E_CHARGE 		1.6021766208E-19 	// C
#define EPSILON_VOID 		8.854187818E-12 	// F.m^{-1}
#define BOHR_RADIUS 		5.2917721056E-11 	// m
#define HARTREE_ENERGY 		4.359744651E-18 	// J
#define ELECTRON_MASS		9.1093837015E-31	// kg
#define KB			1.380649E-23		// J.K^{-1}

// Length conversion

#define BOHR_TO_ANGSTROM	0.529177
#define ANGSTROM_TO_BOHR	1.889727
#define ANGSTROM_TO_M		1E-10
#define M_TO_ANGSTROM		1E+10

// Time conversion
#define	FS_TO_S			1E-15
#define S_TO_FS			1E+15

// Angle conversion
#define RADIAN_TO_DEGREE	57.295779513
#define DEGREE_TO_RADIAN	0.01745329252

// Mass conversion

#define DAL_TO_KG		1.66054E-27
#define KG_TO_DAL		6.02214E+26

// Energy conversions

#define HARTREE_TO_J  		4.35974E-18
#define J_TO_HARTREE  		2.29371E+17
#define HARTREE_TO_EV 		27.2114
#define EV_TO_HARTREE 		3.67493E-2
#define HARTREE_TO_CM 		2.19475E+5
#define CM_TO_HARTREE 		4.55634E-6
#define HARTREE_TO_KJMOL 	2625.50
#define KJMOL_TO_HARTREE 	3.80880E-4
#define HARTREE_TO_KCALMOL 	627.509
#define KCALMOL_TO_HARTREE 	1.59360E-3
#define J_TO_EV 		6.24151E+18
#define EV_TO_J 		1.60218E-19
#define J_TO_CM 		5.03412E+22
#define CM_TO_J 		1.98645E-23
#define J_TO_KJMOL 		96.4853
#define KJMOL_TO_J 		1.66054E-21
#define J_TO_KCALMOL 		1.43933E+20
#define KCALMOL_TO_J 		6.94770E-21
#define EV_TO_CM 		8065.54
#define CM_TO_EV 		1.23984E-4
#define EV_TO_KJMOL 		96.4853
#define KJMOL_TO_EV 		1.03643E-2
#define EV_TO_KCALMOL 		23.0605
#define KCALMOL_TO_EV 		4.33641E-2
#define CM_TO_KJMOL 		1.19627E-2
#define KJMOL_TO_CM 		83.5935
#define CM_TO_KCALMOL 		2.85914E-3
#define KCALMOL_TO_CM 		349.755
#define KJMOL_TO_KCALMOL 	0.239006
#define KCALMOL_TO_KJMOL 	4.18400
#define K_TO_J			1.3806E-23
#define J_TO_K			7.2430E+22
#define HARTREE_TO_K		3.15777E+05
#define K_TO_HARTREE            3.16678E-06
#define K_TO_EV			8.61705E-05
#define EV_TO_K			1.16049E+04
#define K_TO_CM			6.95028E-01
#define CM_TO_K			1.42879
#define K_TO_KCALMOL		1.98717E-03
#define KCALMOL_TO_K		5.03228E+02
#define K_TO_KJMOL		8.31435E-03
#define KJMOL_TO_K		1.20274E+02

#endif //CONSTANTS_H_INCLUDED
