#ifndef SETTINGS_DEFAULTS_cc
#define SETTINGS_DEFAULTS_cc

#include "settings.hh"

/* MARES default parameters.
Default parameters are indicated by a "_" (_parameter).
They include CLHEP-like units in their definition. 
These values can be updated externally (CAREFUL!),
BUT they need to be updated in the correct units described here.
[unit]
*/  

/* --------- Plasma properties  --------------- */

// Mean plasma lifetime [ns] **** can have us, ms, s also
double _lifetime = 10 *ns;      
//double _lifetime = 50 *ns;         

// Collision frequency estimate [THz]
// double _f_coll = 100 *THz;
// RS value for collision frequency
double _f_coll = 64.733 * THz;         

// Effective plasma mass ratio w.r.t electron
// i.e, effective plasma mass in electron masses
double _memp = 1;

/* --------- Medium properties  --------------- */

// Refractive index (Ice)
// double _refindex = 1.78;      //1.78         
double _refindex = 1.78;      //1.78      
// Refractive index (HDPE)
// double _refindex = 1.51;

// Density (Ice) [g/cm^3]
double _rho_ice = 0.92 * g/pow(cm,3); 
// Density (sea-level air)
//double _rho = 1.168e-3;

// Attenuation length [m]
double _att_length = 1450.0 *m;
        
/* --------- Cascade properties  --------------- */

// Moliere Radius (in ice) [cm]
double _r_moliere = 7.0 *cm;

// Electron ionization energy from secondary cascades 
double _E_ionization = 20.0 *eV;
// Electron ionization energy (RadioScatter)
// double _E_ionization = 69 * eV;       

// Electron deposition energy (dE/dX) [MeV/(g/cm^2)]
// Mass stopping power of ice
// i.e. energy loss per ionizing particle (@ 1 GeV)
double _E_deposition = 2.0 *MeV/(g/pow(cm,2));       

// Critical (threshold) energy for electron ionization [MeV]
double _E_c = 78.6 *MeV;

// Radiation "length" (columm density) in ice [g/pow(cm,2)]
double _X_0 = 36.08 *g/pow(cm,2);    
// Radiation length in air
// double _X_0 = 36.7 *g/pow(cm,2);    

/* --------- Simulation properties  --------------- */

// Spatial resolution, effective scatterer size. 
double _dL = 1.0 *cm; // [cm/bin]                                 // default = 1.0*cm   modified test, 0.5cm, 0.5mm
double _dR = 1.0 *mm;                                             // default = 1.0*mm


// default dL dR (resolution and effective scatter)
// double _dL = 1 * cm;
// double _dR = 1 * mm;

// Temporal resolution, sampling ratio.
// We found that best results are computed with 
// a sampling ratio between 10x and 100x freq_obs.
double _sampling = 100;                                            // default = 50

double _c_ice = c_vac/_refindex;  // [m/s] speed of light in ice
double _Z_ice = Z_0/_refindex;            // Impedance of ice
double _L_0 = _X_0/_rho_ice;       // Radiation length = 39.22 cm

// Ltot_factor can be 2 safely for speed, but 3 makes prettier plots. 
double _Ltot_factor = 3.0;
double _Rtot_factor = 2.0;

/*
The cascade values are computed in a frame which dimensions
are determined in terms of the the typical length scale in 
every dimension. Total size = typical* factor
*/

// // Plasma
// double _lifetime = 10 *ns;            // Mean plasma lifetime
// double _f_coll= 100 *THz;         // [Hz] Collision frequency estimate
// //double _f_coll= 64.733 *THz;         // [Hz] RS collision frequency
// double _memp = 1;                    // Plasma to electron mass ratio

//       // Ice  ------------------------------------------------------------------
// double _refindex=1.78;               // refractive index
// // double _refindex=1.51;               // refractive index
// double _att_length=1450 *m;          // Attenuation length
// double _rho_ice = 0.92 *g/pow(cm,3); // Density, from GEANT


// double _r_moliere = 7 *cm;           // Moliere Radius in ice
// // double _E_ionization = 20 *eV;    // e- ionization energy 
// double _E_ionization = 69 *eV;       // e- ionization energy [RS]

// double _E_deposition = 2 *MeV/(g/pow(cm,2));       
// double _E_c = 78.6 * MeV;          // Critical cascade energy for ionization
// double _X_0 = 36.08 *g/pow(cm,2);    // Radiation columm density

//       /*The typical length of a cascade is log(12.72 * fEnergy) * X_0;	
//  	The typical radius of a cascade is the moliere radius;
// */




#endif
