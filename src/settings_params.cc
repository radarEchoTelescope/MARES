#ifndef SETTINGS_PARAMS_cc
#define SETTINGS_PARAMS_cc

#include "settings.hh"

/* These are the mathematical model's variables.
    By making these parameters non-const in the code,
 we allow for them to be changed at the beginning of a run.
 (For example, in a macro file). 

    The numbers here act as a "default", when no specific value is needed. 
    They include CLHEP-like units in their definition. 
    Natural units: mm, ns, MeV.
*/  

/* Default computational parameters */
// These work as the default values for the constructors in MARES.



// Spatial resolution, effective scatterer size. 
double _dL = 1 *cm; // [cm/bin]
double _dR = 1 *mm;
// Unused in a 1D cascade
double _dN = 1 *mm; // This is also a radial direction.

// Temporal resolution, sampling ratio.
// We found that best results are with a sampling between 10x and 100x freq_obs.
double _sampling = 50;


      // Plasma ----------------------------------------------------------------
double _lifetime = 3 *ns;            // Mean plasma lifetime
double _f_coll= 64.733 *THz;         // [Hz] RS collision frequency
double _memp = 1;                    // Plasma to electron mass ratio

      // Ice  ------------------------------------------------------------------
// double _refindex=1.78;               // refractive index
double _refindex=1.51;               // refractive index
double _att_length=1450 *m;          // Attenuation length
double _rho_ice = 0.92 *g/pow(cm,3); // Density, from GEANT


double _r_moliere = 7 *cm;           // Moliere Radius in ice
// double _E_ionization = 20 *eV;    // e- ionization energy 
double _E_ionization = 69 *eV;       // e- ionization energy [RS]

double _E_deposition = 2 *MeV/(g/pow(cm,2));       
double _E_c = 78.6 * MeV;          // Critical cascade energy for ionization
double _X_0 = 36.08 *g/pow(cm,2);    // Radiation columm density

#endif