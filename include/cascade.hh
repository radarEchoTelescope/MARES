/*/////////////////////////////////////////////////////////////////////////////
This file defines the cascade class and all its associated functions
*//////////////////////////////////////////////////////////////////////////////
#ifndef CASCADE_hh
#define CASCADE_hh

#include "settings.hh"

// -----------------------------------------------------------------------------
class Cascade {
public:

  Cascade(int evt, double cenergy, double xpos, double ypos, double zpos,
          double czenith, double cazimuth, int np = 1);

  Cascade(int evt, double cenergy, double xpos, double ypos, double zpos,
          double czenith, double cazimuth, int np, 
          double nenergy, double nzenith, double nazimuth, double oneweight);

// Basic accessors

  int    Evt()    const;
  double Energy() const;
  double Ltot()   const;
  double Rtot()   const;
  double Xtot()   const;
  std::vector<double> Dir() const;
  std::vector<double> Pos() const;
  std::vector<double> Sph() const;
  double* Parent();

// Method's storage accesors
// User should only read and store the values after they are made.

  std::vector<std::vector<double>> Density();
  std::vector<std::vector<double>> PlasmaFreq();
  std::vector<std::vector<double>> Absorption();
  std::vector<std::vector<double>> SkinDepth();

  std::vector<std::vector<double>> Reflectance();
  std::vector<std::vector<double>> Opacity();
  std::vector<std::vector<double>> Transparency();


protected:
    friend class Scatter; // Scatter can access private members.
    // Scatter has access to setter functions. 

  int    fEvent;            // Event number.
  int    fNp;               // Number of primaries; 
  double fEnergy;          // Energy of the cascade.
  
  // Max cascade depth: X_tot = 4* max depth from Heitler model estimate.
  double fLtot;           // [cm]
  double fRtot;           // [cm]
  double fXtot;           // [g/cm^2] Penetration depth

  // Interaction point's position (Shower start, head).
  std::vector<double> fPosition{0,0,0};

  // The cascade development direction is now set from the sph_ang.
  std::vector<double> fDirection{0,0,0}; // Cartesian vector direction

  // Incoming spherical angles (spherical coordinates) w.r.t the Earth frame.
  std::vector<double> fSphericalAngles{0,0};    // {Theta = zenith, phi= azimuth}

  // Neutrino parent values: energy, zenith, azimuth, oneweight.
  double fNeutrino[4];

// Functions and methods

/* (Classical) Particle (electron) number for penetration length X and radius r.
  [g/cm^2, cm, GeV] */
  // Uses 1 particle with fEnergy.
  double Ne(double X, double r, double delta_r);           // [#e-/cm]

  /* Particle (electron) density for penetration length X and radius r.
  [g/cm^2, cm, GeV] */
  // In that case each particle carries E energy, the mean energy per primary.
  double Ne(double X, double r, double delta_r, double Ep, double Np);           // [#e-/cm]

  // (Classical) density, in case the cascade was orginated from 1 particle.
  // Uses 1 particle with fEnergy.
  double Density(double X, double r, double delta_r);      // [#e-/ cm^3]

 /* Particle (electron) density for penetration length X and radius r.
  [g/cm^2, cm, GeV] */
  // Density in case the cascade was orginated from Np particles
  // In that case each particle carries E energy, the mean energy per primary.
  double Density(double X, double r, double delta_r, double Ep, double Np); 

  // TO-DO: 2 species absorption! 

  double PlasmaFreq(const double &dens);         // [Hz]
  // Electron absorption

  // Full electron absorption formula, defaults as free (unbound) 
  // SINGLE SYSTEM OF UNITS! THESE TWO FUNCTIONS NOW RETURN VALUES IN mm! 
  double Absorption(const double &dens, const double &freq_obs, const double &freq_nat = 0);
  double SkinDepth(const double &dens, const double &freq_obs, const double &freq_nat = 0);

  // Cascade's matrix storage and methods

  std::vector<std::vector<double>> fDensity;
  std::vector<std::vector<double>> fPlasmaFrequency;
  std::vector<std::vector<double>> fAbsorption;
  std::vector<std::vector<double>> fSkinDepth;
  std::vector<std::vector<double>> fReflectance;
  std::vector<std::vector<double>> fOpacity;
  std::vector<std::vector<double>> fTransparency;

// TO_DO: Parallelise all the loops from here down with OMP.

// length_vals and radius_vals are 2D matrices with the coordinates to evaluate the density.
  void Density( const std::vector<std::vector<double>> &length_vals,
                const std::vector<std::vector<double>> &radius_vals,
                const double delta_r  );

  void PlasmaFreq(const std::vector<std::vector<double>> &density);
  void Absorption(const std::vector<std::vector<double>> &density, const double & freq);
  void SkinDepth(const std::vector<std::vector<double>> &density, const double & freq);
  // Transparency also computes reflectivity and opacity
  void Transparency(const std::vector<std::vector<double>> &density, const double & freq, const double & delta_r);

  // double fRwaist;
  // std::vector<double> fRcrit; // [cm]
  // void SetRcrit(const std::vector<std::vector<double>> &density = fDensity,
  //               const double &freq = freq_obs);

};

std::vector<Cascade> load_cascade_file(const std::string& cs_filepath);

// Helper functions, used by other functions only
namespace{

  /* Ne, number of particles in the cascade */
  double N(double X, double E);

  /* Shower age */
  double ShowerAge(double X, double E);

  /* Lateral particle distribution for radius r and penetration length X. */
  double wiv1(double r, double s);

  /* Integral of lateral particle distribution between two radii. */
  double intwiv(double r, double delta_r, double s = 1.01);    // [cm, cm, g/cm^2, GeV]
  // "Hard-coded" s = 1 for now. s = 1.01 to avoid divergencies.
}

#endif
// -----------------------------------------------------------------------------


/* You don't want to allow the code to edit these values.
These are needed if you add the default empty constructor.

  // Cascade() {}

  void set_position(double x, double y, double z){ pos = {x,y,z}; }
  void set_direction_w_sph_ang(double theta_ang, double phi_ang){
    theta = theta_ang;
    phi = phi_ang;
    dir[0] = sin(theta)*cos(phi);
    dir[1] = sin(theta)*sin(phi);
    dir[2] = cos(theta);
  }
  void set_direction_cartesian(double x, double y, double z){ dir = {x,y,z}; }
  void set_parent_neutrino(double nenergy, double nazimuth, double nzenith,
    double oneweight) {neutrino = {nenergy, nzenith , nazimuth, oneweight}; }
  void set_max_depth(fEnergy){X_tot = 4*(log(fEnergy/E_c)/log(2) )*X_int;}
*/
