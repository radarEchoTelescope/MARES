/*/////////////////////////////////////////////////////////////////////////////
This file defines the cascade class and all its associated functions
*//////////////////////////////////////////////////////////////////////////////
#ifndef CASCADE_hh
#define CASCADE_hh

#include "macro_settings.hh"

// -----------------------------------------------------------------------------
class Cascade {
public:

  Cascade(int evt, double cenergy, double xpos, double ypos, double zpos,
          double czenith, double cazimuth, double nenergy,
          double nzenith, double nazimuth, double oneweight);

  // Cascade's methods

  /* Particle (electron) density for penetration length X and radius r.
  [g/cm^2, cm, GeV] */
  double Density(double X, double r, double delta_r = 0.01);      // [#e-/ cm^3]
  double PlasmaFreq(const double &dens);         // [Hz]
  double Absorption(const double &dens, const double &freq_obs);
  double SkinDepth(const double &dens, const double &freq_obs);

  // Accesors
  int  Evt() const;
  double  Energy() const;
  std::vector<double> Pos() const;
  std::vector<double> Dir() const;
  std::vector<double> Sph() const;

  double* Parent();

  double Ldiv() const;
  double Xdiv() const;
  double Rdiv() const;

  double Ltot() const;
  double Xtot() const;
  double Rtot() const;

  double Lbins() const;
  double Xbins() const;
  double Rbins() const;

  std::vector<std::vector<double>> Density();


private:
  friend class Scatter; // Scatter can access private members.

  int    fEvent;            // Event number.
  double fEnergy;          // Energy of the cascade.

  // Interaction point's position (Shower start, head).
  std::vector<double> fPosition{0,0,0};

  // The cascade development direction is now set from the sph_ang.
  std::vector<double> fDirection{0,0,0}; // Cartesian vector direction

  // Incoming spherical angles (spherical coordinates) w.r.t the Earth frame.
  std::vector<double> fSphericalAngles{0,0};    // {Theta = zenith, phi= azimuth}

  // Neutrino parent values: energy, zenith, azimuth, oneweight.
  double fNeutrino[4];


  double fLdiv = 1;       // [cm] Length interval, resolution.
  double fXdiv = fLdiv * rho_ice;
  double fRdiv = 1;       // [cm] radial interval, resolution.

  // Max cascade depth: X_tot = 4* max depth from Heitler model estimate.
  double fLtot;           // [cm]
  double fXtot;           // [g/cm^2] Penetration depth
  double fRtot;           // [cm]

  double fLbins;          // Number of bins for cascade length
  double fXbins;          // Number of bins for cascade depth
  double fRbins;          // Number of bins for cascade radius



  //

  // std::vector<std::vector<double>> fDensity;        // [#e-/ cm^3]
  // double fRwaist;
  // std::vector<double> fRcrit; // [cm]
  // std::vector<std::vector<double>> fPlasmaFrequency;
  // std::vector<std::vector<double>> fAbsorption;
  // std::vector<std::vector<double>> fSkinDepth;
  // std::vector<std::vector<double>> fReflectance;
  // std::vector<std::vector<double>> fOpacity;
  //
  // void SetDensity();
  // void Rcrit(const std::vector<std::vector<double>> &density, const double & freq);
  // void PlasmaFreq(const std::vector<std::vector<double>> &density);
  // void Absorption(const std::vector<std::vector<double>> &density, const double & freq);
  // void SkinDepth(const std::vector<std::vector<double>> &density, const double & freq);
  // void Reflectance(const std::vector<std::vector<double>> &absorption);
  // void Opacity(const std::vector<std::vector<double>> &absorption);
  // void RCS(const std::vector<std::vector<double>> &reflectance);




  std::vector<std::vector<double>> fDensity;
  void SetDensity();

  // double fRwaist;
  // std::vector<double> fRcrit; // [cm]
  // void SetRcrit(const std::vector<std::vector<double>> &density = fDensity,
  //               const double &freq = freq_obs);

};

std::vector<Cascade> load_cascade_file(const std::string& cs_filepath);

// Helper functions, used by other functions only
namespace{

  /* Ne, number of particles in the cascade */
  double Ne(double X, double E);

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
