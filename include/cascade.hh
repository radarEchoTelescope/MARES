/*/////////////////////////////////////////////////////////////////////////////
This file defines the cascade class and all its associated functions
*//////////////////////////////////////////////////////////////////////////////
#ifndef CASCADE_hh
#define CASCADE_hh

#include "settings.hh"

// -----------------------------------------------------------------------------
class Cascade {
public:

  Cascade(double xpos, double ypos, double zpos,
          double zenith, double azimuth, 
          double energy, int primaries = 1, 
          int n_num = 0, int t_num = 0);

  Cascade(double xpos, double ypos, double zpos,
          double zenith, double azimuth, 
          double energy, int primaries,
          int n_num, int t_num, int p_id, int i_type, int channel,
          double inelasticity, double oneweight);

// Basic accessors

  int    Evt()    const;
  int    Nprimaries() const;
  double Energy() const;
  double Ltot()   const;
  double Rtot()   const;
  double Xtot()   const;
  std::vector<double> Dir() const;
  std::vector<double> Pos() const;
  std::vector<double> Sph() const;
  int* Info();

// Method's storage accesors
// User should and can only read and store the values after they are made.

  std::vector<std::vector<double>> Ne();
  std::vector<std::vector<double>> Density();
  std::vector<std::vector<double>> PlasmaFreq();

  // 
  std::vector<std::vector<double>> Absorption();
  std::vector<std::vector<double>> SkinDepth();

  std::vector<std::vector<double>> Reflectance();
  std::vector<std::vector<double>> Opacity();
  std::vector<std::vector<double>> Transparency();

  void save_output_files(const std::string& output_path, const std::array<bool, 8>& flags = {0});
  /*In order to save the member functions, you need to 
  pass an array of flags choosing what members to save:

    True = Save this field
    False = Don't save this field

  0 = Ne
  1 = Density
  2 = Plasma freq
  3 = Absorption
  4 = Skin depth
  5 = Reflectance
  6 = Opacity
  7 = Transparency
  */

protected:
    friend class Scatter; // Scatter can access private members.
    // Scatter has access to setter functions. 

  int    fEvent;            // Event number.
  int    fNp;               // Number of primaries; 
  double fEnergy;          // Energy of the cascade.
  double fBy;               // Inelasticity (Bjorken-y)
  double fOneweight;        // oneweight
  
  // Max cascade depth: X_tot = 4* max depth from Heitler model estimate.
  double fXtot;           // Penetration depth
  double fLtot;           // 
  double fRtot;           // 
  double fRnorm;          // Normalization factor for intwiv for this Rtot.

  // Interaction point's position (Shower start, head).
  std::vector<double> fPosition{0,0,0};

  // The cascade development direction is now set from the sph_ang.
  std::vector<double> fDirection{0,0,0}; // Cartesian vector direction

  // Incoming spherical angles (spherical coordinates) w.r.t the Earth frame.
  std::vector<double> fSphericalAngles{0,0};    // {Theta = zenith, phi= azimuth}

  // Neutrino parent values: energy, zenith, azimuth, oneweight.
  int fParent[5];

// Functions and methods

// TODO= Revise Units here

/* (Classical) Particle (electron) number for penetration length X.
  [g/mm^2] */
  // Uses 1 particle with fEnergy.
  double Ne(double X);           // [#e-/mm]

  /* Particle (electron) density for penetration length X.
  [g/mm^2, MeV] */
  // In that case each particle carries E energy, the mean energy per primary.
  double Ne(double X, double Ep, double Np);           // [#e-/cm]


  // Particle electron density evaluated at X for a variable density
  double Ne_vd(double X, double Ep, double Np, double rho_vd);   // [#e-/cm]

  // (Classical) density, in case the cascade was orginated from 1 particle.
  // Uses 1 particle with fEnergy.
  double Density(double X, double r, double delta_r);      // [#e-/ mm^3]

 /* Particle (electron) density for penetration length X and radius r.
  [g/cm^2, cm, GeV] */
  // Density in case the cascade was orginated from Np particles
  // In that case each particle carries E energy, the mean energy per primary.
  double Density(double X, double r, double delta_r, double Ep, double Np); 

  // TO-DO: 2 species absorption! 

  double PlasmaFreq(const double &dens);         // [Hz]
  // Electron absorption

  // Full electron absorption formula, defaults as free (unbound) 
  double Absorption(const double &dens, const double &freq_obs, const double &freq_nat = 0);
  double SkinDepth(const double &dens, const double &freq_obs, const double &freq_nat = 0);

  // Cascade's matrix storage and methods

  std::vector<std::vector<double>> fNe;
  std::vector<std::vector<double>> fDensity;
  std::vector<std::vector<double>> fPlasmaFrequency;
  std::vector<std::vector<double>> fAbsorption;
  std::vector<std::vector<double>> fSkinDepth;
  std::vector<std::vector<double>> fReflectance;
  std::vector<std::vector<double>> fOpacity;
  std::vector<std::vector<double>> fTransparency;

// TO_DO: Parallelise all the loops from here down with OMP.

// length_vals and radius_vals are 2D matrices with the coordinates to evaluate the density.
  void Ne( const std::vector<std::vector<double>> &length_vals,
                const std::vector<std::vector<double>> &radius_vals,
                const double delta_r  );

  void Density( const std::vector<std::vector<double>> &length_vals,
                const std::vector<std::vector<double>> &radius_vals,
                const double delta_r  );

  void PlasmaFreq(const std::vector<std::vector<double>> &density);
  void Absorption(const std::vector<std::vector<double>> &density, const double & freq);
  void SkinDepth(const std::vector<std::vector<double>> &density, const double & freq);
  // Transparency also computes reflectivity and opacity
  void Transparency(const std::vector<std::vector<double>> &density, const double & freq, const double & delta_r);

};

// Load a config or a config file and returns the cascades and 
// outputs the rejected ones into a ostream (file or terminal). 
std::vector<Cascade> load_cascade_file(const std::string& cs_config_filepath, std::ostream& out = std::cout);
std::vector<Cascade> load_cascade_config(libconfig::Config& cs_config, std::ostream& out = std::cout);

// Helper functions, used by other functions only
namespace NKG{

  /* Ne, number of particles in the cascade */
  double N(double X, double E);

  /* Shower age */
  double ShowerAge(double X, double E);

  /* Lateral particle distribution for radius r and penetration length X. */
  double wiv1(double r, double s);

  /* Integral of lateral particle distribution between r and r + dr.
		This integral is performed using the trapezoid rule.

		https://www.whitman.edu/mathematics/calculus_online/section08.06.html

		The error associated with this method is known. Over the the step size

		E(step) = dr^3/ (12* imx^2) *M where M is the maximum value of the second
		derivative of the function over the interval [r, r + dr].
 */
  // double intwiv(double r, double delta_r, double s = 1.01, double imx = 50.0);
  // "Hard-coded" s = 1 for in-ice showers. s = 1.01 to avoid divergencies.
  // Imx is the number of steps that are used for the integral. 
  // double intwiv(double r, double delta_r, double s = 1.01, double imx = 50.0);  //

  //Adjusted LDF function to have a user-set shower age.
  double intwiv(double r, double delta_r, double s, double imx = 50.0);

  // Function to calculate the ice density at a certain depth for a South Pole density model
  double rhoSouthPole(double depth); 

  // Function to calculate the ice density at a certain depth for a Greenland density model
  double rhoGreenland(double depth); 

  // Correction function for the NKG+ cosmic ray method.
  double removeEarlyN(double X, double a, double b);

}

#endif
// -----------------------------------------------------------------------------