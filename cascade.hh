/*/////////////////////////////////////////////////////////////////////////////
This file defines the cascade class and all its associated functions
*//////////////////////////////////////////////////////////////////////////////
#ifndef CASCADE_hh
#define CASCADE_hh

#include "macro_settings.hh"

// ----- Computational constants -----------------------------------------------
const int nbin = 500;                 // # of bins for integration/array filling

// -----------------------------------------------------------------------------
class Cascade {
public:

  Cascade(int evt, double cenergy, double xpos, double ypos, double zpos,
          double czenith, double cazimuth, double nenergy,
          double nzenith, double nazimuth, double oneweight);

  /* Particle (electron) density for penetration length X and radius r. */
  double dens(double X, double r);      // [g/cm^2, cm, GeV]

  std::vector<double> get_rcrit();
  std::vector<std::vector<double>> get_density();
  std::vector<std::vector<double>> get_reflectivty_2D();
  std::vector<std::vector<double>> get_reflectance_2D();

  // Accesors
  double  event();
  double  energy();
  double* sph_angles();
  double* position();
  double* direction();
  double* parent();

  double Xtot();
  double get_X_bin();
  double get_r_tot();
  double get_r_bin();
  double get_L_bin();

private:
  friend class Scatter; // Antena can access private members.

  int    evtnr;            // Event number.
  double E_p;          // Energy of the cascade.
  double sph_ang[2];      // {Theta = zenith, phi= azimuth}
  // Incoming spherical angles (spherical coordinates) w.r.t the detector frame.
  double pos[3];          // Interaction point's position (Shower start, head).
  // The cascade development direction is set from the sph_ang.
  double dir[3];          // Cartesian vector direction
  // Neutrino parent values: energy, zenith, azimuth, oneweight.
  double neutrino[4];

  // Max cascade depth: X_tot = 4* max depth from Heitler model estimate.
  double X_tot;           // [g/cm^2] Penetration depth
  double X_bin;
  double r_tot;           // [cm]
  double r_bin;
  double L_tot;           // [m]
  double L_bin;

  std::vector<std::vector<double>> density;         // [#e-/ cm^3] The variable
  void set_density();   // The function that generates the variable.

  double r_waist;         // [cm]
  std::vector<double> r_crit;
  void set_rcrit(const double & freq_obs = 1E9);    // Default frequency.

  /* Compute the 2D reflectance and reflectivity for varying skin depth.

  Physics note:
  Reflectivity = Property of a material.
  Reflectance = Refers to a specific sample, depends on size and other params.

  Reflectivity is the reflectance value as the object becomes thick.
  Thus, reflectivity is defined here as the integral of reflectance over
  the layers.
  */
  std::vector<std::vector<double>> reflectance2D, reflectivity2D;
  void set_reflectivity_2D();

  // std::vector<double> od_cs_1D;
  // void set_od_cs_1D();

};

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
  void set_max_depth(E_p){X_tot = 4*(log(E_p/E_c)/log(2) )*X_int;}
*/

std::vector<Cascade> load_cascade_file(const std::string& cs_filepath);

// Helper functions, used by other functions only
namespace{

  /* Ne, number of particles in the cascade */
  double Ne(double X, double E_p);

  /* Shower age */
  double sh_age(double X, double E_p);

  /* Lateral particle distribution for radius r and penetration length X. */
  double wiv1(double r, double s);

  /* Integral of lateral particle distribution between two radii. */
  double intwiv(double r, double delta_r, double s = 1.01);    // [cm, cm, g/cm^2, GeV]
  // "Hard-coded" s = 1 for now. s = 1.01 to avoid divergencies.
}

#endif
// -----------------------------------------------------------------------------
