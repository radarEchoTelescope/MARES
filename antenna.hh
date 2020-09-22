/*//////////////////////////////////////////////////////////////////////////////
This header file defines the antenna positions.
*///////////////////////////////////////////////////////////////////////////////
#ifndef ANTENNA_hh
#define ANTENNA_hh

#include "macro_settings.hh"

// ----- Default frequency -----------------------------------------------------
// const double freq_obs = 1E9;
const double freq_obs = 450 * 1E6;

// -----------------------------------------------------------------------------
class Antenna{

public:

  Antenna(double power, double xpos, double ypos, double zpos,
          double theta, double phi, double gain, double frequency);
  Antenna(double power, double xpos, double ypos, double zpos);

  double  power();
  double  gain();
  double* position();
  double* polarization();
  double  distance();
  double  projection();
  double  angle();

  double frequency();
  double lambda();
  double omega();
  double wavenr();


private:
  friend class Scatter;
  double power_;
  double gain_;

  double pos[3];
  double polar[2];    // Polarization

  double dir[3] = {0};         // Vector direction to cs.
  double dist = 0;             // Module of distance to cs.
  double dot = 0;              // Dot (inner) product with cascade direction.
                               // The cosine of the angle between them.
  double ang = 0;            // The angle.

  double IRT_dist[2];
  double IRT_angle[2];

  double f_obs;               // [Hz] Observer frequency = 1E9 , 450 *1E6
  double k_obs;               // Wavenumber.
  double l_obs;               // [m] Detection wavelength
  double w_obs;               // Angular freq_obs

};


class Detector{

public:

// Constructors
  Detector();               // Default, empty?
  Detector(std::string name);   // Model-based constructor
  /* Other models:
    - Bistatic
    - RET_CR
  */

  Detector(int a, int b); // Krijn's old system locations

  std::vector<Antenna> get_transmitters();
  std::vector<Antenna> get_receivers();

  void add_antenna(Antenna& at);

private:

  std::vector<Antenna> Transmitters;
  std::vector<Antenna> Receivers;

};

#endif
