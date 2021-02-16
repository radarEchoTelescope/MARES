/*//////////////////////////////////////////////////////////////////////////////
This header file defines the antenna positions.
*///////////////////////////////////////////////////////////////////////////////
#ifndef ANTENNA_hh
#define ANTENNA_hh

#include "macro_settings.hh"

class Antenna{

public:

  Antenna(double power, double xpos, double ypos, double zpos,
          double xpol, double ypol, double zpol, double gain, double frequency);
  Antenna(double power, double xpos, double ypos, double zpos);

  double  power();
  double  gain();
  // double* position();
  std::vector<double> position();
  std::vector<double> polarization();
  std::vector<double> direction();
  std::vector<double> sph_angles();
  double  distance();

  double  projection();
  double  angle();

  double freq();
  double lambda();
  double omega();
  double wavenr();


private:
  friend class Scatter;
  double _power;
  double _gain;
  double _leff = 1;
  double _efficiency = 1;
  double _load = 50; // pi * [Ohm]

  std::vector<double> _polar{0, 0, 1};    // Polarization
  // Standard antennas are vertically in the ice

  // double pos[3];
  std::vector<double> pos{0, 0, 0};
  std::vector<double> dir{0, 0, 0};         // Vector direction to cs.
  double dist = 0;             // Module of distance to cs.
  std::vector<double> sph_ang{0,0};

  double dot = 0;              // Dot (inner) product with cascade direction.
                               // The cosine of the angle between them.
  double ang = 0;              // The angle.


  double f_obs;               // [Hz] Observer frequency = 1E9 , 450 *1E6
  double k_obs;               // Wavenumber.
  double l_obs;               // [m] Detection wavelength
  double w_obs;               // Angular freq_obs

  double IRT_dist[2];
  double IRT_angle[2];
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

  std::vector<Antenna> transmitters();
  std::vector<Antenna> receivers();

  void add_antenna(Antenna& at);

private:

  std::vector<Antenna> _transmitters;
  std::vector<Antenna> _receivers;

};

#endif
