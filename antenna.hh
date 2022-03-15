/*//////////////////////////////////////////////////////////////////////////////
This header file defines the antenna positions.
*///////////////////////////////////////////////////////////////////////////////
#ifndef ANTENNA_hh
#define ANTENNA_hh

#include "macro_settings.hh"

class Antenna{

public:
  Antenna(double power, double xpos, double ypos, double zpos);

  Antenna(double power, double xpos, double ypos, double zpos,
          double xpol, double ypol, double zpol, double gain, double frequency);

  double GainDipole(double theta);
  double GainDipole(std::vector<double> dir);

  // Accesors

  double Power()  const;
  double Gain()   const;
  double E0()     const;    // [V/m]
  double Leff()   const;
  double Eff()    const;
  double Load()   const;

  double Freq()         const;  // [Hz]
  double Wavelength()   const;  // [m]
  double Wavenumber()   const;  // [rad/m]
  double AngularFreq()  const;  // [Hz * rad]

  std::vector<double> Pos() const;
  std::vector<double> Pol() const;

  double  Dist() const;
  std::vector<double> Dir() const;
  std::vector<double> Sph() const;

  double Dot()   const;
  double Delta() const;

  void SetDirection(std::vector<double> coords);
  void SetAngle(std::vector<double> dir);

private:
  // friend class Scatter;
  friend class Detector;

  double fPower;          // [W]
  double fGain;
  double fLeff = 1;       // m
  double fEff  = 1;
  double fLoad = 50;      // [pi*Ohm]

  // These are all related. Setting or changeing one needs to redo the rest.
  double fFrequency;               // [Hz] Observer frequency = 1E9 , 450 *1E
  double fWavelength;              // [m] Detection wavelength
  double fWavenumber;              // Wavenumber.
  double fAngularFreq;             // Angular freq_obs

  std::vector<double> fPosition{0, 0, 0};
  // Standard antennas are vertically in the ice (Vpol)
  std::vector<double> fPolarization{1, 0, 0};   // Polarization

  double fDistance = 0;                         // Module of distance to cs point.
  std::vector<double> fDirection{0, 0, 0};      // Vector direction to cs point.
  std::vector<double> fSphericalAngles{0,0};

  double fDot = 0;             // Dot (inner) product with reference (cascade) direction.
  // This is the same as the cosine of the angle between them.
  double fDelta = 0;          // The angle.
  // fDelta is the projection angle, defined between 0 and pi only.

  // double IRT_dist[2];
  // double IRT_angle[2];
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

  std::vector<Antenna> Transmitters();
  std::vector<Antenna> Receivers();

  void add_antenna(Antenna& at);
  void add_antenna(Antenna at);

private:

  std::vector<Antenna> fTransmitters;
  std::vector<Antenna> fReceivers;

};

#endif
