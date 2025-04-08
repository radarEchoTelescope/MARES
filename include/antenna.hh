/*//////////////////////////////////////////////////////////////////////////////
This header file defines the antenna positions.
*///////////////////////////////////////////////////////////////////////////////
#ifndef ANTENNA_hh
#define ANTENNA_hh

#include "settings.hh"

class Antenna{

public:
  Antenna(double xpos, double ypos, double zpos,
          double xpol, double ypol, double zpol,
          double power, double frequency, double gaindB);

// Using RS approach now to define the gain in dB, but used as linear.
// 0  dBi = 1
// 1.76 dBi ~= 1.5, ideal small dipole
// 2.15 dBi = 1.64, half-wavelength dipole
// 3  dBi = 2
// 5  dBi = 7   *** wrong way round i think - 7 dBi = 5x
// 10 dBi = 10

  void SetDirection(std::vector<double> coords);
  void SetAngle(std::vector<double> dir);
  double GainDipole(double theta);

  // NOT ADDED YET
  // double GainDipole(std::vector<double> dir);

  // Accesors

  double Power()  const;
  double Gain()   const;
  double E0()     const;    // [V/m]
  double Leff()   const;
  double Load()   const;
  double Efficiency()    const;

  double Freq()         const;  // [Hz]
  double Wavelength()   const;  // [m]
  double Wavenumber()   const;  // [rad/m]
  double AngularFreq()  const;  // [Hz * rad]

  std::vector<double> Pos() const;
  std::vector<double> Pol() const;

// W.r.t the interaction vertex / direction of interest
  double Dot()   const;
  double Delta() const;
  double  Dist() const;
  std::vector<double> Dir() const;
  std::vector<double> Sph() const;

private:
  friend class Detector;

  double fPower;                // [W]
  double fGain;
  double fLeff = 1*m;           // [mm]
  double fEff  = 1;
  double fLoad = 50;            // [Ohm]

  // These are all related. Setting or changeing one needs to redo the rest.
  double fFrequency;               // [Hz] Observer frequency = 1E9 , 450 *1E
  double fWavelength;              // [m] Detection wavelength
  double fWavenumber;              // Wavenumber.
  double fAngularFreq;             // Angular freq_obs

  std::vector<double> fPosition{};
  // Standard antennas are vertically in the ice (Vpol)
  std::vector<double> fPolarization{};   // Polarization

  double fDistance = 0;                         // Module of distance to cs point.
  std::vector<double> fDirection{0, 0, 0};      // Vector direction to cs point.
  std::vector<double> fSphericalAngles{0,0};

  double fDot = 0;             // Dot (inner) product with reference (cascade) direction.
  // This is the same as the cosine of the angle between them.
  double fDelta = 0;          // The angle.
  // fDelta is the projection angle, defined between 0 and pi only.

};

class Detector{

public:

// Constructors
  Detector();               // Default, empty
  Detector(std::string name);   // Model-based constructor
  /* Other models:
    - Bistatic
    - T576
    - Star grid: Plane of receivers @-140 m, half-dipoles, 100 MHz
        You need to add your TX. 
  */

  std::vector<Antenna> Transmitters();
  std::vector<Antenna> Receivers();

  void add_antenna(Antenna& at);
  void add_antenna(Antenna at);


private:

  std::vector<Antenna> fTransmitters;
  std::vector<Antenna> fReceivers;

};

// Loader
Detector load_antenna_list(libconfig::Setting& at_list);
Detector load_detector_config(libconfig::Config& dect_config);
Detector load_detector_file(const std::string& dect_config_filepath);

//  Detector configurations can go in their own namespace,
// to not access the wrong one by mistake. 

namespace {
  /* GEANT 4 Array Style, position of the antennas.

  Each sub vector in the antennas vector hold the x, y and z position of an antenna.
  Remember: the XZ plane is the horizontal plane (ice surface plane),
  Y direction is vertically upwards, origin (0, 0, 0) is in the middle of the
  ice shelf (so shelfSizeY/2 is BELOW the surface).
  e.g. (0, shelfSizeY/2., 0) lies in the middle of the ice shelf surface,
  an Y pos of -20.*m corresponds to 30 m below the surface
  */
  static const std::vector<std::vector<double>> StarGrid =
  {
    {28.2842712475*m, -140.*m, 28.2842712475*m},
    {-28.2842712475*m, -140.*m, -28.2842712475*m},
    {28.2842712475*m, -140.*m, -28.2842712475*m},
    {-28.2842712475*m, -140.*m, 28.2842712475*m},
    {56.5685424949*m, -140.*m, 56.5685424949*m},
    {-56.5685424949*m, -140.*m, -56.5685424949*m},
    {56.5685424949*m, -140.*m, -56.5685424949*m},
    {-56.5685424949*m, -140.*m, 56.5685424949*m},
    {84.8528137424*m, -140.*m, 84.8528137424*m},
    {-84.8528137424*m, -140.*m, -84.8528137424*m},
    {84.8528137424*m, -140.*m, -84.8528137424*m},
    {-84.8528137424*m, -140.*m, 84.8528137424*m},
    {113.13708499*m, -140.*m, 113.13708499*m},
    {-113.13708499*m, -140.*m, -113.13708499*m},
    {113.13708499*m, -140.*m, -113.13708499*m},
    {-113.13708499*m, -140.*m, 113.13708499*m},
    {141.421356237*m, -140.*m, 141.421356237*m},
    {-141.421356237*m, -140.*m, -141.421356237*m},
    {141.421356237*m, -140.*m, -141.421356237*m},
    {-141.421356237*m, -140.*m, 141.421356237*m},
    {169.705627485*m, -140.*m, 169.705627485*m},
    {-169.705627485*m, -140.*m, -169.705627485*m},
    {169.705627485*m, -140.*m, -169.705627485*m},
    {-169.705627485*m, -140.*m, 169.705627485*m},
    {197.989898732*m, -140.*m, 197.989898732*m},
    {-197.989898732*m, -140.*m, -197.989898732*m},
    {197.989898732*m, -140.*m, -197.989898732*m},
    {-197.989898732*m, -140.*m, 197.989898732*m},
    {226.27416998*m, -140.*m, 226.27416998*m},
    {-226.27416998*m, -140.*m, -226.27416998*m},
    {226.27416998*m, -140.*m, -226.27416998*m},
    {-226.27416998*m, -140.*m, 226.27416998*m},
    {254.558441227*m, -140.*m, 254.558441227*m},
    {-254.558441227*m, -140.*m, -254.558441227*m},
    {254.558441227*m, -140.*m, -254.558441227*m},
    {-254.558441227*m, -140.*m, 254.558441227*m},
    {282.842712475*m, -140.*m, 282.842712475*m},
    {-282.842712475*m, -140.*m, -282.842712475*m},
    {282.842712475*m, -140.*m, -282.842712475*m},
    {-282.842712475*m, -140.*m, 282.842712475*m},
    {-300.0*m, -140.*m, 0.0*m},
    {-280.0*m, -140.*m, 0.0*m},
    {-260.0*m, -140.*m, 0.0*m},
    {-240.0*m, -140.*m, 0.0*m},
    {-220.0*m, -140.*m, 0.0*m},
    {-200.0*m, -140.*m, 0.0*m},
    {-180.0*m, -140.*m, 0.0*m},
    {-160.0*m, -140.*m, 0.0*m},
    {-140.0*m, -140.*m, 0.0*m},
    {-120.0*m, -140.*m, 0.0*m},
    {-100.0*m, -140.*m, 0.0*m},
    {-80.0*m, -140.*m, 0.0*m},
    {-60.0*m, -140.*m, 0.0*m},
    {-40.0*m, -140.*m, 0.0*m},
    {-20.0*m, -140.*m, 0.0*m},
    {0.0*m, -140.*m, -300.0*m},
    {0.0*m, -140.*m, -280.0*m},
    {0.0*m, -140.*m, -260.0*m},
    {0.0*m, -140.*m, -240.0*m},
    {0.0*m, -140.*m, -220.0*m},
    {0.0*m, -140.*m, -200.0*m},
    {0.0*m, -140.*m, -180.0*m},
    {0.0*m, -140.*m, -160.0*m},
    {0.0*m, -140.*m, -140.0*m},
    {0.0*m, -140.*m, -120.0*m},
    {0.0*m, -140.*m, -100.0*m},
    {0.0*m, -140.*m, -80.0*m},
    {0.0*m, -140.*m, -60.0*m},
    {0.0*m, -140.*m, -40.0*m},
    {0.0*m, -140.*m, -20.0*m},
    {0.0*m, -140.*m, 0.0*m},
    {0.0*m, -140.*m, 20.0*m},
    {0.0*m, -140.*m, 40.0*m},
    {0.0*m, -140.*m, 60.0*m},
    {0.0*m, -140.*m, 80.0*m},
    {0.0*m, -140.*m, 100.0*m},
    {0.0*m, -140.*m, 120.0*m},
    {0.0*m, -140.*m, 140.0*m},
    {0.0*m, -140.*m, 160.0*m},
    {0.0*m, -140.*m, 180.0*m},
    {0.0*m, -140.*m, 200.0*m},
    {0.0*m, -140.*m, 220.0*m},
    {0.0*m, -140.*m, 240.0*m},
    {0.0*m, -140.*m, 260.0*m},
    {0.0*m, -140.*m, 280.0*m},
    {0.0*m, -140.*m, 300.0*m},
    {20.0*m, -140.*m, 0.0*m},
    {40.0*m, -140.*m, 0.0*m},
    {60.0*m, -140.*m, 0.0*m},
    {80.0*m, -140.*m, 0.0*m},
    {100.0*m, -140.*m, 0.0*m},
    {120.0*m, -140.*m, 0.0*m},
    {140.0*m, -140.*m, 0.0*m},
    {160.0*m, -140.*m, 0.0*m},
    {180.0*m, -140.*m, 0.0*m},
    {200.0*m, -140.*m, 0.0*m},
    {220.0*m, -140.*m, 0.0*m},
    {240.0*m, -140.*m, 0.0*m},
    {260.0*m, -140.*m, 0.0*m},
    {280.0*m, -140.*m, 0.0*m},
    {300.0*m, -140.*m, 0.0*m}
  };
}


#endif
