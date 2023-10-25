// ----- Physical parameters ---------------------------------------------------

  // --- Units -----------------------------------------------------------------

/* Inspired by the CLHEP global system of units.

  use these to keep your numbers in the global system of units defined above:
  ns, GHz, mm, nC
  if you want to write something in terms of MHz, for example, just do
  freq = 1200*MHz
  and then freq will have units of GHz, as it should.

      lengths
  for example, if you wanted to calculate the time it took for a signal
  to propagate 75 feet, you'd do:

  75*ft/c_light

  and it would return the correct time in nanoseconds.
*/

// length
static constexpr double mm = 1;
static constexpr double cm = 10*mm;
static constexpr double m = 1000.*mm;
//static constexpr double mm = .001*m;
// static constexpr double ft = .3047*m;

//energy
static constexpr double MeV = 1.;
static constexpr double EeV = 1e12*MeV;
static constexpr double PeV = 1e9*MeV;
static constexpr double GeV = 1000.*MeV;
static constexpr double KeV = .001*MeV;
static constexpr double eV = 1e-6*MeV;

  //time
static constexpr double ns = 1.;
static constexpr double us = ns*1e3;
static constexpr double ms = ns*1e6;
static constexpr double s = ns*1e9;

  //frequency
static constexpr double GHz = 1./ns;
static constexpr double THz = 1000.*GHz;
static constexpr double MHz = .001*GHz;
static constexpr double kHz = 1e-6*GHz;
static constexpr double Hz = 1e-9*GHz;

  //mass
static constexpr double g = 1.;
static constexpr double kg = 1000.*g;

  //work
static constexpr double W = 1;
// static constexpr double W = kg*pow(m,2)*pow(s,-3);
static constexpr double kW = 1000.*W;

// Universal Constants -------------------------------------------------------------

// ----- Math Constants --------------------------------------------------------
// static constexpr double pi=3.1415926535;            // OLD
static constexpr double pi = 3.14159265358979323846;     // [RS] radians
static constexpr double e = 2.71828;                    // Do not mix with electron charge
static constexpr double deg = pi/180.;     //radians
// static constexpr double degree=pi/180.;  //radians
// static constexpr double rad = 180./pi;     //radians
// static constexpr double twoPi = 2.*pi;   //radians
// static complex<double> I=sqrt(complex<double>(-1));//imaginary unit, used for complex stuff

// ------ Physics Constants --------------------------------
  // static constexpr double c_vac=2.998E8;                        // [m/s]
static constexpr double c_vac=2.9979246E8 *m/s;                    // [mm/ns]
static constexpr double Z_0=377;                                   // [Ohm]
// Thompson e- scattering cs
static constexpr double thomson=6.6524574E-25 *cm*cm;             // [cm^2]

// Not used, available
// static constexpr double m_e=0.510998;                           // [MeV/c^2]
// static constexpr double classic_electr_radius = 2.8179403E-15 *m;
// static constexpr double kelvin=1;                  // [K]
// static constexpr double z_0=50;                    // [Ohm]
// static constexpr double kB=8.617343e-11 *MeV/kelvin;  // [MeV/kelvin]
// static constexpr double kBJoulesKelvin=1.38e-23/kelvin;      // [J/kelvin]
