#ifndef MACRO_SETTINGS_hh
#define MACRO_SETTINGS_hh

#include <iostream>   // C++ only, file read/write
#include <fstream>    // C++ only, screen read/write
#include <stdlib.h>   // exit, EXIT_FAILURE
#include <string.h>
#include <vector>
#include <iterator>
#include <algorithm>  // Copy, max, min, etc.
#include <math.h>     // Math stuff: pow,sqrt,etc.
#include <assert.h>   // Debug purposes
#include <random>
// #define NDEBUG     // Turn off debug.

// Add flags here?
// use_NKG
// use_electron_plasma

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
// static constexpr double kg = 1000*g;

  // "Universal" Constants -------------------------------------------------------------

// You should not need multiple definitions of c_vac
  // static constexpr double c_vac=2.998E8;                        // [m/s]
static constexpr double c_vac=2.9979246E8 *m/s;                    // [mm/ns]
// static constexpr double c_vac_cm=c_vac /cm;                     // [cm/s]
static constexpr double Z0=377;                                   // [Ohm]
// static constexpr double classic_electr_radius = 2.8179403E-15 *m;
// Thompson e- scattering cs
static constexpr double thomson=6.6524574E-25 *cm*cm;             // [cm^2]
// static constexpr double m_e=0.510998;                           // [MeV/c^2]

// static constexpr double kelvin=1;                  // [K]
// static constexpr double z_0=50;                    // [Ohm]
// static constexpr double kB=8.617343e-11 *MeV/kelvin;  // [MeV/kelvin]
// static constexpr double kBJoulesKelvin=1.38e-23/kelvin;      // [J/kelvin]

      // Plasma ----------------------------------------------------------------
static constexpr double f_coll= 64.733 *THz;         // [Hz] RS collision frequency
// static constexpr double f_coll=88 *THz;           // [Hz] collision frequency
static constexpr double memp = 1;                    // Plasma to electron mass ratio

      // Ice  ------------------------------------------------------------------
static constexpr double refindex=1.78;               // refractive index
static constexpr double att_length=1450 *m;          // Attenuation length
static constexpr double c_ice=c_vac/refindex;        // [m/s] speed of light in ice
// static constexpr double c_ice_cm=c_ice /cm;       // [cm/s] c in ice
static constexpr double rho_ice = 0.92 *g/pow(cm,3); // Density, from GEANT
static constexpr double Z_ice=c_vac/refindex;        // Impedance of ice
static constexpr double r_moliere = 7 *cm;           // Moliere Radius in ice

// static constexpr double E_ionization = 20 *eV;    // e- ionization energy 
static constexpr double E_ionization = 69 *eV;       // e- ionization energy [RS]
// Mass stopping power of ice - energy loss per ionizing particle (@ 1 GeV)
static constexpr double E_deposition = 2 *MeV/g*pow(cm,2);       
static constexpr double E_c = 0.0786 * GeV;          // Critical cascade energy for ionization
static constexpr double X_0 = 36.08 *g/pow(cm,2);    // Radiation columm density
static constexpr double L_0 = X_0/rho_ice;           // Radiation length = 39.22 cm

/* Other RS constants

  static constexpr double rho=1.168e-3;//sea level density
  static constexpr double x_0=36.7;//radiation length in air

*/



// ----- Math Constants --------------------------------------------------------
// static constexpr double pi=3.1415926535;            // OLD
static constexpr double pi=3.14159265358979323846;     // [RS] radians
static constexpr double e =2.71828;                    // Do not mix with electron charge
// static constexpr double deg=pi/180.;     //radians
// static constexpr double degree=pi/180.;  //radians
// static constexpr double twoPi = 2.*pi;   //radians
// static complex<double> I=sqrt(complex<double>(-1));//imaginary unit, used for complex stuff

// Math tools ------------------------------------------------------------------

template<typename T>
std::vector<T> arange(T start, T stop, T step ){
  std::vector<T> values;
  for (T value = start; value <= stop; value += step) { values.push_back(value);}
  return values;
}


int sgn(double val);

// template <class T, class Q>
// std::vector <T> operator* (const Q c, std::vector <T> A)
// {
//     std::transform (A.begin (), A.end (), A.begin (),
//                  std::bind1st (std::multiplies <T> () , c)) ;
//     return A ;
// }

// C++ template to print vector container elements
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    os << "[";
    for (int i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i != v.size() - 1)
            os << ", ";
    }
    os << "]\n";
    return os;
}

template<typename N>
N rad2deg(N angle) {
  return angle * 180.0 / pi;
}

template<typename N>
N deg2rad(N radian) {
  return radian / 180.0 * pi;
}

double acos2(double y, double x);
double sinc(double x);

double norm(std::vector<double> a);
std::vector<double> direction(std::vector<double> vec_a, std::vector<double> vec_b);
std::vector<double> direction(double x1, double y1, double z1, double x2, double y2, double z2);

double distance(std::vector<double> u, std::vector<double> v);
double distance(double x1, double y1, double x2, double y2);
double distance(double x1, double y1, double z1, double x2, double y2, double z2);
double projection(std::vector<double> a, std::vector<double> b);

std::vector<double> normalize(std::vector<double> a);
std::vector<double> cross_product(double a1, double a2, double a3, double b1, double b2, double b3);
std::vector<double> cross_product(std::vector<double> a, std::vector<double> b);

  // Random numbers

class RN_uniform{

public:
  RN_uniform (double minval, double maxval, double seed) :
      _generator(seed), _distribution(minval, maxval) {}

  double get() { return _distribution(_generator); }

private:
  std::mt19937 _generator;
  std::uniform_real_distribution<> _distribution;
};


//---- I/O functions -----------------------------------------------------------
void write_1D_array(std::vector<double> array, std::string output_path,
  const bool trigger);

void write_2D_array(std::vector<std::vector<double>> array,
    std::string output_path, const bool trigger);



// Unused ----------------------------------------------------------------------
std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> matrix);



#endif
