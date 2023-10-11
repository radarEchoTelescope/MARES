#ifndef SETTINGS_hh
#define SETTINGS_hh

// #define NDEBUG     // Turn off debug.

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

// "Universal" Constants -------------------------------------------------------------

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

// Simulation parameters --------------------------------------------    
// see "settings_params.hh"

      // Computational parameters  --------------------------------------------

extern double _dL; // [cm/bin]
extern double _dR;
extern double _dN; // This is also a radial direction.

extern double _sampling;


    // Plasma ----------------------------------------------------------------
extern double _lifetime;         // Mean plasma lifetime
extern double _f_coll;           // RS collision frequency
extern double _memp;             // Effective plasma mass in electron masses
// Effective mass plasma to electron mass ratio

    // Ice  ------------------------------------------------------------------
extern double _refindex;         // refractive index
extern double _att_length;       // Attenuation length
extern double _rho_ice;          // Density


// Mass stopping power of ice - energy loss per ionizing particle (@ 1 GeV)
extern double _r_moliere;        // Moliere Radius in ice
extern double _E_ionization;     // e- ionization energy [RS]
extern double _E_deposition;       
extern double _E_c;              // Critical cascade energy for ionization
extern double _X_0;              // Radiation columm density

static double _c_ice = c_vac/_refindex;  // [m/s] speed of light in ice
static double _Z_ice = Z_0/_refindex;            // Impedance of ice
static double _L_0 = _X_0/_rho_ice;       // Radiation length = 39.22 cm

    // Air / Other constants -----------------------------------
// extern double rho=1.168e-3;//sea level density
// extern double x_0=36.7;//radiation length in air

/*"Hiding" the mutable values and accessing through a constant reference
protects the code against accidental changes in a global variable. */
static const double& tau = _lifetime;
static const double& f_coll = _f_coll;
static const double& memp = _memp;
static const double& refindex = _refindex; 
static const double& att_length = _att_length; 
static const double& rho_ice = _rho_ice;
static const double& r_moliere = _r_moliere; 
static const double& E_ionization = _E_ionization; 
static const double& E_deposition = _E_deposition;       
static const double& E_c = _E_c;
static const double& X_0 = _X_0;
static const double& c_ice = _c_ice;
static const double& Z_ice = _Z_ice;
static const double& L_0 = _L_0;

// Math tools ------------------------------------------------------------------

template<typename T>
std::vector<T> arange(T start, T stop, T step ){
  std::vector<T> values;
  for (T value = start; value <= stop; value += step) { values.push_back(value);}
  return values;
}

template<typename T>
std::vector<T>& operator+=(std::vector<T> &lhs, const std::vector<T> &rhs) {
    if (lhs.size() != rhs.size())
        throw std::length_error("vectors must be same size to add");
    for (auto i = 0; i < lhs.size(); ++i) {lhs[i] += rhs[i];}
    return lhs;
}

template<typename T>
std::vector<T> operator+ (std::vector<T> lhs, const std::vector<T> &rhs) {  return lhs += rhs; }

template<typename T>
std::vector<T>& operator*=(std::vector<T> &lhs, const T &rhs) {
    for (auto i = 0; i < lhs.size(); ++i) {lhs[i] *= rhs;}
    return lhs;
}

template<typename T>
std::vector<T> operator* (std::vector<T> lhs, const T &rhs) {  return lhs *= rhs; }
template<typename T>
std::vector<T> operator* (const T &rhs, std::vector<T> lhs) {  return lhs *= rhs; }

int sgn(double val);

double acos2(double y, double x);
double sinc(double x);

double norm(std::vector<double> a);
std::vector<double> direction(std::vector<double> vec_a, std::vector<double> vec_b);
std::vector<double> direction(double x1, double y1, double z1, double x2, double y2, double z2);

double distance(std::vector<double> u, std::vector<double> v);
double distance(double x1, double y1, double x2, double y2);
double distance(double x1, double y1, double z1, double x2, double y2, double z2);
double projection(std::vector<double> a, std::vector<double> b);

double shortest_distance(std::vector<double> plane, std::vector<double> point);

std::vector<double> normalize(std::vector<double> a);
std::vector<double> cross_product(double a1, double a2, double a3, double b1, double b2, double b3);
std::vector<double> cross_product(std::vector<double> a, std::vector<double> b);
std::vector<double> rotate(std::vector<double> vector_in, std::vector<double> rotAxis, double angle);
std::vector<double> find_perpendicular(std::vector<double> vector_in);
std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> matrix);

  // Random number generator
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
};

void write_1D_array(std::vector<double> array, std::string output_path, const bool trigger);

void write_2D_array(std::vector<std::vector<double>> array, std::string output_path, const bool trigger);

#endif
