#ifndef MACRO_SETTINGS_hh
#define MACRO_SETTINGS_hh

#include <iostream>   // C++ only, file read/write
#include <fstream>    // C++ only, screen read/write
#include <string.h>
#include <vector>
#include <iterator>
#include <algorithm>  // Copy, max, min, etc.
#include <math.h>     // Math stuff: pow,sqrt,etc.
#include <assert.h>   // Debug purposes
#include <random>
// #define NDEBUG     // Turn off debug.

// ----- Computational parameters ----------------------------------------------
// const int nbin = 3000;               // # of bins for integration/array filling
// const int nbin = 1E4;             // # of bins for integration/array filling

const double gdPerp = 5.0; // [cm/bin]
const double gdPar  = 0.1;
const double gdNorm = 0.1; // This is also a radial direction.

// const double freq_obs = 450 * 1E6;
const double freq_obs = 5E8;              // [Hz] Default observer frequency
const double freq_sampling = 5E9;
// sampling frequency should be between 10x and 100x freq_obs.

// ----- Physical parameters ---------------------------------------------------

    // Plasma ------------------------------------------------------------------
const double tau = 1E-8;               // [s] Plasma lifetime 10 ns
// const double tau = 20*1E-9;               // [s] Plasma lifetime

// const double f_coll = 0;             // [Hz] collision frequency
// const double f_coll=88E12;              // [Hz] collision frequency
const double f_coll= 64.733E12;              // [Hz] collision frequency
const double memp = 1;                    // Plasma to electron mass ratio

      // Ice  ------------------------------------------------------------------
const double att_length=1450;             // [m] attenuation length
const double refindex=1.78;               // refractive index
const double rho_ice = 0.92;              // [g/cm^3] Density, from GEANT
const double r_moliere = 7;               // [cm] Moliere Radius
const double E_c = 0.0786;                // [GeV] Critical cascade energy
const double X_0 = 36.08;                 // [g/cm^2] radiation columm density
const double L_0 = X_0/rho_ice;           // [cm] = 39.22 radiation length

// ----- Physical constants ----------------------------------------------------
const double pi=3.1415926535;
const double e =2.71828;

    // Macroscopical [IS] ------------------------------------------------------
const double c_vac=2.998E8;               // [m/s]
const double c_ice=c_vac/refindex;        // [m/s] speed of light in ice
const double Z0=119.917 ;                 // [pi * Ohm]

    // Microscopical [cgs] -----------------------------------------------------
const double cvac_cm=c_vac * 100;         // [cm/s]
const double cice_cm=c_ice * 100;         // [cm/s] c in ice
const double thompson=6.6524574E-25;      // [cm^2] Thompson e- scattering cs


// Math tools ------------------------------------------------------------------

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
  RN_uniform (double minval, double maxval, double seed = 42) :
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

#endif

// Unused ----------------------------------------------------------------------
// std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> matrix);
