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

// ----- Computational constants -----------------------------------------------
const int nbin = 3000;                 // # of bins for integration/array filling
// const int nbin = 1E4;                 // # of bins for integration/array filling

//---- I/O functions -----------------------------------------------------------
void write_1D_array(std::vector<double> array, std::string output_path,
  const bool trigger);

void write_2D_array(std::vector<std::vector<double>> array,
    std::string output_path, const bool trigger);

// ----- Physical constants ----------------------------------------------------

// ----- Macroscopical [IS]
const double pi=3.1415926535;
// const double c_vac=3E8;               // [m/s]
const double c_vac=2.998E8;               // [m/s]
// const double refindex=1.7;            // refractive index in ice
const double refindex=1.78;            // refractive index in ice
const double c_med=c_vac/refindex;    // [m/s] speed of light in ice

// ----- Microscopical [ cgs]
const double cvac_cm=c_vac * 100;      // [cm/s]
const double cmed_cm=c_med * 100;      // [cm/s] c in ice
const double thompson=6.6524574E-25;   // [cm^2] Thompson e- scattering cs

      // Ice parameters --------------------------------------------------------
const double rho_ice=0.92;             // [g/cm^3] Density
const double E_c= 0.0786;              // [GeV] Critical cascade energy
const double X_0 = 36.08;              // [g/cm^2] radiation length
const double X_int=25.01;              // [g/cm^2] interaction length
const double w_coll=1E12;             // [rad/s] collision frequency
// const double w_coll=0;             // [rad/s] collision frequency

    // Plasma parameters--------------------------------------------------------
const double r_moliere=7;              // [cm] Moliere Radius
const double tau=20*1E-9;              // [s] Plasma lifetime
const double mme=1;                    // Plasma to electron mass ratio

// Math functions --------------------------------------------------------------

int sgn(double val);
// template <typename T>  int sgn(T val) { return (T(0) < val) - (val < T(0));}

double distance(double x1, double y1, double x2, double y2);
double distance(double x1, double y1, double z1, double x2, double y2, double z2);

std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> matrix);

  // Random numbers

class RN_uniform{

public:
  RN_uniform (double minval, double maxval, double seed = 42) :
      generator(seed), distribution(minval, maxval) {}

  double get() { return distribution(generator); }

private:
  std::mt19937 generator;
  std::uniform_real_distribution<> distribution;
};

#endif
