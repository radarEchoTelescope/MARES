#ifndef SETTINGS_hh
#define SETTINGS_hh

// #define NDEBUG     // Turn off debug.

#include <iostream>   // C++ only, file read/write
#include <fstream>    // C++ only, screen read/write
#include <stdlib.h>   // exit, EXIT_FAILURE
#include <string.h>
#include <array>
#include <vector>
#include <iterator>
#include <algorithm>  // Copy, max, min, etc.
#include <math.h>     // Math stuff: pow,sqrt,etc.
#include <assert.h>   // Debug purposes
#include <random>

#include <libconfig.h++>
#include <settings_units_constants.hh>

// Simulation parameters --------------------------------------------    
// see "settings_params.hh"

      // Computational parameters  --------------------------------------------

extern double _dL; // [cm/bin]
extern double _dR;
extern double _sampling;
extern double _Ltot_factor;
extern double _Rtot_factor;

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


extern double _c_ice;  // [m/s] speed of light in ice
extern double _Z_ice;  // Impedance of ice
extern double _L_0;  // Radiation length = 39.22 cm

    // Air / Other constants -----------------------------------
// extern double rho=1.168e-3;//sea level density
// extern double x_0=36.7;//radiation length in air

// The total cascade dimensions are given in terms of the 
// the typical length scale in every dimension. 
// Total = typical* factor
// 	The typical length of a cascade is log(12.72 * fEnergy) * X_0;	
// 	The typical radius of a cascade is the moliere radius;	

static const double Ltot_factor = _Ltot_factor;
static const double Rtot_factor = _Rtot_factor;

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
double projection_pol(std::vector<double> a, std::vector<double> b);

double shortest_distance(std::vector<double> plane, std::vector<double> point);

std::vector<double> normalize(std::vector<double> a);
std::vector<double> cross_product(double a1, double a2, double a3, double b1, double b2, double b3);
std::vector<double> cross_product(std::vector<double> a, std::vector<double> b);
double dot_product(std::vector<double> a, std::vector<double> b);
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
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v){
    os << "[";
    for (int i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i != v.size() - 1){os << ", ";}
    }
    os << "]\n";
    return os;
};

// C++ template to print vector container elements
template <typename T,std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T,N>& v){
    os << "[";
    for (int i = 0; i < N; ++i) {
        os << v[i];
        if (i != N - 1) {os << ", ";}
    }
    os << "]\n";
    return os;
};

void write_1D_array(std::vector<double> array, std::string output_path, const bool trigger = true);

void write_2D_array(std::vector<std::vector<double>> array, std::string output_path, const bool trigger = true);

// Read a libconfig file and return the configuration. 
// If there is an error, report it and exit early.
void load_config_file(libconfig::Config& cfg, const char* config_file);

#endif
