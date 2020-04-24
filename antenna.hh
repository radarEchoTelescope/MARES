/*//////////////////////////////////////////////////////////////////////////////
This header file defines the antenna positions.
*///////////////////////////////////////////////////////////////////////////////
#ifndef ANTENNA_hh
#define ANTENNA_hh

#include "macro_settings.hh"
#include "cascade.hh"
#include "IceRayTracing.hh"
//using namespace IceRayTracing;  // Call the namespace directly instead
// class Cascade;

class Antenna{

public:

  Antenna(double power, double position, double polarization, double gain, double frequency);
  Antenna(double power, double position);

  void set_direction(const Cascade & cs);
  void IRT_set_direction(const Cascade & cs);

  double get_power();
  double get_gain();
  double* get_position();
  double* get_polarization();

  double get_f_obs();
  double get_l_obs();
  double get_w_obs();
  double get_k_obs();

  double* get_cs_direction();
  double get_cs_distance();
  double get_cs_projection();

  double* get_IRT_distances();
  double* get_IRT_angles();

private:
  double power;
  double gain;

  double pos[3];
  double polar[2];    // Polarization

  double cs_dir[3];           // Vector direction to cs.
  double cs_dist;             // Module of distance to cs.
  double cs_dot;              // Dot (inner) product with cascade direction.
                              // The cosine of the angle between them.
  double IRT_cs_dist[2];
  double IRT_cs_angle[2];

  double f_obs;               // [Hz] Observer frequency = 1E9 , 450 *1E6
  double k_obs;               // Wavenumber.
  double l_obs;               // [m] Detection wavelength
  double w_obs;               // Angular freq_obs

};


class Detector{

public:

// Consturctors
  Detector();         // Defalut, bistatic
  Detector(string det_name); // Model based - RET_CR

  std::vector<Antenna> get_transmitters();
  std::vector<Antenna> get_receivers();

private:

  std::vector<Antenna> Transmitters;
  std::vector<Antenna> Receivers;

};

// // Unnamed namespaces FTW.
// namespace{
// double cosine_rule(double a, double b, double cos_alpha){
//   double rose_hates_this_lack_of_line;
//   rose_hates_this_lack_of_line = sqrt(pow(a,2) + pow(b,2) - 2*a*b*cos_alpha);
//   return rose_hates_this_lack_of_line;
// }
//
// };

#endif

//  Other Detector configurations

/* Old setup locations, Krijn's code.

switch (ii) {  // Position, 5 positions
  case 1: break;                    // Center (default)
  case 2: tx.pos.y = 1E3;  break;   // Top
  case 3: tx.pos.x = -1E3; break;   // Left
  case 4: tx.pos.y = -1E3; break;   // Bottom
  case 5: tx.pos.x = 1E3;  break;   // Right
}

switch (iii) {  // Position, 5 positions
    case 1:  rx.pos.y =   5E2; break;
    case 2:  rx.pos.y =  -5E2; break;
    case 3:  rx.pos.x =   5E2; break;
    case 4:  rx.pos.x =  -5E2; break;
    case 5:  rx.pos.y =  15E2; break;
    case 6:  rx.pos.x =   5E2; rx.pos.y =  10E2; break;
    case 7:  rx.pos.x =  -5E2; rx.pos.y =  10E2; break;
    case 8:  rx.pos.x = -15E2; break;
    case 9:  rx.pos.x = -10E2; rx.pos.y =   5E2; break;
    case 10: rx.pos.x = -10E2; rx.pos.y =  -5E2; break;
    case 11: rx.pos.y = -15E2; break;
    case 12: rx.pos.x =   5E2; rx.pos.y = -10E2; break;
    case 13: rx.pos.x =  -5E2; rx.pos.y = -10E2; break;
    case 14: rx.pos.x =  15E2; break;
    case 15: rx.pos.x =  10E2; rx.pos.y =  5E2; break;
    case 16: rx.pos.x =  10E2; rx.pos.y = -5E2;  break;
    case 17: rx.pos.x =  20E2; break;
  }

switch (jj) {              // Polarization angles, 2 options
  case 1: {
  tx.pos.phi = 0, tx.pos.theta = pi/2;
  rx.pos.phi = 0, rx.pos.theta = pi/2;
  break;}
  case 2: break; // default case, defined as 0
}
*/


/* Simon's GEANT RNOg model

GEANT 4 Array Style, position of the antennas.

Each sub vector in the antennas vector hold the x, y and z position of an antenna.
Remember: the XZ plane is the horizontal plane (ice surface plane),
Y direction is vertically upwards, origin (0, 0, 0) is in the middle of the
 ice shelf (so shelfSizeY/2 is BELOW the surface).
e.g. (0, shelfSizeY/2., 0) lies in the middle of the ice shelf surface,
an Y pos of -20.*m corresponds to 30 m below the surface

static const std::vector<std::vector<G4double>> antennas =
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
*/
