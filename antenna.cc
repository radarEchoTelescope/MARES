#include "antenna.hh"

Antenna::Antenna(double power, double xpos, double ypos, double zpos,
                 double theta, double phi, double gain, double frequency):
  power_(power), pos{xpos, ypos, zpos}, polar{theta, phi}, f_obs(frequency){

  k_obs = w_obs/c_med;
  l_obs = c_med/f_obs;
  w_obs = 2*pi *f_obs;

  if (!gain){ power ? gain_ = 1 : gain_ = pow(l_obs,2.0); }
  /* If gain = 0, default , unphysical value, then set default gains
  Default TX gain (power nonzero) =  1, Isotropic emission.
  Defalut RX gain (power zero): Effective area. */
}

Antenna::Antenna(double power, double xpos, double ypos, double zpos) :
Antenna(power, xpos, ypos, zpos, 0,0, 0, freq_obs){} //default freq_obs

// Accesors

double  Antenna::power(){return power_;}
double  Antenna::gain(){return gain_;}
double* Antenna::position(){return pos;}
double* Antenna::polarization(){return polar;}
double  Antenna::distance(){return dist;}
double  Antenna::projection(){return dot;}
double  Antenna::frequency(){return f_obs;}
double  Antenna::lambda(){return l_obs;}
double  Antenna::omega(){return w_obs;}
double  Antenna::wavenr(){return k_obs;}

// --------------------------------------------

Detector::Detector(){};

Detector::Detector(std::string name){
  if(name == "bistatic" || name == "Bistatic"){
    Transmitters.push_back( Antenna(1, 0,0,0) );
    Receivers.push_back( Antenna(0, 500, 0, 0) );
  } else if (name == "RET_CR" || name == "ret_cr"){
    Transmitters.push_back( Antenna(1, 0,0,0) );
    Receivers.push_back( Antenna(0, 200,0,0) );
    Receivers.push_back( Antenna(0, 100, -100, 0) );
    Receivers.push_back( Antenna(0, -100, -100, 0) );
    Receivers.push_back( Antenna(0, -100,  100, 0) );
    Receivers.push_back( Antenna(0, 200 , 200, 0) );
    Receivers.push_back( Antenna(0, 200 ,-200, 0) );
    Receivers.push_back( Antenna(0, -200 ,-200, 0) );
    Receivers.push_back( Antenna(0, -200 , 200, 0) );
  } else {
    std::cout << " There is no compatible detector configurtion with "
    + name << std::endl;
  }
};
//
// Detector::Detector(int a, int b){
//   switch (a) {  // 5 options
//     case 1: Transmitters.push_back(Antenna(1, (0, 0, 0), (0,0), 0, freq_obs));
//       break;  // Center
//     case 2: Transmitters.push_back(Antenna(1, (0, 1E3, 0), (0,0), 0, freq_obs));
//       break;   // Top
//     case 3: Transmitters.push_back(Antenna(1, (-1E3, 0, 0), (0,0), 0, freq_obs));
//       break;   // Left
//     case 4: Transmitters.push_back(Antenna(1, (0, -1E3, 0), (0,0), 0, freq_obs));
//       break;   // Bottom
//     case 5: Transmitters.push_back(Antenna(1, (1E3, 0, 0), (0,0), 0, freq_obs));
//       break;   // Right
//   }
//
//   switch (b) {  // Position, 5 positions
//     case 1: Receivers.push_back(Antenna(0, (0, 5E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 2: Receivers.push_back(Antenna(0, (0, -5E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 3: Receivers.push_back(Antenna(0, (5E2, 0, 0), (0,0), 0, freq_obs));
//       break;
//     case 4: Receivers.push_back(Antenna(0, (-5E2, 0, 0), (0,0), 0, freq_obs));
//       break;
//     case 5: Receivers.push_back(Antenna(0, (0, 15E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 6: Receivers.push_back(Antenna(0, (5E2, 1E3, 0), (0,0), 0, freq_obs));
//       break;
//     case 7: Receivers.push_back(Antenna(0, (-5E2, 1E3, 0), (0,0), 0, freq_obs));
//       break;
//     case 8: Receivers.push_back(Antenna(0, (-15E2, 0, 0), (0,0), 0, freq_obs));
//       break;
//     case 9: Receivers.push_back(Antenna(0, (-1E3, 5E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 10: Receivers.push_back(Antenna(0, (-1E3, -5E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 11: Receivers.push_back(Antenna(0, (0, -15E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 12: Receivers.push_back(Antenna(0, (5E2, -1E3, 0), (0,0), 0, freq_obs));
//       break;
//     case 13: Receivers.push_back(Antenna(0, (-5E2, -1E3, 0), (0,0), 0, freq_obs));
//       break;
//     case 14: Receivers.push_back(Antenna(0, (15E2, 0, 0), (0,0), 0, freq_obs));
//       break;
//     case 15: Receivers.push_back(Antenna(0, (1E3, 5E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 16: Receivers.push_back(Antenna(0, (1E3, -5E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 17: Receivers.push_back(Antenna(0, (2E3, 0, 0), (0,0), 0, freq_obs));
//       break;
//     }
//
// };

std::vector<Antenna> Detector::get_transmitters(){return Transmitters;}
std::vector<Antenna> Detector::get_receivers(){return Receivers;}

void Detector::add_antenna(Antenna& at){
  if (at.power()) { Transmitters.push_back(at); }
  else { Receivers.push_back(at); }
}

// -----------------------------------------------------------------------------
//  Other Detector configurations

/* Old setup locations, Krijn's code.



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
