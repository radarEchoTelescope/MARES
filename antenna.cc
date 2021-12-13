#include "antenna.hh"

Antenna::Antenna(double power, double xpos, double ypos, double zpos) :
Antenna(power, xpos, ypos, zpos, 0,0,1, 0, freq_obs){} //default freq_obs

Antenna::Antenna(double power, double xpos, double ypos, double zpos,
                 double xpol, double ypol, double zpol, double gain, double frequency):
                fPower(power),
                fFrequency(frequency),
                fWavelength(c_ice/frequency),
                fAngularFreq(2 * pi * frequency),
                fWavenumber(2 * pi * frequency / c_ice),
                fPosition{xpos, ypos, zpos},
                fPolarization{xpol, ypol, zpol}{

  if (!gain){ fPower ? fGain = 1 : fGain = pow(fWavelength,2.0); }
  /* If gain = 0, default , unphysical value, then set default gains
  Default TX gain (nonzero power) =  1, Isotropic emission.
  Defalut RX gain (zero power): Effective area. */
}

void Antenna::SetDirection(std::vector<double> coords){

  // Direction to cascade: Correct orientation.TX -> CS -> RX.
  fPower != 0 ?  fDirection = direction(fPosition, coords) :
                 fDirection = direction(coords, fPosition);

  // Module of distance to cascade
  fDistance = norm(fDirection);

  /* Quick fix, the far field ssumption should be much larger */
  if(fDistance < fWavelength) {fDistance = fWavelength;}

  // Sanity check
  // assert(fDistance != 0 && "Point overlaps antenna");

 // Including correction for sgn(angle);
 fSphericalAngles[0] =  acos( fDirection[2] / fDistance);
 // The zenith angle is defined between 0 and pi so its sign safe.

  // Spherical angles of the line of sight to cascade
  fSphericalAngles[1] = atan2( fDirection[1], fDirection[0]);
  // atan2 has built in corrections for the signs of the angle.

}

void Antenna::SetAngle(std::vector<double> dir){
  // Determine inner product  in l.o.s plane.
  fDot   = projection( fDirection, dir );
  // fDelta is the projection angle, defined between 0 and pi only.
  fDelta = acos(fDot);

  /* The current model breaks down at small angles when
  |L*sin(delta)| < |r*cos(delta)|
  or
  r/L = |tan(delta)|
  so
  delta_critical = arctan(r/L)
  */

  // double delta_crit = atan(cs.Rtot()/cs.Ltot());
  // std::cout << rad2deg(delta_crit) << std::endl;
  // fDelta < delta_crit ? fDelta = delta_crit: 1;
  // This is right now 0.5 degrees

  // Or, something simpler, if delta is smaller than 1 degree, make it 1 degree.
  // (fDelta < pi/180.0) ? fDelta = pi/180.0 : 1;


}

double Antenna::GainDipole(double theta){
  return 1.643 * pow(sin(theta),2.6);
}

// double Antenna::GainDipole(std::vector<double> dir){
//
// };

// Accesors

double  Antenna::Power()  const {return fPower;}
double  Antenna::Gain()   const {return fGain;}
double  Antenna::E0()     const {return fWavenumber * sqrt(Z0*fPower/fGain);}
double  Antenna::Leff()   const {return fLeff;}
double  Antenna::Eff()    const {return fEff;}
double  Antenna::Load()   const {return fLoad;}

double  Antenna::Freq()       const {return fFrequency;}
double  Antenna::Wavelength() const {return fWavelength;}
double  Antenna::Wavenumber() const {return fWavenumber;}
double  Antenna::AngularFreq()const {return fAngularFreq;}

std::vector<double> Antenna::Pos() const {return fPosition;}
std::vector<double> Antenna::Dir() const {return fDirection;}
std::vector<double> Antenna::Sph() const {return fSphericalAngles;}
std::vector<double> Antenna::Pol() const {return fPolarization;}

double  Antenna::Dist()  const {return fDistance;}
double  Antenna::Dot()   const {return fDot;}
double  Antenna::Delta() const {return fDelta;}



// --------------------------------------------

Detector::Detector(){};

Detector::Detector(std::string name){
  if(name == "bistatic" || name == "Bistatic"){
    fTransmitters.push_back( Antenna(1E4, 0,0,0) ); // 1kW
    fReceivers.push_back( Antenna(0, 500, 0, 0) );
  } else if (name == "RET_CR" || name == "ret_cr"){
    fTransmitters.push_back( Antenna(1, 0,0,0) );
    fReceivers.push_back( Antenna(0, 200,0,0) );
    fReceivers.push_back( Antenna(0, 100, -100, 0) );
    fReceivers.push_back( Antenna(0, -100, -100, 0) );
    fReceivers.push_back( Antenna(0, -100,  100, 0) );
    fReceivers.push_back( Antenna(0, 200 , 200, 0) );
    fReceivers.push_back( Antenna(0, 200 ,-200, 0) );
    fReceivers.push_back( Antenna(0, -200 ,-200, 0) );
    fReceivers.push_back( Antenna(0, -200 , 200, 0) );
  } else {
    std::cout << " There is no compatible detector configurtion with "
    + name << std::endl;
  }
};

std::vector<Antenna> Detector::Transmitters(){return fTransmitters;}
std::vector<Antenna> Detector::Receivers()   {return fReceivers;}

void Detector::add_antenna(Antenna& at){
  if (at.fPower) { fTransmitters.push_back(at); }
  else { fReceivers.push_back(at); }
}

void Detector::add_antenna(Antenna at){
  if (at.fPower) { fTransmitters.push_back(at); }
  else { fReceivers.push_back(at); }
}
// -----------------------------------------------------------------------------


//  Other Detector configurations

/* Old setup locations, Krijn's code.

//
// Detector::Detector(int a, int b){
//   switch (a) {  // 5 options
//     case 1: fTransmitters.push_back(Antenna(1, (0, 0, 0), (0,0), 0, freq_obs));
//       break;  // Center
//     case 2: fTransmitters.push_back(Antenna(1, (0, 1E3, 0), (0,0), 0, freq_obs));
//       break;   // Top
//     case 3: fTransmitters.push_back(Antenna(1, (-1E3, 0, 0), (0,0), 0, freq_obs));
//       break;   // Left
//     case 4: fTransmitters.push_back(Antenna(1, (0, -1E3, 0), (0,0), 0, freq_obs));
//       break;   // Bottom
//     case 5: fTransmitters.push_back(Antenna(1, (1E3, 0, 0), (0,0), 0, freq_obs));
//       break;   // Right
//   }
//
//   switch (b) {  // Position, 5 positions
//     case 1: fReceivers.push_back(Antenna(0, (0, 5E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 2: fReceivers.push_back(Antenna(0, (0, -5E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 3: fReceivers.push_back(Antenna(0, (5E2, 0, 0), (0,0), 0, freq_obs));
//       break;
//     case 4: fReceivers.push_back(Antenna(0, (-5E2, 0, 0), (0,0), 0, freq_obs));
//       break;
//     case 5: fReceivers.push_back(Antenna(0, (0, 15E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 6: fReceivers.push_back(Antenna(0, (5E2, 1E3, 0), (0,0), 0, freq_obs));
//       break;
//     case 7: fReceivers.push_back(Antenna(0, (-5E2, 1E3, 0), (0,0), 0, freq_obs));
//       break;
//     case 8: fReceivers.push_back(Antenna(0, (-15E2, 0, 0), (0,0), 0, freq_obs));
//       break;
//     case 9: fReceivers.push_back(Antenna(0, (-1E3, 5E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 10: fReceivers.push_back(Antenna(0, (-1E3, -5E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 11: fReceivers.push_back(Antenna(0, (0, -15E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 12: fReceivers.push_back(Antenna(0, (5E2, -1E3, 0), (0,0), 0, freq_obs));
//       break;
//     case 13: fReceivers.push_back(Antenna(0, (-5E2, -1E3, 0), (0,0), 0, freq_obs));
//       break;
//     case 14: fReceivers.push_back(Antenna(0, (15E2, 0, 0), (0,0), 0, freq_obs));
//       break;
//     case 15: fReceivers.push_back(Antenna(0, (1E3, 5E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 16: fReceivers.push_back(Antenna(0, (1E3, -5E2, 0), (0,0), 0, freq_obs));
//       break;
//     case 17: fReceivers.push_back(Antenna(0, (2E3, 0, 0), (0,0), 0, freq_obs));
//       break;
//     }
//
// };


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
