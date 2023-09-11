#include "antenna.hh"

Antenna::Antenna(double xpos, double ypos, double zpos,
                 double xpol, double ypol, double zpol,
                 double power, double frequency, double gaindB):
                fPower(power),
                fFrequency(frequency),
                fWavelength(c_ice/frequency),
                fAngularFreq(2 * pi * frequency),
                fWavenumber(2 * pi * frequency / c_ice),
                fGain( pow(10., gaindB/10.) ),
                fPosition{xpos, ypos, zpos},
                fPolarization{xpol, ypol, zpol}{

// Antenna factor as sqrt(effective area)
  double factor = (c_ice/frequency)*sqrt(fGain/(4.*pi));

  // RS
  // double factor = sqrt(rx_gain*lambda/(4.*pi));

  fLeff = 1./factor;

}

void Antenna::SetDirection(std::vector<double> coords){

  // Direction to cascade: Correct orientation.TX -> CS -> RX.
  fPower != 0 ?  fDirection = direction(fPosition, coords) :
                 fDirection = direction(coords, fPosition);

  // Module of distance to cascade
  fDistance = norm(fDirection);

  // Sanity check
  assert(fDistance > c_ice/fFrequency && "Point too close to antenna");

  /* Quick fix, the far field ssumption should be much larger */
  // if(fDistance < fWavelength) {fDistance = fWavelength;}

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
}

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

// 1.643
double Antenna::GainDipole(double theta){
  return fGain * pow(sin(theta),2.6);
}

// double Antenna::GainDipole(std::vector<double> dir){
//
// };

// Accesors

double  Antenna::Power()  const {return fPower;}
double  Antenna::Gain()   const {return fGain;}
double  Antenna::E0()         const {
  return fWavenumber * sqrt(fLoad*fPower/fGain);}
double  Antenna::Leff()       const {return fLeff;}
double  Antenna::Efficiency() const {return fEff;}
double  Antenna::Load()       const {return fLoad;}

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

// TO-DO: CORRECT SHORTHAND NOTATION FOR BISTATIC AND RET_CR

Detector::Detector(std::string name){
  if(name == "bistatic" || name == "Bistatic"){
    fTransmitters.push_back( Antenna(0,0,0, 0,0,1, 10, 500*MHz, 2.15) ); // Vertical, 10 kW, 500 Hz
    fReceivers.push_back( Antenna( 0,500*m,0, 0,0,1, 0, 500*MHz, 2.15) );
   } else if (name == "T576" || name == "t576"){
    // A BEAM SETUP USES GEANT COORDINATE SYSTEM:
    // THE BEAM MOVES ALONG Z, Y IS VERTICAL. 
    fTransmitters.push_back( Antenna(-3.877*m, 0, -0.294*m,   // Position
                                      // 0, 1, 0,            // Polarisation
                                      0, 0, 1,            // Polarisation
                                      50, 2.1*GHz, 12) );   // Power, frequency, gain.
    fReceivers.push_back(    Antenna(-3.877*m, 0, 4.294*m,
                                    // 0, 1, 0,            // Polarisation
                                     0, 0, 1,
                                     0, 2.1*GHz, 18) );
    fReceivers.push_back(    Antenna(-6.395*m, 0, 3.553*m,
                                     0, 0, 1,
                                     0, 2.1*GHz, 18) );
    fReceivers.push_back(    Antenna(-6.038*m, 0, 4.536*m,
                                     0, 0, 1,
                                     0, 2.1*GHz, 18) );
  } else if (name == "RET_CR" || name == "ret_cr"){

    // TO-DO UPDATE RETCR VALUES HERE
    
    // fTransmitters.push_back( Antenna(0,0,0, 0,0,1, 50, 200*MHz, 2.15) ); // Vertical, 50 kW, 200 MHz
    // fReceivers.push_back( Antenna( 200*m,      0, 0) );
    // fReceivers.push_back( Antenna( 100*m, -100*m, 0) );
    // fReceivers.push_back( Antenna(-100*m, -100*m, 0) );
    // fReceivers.push_back( Antenna(-100*m,  100*m, 0) );
    // fReceivers.push_back( Antenna( 200*m , 200*m, 0) );
    // fReceivers.push_back( Antenna( 200*m ,-200*m, 0) );
    // fReceivers.push_back( Antenna(-200*m ,-200*m, 0) );
    // fReceivers.push_back( Antenna(-200*m , 200*m, 0) );
  } else if (name == "RNOg"){
    // Simon's GEANT RNOg model
    // GEANT defines vertical direction as y axis, so we need to swap indices.
    // IMPORTANT TO KEEP A RHS rotation frame!
    for (auto& antenna : RNOgAntennas){
      // TO-DO Correct transformation.
      // fReceivers.push_back( Antenna(antenna[1] , antenna[0], antenna[2]) );
    }
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
