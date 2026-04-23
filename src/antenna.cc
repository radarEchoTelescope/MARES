#include "antenna.hh"
Antenna::Antenna(double xpos, double ypos, double zpos,
                 double xpol, double ypol, double zpol,
                 double power, double frequency, double gaindB, int modMode, double modDuration, double modBandwidth):
                fPower(power),
                fFrequency(frequency),
                fWavelength(c_ice/frequency),
                fAngularFreq(2 * pi * frequency),
                fWavenumber(2 * pi * frequency / c_ice),
                fGain( pow(10., gaindB/10.) ),
                fPosition{xpos, ypos, zpos},
                fPolarization{xpol, ypol, zpol},
                fMode(modMode),
                fModDuration(modDuration),
                fModBandwidth(modBandwidth){

  // Hardcode normalised polarisation factors 
  fPolarization = normalize(fPolarization);                
  // Antenna factor as sqrt(effective area)
  double factor = (c_ice/frequency)*sqrt(fGain/(4.*pi));

  // RS
  // double factor = sqrt(rx_gain*lambda/(4.*pi));
  fLeff = 1./factor;

  assert(modDuration != 0.0 && "The total duration of the modulation should be non-zero. If running in the CW mode, the mode should be set to -1 and this parameter to a random non-zero value.");
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

double Antenna::ModBandwidth() const {return fModBandwidth;} 
double Antenna::ModDuration() const {return fModDuration;}
std::string Antenna::Mode() const{
  if (fMode== -1){
    return "CW mode";
  }
  if (fMode== 0){
    return "FMCW mode : Triangular modulation";
  }
  if (fMode== 1){
    return "FMCW mode : Sawtooth modulation";
  } 
  else
  {
    return "Frequency modulation mode was set incorrectly";
  }   
}

double Antenna::PhaseTriangularMod(double time) const {
// The time here is general. This means it can correspond with the retarded TX time, or the observation (RX) time ,... 
  double slopeMod= 0.5* fModBandwidth/(fModDuration/2);
  // std::cout<<"SlopeMod:"<<slopeMod<<std::endl;
  double tmp_phase;
  double eval_time;

  // Note: the eval_time is calculated to make the modulation periodic. It is defined as 0 < eval_time < fModDuration. 
  // The two lines of code below work both in case of negative or positive time values. 
  int multi = floor(time/fModDuration);
  eval_time=time-multi*fModDuration; 
  if(eval_time<fModDuration/2)
  {
  tmp_phase= fFrequency*eval_time + (slopeMod * pow(eval_time,2));
  }
  else {
  tmp_phase= fFrequency*(abs(eval_time-fModDuration)) + (slopeMod *pow(abs(eval_time-fModDuration),2)); // this implementation makes sure the phase oscillates correctly for a triangular modulation on the downward slope
  }
  return 2*pi*tmp_phase; 
}

double Antenna::PhaseSawtoothMod(double time) const {
// The time here is general. This means it can correspond with the retarded TX time, or the observation (RX) time ,... 
  double slopeMod= 0.5*fModBandwidth/fModDuration;
  double tmp_phase;
  double eval_time;

  // Note: the eval_time is calculated to make the modulation periodic. It is defined as 0 < eval_time < fModDuration. 
  // The two lines of code below work both in case of negative or positive time values. 
  int multi = floor(time/fModDuration);
  eval_time=time-multi*fModDuration; 

  tmp_phase= fFrequency*eval_time + (slopeMod * pow(eval_time,2));

  return 2*pi*tmp_phase; 
}

double Antenna::PhaseFMCW(double time) const
 // this function returns the frequency at a specific time given a modulation. 0: Triangular, 1: Sawtooth
{
  double freq_t;
  if(fMode==0){
    freq_t=PhaseTriangularMod(time);
  }
  else if(fMode==1){
    freq_t=PhaseSawtoothMod(time);
  }
  else{
    return std::numeric_limits<double>::quiet_NaN();
  }
  return freq_t ;
}

double Antenna::TriangularMod(double time) const {
// The time here is general. This means it can correspond with the retarded TX time, or the observation (RX) time ,... 
  double slopeMod= fModBandwidth/(fModDuration/2);
  // std::cout<<"SlopeMod:"<<slopeMod<<std::endl;
  double tmp_freq;
  double eval_time;

  // Note: the eval_time is calculated to make the modulation periodic. It is defined as 0 < eval_time < fModDuration. 
  // The two lines of code below work both in case of negative or positive time values. 
  int multi = floor(time/fModDuration);
  eval_time=time-multi*fModDuration; 
  if(eval_time<fModDuration/2)
  {
  tmp_freq= fFrequency + (slopeMod * eval_time);
  }
  else {
  tmp_freq= (fFrequency + fModBandwidth) - slopeMod * (eval_time-fModDuration/2);
  }
  return tmp_freq; 
}

double Antenna::SawtoothMod(double time) const {
// The time here is general. This means it can correspond with the retarded TX time, or the observation (RX) time ,... 
  double slopeMod= fModBandwidth/fModDuration;
  double tmp_freq;
  double eval_time;

  // Note: the eval_time is calculated to make the modulation periodic. It is defined as 0 < eval_time < fModDuration. 
  // The two lines of code below work both in case of negative or positive time values. 
  int multi = floor(time/fModDuration);
  eval_time=time-multi*fModDuration; 

  tmp_freq= fFrequency + (slopeMod * eval_time);

  return tmp_freq; 
}

double Antenna::Freq(double time) const
 // this function returns the frequency at a specific time given a modulation. 0: Triangular, 1: Sawtooth
{
  double freq_t;
  if(fMode==0){
    freq_t=TriangularMod(time);
  }
  else if(fMode==1){
    freq_t=SawtoothMod(time);
  }
  else{
    return std::numeric_limits<double>::quiet_NaN();
  }
  return freq_t ;
}


double Antenna::AngularFreq(double time) const // angular frequency for a triangular modulation. 
{
  return 2* pi * Freq(time); 
}

double Antenna::Wavelength(double time) const // wavelength for a triangular modulation. 
{
  return c_ice/Freq(time); 
}

double Antenna::Wavenumber(double time) const // wavenumber for a triangular modulation. 
{
  return 2* pi* Freq(time)/c_ice; 
}


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
    fTransmitters.push_back( Antenna(0,0,-100*m, 0,1,0, 1E3, 50*MHz, 1) ); // Vertical, 10 W, 500 Hz
    fReceivers.push_back( Antenna( -250*m, 0, 0, 0,1,0, 0, 50*MHz, 1) );
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
  } else if (name == "stargrid" || name == "StarGrid"){
    // Simon's Grid model used in CORSIKA+GEANT simulations. 
    // GEANT defines vertical direction as y axis, so we need to swap indices.
    // IMPORTANT TO KEEP A RHS rotation frame!
    for (auto& antenna : StarGrid){
      fReceivers.push_back( Antenna(
        antenna[0] , antenna[2], antenna[1],
        0, 0, 1, 
        0, 100*MHz, 2.15));
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

void load_antenna_list(const libconfig::Setting& at_list, Detector& dect){
  // Antenna placeholder variables.
  double xpos, ypos, zpos, xpol, ypol, zpol, power, freq, gaindB,modBandwidth,modDuration;
  int modMode;
   for (int i = 0; i < at_list.getLength(); i++){
    power = NAN;
    // Grab the next antenna
    const libconfig::Setting &at = at_list[i];

    // Try to load the values from the antenna
    if( !(
          at["position"].lookupValue("x", xpos)     &&
          at["position"].lookupValue("y", ypos)     &&
          at["position"].lookupValue("z", zpos)     &&
          at["polarization"].lookupValue("x", xpol) &&
          at["polarization"].lookupValue("y", ypol) &&
          at["polarization"].lookupValue("z", zpol) &&
          at.lookupValue("power", power)            &&
          at.lookupValue("frequency", freq)         &&
          at.lookupValue("gaindB", gaindB)          &&
          at.lookupValue("modBandwidth",modBandwidth)&&
          at.lookupValue("modDuration",modDuration) &&
          at.lookupValue("modMode",modMode)

        )
    ){
      if(power != 0){
        std::cerr << " Transmitter " << i << " was skipped" << std::endl;
      } else if (power == 0){
        std::cerr << " Receiver " << i << " was skipped" << std::endl;
      } else{
        std::cerr << " Unknown antenna " << i << " was skipped" << std::endl;
      } 
      continue;
    }
    dect.add_antenna( Antenna(
              xpos *m, ypos *m, zpos *m,
              xpol, ypol, zpol,
              power *W, freq *GHz, gaindB, modMode, modDuration*ns, modBandwidth *GHz
      )
    );
  }
}


Detector load_detector_config(libconfig::Config& dect_config){
  Detector dect;
  libconfig::Setting &root = dect_config.getRoot();

  try{
    const libconfig::Setting& tx_list = root["detector"]["transmitter"]; 
  } catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Your config file is missing transmitter group" << std::endl;
    exit(EXIT_FAILURE);
  }
  const libconfig::Setting& tx_list = root["detector"]["transmitter"];
  load_antenna_list(tx_list, dect);

  try{
    const libconfig::Setting& rx_list = root["detector"]["receiver"]; 
  } catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Your config file is missing receiver group" << std::endl;
    exit(EXIT_FAILURE);
  }
  const libconfig::Setting& rx_list = root["detector"]["receiver"]; 
  load_antenna_list(rx_list, dect);

  return dect;


}

Detector load_detector_file(const std::string& dect_config_filepath){
  libconfig::Config cfg;
  load_config_file(cfg, dect_config_filepath.c_str());
  Detector dect = load_detector_config(cfg);
  return dect;
}

// -----------------------------------------------------------------------------
