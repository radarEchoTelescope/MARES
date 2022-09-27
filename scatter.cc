#include "scatter.hh"

Scatter::Scatter(Antenna& tx, Antenna& rx, std::vector<ScatterPoint> points,
                  const double& lifetime, const double& sampling):
  fTX(tx), fRX(rx), fPoints(points), tau(lifetime), sampling_ratio(sampling){}

Scatter::Scatter(Antenna& tx, Antenna& rx, const double& lifetime, const double& sampling):
  Scatter(tx, rx, std::vector<ScatterPoint> {}, lifetime, sampling){}


void Scatter::UpdateRX(const Antenna& new_RX){
  fRX = new_RX;
}

// TO-D0
    /* Updates the RX and all the RX information.
  Automatically reruns SetSegments() and RunScatter() to get updated information.
  This avoids conflicting data in the same object.
  */
  // fRX.SetDirection( fCS.Pos() );
  // SetSegments(nL);
  // RunScatter();

void Scatter::AddPoint(const ScatterPoint& p){fPoints.push_back(p);}

void Scatter::AddPoints(const std::vector<ScatterPoint> new_points){
  fPoints.insert(fPoints.end(), new_points.begin(), new_points.end());
}

// TODO: CHECK THAT POLARISATION EFF IS DEFINED AT THE ELECTRIC FIELD LEVEL

void Scatter::SetInConstIce(){ for(auto& p: fPoints) {SetInConstIce(p);} }

// Simplest model, straight lines, constant n.
void Scatter::SetInConstIce(ScatterPoint& p){
  // Set distances
  p.TXDir = direction( fTX.Pos(), p.Position) ;
  p.RTX = norm(p.TXDir);

  p.RXDir = direction(p.Position, fRX.Pos());
  p.RRX = norm(p.RXDir);

  // Set ArrivalTime
  p.ArrivalTime = p.StartTime + p.RRX/c_ice;

  // Set phase
  p.Phase = fTX.Wavenumber()*(p.RTX + p.RRX) - pi/2;

  /*
  tan^-1(f_coll/f_TX) is the phase shift in the oscillation caused by collisions.
  This would a term: 
   + atan2( f_coll, - fTX.Freq() );
  in our regime (f_Coll >>> f_TX), it's almost perfeclty -pi/2.
  */ 

  // E field attenuation at reciever from position dependant factors:
  // e^-r/Latt from medium attenuation
  p.Attenuation =  pow(e, -(p.RTX + p.RRX)/(2*att_length) ) ;

  /* Attenuation model example
  // Parametrized attenuation length for the Ross Ice Shelf, South Pole.
  double att(double freq_obs){
    double a1=469;                  // [m] Attanuation length parameter
    double a2=-0.205;               // Attanuation length parameter
    double a3=4.87E-5;              // Attanuation length parameter
    return a1+a2*freq_obs/1E6+a3*pow(freq_obs/1E6,2);
  }
  */

  // Segment's polarization direction == E field at segment.
  p.Polarization = cross_product(normalize(p.TXDir), cross_product(normalize(p.TXDir), fTX.Pol() ) );

  // Direction of electric field at the receiver.
  p.EFieldAtRX   = cross_product(normalize(p.RXDir), cross_product(normalize(p.RXDir), p.Polarization) );

  // Set polarization efficiency
  p.PolEff = abs(projection(p.EFieldAtRX, fRX.Pol()));
  //In the thin-wire theory, only the  component of  the electric- field vector parallel to the wire
  // axis can interact to form a scattered  wave. That is not our case, our layers will scatter as a free charge

  // Directivty is hardcoded as a small Herztian dipole.
  p.Directivity = 1.5*norm(p.EFieldAtRX);
   // 3/2*sin(theta_R)*sin(theta_T)
}

void Scatter::SetInBeam( const double& yb, const double& na, const double& nb){
  // If you wish to make a lookup table, it should be made here.
   for(auto& p: fPoints) {SetInBeam(p, yb, na, nb);}
 }

// TO DOUBLE CHECK, EVERYTHING HERE SHOULD BE DONE IN METERS
// TO DOUBLE CHECK, THE ANTENNAS SHOULD BE PLACED CORRECTLY ALONG THE CONVENTION
// COORDINATES
void Scatter::SetInBeam(ScatterPoint& p, const double& yb,
                        const double& na, const double& nb){
  double da, db, dt, ra, rb, ya, yt;
  double n  = nb/na;
  p.TXDir = direction(fTX.Pos(), p.Position);
  p.RXDir = direction(p.Position, fRX.Pos() );

  // Careful, because the medium change and refraction means that the relative
  // position of the antennas does not correlate with the distance that the light travels.

  /* Randomly picking coordinates. WARNING: This does not match IRT coordinates.
      x as the direction of the beam.
      y as the perpendicular to the beam in the plane where RX lives.
      z as the vertical.

      (P) - - - - - - - - (xt) - - - - - > beam (x direction = d direction)
      |                    |
      |            (yb)    |
      |------------------------------- (interface)
      |            (ya)    |
      |                    |
      (yt) ------------ (RX)
      |
      |
      (y direction)

  yb is the distance between the beam line and the interface.
  ya is the distance between the interface and the receiver.
  yt is the total distance between the point and the antenna
    in the "vertical" direction of reflection.
  dt is the total distance between P and RX out of the y direction.
    (dt = xt in 2D, dt incorporates x and z in 3D).
  We need to figure out the da, db values.
  */

  dt = sqrt( pow(p.TXDir[0],2) + pow(p.TXDir[2],2) );
  yt = p.TXDir[1];
  ya = yt - yb;

  db = estimate_interface_x( ya, yb, dt, n);
  da = dt - db;

  // ra = sqrt(pow(da,2) + pow(ya,2));
  // rb = sqrt(pow(db,2) + pow(yb,2));

  // Total travelled distance by the TX ray of light
  p.RTX = sqrt(pow(da,2) + pow(ya,2)) + sqrt(pow(db,2) + pow(yb,2));

  dt = sqrt( pow(p.RXDir[0],2) + pow(p.RXDir[2],2) );
  yt = p.RXDir[1];
  ya = yt - yb;

  db = estimate_interface_x( ya, yb, dt, n);
  da = dt - db;

  // Total travelled distance by the ray of light
  p.RRX = sqrt(pow(da,2) + pow(ya,2)) + sqrt(pow(db,2) + pow(yb,2));

  // TODO-TRIPLE CHECK THIS

  // Set ArrivalTime
  p.ArrivalTime = p.StartTime + p.RRX/c_ice;

  // Set phase
  p.Phase = fTX.Wavenumber()*(p.RTX + p.RRX)
                + atan2( f_coll, - fTX.Freq() );
  // tan^-1(f_coll/f_TX) is the phase shift in the oscillation caused by collisions.
  // in our regime, it could be almost fixed to -Pi/2.

  // DOUBLE CHECK MINUS SIGN

  /* E field attenuation at reciever from long-distance attenuation is neglected in 
  a beam setup */
  p.Attenuation = 1;
  // p.Attenuation =  pow(e, -(p.RTX + p.RRX)/(2*att_length) ) ;


  // Segment's polarization direction == E field at segment.
  p.Polarization = cross_product(normalize(p.TXDir), cross_product(normalize(p.TXDir), fTX.Pol() ) );

  // Direction of electric field at the receiver.
  p.EFieldAtRX   = cross_product(normalize(p.RXDir), cross_product(normalize(p.RXDir), p.Polarization) );

  // Set polarization efficiency
  p.PolEff = abs(projection(p.EFieldAtRX, fRX.Pol()));
  //In the thin-wire theory, only the  component of  the electric- field vector parallel to the wire
  // axis can interact to form a scattered  wave. That is not our case, our layers will scatter as a free charge

  // Directivty is hardcoded as a small Herztian dipole.
  p.Directivity = 1.5*norm(p.EFieldAtRX);
   // 3/2*sin(theta_R)*sin(theta_T)
}

double Scatter::estimate_interface_x( const double& ya, const double& yb,
                             const double& dt, const double& n,
                             const double& tolerance ){
  int steps = 0, total_tries = 1E3;
  double dt_est, x_est, tmp, x_start = 0, x_end = dt;
  while (steps < total_tries){
    x_est = (x_start + x_end)/2.;
    // rb = sqrt(pow(x_est,2)+pow(yb,2));
    tmp = n/sqrt(pow(x_est,2)+pow(yb,2));
    dt_est = x_est + x_est*ya*tmp  / sqrt(1 - pow(x_est*tmp,2));
    if (abs(dt-dt_est) < tolerance){return x_est;}
    steps++;
    dt > dt_est ? x_start = x_est : x_end = x_est;
  }
  return nan("");
}

void Scatter::SetWithIRT(){ for(auto& p: fPoints) {SetWithIRT(p);} }

void Scatter::SetWithIRT(ScatterPoint& p) {
  /* code call to IRT */

  // TODO TRIPLE CHECK THAT THIS IS STILL VALID
  // Set phase
  p.Phase = fTX.Wavenumber()*(p.RTX + p.RRX)
                + atan2( f_coll, - fTX.Freq() );
  // tan^-1(f_coll/f_TX) is the phase shift in the oscillation caused by collisions.
  // in our regime, it could be almost fixed to -Pi/2.

  // DOUBLE CHECK MINUS SIGN

// DOES IRT INCLUDE ICE ATTENUATION?
  // E field attenuation at reciever from position dependant factors:
  // e^-r/Latt from medium attenuation
  p.Attenuation =  pow(e, -(p.RTX + p.RRX)/(2*att_length) ) ;

  /* Attenuation model example
  // Parametrized attenuation length for the Ross Ice Shelf, South Pole.
  double att(double freq_obs){
    double a1=469;                  // [m] Attanuation length parameter
    double a2=-0.205;               // Attanuation length parameter
    double a3=4.87E-5;              // Attanuation length parameter
    return a1+a2*freq_obs/1E6+a3*pow(freq_obs/1E6,2);
  }
  */

  // Segment's polarization direction == E field at segment.
  p.Polarization = cross_product(normalize(p.TXDir), cross_product(normalize(p.TXDir), fTX.Pol() ) );

  // Direction of electric field at the receiver.
  p.EFieldAtRX   = cross_product(normalize(p.RXDir), cross_product(normalize(p.RXDir), p.Polarization) );

  // Set polarization efficiency
  p.PolEff = abs(projection(p.EFieldAtRX, fRX.Pol()));
  //In the thin-wire theory, only the  component of  the electric- field vector parallel to the wire
  // axis can interact to form a scattered  wave. That is not our case, our layers will scatter as a free charge

  // Directivty is hardcoded as a small Herztian dipole.
  p.Directivity = 1.5*norm(p.EFieldAtRX);
   // 3/2*sin(theta_R)*sin(theta_T)
}

/* Run time loop */
// Requires set_segments and set_radar_cs();
void Scatter::RunScatter(const bool save2Dmatrices){
  // Temporary variables
  int steps;
  double t, t_start, t_end, freq_sampling;
  ScatterPoint p;
	std::vector<double> phase_time    (nP, 0.0);
	std::vector<double> sqrt_rcs_time (nP, 0.0);
	std::vector<double> voltage_time  (nP, 0.0);

  std::vector<double> fArrivalTime = ArrivalTime();
  t_start  = *std::min_element(fArrivalTime.begin(), fArrivalTime.end()) - 5E-9 ;
  t_end    = *std::max_element(fArrivalTime.begin(), fArrivalTime.end()) + tau +5E-9;

  freq_sampling = fTX.Freq()*sampling_ratio;
  steps = (t_end - t_start)*sampling_ratio;

  // Memory allocation
  fDuration     = std::vector<double>(steps, 0);    // The time
  fRCS          = std::vector<double>(steps, 0);    // The RCS
  fWaveform     = std::vector<double>(steps, 0);    // The electric field
  fPower        = std::vector<double>(steps, 0);    // The power

  if (save2Dmatrices){
    fPhaseTime    = std::vector<std::vector<double>>(steps, std::vector<double> (nP, 0));
    fRCSTime      = std::vector<std::vector<double>>(steps, std::vector<double> (nP, 0));
    fWaveformTime = std::vector<std::vector<double>>(steps, std::vector<double> (nP, 0));
  }

  // Radar scatter constants
  double V0 = fTX.Wavelength()/ pow(4*pi,1.5) *
                        sqrt( fTX.Power() * fTX.Gain() *
                        fRX.Load() * fRX.Gain() ) ;

  // Time Loop!
  for (int ts = 0; ts < steps; ts++){
    t = ts/freq_sampling + t_start;
    fDuration[ts] = t;

  	std::fill(sqrt_rcs_time.begin(), sqrt_rcs_time.end(), 0.0);
		std::fill(voltage_time.begin(), voltage_time.end(), 0.0);
		std::fill(phase_time.begin(), phase_time.end(), 0.0);

    for (int i = 0; i < nP; i++){
      p = fPoints[i];

      // If active, add its contribution.
      if(t>p.ArrivalTime){

      // You can also add an arbitrary cutoff (no smaller than 5*tau)
      // if(t>p.ArrivalTime && t<=(p.fArrivalTime + 5*tau)){

        phase_time[i] = cos(p.Phase - fTX.AngularFreq()*t);

        sqrt_rcs_time[i] = p.TCS *                          // Th * transparency * damping *N_e^2
                        pow(e,-(t-p.ArrivalTime)/tau) *  // Lifetime decay
                        phase_time[i];

        // Waveform == Voltage now
        // The voltage has to also include polarization and attenuation effects. 
        voltage_time[i] =  sqrt_rcs_time[i] * 1.0/(p.RTX*p.RRX) *
                            p.PolEff * p.Attenuation;
      }
    }

    // The final RCS, E field value for a given timestep is the sum of the effects of all segments.
    fRCS[ts] = std::accumulate(std::begin(sqrt_rcs_time), std::end(sqrt_rcs_time), 0.0);
    fRCS[ts] = pow(fRCS[ts],2);

    fWaveform[ts] = V0* std::accumulate(std::begin(voltage_time), std::end(voltage_time), 0.0);
    fPower[ts] = pow(fWaveform[ts],2)/fRX.Load();

    if(save2Dmatrices) {
      fPhaseTime[ts] = phase_time;
      fRCSTime[ts] = sqrt_rcs_time;
      fWaveformTime[ts] = voltage_time;
    }
  }
}

// void Scatter::SetAtDirection(Antenna &at){ at.SetDirection( cs.Pos() ); }

// /* Set the antennas directions, module and dot product with cs */
// void Scatter::SetAtDirCenter(Antenna &at){
//   // Find the cascade center
//   std::vector<double> center_pos = { cs.Pos()[0] + cs.fLtot/2*cs.Dir()[0],
//                                      cs.Pos()[1] + cs.fLtot/2*cs.Dir()[1],
//                                      cs.Pos()[2] + cs.fLtot/2*cs.Dir()[2]
//                                    };
//   at.SetDirection( center_pos );
// }

/*  --------------------------------------------------------------------------*/


// Accesors
Antenna Scatter::TX(){return fTX;}
Antenna Scatter::RX(){return fRX;}
// TODO ADD ACCEESOR FOR POINTS!

// TO-DO ADD ACCESOR FOR COORDINATES

std::vector<double> Scatter::Phase(){
  std::vector<double> phase (fPoints.size(), 0);
  std::transform(fPoints.begin(), fPoints.end(), phase.begin(),
                  [](ScatterPoint p){return p.Phase;});
  return phase;
}

std::vector<double> Scatter::ArrivalTime(){
  std::vector<double> arrivals (fPoints.size(), 0);
  std::transform(fPoints.begin(), fPoints.end(), arrivals.begin(),
                  [](ScatterPoint p){return p.ArrivalTime;});
 return arrivals;
 }

std::vector<double> Scatter::Attenuation(){
  std::vector<double> attenuation (fPoints.size(), 0);
  std::transform(fPoints.begin(), fPoints.end(), attenuation.begin(),
                  [](ScatterPoint p){return p.Attenuation;});
  return attenuation;
 }

std::vector<double> Scatter::Directivity(){
  std::vector<double> directivity (fPoints.size(), 0);
  std::transform(fPoints.begin(), fPoints.end(), directivity.begin(),
                  [](ScatterPoint p){return p.Directivity;});
  return directivity;
 }

std::vector<double> Scatter::Polarization(){
  std::vector<double> polarization (fPoints.size(), 0);
  std::transform(fPoints.begin(), fPoints.end(), polarization.begin(),
                  [](ScatterPoint p){return p.PolEff;});
  return polarization;
}

std::vector<double> Scatter::TCS(){
  std::vector<double> tcs (fPoints.size(), 0);
  std::transform(fPoints.begin(), fPoints.end(), tcs.begin(),
                  [](ScatterPoint p){return p.TCS;});
  return tcs;
 }

// std::vector<std::vector<double>> Scatter::Coordinates(){return fSegmentCoord;}


std::vector<double> Scatter::Duration(){ return fDuration; }
std::vector<double> Scatter::Waveform(){ return fWaveform; }
std::vector<double> Scatter::Power(){return fPower;}
std::vector<double> Scatter::RCS(){return fRCS;}
std::vector<std::vector<double>> Scatter::RCS_time(){ return fRCSTime; }
std::vector<std::vector<double>> Scatter::Phase_time(){ return fPhaseTime; }
std::vector<std::vector<double>> Scatter::E_time(){ return fWaveformTime; }
