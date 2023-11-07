#include "scatter.hh"

Scatter::Scatter(Antenna& tx, Antenna& rx, std::vector<ScatterPoint> points,
                  const double& sampling):
  fTX(tx), fRX(rx), fPoints(points), sampling_ratio(sampling){}

Scatter::Scatter(Antenna& tx, Antenna& rx, const double& sampling):
  Scatter(tx, rx, std::vector<ScatterPoint> {}, sampling){}


void Scatter::UpdateRX(const Antenna& new_RX){
  fRX = new_RX;
}

void Scatter::AddPoint(const ScatterPoint& p){
  fPoints.push_back(p);
  // Recalculate the number of points
  nP = fPoints.size();
}

void Scatter::AddPoints(const std::vector<ScatterPoint> new_points){
  fPoints.insert(fPoints.end(), new_points.begin(), new_points.end());
    // Recalculate the number of points
  nP = fPoints.size();
}


void Scatter::SetInDirection(std::vector<double> vertex,
                        std::vector<double> direction,
                        double length, double rand_seed){

  double dSeg = length/nP;
  
  RN_uniform rand_line(-0.5, 0.5, rand_seed);

  for(int i = 0 ; i < nP; i++){
    ScatterPoint& p = fPoints[i];
    
    if(rand_seed){ 
      p.L = (i + rand_line.get()) * dSeg;
    } else{ 
      p.L = i * dSeg; 
    }

    p.StartTime = p.L/c_vac;

    // Set position
    p.Position[0] = p.L*direction[0] + vertex[0];
    p.Position[1] = p.L*direction[1] + vertex[1];
    p.Position[2] = p.L*direction[2] + vertex[2];

  }
}



void Scatter::SetInCSPlane(std::vector<double> vertex,
                    std::vector<double> direction,
                    std::vector<double> & l_vals,
                    const std::vector<double>& r_vals){


  std::vector<double> R_TX_dir{0,0,0}, L_TX_dir{0,0,0},
                      seg_pos{0,0,0}, tx_seg{0,0,0}, rx_seg{0,0,0};
  R_TX_dir = normalize(fTX.Dir()) ;

// This only breaks for R_TX_dir == this->fDirection, exactly, we can guard against it easily
// A deviation of 10^-20 in the direction is well within numerical errors, and avoids
// making the cross product exactly zero.

  if( ( R_TX_dir == normalize( direction ) ) ) {
    R_TX_dir[0] += 1E-20;
    R_TX_dir[1] += 1E-20;
    R_TX_dir[2] += 1E-20;
  }

  L_TX_dir = cross_product( normalize( cross_product(R_TX_dir, normalize(direction) ) ), R_TX_dir);
  L_TX_dir = normalize(L_TX_dir);

  // To separate between the cascade pointing towards TX and away from TX
  double s = sgn( projection(R_TX_dir, direction) ); // cos(delta)
  if(s < 0){  std::reverse( l_vals.begin(), l_vals.end() ); }
  
  // Segment loop
  // The length of l_vals and r_vals must be nL (nP), as that is the number of
  // segments to be placed.
  ScatterPoint p;
  for(int i = 0 ; i < l_vals.size(); i++){
    p = fPoints[i];

    // Set distance from the shower head (starting point)
    p.Position[0]  = r_vals[i]*R_TX_dir[0]*s + l_vals[i]*L_TX_dir[0];
    p.Position[1]  = r_vals[i]*R_TX_dir[1]*s + l_vals[i]*L_TX_dir[1];
    p.Position[2]  = r_vals[i]*R_TX_dir[2]*s + l_vals[i]*L_TX_dir[2];
    p.L = norm(p.Position);

    p.StartTime = p.L/c_vac;

    // Set Position in lab frame (TX)
    p.Position[0] += vertex[0];
    p.Position[1] += vertex[1];
    p.Position[2] += vertex[2];
    }
}


/*  For each segment in the L_TX direction, we can find the coordinates (r_max,l_max)
  of the point where the reflectivity is the highest (it should correlate with the
  peak density for the segment, too), and place the segment at that point.
*/
void Scatter::SetInMax(std::vector<double> vertex, std::vector<double> direction, 
                          const double& dL, const double& dR,
                          const std::vector<std::vector<double>> &variable){
// i lives in L_TX, j lives in R_TX
// (The integral of variable happens along j)
    int i, j;
    double l,r;
    std::vector<double> fTXMaxLength;
    std::vector<double> fTXMaxRadius;
    double nL = variable.size();
    fTXMaxLength = std::vector<double> (nL);
    fTXMaxRadius = std::vector<double> (nL);
  for (i = 0; i < nL; i++){
    // We want the index, not the value, therefore distance()
    j = std::distance(variable[i].begin(),
            max_element(variable[i].begin(), variable[i].end() ) );
    fTXMaxLength[i] = i*dL;   // This is "b" before the rotation
    fTXMaxRadius[i] = j*dR;   // This is "a" before the rotation
  }
  
  SetInCSPlane(vertex, direction, fTXMaxLength, fTXMaxRadius);
}


/* Run time loop */
void Scatter::RunScatter(const bool save2Dmatrices){
  int steps;
  double t, t_start, t_end, freq_sampling;
  
	std::vector<double> phase_time    (nP, 0.0);
	std::vector<double> sqrt_rcs_time (nP, 0.0);
	std::vector<double> voltage_time  (nP, 0.0);

  std::vector<double> fArrivalTime = ArrivalTime();
  t_start  = *std::min_element(fArrivalTime.begin(), fArrivalTime.end()) - 5*ns ;
  t_end    = *std::max_element(fArrivalTime.begin(), fArrivalTime.end()) + 5*tau + 5*ns;

  freq_sampling = fTX.Freq()*sampling_ratio;
  steps = (t_end - t_start)*freq_sampling;

  // Memory allocation
  fDuration     = std::vector<double>(steps, 0);    // The time
  fRCS          = std::vector<double>(steps, 0);    // The RCS
  fVoltage     = std::vector<double>(steps, 0);    // The electric field
  fPower        = std::vector<double>(steps, 0);    // The power

  if (save2Dmatrices){
    fPhaseTime    = std::vector<std::vector<double>>(steps, std::vector<double> (nP, 0));
    fRCSTime      = std::vector<std::vector<double>>(steps, std::vector<double> (nP, 0));
    fVoltageTime = std::vector<std::vector<double>>(steps, std::vector<double> (nP, 0));
  }

  // Radar scatter constants
  double V0 = fTX.Wavelength()/ pow(4*pi,1.5) *
                        sqrt( 
                        fTX.Power() * fTX.Gain() *
                        fRX.Load() * fRX.Gain() ) ;

  // Time Loop!
  for (int ts = 0; ts < steps; ts++){
    t = ts/freq_sampling + t_start;
    fDuration[ts] = t;

  	std::fill(sqrt_rcs_time.begin(), sqrt_rcs_time.end(), 0.0);
		std::fill(voltage_time.begin(), voltage_time.end(), 0.0);
		std::fill(phase_time.begin(), phase_time.end(), 0.0);
    for (int i = 0; i < nP; i++){
      ScatterPoint& p = fPoints[i];

      // If active, add its contribution.
      // if(t > p.ArrivalTime){
    
      // You can also add an arbitrary cutoff (no smaller than 5*tau)
      if(t>p.ArrivalTime && t<=(p.ArrivalTime + 5*tau)){

        phase_time[i] = cos(p.Phase - fTX.AngularFreq()*t);
  
        sqrt_rcs_time[i] = sqrt(p.TCS) * // TCS =  Th * transparency * damping *N_e^2
                        pow(e,-(t-p.ArrivalTime)/tau) *  // Lifetime decay
                        phase_time[i];

        // The voltage has to also include polarization and attenuation effects. 
        voltage_time[i] =  sqrt_rcs_time[i] * 1.0/(p.RTX*p.RRX) *
                            p.PolEff * p.Attenuation;
                            // 1;
      }
    }

    // The final RCS, E field value for a given timestep is the sum of the effects of all segments.
    fRCS[ts] = std::accumulate(std::begin(sqrt_rcs_time), std::end(sqrt_rcs_time), 0.0);
    fRCS[ts] = pow(fRCS[ts],2);

    fVoltage[ts] = V0* std::accumulate(std::begin(voltage_time), std::end(voltage_time), 0.0);
    fPower[ts] = pow(fVoltage[ts],2)/fRX.Load();

    if(save2Dmatrices) {
      fPhaseTime[ts] = phase_time;
      fRCSTime[ts] = sqrt_rcs_time;
      fVoltageTime[ts] = V0*voltage_time;
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
std::vector<ScatterPoint> Scatter::Points(){ return fPoints; }

std::vector<std::vector<double>> Scatter::Position(){
  std::vector<std::vector<double>> position;
  std::transform(fPoints.begin(), fPoints.end(), position.begin(),
                  [](ScatterPoint p){return p.Position;});
  return position;
}

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




std::vector<double> Scatter::Duration(){ return fDuration; }
std::vector<double> Scatter::Voltage(){ return fVoltage; }
std::vector<double> Scatter::Power(){return fPower;}
std::vector<double> Scatter::RCS(){return fRCS;}
std::vector<std::vector<double>> Scatter::RCS_time(){ return fRCSTime; }
std::vector<std::vector<double>> Scatter::Phase_time(){ return fPhaseTime; }
std::vector<std::vector<double>> Scatter::E_time(){ return fVoltageTime; }

void Scatter::save_output_files(const std::string& output_path, const std::array<bool, 13>& flags){

  if(flags[0]){ write_1D_array(Duration(),    output_path + "_duration.txt");}
  if(flags[1]){ write_1D_array(Voltage(),    output_path + "_voltage.txt");}
  if(flags[2]){ write_1D_array(Power(),       output_path + "_power.txt");}
  if(flags[3]){ write_1D_array(TCS(),         output_path + "_TCS.txt");}
  if(flags[4]){ write_1D_array(RCS(),         output_path + "_RCS.txt");}

  if(flags[5]){ write_2D_array(Position(),    output_path + "_position.txt");}
  if(flags[6]){ write_1D_array(Phase(),       output_path + "_phase.txt");}
  if(flags[7]){ write_1D_array(ArrivalTime(), output_path + "_arrival_t.txt");}
  if(flags[8]){ write_1D_array(Attenuation(), output_path + "_attenuation.txt");}
  if(flags[9]){ write_1D_array(Polarization(),output_path + "_polarization.txt");}

  if(flags[10]){ write_2D_array(Phase_time(),  output_path + "_phase_time.txt");}
  if(flags[11]){ write_2D_array(RCS_time(),    output_path + "_RCS_time.txt");}
  if(flags[12]){ write_2D_array(E_time(),      output_path + "_E_time.txt");}
  
}