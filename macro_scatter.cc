#include "macro_scatter.hh"

Scatter::Scatter(Antenna& tx1, Antenna& rx1, Cascade& cs1):
  tx(tx1), rx(rx1), cs(cs1){
  set_direction(tx);
  set_direction(rx);

  // Attenuation model goes here
}

Antenna Scatter::transmitter(){return tx;}
Antenna Scatter::receiver(){return rx;}
Cascade Scatter::cascade(){return cs;}

/* Set the antennas directions, module and dot product with cs */
void Scatter::set_direction(Antenna& at){

  // Direction to cascade
  at.dir[0] = (cs.pos[0] - at.pos[0]);
  at.dir[1] = (cs.pos[1] - at.pos[1]);
  at.dir[2] = (cs.pos[2] - at.pos[2]);

  // Correct orientation for the Rx case (cs_to_at).
  if (at.power_ == 0) {   // Rx
    at.dir[0] = - at.dir[0];
    at.dir[1] = - at.dir[1];
    at.dir[2] = - at.dir[2];
  }

  // Module of distance to cascade
  at.dist = sqrt(pow(at.dir[0],2.0)+pow(at.dir[1],2.0)+pow(at.dir[2],2.0));

  // Sanity check
  assert(at.dist!= 0 && "Cascade overlaps antenna");
  /* THIS SHOULD THROW AND EXCEPTION SO YOU CAN CHOOSE HOW TO SOLVE IT
      FOR NOW:
 */
 if(at.dist== 0) {at.dist = 10*at.l_obs;}

  //Determine inner product between point_tc and cascade direction
  at.dot =  (cs.dir[0]*at.dir[0] + cs.dir[1]*at.dir[1] + cs.dir[2]*at.dir[2])
            /at.dist; //Normalize
}


/* Set the antennas positions and directions with IceRayTracing corrections */

/*To add, how to use IRT*/

// double * at_to_cs = IceRayTracing::IceRayTracing(x0,z0,x1,z1);

// Get the 2 rays between tx and cascade

// Get the proper distances from the solutions.
// Module of distance to cascade is now the optical path length.
// The cs_angle is the IRT recieved angle for the at_to_cs and the emitted
// angle for cs_to_rx case.

void Scatter::set_IRT_direction(Antenna& at){
  double at_r,at_z,cs_r,cs_z;
  double * at_to_cs;

  at_r = sqrt(pow(at.pos[0],2) + pow(at.pos[1],2));
  at_z = at.pos[2];
  cs_r = sqrt(pow(cs.pos[0],2)+pow(cs.pos[1],2));
  cs_z = cs.pos[2];

  if(at.power()) { // TX case
    // The antenna emits, the cascade is the one that recieves the signal.
    at_to_cs = IceRayTracing::IceRayTracing(at_r,at_z,cs_r,cs_z);

  } else { // RX case
    // The "emitter" now is the cascade and the "antenna" the receiver.
    at_to_cs = IceRayTracing::IceRayTracing(cs_r,cs_z,at_r,at_z);

  }

  // Now, let's unpack the values that we recieve from IRT.

  if(at_to_cs[6] == 0){               // No direct ray
    if(at_to_cs[4] < at_to_cs[5]){    // Reflected is shortest
      at.IRT_dist[0] = at_to_cs[4]*c_med;        // [m]
      at.IRT_dist[1] = at_to_cs[5]*c_med;   // [m]
      at.IRT_angle[0] = at_to_cs[7];      // degrees
      at.IRT_angle[1] = at_to_cs[8];      // degrees

    } else {
      at.IRT_dist[0] = at_to_cs[5]*c_med;        // [m]
      at.IRT_dist[1] = at_to_cs[4]*c_med;
      at.IRT_angle[0] = at_to_cs[8];
      at.IRT_angle[1] = at_to_cs[7];      // degrees

    }
  }

  if(at_to_cs[7] == 0){               // No reflected ray
    if(at_to_cs[3] < at_to_cs[5]){
      at.IRT_dist[0] = at_to_cs[3]*c_med;        // [m]
      at.IRT_dist[1] = at_to_cs[5]*c_med;
      at.IRT_angle[0] = at_to_cs[6];      // degrees
      at.IRT_angle[1] = at_to_cs[8];      // degrees

    } else {
      at.IRT_dist[0] = at_to_cs[5]*c_med;        // [m]
      at.IRT_dist[1] = at_to_cs[3]*c_med;
      at.IRT_angle[0] = at_to_cs[8];
      at.IRT_angle[1] = at_to_cs[6];      // degrees

    }
  }

  if(at_to_cs[8] == 0){               // No refracted ray
    if(at_to_cs[3] < at_to_cs[4]){
      at.IRT_dist[0] = at_to_cs[3]*c_med;        // [m]
      at.IRT_dist[1] = at_to_cs[4]*c_med;
      at.IRT_angle[0] = at_to_cs[6];      // degrees
      at.IRT_angle[1] = at_to_cs[7];      // degrees

    } else {
      at.IRT_dist[0] = at_to_cs[4]*c_med;        // [m]
      at.IRT_dist[1] = at_to_cs[3]*c_med;
      at.IRT_angle[0] = at_to_cs[7];      // degrees
      at.IRT_angle[1] = at_to_cs[6];      // degrees

    }
  }

  // Sanity check
  // In theory, only the first value (shortest distance) needs to be evaluated.
  assert(at.IRT_dist[0] != 0 && "Cascade overlaps antenna");
  /* THIS SHOULD THROW AND EXCEPTION SO YOU CAN CHOOSE HOW TO SOLVE IT
      FOR NOW:
 */

 // But both need to be fixed.
 if(at.IRT_dist[0] == 0) {at.IRT_dist[0] = 10*at.l_obs;}
 if(at.IRT_dist[1] == 0) {at.IRT_dist[1] = 10*at.l_obs;}

}

// ----------------------------------------------------------------------------

Scatter1D::Scatter1D(Antenna& tx, Antenna& rx, Cascade& cs):
  Scatter(tx, rx, cs){
    set_1D_values();
    get_od_cs_1D();
    run_time_loop();
  }

// Set distances, times and E field for segments.
/* Loop over the 1D segments  */
void Scatter1D::set_1D_values(){

  double l,xpos, ypos, rt, rr;
  Length.reserve(nbin), Xpos.reserve(nbin), Ypos.reserve(nbin),
  Rt.reserve(nbin), Rr.reserve(nbin), Phase.reserve(nbin),
  Arrivals.reserve(nbin), Amplitude.reserve(nbin);

  double E0 = 100;
  // TBD proper fomula for Transmitter field.

  for(int i = 0 ; i < nbin; i++){

    // Distances
    l= i*cs.get_L_bin();      // [m] distance from the shower head (starting point)
    xpos = cs.position()[0] + l*cos(cs.sph_angles()[1]);
    ypos = cs.position()[1] + l*sin(cs.sph_angles()[1]);
    rt = sqrt(pow(xpos-tx.position()[0],2)+pow(ypos-tx.position()[1],2));
    rr = sqrt(pow(xpos-rx.position()[0],2)+pow(ypos-tx.position()[1],2));

    Length.push_back(l);
    Xpos.push_back(xpos);
    Ypos.push_back(ypos);
    Rt.push_back(rt);
    Rr.push_back(rr);

    // Time evaluation
    // T0 = 0 by definition when the cascade begins (first element = head).
    //birth = l/c_vac; // time where the i'th segment starts scattering.

    // (Retarded) time where the scattered radio signal by the segment is produced.
    // production = bith - Rt/c_med;

    // (Advanced) time where the scattered signal by the segments arrives in the receiver.
    Arrivals.push_back(l/c_vac + rr/c_med);

    // E field amplitude at reciever from constant and distance dependant factors.
    // TO ADD ATTENUATION

      // Original
    // Er = E0 / (4*pi*Rt*Rr) * sqrt(rx.area);

      // Dieder Test
    Amplitude.push_back(E0 / rr);

    // Phase
      // Orignal
    //phase = omega*l/cvac + k_obs*(Rr + Rt);
      // Dieder test
    Phase.push_back(- tx.wavenr()*(rr + rt));
  }


}


/* Compute the overdense area of the cascade in slices.
The reflectance corrections are only applied to the od sections. */
void Scatter1D::get_od_cs_1D(){
  double X, r_bin, L_bin, k_mid, od_cs;
  od_cs_1D.reserve(nbin);

  std::vector<double> r_crit = cs.get_rcrit();
  std::vector<std::vector<double>> reflectance2D = cs.get_reflectance_2D();

  // Loop over regions.
  for (int i = 0; i < nbin; i++){

    // Number of layers per segment.
    int k_max = (int) r_crit[i]/r_bin;
    od_cs = 0;
    for (int k = 0; k < nbin; k++){
      // Only layers from od_region are looked, outside-in.
      k >= (nbin - k_max) ? k_mid = nbin - (k + 0.5) : k_mid = 0;
      // To evaluate for the full space as od
      // k_mid = nbin - (k + 0.5);

      od_cs += 2*k_mid*r_bin*reflectance2D[i][k]*L_bin;      // [cm^2]
      // Factor of 2 required to account for plasma tube symmetry!
    }
    od_cs_1D.push_back(od_cs);
  }
}

/* Run time loop */
void Scatter1D::run_time_loop(){

  od_cs_1D_time.reserve(nbin);
  Er_time.reserve(nbin);
  phase_time.reserve(nbin);

  vector<double> od_time_row;
  vector<double> Er_time_row;
  vector<double> phase_time_row;
  od_time_row.resize(nbin);
  Er_time_row.resize(nbin);
  phase_time_row.resize(nbin);

  double timestep= 1.0/(100*freq_obs);    // sampling frequency
  double t_start = *min_element(Arrivals.begin(), Arrivals.end()) - 5E-9;
  double t_end   = *max_element(Arrivals.begin(), Arrivals.end()) + 5*tau;
  vector<double> time;
  time.reserve((t_end - t_start)/timestep);

  double dieder_phase;

  // Loop over time.
  for (double t = t_start; t<t_end;t += timestep){
    time.push_back(t);
    // Loop over depth segments.
    for (int i = 0; i < nbin; i++){

      // Select some length values
       if(i % 2 == 0){
        // If active, add its contribution.
        if(t>Arrivals[i] && t<(Arrivals[i]+tau)){

          //phase = omega*t + k_obs*Rt;
          dieder_phase = tx.omega()*t + Phase[i];

          od_time_row.push_back(od_cs_1D[i]);
          Er_time_row.push_back(Amplitude[i]* cos(dieder_phase));
          phase_time_row.push_back(tx.omega()*(t - Arrivals[i]) + dieder_phase);
        }
        // else {
          //   Er_time_row.push_back(0);
          //   od_time_row.push_back(0);
          //   phase_time_row.push_back(0);
          // }

      }

    }

    Er_time.push_back(Er_time_row);
    od_cs_1D_time.push_back(od_time_row);
    phase_time.push_back(phase_time_row);

    Er_time_row.clear();
    od_time_row.clear();
    phase_time_row.clear();
  }

}
