#include "macro_scatter.hh"

Scatter::Scatter(Antenna& tx1, Antenna& rx1, Cascade& cs1):
  tx(tx1), rx(rx1), cs(cs1){
    // set_direction_center(tx);
  set_direction(tx);
  set_direction(rx);
  set_segments();

  // Attenuation model goes here

  /*
  // Parametrized attenuation length for the Ross Ice Shelf, South Pole.
  double att(double freq_obs){
  double a1=469;                  // [m] Attanuation length parameter
  double a2=-0.205;               // Attanuation length parameter
  double a3=4.87E-5;              // Attanuation length parameter
  return a1+a2*freq_obs/1E6+a3*pow(freq_obs/1E6,2);
  }
  */


}

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

  cs.pos[1] > 0 ? at.ang = acos(at.dot) : at.ang = - acos(at.dot);
}

/* Set the antennas directions, module and dot product with cs */
void Scatter::set_direction_center(Antenna& at){

  // 2D FOR NOW

  // Find the cascade center
  double cs_xpos = cs.pos[0] + cs.L_tot/2*cos(cs.sph_ang[1]);
  double cs_ypos = cs.pos[1] + cs.L_tot/2*sin(cs.sph_ang[1]);

  // Direction to cascade
  at.dir[0] = (cs_xpos - at.pos[0]);
  at.dir[1] = (cs_ypos - at.pos[1]);
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

  // std::cout << cs.sph_ang[0] << '\t' << cs.sph_ang[1] << '\n'
  //           << cs.dir[0] << '\t' << cs.dir[1] << '\t'<< cs.dir[2] << '\n'
  //           << at.dir[0] << '\t' << at.dir[1] << '\t'<< at.dir[2] << '\n'
  //           << at.dot << std::endl;

  cs.pos[1] > 0 ? at.ang = acos(at.dot) : at.ang = - acos(at.dot);

}

// Set distances, times and E field for segments.
/* Loop over the 1D segments  */
void Scatter::set_segments(){
  // RN_uniform rand_line(-0.5, 0.5, 42);          // Same seed for debugging.
  RN_uniform rand_line(-1, 1, 42);          // Same seed for debugging.
  // RN_uniform rand_line(-1, 1, time(0));  // Different seed for random, independent runs.

  double l,xpos, ypos, rt, rr, rt0;
  _phase.reserve(nbin), _arrival.reserve(nbin), _amplitude.reserve(nbin);
  _segment_coords = std::vector<std::vector<double>>(nbin, vector<double> (5, 0));

  // TO DO proper fomula for Transmitter field.
  double E0 = 100;

  for(int i = 0 ; i < nbin; i++){

    // Distances
    l    = (i + rand_line.get()) * cs.get_L_bin();      // [m] distance from the shower head (starting point)
    // l    = i*cs.get_L_bin();      // [m] distance from the shower head (starting point)
    xpos = cs.position()[0] + l*cos(cs.sph_angles()[1]);
    ypos = cs.position()[1] + l*sin(cs.sph_angles()[1]);
    rt   = distance(xpos, ypos, tx.position()[0], tx.position()[1]);
    rr   = distance(xpos, ypos, rx.position()[0], rx.position()[1]);

    if(i == 0){rt0 = rt;}

    _segment_coords[i][0] = l;
    _segment_coords[i][1] = xpos;
    _segment_coords[i][2] = ypos;
    _segment_coords[i][3] = rt;
    _segment_coords[i][4] = rr;

    // Time evaluation
    // T0 = 0 by definition when the cascade begins (first element = head).
    //birth = l/c_vac; // time where the i'th segment starts scattering.

    // (Retarded) time where the scattered radio signal by the segment is produced.
    // production = birth - Rt/c_med;

    // (Advanced) time where the scattered signal by the segments arrives in the receiver.
    _arrival.push_back(l/c_vac + rr/c_med);

    // E field amplitude at reciever from constant and distance dependant factors.

    // Dieder's Amplitude

    _amplitude.push_back(E0 / (rt * rr) ); // Correct
    // _amplitude.push_back(E0 / (rr * rr) ); // Dieder, with mistake

      // Original
    // Er = E0 / (4*pi*Rt*Rr) * sqrt(rx.area);

    /* Dieder's phase
    phase = k*(rt + rr) - omega(t + l/c_vac + (rt(0) - rt)/c_med - rr/c_med)

    where:

      k*(rt + rr) = tx.wavenr()*(rt + rr); spatial phase
      l/c_vac; cascade propagation term
      (rt(0) - rt)/c_med; change in phase for each point w.r.t starting point (i = 0).
      rt(0) = | cs_pos - tx.pos | = cs.distance if computed from the head.
      rr/c_med; retardation effects ( = time@receiver - t@cascade)

      This becomes ( with omega/c_med = k):
      // _phase.push_back( tx.wavenr()*(2*(rt + rr) - rt0 - l/refindex) );
    */

    // _phase.push_back( tx.wavenr()*(2*(rt + rr) - rt0 - l/refindex) ); // old
    // _phase.push_back( tx.wavenr()*(rt + rr - rt0 - l/refindex) ); // Simplified, supposedly
      _phase.push_back( tx.wavenr()*(2*rt + rr - rt0) ); // Simplified from Dieder code
  }
}

/* Run time loop */
void Scatter::run_time_loop(){
// Requires set_segments and set_radar_cs();

  double t;
  double t_start  = *min_element(_arrival.begin(), _arrival.end()) - 5E-9;
  double t_end    = *max_element(_arrival.begin(), _arrival.end()) + 5*tau;
  // std::cout << t_start << '\t' << t_end << std::endl;

  double sampling = (100*tx.freq());
  int steps = (t_end - t_start)*sampling;

  _time       = std::vector<double>(steps, 0);    // The time
  _waveform   = std::vector<double>(steps, 0);    // The electric field

  _er_time    = std::vector<std::vector<double>>(steps, vector<double> (nbin, 0));
  _phase_time = std::vector<std::vector<double>>(steps, vector<double> (nbin, 0));
  _rcs_time   = std::vector<std::vector<double>>(steps, vector<double> (nbin, 0));

  // Loop over time.
  for (int j = 0; j < steps; j ++){
    t = j/sampling + t_start;
    _time[j] = t;
    for (int i = 0; i < nbin; i++){

      // Select some length values?
       // if(i == 0){

        // If active, add its contribution.
        if(t>_arrival[i] && t<(_arrival[i]+tau)){

          _rcs_time[j][i] = _rcs[i];
          _phase_time[j][i] = _phase[i] - tx.omega()*t;
          //_phase_time[j][i] = _phase[i] - tx.omega()*(t-t_start);

          _er_time[j][i] = _amplitude[i]* _rcs[i] * cos(_phase[i] - tx.omega()*t);
        }

      // }  // Selection closing

    }

    // The final E field value is the sum of all particles for the timestep.
    _waveform[j] = std::accumulate(std::begin(_er_time[j]), std::end(_er_time[j]), 0.0);
  }
}

// Accesors
Antenna Scatter::transmitter(){return tx;}
Antenna Scatter::receiver(){return rx;}
Cascade Scatter::cascade(){return cs;}
std::vector<double> Scatter::radar_cs(){ return _rcs; }

  // Not set before set_segments()
std::vector<std::vector<double>> Scatter::segement_coords(){return _segment_coords;}
std::vector<double> Scatter::amplitude(){ return _amplitude; }
std::vector<double> Scatter::arrivals(){ return _arrival ; }
std::vector<double> Scatter::phase(){return _phase;}


  // Not set before run_time_loop()
std::vector<double> Scatter::time(){ return _time ; }
std::vector<double> Scatter::waveform(){ return _waveform ; }
std::vector<std::vector<double>> Scatter::rcs_time(){ return _rcs_time; }
std::vector<std::vector<double>> Scatter::phase_time(){ return _phase_time; }
std::vector<std::vector<double>> Scatter::wave_time(){ return _er_time; }

// // --------------------------------------------------------------------------

//==============================================================================
// Line1D is Dieder's model

Line1D::Line1D(Antenna& tx, Antenna& rx, Cascade& cs): Scatter(tx, rx, cs){
  _rcs = std::vector<double>(nbin, 1); // The line's segments have rcs untiy.
  run_time_loop();
  }

// -----------------------------------------------------------------------------
// THE ROTATED VERSION

Cascade1D::Cascade1D(Antenna& tx, Antenna& rx, Cascade& cs):
  Scatter(tx, rx, cs){
    set_rotated_density();

    // write_2D_array(get_density_cs(), "_density_cs.txt", 1);
    // write_2D_array(get_density_tx(), "_density_tx.txt", 1);

    set_reflectivity();
    set_radar_cs();
    run_time_loop();
  }


/* Make 2D-array of density profile */
void Cascade1D::set_rotated_density(){
	double l, r, x, y;

  coords = std::vector<std::vector<double>> (nbin, vector<double> (4, 0));
  density_cs = std::vector<std::vector<double>> (nbin, vector<double> (nbin, 0));
  density_tx = std::vector<std::vector<double>> (nbin, vector<double> (nbin, 0));

  /* The rotation we have to perform over the angle is: */
  // Alpha is the rotation of the r lines over the centers in the vertical or horizontal;
  alpha = pi/2 + tx.angle() - cs.sph_angles()[1];
  double L_2 = cs.get_L_tot()/2.0 * 100.0; // [cm]

//  std::cout
//            << cs.get_L_tot()*100/2.0 << '\t'  <<cs.get_r_tot() << '\n'
//            << tx.angle() << '\t' << cs.sph_angles()[1] << '\t' << alpha << '\n'
//            << cos(alpha) <<'\t' << sin(alpha) << '\n'
//            << abs(cs.get_L_tot()*cos(alpha) *100.0) << '\t'
//            << abs(cs.get_r_tot()*sin(alpha)) << '\n'
// // << max(cs.get_L_tot()/2.0 * 100.0 * cos(alpha), cs.get_r_tot()*cos(alpha)) << '\t'
// // <<  max(cs.get_L_tot()/2.0 * 100.0 * sin(alpha), cs.get_r_tot()*sin(alpha)) << '\n'
// <<  std::endl;

if(abs(cs.get_L_tot()*cos(alpha) *100.0) > abs(cs.get_r_tot()*sin(alpha))){
  // EDGE-ON MODE: Centers along the l dimension.

  for (int i = 0; i < nbin; i++){
    l = (2.0*i/nbin - 1) * L_2 ;                // Length from -L/2 to L/2.
    for (int j = 0; j < nbin; j++){
      r = (1 - 2.0*j/nbin) * cs.get_r_tot();    // Radius between +r and -r.

        // The minus sign is to correct for the loop order
        // which is fixed because we need to loop "outside-in".
        x = l - r*sin(alpha);
        y = r*sgn(cos(alpha));

        coords[i][0] = l;
        coords[i][1] = r;
        coords[i][2] = x;
        coords[i][3] = y;
        density_cs[j][i] = cs.dens((l + L_2)*rho_ice,r); // Cascade frame (non-rotated) "green"
        density_tx[j][i] = cs.dens((x + L_2)*rho_ice,y);
        // density_cs[i][j] = cs.dens((l + L_2)*rho_ice,r); // Cascade frame (non-rotated) "green"
        // density_tx[i][j] = cs.dens((x + L_2)*rho_ice,y);
      }
    }
  } else {
  // FACE-ON: // Centers along the r dimension.
    std::cout << "WARNING: Face-on" << std::endl;
    for (int j = 0; j < nbin; j++){
      r = (1 - 2.0*j/nbin) * cs.get_r_tot();    // Radius between +r and -r.
      for (int i = 0; i < nbin; i++){
        l = (2.0*i/nbin - 1) * L_2 ;                // Length from -L/2 to L/2.

        y = r - l*cos(alpha);  //
        x = l;//*sin(alpha);      //
        // Due to the limted range, there is not need for sgn() function.

        coords[i][0] = l;
        coords[i][1] = r;
        coords[i][2] = x;
        coords[i][3] = y;
        density_cs[i][j] = cs.dens((l + L_2)*rho_ice,r); // Cascade frame (non-rotated) "green"
        density_tx[i][j] = cs.dens((x + L_2)*rho_ice,y);
      }
    }
    // NO NEED FOR TRANSPOSING, THE SWAPPED INDICES ABOVE SHOULLD DO THE TRICK.
    // density_tx = transpose(density_tx);
    // density_cs = transpose(density_cs);
  }

}


double Cascade1D::absorption(double& dens){
  double f_plasma = cs.fplasma(dens);
  double w, a, b, q;

  // Exact solution from dispersion relation with collisions.
  if (f_coll != 0){
    w = pow(f_plasma,2)/( pow(tx.freq(),2) + pow(f_coll,2) );
    a = 1 - w;
    b = (f_coll/ tx.freq()) * w;

    // double p = (w_obs/c_med)*np.sqrt((np.sqrt(a**2 + b**2) + a) /2 )
    q = (tx.freq()/c_med)*sqrt((sqrt(pow(a,2) + pow(b,2)) - a) /2 );
  } else {
  // Collisionless model
    f_plasma > tx.freq() ? q = cmed_cm/(f_plasma) : q = 0;
  }
  return q;
}

double Cascade1D::skin_depth(double& dens){
  return (1/absorption(dens));
}

// TODO FIx the orientation mode.
void Cascade1D::set_reflectance(){
	double reflectance, reflectivity;      // Unitless
  reflectance_matrix = std::vector<std::vector<double>> (nbin, vector<double> (nbin, 0));
  reflectivity_matrix = std::vector<std::vector<double>> (nbin, vector<double> (nbin, 0));

	// Loop over the density matrix
	for (int i = 0; i < nbin; i++){
		reflectance = 0, reflectivity = 0;
		for (int j = 0; j < nbin; j++){

// TO CHECK: dr = r_bin is valid? exact ????
// Previous definition
      // reflectance = (1-reflectivity)*(1-exp(-1*cs.get_r_bin()/skin_depth(density_tx[j][i])));
      // reflectivity += reflectance;

// New definition
      reflectance = (1-exp(-1*cs.get_r_bin()*absorption( density_tx[j][i] ) ));
			reflectivity += (1-reflectivity)*reflectance;

      // std::cout << i << std::endl;
      // std::cout << j << std::endl;

			assert(reflectivity < 1 && "Reflectivity larger than 1!");

      reflectance_matrix[j][i] = reflectance;
      reflectivity_matrix[j][i] = reflectivity;

		}
	}
}

void Cascade1D::set_reflectivity(){set_reflectance(); }

/* Compute the overdense area of the cascade in slices.
Integrate the reflectance matrix over the l.o.s. (sum over all the segments)
The cell (parallelogram) has constant size idependently of the tilt, but it
changes with r!
*/

void Cascade1D::set_radar_cs(){
  double r;
  if (reflectivity_matrix.empty()) {set_reflectance();}
  _rcs = std::vector<double> (nbin);
  for (int i = 0; i < nbin; i++){
    for (int j = 0; j < nbin; j++){
      r = (1 - abs(2.0*j/nbin - 1)) * cs.get_r_tot();    // Radius between +r and -r.
      // CHANGE THIS BY A SINE FUNCTION
      _rcs[i] += r*cs.get_L_bin()*100 *reflectance_matrix[i][j];      // [cm^2]
    }
  }
}
// VALID FOR EDGE ON AND FACE ON?! (If we correct the loop order)

// The reflectance corrections are only applied to the od sections.
// std::vector<double> r_crit = cs.get_rcrit();
// // Number of layers per segment.
// int k_max = (int) r_crit[i]/r_bin;
 // Only layers from od_region are looked, outside-in.
// k >= (nbin - k_max) ? k_mid = nbin - (k + 0.5) : k_mid = 0;

// for(int i=0; i < rcs.size(); i++) {std::cout << rcs.at(i) << ' ';}
// std::remove_copy(reflectivity_matrix.back().begin(), reflectivity_matrix.back().end(), v2.begin(), 0);
// std::cout << v2[0] << std::endl;
// std::cout << rcs.empty() << std::endl;
// std::cout << cell_size << std::endl;
// std::cout << rcs[0] << " " << rcs[nbin - 1] << std::endl;


std::vector<std::vector<double>> Cascade1D::get_density_cs(){ return density_cs; }
std::vector<std::vector<double>> Cascade1D::get_density_tx(){ return density_tx; }

std::vector<std::vector<double>> Cascade1D::get_reflectance(){ return reflectance_matrix; }
std::vector<std::vector<double>> Cascade1D::get_reflectivity(){ return reflectivity_matrix; }

std::vector<double> Cascade1D::radar_cs(){
  if (_rcs.empty()) {set_radar_cs();}
  return _rcs;
}
// =============================================================================
