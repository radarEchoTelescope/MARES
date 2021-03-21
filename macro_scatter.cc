#include "macro_scatter.hh"

Scatter::Scatter(Antenna& tx, Antenna& rx, Cascade& cs):
  fTX(tx), fRX(rx), fCS(cs){

  /* Set the antennas directions, module and dot product with cs */
  fTX.SetDirection( fCS.Pos() );
  fRX.SetDirection( fCS.Pos() );
  /* Attenuation model example
  // Parametrized attenuation length for the Ross Ice Shelf, South Pole.
  double att(double freq_obs){
  double a1=469;                  // [m] Attanuation length parameter
  double a2=-0.205;               // Attanuation length parameter
  double a3=4.87E-5;              // Attanuation length parameter
  return a1+a2*freq_obs/1E6+a3*pow(freq_obs/1E6,2);
  }
  */


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


/* Loop over the 1D segments  */
void Scatter::SetSegments( const double& nSeg){
  double l, rt, rr, rt0, dSeg;
  std::vector<double> seg_pos{0,0,0}, tx_seg{0,0,0}, rx_seg{0,0,0};

  fPhase        = std::vector<double> ( nSeg, 0 );
  fArrivalTime  = std::vector<double> ( nSeg, 0 );
  fAttenuation  = std::vector<double> ( nSeg, 0 );
  fPolarization = std::vector<double> ( nSeg, 0 );
  fSegmentCoord = std::vector<std::vector<double>>(nSeg, vector<double> (6, 0));

  // Giving a non-uniform position in the cascade gets rid of artifacts in the FFT.
  RN_uniform rand_line(-0.5, 0.5, 42);          // Same seed for debugging.
  // RN_uniform rand_line(-0.5, 0.5, time(0));  // Different seed for random, independent runs.

  // This way we are making sure that we are covering the whole cascade length
  dSeg = fCS.Ltot()/100 / nSeg; // [m]

  // Sanity check! Your sections are not unphysical due to lifetime constraint.
  assert(dSeg <= cvac_cm*tau && "Cascade segments are too large!");
  // This should be further checked against the probing wavelength.

  // Segment loop
  for(int i = 0 ; i < nSeg; i++){
    // Set position
    l = (i + rand_line.get()) * dSeg;
    // [m] distance from the shower head (starting point)

    seg_pos[0]  = fCS.Pos()[0] + l*fCS.Dir()[0];
    seg_pos[1]  = fCS.Pos()[1] + l*fCS.Dir()[1];
    seg_pos[2]  = fCS.Pos()[2] + l*fCS.Dir()[2];

    // Set distances
    tx_seg = direction(fTX.Pos(), seg_pos);
    rt = norm(tx_seg);

    if(i == 0){rt0 = rt;}

    rx_seg = direction(seg_pos, fRX.Pos());
    rr = norm(rx_seg);

    fSegmentCoord[i][0] = l;
    fSegmentCoord[i][1] = seg_pos[0];
    fSegmentCoord[i][2] = seg_pos[1];
    fSegmentCoord[i][3] = seg_pos[2];
    fSegmentCoord[i][4] = rt;
    fSegmentCoord[i][5] = rr;

    // Set times
    // T0 = 0 by definition when the cascade begins (first element = head).
    //birth = l/c_vac; // time where the i'th segment starts scattering.

    // (Retarded) time where the scattered radio signal by the segment is produced.
    // production = birth - Rt/c_ice;

    // (Advanced) time where the scattered signal by the segments arrives in the receiver.
    fArrivalTime[i] = l/c_vac + rr/c_ice;

    // Set phase
    fPhase[i] = fTX.Wavenumber()*(2*rt + rr - rt0) ; // Simplified from Dieder code
    // fPhase.push_back( tx.wavenumber()*(2*(rt + rr) - rt0 - l/refindex) ); // Simplified form Dieder's thesis.

    /* Dieder's phase from the thesis
    phase = k*(rt + rr) - AngularFreq(l/c_vac + (rt(0) - rt)/c_ice - rr/c_ice) (- AngularFreq*t)
    // fPhase.push_back(tx.wavenumber()*(rt + rr) - tx.AngularFreq()*(l/c_vac + (rt0 - rt)/c_ice - rr/c_ice));

    where:

      k*(rt + rr) = tx.wavenumber()*(rt + rr); spatial phase
      l/c_vac; cascade propagation term
      (rt(0) - rt)/c_ice; change in phase for each point w.r.t starting point (i = 0).
      rt(0) = | cs_pos - tx.pos | = cs.distance if computed from the head.
      rr/c_ice; retardation effects ( = time@receiver - t@cascade)

      This becomes ( with AngularFreq/c_ice = k):
      // fPhase.push_back( tx.wavenumber()*(2*(rt + rr) - rt0 - l/refindex) );
    */


    // E field attenuation at reciever from position dependant factors:
    // 1/R from electric field.
    // e^-r/Latt from medium attenuation
    // p from signal polarization

    // Set polarization
    fPolarization[i] = projection(
                        cross_product(normalize(rx_seg),
                          cross_product(normalize(rx_seg),
                            cross_product(normalize(tx_seg),
                              cross_product(normalize(tx_seg), fTX.Pol()
                                           )
                                         )
                                       )
                                     )
                      , fRX.Pol() );


    fAttenuation[i] = ( 1 / ( rt*rr ) * pow(e, -(rt + rr)/(2*att_length) )
                      * fPolarization[i] );

  }
}

/* Run time loop */
// Requires set_segments and set_radar_cs();
void Scatter::run_time_loop(){

  double nSeg = fRCS.size();

  double t, E0;
  double t_start  = *min_element(fArrivalTime.begin(), fArrivalTime.end()) - 5E-9;
  double t_end    = *max_element(fArrivalTime.begin(), fArrivalTime.end()) + tau + 5E-9;

  double sampling = (100*fTX.Freq());
  int steps = (t_end - t_start)*sampling;

  fDuration     = std::vector<double>(steps, 0);    // The time
  fWaveform     = std::vector<double>(steps, 0);    // The electric field
  fPhaseTime    = std::vector<std::vector<double>>(steps, vector<double> (nSeg, 0));
  fRCSTime      = std::vector<std::vector<double>>(steps, vector<double> (nSeg, 0));
  fWaveformTime = std::vector<std::vector<double>>(steps, vector<double> (nSeg, 0));

  /* THE RADAR RETURN EQUATION AT ELECTRIC FIELD LEVEL
  Er = E0 / tx.wavelength() * (leff_tx)*(leff_rx) * Z0/Zload * sqrt(eta_tx/(4*pi) *        (Antenna/s constants)
      1/ ( rt*rr) * e^(abs(rt + rr)/2*att_length )              (Position dependant variables, atenuattion)
      rx_seg X (rx_seg X [tx_seg X (tx_seg X tx_pol)]) · rx_pol      (Polarization = sin * sin * cos)
      sqrt(sigmafRCS);                                            (The cross section, computed later)
  */

  // Position independent factors
  E0 =  fTX.E0()/fTX.Wavelength() * (fTX.Leff()) * (fRX.Leff()) * Z0/fRX.Load() *
        sqrt(fTX.Eff()/(4*pi));

  // Loop over time.
  for (int s = 0; s < steps; s ++){
    t = s/sampling + t_start;
    fDuration[s] = t;
    for (int i = 0; i < nSeg; i++){

      // Select some length values?
       // if(i == 0){

        // If active, add its contribution.
        if(t>fArrivalTime[i] && t<(fArrivalTime[i]+tau)){

          fRCSTime[s][i] = fRCS[i];
          fPhaseTime[s][i] = fPhase[i] - fTX.AngularFreq()*t;
          //fPhaseTime[s][i] = fPhase[i] - fTX.AngularFreq()*(t-t_start);

          fWaveformTime[s][i] = E0 * fAttenuation[i] * sqrt(fRCS[i])/100 *
                                cos(fPhase[i] - fTX.AngularFreq()*t);
          // RCS is computed in cm^2 and we are moving now to m
        }

      // }  // Selection closing

    }

    // The final E field value is the sum of all particles for the timestep.
    fWaveform[s] = std::accumulate(std::begin(fWaveformTime[s]), std::end(fWaveformTime[s]), 0.0);
  }
}

// Accesors
Antenna Scatter::transmitter(){return fTX;}
Antenna Scatter::receiver(){return fRX;}
Cascade Scatter::cascade(){return fCS;}
std::vector<double> Scatter::RCS(){ return fRCS; }

  // Not set before set_segments()
std::vector<double> Scatter::phase(){return fPhase;}
std::vector<double> Scatter::Attenuation(){ return fAttenuation; }
std::vector<double> Scatter::arrivals(){ return fArrivalTime ; }
std::vector<std::vector<double>> Scatter::Coordinates(){return fSegmentCoord;}

  // Not set before run_time_loop()
std::vector<double> Scatter::duration(){ return fDuration ; }
std::vector<double> Scatter::waveform(){ return fWaveform ; }
std::vector<std::vector<double>> Scatter::rcs_time(){ return fRCSTime; }
std::vector<std::vector<double>> Scatter::phase_time(){ return fPhaseTime; }
std::vector<std::vector<double>> Scatter::wave_time(){ return fWaveformTime; }

// // --------------------------------------------------------------------------

//==============================================================================
// Line1D is Dieder's model

Line1D::Line1D(Antenna& tx, Antenna& rx, Cascade& cs): Scatter(tx, rx, cs){

  fRCS = std::vector<double>(fCS.Lbins(), 1); // The line's segments have rcs unity.
  SetSegments( fCS.Lbins() );
  run_time_loop();
  }

// -----------------------------------------------------------------------------
// THE ROTATED VERSION

Cascade1D::Cascade1D(Antenna& tx, Antenna& rx, Cascade& cs): Scatter(tx, rx, cs){
    /*
    The density is computed not in  the lab frame, but in the plane between the
    tx.direction vector and the cs direction vector. Each one defines a frame
    with their own perpendicular direction.

    The two frames are separated by an angle delta. This angle delta behaves like
    the declination, is only defined between 0 and pi. Due to the cascade radial
    symmetry, and the choice of plane (that cuts the cascade longitudinally)
    the solutions for delta and minus delta in the plane frame are equivalent.

    Delta is actually the angle
    */
    //
    // std::cout << fTX.Dir()[0] << " " << fTX.Dir()[1] << " " << fTX.Dir()[2] << endl;
    // std::cout << fCS.Dir()[0] << " " << fCS.Dir()[1] << " " << fCS.Dir()[2] << endl;
    // // std::cout << cD << " " << sD << endl;

    // First: Find the angle between the direction to the cascade and the cascade's direction.
    // Determine inner product  in l.o.s plane.
    fDot   = projection(fTX.Dir(), fCS.Dir());
    fDelta = acos(fDot);
    // fDelta is defined between 0 and pi only.
    // We can avoid computing the same values thousands of times.
    cD = cos(fDelta);
    sD = sin(fDelta);

    // std::cout << fDot << " " << fDelta << endl;
    // std::cout << cD << " " << sD << endl;

    // These are the dimensions of the axis in the projection into the incidence frame.
    fPar  = fCS.Ltot()*abs(cD) + 2*fCS.Rtot()*abs(sD);
    fPerp = fCS.Ltot()*abs(sD) + 2*fCS.Rtot()*abs(cD);

    // nbins  = size / division.
    nPar  = (int) ceil(fPar  / dPar);
    nPerp = (int) ceil(fPerp / dPerp);

    // std::cout << fPerp << " " << fPar << endl;
    // std::cout << nPerp << " " << nPar << endl;
    // From Density to RCS.
    SetDensity();
    PlasmaFreq(fDensity); // Not necessary
    Absorption(fDensity, fTX.Freq() );
    SkinDepth(fDensity, fTX.Freq() ); // Not necessary
    // Reflectance();
    Opacity(fAbsorption);
    // std::cout << "Check" << endl;
    RCS(fOpacity);

    SetSegments(nPerp);
    // run_time_loop();
  }


/* Make 2D-array of density profile in the TX frame
Transmitter frame (a, b) goes from (0 -> A , 0 -> B)
Cascade frame (L,R) goes to (0 -> L, -r/2 -> r/2)
The rotation between the frames needs to happen at cascade's center.
The formula below is translation -> rotation -> translation back.
*/
void Cascade1D::SetDensity(){
	double a, b, l, r;
  fDensity = std::vector<std::vector<double>> (nPerp, std::vector<double> (nPar, 0));

  for (int i = 0; i < nPerp; i++){
    b = i* dPerp;
    for (int j = 0; j < nPar; j++){
      a = j * dPar;

      l =  (a - fPar/2)*cD + (b - fPerp/2)*sD + fCS.Ltot()/2;
      r = -(a - fPar/2)*sD + (b - fPerp/2)*cD;
      // std::cout << a << " " << b << endl;
      // std::cout << l << " " << r << endl;

      fDensity[i][j] = fCS.Density(l*rho_ice, r);
    }
  }
}
// NOTE: l,r in the function are meant to be positions in the cascade frame (after rotation)
// but still move along TX frame.

/* Computes r_crit, the radial values of the line where the plasma frequency
equals the detection frequency (overdense line).
The shower waist is by definition the maximum radius for the same plasma density.
*/
/* [(#e-) cm^-3]  Critical e- density for overdense scattering condition
wp > w > 8980*sqrt(ne), w is frequency [Hz]!! */
void Cascade1D::Rcrit(const std::vector<std::vector<double>> &density, const double & freq){
  double r, dens_crit = memp*pow(freq/8980,2);
  fRcrit = std::vector<double>(density.size(), 0);

  // Loop over density matrix rows and then columms.
  for (int i = 0; i < density.size(); i++){
    r = 0;
    for (int j = 0; j < density[i].size(); j++){
      // If overdense, remember the layer.
      if (density[i][j] > dens_crit) { r = j;}   // Layer
    }
    // fRcrit[i] = r*fRdiv; // [cm]
  }

  fRwaist = *max_element(fRcrit.begin(), fRcrit.end());
}

void Cascade1D::PlasmaFreq(const std::vector<std::vector<double>> &density){
  fPlasmaFrequency = std::vector<std::vector<double>> (density.size(),
                      vector<double> (density[0].size(), 0));

  for (int i = 0; i < density.size(); i++){
    for (int j = 0; j < density[i].size(); j++){
        fPlasmaFrequency[i][j] = fCS.PlasmaFreq(density[i][j]);
    }
  }
}

void Cascade1D::Absorption(const std::vector<std::vector<double>> &density, const double & freq){
  fAbsorption = std::vector<std::vector<double>> (density.size(),
                  vector<double> (density[0].size(), 0));

  for (int i = 0; i < density.size(); i++){
    for (int j = 0; j < density[i].size(); j++){
        fAbsorption[i][j] = fCS.Absorption(density[i][j], freq);
    }
  }
}

void Cascade1D::SkinDepth(const std::vector<std::vector<double>> &density, const double & freq){
  fSkinDepth = std::vector<std::vector<double>> (density.size(),
                vector<double> (density[0].size(), 0));

  for (int i = 0; i < density.size(); i++){
    for (int j = 0; j < density[i].size(); j++){
        fSkinDepth[i][j] = fCS.SkinDepth(density[i][j], freq);
    }
  }
}

/* Compute the 2D reflectance and opacity for varying skin depth.

Physics note:
Opacity = Property of a material.
Reflectance = Refers to a specific sample, depends on size and other params.

Opacity is the reflectance value as the object becomes thick.
Thus, opacity is defined here as the integral of reflectance over
the layers.
*/

void Cascade1D::Reflectance(const std::vector<std::vector<double>> &absorption){
	double reflectance, opacity;      // Unitless, tmp variables.
  fReflectance = std::vector<std::vector<double>> (absorption.size(),
                  vector<double> (absorption[0].size(), 0));
  fOpacity = std::vector<std::vector<double>> (absorption.size(),
                  vector<double> (absorption[0].size(), 0));

	// Loop over the absorption matrix
  for (int i = 0; i < absorption.size(); i++){
		reflectance = 0, opacity = 0;
    for (int j = 0; j < absorption[0].size(); j++){
      // The bin size [in cm] is the parallel dimension resolution.
      reflectance = (1-opacity)*(1-exp(-1.0 * dPar * absorption[i][j]));
      opacity += reflectance;
			assert(opacity < 1 && "Opacity larger than 1!");   // sanity check
      fReflectance[i][j] = reflectance;
      fOpacity[i][j] = opacity;
		}
	}
}

void Cascade1D::Opacity(const std::vector<std::vector<double>> &absorption){ Reflectance(absorption); }

/* Compute the RCS of the cascade from the slices.
For a constant cell size with layers, we just need the opacity (accumulated reflectivity)
and multiply by the cell area.
*/
void Cascade1D::RCS(const std::vector<std::vector<double>> &opacity){
  fRCS = std::vector<double> ( opacity.size() );
  std::vector<double> op = opacity.back();
  fRCS = dPerp * dNorm * op;
}


/* Compute the RCS of the cascade from the slices.
For a NON constant cell size w.r.t layers, need to loop over the reflectivity
and weigh the cell by its respective area before adding.
*/
// void Cascade1D::RCS(const std::vector<std::vector<double>> &reflectance){
//   double r;     // No longer constant, dNorm.
//   fRCS = std::vector<double> ( reflectance.size() );
//   for (int i = 0; i < reflectance.size(); i++){
//     for (int j = 0; j < reflectance[i].size(); j++){
//       // Whatever parametrization for dNorm that we want goes here [to be investigated]
//
//       // r = -1.0 *- sin ( pi* (2.0*j/nbin - 1) ) * dNorm; // e.g.
//       // r = abs(r);        // Radius evaluated in [+r, -r] but we need the module
//       // fRCS[i] += r * dPerp * reflectance[i][j];      // [cm^2]
//     }
//   }
// }

// ----------------------------------------------------------------------------
// Accesors
std::vector<std::vector<double>> Cascade1D::Density()  { return fDensity; }
double  Cascade1D::Rwaist(){ return fRwaist; }
std::vector<double> Cascade1D::Rcrit(){ return fRcrit; }
std::vector<std::vector<double>> Cascade1D::PlasmaFreq(){ return fPlasmaFrequency; }
std::vector<std::vector<double>> Cascade1D::Absorption(){ return fAbsorption; }
std::vector<std::vector<double>> Cascade1D::SkinDepth(){ return fSkinDepth; }
std::vector<std::vector<double>> Cascade1D::Reflectance(){ return fReflectance; }
std::vector<std::vector<double>> Cascade1D::Opacity(){ return fOpacity; }
std::vector<double> Cascade1D::RCS(){ return Scatter::RCS(); };

// std::vector<double> Cascade1D::radar_cs(){
//   if (fRCS.empty()) {set_radar_cs();}
//   return fRCS;
// }
// =============================================================================
