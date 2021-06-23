#include "macro_scatter.hh"

Scatter::Scatter(Antenna& tx, Antenna& rx, Cascade& cs):
  fTX(tx), fRX(rx), fCS(cs){

  /* First: Set the antennas directions, module and dot product with cs */
  fTX.SetDirection( fCS.Pos() );
  fRX.SetDirection( fCS.Pos() );

  /* Second: Build the TX frame.
  Find the angle delta between the direction to the cascade and the cascade's direction.

  The density is computed not in  the lab frame, but in the plane between the
  tx.direction vector and the cs direction vector. Each one defines a frame
  with their own perpendicular direction.

  The two frames are separated by an angle delta. This angle delta behaves like
  the declination, is only defined between 0 and pi. Due to the cascade radial
  symmetry, and the choice of plane (that cuts the cascade longitudinally)
  the solutions for delta and minus delta in the plane frame are equivalent.
  */

  // Determine inner product  in l.o.s plane.
  fDot   = projection(fTX.Dir(), fCS.Dir());
  fDelta = acos(fDot);
  // fDelta is defined between 0 and pi only.

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
  (fDelta < pi/180.0) ? fDelta = pi/180.0 : 1;

  // We can avoid computing the same values thousands of times.
  cD = cos(fDelta);
  sD = sin(fDelta);

  // These are the dimensions of the TX frame axis. [cm]
  fPar  = fCS.Ltot()*abs(cD) + 2*fCS.Rtot()*abs(sD);
  fPerp = fCS.Ltot()*abs(sD) + 2*fCS.Rtot()*abs(cD);

  // nbins  = size / division.
  nPar  = (int) ceil(fPar  / gdPar);
  nPerp = (int) ceil(fPerp / gdPerp);


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
    // l = i * dSeg;
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
                              cross_product(normalize(tx_seg), fTX.Pol() )
                                         )
                                       )
                                     )
                      , fRX.Pol() );

    fAttenuation[i] = ( 1 / ( rt*rr ) * pow(e, -(rt + rr)/(2*att_length) ) * fPolarization[i] );

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
}

/* Run time loop */
// Requires set_segments and set_radar_cs();
void Scatter::run_time_loop(){

  double nSeg = fRCS.size();
  // cout << nSeg << endl;

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
  for (int s = 0; s < steps; s++){
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


//        if(flag = 0)skip{}
          fWaveformTime[s][i] = E0 * fAttenuation[i] * sqrt(fRCS[i]) / 100.0 *
                                cos(fPhase[i] - fTX.AngularFreq()*t);
          // RCS is computed in cm^2 and we are moving now to m
          // if (fRCS[0]!= 1) { std::cout << fRCS[i] << std::endl;}
        }

      // }  // Selection closing

    }

    // The final E field value for a given timestep is the sum of the effects of all segments.
    fWaveform[s] = std::accumulate(std::begin(fWaveformTime[s]), std::end(fWaveformTime[s]), 0.0);
    // if (fRCS[0]!= 1 && fWaveform[s] != 0) { std::cout << fWaveform[s] << std::endl;}
  }
}

// Accesors
Antenna Scatter::transmitter(){return fTX;}
Antenna Scatter::receiver(){return fRX;}
Cascade Scatter::cascade(){return fCS;}
std::vector<double> Scatter::ESA(){ return fRCS; }

  // Not set before set_segments()
std::vector<double> Scatter::phase(){return fPhase;}
std::vector<double> Scatter::Polarization(){ return fPolarization; }
std::vector<double> Scatter::Attenuation(){ return fAttenuation; }
std::vector<double> Scatter::arrivals(){ return fArrivalTime ; }
std::vector<std::vector<double>> Scatter::Coordinates(){return fSegmentCoord;}

  // Not set before run_time_loop()
std::vector<double> Scatter::duration(){ return fDuration; }
std::vector<double> Scatter::Waveform(){ return fWaveform; }
std::vector<std::vector<double>> Scatter::rcs_time(){ return fRCSTime; }
std::vector<std::vector<double>> Scatter::phase_time(){ return fPhaseTime; }
std::vector<std::vector<double>> Scatter::wave_time(){ return fWaveformTime; }

//==============================================================================

// -----------------------------------------------------------------------------
// THE ROTATED VERSION

Cascade1D::Cascade1D(Antenna& tx, Antenna& rx, Cascade& cs): Scatter(tx, rx, cs){

    // check first if RCS is in folder, if not, run all of This

    // From Density to RCS.
    SetCascadeCoodinates();
    Density(fCSLength, fCSRadius);
    PlasmaFreq(fDensity); // Not necessary
    Absorption(fDensity, fTX.Freq() );
    SkinDepth(fDensity, fTX.Freq() ); // Not necessary
    // Reflectance();
    Opacity(fAbsorption);
    // std::cout << "Check" << endl;
    ESA(fReflectance);

    //  iF THE ecs is in folder, load and jump here

    SetSegments(nPerp);
    run_time_loop();
  }

/* Make 2D-array of density profile in the TX frame
Transmitter frame (a, b) goes from (0 -> A , 0 -> B)
Cascade frame (L,R) goes to (0 -> L, -r/2 -> r/2)
The rotation between the frames needs to happen at cascade's center.
The formula below is translation -> rotation -> translation back.
NOTE: l,r in the function are meant to be positions in the cascade frame (after rotation)
but still move along TX frame.
*/
void Cascade1D::SetCascadeCoodinates(){
  double a, b;
  fCSLength = std::vector<std::vector<double>> (nPerp, std::vector<double> (nPar, 0));
  fCSRadius = std::vector<std::vector<double>> (nPerp, std::vector<double> (nPar, 0));

  for (int i = 0; i < nPerp; i++){
    b = i* gdPerp;
    for (int j = 0; j < nPar; j++){
      a = j * gdPar;

      fCSLength[i][j] = (a - fPar/2)*cD + (b - fPerp/2)*sD + fCS.Ltot()/2;
      fCSRadius[i][j] = -(a - fPar/2)*sD + (b - fPerp/2)*cD;

    }
  }
}


void Cascade1D::Density( const std::vector<std::vector<double>> &fCSLength,
                            const std::vector<std::vector<double>> &fCSRadius ){
  fDensity = std::vector<std::vector<double>> (fCSLength.size(),
              std::vector<double> (fCSLength[0].size(), 0.0));

  for (int i = 0; i < fDensity.size(); i++){
    for (int j = 0; j < fDensity[i].size(); j++){
      // Check if we are not too far out from the cascade direction
      if ( abs(fCSRadius[i][j]) <= fCS.Rtot() ) {
        fDensity[i][j] = fCS.Density(rho_ice*fCSLength[i][j], fCSRadius[i][j]);
      }
    }
  }
}

/* Computes r_crit, the radial values of the line where the plasma frequency
equals the detection frequency (overdense line).
The shower waist is by definition the maximum radius for the same plasma density.
*/
/* [(#e-) cm^-3]  Critical e- density for overdense scattering condition
wp > w > 8980*sqrt(ne), w is frequency [Hz]!! */
// void Cascade1D::Rcrit(const std::vector<std::vector<double>> &density, const double & freq){
//   double r, dens_crit = memp*pow(freq/8980,2);
//   fRcrit = std::vector<double>(density.size(), 0);
//
//   // Loop over density matrix rows and then columms.
//   for (int i = 0; i < density.size(); i++){
//     r = 0;
//     for (int j = 0; j < density[i].size(); j++){
//       // If overdense, remember the layer.
//       if (density[i][j] > dens_crit) { r = j;}   // Layer
//     }
//     // fRcrit[i] = r*fRdiv; // [cm]
//   }
//
//   fRwaist = *max_element(fRcrit.begin(), fRcrit.end());
// }

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
      reflectance = (1-opacity)*(1-exp(-1.0 * gdPar * absorption[i][j]));
      opacity += reflectance;
			assert(opacity < 1 && "Opacity larger than 1!");   // sanity check
      fReflectance[i][j] = reflectance;
      fOpacity[i][j] = opacity;
		}
	}
}

void Cascade1D::Opacity(const std::vector<std::vector<double>> &absorption){ Reflectance(absorption); }

/* Compute the RCS of the cascade from the slices.
NON constant layer size, need to loop over the radius.
Per radius, we consider scattering over a cylinder shell of area:

 Shell area = 2pi*base * height_segment = 2*pi*r * gdPerp;

Each radius r actually contributes to half a cylinder (r goes from -r_max ro r_max):

  Half-shell area = pi* r * gdPerp;

Each shell is weighted by its reflectivity.
*/
void Cascade1D::ESA(const std::vector<std::vector<double>> &reflectance){
  // double r;      // No longer constant, normal dimension.
  fRCS = std::vector<double> ( reflectance.size(), 0.0 );
  for (int i = 0; i < reflectance.size(); i++){
    for (int j = 0; j < reflectance[i].size(); j++){
      fRCS[i] += abs(fCSRadius[i][j]) * reflectance[i][j];      // [cm^2]
    }
    fRCS[i] *= pi * gdPerp;
  }
}

// Accesors
std::vector<std::vector<double>> Cascade1D::Radius() {return fCSRadius;}
std::vector<std::vector<double>> Cascade1D::Density()  { return fDensity; }
// double  Cascade1D::Rwaist(){ return fRwaist; }
std::vector<double> Cascade1D::Rcrit(){ return fRcrit; }
std::vector<std::vector<double>> Cascade1D::PlasmaFreq(){ return fPlasmaFrequency; }
std::vector<std::vector<double>> Cascade1D::Absorption(){ return fAbsorption; }
std::vector<std::vector<double>> Cascade1D::SkinDepth(){ return fSkinDepth; }
std::vector<std::vector<double>> Cascade1D::Reflectance(){ return fReflectance; }
std::vector<std::vector<double>> Cascade1D::Opacity(){ return fOpacity; }
std::vector<double> Cascade1D::ESA(){ return Scatter::ESA(); };

// std::vector<double> Cascade1D::radar_cs(){
//   if (fRCS.empty()) {set_radar_cs();}
//   return fRCS;
// }

// -----------------------------------------------------------------------------

/* Cylinder1D
This should be a close approximation to the ideal cascade. The radial profile
r(L) is not trivial for even simple geometries, so we have a cylinder instead.
A cylinder that has perfect reflectivity only scatters out of his outer shell.
The typical size for the cylinder is the moliere radius.

Outer cylinder shell area per segment = 2pibase* height = pi* r_moliere* gdPerp
Opacity is 1;
*/

Cylinder1D::Cylinder1D(Antenna& tx, Antenna& rx, Cascade& cs): Scatter(tx, rx, cs){
  fRCS = std::vector<double>(nPerp, gdPerp * r_moliere * pi); // [cm^2]
  SetSegments( nPerp );
  run_time_loop();
}
// -----------------------------------------------------------------------------

/* Line1D is Dieder's model.
This should be the closest representation to the ideal thin-wire solution
(there is an analytical solution that does not contain retardation effects)
*/
Line1D::Line1D(Antenna& tx, Antenna& rx, Cascade& cs): Scatter(tx, rx, cs){
  fRCS = std::vector<double>(fCS.Lbins(), gdNorm*gdPerp); // The line's segments have rcs unity.
  SetSegments( fCS.Lbins() );
  run_time_loop();
}
// =============================================================================
