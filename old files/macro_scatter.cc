#include "macro_scatter.hh"

Scatter::Scatter(Antenna& tx, Antenna& rx, Cascade& cs,
        const double deltaL, const double deltaR,  const double deltaN,
        const double lifetime, const double sampling):

  fTX(tx), fRX(rx), fCS(cs), dL(deltaL), dR(deltaR), dN(deltaN),
  tau(lifetime), sampling_ratio(sampling){

    // Sanity check! Your sections are not unphysical due to lifetime constraint.
 	 assert(deltaL <= c_vac*tau && "Cascade resolution is too large!");
 	 // This should be further checked against the probing wavelength.

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

  fTX.SetAngle( fCS.Dir() );
  fRX.SetAngle( fCS.Dir() );


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


  // We can avoid computing the same values thousands of times.
  cD = cos( fTX.Delta() );
  sD = sin( fTX.Delta() );

  // These are the dimensions of the TX frame axis. [cm]
  R = fCS.Ltot()*abs(cD) + 2*fCS.Rtot()*abs(sD);
  L = fCS.Ltot()*abs(sD) + 2*fCS.Rtot()*abs(cD);

  /*
  If delta = pi/2: Perpendicular incidence
    - L (TX frame) = L (CS frame)
    - R (TX frame) = 2R (CS frame)

  If delta = 0:
    - L (TX frame) = 2R (CS fame)
    - R (TX frame) = L (CS frame)
  */

  // dL = fCS.Ldiv();
  // dR = fCS.Rdiv();
  // dN = fCS.Rdiv();

  // nbins  = size / division.
  nR  = (int) ceil(R  / dR);
  nL  = (int) ceil(L / dL);
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
void Scatter::SetSegments( const int& nSeg){
  double l, rt, rr, rt0, dSeg;
  std::vector<double> R_TX_dir{0,0,0}, L_TX_dir{0,0,0}, seg_pos{0,0,0}, tx_seg{0,0,0}, rx_seg{0,0,0};

  fPhase        = std::vector<double> ( nSeg, 0 );
  fArrivalTime  = std::vector<double> ( nSeg, 0 );
  fAttenuation  = std::vector<double> ( nSeg, 0 );
  fPolarization = std::vector<double> ( nSeg, 0 );
  fDirectivity  = std::vector<double> ( nSeg, 0 );
  fSegmentCoord = std::vector<std::vector<double>>(nSeg, std::vector<double> (6, 0));

  // Giving a non-uniform position in the cascade gets rid of artifacts in the FFT.
  RN_uniform rand_line(-0.5, 0.5, 42);          // Same seed for debugging.
  // RN_uniform rand_line(-0.5, 0.5, time(0));  // Different seed for random, independent runs.

  // This way we are making sure that we are covering the whole cascade length
  // dSeg = fCS.Ltot()/ nSeg / 100.0; // [m]

  // Sanity check! Your sections are not unphysical due to lifetime constraint.
  // assert(dSeg <= c_vac*tau && "Cascade segments are too large!");
  // This should be further checked against the probing wavelength.



  // Segment loop
  for(int i = 0 ; i < nSeg; i++){
    //     // Set position
//     l = (i + rand_line.get()) * dSeg;
//     // l = i * dSeg;
//     // [m] distance from the shower head (starting point)
//
// // TEST = Small randomisation (0.1 mm) in all directions.
//     seg_pos[0]  = fCS.Pos()[0] + l*fCS.Dir()[0] + 1E-2*rand_line.get();
//     seg_pos[1]  = fCS.Pos()[1] + l*fCS.Dir()[1] + 1E-2*rand_line.get();
//     seg_pos[2]  = fCS.Pos()[2] + l*fCS.Dir()[2] + 1E-2*rand_line.get();
//
// //    seg_pos[0]  = fCS.Pos()[0] + l*fCS.Dir()[0];
// //    seg_pos[1]  = fCS.Pos()[1] + l*fCS.Dir()[1];
// //    seg_pos[2]  = fCS.Pos()[2] + l*fCS.Dir()[2];

    // Set position
    // [m] distance from the shower head (starting point)


/* So far, the plane R_TX, L_TX was created to cover the cascade at an angle and
the values within were rotated to find the values in the cascade frame (density, etc).
For geometry, a point P with coordinates (r,l) in the TX frame is placed in the
3D lab frame at:

P = CS_vertex + r*e(R_TX) +l*e(L_TX) (capitals == vector)

For each segment in the L_TX direction, we can find the coordinates (r_max,l_max)
of the point where the reflectivity is the highest (it should correlate with the
peak density for the segment, too), and place the segment at that point.

The choice of L_TX is the direction normal to R_TX that is consistent with the
rotation matrix used in SetCascadeCoodinates(), which is the clockwise rotation
of the TX frame with positive angle, or the cascade rotating anti-clockwise from
the frame with positive delta.

*/

    R_TX_dir = normalize(fTX.Dir()) ;
    L_TX_dir = cross_product( cross_product(R_TX_dir, normalize(fCS.Dir()) ), R_TX_dir);
// Is this last step needed? I dont think so but just in case.
    L_TX_dir = normalize(L_TX_dir);

    // fTXMaxLength, fTXMaxRadius should have been computed in cm
    // fCSPos is given in m

    seg_pos[0]  = fTXMaxRadius[i]*R_TX_dir[0] + fTXMaxLength[i]*L_TX_dir[0];
    seg_pos[1]  = fTXMaxRadius[i]*R_TX_dir[1] + fTXMaxLength[i]*L_TX_dir[1];
    seg_pos[2]  = fTXMaxRadius[i]*R_TX_dir[2] + fTXMaxLength[i]*L_TX_dir[2];

    // fTXMaxLength, fTXMaxRadius should have been computed in cm
    // fCSPos is given in m
    seg_pos[0]  = fCS.Pos()[0] + seg_pos[0]*cm/m;
    seg_pos[1]  = fCS.Pos()[1] + seg_pos[1]*cm/m;
    seg_pos[2]  = fCS.Pos()[2] + seg_pos[2]*cm/m;
    // Now EVERYTHING is in m.


    // Set distances
    tx_seg = direction(fTX.Pos(), seg_pos);
    rt = norm(tx_seg);

    // if(i == 0){rt0 = rt;}

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
    // fArrivalTime[i] = rr/c_ice;

    // Set phase
    fPhase[i] = fTX.Wavenumber()*(rt + rr) + atan2( f_coll, - fTX.Freq() );
    // tan^-1(f_coll/f_TX) is the phase shift in the oscillation caused by collisions.
    // in our regime, it could be almost fixed to -Pi/2.

    // E field attenuation at reciever from position dependant factors:
    // 1/R from electric field.
    // e^-r/Latt from medium attenuation
    fAttenuation[i] =  1 / ( rt*rr ) * pow(e, -(rt + rr)/(2*att_length) ) ;

    /* Attenuation model example
    // Parametrized attenuation length for the Ross Ice Shelf, South Pole.
    double att(double freq_obs){
      double a1=469;                  // [m] Attanuation length parameter
      double a2=-0.205;               // Attanuation length parameter
      double a3=4.87E-5;              // Attanuation length parameter
      return a1+a2*freq_obs/1E6+a3*pow(freq_obs/1E6,2);
    }
    */

    // Set segment directivity for a square aperture
    // Total directivity = Max directivty * normalised gain;

    // theta and phi are the spherical coordinates for the rx_seg
    // (orfanidi, fig 18.4.1)
    // double theta, phi;
    // theta =  acos( rx_seg[2] / rr);
    // phi = atan2( rx_seg[1], rx_seg[0]);
    // The zenith angle is defined between 0 and pi so its sign safe.
    // atan2 has built in corrections for the signs of the angle.

    /* The segment normal is oriented parallel to tx_seg
    ( (1 + cos( delta ) ) / 2 )^2 is the obliquity factor,
    delta is the angle between tx and rx directions
    cos (delta) = cos(acos(projection(tx_seg, rx_seg))).
    acos is really unstable around -1, so we just avoid it. */

    // std::vector<double> cs_normal = cross_product(fCS.Dir(), cross_product(normalize(tx_seg), fCS.Dir() ) );

    /* Aperture-like directivty
      Aperture peak approx:  4*pi*dL*dR/pow((fTX.Wavelength()*100),2)

      Full aperture directivty:
        pow( sinc( pi * dL/ (fTX.Wavenumber()*100) * sin(theta)*sin(phi) ), 2)*
        pow( sinc( pi * dN/ (fTX.Wavelength()*100)  * sin(theta)*cos(phi) ), 2)

      Obliquity factor (3 definitions):
        pow( ((projection(cs_normal, tx_seg) ) - projection(cs_normal, rx_seg) ) / 2.0, 2) *
        pow((1 - projection(rx_seg, tx_seg) ) / 2.0, 2) *
        pow((1 + projection(rx_seg, tx_seg) ) / 2.0, 2) *
    */

    // Segment's polarization direction == E field at segment.
    std::vector<double> seg_pol = cross_product(normalize(tx_seg), cross_product(normalize(tx_seg), fTX.Pol() ) );
    // Direction of electric field at the receiver.
    std::vector<double> rx_Efield = cross_product(normalize(rx_seg), cross_product(normalize(rx_seg), seg_pol ) );

    // Directivty is hardcoded as a small Herztian dipole.
    fDirectivity[i] = 1.5*norm(rx_Efield); // 3/2*sin(theta_R)*sin(theta_T)

    // Set polarization efficiency
    fPolarization[i] = abs(projection(rx_Efield, fRX.Pol()));
    //In the thin-wire theory, only the  component of  the electric- field vector parallel to the wire
    // axis can interact to form a scattered  wave. That is not our case, our layers will scatter as a free charge
  }
}






/* Run time loop */
// Requires set_segments and set_radar_cs();
void Scatter::RunScatter(){
  double nSeg = fTCS.size(); // nL

  double t;
  double t_start  = *min_element(fArrivalTime.begin(), fArrivalTime.end()) - 5E-9 ;
  double t_end    = *max_element(fArrivalTime.begin(), fArrivalTime.end()) + tau +5E-9;

  double freq_sampling = fTX.Freq()*sampling_ratio;
  int steps = (t_end - t_start)*freq_sampling;


  fDuration     = std::vector<double>(steps, 0);    // The time
  fWaveform     = std::vector<double>(steps, 0);    // The electric field
  fPower        = std::vector<double>(steps, 0);    // The electric field
  fRCS          = std::vector<double>(steps, 0);    // The electric field
  fPhaseTime    = std::vector<std::vector<double>>(steps, std::vector<double> (nSeg, 0));
  fRCSTime      = std::vector<std::vector<double>>(steps, std::vector<double> (nSeg, 0));
  fTCSTime      = std::vector<std::vector<double>>(steps, std::vector<double> (nSeg, 0));
  fWaveformTime = std::vector<std::vector<double>>(steps, std::vector<double> (nSeg, 0));

  /* THE RADAR RETURN EQUATION AT ELECTRIC FIELD LEVEL
  Er = E0 (Antenna/s constants)
      1/ ( rt*rr) * e^(abs(rt + rr)/2*att_length )              (Position dependant variables, atenuattion)
      rx_seg X (rx_seg X [tx_seg X (tx_seg X tx_pol)]) · rx_pol      (Polarization = sin * sin * cos)
      sqrt(sigmafTCS);                                            (The cross section, computed later)
  */

  // Flag variable to select just some particular wavefrom during debugging.
  // int flag = 0;

  // Time Loop!
  for (int ts = 0; ts < steps; ts++){
    t = ts/freq_sampling + t_start;

    fDuration[ts] = t;
    for (int i = 0; i < nSeg; i++){

      // Select some length values?
      // if(i == 0){


        // If active, add its contribution.
        if(t>fArrivalTime[i] && t<=(fArrivalTime[i]+3*tau)){

          fTCSTime[ts][i] = fTCS[i];
          fPhaseTime[ts][i] = cos(fPhase[i] - fTX.AngularFreq()*t);

          // Compute the RCS here
          fRCSTime[ts][i] = sqrt(fTCS[i]  *      // A * rho
                                fDirectivity[i]       // Directivity of a segment
                                ) / 100.0 *           // [m^2]
                            pow(e,-(t-fArrivalTime[i])/tau)   // Lifetime decay
                            * fPhaseTime[ts][i];

  // TCS was computed in cm^2 and we are moving now to m^2 outside of the sqrt)

          // If (flag), find the waveform to store it
          fWaveformTime[ts][i] = 1.0/(4*pi) *
                                sqrt(2 * Z_ice * fTX.Power()  *
                                  fTX.GainDipole( fTX.Delta() ) *
                                  fRX.GainDipole( fRX.Delta() )
                                ) *
                                fAttenuation[i] *
                                fPolarization[i] *
                                fRCSTime[ts][i];
                              }

        // }  // Selection closing

      // }

      // The final RCS, E field value for a given timestep is the sum of the effects of all segments.
      fRCS[ts] = std::accumulate(std::begin(fRCSTime[ts]), std::end(fRCSTime[ts]), 0.0);
      fRCS[ts] = pow(fRCS[ts],2);

      fWaveform[ts] = std::accumulate(std::begin(fWaveformTime[ts]), std::end(fWaveformTime[ts]), 0.0);
      // if (fTCS[0]!= 1 && fWaveform[ts] != 0) { std::cout << fWaveform[ts] << std::endl;}

      //Alternative E field definition, from RCS_time

      // fWaveform[ts] = 1/(4*pi) * fAttenuation[i] * fPolarization[i]:
      //                sqrt(2 * Z_ice * fTX.Power() * fTX.GainDipole( fTX.Delta() ) * fRX.GainDipole( fRX.Delta() ) ) *
                     // fRCSTime[ts];

      fPower[ts] = 1.0/(2*Z_ice) * pow(fWaveform[ts],2)* pow(fTX.Wavelength(),2) /(4*pi) ;

      // fRCS[ts] = pow(4 * pi , 2)/
      // (fTX.Power() * fTX.GainDipole( fTX.Delta() ) * fRX.GainDipole( fRX.Delta() )
      // * pow(fAttenuation[0],2) // Already includes RT, RR at E level
      // * fPolarization[0]
      // * pow(fTX.Wavelength(),2) /(4*pi) ) * fPower[ts];

    }
  }

}
// Accesors
Antenna Scatter::TX(){return fTX;}
Antenna Scatter::RX(){return fRX;}
Cascade Scatter::CS(){return fCS;}
std::vector<double> Scatter::TCS(){ return fTCS; }

  // Defined at SetSegments()
std::vector<double> Scatter::Phase(){return fPhase;}
std::vector<double> Scatter::ArrivalTime(){ return fArrivalTime ; }
std::vector<double> Scatter::Attenuation(){ return fAttenuation; }
std::vector<double> Scatter::Directivity(){ return fDirectivity;}
std::vector<double> Scatter::Polarization(){ return fPolarization; }
std::vector<std::vector<double>> Scatter::Coordinates(){return fSegmentCoord;}

  // Not set before RunScatter()
std::vector<double> Scatter::Duration(){ return fDuration; }
std::vector<double> Scatter::Waveform(){ return fWaveform; }
std::vector<double> Scatter::Power(){return fPower;}
std::vector<double> Scatter::RCS(){return fRCS;}
std::vector<std::vector<double>> Scatter::TCS_time(){return fTCSTime;}
std::vector<std::vector<double>> Scatter::RCS_time(){ return fRCSTime; }
std::vector<std::vector<double>> Scatter::Phase_time(){ return fPhaseTime; }
std::vector<std::vector<double>> Scatter::E_time(){ return fWaveformTime; }

//==============================================================================
// -----------------------------------------------------------------------------

/* Line1D is Dieder's model.
This should be the closest representation to the ideal thin-wire solution
(there is an analytical solution that does not contain retardation effects)
*/
Line1D::Line1D(Antenna& tx, Antenna& rx, Cascade& cs): Scatter(tx, rx, cs){

  // std::cout << 2* 1.23 * fCS.Rtot()*fCS.Rtot()*2.0 <<'\t' << '\n';

  fTCS = std::vector<double>(nL, dL*dN*2.0); // [cm^2]
  //The line's segements are is 2a*a, with perfect reflectivity.
  // A line is expected to have Different polarization than a cascade
  SetSegments( nL );
  RunScatter();
}

// -----------------------------------------------------------------------------
Cascade1D::Cascade1D(Antenna& tx, Antenna& rx, Cascade& cs,
        const double deltaL, const double deltaR, const double deltaN,
        const double lifetime, const double sampling):
        Scatter(tx, rx, cs, deltaL, deltaR, deltaN, lifetime, sampling){

    // Missing for multi-receiver setups:
    //  Figure out if TCS has been computed already.

    fDamping = 1.0 / sqrt( pow(fTX.Freq(), 4) + pow(fTX.Freq()*f_coll,2) );

    // From Density to TCS.
    SetCascadeCoodinates();


    // Density_beam TEST!!!


    Density2(fCSLength, fCSRadius);
    PlasmaFreq(fDensity);                         // Not necessary for TCS
    Absorption(fDensity, fTX.Freq() );
    SkinDepth(fDensity, fTX.Freq() );             // Not necessary for TCS
    // Reflectance();
    Opacity(fAbsorption);
    TCS(fReflectance);

    //  Once we have the TCS,

    SegmentMaxCoords(fReflectance);


    SetSegments(nL);
    RunScatter();

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
  fCSLength = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  fCSRadius = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
// dL, dR should be given in cm
// fCSltot should be given in cm for now


  for (int i = 0; i < nL; i++){
    b = i* dL;
    for (int j = 0; j < nR; j++){
      a = j * dR;

      // What it should be ??
      // fCSLength[i][j] = (a - R/2)*cD + (b )*sD;
      fCSLength[i][j] = (a - R/2)*cD + (b - L/2)*sD + fCS.Ltot()/2;
      fCSRadius[i][j] = -(a - R/2)*sD + (b - L/2)*cD;

    }
  }
}


void Cascade1D::Density( const std::vector<std::vector<double>> &fCSLength,
                            const std::vector<std::vector<double>> &fCSRadius ){
  fDensity = std::vector<std::vector<double>> (fCSLength.size(),
              std::vector<double> (fCSLength[0].size(), 0.0));

  for (int i = 0; i < fDensity.size(); i++){
    for (int j = 0; j < fDensity[i].size(); j++){
      // Simple check to avoid computing values too far out from the cascade direction
      if ( abs(fCSRadius[i][j]) <= fCS.Rtot() ) {
        fDensity[i][j] = fCS.Density(rho_ice*fCSLength[i][j], fCSRadius[i][j]);
      }
    }
  }
}

void Cascade1D::Density2( const std::vector<std::vector<double>> &fCSLength,
                            const std::vector<std::vector<double>> &fCSRadius ){
  fDensity = std::vector<std::vector<double>> (fCSLength.size(),
              std::vector<double> (fCSLength[0].size(), 0.0));

  for (int i = 0; i < fDensity.size(); i++){
    for (int j = 0; j < fDensity[i].size(); j++){
      // Simple check to avoid computing values too far out from the cascade direction
      if ( abs(fCSRadius[i][j]) <= fCS.Rtot() ) {
        fDensity[i][j] = fCS.Density_beam(rho_ice*fCSLength[i][j], fCSRadius[i][j], 14.4, 1E9) ;
      }
    }
  }
}

std::vector<std::vector<double>> Density();

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
                      std::vector<double> (density[0].size(), 0));

  for (int i = 0; i < density.size(); i++){
    for (int j = 0; j < density[i].size(); j++){
        fPlasmaFrequency[i][j] = fCS.PlasmaFreq(density[i][j]);
    }
  }
}

void Cascade1D::Absorption(const std::vector<std::vector<double>> &density, const double & freq){
  fAbsorption = std::vector<std::vector<double>> (density.size(),
                  std::vector<double> (density[0].size(), 0));

  for (int i = 0; i < density.size(); i++){
    for (int j = 0; j < density[i].size(); j++){
      fAbsorption[i][j] = fCS.Absorption(density[i][j], freq);
        // fAbsorption[i][j] = fCS.Absorption2(density[i][j], freq);
    }
  }
}

void Cascade1D::SkinDepth(const std::vector<std::vector<double>> &density, const double & freq){
  fSkinDepth = std::vector<std::vector<double>> (density.size(),
                std::vector<double> (density[0].size(), 0));

  for (int i = 0; i < density.size(); i++){
    for (int j = 0; j < density[i].size(); j++){
      fSkinDepth[i][j] = fCS.SkinDepth(density[i][j], freq);
        // fSkinDepth[i][j] = fCS.SkinDepth2(density[i][j], freq);
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
	//OLD double reflectance, opacity;      // Unitless, tmp variables.
  double transparency, transmitivity;
  fTransparency = std::vector<std::vector<double>> (absorption.size(),
                  std::vector<double> (absorption[0].size(), 0));
  fOpacity = std::vector<std::vector<double>> (absorption.size(),
                  std::vector<double> (absorption[0].size(), 0));
  fReflectance = std::vector<std::vector<double>> (absorption.size(),
                  std::vector<double> (absorption[0].size(), 0));

	// Loop over the absorption matrix
  for (int i = 0; i < absorption.size(); i++){
		//OLD  reflectance = 0, opacity = 0;
    transmitivity = 1; transparency = 1;
    for (int j = 0; j < absorption[0].size(); j++){
      // NEW
      fTransparency[i][j] = transparency;
      fOpacity[i][j] = 1 - transparency;

/* For consistency and usefulness, the transparency and opacity are defined
___before___ reaching a certain layer.

      transparency after layer = transparency before layer * transmitivity
      power scattered per later = transparency before layer * reflectivity
      This was the old "reflectivity" definition.
*/
      transmitivity = exp(-1.0 * dR * absorption[i][j]);
      transparency *= transmitivity;
      assert(transparency > 0 && "Opacity larger than 1!");   // sanity check
      fReflectance[i][j] = 1 - transmitivity;


// OLD
    //   // The bin size [in cm] is the parallel dimension resolution.
    //   reflectance = (1-opacity)*(1-exp(-1.0 * dR * absorption[i][j]));
    //   opacity += reflectance;
		// 	assert(opacity < 1 && "Opacity larger than 1!");   // sanity check
    //   fReflectance[i][j] = reflectance;
    //   fOpacity[i][j] = opacity;
		}
	}
}

void Cascade1D::Opacity(const std::vector<std::vector<double>> &absorption){ Reflectance(absorption); }

/* Compute the TCS of the cascade from the slices.
For every cylindircal shell of thickness dR, we consider scattering over an area:

 Shell area = 2pi*base * height_segment = 2*pi* dL * dN;

NON constant layer size, dN is the normal direction, also over the radius.
Each radii dN actually  contributes to half a cylinder (dN goes from -r_max ro r_max):

  Half-shell area = pi* dN * dL;

This should work even if the frame of the shells is not perpendicular to the
cascade direction, as the shells would still form a cylinder.

Each shell is weighted by its reflectivity and the coherent scattering factor,
the number of electrons N_e in each shell.
*/
void Cascade1D::TCS(const std::vector<std::vector<double>> &reflectance){
  fTCS = std::vector<double> ( reflectance.size(), 0.0 );
  for (int i = 0; i < reflectance.size(); i++){
    for (int j = 0; j < reflectance[i].size(); j++){
      fTCS[i] += fTransparency[i][j] * pow(fDensity[i][j] * abs(fCSRadius[i][j]),2);
    }
    fTCS[i] *= 6.6524587*1E-25 * fDamping* pow( fTX.Freq(),2) *pow(dR * pi* dL,2); // [cm^2]
    // Thomson scattering cross section in cm^2
    // N_e = n_e*dV = n_e * dR * dA = n_e * dR* pi*dL*dN
  }
}

void Cascade1D::SegmentMaxCoords(const std::vector<std::vector<double>> &reflectance){
  double l,r;
  double nL = reflectance.size();
  fTXMaxLength = std::vector<double> (nL);
  fTXMaxRadius = std::vector<double> (nL);
  // i lives in L_TX, j lives in R_TX
  // (The integral of reflectance happens along j)
  for (int i = 0; i < nL; i++){
    // We want the index, not the value, therefore distance()
    int j = distance(reflectance[i].begin(),
            max_element(reflectance[i].begin(), reflectance[i].end() ) );
    fTXMaxLength[i] = i*dL;   // This is b before the rotation
    fTXMaxRadius[i] = j*dR;   // This is a before the rotation
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
std::vector<double> Cascade1D::TCS(){ return Scatter::TCS(); };

// std::vector<double> Cascade1D::radar_cs(){
//   if (fTCS.empty()) {set_radar_cs();}
//   return fTCS;
// }

// -----------------------------------------------------------------------------

/* Cylinder1D
This should be a close approximation to the ideal cascade. The radial profile
r(L) is not trivial for even simple geometries, so we have a cylinder instead.
A cylinder that has perfect reflectivity only scatters out of his outer shell.
The typical size for the cylinder is the moliere radius.

Outer cylinder shell area per segment = 2pibase* height = pi* r_moliere* dL
Opacity is 1;
*/

Cylinder1D::Cylinder1D(Antenna& tx, Antenna& rx, Cascade& cs): Scatter(tx, rx, cs){
  fTCS = std::vector<double>(nL, dL * r_moliere * pi); // [cm^2]
  SetSegments( nL );
  RunScatter();
}

// =============================================================================
