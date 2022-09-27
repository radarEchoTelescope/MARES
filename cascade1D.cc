#include "cascade1D.hh"

Cascade1D::Cascade1D(Antenna& tx, Antenna& rx, Cascade& cs,
                    const double deltaL, const double deltaR, const double deltaN,
                    const double lifetime, const double sampling):
    Scatter(tx, rx, lifetime, sampling), Cascade(cs){

/* First: Set the antennas directions, module and dot product with cs */
    fTX.SetDirection( this->fPosition );
    fRX.SetDirection( this->fPosition );

/* Second: Build the TX frame. */

    fTX.SetAngle( this->fDirection );
    fRX.SetAngle( this->fDirection );
    SetTXFrame();

    TCS(fCSRadius, fDensity);

  // From Density to TCS.
    // fDensity = Density(fCSLength, fCSRadius, dR);
    // Transparency(fDensity);
    // TCS(fCSRadius, fDensity, fTransparency);

    //  Once we have the TCS,
    // Basic option
    // SetWRTCSDirection();
    // This is a call to Scatter::SetInLine();

    // or

    // SegmentMaxCoords(fReflectance);
    // SetWRTCSPosition(fTXMaxLength, fTXMaxRadius);

    SetInConstIce();
    // SetInBeam();
    // SetWithIRT();

    RunScatter();


    // Set times
    // (Retarded) time where the scattered radio signal by the segment is produced.
    // production = birth - fRTX/c_ice;

    // (Advanced) time where the scattered signal by the segments arrives in the receiver.
    // arrival = birth + Rr/c_ice;

    // T0 = 0 by definition when the cascade begins (first element = head).
    // p.ArrivalTime = p.L/c_vac + p.RRX/c_ice;
    // time where the segment starts scattering.

  }

/* Make 2D-array of density profile in the TX frame
Transmitter frame (a, b) goes from (0 -> A , 0 -> B)
Cascade frame (L,R) goes to (0 -> L, -r/2 -> r/2)
The rotation between the frames needs to happen at cascade's center.
The formula below is translation -> rotation -> translation back.
NOTE: l,r in the function are meant to be positions in the cascade frame (after rotation)
but still move along TX frame.
*/
void Cascade1D::SetTXFrame(){
  // cD = cosine delta, sD, sine delta, the angle between TX_dir and CS_dir.
  double cD, sD;
  cD = fTX.Dot();
  sD = sqrt(1 - pow(cD,2));
  // We can avoid computing the same values thousands of times.

// R,L are in cm since everything in the cascade should be in cm.
  R = fLtot*abs(cD) + 2*fRtot*abs(sD);
  L = fLtot*abs(sD) + 2*fRtot*abs(cD);

  nR = (int) ceil(R  / dR);
  nL = (int) ceil(L / dL);

// Sanity check! Your sections are not unphysical due to lifetime constraint.
  assert(fLtot / nL  <= c_vac*tau && "Cascade segments are too large!");

// This can be further checked against the probing wavelength.
// Let have a WARNING if we are not in the Fraunhofer regime.

  fCSLength = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  fCSRadius = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  fDensity = std::vector<std::vector<double>> (fCSLength.size(),
              std::vector<double> (fCSLength[0].size(), 0));

// Now we find the length and radial values in the transmitter frame, and
double r, l, l_tmp, r_tmp;
  for (int i = 0; i < nL; i++){
    l = i* dL;
    for (int j = 0; j < nR; j++){
      r = j * dR;
      l_tmp = (r - R/2)*cD + (l - L/2)*sD + fLtot/2;
      r_tmp = -(r - R/2)*sD + (l - L/2)*cD;
      // Simple check to avoid computing values too far out from the cascade direction
      if ( abs(r_tmp) <= fRtot ) { fDensity[i][j] = Cascade::Density(rho_ice*l_tmp, r_tmp, dR); }
      fCSLength[i][j] = l_tmp;
      fCSRadius[i][j] = r_tmp;
    }
  }
}


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

// N_e = n_e*dV = n_e * dR * dA = n_e * dR* pi*dL*dN
// N_e[i] += fDensity[i][j] * pi * dL * dR * abs(fCSRadius[i][j]);
*/
void Cascade1D::TCS(const std::vector<std::vector<double>> &radius,
                    const std::vector<std::vector<double>> &density){

  double transparency, transmitivity, tmp_tcs;
  double fDamping = 1.0 / sqrt( pow(fTX.Freq(), 4) + pow(fTX.Freq()*f_coll,2) );
  fTCS = std::vector<double> ( density.size(), 0.0 );

  for (int i = 0; i < density.size(); i++){
    transmitivity = 1; transparency = 1, tmp_tcs = 0;
    for (int j = 0; j < density[i].size(); j++){
      if(density[i][j] == 0){continue;}

      transmitivity = exp(-1.0 * dR * Absorption( density[i][j],fTX.Freq() ) );
      transparency *= transmitivity;
      assert(transparency > 0 && "Opacity larger than 1!");   // sanity check

      tmp_tcs += density[i][j] * abs(radius[i][j]) * transparency;
    }
    fTCS[i] = pow(tmp_tcs,2) * pi * dL * dR * thomson * fDamping; // [cm^2]
    // Thomson scattering cross section in cm^2 converts from number of coherent electrons to area
  }
}

// Set the position and arrival time of your segments

void Cascade1D::SetInLine(std::vector<double> vertex,
                        std::vector<double> direction,
                        double length, double rand_seed){

  ScatterPoint p;
  double dSeg = length/nP;
  // This does not work
  // if(rand_seed){ RN_uniform rand_line(-0.5, 0.5, rand_seed);}
  RN_uniform rand_line(-0.5, 0.5, rand_seed);

  for(int i = 0 ; i < nP; i++){
    p = fPoints[i];

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


void Cascade1D::SetWRTCSPosition(std::vector<double> & CSLength,
  const std::vector<double>& CSRadius){
  ScatterPoint seg;
  // fPoints = std::vector<ScatterPoint> (nL);
  // The length of CSLength and CS Radius must be nL, as that is the number of
  // segments to be placed.

  // This is not used here, but we could insert randomisation too.

  // // Giving a non-uniform position in the cascade gets rid of artifacts in the FFT.
  // RN_uniform rand_line(-0.5, 0.5, 42);          // Same seed for debugging.
  // // RN_uniform rand_line(-0.5, 0.5, time(0));  // Different seed for random, independent runs.

  std::vector<double> R_TX_dir{0,0,0}, L_TX_dir{0,0,0},
                      seg_pos{0,0,0}, tx_seg{0,0,0}, rx_seg{0,0,0};

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

// This only breaks for R_TX_dir == this->fDirection, exactly, we can guard against it easily
// A deviation of 10^-20 in the direction is well within numerical errors, and avoids
// making the cross product exactly zero.

  if( ( R_TX_dir == normalize( this->fDirection ) ) ) {
    R_TX_dir[0] += 1E-20;
    R_TX_dir[1] += 1E-20;
    R_TX_dir[2] += 1E-20;
  }

  L_TX_dir = cross_product( normalize( cross_product(R_TX_dir, normalize(this->fDirection) ) ), R_TX_dir);
  L_TX_dir = normalize(L_TX_dir);

  // To separate between the cascade pointing towards TX and away from TX
  double s = sgn( projection(R_TX_dir, this->fDirection) ); // cos(delta)
  if(s < 0){
    std::reverse( CSLength.begin(), CSLength.end() );
  }

  // Segment loop
  for(int i = 0 ; i < nL; i++){
    seg = fPoints[i];

    // Set position
    // [m] distance from the shower head (starting point)

    seg.Position[0]  = CSRadius[i]*R_TX_dir[0]*s + CSLength[i]*L_TX_dir[0];
    seg.Position[1]  = CSRadius[i]*R_TX_dir[1]*s + CSLength[i]*L_TX_dir[1];
    seg.Position[2]  = CSRadius[i]*R_TX_dir[2]*s + CSLength[i]*L_TX_dir[2];

    // CSLength, CSRadius should have been computed in cm
    seg.Position[0] /= 100.;
    seg.Position[1] /= 100.;
    seg.Position[2] /= 100.;

    // fCSPos is given in m
    seg.Position[0] += this->fPosition[0];
    seg.Position[1] += this->fPosition[1];
    seg.Position[2] += this->fPosition[2];

    seg.L = distance(seg.Position[0], seg.Position[1], seg.Position[2],
            fPoints[0].Position[0],  fPoints[0].Position[1], fPoints[0].Position[2]);
    }
}


// Accesors

std::vector<double> Cascade1D::TCS(){return fTCS; }
std::vector<std::vector<double>> Cascade1D::Length() {return fCSLength;}
std::vector<std::vector<double>> Cascade1D::Radius() {return fCSRadius;}



// std::vector<double> Cascade1D::TCS(){Scatter::TCS(); };


// double  Cascade1D::Rwaist(){ return fRwaist; }
// std::vector<double> Cascade1D::Rcrit(){ return fRcrit; }

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


// -----------------------------------------------------------------------------

// /* Line1D is Dieder's model.
// This should be the closest representation to the ideal thin-wire solution
// (there is an analytical solution that does not contain retardation effects)
// */
// Line1D::Line1D(Antenna& tx, Antenna& rx, Cascade& cs): Scatter(tx, rx, cs){
//
//   // std::cout << 2* 1.23 * fCS.Rtot()*fCS.Rtot()*2.0 <<'\t' << '\n';
//
//   fTCS = std::vector<double>(nL, dL*dN*2.0); // [cm^2]
//   //The line's segements are is 2a*a, with perfect reflectivity.
//   // A line is expected to have Different polarization than a cascade
//   SetSegments( nL );
//   RunScatter();
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

// Cylinder1D::Cylinder1D(Antenna& tx, Antenna& rx, Cascade& cs): Scatter(tx, rx, cs){
//   fTCS = std::vector<double>(nL, dL * r_moliere * pi); // [cm^2]
//   SetSegments( nL );
//   RunScatter();
// }

// =============================================================================
