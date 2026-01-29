#include "cascade1D.hh"
#include <sstream>

Cascade1D::Cascade1D(Antenna& tx, Antenna& rx, Cascade& cs,
                    const double deltaL, const double deltaR,
                    const double sampling):
  Scatter(tx, rx, sampling), Cascade(cs), dL(deltaL), dR(deltaR){

/* First: Set the antennas directions, module and dot product with cs */
  fTX.SetDirection( cs.Pos() );
  fRX.SetDirection( cs.Pos() );

  fTX.SetAngle( cs.Dir() );
  fRX.SetAngle( cs.Dir() );

/* Second: Build the TX frame. Compute the density matrix */
// Choose type of cascade to run:
if (cascade_type == "neutrino"){
  SetTXFrame();     // Default cascade generation method (neutrino NKG)
} else if (cascade_type == "cosmic_ray"){
  SetTXFrameCR();  // Cosmic ray cascades (binned NKG)
} else if (cascade_type == "faerie"){
  SetTXFrameFaerie();  // Read in Faerie profile (in progress, not well tested)
} else {
  std::cerr<< "Cascade generation type not recognised. Using a neutrino NKG cascade."<<std::endl;
  SetTXFrame();
}
// From Density to TCS.
    // Density = Density(Length, Radius, dR);
    // Transparency = Transparency(Attenuation(Density), dR);
    // TCS(SRadius, Density, Transparency, damping, ...);


// UNCOMMENT HERE AS NEEDED FOR NOW
    // PlasmaFreq(fDensity);
    // Absorption(fDensity, fTX.Freq());
    // SkinDepth(fDensity, fTX.Freq());
    // Transparency(fDensity,fTX.Freq(), dR);

// Third, compute the transparency and the TCS in one go
  SetTCS();

// Lastly, start filling the scattering centres
  ScatterPoint p; 
  for(int i = 0 ; i < nL; i++ ){
      p = {};
      p.ID = i;
      p.TCS = fTCS[i];
      AddPoint(p);
  }

}

void Cascade1D::UpdateRX(const Antenna& new_RX){
  Scatter::UpdateRX(new_RX);
  fRX.SetDirection( Cascade::fPosition );
  fRX.SetAngle( Cascade::fDirection);
}

void Cascade1D::SetInDirection(double rand_seed){
  Scatter::SetInDirection( Cascade::fPosition, Cascade::fDirection,
                Cascade::fLtot, rand_seed);
}


/* Make 2D-array of density profile in the TX frame
Transmitter frame (a, b) goes from (0 -> A , 0 -> B)
Cascade frame (L,R) goes to (0 -> L, -r/2 -> r/2)
The rotation between the frames needs to happen at cascade's center.
The formula below is translation -> rotation -> translation back.
NOTE: l,r in the function are positions in the cascade frame (after rotation)
but move along the TX frame (are defined there.).
*/
void Cascade1D::SetTXFrame(){
  // We can avoid computing the same values thousands of times later.
  cD = fTX.Dot();
  sD = sqrt(1 - pow(cD,2));

  // cout<<sD<<endl;    //debug cout
  // cout<<cD<<endl;    //debug cout

  R = fLtot*abs(cD) + 2*fRtot*abs(sD);
  L = fLtot*abs(sD) + 2*fRtot*abs(cD);

  nR = (int) ceil(R  / dR);
  nL = (int) ceil(L / dL);

// Sanity check! Your sections are not unphysical due to lifetime constraint.
  // assert(fLtot / nL  <= c_vac*tau && "Cascade segments are too large!");

// This can be further checked against the probing wavelength.
// Let have a WARNING if we are not in the Fraunhofer regime.

  fCSLength = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  fCSRadius = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  fNe       = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  fDensity  = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  std::vector<std::vector<double>> Iwr  = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));

// Fill the length and radial values in the transmitter frame. 
double r, l, l_tmp, r_tmp, norm;

  for (int i = 0; i < nL; i++){
    l = i* dL;
    norm = 0.0;
    for (int j = 0; j < nR; j++){
      r = j * dR;
      l_tmp = (r - R/2)*cD + (l - L/2)*sD + fLtot/2;
      r_tmp = -(r - R/2)*sD + (l - L/2)*cD;

      // Simple check to avoid computing values too far out from the cascade direction
      if ( abs(r_tmp) <= fRtot ) {
          Iwr[i][j] = NKG::intwiv(abs(r_tmp),dR,shower_age);
          fNe[i][j] = Iwr[i][j]*Cascade::Ne(rho_ice*l_tmp, Cascade::fEnergy, Cascade::fNp)*dL;
          fDensity[i][j] = fNe[i][j]/ (pi*dL*(pow(dR,2) + 2*abs(r_tmp)*dR));
          // or,
          // fDensity[i][i] = Cascade::Density(rho_ice*l_tmp, r_tmp, dR, Cascade::fEnergy, Cascade::fNp  );
      }
      fCSLength[i][j] = l_tmp;
      fCSRadius[i][j] = r_tmp;
    }
    // Normalisation in R!
    norm = std::accumulate(std::begin(Iwr[i]), std::end(Iwr[i]), 0.0);
    std::transform(fNe[i].begin(), fNe[i].end(), fNe[i].begin(), [norm](double &c){ return c/norm; });
    std::transform(fDensity[i].begin(), fDensity[i].end(), fDensity[i].begin(), [norm](double &c){ return c/norm; });

  }
}

// New SetTXFrame function with variable density:
void Cascade1D::SetTXFrameCR(){

  // Read in the Ep and Np text files for the cascade generation
  ifstream EpFile  (epfilepath);
  ifstream NpFile  (npfilepath);

  fEnergyVals  = read_1D_array( EpFile ) * GeV ; // convert from GeV -> MeV (MARES energy unit.)
  fNumVals     = read_1D_array( NpFile ) ;
  int nEl      = fEnergyVals.size(); // number of bins in the file

  // We can avoid computing the same values thousands of times later.
  cD = fTX.Dot();
  sD = sqrt(1 - pow(cD,2));

  R = fLtot*abs(cD) + 2*fRtot*abs(sD); // = B     (the cascade in LOS frame, rectangle of A and B)
  L = fLtot*abs(sD) + 2*fRtot*abs(cD); // = A 

  nR = (int) ceil(R  / dR);
  nL = (int) ceil(L / dL);

  // cout<<fEnergyVals[0]*GeV<<endl;    //debug cout
  // cout<<fNumVals[0]<<endl;           //debug cout
  // cout<<E_c<<endl;                   //debug cout
  // cout<<fDirection<<endl;            //debug cout
  // cout<<normalize(fTX.Dir())<<endl;  //debug cout
  // cout<<nEl<<endl;                   //debug cout
  // cout<<"sD: "<<sD<<endl;            //debug cout
  // cout<<"cD: "<<cD<<endl;            //debug cout
  // cout<<"L: "<<L<<endl;              //debug cout
  // cout<<"R: "<<R<<endl;              //debug cout
  // cout<<nL<<endl;                    //debug cout
  // cout<<nR<<"\n"<<endl;              //debug cout
  // cout<<fDirection<<endl;            //debug cout
  // cout<<fPosition<<endl;             //debug cout
  // cout<<rho_ice<<endl;               //debug cout

// Sanity check! Your sections are not unphysical due to lifetime constraint.
// assert(fLtot / nL  <= c_vac*tau && "Cascade segments are too large!");

// This can be further checked against the probing wavelength.
// Let have a WARNING if we are not in the Fraunhofer regime.

  fCSLength = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  fCSRadius = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  fNe       = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  fDensity  = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  std::vector<std::vector<double>> Iwr  = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));

  // basis vector in CS frame, used to find the depth of each point
  e_R = cross_product(fDirection,normalize(cross_product(normalize(fTX.Dir()),fDirection))) ;
  e_R = normalize(e_R);

// Fill the length and radial values in the transmitter frame. 
double r, l, l_tmp, r_tmp, norm;

// new variables/arrays needed for cosmic ray method.
double z_depth, rho_ice_d, NEval, NumParticles,dl_;
double l_tmp_prev = 0.0;
double columnDepth[nR];
double lvals[nR];
for (int i = 0; i < nR; i++){
  lvals[i] = 0.;
  columnDepth[i] = 0.;
}

double X1 = 18.0 *(g/pow(cm,2)) ; // correction factor for the NKG formula. Currently hardcoded, shouldn't be changed.
double ELowerLim = 0.2*GeV ;  // energy lower limit for NKG. Currently hardcoded, shouldn't be changed.


  // start the cascade generation.
  for (int i = 0; i < nL; i++){
    l = i* dL;
    norm = 0.0;

    for (int j = 0; j < nR; j++){
      r = j * dR;
      // obtaining the CS frame ([0,L][-R/2,R/2]) coordinates of the points as defined in the LOS frame ([0,nA][0,nB])
      l_tmp = (r - R/2)*cD + (l - L/2)*sD + fLtot/2;  
      r_tmp = -(r - R/2)*sD + (l - L/2)*cD;

      // Calculate the depth at each successive l and r value in the lab frame
      z_depth   = abs((fPosition + fDirection*l_tmp + e_R*r_tmp )[2]) ;
      // Return the ice density at the given depth
      rho_ice_d = NKG::rhoGreenland(z_depth) ;

      // calculating the column density for variable density ice (cumulative sum)
      if (l_tmp > 0.){ // l_tmp may be negative at the start, as the rectangles are misaligned,
                        // so we make sure to skip those values
        dl_ = l_tmp - lvals[j]; //obtaining bin width (length travelled by that bin) in the cascade frame 
        columnDepth[j] += rho_ice_d * dl_; // accumulate column depth with new bin added
        lvals[j] = l_tmp; // save new l val to use in next row
      }
    
      // Simple check to avoid computing values too far out from the cascade direction
      if ( abs(r_tmp) <= fRtot ) {
          Iwr[i][j] = NKG::intwiv(abs(r_tmp),dR,shower_age);

          NEval = 0.0;
          for (int k = 0; k < nEl; k++){
            if (fEnergyVals[k] > ELowerLim){
              // adjust the number of particles with the early correction formula
              NumParticles = fNumVals[k] * NKG::removeEarlyN( columnDepth[j]/(g/pow(cm,2)), 9e-2, 0.407 );
              // accumulate the NKG number for each energy bin
              NEval += Cascade::Ne_vd(columnDepth[j]+(X1), fEnergyVals[k], NumParticles, rho_ice_d);
            }
          }

          fNe[i][j] = Iwr[i][j]*NEval*dL;
          fDensity[i][j] = fNe[i][j]/ (pi*dL*(pow(dR,2) + 2*abs(r_tmp)*dR));

      }
      fCSLength[i][j] = l_tmp;
      fCSRadius[i][j] = r_tmp;
    }
    // Normalisation in R!
    norm = std::accumulate(std::begin(Iwr[i]), std::end(Iwr[i]), 0.0);
    std::transform(fNe[i].begin(), fNe[i].end(), fNe[i].begin(), [norm](double &c){ return c/norm; });
    std::transform(fDensity[i].begin(), fDensity[i].end(), fDensity[i].begin(), [norm](double &c){ return c/norm; });

  }
}

// Another new SetTXFrame function to read in Faerie profiles
// Here, the density and number cascades are already generated by Faerie, and we only need to read them in. 
// This requires that the Faerie data is interpolated to match MARES binning for each geometry.
// CURRENTLY ONLY WORKS FOR SIDE-ON GEOMETRY, WHERE THE CASCADE AND LOS FRAMES ARE ALIGNED - needs further testing.
void Cascade1D::SetTXFrameFaerie(){

  cD = fTX.Dot();
  sD = sqrt(1 - pow(cD,2));

  R = fLtot*abs(cD) + 2*fRtot*abs(sD);
  L = fLtot*abs(sD) + 2*fRtot*abs(cD);

  nR = (int) ceil(R  / dR);
  nL = (int) ceil(L / dL);


// Sanity check! Your sections are not unphysical due to lifetime constraint.
  // assert(fLtot / nL  <= c_vac*tau && "Cascade segments are too large!");

// This can be further checked against the probing wavelength.
// Let have a WARNING if we are not in the Fraunhofer regime.

  fCSLength = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  fCSRadius = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  fNe       = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  fDensity  = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));
  std::vector<std::vector<double>> Iwr  = std::vector<std::vector<double>> (nL, std::vector<double> (nR, 0));


  // Fill the length and radial values in the transmitter frame. 
  double r, l, l_tmp, r_tmp, norm;

  // Faerie file to read in
  std::ifstream file(faeriefilepath); // density file hardcoded, to be removed.
  std::ifstream densityfile("/Users/ishaloudon/local_code/inIceCR/c5e16.5z0a0v0_2dprofile_dens.txt");
  string line  = "";
  string line2 = "";

  int i = 0;
  while(getline(file, line)){  // going line by line (same as looping over i in standard cascade function)
    l = i * dL;

    // Splitting the line into parts by spaces
    std::stringstream ss(line);
    string entry = "";

    getline(densityfile, line2);
    std::stringstream dd(line2);
    string entry2 = "";

    int j = 0;
    while(std::getline(ss, entry, ' ')){  // going across each line (looping over j)
      double val = (double)std::stod(entry);

      std::getline(dd, entry2, ' ');
      double den = (double)std::stod(entry2);

      r = j * dR;
      l_tmp = (r - R/2)*cD + (l - L/2)*sD + fLtot/2;
      r_tmp = -(r - R/2)*sD + (l - L/2)*cD;

      fNe[i][j] = val;
      fDensity[i][j] = den;

      fCSLength[i][j] = l_tmp;
      fCSRadius[i][j] = r_tmp;

      j++;
      }

    i++;
    }

    file.close();

    // }

    // ToDo before using - determine effects of normalisation in this mode! Currently disabled.

    // norm = std::accumulate(std::begin(Iwr[i]), std::end(Iwr[i]), 0.0);
    // std::transform(fNe[i].begin(), fNe[i].end(), fNe[i].begin(), [norm](double &c){ return c/norm; });
    // std::transform(fDensity[i].begin(), fDensity[i].end(), fDensity[i].begin(), [norm](double &c){ return c/norm; });

  // }

}

/* Compute the TCS of the cascade from the slices.
For every cylindrical shell of thickness dR, we consider scattering over an area:

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
void Cascade1D::SetTCS(){

  double transparency, transmissivity, tmp_tcs;
  // This is the damping factor, \omega*W. 
  double fDamping = 1.0 / (1 + pow(f_coll/fTX.AngularFreq(), 2) );

  fTCS = std::vector<double> ( fDensity.size(), 0.0 );
  for (int i = 0; i < fDensity.size(); i++){
    transmissivity = 1; transparency = 1, tmp_tcs = 0;
    for (int j = 0; j < fDensity[i].size(); j++){
      if(fDensity[i][j] == 0){continue;}
      transmissivity = exp(-1.0 * dR * Absorption( fDensity[i][j],fTX.Freq() ) );
      transparency *= transmissivity;
      // Very important sanity check
      // Reveals is something goes wrong with the density calc.
      assert(transparency > 0 && "Opacity larger than 1!");
      
      // We sum into one TCS the electrons from the TCS from the separate shells.
      tmp_tcs += fNe[i][j];
    }
    
    fTCS[i] = pow(tmp_tcs,2)*transparency*fDamping * thomson *1.5  /dR;
    //1.5 is (the gain of) the Herzian dipole factor
    // dR is necessary to normaise here the number of steps/iterations that we do in this loop. 
  }
}


void Cascade1D::SetInCSPlane(std::vector<double> vertex,
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

  SetInPlane(vertex, L_TX_dir, R_TX_dir*s, l_vals, r_vals);

}

/*  For each segment in the L_TX direction, we can find the coordinates (r_max,l_max)
  of the point where the reflectivity is the highest (it should correlate with the
  peak density for the segment, too), and place the segment at that point.
*/
void Cascade1D::SetInMax(std::vector<double> vertex, std::vector<double> direction, 
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

void Cascade1D::SetInMaxReflectivity(){
  // Compute the reflectivty matrix first.
  Cascade::Transparency(Cascade::fDensity, Scatter::fTX.Freq(), dR );
  // Place the segments along the values. 
  SetInMax( Cascade::fPosition, Cascade::fDirection,
                dL, dR, Cascade::fReflectance);
}


// Accessors
std::vector<double> Cascade1D::TCS(){return fTCS; }
std::vector<std::vector<double>> Cascade1D::Length() {return fCSLength;}
std::vector<std::vector<double>> Cascade1D::Radius() {return fCSRadius;}
