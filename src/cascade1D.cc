#include "cascade1D.hh"

Cascade1D::Cascade1D(Antenna& tx, Antenna& rx, Cascade& cs,
                    const double deltaL, const double deltaR, const double deltaN,
                    const double sampling):
  Scatter1D(tx, rx, sampling), Cascade(cs), dL(deltaL), dR(deltaR), dN(deltaN){

/* First: Set the antennas directions, module and dot product with cs */
  fTX.SetDirection( cs.Pos() );
  fRX.SetDirection( cs.Pos() );
  fTX.SetAngle( cs.Dir() );
  fRX.SetAngle( cs.Dir() );

/* Second: Build the TX frame. Compute the density matrix */
  SetTXFrame();

// From Density to TCS.
    // Density = Density(Length, Radius, dR);
    // Transparency = Transparency(Attenuation(Density), dR);
    // TCS(SRadius, Density, Transparency, damping, ...);

    // Absorption(fDensity, fTX.Freq());

// Third, compute the transparency and the TCS in one go
  SetTCS();

// Lastly, start filling the scattering centres
  ScatterPoint p; 
  for(int i = 0 ; i < nL; i++ ){
      p = {};
      p.ID = i;
      p.TCS = fTCS[i];
      // p.Size = {dL, dR, dN};
      AddPoint(p);
  }

  }

void Cascade1D::SetInDirection(double rand_seed){
  Scatter1D::SetInDirection( Cascade::fPosition, Cascade::fDirection,
                Cascade::fLtot, rand_seed);
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
  // We can avoid computing the same values thousands of times later.
  cD = fTX.Dot();
  sD = sqrt(1 - pow(cD,2));

  R = fLtot*abs(cD) + 2*fRtot*abs(sD);
  L = fLtot*abs(sD) + 2*fRtot*abs(cD);

  nR = (int) ceil(R  / dR);
  nL = (int) ceil(L / dL);
  std::cout << dL << std::endl;

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
          Iwr[i][j] = NKG::intwiv(abs(r_tmp),dR);
          fNe[i][j] = Iwr[i][j]*Cascade::Ne(rho_ice*l_tmp, Cascade::fEnergy, Cascade::fNp)*dL;
          fDensity[i][j] = fNe[i][j]/ (pi*dL*(pow(dR,2) + 2*abs(r_tmp)*dR));
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
void Cascade1D::SetTCS(){

  double transparency, transmitivity, tmp_tcs;
  // This is the damping factor, already squared. 
  double fDamping = 1.0 / sqrt( pow(fTX.Freq(), 4) + pow(fTX.Freq()*f_coll,2) );
  fTCS = std::vector<double> ( fDensity.size(), 0.0 );
  for (int i = 0; i < fDensity.size(); i++){
    transmitivity = 1; transparency = 1, tmp_tcs = 0;
    for (int j = 0; j < fDensity[i].size(); j++){
      if(fDensity[i][j] == 0){continue;}
      transmitivity = exp(-1.0 * dR * Absorption( fDensity[i][j],fTX.Freq() ) );
      transparency *= transmitivity;
      // Very important sanity check
      // Reveals is something goes wrong with the density calc.
      assert(transparency > 0 && "Opacity larger than 1!");
      
      // We sum into one TCS the contributions from the TCS from the 
      // separate shells.
      tmp_tcs += pow(fNe[i][j],2) * transparency;
      // CAREFUL! Sum(n)^2 != Sum(n^2)!!!
    }
    
    fTCS[i] = tmp_tcs*fDamping * thomson * 1.5 /dR;
    //1.5 is (the gain of) the Herzian dipole factor
    // dR is necessary to normaise here the number of steps/iterations that we do in this loop. 
  }
}



// Accesors
std::vector<double> Cascade1D::TCS(){return fTCS; }
std::vector<std::vector<double>> Cascade1D::Length() {return fCSLength;}
std::vector<std::vector<double>> Cascade1D::Radius() {return fCSRadius;}