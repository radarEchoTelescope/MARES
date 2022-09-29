#include "cascade1D.hh"

Cascade1D::Cascade1D(Antenna& tx, Antenna& rx, Cascade& cs,
                    const double deltaL, const double deltaR, const double deltaN,
                    const double sampling):
  Scatter1D(tx, rx, sampling), Cascade(cs){

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

// Third, compute the transparency and the TCS in one go
  TCS(fCSRadius, fDensity);

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
  fDensity = std::vector<std::vector<double>>  (nL, std::vector<double> (nR, 0));

// Fill the length and radial values in the transmitter frame. 
double r, l, l_tmp, r_tmp;
  for (int i = 0; i < nL; i++){
    l = i* dL;
    for (int j = 0; j < nR; j++){
      r = j * dR;
      l_tmp = (r - R/2)*cD + (l - L/2)*sD + fLtot/2;
      r_tmp = -(r - R/2)*sD + (l - L/2)*cD;
      // Simple check to avoid computing values too far out from the cascade direction
      if ( abs(r_tmp) <= fRtot ) {
         fDensity[i][j] = Cascade::Density(rho_ice*l_tmp, r_tmp, dR); 
      }
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
    fTCS[i] = pow(tmp_tcs,2) * pi * dL * dR * thomson * fDamping;
  }
}

// Accesors
std::vector<double> Cascade1D::TCS(){return fTCS; }
std::vector<std::vector<double>> Cascade1D::Length() {return fCSLength;}
std::vector<std::vector<double>> Cascade1D::Radius() {return fCSRadius;}