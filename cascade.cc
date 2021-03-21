#include "cascade.hh"

Cascade::Cascade(int event, double cenergy, double xpos, double ypos, double zpos,
								 double czenith, double cazimuth,  double nenergy,
								 double nzenith, double nazimuth, double oneweight):
				 			   fEvent(event),
								 fEnergy(cenergy),
								 fPosition{xpos, ypos, zpos},
								 fSphericalAngles{czenith,cazimuth},
								 fNeutrino{nenergy, nzenith, nazimuth,oneweight} {

	 // Sanity check! Your sections are not unphysical due to lifetime constraint.
	 assert(fLdiv <= cvac_cm*tau && "Cascade resolution is too large!");
	 // This should be further checked against the probing wavelength.

	fDirection[0] = sin(fSphericalAngles[0])*cos(fSphericalAngles[1]);
	fDirection[1] = sin(fSphericalAngles[0])*sin(fSphericalAngles[1]);
	fDirection[2] = cos(fSphericalAngles[0]);

	/* X_tot 	= 4* number of divisions for shower_max * interaction length
				  = 4*(log(fEnergy/E_c)/log(2) )*X_int;
		 			= 4*(log(fEnergy/E_c)) * X_0
 */
	fXtot = 4 * log(12.72 * fEnergy) * X_0;			// [g/cm^2]
	fLtot = fXtot/rho_ice; 															// [cm]
	fRtot = 5*r_moliere;//*log(fEnergy/1E9);						// [cm]
	// THIS WORKS FOR E 1E9, switch to E MIN 10^6 GeV AND VERIFY FOR ALL ENERGIES.

	fLbins = (int) ceil( fLtot/fLdiv );
	fXbins = (int) ceil( fXtot/fXdiv );
	fRbins = (int) ceil( fRtot/fRdiv );

}

// Cascade's Methods

double Cascade::Density(double X, double r){
	double dens, delta_r = 0.05;
	// s = ShowerAge(X,fEnergy);*step
	if (r < 0) {r = -r;}
	if (X < 0) {X = 0;}
	dens = Ne(X,fEnergy) * intwiv(r, delta_r) / (pi*(pow(r + delta_r,2) - pow(r,2)));
	assert(dens >= 0 && "Negative density value"); // Sanity check
	return dens;
}

double Cascade::PlasmaFreq(const double &dens){ return (8980 * sqrt(memp) * sqrt(dens)); }

double Cascade::Absorption(const double &dens, const double &freq_obs){
	double w, a, b, q;     // Some temp varables.
  double fplasma = PlasmaFreq(dens);

  // Exact solution from dispersion relation with collisions.
  if (f_coll != 0){
    w = pow(fplasma,2)/( pow(freq_obs,2) + pow(f_coll,2) );
    a = 1 - w;
    b = (f_coll / freq_obs) * w;

    // double p = (fAngularFreq/cvac_cm)*np.sqrt((np.sqrt(a**2 + b**2) + a) /2 )
    q = (freq_obs/cvac_cm)*sqrt((sqrt(pow(a,2) + pow(b,2)) - a) /2 ); // [1/cm]
  } else {
  // Collisionless model
    fplasma > freq_obs ? q = fplasma/cice_cm : q = 0;               // [1/cm]
  }
  return q;
}

double Cascade::SkinDepth(const double &dens, const double &freq_obs){
	return 1/Absorption(dens, freq_obs);
}

/* Compute the 2D density profile in the cascade frame */
void Cascade::SetDensity(){
	double l, r;
	fDensity = std::vector<std::vector<double>> (fLbins, std::vector<double> (fRbins, 0));

	for (int i = 0; i < fLbins; i++){
		l = i*fLdiv;
		for (int j = 0; j < fRbins; j++){
			r = j*fRdiv;
			fDensity[i][j] = Density(l*rho_ice,r);	// dens(X,r) takes columm density.
		}
	}
}
//
// /* Computes crital radius, the radial values of the line where the plasma frequency
// equals the detection frequency (overdense line).
// The density associated with the plasma frequency is the critical density [(#e-) cm^-3]
// If the local density is higher than that then the medium is overdense.
//
// The shower waist is by definition the maximum radius for the same plasma density.
// */
// void Cascade::SetRcrit(std::vector<std::vector<double>> &density, const double & freq){
// 	double r, dens_crit = memp*pow(freq/8980,2);
// 	fRcrit = std::vector<double>(density.size(), 0);
//
// 	// Loop over density matrix rows and then columms.
// 	for (int i = 0; i < density.size(); i++){
// 		r = 0;
// 		for (int j = 0; j < density[0].size(); j++){
// 			// If overdense, remember the layer.
// 			if (density[i][j] > dens_crit) { r = j;}   // Layers
// 		}
// 		fRcrit[i] = r*fRdiv; // [cm]
// 	}
//
// 	fRwaist = *max_element(fRcrit.begin(), fRcrit.end());
// }
// ----------------------------------------------------------------------------
// Accesors
int  Cascade::Evt()	const {return fEvent;}
double  Cascade::Energy() const {return fEnergy;}
std::vector<double> Cascade::Sph() const {return fSphericalAngles;};
std::vector<double> Cascade::Pos() const {return fPosition;}
std::vector<double> Cascade::Dir() const {return fDirection;}
double* Cascade::Parent() {return fNeutrino;}

double Cascade::Ldiv() const {return fLdiv;}
double Cascade::Xdiv() const {return fXdiv;}
double Cascade::Rdiv() const {return fRdiv;}

double Cascade::Ltot() const {return fLtot;}
double Cascade::Xtot() const {return fXtot;}
double Cascade::Rtot() const {return fRtot;}

double Cascade::Lbins() const {return fLbins;}
double Cascade::Xbins() const {return fXbins;}
double Cascade::Rbins() const {return fRbins;}

std::vector<std::vector<double>> Cascade::Density()  {
  if ( fDensity.empty() ) {SetDensity();}
  return fDensity;
}

// -----------------------------------------------------------------------------
std::vector<Cascade> load_cascade_file(const std::string& cs_filepath){
  std::ifstream cs_file(cs_filepath);

  if (!cs_file.is_open()){
    std::cout << "The input file is not opening" << std::endl;
    throw 1;
  }

  int eventnr;
  double nzenith,nazimuth,nenergy,czenith,cazimuth,cenergy,xpos,ypos,zpos,oneweight;

  std::vector<Cascade> Cascades;
  std::vector<double> rejected;
  std::vector<std::vector<double>> rejected_particles;

  while (cs_file  >> eventnr >> nzenith >> nazimuth >> nenergy >> czenith >> cazimuth >>
  cenergy >> xpos >> ypos >> zpos >> oneweight) {
  // Read line from file

    if (cenergy < 1E6) {              // Energy check
       std::cout << "Skipping cascade " << eventnr <<
       " since it has energy below the 1 PeV threshold" << std::endl;
       rejected = {eventnr, nzenith, nazimuth, nenergy, czenith, cazimuth,
       cenergy, xpos, ypos, zpos, oneweight};
    }

/*  This is omitted for testing purposes
     else if((zpos - 1389) > 0){   // Position check
      // 1389 = 2778/2, is the middle of the ice shelf [m]
      cout << "The cascade "<< eventnr <<
        " has a positive z-position (out-of-ice)" << endl;
      rejected = {eventnr, nzenith, nazimuth, energy, czenith, cazimuth,
      cenergy, xpos, ypos, zpos, oneweight};
    }
*/

    if (rejected.empty()){
        Cascade cs(eventnr, cenergy, xpos,ypos,zpos,
                   czenith, cazimuth, nenergy, nzenith, nazimuth, oneweight);
        Cascades.push_back(cs);
    } else {
      rejected_particles.push_back(rejected);
      rejected.clear();
    }

  }

  if (!rejected_particles.empty()){
    std::string path_out = cs_filepath.substr(0,cs_filepath.find_last_of("."))
      + "_rejected.out";
    write_2D_array(rejected_particles, path_out, 1);
  }
  return Cascades;
}


// -----------------------------------------------------------------------------

/* Helper functions, used only by other functions*/
namespace {

	/* Shower age */
	double ShowerAge(double X, double E){            // [g/cm^2, GeV]
		return (3*X/X_0) / ((X/X_0)+2*log(E/E_c));  // Unitless
	}

	/* Ne, number of ionized electrons per unit length */
	// N(X,fEnergy)					Number of particles in the cascade
	// 2E6           		[eV/ g/cm^2]  Energy loss per ionizing particle (@ 1 GeV).
	// 20            		[eV] Electron ionization energy (per electron).
	// 1E5 = 2E6/20 		[# ionized electrons/ # HE ionizing particle/ g/cm^2]
	// N*1E5  					[# electrons/ g/cm^2]
	// N*1E5*rho_ice	 	[#e / cm]
	double Ne(double X, double E){                // [g/cm^2, GeV]
		assert(X >= 0 && "Ne's X is ill defined");
		double N, Ne = 0;
	  if (X > 0){
	    N = 0.31*exp((X/X_0)*(1-1.5*log(ShowerAge(X,E))))/sqrt(log(E/E_c));
			Ne = N * 1E5 * rho_ice;
	  }
	  return Ne;       // [#e / cm]
	}

	/* Lateral particle distribution for radius r and shower age s. */
	// Last division is for normalization, missing in paper.
	double wiv1(double r, double s){      // [cm, Unitless]
	  double wiv1;
    wiv1 = exp(lgamma(4.5-s)-lgamma(s)-lgamma(4.5-2*s))
							*pow(r/r_moliere,s-1) *pow(r/r_moliere+1,s-4.5)/r_moliere;
	  return wiv1;     // [1/cm] Differential.
	}

	/* Integral of lateral particle distribution between r and r + dr. */
	double intwiv(double r, double delta_r, double s){    // [cm, cm, unitless]
		double intwiv = 0, imx = 50.0, step = delta_r/imx;
		assert(r >= 0 && "Intwiv's r < 0");
	  for (int i = 0; i < imx; i++){ intwiv += wiv1(r + i*step, s);}
	  return intwiv*step;         // [Unitless]
	}

}
