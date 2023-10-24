#include "cascade.hh"

Cascade::Cascade(int event, double cenergy, double xpos, double ypos, double zpos,
								 double czenith, double cazimuth, int np):
				 Cascade(event, cenergy, xpos, ypos, zpos, czenith,cazimuth, np,
								 cenergy, czenith, cazimuth, 1){}

Cascade::Cascade(int event, double cenergy, double xpos, double ypos, double zpos,
								 double czenith, double cazimuth,  int np, 
                 double nenergy, double nzenith, double nazimuth, double oneweight):
				 			   fEvent(event),
								 fEnergy(cenergy),
                 fNp(np),
								 fPosition{xpos, ypos, zpos},
								 fSphericalAngles{czenith,cazimuth},
								 fNeutrino{nenergy, nzenith, nazimuth,oneweight} {

	fDirection[0] = sin(fSphericalAngles[0])*cos(fSphericalAngles[1]);
	fDirection[1] = sin(fSphericalAngles[0])*sin(fSphericalAngles[1]);
	fDirection[2] = cos(fSphericalAngles[0]);

	/* X_tot 	= 4* number of divisions for shower_max * interaction length
				  = 4*(log(fEnergy/E_c)/log(2) )*X_int;
		 			= 4*(log(fEnergy/E_c)) * X_0
 */

	// fXtot = 4 * log(12.72 * fEnergy) * X_0;	
 	fXtot = 3 * log(12.72 * fEnergy) * X_0;			
	fLtot = fXtot/rho_ice; 											
	fRtot = 2*r_moliere;
	// fRtot = 3*r_moliere;											

	/* If you want your cascade to scale with energy in both dimensions,
	this should work for any energy above 1 PeV / 10^6 GeV/ 10^15 eV.
	*/
	// fRtot = r_moliere*(log(fEnergy/1E5));	

// We calculate the integral of the radial profile for a given total radius.
// This should used to normalize all other integrals using intwiv. 
  // fRnorm = intwiv(0, fRtot, 1.01, 1000);
  //   std::cout << fRnorm << std::endl; 


}

// Cascade's Methods

/* Ne, number of ionized electrons per unit length */
  // Usual N(e), 1 particle of energy Ep.
double Cascade::Ne(double X){
	return Ne(X, fEnergy, fNp);
}

  // N(e)/mm, applies to a cascade of Np particles with energy Ep per particle. 
double Cascade::Ne(double X, double Ep, double Np){
	return Np* NKG::N(X,Ep)*         
          // Number of ionizing primaries
          E_deposition/E_ionization * rho_ice;
          // # e/mm per ionising primary
}

  // Usual density function, 1 particle of energy Ep.
double Cascade::Density(double X, double r, double delta_r){
	 return Density(X, r, delta_r, fEnergy, 1.);
             // Weight for the fraction inside the r+dr ring.
 }

 // Density for beam containing N particles of energy E (per particle).
double Cascade::Density(double X, double r, double delta_r, double E, double Np){
	double dens;
  if(r<0.0){r = -r;}
	// s = ShowerAge(X,fEnergy);*step
  dens = Ne(X, E, Np)*
        NKG::intwiv(r,delta_r) 
        / (pi*(pow(delta_r,2) + 2*r*delta_r));
	// Equivalent to:
  //  / (pi*(pow(r + delta_r,2) - pow(r,2)));
	assert(dens >= 0 && "Negative density value"); // Sanity check
	return dens;
}

// This formula only works if the electron density is in cm^3!!
// Also, this formula returns frequency in Hz.
double Cascade::PlasmaFreq(const double &dens){ 
  return (8980 * sqrt(memp) * sqrt(dens*pow(cm,3)))*Hz ; 
}

double Cascade::Absorption(const double &dens, const double &freq_obs, const double &freq_nat){
	double w_sq, a, b, q, fnat_sq, fobs_sq, fplasma_sq;     // Some temp varables.
  double fplasma = PlasmaFreq(dens);

	/* Collisionless model
      This is meant to protect against error, not to actually run in 
      a unphysical setting. Might add error code instead.
  */
  if (f_coll == 0){
    std::cout << "You are running an unphysical collsionless model" << std::endl;
		fplasma > freq_obs ? q = 2*pi*fplasma/c_ice : q = 0;               // [1/mm]
		return q;
	}

	fnat_sq    = pow(freq_nat,2);
	fobs_sq    = pow(freq_obs,2);
	fplasma_sq = pow(fplasma,2);

	// Exact solution from dispersion relation with collisions.
  w_sq = 1.0/( pow( fnat_sq - fobs_sq,2) + fobs_sq*pow(f_coll,2) );
	// a = 1 + wp^2 * Re(W)
  a = 1 + fplasma_sq * (fnat_sq - fobs_sq)*w_sq ;
	// b = w_p^2 * Im(W)
  b = fplasma_sq * freq_obs * f_coll * w_sq;

  // double p = (freq_obs/c_vac)*np.sqrt((np.sqrt(a**2 + b**2) + a) /2 )
  q = (2*pi*freq_obs/c_vac)*sqrt((sqrt(pow(a,2) + pow(b,2)) - a) /2 ); // [1/mm]
  
  return q;
}

double Cascade::SkinDepth(const double &dens, const double &freq_obs, const double &freq_nat){
	return 1./Absorption(dens, freq_obs, freq_nat);
}

//-------- These are the matrix methods for the 2D storage variables ----------

void Cascade::Ne( const std::vector<std::vector<double>> &length_vals,
                       const std::vector<std::vector<double>> &radius_vals,
										   const double delta_r){
  fNe = std::vector<std::vector<double>> (length_vals.size(),
              std::vector<double> (length_vals[0].size(), 0));

  for (int i = 0; i < fNe.size(); i++){
    for (int j = 0; j < fNe[i].size(); j++){
      // Simple check to avoid computing values too far out from the cascade direction
      if ( abs(radius_vals[i][j]) <= fRtot ) {
        fNe[i][j] = Ne(rho_ice*length_vals[i][j], radius_vals[i][j], delta_r);
      }
    }
  }
}

void Cascade::Density( const std::vector<std::vector<double>> &length_vals,
                       const std::vector<std::vector<double>> &radius_vals,
										   const double delta_r){
  fDensity = std::vector<std::vector<double>> (length_vals.size(),
              std::vector<double> (length_vals[0].size(), 0));

  for (int i = 0; i < fDensity.size(); i++){
    for (int j = 0; j < fDensity[i].size(); j++){
      // Simple check to avoid computing values too far out from the cascade direction
      if ( abs(radius_vals[i][j]) <= fRtot ) {
        fDensity[i][j] = Density(rho_ice*length_vals[i][j], radius_vals[i][j], delta_r);
      }
    }
  }
}

void Cascade::PlasmaFreq(const std::vector<std::vector<double>> &density){
  fPlasmaFrequency = std::vector<std::vector<double>> (density.size(),
                      std::vector<double> (density[0].size()));

  for (int i = 0; i < density.size(); i++){
    for (int j = 0; j < density[i].size(); j++){
      if (density[i][j] != 0) {fPlasmaFrequency[i][j] = Cascade::PlasmaFreq(density[i][j]);}
    }
  }
}

void Cascade::Absorption(const std::vector<std::vector<double>> &density, const double & freq){
  fAbsorption = std::vector<std::vector<double>> (density.size(),
                 std::vector<double> (density[0].size()));

  for (int i = 0; i < density.size(); i++){
    for (int j = 0; j < density[i].size(); j++){
      if (density[i][j] != 0) {fAbsorption[i][j] = Cascade::Absorption(density[i][j], freq);}
    }
  }
}

void Cascade::SkinDepth(const std::vector<std::vector<double>> &density, const double & freq){
  fSkinDepth = std::vector<std::vector<double>> (density.size(),
                std::vector<double> (density[0].size()));

  for (int i = 0; i < density.size(); i++){
    for (int j = 0; j < density[i].size(); j++){
      if (density[i][j] != 0) {fSkinDepth[i][j] = Cascade::SkinDepth(density[i][j], freq);}
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
/* For consistency and usefulness, the transparency and opacity are defined
___before___ reaching a certain layer.

transparency after layer = transparency before layer * transmitivity
power scattered per later = transparency before layer * reflectivity
This was the old "reflectivity" definition.
*/

void Cascade::Transparency(const std::vector<std::vector<double>> &density,
                            const double & freq, const double & delta_r){
  double transparency, transmitivity;
  fTransparency = std::vector<std::vector<double>> (density.size(),
                  std::vector<double> (density[0].size(), 1));
  fOpacity = std::vector<std::vector<double>> (density.size(),
                  std::vector<double> (density[0].size(), 0));
  fReflectance = std::vector<std::vector<double>> (density.size(),
                  std::vector<double> (density[0].size(), 0));

	// Loop over the density matrix
  for (int i = 0; i < density.size(); i++){
    transmitivity = 1; transparency = 1;
    for (int j = 0; j < density[0].size(); j++){
      if(density[i][j] == 0){continue;}

      fTransparency[i][j] = transparency;
      fOpacity[i][j] = 1 - transparency;

      transmitivity = exp(-1.0 * delta_r * Absorption( density[i][j], freq) );
      transparency *= transmitivity;
      assert(transparency > 0 && "Opacity larger than 1!");   // sanity check
      fReflectance[i][j] = 1 - transmitivity;
		}
	}
}


// ----------------------------------------------------------------------------
// Accesors
int  Cascade::Evt()	const {return fEvent;}
double  Cascade::Energy() const {return fEnergy;}
double Cascade::Ltot() const {return fLtot;}
double Cascade::Rtot() const {return fRtot;}
double Cascade::Xtot() const {return fXtot;}
std::vector<double> Cascade::Pos() const {return fPosition;}
std::vector<double> Cascade::Dir() const {return fDirection;}
std::vector<double> Cascade::Sph() const {return fSphericalAngles;};
double* Cascade::Parent() {return fNeutrino;}

std::vector<std::vector<double>> Cascade::Ne()  { return fNe; }
std::vector<std::vector<double>> Cascade::Density()  { return fDensity; }

// If the values are not necessary to be stored for TCS,
// then the values are only stored upon demand
std::vector<std::vector<double>> Cascade::PlasmaFreq(){
  if (fPlasmaFrequency.empty()) {PlasmaFreq(fDensity);}
  return fPlasmaFrequency;
}

std::vector<std::vector<double>> Cascade::Absorption(){
  if (fAbsorption.empty()) {
    std::cerr << "Absorption has not been set!" << std::endl;
    exit(EXIT_FAILURE);}
  return fAbsorption;
}

std::vector<std::vector<double>> Cascade::SkinDepth(){
  if (fSkinDepth.empty()) {
    std::cerr << "Absorption/Skin depth has not been set!" << std::endl;
    exit(EXIT_FAILURE);
    }
  return fSkinDepth;
}

std::vector<std::vector<double>> Cascade::Reflectance(){ 
  if (fReflectance.empty()){
    std::cerr << "Reflectance (Transparency) has not been set!" << std::endl;
    exit(EXIT_FAILURE);
    }
  return fReflectance; 
}
std::vector<std::vector<double>> Cascade::Opacity(){
  if (fOpacity.empty()){
    std::cerr << "Opacity (Transparency) has not been set!" << std::endl;
    exit(EXIT_FAILURE);
    } 
  return fOpacity; 
}
std::vector<std::vector<double>> Cascade::Transparency(){
  if (fTransparency.empty()){
    std::cerr << "Transparency has not been set!" << std::endl;
    exit(EXIT_FAILURE);
  }
   return fTransparency; 
}
// ----------------------------------------------------------------------------

void Cascade::save_output_files(const std::string& output_path, const std::array<bool, 8>& flags){

  if(flags[0]){ write_2D_array(Ne(),          output_path + "_e_number.txt");}
  if(flags[1]){ write_2D_array(Density(),     output_path + "_density_tx.txt");}
  if(flags[2]){ write_2D_array(PlasmaFreq(),  output_path + "_plasma_freq.txt");}
  if(flags[3]){ write_2D_array(Absorption(),  output_path + "_absorption.txt");}
  if(flags[4]){ write_2D_array(SkinDepth(),   output_path + "_skin_depth.txt");}
  if(flags[5]){ write_2D_array(Reflectance(), output_path + "_reflectance.txt");}
  if(flags[6]){ write_2D_array(Opacity(),     output_path + "_opacity.txt");}
  if(flags[7]){ write_2D_array(Transparency(),output_path + "_transparency.txt");}
}

// -----------------------------------------------------------------------------
std::vector<Cascade> load_cascade_file(const std::string& cs_filepath){
  std::ifstream cs_file(cs_filepath);

  if (!cs_file.is_open()){ 
    std::cout << "The input file is not opening" << std::endl;
    throw 1;
  }

  std::string path_out = cs_filepath.substr(0,cs_filepath.find_last_of("."))
    + "_rejected.out";

  int eventnr, np;
  bool rejected = false; 
  double nzenith,nazimuth,nenergy,czenith,cazimuth,cenergy,xpos,ypos,zpos,oneweight;
  std::vector<Cascade> Cascades;

  while (cs_file  >> eventnr >> nzenith >> nazimuth >> nenergy >> np >> czenith >> cazimuth >>
  cenergy >> xpos >> ypos >> zpos >> oneweight) {
  // Read line from file

    if (cenergy < 1E6) {              // Energy check
       std::cout << "Skipping cascade " << eventnr <<
       " since it has energy below the 1 PeV threshold" << std::endl;
       rejected = true;
    }

/*  This is omitted for testing purposes
     else if((zpos - 1389) > 0){   // Position check
      // 1389 = 2778/2, is the middle of the ice shelf [m]
      cout << "The cascade "<< eventnr <<
        " has a positive z-position (out-of-ice)" << endl;
      rejected = true;
    }
*/

    if (!rejected){
        Cascade cs(eventnr, cenergy, xpos,ypos,zpos, czenith, cazimuth,
          np, nenergy, nzenith, nazimuth, oneweight);
        Cascades.push_back(cs);
    } else {
      std::ofstream log(path_out, std::ios_base::app | std::ios_base::out);
      log << eventnr << nzenith  << nazimuth  << nenergy << np << czenith << cazimuth
       << cenergy << xpos << ypos << zpos << oneweight << std::endl;
    }

  }
  return Cascades;
}


// -----------------------------------------------------------------------------

// TODO Add SImon's param in different namespace
// TODO Take intwiv out, add the hypergeomtric function first. 

/* Helper functions, used only by other functions*/
namespace NKG{

	/* Shower age */
	double ShowerAge(double X, double E){            // [g/mm^2, MeV]
		return (3*X/X_0) / ((X/X_0)+2*log(E/E_c));  // Unitless
	}

	/* Number of charged particles in the cascade */
	double N(double X, double E){                // [g/cm^2, GeV]
		double N = 0;
		if (X > 0) {
	    N = 0.31*exp((X/X_0)*(1-1.5*log(ShowerAge(X,E))))/sqrt(log(E/E_c));
		}
	  return N;
	}

	/* Lateral particle distribution for radius r and shower age s. */
	// Last r_moliere is for normalization, missing in paper.
	double wiv1(double r, double s){      // [mm, Unitless]
	  double wiv1;
    wiv1 = exp(lgamma(4.5-s)-lgamma(s)-lgamma(4.5-2*s))
							*pow(r/r_moliere,s-1) *pow(r/r_moliere+1,s-4.5) /r_moliere;
	  return wiv1;     // [1/mm] Differential.
	}

	/* Integral of lateral particle distribution between r and r + dr.	*/
	double intwiv(double r, double delta_r, double s, double imx){
		double intwiv = 0, step = delta_r/imx;
		assert(r >= 0 && "Intwiv's r < 0");
	  for (int i = 0; i <= imx; i++){
			if(i == 0 || i == imx){
				intwiv += wiv1(r + i*step, s) / 2.0 ;
			} else {
				intwiv += wiv1(r + i*step, s);
			}
		}
	  return intwiv*step;         // [Unitless]

	}

}
