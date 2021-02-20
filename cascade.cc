#include "cascade.hh"

namespace {
	double X, r, k_mid; 		// temp variables
}

Cascade::Cascade(int event, double cenergy, double xpos, double ypos, double zpos,
								 double czenith, double cazimuth,  double nenergy,
								 double nzenith, double nazimuth, double oneweight):
				 evtnr(event), E_p(cenergy), pos{xpos, ypos, zpos}, sph_ang{czenith,cazimuth},
				 neutrino{nenergy, nzenith, nazimuth,oneweight} {

	dir[0] = sin(sph_ang[0])*cos(sph_ang[1]);
	dir[1] = sin(sph_ang[0])*sin(sph_ang[1]);
	dir[2] = cos(sph_ang[0]);

	X_tot = 4*(log(E_p/E_c)/log(2) )*X_int;
	X_bin = X_tot/nbin;

	// CORRECT FOR E 1E9, switch to E MIN 10^6 GeV AND VERIFY FOR ALL ENERGIES.
	r_tot = 5*r_moliere;//*log(E_p/1E9);
	r_bin = r_tot/nbin;

	L_tot = X_tot/rho_ice/100; 			//L in meters m.
	L_bin = L_tot/nbin;
	// Sanity check! Your sections are not unphysical due to lifetime constraint.
	assert(L_bin <= c_vac*tau && "Cascade segments are too large!");
}

/* Particle (electron) density for penetration length X and radius r. */
double Cascade::dens(double X, double r){      // [g/cm^2, cm, GeV]
	double dens, delta_r = 0.05;
	// s = sh_age(X,E_p);*step
	if (r < 0 ) {r = -r;}
	if (X < 0){X = 0;}
	dens = Ne(X,E_p) * intwiv(r, delta_r) / (pi*(pow(r + delta_r,2) - pow(r,2)));
	assert(dens >= 0 && "Negative density value");
	return dens;        // [#e-/ cm^3]
}

double Cascade::fplasma(double& dens){return (8980 * sqrt(memp) * sqrt(dens));} // [Hz]


/* Make 2D-array of density profile */
void Cascade::set_density(){
	//double X, r;
  std::vector <double> row;

  for (int i = 0; i < nbin; i++){
    X=i*X_bin;
    for (int j = 0; j < nbin; j++){
      r = j*r_bin;
      row.push_back(dens(X,r));
    }
    density.push_back(row);
    row.clear();
  }
}




/* Get the radial values of the od/ud line (r_crit).
This is where the plasma frequency equals the detection frequency. */
void Cascade::set_rcrit(const double & freq_obs){
    // double X, r;

		get_density();
		double dens_crit = memp*pow(freq_obs/8980,2);
		/* [(#e-) cm^-3]  Critical e- density for overdense scattering condition
		 wp > w > 8980*sqrt(ne), w is frequency [Hz]!! */

    // Loop over cascade depth
    for (int i = 0; i < nbin; i++){
      X = i*X_bin;
      r = 0;
      // Loop over cascade radius
      for (int j = 0; j < nbin; j++){
        // Update value of plasma radius per depth step
        if (density[i][j] > dens_crit) { r = j*r_bin;}   // [cm]
      }
      r_crit.push_back(r);
    }

    // The shower waist is located by definition at the maximum 	plasma radius.
    r_waist = *max_element(r_crit.begin(), r_crit.end());
}

// ----------------------------------------------------------------------------
// Accesors
int  Cascade::event(){return evtnr;}
double  Cascade::energy(){return E_p;}
// double* Cascade::sph_angles(){return sph_ang;}
// double* Cascade::position(){return pos;}
// double* Cascade::direction(){return dir;}
std::vector<double> Cascade::sph_angles(){return sph_ang;};
std::vector<double> Cascade::position(){return pos;}
std::vector<double> Cascade::direction(){return dir;}
double* Cascade::parent(){return neutrino;}
double  Cascade::Xtot(){return X_tot;}
double  Cascade::get_X_bin(){return X_bin;}
double 	Cascade::get_r_tot(){return r_tot;}
double  Cascade::get_r_bin(){return r_bin;}
double  Cascade::get_L_tot(){return L_tot;}
double  Cascade::get_L_bin(){return L_bin;}
double  Cascade::get_rwaist(){return r_waist;}

std::vector<std::vector<double>> Cascade::get_density(){
	// if(density.empty()){ set_density(); }
	return density;
}

std::vector<double> Cascade::get_rcrit(){
	if(r_crit.empty()){ set_rcrit(); }
	return r_crit;
}
//
// std::vector<std::vector<double>> Cascade::get_reflectivty_2D(){
// 	if(reflectivity2D.empty()) {set_reflectivity_2D();}
// 	return reflectivity2D;
// }
//
// std::vector<std::vector<double>> Cascade::get_reflectance_2D(){
// 	if(reflectance2D.empty()) {set_reflectivity_2D();}
// 	return reflectance2D;
// }

// -----------------------------------------------------------------------------
std::vector<Cascade> load_cascade_file(const std::string& cs_filepath){
  std::ifstream cs_file(cs_filepath);

  if (!cs_file.is_open()){
    std::cout << "The input file is not opening" << std::endl;
    //outfile << "The input file is not opening" << endl;
    //exit (EXIT_FAILURE);
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

namespace {

	/* Helper functions, used only by other functions*/

	/* Shower age */
	double sh_age(double X, double E_p){            // [g/cm^2, GeV]
		return (3*X/X_0) / ((X/X_0)+2*log(E_p/E_c));  // Unitless
	}

	/* Ne, number of ionized electrons per unit length */
	// N(X,E_p)					Number of particles in the cascade
	// 2E6           		[eV/ g/cm^2]  Energy loss per ionizing particle (@ 1 GeV).
	// 20            		[eV] Electron ionization energy (per electron).
	// 1E5 = 2E6/20 		[# ionized electrons/ # HE ionizing particle/ g/cm^2]
	// N*1E5  					[# electrons/ g/cm^2]
	// N*1E5*rho_ice	 	[#e / cm]
	double Ne(double X, double E_p){                // [g/cm^2, GeV]
		assert(X >= 0 && "Ne's X is ill defined");
		double N, Ne = 0;
	  if (X > 0){
	    N = 0.31*exp((X/X_0)*(1-1.5*log(sh_age(X,E_p))))/sqrt(log(E_p/E_c));
			Ne = N * 1E5 * rho_ice;
	  }
	  return Ne;       // [#e / cm]
	}

	/* Lateral particle distribution for radius r and shower age s. */
	// Last division is for Normalization, missing in paper.
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
