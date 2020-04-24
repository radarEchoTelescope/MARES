#include "cascade.hh"

namespace {
	double X, r, k_mid; 		// temp variables
}

Cascade::Cascade(int evt, double position, double cenergy, double cazimuth, double czenith,
				 double nenergy, double nzenith, double nazimuth, double oneweight):
				 eventnr(evt), energy(cenergy), pos{position}, sph_ang{czenith,cazimuth},
				 neutrino{nenergy, nzenith, nazimuth,oneweight} {

	dir[0] = sin(sph_ang[0])*cos(sph_ang[1]);
	dir[1] = sin(sph_ang[0])*sin(sph_ang[1]);
	dir[2] = cos(sph_ang[0]);

	X_tot = 4*(log(energy/E_c)/log(2) )*X_int;
	X_bin = X_tot/nbin;

	// CORRECT FOR E 1E9, switch to E MIN 10^6 GeV AND VERIFY FOR ALL ENERGIES.
	r_tot = 5*r_moliere;//*log(energy/1E9);
	r_bin = r_tot/nbin;

	L_tot = X_tot/rho_ice/100; 			//L in meters m.
	L_bin = L_tot/nbin;
	// Sanity check! Your sections are not unphysical due to lifetime constraint.
	assert(L_bin <= c_vac*tau && "Overdense sections are too large!");
}

/* Make 2D-array of density profile */
void Cascade::set_density(){
	//double X, r;
  std::vector <double> row;

  for (int i = 0; i < nbin; i++){
    X=i*X_bin;
    for (int j = 0; j < nbin; j++){
      r = j*r_bin;
      row.push_back(dens(X,r,energy));
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
		double dens_crit = mme*pow(freq_obs/8980,2);
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

void Cascade::set_reflectivity_2D(){
	// double X, X_bin, r, r_bin, k_mid;
	double skin, dr;                      // [cm]
	double wplasma;                       // [Hz]

	double reflectance, reflectivity;      // Unitless
	std::vector<double> reflectance_row, reflectivity_row;

	// Loop over depth
	for (int i = 0; i < nbin; i++){
		X=i*X_bin;
		reflectance = 0, reflectivity = 0;
		// Loop over layers
		for (int k = 0; k < nbin; k++){
			// We actually need to move from the outside - in!
			k_mid = nbin - (k + 0.5);
			// Fix for the last layer variable size (since we shift them).
			k == nbin ? dr = 0.5*r_bin : dr = r_bin;

			wplasma=8980*sqrt(dens(X,k_mid*r_bin,energy))*sqrt(1/mme);
			skin=cmed_cm/(2*wplasma);

			reflectance = (1-reflectivity)*(1-exp(-1*dr/skin));
			reflectivity += reflectance;
			assert(reflectivity < 1 && "Reflectivity larger than 1!");

			reflectance_row.push_back(reflectance);
			reflectivity_row.push_back(reflectivity);
		}
		reflectance2D.push_back(reflectance_row);
		reflectivity2D.push_back(reflectivity_row);
		reflectance_row.clear();
		reflectivity_row.clear();
	}
}

/* Compute the overdense area of the cascade in slices.
The reflectance corrections are only applied to the od sections. */
void Cascade::set_od_cs_1D(){
	// double X, r_bin, L_bin, k_mid;
	double od_cs;

	// The arrays of k and ref are no longer needed.
	// std::vector<double> row_k;
	// std::vector<double> row_ref;
	// std::vector<std::vector<double>> array_k;
	// std::vector<std::vector<double>> array_ref;

	get_rcrit();
	get_reflectance_2D();

	// Loop over regions.
	for (int i = 0; i < nbin; i++){

		// Number of layers per segment.
		int k_max = (int) r_crit[i]/r_bin;
		od_cs = 0;
		for (int k = 0; k < nbin; k++){
			// Only layers from od_region are looked, outside-in.
			k >= (nbin - k_max) ? k_mid = nbin - (k + 0.5) : k_mid = 0;
			// To evaluate for the full space as od
			// k_mid = nbin - (k + 0.5);

			od_cs += 2*k_mid*r_bin*reflectance2D[i][k]*L_bin;      // [cm^2]
			// Factor of 2 required to account for plasma tube symmetry!

			// row_k.push_back(k_mid);
			// row_ref.push_back(reflectance[i][k]);
		}
		od_cs_1D.push_back(od_cs);
		// array_k.push_back(row_k);
		// array_ref.push_back(row_ref);
		// row_k.clear();
		// row_ref.clear();
	}
}

/* Compute the total overdense surface of the cascade */
void Cascade::set_od_cs_0D(double tx_dot_cs, double cs_dot_rx, double lambda){
	//double X,r;
	get_rcrit();

	/*First, get some critical shower values*/
	double dens_max = 0, X_max;
	int i_max;
	// Loop over cascade depth
	for (int i = 0; i < nbin; i++){
		X = i*X_bin;
		r = 0;
		// Loop over cascade radius
		for (int j = 0; j < nbin; j++){

			// Update the maximum values
			if (density[i][j] > dens_max){
				dens_max = density[i][j];
				X_max = X;                // Penetration depth at max density.
				i_max = i;
			}
		}
	}

	// If shower age = 1, r_waist is at the maximum density.
	assert(r_waist == r_crit[i_max] && "Ill-defined waist radius.");

	/* Second, find the range of layers of "constant" density.
		These layers need to keep the radial spacing r_bin.
	*/
	// double k_mid;              // [unitless]
	double len;										// [cm]
	double dens_layer;            // [g/cm^3]
	double X_low, X_high;					// [g/cm^2]
	std::vector<double> length;
	int k_max = (int) r_waist/r_bin;

	// Loop over od layers
	for(int k = 0; k < nbin; k++){
		k_mid = k + 0.5;
		// Density of midpoint of k layer at waist
		dens_layer=dens(X_max,k_mid*r_bin, energy);
		X_low = X_tot;
		X_high = 0;
		// Loop over cascade depth
		for (int i = 0; i < nbin; i++){
			X=i*X_bin;
			// Loop over cascade radius until waist, where the od values are defined.
			for (int j = 0; j < k_max; j++){
				// Update size of the od layer
				if(density[i][j] > dens_layer){
					if(X<X_low) {X_low=X;}
					if(X>X_high) {X_high=X;}
				}
			}
		}
		len = (X_high - X_low)/rho_ice;
		assert(len > 0 && "Layer length below zero!");
		// If length is bigger than lifetime constraint, update else 1 (Do nothing).
		len > cvac_cm*tau ? len = cvac_cm*tau: 1;
		length.push_back(len);
	}

	/* Thirdly, we find the total surface by integrating the reflectance:
	We try to  save some time if the reflectance2D values are already computed. */
	double od_tot;
	// if (!reflectance2D.empty()){
	// 	for (int k = 0; k < nbin; k++){
	// 		od_tot += 2*k_mid*r_bin*length[k]*reflectance2D[i_max][k];      // [cm^2]
	// 		/* Factor of 2 needed for shower symmetry. */
	// 	}
	// } else {
	/* We compute the reflectance layers along the waist (max particle density) */
		double skin, dr;                      // [cm]
		double wplasma;                       // [Hz]
		double reflectivity, reflectance; // Unitless

		for (int k = 0; k < nbin; k++){
			// k goes from outside to inside.
			k_mid = nbin - k + 0.5;
			// k is placed between radius layers, so the last layer is half size.
			k == nbin ? dr = 0.5*r_bin : dr = r_bin;              // [cm]

			//Directly take into account for skin effects in radar cross-section determination
			wplasma=8980*sqrt(dens(X_max,k_mid*r_bin,energy))*sqrt(1/mme);  // [Hz]
			skin=cmed_cm/(2*wplasma);                             // [cm] skin-depth

			reflectance = (1-reflectivity)*(1-exp(-1*dr/skin));
			reflectivity += reflectance;
			assert(reflectivity < 1 && "Reflectivity larger than 1!");

			od_tot += 2*k_mid*r_bin*length[k]*reflectance;      // [cm^2]
		}
	// }

	/* Finally, we factor the geometrical corrections for the od surface.
	f_geometry: norm_in * G3  * Polarization
		- norm_in = Correction factor for angle of normal incidence
			(line of sight/apparent size).
		- G3 = shower boost of radio signal "directional gain"
			Right now, it is obtained by f_diff, the single slit diffraction pattern.
		-	Polarization = Missing here!
	*/
	double f_geometry, norm_in, f_diff = 0;

		// Compute norm_in
	norm_in  = 1 - abs(tx_dot_cs);

		// Compute f_diff
	double l_max, gam, av_int = 0, ang, steps = 200000.0;

	// Angle: Distance to relative to full internal reflection
	// "gam" = "alpha" - "beta" in the paper
	gam = acos(tx_dot_cs) - acos(cs_dot_rx);
// BIG RED FLAG WITH ACOS FUNCTION, IT HAS FAILED IN THE PAST.
// cos(a) = cos(-a) -> acos(cos(-a)) = a

	l_max = *max_element(length.begin(), length.end());
	// Average intensity over the -pi/2 to pi/2 region
	for (int i = 0; i < steps ; i++){
		ang = (-1/2.0 + (i)/steps)*pi; // Some integration value
		av_int += pow(sinc((pi*l_max/lambda)*sin(ang)),2);     // Add contibutions
	}
	av_int /= steps;
	// Intentsity at gamma angle.
	f_diff =  pow(sinc((pi*l_max/lambda)*sin(gam)),2)/av_int;

	f_geometry=abs(norm_in*f_diff);

	od_cs_0D = od_tot* f_geometry;

}
/* Compute under-dense cross-section */
void Cascade::set_ud_cs_0D(){
	/*
	r_crit[i] is the overdense limit for radius i
	r_tot is the radious of the plasma rod that we consider
	Ne(X) is the number of ionized electrons per charge per unit depth.
	intwiv(rlow,rup,X) is the fraction of under-dense charges for depth X
	*/
	double N_ch=0;
	get_rcrit();
	for (int i = 0; i < nbin; i++){
		X=i*X_bin;
		N_ch += Ne(X,energy)*intwiv(r_crit[i],r_tot-r_crit[i],1.01);
	}
	// Number of ud-scattering electrons for the full length
	ud_cs_0D = N_ch*X_bin*pow((1/mme),2)*thompson;               // [cm^2]
	// Where is alpha? Why == 1?? Should not be set to 2 by definition?
}

// -------------------------------------------
// abs((pi*pow(r,2))-(pi*pow(r + delta_r,2)));
// Acess values
double  Cascade::get_event(){return eventnr;}
double  Cascade::get_energy(){return energy;}
double* Cascade::get_sph_ang(){return sph_ang;}
double* Cascade::get_position(){return pos;}
double* Cascade::get_direction(){return dir;}
double* Cascade::get_parent_neutrino(){return neutrino;}
double  Cascade::get_max_depth(){return X_tot;}
double  Cascade::get_X_bin(){return X_bin;}
double 	Cascade::get_r_tot(){return r_tot;}
double  Cascade::get_r_bin(){return r_bin;} 
double  Cascade::get_L_bin(){return L_bin;}

std::vector<std::vector<double>> Cascade::get_density(){
	if(density.empty()){ set_density(); }
	return density;
}

std::vector<double> Cascade::get_rcrit(){
	if(r_crit.empty()){ set_rcrit(); }
	return r_crit;
}

std::vector<std::vector<double>> Cascade::get_reflectivty_2D(){
	if(reflectivity2D.empty()) {set_reflectivity_2D();}
	return reflectivity2D;
}

std::vector<std::vector<double>> Cascade::get_reflectance_2D(){
	if(reflectance2D.empty()) {set_reflectivity_2D();}
	return reflectance2D;
}

std::vector<double> Cascade::get_od_cs_1D(){
	if(od_cs_1D.empty()) {set_od_cs_1D();}
	return od_cs_1D;
}

double Cascade::get_od_cs_0D(double tx_dot_cs, double cs_dot_rx, double lambda){
	if(od_cs_0D == 0) {set_od_cs_0D(tx_dot_cs, cs_dot_rx, lambda);}
	return od_cs_0D;
}

double Cascade::get_ud_cs_0D(){
	if(ud_cs_0D == 0) {set_ud_cs_0D();}
	return ud_cs_0D;
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
	double wiv1(double r, double s){      // [cm, Unitless]
	  double wiv1;
	  if(s == 1.01) { // default value for s
			wiv1 = 2.5217 * pow(r/r_moliere,s-1) *pow(r/r_moliere+1,s-4.5) /r_moliere;
		} else {
	    wiv1 = exp(lgamma(4.5-s)-lgamma(s)-lgamma(4.5-2*s))
								*pow(r/r_moliere,s-1) *pow(r/r_moliere+1,s-4.5)/r_moliere;
	  }
	  return wiv1;     // [1/cm] Differential.
	}
	// Last division is for Normalization, missing in paper ?

	/* Integral of lateral particle distribution between r and r + dr. */
	double intwiv(double r, double delta_r, double s){    // [cm, cm, unitless]
		assert(r >= 0 && "Intwiv's r < 0");
	  double intwiv = 0, imx = 50.0, step = delta_r/imx;
	  for (int i = 0; i < imx; i++){ intwiv += wiv1(r + i*step, s);}
	  return intwiv*step;         // [Unitless]
	}

	/* Particle (electron) density for penetration length X and radius r. */
	double dens(double X, double r, double E_p){      // [g/cm^2, cm, GeV]
		double dens, delta_r = 0.05;
		// s = sh_age(X,E_p);*step
		dens = Ne(X,E_p) * intwiv(r, delta_r) / (pi*(pow(r + delta_r,2) - pow(r,2)));
		assert(dens >= 0 && "Negative density value");
		return dens;        // [#e-/ cm^3]
	}

	/* Sinc function*/
	double sinc(double x){
	  double tmp = 0;
	  abs(x) > 10E-8 ? tmp = sin(x)/x : tmp = 1;
	  return tmp;
	}

}
