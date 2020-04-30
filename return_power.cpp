/* Radar computation - return electric field */
#include "macro_scatter.hh"

//---- Free (user) parameters---------------------------------------------------
// const std::string identifier = "return_power_0D";


double sinc(double x);
double get_od_cs_0D(Antenna& tx, Antenna& rx, Cascade& cs);
double get_ud_cs_0D(Cascade & cs);

int main(int argc, char** argv){

  // I/O Setup
  string cs_filepath = argv[1];
  // Read data from file, get separate unphysical events
  vector<Cascade> Cascades = load_cascade_file(cs_filepath);

  // Make your detector
  string det_type = "bistatic";
  Detector bistatic(det_type);

  // Compute the return power
  double Rt, Rr;
  double eta, Latt, rad_cs, P_r, P_bkg;
  int triggered;

  // For every cascade
  for (auto& cs: Cascades){
    // For every transmitter
    for (auto& tx : bistatic.get_transmitters()){
      // For every receiver
      for (auto& rx : bistatic.get_receivers()){

        // cout <<
        //       rx.power() << " " <<
        //       rx.position()[0] << " " <<
        //       rx.position()[1] << " " <<
        //       rx.position()[2] << " " <<
        //       rx.distance() << " " <<
        // endl;
        //
        // cout <<
        //       rx.power << " " <<
        //       rx.pos[0] << " " <<
        //       rx.pos[1] << " " <<
        //       rx.pos[2] << " " <<
        //       rx.dist << " " <<
        // endl;

        Scatter event(tx,rx,cs);

        Antenna new_tx =  event.transmitter();
        Antenna new_rx =  event.receiver();

        // set_direction(cs.get_position(), cs.get_direction());
        // rx.set_direction(cs.get_position(), cs.get_direction());
        Rt = new_tx.distance();
        Rr = new_rx.distance();

        cout  << Rt << '\t'
              << Rr << '\t'
              << tx.distance() << '\t'
              << rx.distance() << '\t'
              << endl;

        double od_cs = get_od_cs_0D(new_tx, new_rx, cs);
        double ud_cs = get_ud_cs_0D(cs);

        vector<vector<double>> dens = cs.get_density();
      	rad_cs = (od_cs + ud_cs)/1E4;   	     // [m^2] !! Unit change happens here!

	// Inefficiency --------------------------------------------------------------
      	eta = 1;                     // For testing purposes
      	// eta=4.4E-3;//*1.d-2          // @450 MHz; 20ns
      	//eta=0.18d0*1.d-2              // @50 MHz; 1000ns
      	//eta=2.64d-5                   // @10 MHz superdarn detected @ 200-500 MHz

      	Latt = 1000;                 // [m] attenuation length
        //Latt=att(freq_obs);        // It should be parametrized depeding on location

      	P_r = eta*tx.power()*tx.gain()*rx.gain()*      // Constants
      	exp(-2.0*(Rt+Rr)/Latt)/pow(4*pi*Rt*Rr,2)*    // Distance-dependant term
      	rad_cs;

        cout << tx.power() << '\t'
              << tx.gain() << '\t'
              << rx.gain() << '\t'
              << endl;

        // Check if the event is over thermal background
      	double k_b=1.3806503E-23;       // [SI] Boltzmann's constant
      	double T_sys=325;               // [K] system temperature from ARA paper
      	double deltaf=1E5;              // [Hz]

      	P_bkg=k_b*T_sys*deltaf/1E3;         //[kW]
      	(P_r > P_bkg) ? triggered = 1 : triggered = 0;

      } // Closes antenna loop
    } // Closes transmitter loop

    cout << rad_cs << '\t' << P_r << '\t' << triggered << endl;

/*
    cout <<
      "Event:" << scientific << '\t' <<
      "R_t = " << tx.R + rx.R << '\t' <<
      //"ud_cs = " << ud_cs << '\t' <<
      "od_tot = " << od_tot << '\t' <<
      "rad_cs = " << rad_cs << '\t' <<
      "P_r = " << P_r << defaultfloat <<
    endl;
*/

  } // Closes cascade loop
  cout << "END" << endl;
}
// End of main

/* Compute the total overdense surface of the cascade */
double get_od_cs_0D(Antenna& tx, Antenna& rx, Cascade& cs){
  double X, r, k_mid;       // Variable
  double od_cs_0D;
  double X_tot  = cs.Xtot();
  double X_bin  = cs.get_X_bin();
  double r_bin  = cs.get_r_bin();
  double lambda = tx.lambda();
  double tx_dot_cs = tx.projection();
  double cs_dot_rx = rx.projection();
	std::vector<double> r_crit = cs.get_rcrit();
  std::vector<std::vector<double>> density = cs.get_density();

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

  double r_waist = r_crit[i_max];
	// // If shower age = 1, r_waist is at the maximum density.
	// assert(r_waist == r_crit[i_max] && "Ill-defined waist radius.");

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
		dens_layer= cs.dens(X_max,k_mid*r_bin);
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
	We compute the reflectance layers along the waist (max particle density) */
	double od_tot;
	double skin, dr;                      // [cm]
	double wplasma;                       // [Hz]
	double reflectivity = 0, reflectance = 0; // Unitless

	for (int k = 0; k < nbin; k++){
		// k goes from outside to inside.
		k_mid = nbin - k + 0.5;
		// k is placed between radius layers, so the last layer is half size.
		k == nbin ? dr = 0.5*r_bin : dr = r_bin;              // [cm]

		//Directly take into account for skin effects in radar cross-section determination
		wplasma=8980*sqrt(cs.dens(X_max,k_mid*r_bin))*sqrt(1/mme);  // [Hz]
		skin=cmed_cm/(2*wplasma);                             // [cm] skin-depth

		reflectance = (1-reflectivity)*(1-exp(-1*dr/skin));
		reflectivity += reflectance;
		assert(reflectivity < 1 && "Reflectivity larger than 1!");

		od_tot += 2*k_mid*r_bin*length[k]*reflectance;      // [cm^2]
	}

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

  return od_cs_0D;
}


/* Compute under-dense cross-section */
double get_ud_cs_0D(Cascade & cs){
	double ud_cs, N_ch=0, r_bin = cs.get_r_bin();
	std::vector<double> r_crit = cs.get_rcrit();
  std::vector<std::vector<double>> density = cs.get_density();

	for (int i = 0; i < nbin; i++){
    int r_min = (int) r_crit[i]/r_bin ;
    for (int j = r_min; j < nbin; j++){
      N_ch += density[i][j]* pi*(pow((j+1)*r_bin,2) - pow(j*r_bin,2));
    }
	}
	// Number of ud-scattering electrons for the full length
	ud_cs = N_ch*cs.get_X_bin()*pow((1/mme),2)*thompson;               // [cm^2]
	// Where is alpha? Why == 1?? Should not be set to 2 by definition?
  return ud_cs;
}

/* Sinc function*/
double sinc(double x){
  double tmp = 0;
  abs(x) > 10E-8 ? tmp = sin(x)/x : tmp = 1;
  return tmp;
}
//------------------------------------------------------------------------------
/* NO LONGER IN USE */

/*
// Parametrized attenuation length for the Ross Ice Shelf, South Pole.
double att(double freq_obs){
double a1=469;                  // [m] Attanuation length parameter
double a2=-0.205;               // Attanuation length parameter
double a3=4.87E-5;              // Attanuation length parameter
return a1+a2*freq_obs/1E6+a3*pow(freq_obs/1E6,2);
}
*/
