/* Radar computation - return electric field */
#include "macro_scatter.hh"

// ----- Computational constants -----------------------------------------------
const int nbin = 1E4;                 // # of bins for integration/array filling


//---- Free (user) parameters---------------------------------------------------
// const std::string identifier = "return_power_0D";

class Scatter0D: public Scatter {
public:

  Scatter0D(Antenna& tx, Antenna& rx, Cascade& cs);
  double get_od_cs(){ if (od_cs == 0) {set_od_cs();} return od_cs; }
  double get_ud_cs(){ if (ud_cs == 0) {set_ud_cs();} return ud_cs; }
  double get_power(){ if (power == 0) {set_power();} return power; }
  int accepted()    {set_power(); return triggered;}

private:

  int triggered;
  double od_cs;
  double ud_cs;
  double power;
  void set_od_cs();
  void set_ud_cs();
  void set_power();

  double sinc(double x);

};

int main(int argc, char** argv){

  // I/O Setup
  string cs_filepath = argv[1];
  // Read data from file, get separate unphysical events
  vector<Cascade> Cascades = load_cascade_file(cs_filepath);

  // Make your detector
  string det_type = "bistatic";
  Detector bistatic(det_type);

  // For every cascade
  for (auto& cs: Cascades){
    // For every transmitter
    for (auto& tx : bistatic.get_transmitters()){
      // For every receiver
      for (auto& rx : bistatic.get_receivers()){

        Scatter0D event(tx,rx,cs);

        Antenna new_tx =  event.transmitter();
        Antenna new_rx =  event.receiver();

        // cout  << new_tx.power() << '\t'
        //       << new_tx.gain() << '\t'
        //       << new_rx.gain() << '\t'
        //       << endl;
        //
        // cout  << tx.distance() << '\t'
        //       << rx.distance() << '\t'
        //       << new_tx.distance() << '\t'
        //       << new_rx.distance() << '\t'
        //       << endl;

        cout <<
          "Event:" << scientific << '\t' <<
          "R_t = " << new_tx.distance() << '\t' <<
          "R_r = " << new_rx.distance() << '\t' <<
          "od_tot = " << event.get_od_cs()/1E4 << '\t' <<
          "ud_tot = " << event.get_ud_cs() << '\t' <<
          // "rad_cs = " << rad_cs << '\t' <<
          // "P_r = " << P_r << defaultfloat <<
          "accepted = " << event.accepted() <<
        endl;

        // cout << rad_cs << '\t'  << P_r << '\t' << triggered << endl;

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

      } // Closes antenna loop
    } // Closes transmitter loop
  } // Closes cascade loop
  cout << "END" << endl;
}
// End of main

Scatter0D::Scatter0D(Antenna& tx, Antenna& rx, Cascade& cs):
  Scatter(tx, rx, cs){
    set_od_cs();
    set_ud_cs();
    set_power();
  }

/* Compute the total overdense surface of the cascade */
void Scatter0D::set_od_cs(){
  double X, r, k_mid;       // Variable
  double r_bin  = cs.get_r_bin();
  double r_waist = cs.get_rwaist();

	std::vector<double> r_crit = cs.get_rcrit();
  std::vector<std::vector<double>> density = cs.get_density();

	/*First, get some critical shower values*/
	double dens_max = 0, X_max;
	int i_max;
	// Loop over cascade depth
	for (int i = 0; i < nbin; i++){
		X = i*cs.get_X_bin();
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

  // The shower waist is located by definition at the maximum 	plasma radius.
  double r_waist2 = *max_element(r_crit.begin(), r_crit.end());

	// // If shower age = 1, r_waist is at the maximum density.
	assert(r_waist2 == r_crit[i_max] && "Ill-defined waist radius.");

	/* Second, find the range of layers of "constant" density.
		These layers need to keep the radial spacing r_bin.
	*/
	// double k_mid;              // [unitless]
	double len;										// [cm]
	double dens_layer;            // [g/cm^3]
	double X_low, X_high;					// [g/cm^2]
	std::vector<double> length;
	int k_max = (int) r_waist2/r_bin;

	// Loop over od layers
	for(int k = 0; k < nbin; k++){
		k_mid = k + 0.5;
		// Density of midpoint of k layer at waist
		dens_layer= cs.dens(X_max,k_mid*r_bin);
		X_low = cs.Xtot();
		X_high = 0;
		// Loop over cascade depth
		for (int i = 0; i < nbin; i++){
			X=i*cs.get_X_bin();
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
	norm_in  = 1 - abs(tx.projection());

		// Compute f_diff
	double l_max, gam, av_int = 0, ang, steps = 200000.0;

	// Angle: Distance to relative to full internal reflection
	// "gam" = "alpha" - "beta" in the paper
	gam = acos(tx.projection()) - acos(rx.projection());
// BIG RED FLAG WITH ACOS FUNCTION, IT HAS FAILED IN THE PAST.
// cos(a) = cos(-a) -> acos(cos(-a)) = a

	l_max = *max_element(length.begin(), length.end());
	// Average intensity over the -pi/2 to pi/2 region
	for (int i = 0; i < steps ; i++){
		ang = (-1/2.0 + (i)/steps)*pi; // Some integration value
		av_int += pow(sinc((pi*l_max/tx.lambda())*sin(ang)),2);     // Add contibutions
	}
	av_int /= steps;
	// Intentsity at gamma angle.
	f_diff =  pow(sinc((pi*l_max/tx.lambda())*sin(gam)),2)/av_int;

	f_geometry=abs(norm_in*f_diff);

	od_cs = od_tot* f_geometry;

}

/* Compute under-dense cross-section */
void Scatter0D::set_ud_cs(){
	double N_ch=0, r_bin = cs.get_r_bin();
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
}

// Compute the return power
void Scatter0D::set_power(){
  double Rt, Rr, eta, Latt, rad_cs, P_r, P_bkg;

// Inefficiency --------------------------------------------------------------
  eta = 1;                     // For testing purposes
  // eta=4.4E-3;//*1.d-2          // @450 MHz; 20ns
  //eta=0.18d0*1.d-2              // @50 MHz; 1000ns
  //eta=2.64d-5                   // @10 MHz superdarn detected @ 200-500 MHz

  Latt = 1000;                 // [m] attenuation length
  //Latt=att(freq_obs);        // It should be parametrized depeding on location

  Rt = tx.distance();
  Rr = rx.distance();

  P_r = eta*tx.power()*tx.gain()*rx.gain()*   // Constants
  exp(-2.0*(Rt+Rr)/Latt)/pow(4*pi*Rt*Rr,2)*   // Distance-dependant term
  (od_cs + ud_cs)/1E4;   	                    // [m^2] rad_cs!!
  // Unit change happens here!;

  // Check if the event is over thermal background
  double k_b=1.3806503E-23;       // [SI] Boltzmann's constant
  double T_sys=325;               // [K] system temperature from ARA paper
  double deltaf=1E5;              // [Hz]

  P_bkg=k_b*T_sys*deltaf/1E3;         //[kW]
  (P_r > P_bkg) ? triggered = 1 : triggered = 0;
}

/* Sinc function*/
double Scatter0D::sinc(double x){
  double tmp = 0;
  abs(x) > 10E-8 ? tmp = sin(x)/x : tmp = 1;
  return tmp;
}
