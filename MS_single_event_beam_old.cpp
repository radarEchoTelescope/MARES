/* Bistatic Radar Scatter computation 

*/

#include "settings.hh"
#include "cascade1D.hh"

/* HOW TO MAKE A MARES SCATTER */ 

int main(int argc, char** argv){

  // !!! Include CLHEP-like global units at first definition fo your parameters.  


  // For a single event, the result is stored in the same folder as the executable.
  // std::string path_out = "";
  std::string identifier  = argv[1];
  std::string identifier_c = "Beam1D_" + identifier;





// If you want to save the traces for every time step in the simulation. 
  const bool save_time_profiles = false;

// Other model parameters' defaults can be found in file "settings_params.hh"
// These are optional to set, and override the default. 


// Load a predefined detector configuration (defined in antenna.hh and antenna.cc)
  Detector slac("t576");
  Antenna tx = slac.Transmitters()[0];
  Antenna rx = slac.Receivers()[0];

// Or define your own antennas manually:

  // All the free parameters of the antennas

  // double txxpos    = 0 *m;     // 8
  // double txypos    = 0  *m;      // 9
  // double txzpos    = -2  *m;      // 10
  // double txxpol    = 1;         // 11
  // double txypol    = 0;         // 12
  // double txzpol    = 0;         // 13
  // double txpower   = 50;        // 14
  // double txfreq    = 2  *GHz;    // 15

  // double rxxpos    = 0  *m;       // 16
  // double rxypos    = 0  *m;         // 17
  // double rxzpos    = -2  *m;         // 18
  // double rxxpol    = 1;         // 19
  // double rxypol    = 0;         // 20
  // double rxzpol    = 0;         // 21

  // Then, create the antenna objects. 

  // Antenna transmitter(txxpos, txypos, txzpos,
  //                     txxpol, txypol, txzpol,
  //                     txpower, txfreq);
  // Antenna receiver(   rxxpos, rxypos, rxzpos,
  //                     rxxpol, rxypol, rxzpol);

// Next, define your cascade

  double csenergy   = 10 *GeV;    // Energy per particle
  double csxpos     = 0 *m;       // 3
  double csypos     = 0 *m;       // 4
  double cszpos     = 0 *m;       // 5
  // To put the T576 cascade in front of the TX (sanity check). 
  // double cszpos     = -0.294 *m;       // 5
  double cszenith   = 0;          // 6
  double csazimuth  = 0;          // 7
  double csnumber   = 1E9;          // Number of primaries


// (90,90) = [0,1,0]
// (0,--) = [0,0,1]

  Cascade cascade( 0, csenergy , csxpos, csypos, cszpos,
                  cszenith*deg, csazimuth*deg, csnumber) ;

  

// Create pointers
  // Antenna tx = transmitter;
  // Antenna rx = receiver;

  Cascade& cs = cascade;

// Cascade1D and ToyCascade1D are two types of one-dimensional cascade scatter models. 
// Making a scatter model from the cascade turns the cascade into a vector of scattering points (segments). 
    // std::cout << "What" << std::endl;
  // Cascade1D  nu_cascade(tx,rx,cs, dL, dR, dN, sampling );
  Cascade1D  nu_cascade(tx,rx,cs);

        //               nu_cascade.Length()[150] << '\n' <<
    //               nu_cascade.Radius()[150] << std::endl;
  // ToyCascade1D  nu_cascade(tx,rx,cs, dL, dR, dN, sampling);

// First, choose how to place the scattering points 
// There are 3 possible choices:

  // Along the direction of propagation of the cascade
  // std::cout << nu_cascade.Points()[156].Position << std::endl; 
  nu_cascade.SetInDirection();
  // std::cout << nu_cascade.Points()[156].Position << std::endl; 
  // At the point of maximun reflectivity within the cascade that this segment represents. 
  // nu_cascade.SetInMax()


// Next, choose the propagation mode of the radio waves.
// The options are:
  
  // 1 - Uniform constant medium
  // nu_cascade.SetInConstIce();
    // std::cout << rho_ice << '\n' ;
  // std::cout << nu_cascade.Points()[180].Position << std::endl; 
  
  // 2 - Beam-like setting, where the points and the antennas are
  // in two different media, separated by an interface.

  // The interface is a plane with equation Ax + By + Cz + D = 0
  // na, nb are the refractive indices of the two media at both sides of the interface.

  // The T576 plane is  x = -0.6 m, so
  std::vector<double> interface {1,0,0,0.6*m};
  nu_cascade.SetInBeam(interface, 1, 1.51);

  // 3 - Use Uzair's IceRayTracing for non-constant media. 
   // SetWithIRT();


// Finally, run the scatter. 
  nu_cascade.RunScatter(save_time_profiles);
    // std::cout << "What" << std::endl;


 
  /* Choose what to write out by uncommenting the lines. */
  // write_2D_array(nu_cascade.Coordinates(),    identifier_c + "_coords.txt", 1);
    // write_2D_array(nu_cascade.Length(),          identifier_c + "_length.txt", 1);
    // write_2D_array(nu_cascade.Radius(),          identifier_c + "_radius.txt", 1);

  // write_1D_array(nu_cascade.Phase(),          identifier_c + "_phase.txt", 1);
  write_1D_array(nu_cascade.Attenuation(),    identifier_c + "_attenuation.txt", 1);
  // write_1D_array(nu_cascade.ArrivalTime(),    identifier_c + "_arrival_t.txt", 1);
  write_1D_array(nu_cascade.Directivity(),    identifier_c + "_directivity.txt", 1);
  write_1D_array(nu_cascade.Polarization(),   identifier_c + "_polarization.txt", 1);

  write_1D_array(nu_cascade.Duration(),       identifier_c + "_duration.txt", 1);
  write_1D_array(nu_cascade.Waveform(),       identifier_c + "_waveform.txt", 1);
  // write_1D_array(nu_cascade.Power(),          identifier_c + "_power.txt", 1);
  write_1D_array(nu_cascade.TCS(),            identifier_c + "_TCS.txt", 1);
  // write_1D_array(nu_cascade.RCS(),            identifier_c + "_RCS.txt", 1);
  // write_2D_array(nu_cascade.Phase_time(),     identifier_c + "_phase_time.txt", 1);
  // write_2D_array(nu_cascade.RCS_time(),       identifier_c + "_RCS_time.txt", 1);
  // write_2D_array(nu_cascade.E_time(),         identifier_c + "_E_time.txt", 1);

  //   // Cascade-specific paramters
  // write_2D_array(nu_cascade.Radius(),         identifier_c + "_radial_values.txt", 1);
  write_2D_array(nu_cascade.Ne(),        identifier_c + "_enumber.txt", 1);

  write_2D_array(nu_cascade.Density(),        identifier_c + "_density_tx.txt", 1);
  write_2D_array(nu_cascade.PlasmaFreq(),     identifier_c + "_plasma_freq.txt", 1);
  
  // write_2D_array(nu_cascade.Absorption(),     identifier_c + "_absorption.txt", 1);
  // write_2D_array(nu_cascade.SkinDepth(),      identifier_c + "_skin_depth.txt", 1);
  // write_2D_array(nu_cascade.Reflectance(),    identifier_c + "_reflectance.txt", 1);
  // write_2D_array(nu_cascade.Opacity(),        identifier_c + "_opacity.txt", 1);

}
// End
