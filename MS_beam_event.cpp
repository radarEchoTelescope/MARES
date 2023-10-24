/* Bistatic Radar Scatter executable for MARES
This executable needs a single string as an argument to run
E.g.,

./MS_beam_event test_run

This executable was designed to replice a radar detection within a 
test beam, like the SLAC T576 experiment, but it works for any interface
between two media. 

Here, as an example, the T576 paramters are used. 
*/

#include "cascade1D.hh"

int main(int argc, char** argv){

  // 1- Unique event identifier
  /* This is the name of this event, used in the output files. 
  These files are produced in the same folder as the executable.
  Usually the cluster takes care of the final output directory. 
  */
  std::string identifier  = argv[1];
 
 // Full identifier (in case you have similar executables). 
  std::string identifier_c = "Beam1D_" + identifier;

  // 2- Cascade's parameters 
  // [!NOTICE] Units are not assumed in MARES.
  // Always include CLHEP-like global units at first definition. 
  // List of available units at settings.hh

    // Energy per particle
  double csenergy   = 10  *GeV;    
  double csxpos     = 0   *m;       
  double csypos     = 0   *m;
  double cszpos     = 0   *m;
  // To put the T576 cascade in front of the TX (sanity check). 
  // double cszpos     = -0.294 *m;       // 5
  double cszenith   = 0 *deg;
  double csazimuth  = 0 *deg;
    // Number of primaries
  // A beam usually starts with a bunch of electrons 
  // That's why the energy above is defined per primary particle.
  double csnumber   = 1E9;

  // Now we make the cascade's object
  Cascade cascade( 0, csenergy, csxpos, csypos, cszpos, 
                      cszenith, csazimuth, csnumber);
  // And its pointer. 
  Cascade& cs = cascade;

  // 3 - Antennas
  // There are pre-defined detector configurations available in antenna.hh and antenna.cc
  // You can use a detector configuration and pick a TX and a RX. 

  Detector slac("t576");
  Antenna tx = slac.Transmitters()[0];
  // Using only one receiver antenna
  Antenna rx = slac.Receivers()[0];

  // 4 - With TX, RX and CS, we make a Cascade1D object.
  // This generates the density frame and a list of scattering points. 
  Cascade1D nu_cascade(tx,rx,cs);

  // We are using the default simulation resolution variables. 
  // You can specify your own at the moment of construction. 

  // The spatial resolution is the effective scatterer size. 
  // const double dL = 1 *cm; // [cm/bin]
  // const double dR = 1 *mm;
  // const double dN = 1 *mm; // This is also a radial direction.

  // Sampling ratio, should be between 10x and 100x freq_obs.
  // const double sampling = 50;

  // Cascade1D nu_cascade(tx,rx,cs, dL, dR, dN, sampling);

  // The defaults for the resolution variables and other model 
  // parameters are in file "settings_params.hh"

  // The T576 beam was assumed to have a free electron lifetime of 3 ns, 
  // this value needs to be changed in settings_params.cc

  // 5 - Choose how to place the segments. 

    // Along the direction of propagation of the cascade.
  nu_cascade.SetInDirection();

// OR

    // Set the scatterer in the coordinates of maximum reflectivity for that segment.
  // nu_cascade.SetInMax();

  // 6 - Choose the propagation mode.

  /* Propagate with an interface. 
  The scattering points and the antennas are in two different media.

  The interface between the media is a plane with equation Ax + By + Cz + D = 0
  na = 1, nb =1.51 are the refractive indices of the two media at both sides of the interface.
  */

  // The T576 plane is  x = -0.6 m, so
  std::vector<double> interface {1,0,0,0.6*m};
  nu_cascade.SetInBeam(interface, 1, 1.51);

  /* IMPORANT! If you want to disable the Fresnel coefficients effects while 
  running in beam mode, you need to edit out manually one line in scatter.cc
  
  -------------------------------------
                              p.PolEff * p.Attenuation;
  -------------------------------------

  Right now, it is line 649 in scatter.cc

  You need to re-enable it if you want to run a different mode.

  TO-DO: Make this a flag-enabled option.
  */

  // 7 - Run the scatter proper. 

    // If you want to save the individual traces for every time step in the simulation. 
  const bool save_time_profiles = false;
  nu_cascade.RunScatter(save_time_profiles);

  // 8 - Choose to save to disk.
  
    /*
      There is a 1 at the end that acts as a flag. 
      0 will not write the array.
    */
  
      // Cascade (cs) products.
    /* Some of these will only be enerated on demand. They are not needed for the scatter.*/

  // write_2D_array(nu_cascade.Radius(),         identifier_c + "_radial_values.txt", 1);
  // write_2D_array(nu_cascade.Ne(),             identifier_c + "_e_number.txt", 1);
  // write_2D_array(nu_cascade.Density(),        identifier_c + "_density_tx.txt", 1);
  // write_2D_array(nu_cascade.PlasmaFreq(),     identifier_c + "_plasma_freq.txt", 1);
  // write_2D_array(nu_cascade.Absorption(),     identifier_c + "_absorption.txt", 1);
  // write_2D_array(nu_cascade.SkinDepth(),      identifier_c + "_skin_depth.txt", 1);
  // write_2D_array(nu_cascade.Reflectance(),    identifier_c + "_reflectance.txt", 1);
  // write_2D_array(nu_cascade.Opacity(),        identifier_c + "_opacity.txt", 1);

      // Scatter object products: Positions and propagation 
  // write_2D_array(nu_cascade.Coordinates(),    identifier_c + "_coords.txt", 1);
  // write_1D_array(nu_cascade.Phase(),          identifier_c + "_phase.txt", 1);
  // write_1D_array(nu_cascade.Attenuation(),    identifier_c + "_attenuation.txt", 1);
  // write_1D_array(nu_cascade.ArrivalTime(),    identifier_c + "_arrival_t.txt", 1);
  // write_1D_array(nu_cascade.Directivity(),    identifier_c + "_directivity.txt", 1);
  // write_1D_array(nu_cascade.Polarization(),   identifier_c + "_polarization.txt", 1);

      // Scatter products: Time-dependent variables. 
  write_1D_array(nu_cascade.Duration(),       identifier_c + "_duration.txt", 1);
  write_1D_array(nu_cascade.Voltage(),       identifier_c + "_voltage.txt", 1);
  // write_1D_array(nu_cascade.Power(),          identifier_c + "_power.txt", 1);
  write_1D_array(nu_cascade.TCS(),            identifier_c + "_TCS.txt", 1);
  write_1D_array(nu_cascade.RCS(),            identifier_c + "_RCS.txt", 1);

      // These are the time_profiles. 
  // write_2D_array(nu_cascade.Phase_time(),     identifier_c + "_phase_time.txt", 1);
  // write_2D_array(nu_cascade.RCS_time(),       identifier_c + "_RCS_time.txt", 1);
  // write_2D_array(nu_cascade.E_time(),         identifier_c + "_E_time.txt", 1);

}
// End
