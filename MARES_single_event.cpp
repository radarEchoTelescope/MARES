/* Bistatic Radar Scatter executable for MARES
This executable needs a single string as an argument to run
E.g.,

./MS_single_event test_run
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
  std::string identifier_c = "Cascade1D_" + identifier;

  // 2- Cascade's parameters 
  // [!NOTICE] Units are not assumed in MARES.
  // Always include CLHEP-like global units at first definition. 
  // List of available units at settings.hh

  double csxpos     = 0   *m;       
  double csypos     = 0   *m;
  double cszpos     = 0   *m;
  double cszenith   = 90  *deg;
  double csazimuth  = 90  *deg;
    // Energy per particle
  double csenergy   = 1E7  *GeV;    
    // Number of primaries
  double csnumber   = 1;

  /* CASCADE DIRECTION REMINDER
    (90,90) = [0,1,0]
    (0,--) = [0,0,1]
  */

  // Now we make the cascade's object
  Cascade cascade( csxpos, csypos, cszpos, 
                   cszenith, csazimuth, csenergy, csnumber);
  // And its pointer. 
  Cascade& cs = cascade;

  // 3 - Antennas
  // There are pre-defined detector configurations available in antenna.hh and antenna.cc
  // You can use a detector configuration and pick a TX and a RX. 

  // Detector dect("bistatic");
  // Antenna tx = dect.Transmitters()[0];
  // Antenna rx = dect.Receivers()[0];
  
  // Or define your own antennas manually:
  
    // Transmitter, TX
  double txxpos   = 0   *m;
  double txypos   = 0   *m;
  double txzpos   = -100  *m;
  double txxpol   = 0;
  double txypol   = 1;
  double txzpol   = 0;
  double txpower  = 1E3;
  double txfreq   = 50 *MHz;
  double txgaindB = 0;
  // Half-dipole gaindB = 2.15 dBi

  Antenna transmitter(txxpos, txypos, txzpos,
                      txxpol, txypol, txzpol,
                      txpower, txfreq, txgaindB);
  Antenna& tx = transmitter;

  double rxxpos   = -250  *m;
  double rxypos   = 0  *m;
  double rxzpos   = 0  *m;
  double rxxpol   = 0;
  double rxypol   = 1;
  double rxzpol   = 0;

  // The frequency is used to determine the antenna's effective area
  // So far, we have taken rx and tx to operate at the same freq.
  
  double rxfreq = txfreq;
  double rxgaindB = 0;

  Antenna receiver(   rxxpos, rxypos, rxzpos,
                      rxxpol, rxypol, rxzpol,
                      0 , txfreq, rxgaindB);
  Antenna& rx = receiver;

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

  // 5 - Choose how to place the segments. 

    // Along the direction of propagation of the cascade.
  nu_cascade.SetInDirection();

// OR

    // Set the scatterer in the coordinates of maximum reflectivity for that segment.
  // nu_cascade.SetInMax();

  // 6 - Choose the propagation mode.

    // Constant and uniform medium of density n, large-scale attenuation given by att_length. 
  nu_cascade.SetInConstMedium();

// OR
    // Use IceRayTracing to propagate through non-constant, realistic ice media. 
  
  /* WARNING */
  /* 
      This propagation mode has been implemented but not tested. 
      Use at your own risk and apply appropiate sanity checks.
   */

  // nu_cascade.SetWithIRT();

// OR
    // Propagate with an interface. See MS_single_event_beam.cpp.  
  // nu_cascade.SetInBeam();

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
 

      // Scatter object products: Positions and propagation 
  // write_2D_array(nu_cascade.Coordinates(),    identifier_c + "_coords.txt", 1);
  // write_1D_array(nu_cascade.Phase(),          identifier_c + "_phase.txt", 1);
  // write_1D_array(nu_cascade.Attenuation(),    identifier_c + "_attenuation.txt", 1);
  // write_1D_array(nu_cascade.ArrivalTime(),    identifier_c + "_arrival_t.txt", 1);
  // write_1D_array(nu_cascade.Directivity(),    identifier_c + "_directivity.txt", 1);
  // write_1D_array(nu_cascade.Geometry(),       identifier_c + "_geom_efficiency.txt", 1);

      // Scatter products: Time-dependent variables. 
  write_1D_array(nu_cascade.Duration(),       identifier_c + "_duration.txt", 1);
  write_1D_array(nu_cascade.Voltage(),        identifier_c + "_voltage.txt", 1);
  // write_1D_array(nu_cascade.Power(),          identifier_c + "_power.txt", 1);
  // write_1D_array(nu_cascade.TCS(),            identifier_c + "_TCS.txt", 1);
  // write_1D_array(nu_cascade.RCS(),            identifier_c + "_RCS.txt", 1);

      // These are the time_profiles. 
  // write_2D_array(nu_cascade.Phase_time(),     identifier_c + "_phase_time.txt", 1);
  // write_2D_array(nu_cascade.RCS_time(),       identifier_c + "_RCS_time.txt", 1);
  // write_2D_array(nu_cascade.E_time(),         identifier_c + "_E_time.txt", 1);

}
// End
