/* Bistatic Radar Scatter executable for MARES
This executable takes most input variables needed for the event through arguments.
Useful in a cluster enviroment to run events in parallel. 
*/

#include "cascade1D.hh"

int main(int argc, char** argv){
  
/* PASSED BY ARGUMENTS */

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

  double csxpos     = atof(argv[2]) *m;
  double csypos     = atof(argv[3]) *m;
  double rxxpos     = atof(argv[4])  *m;
  double rxypos     = atof(argv[5])  *m;
  double csenergy   = atof(argv[6]) *PeV;
  double cszenith   = atof(argv[7]) *deg;
  double csazimuth  = atof(argv[8]) *deg;
  double txfreq     = atof(argv[9]) *GHz;
  double txgaindB   = atof(argv[10]);
  double txzpos     = atof(argv[11])  *m;
  double rxzpos     = atof(argv[12])  *m;

  double cszpos     = atof(argv[13]) *m;
  double csnumber   = atof(argv[14]);

  // Now we make the cascade's object
  Cascade cascade( csxpos, csypos, cszpos, cszenith, csazimuth, csenergy, csnumber);
  // And pointer. 
  Cascade& cs = cascade;

  // 3 - Antennas
  // There are pre-defined detector configurations available in antenna.hh and antenna.cc
  // We assume that you want to make the antennas too 
  
    // Transmitter, TX
  double txxpos   = atof(argv[15])   *m;
  double txypos   = atof(argv[16])   *m;
  double txxpol   = atof(argv[17]);
  double txypol   = atof(argv[18]);
  double txzpol   = atof(argv[19]);
  double txpower  = atof(argv[20]);

  Antenna transmitter(txxpos, txypos, txzpos,
                      txxpol, txypol, txzpol,
                      txpower, txfreq, txgaindB);
  Antenna& tx = transmitter;

  double rxxpol   = atof(argv[21]);
  double rxypol   = atof(argv[22]);
  double rxzpol   = atof(argv[23]);
  // The frequency is used in the receiverto determine the antenna's effective area
  // So far, we have taken rx and tx to operate at the same freq.
  // double rxfreq   = atof(argv[24])   *GHz;
  double rxfreq = txfreq;
  // double rxgaindB = atof(argv[25]);
  double rxgaindB = txgaindB;

  Antenna receiver(   rxxpos, rxypos, rxzpos,
                      rxxpol, rxypol, rxzpol,
                      0 , txfreq, rxgaindB);    // tx frequency not rx frequency
  Antenna& rx = receiver;

  // 4 - With TX, RX and CS, we make a Cascade1D object.
  // This generates the density frame and a list of scattering points. 
  Cascade1D nu_cascade(tx,rx,cs);

  // We are using the default simulation resolution variables. 
  // You can specify your own at the moment of construction. 
  // Cascade1D nu_cascade(tx,rx,cs, dL, dR, dN, sampling);

  // The defaults for resolution variables and other model 
  // parameters are in file "settings_params.hh"

  // 5 - Choose how to place the segments. 

    // Along the direction of propagation of the cascade.
  nu_cascade.SetInDirection();

    // Set the scatterer in the coordinates of maximum reflectivity for that segment.
  // nu_cascade.SetInMax();

  // 6 - Choose the propagation mode.

    // Constant and uniform medium of density n, large-scale attenuation given by att_length. 
  nu_cascade.SetInConstMedium();

    // Use IceRayTracing to propagate through non-constant, realistic ice media. 
  
  /* WARNING */
  /* 
      This propagation mode has been implemented but not tested. 
      Use at your own risk and apply appropiate sanity checks.
   */

  // nu_cascade.SetWithIRT();

    // Propagate with an interface. See MS_single_event_beam.cpp.  
  // nu_cascade.SetInBeam();

  // 7 - Run the scatter proper. 

    // If you want to save the individual traces for every time step in the simulation. 
  const bool save_time_profiles = false;      //default MARES
  // const bool save_time_profiles = true;

  nu_cascade.RunScatter(save_time_profiles);

  // 8 - Choose to save to disk.
  
    /*
      There is a 1 at the end that acts as a flag. 
      0 will not write the array.
    */
  
      // Cascade (cs) products.
    /* Some of these will only be enerated on demand. They are not needed for the scatter.*/

  // write_2D_array(nu_cascade.Ne(),             identifier_c + "_e_number.txt", 1); //****      // cascade

  //write_2D_array(nu_cascade.Radius(),         identifier_c + "_radial_values.txt", 1);        // cascade1D  ?
  // CSlength?
  //write_2D_array(nu_cascade.Density(),        identifier_c + "_density_tx.txt", 1);           // cascade
  // write_2D_array(nu_cascade.PlasmaFreq(),     identifier_c + "_plasma_freq.txt", 1);         // cascade
  // write_2D_array(nu_cascade.Absorption(),     identifier_c + "_absorption.txt", 1);          // cascade
  // write_2D_array(nu_cascade.SkinDepth(),      identifier_c + "_skin_depth.txt", 1);          // cascade
  // write_2D_array(nu_cascade.Reflectance(),    identifier_c + "_reflectance.txt", 1);         // cascade
  // write_2D_array(nu_cascade.Opacity(),        identifier_c + "_opacity.txt", 1);             // cascade
  // transparency?                                                                              // cascade

      // Scatter object products: Positions and propagation 
  // write_2D_array(nu_cascade.Position(),          identifier_c + "_position.txt",1);  //**** 

  // write_2D_array(nu_cascade.Coordinates(),    identifier_c + "_coords.txt", 1);
  // write_1D_array(nu_cascade.Phase(),          identifier_c + "_phase.txt", 1);
  // write_1D_array(nu_cascade.Attenuation(),    identifier_c + "_attenuation.txt", 1);
  // write_1D_array(nu_cascade.ArrivalTime(),    identifier_c + "_arrival_t.txt", 1);
  // write_1D_array(nu_cascade.Directivity(),    identifier_c + "_directivity.txt", 1);
  // write_1D_array(nu_cascade.Geometry(),       identifier_c + "_geom_efficiency.txt", 1);

      // Scatter products: Time-dependent variables. 
  write_1D_array(nu_cascade.Duration(),       identifier_c + "_duration.txt", 1);
  write_1D_array(nu_cascade.Voltage(),       identifier_c + "_voltage.txt", 1);
  // write_1D_array(nu_cascade.Power(),          identifier_c + "_power.txt", 1);
  // write_1D_array(nu_cascade.TCS(),            identifier_c + "_TCS.txt", 1);               // cascade1D
  // write_1D_array(nu_cascade.RCS(),            identifier_c + "_RCS.txt", 1);

      // These are the time_profiles. 
  // write_2D_array(nu_cascade.Phase_time(),     identifier_c + "_phase_time.txt", 1);
  // write_2D_array(nu_cascade.RCS_time(),       identifier_c + "_RCS_time.txt", 1);
  // write_2D_array(nu_cascade.E_time(),         identifier_c + "_E_time.txt", 1);

}
// End
