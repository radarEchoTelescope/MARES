/* Bistatic Radar Scatter executable for MARES
This executable needs a single string as an argument to run
E.g.,

./MS_single_point point_example

THIS IS A BASIC EXAMPLE FOR A SINGLE POINT IN SPACE. 
For a full cascade simulation, see MS_single_event.cpp
*/

#include "cascade1D.hh"

int main(int argc, char** argv){

  // 1- Unique event identifier -----------------------------------
  /* This is the name of this event, used in the output files. 
  These files are produced in the same folder as the executable.
  Usually the cluster takes care of the final output directory. 
  */
  std::string identifier  = argv[1];
 
 // Full identifier (in case you have similar executables). 
  std::string identifier_c = "SinglePoint_" + identifier;

  // --------------------------------------------------------------- 
  // 2 - Antennas
  // There are pre-defined detector configurations available in antenna.hh and antenna.cc
  // You can use a detector configuration and pick a TX and a RX. 

  Detector dect("bistatic");
  Antenna tx = dect.Transmitters()[0];
  Antenna rx = dect.Receivers()[0];
  
  // Or define your own antennas manually:
  
    // Transmitter, TX
  // double txxpos   = 0   *m;
  // double txypos   = 0   *m;
  // double txzpos   = -100  *m;
  // double txxpol   = 0;
  // double txypol   = 1;
  // double txzpol   = 0;
  // double txpower  = 1E3;
  // double txfreq   = 50 *Hz;
  // double txgaindB = 0;
  // // Half-dipole gaindB = 2.15 dBi

  // Antenna transmitter(txxpos, txypos, txzpos,
  //                     txxpol, txypol, txzpol,
  //                     txpower, txfreq, txgaindB);
  // Antenna& tx = transmitter;

  // double rxxpos   = -250  *m;
  // double rxypos   = 0  *m;
  // double rxzpos   = 0  *m;
  // double rxxpol   = 0;
  // double rxypol   = 1;
  // double rxzpol   = 0;

  // // The frequency is used to determine the antenna's effective area
  // // So far, we have taken rx and tx to operate at the same freq.
  
  // double rxfreq = txfreq;
  // double rxgaindB = 0;

  // Antenna receiver(   rxxpos, rxypos, rxzpos,
  //                     rxxpol, rxypol, rxzpol,
  //                     0 , txfreq, rxgaindB);
  // Antenna& rx = receiver;

  // ------------------------------------------------------------
  // 3 - We create a generic scatter object devoid of scattering points. 
  
  /* 
    The cascade classes that are (or can be) implemented later  
    - so far, only Cascade1D - are derived from this scatter class.

    If you find another way to represent your case as a collection of points,
    Scatter will run the event. 
  */

  Scatter test_points(tx,rx);
  
  // We can specify here the sampling ratio for the time-dependent part
  // of the simulation (generating the voltage waveform). 
  // It represents how much higer is the sampling frequency w.r.t TX's freq.
  // The sampling ratio should be between 10 and 100.

  // The defaults for the model variables in the file "src/settings_params.hh"

  // const double sampling = 100;
  // Scatter test_points(tx,rx,sampling);

  // Let's create an example point now. 
  // The information about the scattering element is stored in the ScatterPoint struct. 

  ScatterPoint p;

    // Add position
  p.Position = {250,0,-250*m};
    // Physical dimensions. 
  p.L = 1*cm;
    // Cross section 0.5 m^2
  p.TCS = 0.5*pow(m,2);

  test_points.AddPoint(p);
  // ------------------------------------------------------------

  // 4 - Choose the propagation mode.

    // Constant and uniform medium of density n, large-scale attenuation given by att_length. 
  test_points.SetInConstMedium();

// OR
    // Use IceRayTracing to propagate through non-constant, realistic ice media. 
  
  /* WARNING */
  /* 
      This propagation mode has been implemented but not tested. 
      Use at your own risk and apply appropiate sanity checks.
   */

  // test_points.SetWithIRT();

// OR
    // Propagate with an interface. See MS_single_event_beam.cpp.  
  // test_points.SetInBeam();

  // 7 - Run the scatter proper. 

    // If you want to save the individual traces for every time step in the simulation. 
  const bool save_time_profiles = false;

  test_points.RunScatter(save_time_profiles);

  // 8 - Choose to save to disk.
  
    /*
      There is a 1 at the end that acts as a flag. 
      0 will not write the array.
    */
  
      
      // Scatter object products: Positions and propagation 
  // write_2D_array(test_points.Coordinates(),    identifier_c + "_coords.txt", 1);
  // write_1D_array(test_points.Phase(),          identifier_c + "_phase.txt", 1);
  // write_1D_array(test_points.Attenuation(),    identifier_c + "_attenuation.txt", 1);
  // write_1D_array(test_points.ArrivalTime(),    identifier_c + "_arrival_t.txt", 1);
  // write_1D_array(test_points.Directivity(),    identifier_c + "_directivity.txt", 1);
  // write_1D_array(test_points.Geometry(),   identifier_c + "_geom_efficiency.txt", 1);

      // Scatter products: Time-dependent variables. 
  write_1D_array(test_points.Duration(),       identifier_c + "_duration.txt", 1);
  write_1D_array(test_points.Voltage(),       identifier_c + "_voltage.txt", 1);
  // write_1D_array(test_points.Power(),          identifier_c + "_power.txt", 1);
  // write_1D_array(test_points.TCS(),            identifier_c + "_TCS.txt", 1);
  // write_1D_array(test_points.RCS(),            identifier_c + "_RCS.txt", 1);

      // These are the time_profiles. 
  // write_2D_array(test_points.Phase_time(),     identifier_c + "_phase_time.txt", 1);
  // write_2D_array(test_points.RCS_time(),       identifier_c + "_RCS_time.txt", 1);
  // write_2D_array(test_points.E_time(),         identifier_c + "_E_time.txt", 1);

}
// End
