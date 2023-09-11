/* Radar scatter computation - return electric field of a list of events
All variables need to be passed through arguments.
Written to be used in a parallel computation (cluster) enviroment.
*/

#include "cascade1D.hh"

// !!! Include CLHEP-like global units at first definition. 
// !!! UNITS ARE NO LONGER ASSUMED. 

int main(int argc, char** argv){
  // Computational parameters ----------------------------------------------

  // For a single event, the result is stored in the same folder as the executable.
  // std::string path_out = "";

  std::string identifier  = argv[1];

  // sampling frequency should be between 10x and 100x freq_obs.
  const double sampling = 100;

  // If you want electric field vs time, enable this. 
  // const bool save_time_profiles = false;

  // Default model variables are in file "settings_params.hh"
  // If set here, they will override the defaults.

// For a single bistatic event, we can use a detector configuration

  Detector slac("t576");

// Or just make two antennas

  // Antenna transmitter(txxpos, txypos, txzpos,
  //                     txxpol, txypol, txzpol,
  //                     txpower, txfreq);
  // Antenna receiver(   rxxpos, rxypos, rxzpos,
  //                     rxxpol, rxypol, rxzpol);

  // Antenna& tx = transmitter;
  // Antenna& rx = receiver;

  // You can create an "empty" scatter .
  Scatter testbeam(slac.Transmitters()[0],slac.Receivers()[0],sampling);

  // This way, you can create and add your own points later. 
  ScatterPoint p;
 
  p.L = 1*cm;
  p.Position = {0,0,2*m};

  // SLAC aligment test, point in front of TX. 
  // p.Position = {0,0,-0.294*m};

  p.TCS = 0.4*pow(m,2);

  testbeam.AddPoint(p);
  // testbeam.AddPoint(p);
  // testbeam.AddPoint(p);
  // testbeam.AddPoint(p);
  // testbeam.AddPoint(p);

  // testbeam.AddPoint(p);
  // testbeam.AddPoint(p);
  // testbeam.AddPoint(p);
  // testbeam.AddPoint(p);
  // testbeam.AddPoint(p);

// Next, choose the propagation mode: 
  
  // Uniform constant medium
  testbeam.SetInConstIce();
  std::cout << testbeam.Points()[1].Position<< std::endl; 

  
  // Beam-like setting, where the points and the antennas are
  // in two different media, separated by an interface.
  
  // The interface is a plane with equation Ax + By + Cz + D = 0
  // na, nb are the refractive indices of the two media at both sides of the interface.

  // The T576 plane is  x = -0.6 m, so
  // std::vector<double> interface {1,0,0,0.6*m};
  // testbeam.SetInBeam(interface, 1, 1.51);

  // Use Uzair's IceRayTracing for non-constant media. 
   // SetWithIRT();


// Finally, run the scatter. 
  testbeam.RunScatter();

  std::string identifier_c = "SinglePoint_" + identifier;

 
  /* Choose what to write out by uncommenting the lines. */
  // write_2D_array(testbeam.Coordinates(),    identifier_c + "_coords.txt", 1);
  write_1D_array(testbeam.Phase(),          identifier_c + "_phase.txt", 1);
  write_1D_array(testbeam.Attenuation(),    identifier_c + "_attenuation.txt", 1);
  write_1D_array(testbeam.ArrivalTime(),    identifier_c + "_arrival_t.txt", 1);
  write_1D_array(testbeam.Directivity(),    identifier_c + "_directivity.txt", 1);
  write_1D_array(testbeam.Polarization(),   identifier_c + "_polarization.txt", 1);

  write_1D_array(testbeam.Duration(),       identifier_c + "_duration.txt", 1);
  write_1D_array(testbeam.Waveform(),       identifier_c + "_waveform.txt", 1);
  // write_1D_array(testbeam.Power(),          identifier_c + "_power.txt", 1);
  write_1D_array(testbeam.TCS(),            identifier_c + "_TCS.txt", 1);
  // write_1D_array(testbeam.RCS(),            identifier_c + "_RCS.txt", 1);

  // write_2D_array(testbeam.Phase_time(),     identifier_c + "_phase_time.txt", 1);
  // write_2D_array(testbeam.RCS_time(),       identifier_c + "_RCS_time.txt", 1);
  // write_2D_array(testbeam.E_time(),         identifier_c + "_E_time.txt", 1);

}
// End
