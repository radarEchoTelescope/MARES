/* Radar scatter computation - return electric field of a list of events
All variables need to be passed through arguments.
Written to be used in a parallel computation (cluster) enviroment.
*/

#include "cascade1D.hh"

// !!! Include CLHEP-like global units at first definition. 

int main(int argc, char** argv){
  // Computational parameters ----------------------------------------------

  // For a single event, the result is stored in the same folder as the executable.
  // std::string path_out = "";

  std::string identifier  = argv[1];
  std::string identifier_c = "Cascade1D_" + identifier;

  const double dL = 1   *cm; // [cm/bin]
  const double dR = 0.1 *cm;
  const double dN = 1   *cm; // This is also a radial direction.

  const double sampling = 100;
  // sampling frequency should be between 10x and 100x freq_obs.

  // If you want to save the traces for every time step in the simulation. 
  const bool save_time_profiles = false;

  // Other model variables are in file "settings_params.hh"
  // These are optional to set. 


// NOW WE HAVE TO PASS BY UNITS WHEN PASSING THE VALUES, UNITS ARE NO LONGER ASSUMED. 

/* VALUES FOR REFERENCE AND FOR TEST */
//   double csenergy   = 10 *PeV;  // 2
//   double csxpos     = 0 *m;        // 3
//   double csypos     = 0 *m;       // 4
//   double cszpos     = 0 *m;        // 5
//   double cszenith   = 90;       // 6
//   double csazimuth  = 90;       // 7
// // (90,90) = [0,1,0]

//   double txxpos    = -2 *m;     // 8
//   double txypos    = 6  *m;      // 9
//   double txzpos    = 0  *m;      // 10
//   double txxpol    = 1;         // 11
//   double txypol    = 0;         // 12
//   double txzpol    = 0;         // 13
//   double txpower   = 50;        // 14
//   double txfreq    = 2  *GHz;    // 15

//   double rxxpos    = 4  *m;       // 16
//   double rxypos    = 6  *m;         // 17
//   double rxzpos    = 0  *m;         // 18
//   double rxxpol    = 1;         // 19
//   double rxypol    = 0;         // 20
//   double rxzpol    = 0;         // 21
  
/* ACTUAL PASSED BY VALUES */

  double csenergy = atof(argv[2]);
  double csxpos = atof(argv[3]);
  double csypos = atof(argv[4]);
  double cszpos = atof(argv[5]);
  double cszenith = atof(argv[6]);
  double csazimuth = atof(argv[7]);
  
  double txxpos = atof(argv[8]);
  double txypos = atof(argv[9]);
  double txzpos = atof(argv[10]);
  double txxpol = atof(argv[11]);
  double txypol = atof(argv[12]);
  double txzpol = atof(argv[13]);
  double txpower = atof(argv[14]);
  double txfreq = atof(argv[15]);
  
  double rxxpos = atof(argv[16]);
  double rxypos = atof(argv[17]);
  double rxzpos = atof(argv[18]);
  double rxxpol = atof(argv[19]);
  double rxypol = atof(argv[20]);
  double rxzpol = atof(argv[21]);

  // double lifetime = atof(argv[22]);

  // For a single bistatic event, no need to load a detector, just make two antennas

  Antenna transmitter(txxpos *m, txypos *m, txzpos *m,
                      txxpol, txypol, txzpol,
                      txpower, txfreq *Hz);
  Antenna receiver(   rxxpos *m, rxypos *m, rxzpos *m,
                      rxxpol, rxypol, rxzpol);

  Cascade cascade( 0, csenergy * GeV , csxpos *m, csypos *m, cszpos *m,
                  cszenith*deg, csazimuth*deg) ;

  Antenna& tx = transmitter;
  Antenna& rx = receiver;
  Cascade& cs = cascade;


// Making a cascade object generates a list of scattering points. 
  Cascade1D     nu_cascade(tx,rx,cs, dL, dR, dN, sampling);
  // ToyCascade1D  nu_cascade(tx,rx,cs, dL, dR, dN, sampling);

  // Scatter testbeam(slac.Transmitters()[0],slac.Receivers()[0],sampling);
  // ScatterPoint p;
  // p.Position = {0,0,2*m};
  // p.L = 1*cm;
  // p.TCS = 1;

  // testbeam.AddPoint(p);



// First, choose how to place the segments. 
// There are 3 possible choices
  nu_cascade.SetInDirection();
  // SetInPosition()
  // SetInMax()

// Second, choose the propagation mode: 
  nu_cascade.SetInConstIce();
  // SetWithIRT();
  // SetInBeam();
  // The interface is a plane with equation Ax + By + Cz + D = 0

  // The T576 plane is  x = -0.6 m, so
  // std::vector<double> interface {1,0,0,0.6*m};
  // The other two numbers are the refractive indices at the two sides of the plane. 
  // testbeam.SetInBeam(interface, 1, 1.51);


// Finally, run the scatter. 
  // testbeam.RunScatter();
  nu_cascade.RunScatter(save_time_profiles);

 
  /* Choose what to write out by uncommenting the lines. */
  // write_2D_array(nu_cascade.Coordinates(),    identifier_c + "_coords.txt", 1);
  write_1D_array(nu_cascade.Phase(),          identifier_c + "_phase.txt", 1);
  write_1D_array(nu_cascade.Attenuation(),    identifier_c + "_attenuation.txt", 1);
  write_1D_array(nu_cascade.ArrivalTime(),    identifier_c + "_arrival_t.txt", 1);
  write_1D_array(nu_cascade.Directivity(),    identifier_c + "_directivity.txt", 1);
  write_1D_array(nu_cascade.Polarization(),   identifier_c + "_polarization.txt", 1);

  write_1D_array(nu_cascade.Duration(),       identifier_c + "_duration.txt", 1);
  write_1D_array(nu_cascade.Waveform(),       identifier_c + "_waveform.txt", 1);
  // write_1D_array(nu_cascade.Power(),          identifier_c + "_power.txt", 1);
  write_1D_array(nu_cascade.TCS(),            identifier_c + "_TCS.txt", 1);
  write_1D_array(nu_cascade.RCS(),            identifier_c + "_RCS.txt", 1);
  // write_2D_array(nu_cascade.Phase_time(),     identifier_c + "_phase_time.txt", 1);
  // write_2D_array(nu_cascade.RCS_time(),       identifier_c + "_RCS_time.txt", 1);
  // write_2D_array(nu_cascade.E_time(),         identifier_c + "_E_time.txt", 1);

  //   // Cascade-specific paramters
  // write_2D_array(nu_cascade.Radius(),         identifier_c + "_radial_values.txt", 1);
  // write_2D_array(nu_cascade.Density(),        identifier_c + "_density_tx.txt", 1);
  // write_2D_array(nu_cascade.PlasmaFreq(),     identifier_c + "_plasma_freq.txt", 1);
  // write_2D_array(nu_cascade.Absorption(),     identifier_c + "_absorption.txt", 1);
  // write_2D_array(nu_cascade.SkinDepth(),      identifier_c + "_skin_depth.txt", 1);
  // write_2D_array(nu_cascade.Reflectance(),    identifier_c + "_reflectance.txt", 1);
  // write_2D_array(nu_cascade.Opacity(),        identifier_c + "_opacity.txt", 1);

}
// End
