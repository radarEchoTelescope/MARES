/* Radar scatter computation - return electric field of a list of events
All variables need to be passed through arguments.
Written to be used in a parallel computation (cluster) enviroment.
*/

#include "macro_scatter.hh"
int main(int argc, char** argv){

  // I/O Setup
  // For a single event, the result is stored in the same folder as the executable.
  std::string path_out = "";

  // Cascade parameters

  std::string identifier  = argv[1];

  double cenergy = atof(argv[2]);
  double xpos = atof(argv[3]);
  double ypos = atof(argv[4]);
  double zpos = atof(argv[5]);
  double czenith = atof(argv[6]);
  double cazimuth = atof(argv[7]);

  double txpower = atof(argv[8]);
  double txxpos = atof(argv[9]);
  double txypos = atof(argv[10]);
  double txzpos = atof(argv[11]);

  double rxxpos = atof(argv[12]);
  double rxypos = atof(argv[13]);
  double rxzpos = atof(argv[14]);

/* REFERENCE VALUES FOR A SIMPLE TEST*/

  // double cenergy   = 1E9;
  // double xpos      = 0;
  // double ypos      = -250 ;
  // double zpos      = 0;
  // double czenith   = 90;
  // double cazimuth  = 45;

  // double txpower  = 1E4;
  // double txxpos    = 0;
  // double txypos    = 0;
  // double txzpos    = 0;

  // double rxxpos    = 500;
  // double rxypos    = 0;
  // double rxzpos    = 0;

  czenith = deg2rad(czenith);
  cazimuth = deg2rad(cazimuth);


  // No need to load a detector, just make two antennas
  Antenna transmitter(txpower, txxpos, txypos, txzpos);
  Antenna receiver(0, rxxpos, rxypos, rxzpos);

  Cascade cascade( 0, cenergy , xpos, ypos, zpos, czenith, cazimuth) ;

  Antenna& tx = transmitter;
  Antenna& rx = receiver;
  Cascade& cs = cascade;

  std::string identifier_c = path_out + "Cascade1D_" + identifier;
  Cascade1D nu_cascade(tx,rx,cs);

  /* Choose what to write out by uncommenting the lines. */
  write_2D_array(nu_cascade.Coordinates(),    identifier_c + "_coords.txt", 1);
  // write_1D_array(nu_cascade.Phase(),       identifier_c + "_phases.txt", 1);
  write_1D_array(nu_cascade.Attenuation(),    identifier_c + "_attenuation.txt", 1);
  // write_1D_array(nu_cascade.ArrivalTime(),    identifier_c + "_arrival_times.txt", 1);
  write_1D_array(nu_cascade.Directivity(),     identifier_c + "_directivity.txt", 1);
  write_1D_array(nu_cascade.Polarization(),   identifier_c + "_polarization.txt", 1);

  write_1D_array(nu_cascade.Duration(),       identifier_c + "_duration.txt", 1);
  write_1D_array(nu_cascade.Waveform(),       identifier_c + "_waveform.txt", 1);
  write_1D_array(nu_cascade.Power(),         identifier_c + "_power.txt", 1);
  write_1D_array(nu_cascade.TCS(),            identifier_c + "_TCS.txt", 1);
  write_1D_array(nu_cascade.RCS(),            identifier_c + "_RCS.txt", 1);
  // write_2D_array(nu_cascade.Phase_time(),  identifier_c + "_phase_time.txt", 1);
  write_2D_array(nu_cascade.TCS_time(),       identifier_c + "_TCS_time.txt", 1);
  // write_2D_array(nu_cascade.RCS_time(),    identifier_c + "_RCS_time.txt", 1);
  // write_2D_array(nu_cascade.E_time(),   identifier_c + "_E_time.txt", 1);

    // Cascade-specific paramters
  write_2D_array(nu_cascade.Radius(),         identifier_c + "_radial_values.txt", 1);
  write_2D_array(nu_cascade.Density(),        identifier_c + "_density_tx.txt", 1);
  write_2D_array(nu_cascade.PlasmaFreq(),     identifier_c + "_plasma_freq.txt", 1);
  write_2D_array(nu_cascade.Absorption(),     identifier_c + "_absorption.txt", 1);
  write_2D_array(nu_cascade.SkinDepth(),      identifier_c + "_skin_depth.txt", 1);
  write_2D_array(nu_cascade.Reflectance(),    identifier_c + "_reflectance.txt", 1);
  write_2D_array(nu_cascade.Opacity(),        identifier_c + "_opacity.txt", 1);

  // std::string identifier_l = path_out + "Line1D_" + identifier;
  // Line1D thinwire(tx,rx,cs);

  /* Choose what to write out by uncommenting the lines. */
  // write_2D_array(thinwire.Coordinates(),    indentifier_l + "_coords.txt", 1);
  // // write_1D_array(thinwire.Phase(),       indentifier_l + "_phases.txt", 1);
  // write_1D_array(thinwire.Attenuation(),    indentifier_l + "_attenuation.txt", 1);
  // // write_1D_array(thinwire.ArrivalTime(),    indentifier_l + "_arrival_times.txt", 1);
  // write_1D_array(thinwire.Directivty(),     indentifier_l + "_directivity.txt", 1);
  // write_1D_array(thinwire.Polarization(),   indentifier_l + "_polarization.txt", 1);
  //
  // write_1D_array(thinwire.Duration(),       indentifier_l + "_duration.txt", 1);
  // write_1D_array(thinwire.Waveform(),       indentifier_l + "_waveform.txt", 1);
  // write_1D_array(thinwire.Power()),         indentifier_l + "_power.txt", 1);
  // write_1D_array(thinwire.TCS(),            indentifier_l + "_target_cs.txt", 1);
  // write_1D_array(thinwire.RCS(),            indentifier_l + "_target_cs.txt", 1);
  // // write_2D_array(thinwire.Phase_time(),  indentifier_l + "_phase_time.txt", 1);
  // write_2D_array(thinwire.TCS_time(),       indentifier_l + "_TCS_time.txt", 1);
  // // write_2D_array(thinwire.RCS_time(),    indentifier_l + "_RCS_time.txt", 1);
  // // write_2D_array(thinwire.E_time(),   indentifier_l + "_E_time.txt", 1);
}
// End
