/* Radar scatter computation - return electric field of a list of events
All variables need to be passed through arguments.
Written to be used in a parallel computation (cluster) enviroment.
*/

#include "macro_scatter.hh"
int main(int argc, char** argv){
  // I/O Setup
  // For a single event, the result is stored in the same folder as the executable.
  // std::string path_out = "";

  // ----- Computational parameters ----------------------------------------------
  std::string identifier  = argv[1];

  // const double freq_obs = 450 * 1E6;
  // const double freq_obs = 5E7;              // [Hz] Default observer frequency

  const double dL = 1; // [cm/bin]
  const double dR = 0.1;
  const double dN = 1; // This is also a radial direction.

  const double sampling = 100;
  // sampling frequency should be between 10x and 100x freq_obs.

  // Cascade parameters

/* VALUES FOR REFERENCE AND FOR TEST */
  // double csenergy   = 1E8;       // 2
  // double csxpos     = 0;         // 3
  // double csypos     = -250 ;     // 4
  // double cszpos     = 0;         // 5
  // double cszenith   = 90;        // 6
  // double csazimuth  = 45;        // 7
  //
  // double txxpos     = 0;         // 8
  // double txypos    = 0;         // 9
  // double txzpos    = 0;         // 10
  // double txxpol    = 0;         // 11
  // double txypol    = 0;         // 12
  // double txzpol    = 1;         // 13
  // double txpower   = 1E3;       // 14
  // double txfreq    = 5E8;       // 15
  //
  // double rxxpos    = 500;       // 16
  // double rxypos    = 0;         // 17
  // double rxzpos    = 0;         // 18
  // double rxxpol    = 0;         // 19
  // double rxypol    = 0;         // 20
  // double rxzpol    = 1;         // 21
  //
  // const double lifetime = 1E-8;               // [s] Plasma lifetime 10 ns
  // const double lifetime = 1E-9;               // [s] Plasma lifetime 1 ms

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

  double lifetime = atof(argv[22]);

  // No need to load a detector, just make two antennas

  Antenna transmitter(txxpos, txypos, txzpos,
                      txxpol, txypol, txzpol,
                      txpower, txfreq); // txgain = 0
  Antenna receiver(   rxxpos, rxypos, rxzpos,
                      rxxpol, rxypos, rxzpos);

  Cascade cascade( 0, csenergy , csxpos, csypos, cszpos,
                  deg2rad(cszenith), deg2rad(csazimuth)) ;

  Antenna& tx = transmitter;
  Antenna& rx = receiver;
  Cascade& cs = cascade;

  std::string identifier_c = "Cascade1D_" + identifier;
  Cascade1D nu_cascade(tx,rx,cs, dL, dR, dN, lifetime, sampling);

  /* Choose what to write out by uncommenting the lines. */
  write_2D_array(nu_cascade.Coordinates(),    identifier_c + "_coords.txt", 1);
  write_1D_array(nu_cascade.Phase(),          identifier_c + "_phase.txt", 1);
  write_1D_array(nu_cascade.Attenuation(),    identifier_c + "_attenuation.txt", 1);
  write_1D_array(nu_cascade.ArrivalTime(),    identifier_c + "_arrival_t.txt", 1);
  write_1D_array(nu_cascade.Directivity(),    identifier_c + "_directivity.txt", 1);
  write_1D_array(nu_cascade.Polarization(),   identifier_c + "_polarization.txt", 1);

  write_1D_array(nu_cascade.Duration(),       identifier_c + "_duration.txt", 1);
  write_1D_array(nu_cascade.Waveform(),       identifier_c + "_waveform.txt", 1);
  write_1D_array(nu_cascade.Power(),          identifier_c + "_power.txt", 1);
  write_1D_array(nu_cascade.TCS(),            identifier_c + "_TCS.txt", 1);
  write_1D_array(nu_cascade.RCS(),            identifier_c + "_RCS.txt", 1);
  write_2D_array(nu_cascade.Phase_time(),     identifier_c + "_phase_time.txt", 1);
  write_2D_array(nu_cascade.TCS_time(),       identifier_c + "_TCS_time.txt", 1);
  write_2D_array(nu_cascade.RCS_time(),       identifier_c + "_RCS_time.txt", 1);
  write_2D_array(nu_cascade.E_time(),         identifier_c + "_E_time.txt", 1);

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
  // // write_1D_array(thinwire.Phase(),       indentifier_l + "_phase.txt", 1);
  // write_1D_array(thinwire.Attenuation(),    indentifier_l + "_attenuation.txt", 1);
  // // write_1D_array(thinwire.ArrivalTime(),    indentifier_l + "_arrival_time.txt", 1);
  // write_1D_array(thinwire.Directivty(),     indentifier_l + "_directivity.txt", 1);
  // write_1D_array(thinwire.Polarization(),   indentifier_l + "_polarization.txt", 1);
  //
  // write_1D_array(thinwire.Duration(),       indentifier_l + "_duration.txt", 1);
  // write_1D_array(thinwire.Waveform(),       indentifier_l + "_waveform.txt", 1);
  // write_1D_array(thinwire.Power()),         indentifier_l + "_power.txt", 1);
  // write_1D_array(thinwire.TCS(),            indentifier_l + "_TCS.txt", 1);
  // write_1D_array(thinwire.RCS(),            indentifier_l + "_RCS.txt", 1);
  // // write_2D_array(thinwire.Phase_time(),  indentifier_l + "_phase_time.txt", 1);
  // write_2D_array(thinwire.TCS_time(),       indentifier_l + "_TCS_time.txt", 1);
  // // write_2D_array(thinwire.RCS_time(),    indentifier_l + "_RCS_time.txt", 1);
  // // write_2D_array(thinwire.E_time(),   indentifier_l + "_E_time.txt", 1);
}
// End
