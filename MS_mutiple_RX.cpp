/* Radar computation - return electric field of a list of events*/
#include "macro_scatter.hh"
int main(int argc, char** argv){

// For now, we stick with the variables being defined either as write-in or
// fed through the command line.

  // I/O Setup

  // As default, the result is stored in the same folder as the executable.
  std::string path_out = "";
  // std::string path_out = cs_filepath.substr(0,cs_filepath.find_last_of("."));

  std::string identifier  = argv[1];

  /* Make the cascade OR parse through load_cascade_file() */
  int evt = 0;
  double cenergy = 1E9; //[GeV]
  // double xpos = 0;
  double xpos = 0;
  double ypos = 0 ;
  double zpos = 0;
  double czenith = 90;
  double cazimuth = 0;
  czenith = deg2rad(czenith);
  cazimuth = deg2rad(cazimuth);

  double nenergy = 1E9;
  double nzenith = 0;
  double nazimuth = 0;
  double oneweight = 1;


  Cascade cascade(evt, cenergy, xpos, ypos, zpos, czenith, cazimuth,
                  nenergy, nzenith, nazimuth, oneweight);

  Cascade& cs = cascade;

  /* arange() works as in python: Start,stop, steps
    For a single value, start = stop and step > 0.
    E.g: 90 = arange<double>(90, 90, 1)
  */

  double R = 200; // Antenna - Wire distance [m]
  std::vector<double> TX_angles = arange<double>(90, 90, 1);
  std::vector<double> RX_angles = arange<double>(0, 360, 6);

  // Make your detector
  Detector lab;
  for(auto& i : TX_angles){
    double ii = deg2rad( i );
    lab.add_antenna( Antenna(1E4,  R*cos(ii), R*sin(ii),0 ) );
  }

  for(auto& j : RX_angles){
    double jj = deg2rad(j);
    lab.add_antenna( Antenna(0,  R*cos(jj) , R*sin(jj), 0 ) );
  }

  int i = 0;
  int j;

  // For every transmitter-receiver pair
  for (auto& tx : lab.Transmitters()){
    j = 0;
    for (auto& rx : lab.Receivers()){
      std::string identifier_c = path_out + identifier +
                        + "_TX_" + std::to_string((int) TX_angles[i]) + "_deg"
                        + "_RX_" + std::to_string((int) RX_angles[j]) + "_deg";

      Cascade1D nu_cascade(tx,rx,cs);

      // write_2D_array(nu_cascade.Coordinates(),    identifier_c + "_coords.txt", 1);
      // write_2D_array(nu_cascade.Radius(),         identifier_c + "_radial_values.txt", 1);
      // write_2D_array(nu_cascade.Density(),        identifier_c + "_density_tx.txt", 1);
      // write_2D_array(nu_cascade.PlasmaFreq(),     identifier_c + "_plasma_freq.txt", 1);
      // write_2D_array(nu_cascade.Absorption(),     identifier_c + "_absorption.txt", 1);
      // write_2D_array(nu_cascade.SkinDepth(),      identifier_c + "_skin_depth.txt", 1);
      // write_2D_array(nu_cascade.Reflectance(),    identifier_c + "_reflectance.txt", 1);
      // write_2D_array(nu_cascade.Opacity(),        identifier_c + "_opacity.txt", 1);
      write_1D_array(nu_cascade.TCS(),            identifier_c + "_TCS.txt", 1);
      write_1D_array(nu_cascade.RCS(),            identifier_c + "_RCS.txt", 1);
      write_2D_array(nu_cascade.Coordinates(),    identifier_c + "_coords.txt", 1);
      write_1D_array(nu_cascade.Polarization(),   identifier_c + "_polarization.txt", 1);
      write_1D_array(nu_cascade.Attenuation(),    identifier_c + "_attenuation.txt", 1);
      write_1D_array(nu_cascade.Directivity(),      identifier_c + "_directivity.txt", 1);

      // // write_1D_array(nu_cascade.arrivals(),    identifier_c + "_t_arrivals.txt", 1);
      // // write_1D_array(nu_cascade.Phase(),       identifier_c + "_phase.txt", 1);
      write_1D_array(nu_cascade.Duration(),       identifier_c + "_duration.txt", 1);
      write_1D_array(nu_cascade.Waveform(),       identifier_c + "_waveform.txt", 1);
      write_2D_array(nu_cascade.TCS_time(),   identifier_c + "_tcs_time_profile.txt", 1);
      write_2D_array(nu_cascade.RCS_time(),   identifier_c + "_rcs_time_profile.txt", 1);
      write_2D_array(nu_cascade.E_time(),   identifier_c + "_Er_time_profile.txt", 1);

      // /* Choose what to write out by uncommenting the lines. */
      // std::string identifier_l = path_out + identifier +
      //                             + "_TX_" + std::to_string((int) TX_angles[i]) + "_deg"
      //                             + "_RX_" + std::to_string((int) RX_angles[j]) + "_deg";
      // Line1D thinwire(tx,rx,cs);
      //
      //
      // // write_2D_array(thinwire.Coordinates(),      identifier_l + "_coords.txt", 1);
      // // // write_1D_array(thinwire.Attenuation(),   identifier_l + "_attenuation.txt", 1);
      // // // // write_1D_array(thinwire.arrivals(),   identifier_l + "_t_arrivals.txt", 1);
      // write_1D_array(thinwire.TCS(),              identifier_l + "_TCS.txt", 1);
      // write_1D_array(thinwire.Duration(),         identifier_l + "_duration.txt", 1);
      // write_1D_array(thinwire.Directivity(),      identifier_l + "_directivity.txt", 1);
      // write_1D_array(thinwire.Waveform(),         identifier_l + "_waveform.txt", 1);
      // write_1D_array(thinwire.Power(),            identifier_l + "_power.txt", 1);
      // write_1D_array(thinwire.RCS(),              identifier_l + "_RCS.txt", 1);
      // // write_2D_array(thinwire.Phase_time(),    identifier_l + "_phase_time_profile.txt", 1);
      // // write_2D_array(thinwire.E_time(),     identifier_l + "_Er_time_profile.txt", 1);

      j++;
    }
  i++;
  }
}
// End
