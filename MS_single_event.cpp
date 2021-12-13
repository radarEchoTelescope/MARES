/* Radar computation - return electric field of a list of events*/
#include "macro_scatter.hh"
int main(int argc, char** argv){

  // I/O Setup

  // For a single event, we don't know a reference path, the result is stored
  // in the same folder as the executable.
  std::string path_out = "";
  // std::string path_out = cs_filepath.substr(0,cs_filepath.find_last_of("."));

  // const std::string identifier = "Particle1to500";
  // std::string identifier = "../Visualization/testing/testing_500_particles";

  // Cascade params

  // int evt = atof(argv[1]);
  // double cenergy = atof(argv[2]);
  // double xpos = atof(argv[3]);
  // double ypos = atof(argv[4]);
  // double zpos = atof(argv[5]);
  // double czenith = atof(argv[6]);
  // double cazimuth = atof(argv[7]);
  // std::string identifier  = argv[8];

  int evt = 0;
  double cenergy = 1E10;
  double xpos = 0;
  double ypos = -250 ;
  double zpos = 0;
  double czenith = 90;
  double cazimuth = 0;
  std::string identifier  = argv[1];

  double nenergy = 1E9;
  double nzenith = 0;
  double nazimuth = 0;
  double oneweight = 1;

  czenith = deg2rad(czenith);
  cazimuth = deg2rad(cazimuth);

  std::string identifier_l = path_out + "Line1D_" + identifier;
  std::string identifier_c = path_out + "Cascade1D_" + identifier;


  Cascade cascade(evt, cenergy, xpos, ypos, zpos, czenith, cazimuth,
                  nenergy, nzenith, nazimuth, oneweight);

  Cascade& cs = cascade;
  // write_2D_array(cs.Density(),    identifier_c + "_density_cs.txt", 1);

  // Run here the validity check for the cascasde?

  // Make your detector
  string det_type = "bistatic";
  Detector bistatic(det_type);


  // For every transmitter-receiver pair
  for (auto& tx : bistatic.Transmitters()){
    for (auto& rx : bistatic.Receivers()){

      // Chose the type of scatter event: Cascade, Line, Cylinder?
      Line1D thinwire(tx,rx,cs);

      Cascade1D nu_cascade(tx,rx,cs);

      // Is a specific name needed per every event?
      // std::string identifier_c = path_out + "Cascade1D_" + std::to_string((int) rad2deg( cs.sph_angles()[1]) ) + "_deg";

      /* Choose what to write out by uncommenting the lines. */
      write_2D_array(thinwire.Coordinates(),      identifier_l + "_coords.txt", 1);
      // write_1D_array(thinwire.Attenuation(),   identifier_l + "_attenuation.txt", 1);
      // // write_1D_array(thinwire.arrivals(),   identifier_l + "_t_arrivals.txt", 1);
      // // write_1D_array(thinwire.phase(),      identifier_l + "_phases.txt", 1);
      write_1D_array(thinwire.Duration(),         identifier_l + "_duration.txt", 1);
      write_1D_array(thinwire.Waveform(),         identifier_l + "_waveform.txt", 1);
      // write_2D_array(thinwire.phase_time(),    identifier_l + "_phase_time_profile.txt", 1);
      // write_2D_array(thinwire.wave_time(),     identifier_l + "_Er_time_profile.txt", 1);


      write_2D_array(nu_cascade.Radius(),         identifier_c + "_radial_values.txt", 1);
      write_2D_array(nu_cascade.Density(),        identifier_c + "_density_tx.txt", 1);
      write_2D_array(nu_cascade.PlasmaFreq(),     identifier_c + "_plasma_freq.txt", 1);
      write_2D_array(nu_cascade.Absorption(),     identifier_c + "_absorption.txt", 1);
      write_2D_array(nu_cascade.SkinDepth(),      identifier_c + "_skin_depth.txt", 1);
      write_2D_array(nu_cascade.Reflectance(),    identifier_c + "_reflectance.txt", 1);
      write_2D_array(nu_cascade.Opacity(),        identifier_c + "_opacity.txt", 1);
      write_1D_array(nu_cascade.TCS(),            identifier_c + "_target_cs.txt", 1);
      write_2D_array(nu_cascade.Coordinates(),    identifier_c + "_coords.txt", 1);
      write_1D_array(nu_cascade.Polarization(),   identifier_c + "_polarization.txt", 1);
      write_1D_array(nu_cascade.Attenuation(),    identifier_c + "_attenuation.txt", 1);
      // write_1D_array(nu_cascade.arrivals(),    identifier_c + "_t_arrivals.txt", 1);
      // write_1D_array(nu_cascade.phase(),       identifier_c + "_phases.txt", 1);
      write_1D_array(nu_cascade.Duration(),       identifier_c + "_duration.txt", 1);
      write_1D_array(nu_cascade.Waveform(),       identifier_c + "_waveform.txt", 1);
      // write_2D_array(nu_cascade.wave_time(),   identifier_c + "_Er_time_profile.txt", 1);
      // write_2D_array(nu_cascade.rcs_time(),    identifier_c + "_rcs_time_profile.txt", 1);
      // write_2D_array(nu_cascade.phase_time(),  identifier_c + "_phase_time_profile.txt", 1);

    }
  }
}
// End
