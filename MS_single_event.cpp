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

  int evt = atof(argv[1]);
  double cenergy = atof(argv[2]);
  double xpos = atof(argv[3]);
  double ypos = atof(argv[4]);
  double zpos = atof(argv[5]);
  double czenith = atof(argv[6]);
  double cazimuth = atof(argv[7]);
  std::string identifier  = argv[8];

  //std::string identifier = argv[1];
  //int evt = 0;
  //double cenergy = 1E9;
  //double xpos = 0;
  //double ypos = 0;
  //double zpos = 250;
  //double czenith = 90;
  //double cazimuth = 0;

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

  // For every transmitter
  for (auto& tx : bistatic.Transmitters()){
    // For every receiver
    for (auto& rx : bistatic.Receivers()){
      // Make the bistatic event.
      // std::cout << "Event: " << std::to_string((int)cs.event()) << endl;

      Line1D event0(tx,rx,cs);

      /* Write out all the information. */
      // write_2D_array(event0.segement_coords(),  identifier_l + "_coords.txt", 1);
      write_1D_array(event0.Attenuation(),        identifier_l + "_attenuation.txt", 1);
      // write_1D_array(event0.arrivals(),         identifier_l + "_t_arrivals.txt", 1);
      // write_1D_array(event0.phase(),            identifier_l + "_phases.txt", 1);
      write_1D_array(event0.duration(),         identifier_l + "_duration.txt", 1);
      write_1D_array(event0.Waveform(),         identifier_l + "_waveform.txt", 1);
      // write_2D_array(event0.phase_time(),       identifier_l + "_phase_time_profile.txt", 1);
      // write_2D_array(event0.wave_time(),        identifier_l + "_Er_time_profile.txt", 1);

      // std::string identifier_c = path_out + "Cascade1D_" + std::to_string((int) rad2deg( cs.sph_angles()[1]) ) + "_deg";
      //
      Cascade1D event(tx,rx,cs);
      // std::cout << cazimuth << std::endl;


      write_2D_array(event.Density(),     identifier_c + "_density_tx.txt", 1);
      write_2D_array(event.PlasmaFreq(),  identifier_c + "_plasma_freq.txt", 1);
      write_2D_array(event.Absorption(),  identifier_c + "_absorption.txt", 1);
      write_2D_array(event.SkinDepth(),   identifier_c + "_skin_depth.txt", 1);
      write_2D_array(event.Reflectance(), identifier_c + "_reflectance_matrix.txt", 1);
      write_2D_array(event.Opacity(),     identifier_c + "_opacity_matrix.txt", 1);
      write_1D_array(event.RCS(),         identifier_c + "_radar_cs.txt", 1);
      write_2D_array(event.Coordinates(),   identifier_c + "_coords.txt", 1);
      write_1D_array(event.Attenuation(),         identifier_c + "_attenuation.txt", 1);
      // write_1D_array(event.arrivals(),          identifier_c + "_t_arrivals.txt", 1);
      // write_1D_array(event.phase(),             identifier_c + "_phases.txt", 1);
      write_1D_array(event.duration(),          identifier_c + "_duration.txt", 1);
      write_1D_array(event.Waveform(),          identifier_c + "_waveform.txt", 1);
      // write_2D_array(event.wave_time(),         identifier_c + "_Er_time_profile.txt", 1);
      // write_2D_array(event.rcs_time(),          identifier_c + "_rcs_time_profile.txt", 1);
      // write_2D_array(event.phase_time(),        identifier_c + "_phase_time_profile.txt", 1);

    }
  }
}
// End
