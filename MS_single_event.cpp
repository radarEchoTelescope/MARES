/* Radar computation - return electric field of a list of events*/
#include "macro_scatter.hh"
int main(int argc, char** argv){

  // I/O Setup

  // For a single event, we don't know a reference path, the result is stored
  // in the same folder as the executable.
  std::string path_out = "";
  // std::string cs_filepath = argv[1];
  // std::string path_out = cs_filepath.substr(0,cs_filepath.find_last_of("."));

  std::string identifier  = argv[2];
  // const std::string identifier = "Particle1to500";
  // std::string identifier = "../Visualization/testing/testing_500_particles";



  // Cascade params
  int evt = 0;
  double cenergy = 1E9;
  double xpos = 250;
  double ypos = -250;
  double zpos = 0;
  double czenith = pi/2;
  double cazimuth = atof(argv[1]);
  double nenergy = 1E9;
  double nzenith = 0;
  double nazimuth = 0;
  double oneweight = 1;

  Cascade cascade(evt, cenergy, xpos, ypos, zpos, czenith, cazimuth,
                  nenergy, nzenith, nazimuth, oneweight);

  Cascade& cs = cascade;

  // Run here the validity check for the cascasde?

  // Make your detector
  string det_type = "bistatic";
  Detector bistatic(det_type);

  // For every transmitter
  for (auto& tx : bistatic.transmitters()){
    // For every receiver
    for (auto& rx : bistatic.receivers()){
      // Make the bistatic event.
      // std::cout << "Event: " << std::to_string((int)cs.event()) << endl;

      std::string identifier_l = path_out + "Line1D_" + identifier;
      Line1D event0(tx,rx,cs);

      // std::cout << cazimuth << std::endl;

      /* Write out all the information. */
      write_2D_array(event0.segement_coords(),  identifier_l + "_coords.txt", 1);
      write_1D_array(event0.amplitude(),        identifier_l + "_Er_slices.txt", 1);
      write_1D_array(event0.arrivals(),         identifier_l + "_t_arrivals.txt", 1);
      write_1D_array(event0.phase(),            identifier_l + "_phases.txt", 1);
      write_1D_array(event0.duration(),         identifier_l + "_duration.txt", 1);
      write_1D_array(event0.waveform(),         identifier_l + "_waveform.txt", 1);
      write_2D_array(event0.phase_time(),       identifier_l + "_phase_time_profile.txt", 1);
      write_2D_array(event0.wave_time(),        identifier_l + "_Er_time_profile.txt", 1);


      // std::string identifier_c = path_out + "Cascade1D_" + std::to_string((int) rad2deg( cs.sph_angles()[1]) ) + "_deg";
      std::string identifier_c = path_out + "Cascade1D_" + identifier;
      Cascade1D event(tx,rx,cs);
      // std::cout << cazimuth << std::endl;

      write_2D_array(event.get_density_cs(),    identifier_c + "_density_cs.txt", 1);
      write_2D_array(event.get_density_tx(),    identifier_c + "_density_tx.txt", 1);
      write_2D_array(event.get_plasma_freq(),   identifier_c + "_plasma_freq.txt", 1);
      write_2D_array(event.get_skin_depth(),    identifier_c + "_skin_depth.txt", 1);
      write_2D_array(event.get_reflectance(),   identifier_c + "_reflectance_matrix.txt", 1);
      write_2D_array(event.get_reflectivity(),  identifier_c + "_reflectivity_matrix.txt", 1);
      write_1D_array(event.radar_cs(),          identifier_c + "_radar_cs.txt", 1);
      write_2D_array(event.segement_coords(),   identifier_c + "_coords.txt", 1);
      write_1D_array(event.amplitude(),         identifier_c + "_Er_slices.txt", 1);
      write_1D_array(event.arrivals(),          identifier_c + "_t_arrivals.txt", 1);
      write_1D_array(event.phase(),             identifier_c + "_phases.txt", 1);
      write_1D_array(event.duration(),          identifier_c + "_duration.txt", 1);
      write_1D_array(event.waveform(),          identifier_c + "_waveform.txt", 1);
      write_2D_array(event.wave_time(),         identifier_c + "_Er_time_profile.txt", 1);
      write_2D_array(event.rcs_time(),          identifier_c + "_rcs_time_profile.txt", 1);
      write_2D_array(event.phase_time(),        identifier_c + "_phase_time_profile.txt", 1);

    }
  }
}
// End
