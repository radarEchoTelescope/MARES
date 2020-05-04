/* Radar computation - return electric field */
#include "macro_scatter.hh"
const std::string identifier = "Particle1to500";

int main(int argc, char** argv){

  // I/O Setup
  string cs_filepath = argv[1];
  // Read data from file, get separate unphysical events
  vector<Cascade> cascade_list = load_cascade_file(cs_filepath);

  // Make your detector
  string det_type = "bistatic";
  Detector bistatic(det_type);

  // For every cascade
  for (auto& cs: cascade_list){
    // For every transmitter
    for (auto& tx : bistatic.get_transmitters()){
      // For every receiver
      for (auto& rx : bistatic.get_receivers()){
        // Make the bistatic event.
        Scatter1D event(tx,rx,cs);

        /* Write out all the information. */

        // On screen
        // cout << event.transmitter().power() << endl;
        // cout << t_start << '\t' << t_end << '\t'<< loops <<endl;
        // This is needed for plotting purposes

        // On file.

        write_1D_array(event.get_length(), identifier + "_length.txt", 1);
        write_1D_array(event.get_xpos(), identifier + "_xpos.txt", 1);
        write_1D_array(event.get_ypos(), identifier + "_ypos.txt", 1);
        write_1D_array(event.get_Rt(), identifier + "_Rt.txt", 1);
        write_1D_array(event.get_Rr(), identifier + "_Rr.txt", 1);
        write_1D_array(event.get_time(), identifier + "_time.txt", 1);
        write_1D_array(event.get_amplitude(), identifier + "_Er_slices.txt", 1);
        // write_1D_array(od_segmented, identifier + "_od_segmented_kmax.txt");

        std::vector<std::vector<double>> Er_time = event.get_Er_time();
        write_2D_array(Er_time, identifier + "_Er_time_profile.txt", 1);
        // write_2D_array(od_cs_1D_time, identifier + "_od_time_profile.txt", 1);
        write_2D_array(event.get_phase_time(), identifier + "_phase.txt", 1);

        // vector<vector<double>> density = get_density_profile(cs);
        // write_2D_array(density, identifier + "_density_profile.txt");

      }
    }
  }
  cout << "END" << endl;
}
// End of main
