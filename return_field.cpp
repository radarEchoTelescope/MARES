/* Radar computation - return electric field */
#include "macro_scatter.hh"

int main(int argc, char** argv){
  const std::string identifier = "test_time_short_4";

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

        // tx.set_direction(cs);
        // rx.set_direction(cs);
        // tx_pos = tx.position();
        // rx_pos = rx.position();
        // Compute RCS = 1D OD for now.
        // vector<double> od_segmented = get_od_cs_1D();



        /* Write out all the information. */

        // On screen
        // cout << t_start << '\t' << t_end << '\t'<< loops <<endl;
        // This is needed for plotting purposes

        // On file.
        // write_2D_array(position, identifier + "_positions.txt");
        // write_1D_array(Er_slices, identifier + "_Er_slices.txt");
        // write_2D_array(phase_array, identifier + "_phase.txt");

        // write_2D_array(od_cs_1D_time, identifier + "_od_time_profile.txt", 1);
        // write_2D_array(Er_time, identifier + "_Er_time_profile.txt", 1);

        // vector<vector<double>> density = get_density_profile(cs);
        // write_2D_array(density, identifier + "_density_profile.txt");
        // write_1D_array(od_segmented, identifier + "_od_segmented_kmax.txt");
      }
    }
  }
  cout << "END" << endl;
}
// End of main
