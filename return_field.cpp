/* Radar computation - return electric field */
#include "macro_scatter.hh"
// const std::string identifier = "Particle1to500";

int main(int argc, char** argv){

  // I/O Setup
  string cs_filepath = argv[1];
  // Read data from file, get separate unphysical events
  vector<Cascade> cascade_list = load_cascade_file(cs_filepath);

  // Make your detector
  string det_type = "bistatic";
  Detector bistatic(det_type);

  std::string path_out = cs_filepath.substr(0,cs_filepath.find_last_of("."));

  // std::string identifier = "../Visualization/testing/testing_500_particles";
  // std::string identifier = "case_2b_new_skin";

  // For every cascade
  for (auto& cs: cascade_list){
    // For every transmitter
    for (auto& tx : bistatic.get_transmitters()){
      // For every receiver
      for (auto& rx : bistatic.get_receivers()){
        // Make the bistatic event.

        Line1D event0(tx,rx,cs);

        std::string identifier0 = path_out + "_Line1D_" + std::to_string((int)cs.event())  + "_pi_over_24";
        write_2D_array(event0.segement_coords(), identifier0 + "_coords.txt", 1);
        write_1D_array(event0.amplitude(), identifier0 + "_Er_slices.txt", 1);
        write_1D_array(event0.arrivals(), identifier0 + "_t_arrivals.txt", 1);
        write_1D_array(event0.phase(), identifier0 + "_phases.txt", 1);
        write_1D_array(event0.time(), identifier0 + "_time.txt", 1);
        write_1D_array(event0.waveform(), identifier0 + "_waveform.txt", 1);
        write_2D_array(event0.phase_time(), identifier0 + "_phase_time_profile.txt", 1);
        write_2D_array(event0.wave_time(), identifier0 + "_Er_time_profile.txt", 1);

        Cascade1D event(tx,rx,cs);

        std::string identifier = path_out + "_Cascade1D_" + std::to_string((int)cs.event()) + "_pi_over_24";
        write_2D_array(event.get_density_cs(), identifier + "_density_cs.txt", 1);
        write_2D_array(event.get_density_tx(), identifier + "_density_tx.txt", 1);
        write_2D_array(event.get_reflectance(), identifier + "_reflectance_matrix.txt", 1);
        write_2D_array(event.get_reflectivity(), identifier + "_reflectivity_matrix.txt", 1);
        write_1D_array(event.radar_cs(), identifier + "_radar_cs.txt", 1);
        write_2D_array(event.segement_coords(), identifier + "_coords.txt", 1);
        write_1D_array(event.amplitude(), identifier + "_Er_slices.txt", 1);
        write_1D_array(event.arrivals(), identifier + "_t_arrivals.txt", 1);
        write_1D_array(event.phase(), identifier + "_phases.txt", 1);
        write_1D_array(event.time(), identifier + "_time.txt", 1);
        write_1D_array(event.waveform(), identifier + "_waveform.txt", 1);
        write_2D_array(event.wave_time(), identifier + "_Er_time_profile.txt", 1);
        write_2D_array(event.rcs_time(), identifier + "_rcs_time_profile.txt", 1);
        write_2D_array(event.phase_time(), identifier + "_phase_time_profile.txt", 1);

        /* Write out all the information. */

        // On screen
        // std::cout << "Check" << std::endl;
        // std::cout << densit[10][10] << std::endl;
        // cout << t_start << '\t' << t_end << '\t'<< loops <<endl;
        // std::vector<double> rcs = event.get_reflectivity().back();
        // for(int i=0; i < rcs.size(); i++) {std::cout << rcs.at(i) << ' ';}
        // std::cout << event.get_reflectivity().back()[250] << '\t' << event.get_reflectivity()[nbin - 1][250] << std::endl;
      }
    }
  }
  // cout << "END" << endl;
}
// End of main
