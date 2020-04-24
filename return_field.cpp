/* Radar computation - return electric field */
#include "antenna.hh"
#include "cascade.hh"
const string identifier = "test_time_short_4";

using namespace std;

int main(int argc, char** argv){

  // I/O Setup
  string cs_filepath = argv[1];
  // Read data from file, get separate unphysical events
  vector<Cascade> Cascades = load_cascade_file(cs_filepath);
  vector<Cascade>::iterator cs_ptr;

  // Make your detector
  Detector bistatic;
  vector<Antenna> Transmitters = bistatic.get_transmitters();
  vector<Antenna> Receivers = bistatic.get_receivers();
  vector<Antenna>::iterator tx_ptr, rx_ptr;
  double* tx_pos;
  double* rx_pos;
  double* cs_pos;
  double* cs_ang;

  // For every cascade
  for (cs_ptr = Cascades.begin(); cs_ptr < Cascades.end(); cs_ptr++){
    Cascade cs = *cs_ptr;
    cs_pos = cs.get_position();
    cs_ang = cs.get_sph_ang();

    // Compute RCS = 1D OD for now.
    vector<double> od_segmented = cs.get_od_cs_1D();

    // For every transmitter
    for (tx_ptr = Transmitters.begin(); tx_ptr < Transmitters.end(); tx_ptr++){
      Antenna tx = *tx_ptr;
      tx.set_direction(cs);
      tx_pos = tx.get_position();
      double omega = tx.get_w_obs();

      // For every receiver
      for (rx_ptr = Receivers.begin(); rx_ptr < Receivers.end(); rx_ptr++){
        Antenna rx = *rx_ptr;
        rx.set_direction(cs);
        rx_pos = rx.get_position();

        // Set distances, times and E field for segments.

        double l, xpos,ypos, Rt,Rr, birth, arrival, E0 = 100, Er,rad_cs, phase;
        vector<double> pos_row, phase_row, Births, Arrivals, Er_slices;
        vector<vector<double>> position;
        // TBD proper fomula for Transmitter field.

        /* Loop over the 1D segments  */
        for(int i = 0 ; i < nbin; i++){
          l= i*cs.get_L_bin();      // [m] distance from the shower head (starting point)

          // Distances
          xpos = cs_pos[0] + l*cos(cs_ang[1]);
          ypos = cs_pos[1] + l*sin(cs_ang[1]);

          Rt = sqrt(pow(xpos-tx_pos[0],2)+pow(ypos-tx_pos[1],2));
          Rr = sqrt(pow(xpos-rx_pos[0],2)+pow(ypos-rx_pos[1],2));

          pos_row.push_back(l);
          pos_row.push_back(xpos);
          pos_row.push_back(ypos);
          pos_row.push_back(Rt);
          pos_row.push_back(Rr);

          // Times

          //Birth = time where the segment starts scattering.
          // T0 = 0 by definition when the cascade begins (first element = head).
          birth = l/c_vac;
          Births.push_back(birth);

          // (Retarded) time where the scattered radio signal by the segment is produced.
          //production = bith - Rt/c_med;

          // (Advanced) time where the scattered signal by the segments arrives in the receiver.
          arrival= birth + Rr/c_med;
          Arrivals.push_back(arrival);

          // E field
          rad_cs = od_segmented[i];
            // Original
          // Er = E0 / (4*pi*Rt*Rr) * sqrt(rad_cs * rx.area);
            // Dieder Test
          Er = E0 / (Rr) ;
          Er_slices.push_back(Er);

          // Phase
            //, orignal
          //phase = omega*l/cvac + k_obs*(Rr + Rt);
            // Dieder test
          phase = - tx.get_k_obs()*(Rr + Rt);
          phase_row.push_back(phase);
        }

        /* Run time loop */
        vector<double> od_time_row;
        vector<double> Er_time_row;
        vector<double> phase_array_row;
        vector<double> time;
        vector<vector<double>> od_time_profile;
        vector<vector<double>> Er_time_profile;
        vector<vector<double>> phase_array;

        double timestep= 1.0/(100*freq_obs);    // sampling frequency
        double t_start = *min_element(Arrivals.begin(), Arrivals.end()) - 5E-9;
        double t_end   = *max_element(Arrivals.begin(), Arrivals.end()) + 5*tau;
        double dieder_phase, loops = 0;

        // Loop over time.
        for (double t = t_start; t<t_end;t += timestep){
          // Loop over depth segments.
          for (int i = 0; i < nbin; i++){

            // Select some length values
             if(i % 2 == 0){
              // If active, add its contribution.
              if(t>Arrivals[i] && t<(Arrivals[i]+tau)){
                od_time_row.push_back(od_segmented[i]);

                //phase = omega*t + k_obs*Rt;
                dieder_phase = omega*(t) + phase_row[i];

                Er_time_row.push_back(Er_slices[i]* cos(dieder_phase));

                phase_array_row.push_back(omega*(t - Arrivals[i]) + phase);
              } else {
                Er_time_row.push_back(0);
                od_time_row.push_back(0);
                phase_row.push_back(0);
              }

            }

          }
          time.push_back(t);

          Er_time_profile.push_back(Er_time_row);
          od_time_profile.push_back(od_time_row);

          Er_time_row.clear();
          od_time_row.clear();

          phase_array.push_back(phase_array_row);
          phase_row.clear();
        }

        /* Write out all the information. */

        // On screen
        cout << t_start << '\t' << t_end << '\t'<< loops <<endl;
        // This is needed for plotting purposes

        // On file.
        // write_1D_array(Births, "./births.txt");
        // write_1D_array(Arrivals, "./arrivals.txt");
        // write_2D_array(position, identifier + "_positions.txt");
        // write_1D_array(Er_slices, identifier + "_Er_slices.txt");
        // write_2D_array(phase_array, identifier + "_phase.txt");

        write_2D_array(od_time_profile, identifier + "_od_time_profile.txt", 1);
        write_2D_array(Er_time_profile, identifier + "_Er_time_profile.txt", 1);

      } // Closes antenna loop
    } // Closes transmitter loop

  // vector<vector<double>> density = get_density_profile(cs);
  // write_2D_array(density, identifier + "_density_profile.txt");
  write_1D_array(od_segmented, identifier + "_od_segmented_kmax.txt");

  } // Closes cascade loop
  cout << "END" << endl;
}
// End of main
