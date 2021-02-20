// Calculates the return power for a bi-static radat setup
// Branched of power.C at 23/1/2020
// Enrique Huesca Santiago, 10-2019
// Original FORTRAN code by Krijn D. de Vries 20-10-2019

#ifndef MACRO_SCATTER
#define MACRO_SCATTER

// #include "macro_settings.hh"
#include "cascade.hh"
#include "antenna.hh"
#include "IceRayTracing.hh"
// #define NDEBUG     // Turn off debug.

/* Simple scatter event with bistatic configuration*/
class Scatter {
public:

  Scatter(Antenna& tx, Antenna& rx, Cascade& cs);

  Antenna transmitter();
  Antenna receiver();
  Cascade cascade();

  std::vector<std::vector<double>> segement_coords();
  std::vector<double> amplitude();
  std::vector<double> arrivals();
  std::vector<double> phase();

  std::vector<double> duration();
  std::vector<double> waveform();
  std::vector<double> radar_cs();
  std::vector<std::vector<double>> phase_time();
  std::vector<std::vector<double>> rcs_time();
  std::vector<std::vector<double>> wave_time();


protected:

  Antenna tx;
  Antenna rx;
  Cascade cs;

  void set_direction(Antenna& at);
  void set_direction_center(Antenna& at);

  void set_segments();
  void run_time_loop();

  std::vector<double> _rcs;

private:

  std::vector<double> _phase, _arrival, _amplitude;
  std::vector<double> _duration; // [ns]
  std::vector<double> _waveform; // [V/m]

  std::vector<std::vector<double>> _segment_coords; // Length, xpos, ypos, zpos, Rt/D1i, Rr/D2i
  std::vector<std::vector<double>> _rcs_time;
  std::vector<std::vector<double>> _phase_time;
  std::vector<std::vector<double>> _er_time;


};

class Line1D:public Scatter{
public:

  Line1D(Antenna& tx, Antenna& rx, Cascade& cs);
};

class Cascade1D: public Scatter {
public:

  Cascade1D(Antenna& tx, Antenna& rx, Cascade& cs);

  std::vector<std::vector<double>> get_density_cs();
  std::vector<std::vector<double>> get_density_tx();

  std::vector<std::vector<double>> get_plasma_freq();
  std::vector<std::vector<double>> get_absorption();
  std::vector<std::vector<double>> get_skin_depth();

  std::vector<std::vector<double>> get_reflectance();
  std::vector<std::vector<double>> get_reflectivity();

  std::vector<double> radar_cs();


private:
  double alpha;

  std::vector<std::vector<double>> coords; // To get rid of?
  std::vector<std::vector<double>> density_cs;
  std::vector<std::vector<double>> density_tx;
  void set_rotated_density();

  double absorption(double& dens);
  double skin_depth(double& dens);

  std::vector<std::vector<double>> fplasma_matrix;
  std::vector<std::vector<double>> absorption_matrix;
  std::vector<std::vector<double>> skin_depth_matrix;
  void set_fplasma();
  void set_absorption();
  void set_skin_depth();

  std::vector<std::vector<double>> reflectance_matrix;
  std::vector<std::vector<double>> reflectivity_matrix;
  void set_reflectance();
  void set_reflectivity();

  void set_radar_cs();
};



// std::vector<Scatter> run_scatter_events(Detector det, std::vector<Cascade> cascade_list){
//   std::vector<Scatter> event_list;
//   // For every cascade
//   for (auto& cs: cascade_list){
//     // For every transmitter
//     for (auto& tx : det.get_transmitters()){
//       // For every receiver
//       for (auto& rx : det.get_receivers()){
//         // Make the bistatic event.
//         event_list.push_back(Scatter(tx,rx,cs));
//       }
//     }
//   }
// }


#endif
