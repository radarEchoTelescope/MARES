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

protected:

  Antenna tx;
  Antenna rx;
  Cascade cs;

  void set_direction(Antenna& at);
  void set_IRT_direction(Antenna& at);

};


class Scatter1D: public Scatter {
public:

  Scatter1D(Antenna& tx, Antenna& rx, Cascade& cs);

private:


    // Make these things arrays?
  std::vector<double> Length, Xpos, Ypos, Rt, Rr, Phase, Arrivals, Amplitude;
  void set_1D_values();

  std::vector<double> od_cs_1D;
  void get_od_cs_1D();

  vector<vector<double>> od_cs_1D_time;
  vector<vector<double>> Er_time;
  vector<vector<double>> phase_time;

  void run_time_loop();
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
