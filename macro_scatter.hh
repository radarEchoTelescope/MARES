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


  std::vector<double> get_length(){ return Length ; }
  std::vector<double> get_xpos(){ return Xpos ; }
  std::vector<double> get_ypos(){ return Ypos ; }
  std::vector<double> get_Rt(){ return Rt ; }
  std::vector<double> get_Rr(){ return Rr ; }
  std::vector<double> get_arrivals(){ return Arrivals ; }


  std::vector<double> get_amplitude();





  std::vector<std::vector<double>> get_od_cs_time();
  std::vector<std::vector<double>> get_Er_time();
  std::vector<std::vector<double>> get_phase_time(){return phase_time;}
  std::vector<double> test_array;
  std::vector<double> get_time(){ return test_array ; }


private:

  std::vector<double> Length, Xpos, Ypos, Rt, Rr, Phase, Arrivals, Amplitude;
  void set_1D_values();

  std::vector<double> od_cs;
  void get_od_cs();


  vector<vector<double>> od_cs_time;
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


//------------------------------------------------------------------------------
/* NO LONGER IN USE */

/*
// Parametrized attenuation length for the Ross Ice Shelf, South Pole.
double att(double freq_obs){
double a1=469;                  // [m] Attanuation length parameter
double a2=-0.205;               // Attanuation length parameter
double a3=4.87E-5;              // Attanuation length parameter
return a1+a2*freq_obs/1E6+a3*pow(freq_obs/1E6,2);
}
*/
