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

  std::vector<std::vector<double>> Coordinates();
  std::vector<double> Attenuation();
  std::vector<double> arrivals();
  std::vector<double> phase();

  std::vector<double> duration();
  std::vector<double> waveform();
  std::vector<double> RCS();
  std::vector<std::vector<double>> phase_time();
  std::vector<std::vector<double>> rcs_time();
  std::vector<std::vector<double>> wave_time();

  // Attenuation model goes here

protected:

  Antenna fTX;
  Antenna fRX;
  Cascade fCS;
  std::vector<double> fRCS; // currently given in cm^2

  // Always done at construction, no reason to be changed.
  // void SetAtDirection(Antenna &at);
  // void SetAtDirCenter(Antenna &at);

  /* Segment the cascade in N segments and compute the properties for each segment:
    Positions, distances, times and E fields */
  void SetSegments(const double& nSeg);
  void run_time_loop();


private:

  std::vector<double> fPhase;
  std::vector<double> fArrivalTime;
  std::vector<double> fAttenuation;
  std::vector<double> fPolarization;
  std::vector<std::vector<double>> fSegmentCoord;
  // Length, xpos, ypos, zpos, R_TX, R_RX

  std::vector<double> fDuration; // [ns]
  std::vector<double> fWaveform; // [V/m]


  std::vector<std::vector<double>> fRCSTime;
  std::vector<std::vector<double>> fPhaseTime;
  std::vector<std::vector<double>> fWaveformTime;


};

class Line1D: public Scatter{
public:

  Line1D(Antenna& tx, Antenna& rx, Cascade& cs);
};

class Cascade1D: public Scatter {
public:

  Cascade1D(Antenna& tx, Antenna& rx, Cascade& cs);


  // Accesors
  double  Delta();
  double  Dot();

  std::vector<std::vector<double>> Density();
  double Rwaist();
  std::vector<double> Rcrit();
  std::vector<std::vector<double>> PlasmaFreq();
  std::vector<std::vector<double>> Absorption();
  std::vector<std::vector<double>> SkinDepth();
  std::vector<std::vector<double>> Reflectance();
  std::vector<std::vector<double>> Opacity();
  std::vector<double> RCS();

  // std::vector<double> radar_cs();


private:

  double fDot = 0;             // Dot (inner) product with cascade direction.
  // The cosine of the angle between them.
  double fDelta = 0;          // The angle.
  double cD, sD;              // cosine and sine of the angle

  int nPerp, nPar;
  double fPerp, fPar;
  double dPerp = 1.0;
  double dPar  = 1.0;
  double dNorm = 1.0;

  std::vector<std::vector<double>> fDensity;
  double fRwaist;
  std::vector<double> fRcrit; // [cm]
  std::vector<std::vector<double>> fPlasmaFrequency;
  std::vector<std::vector<double>> fAbsorption;
  std::vector<std::vector<double>> fSkinDepth;
  std::vector<std::vector<double>> fReflectance;
  std::vector<std::vector<double>> fOpacity;

  void SetDensity();
  void Rcrit(const std::vector<std::vector<double>> &density, const double & freq);
  void PlasmaFreq(const std::vector<std::vector<double>> &density);
  void Absorption(const std::vector<std::vector<double>> &density, const double & freq);
  void SkinDepth(const std::vector<std::vector<double>> &density, const double & freq);
  void Reflectance(const std::vector<std::vector<double>> &absorption);
  void Opacity(const std::vector<std::vector<double>> &absorption);
  void RCS(const std::vector<std::vector<double>> &opacity);          // Const layers
  // void RCS(const std::vector<std::vector<double>> &reflectance);   // Non const layers

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
