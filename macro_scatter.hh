// Calculates the return power for a bi-static radar setup
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

  Antenna TX();
  Antenna RX();
  Cascade CS();

  // Inner product between CS direction and TX-CS (line of sight).
  double  Dot();
  // Angle of the inner product (between 0 and pi).
  double  Delta();

  std::vector<std::vector<double>> Coordinates();
  std::vector<double> Phase();
  std::vector<double> Attenuation();
  std::vector<double> ArrivalTime();
  std::vector<double> Directivity();
  std::vector<double> Polarization();

  std::vector<double> Duration();
  std::vector<double> Waveform();
  std::vector<double> Power();
  std::vector<double> TCS();
  std::vector<double> RCS();
  std::vector<std::vector<double>> Phase_time();
  std::vector<std::vector<double>> TCS_time();
  std::vector<std::vector<double>> RCS_time();
  std::vector<std::vector<double>> E_time();

  // Attenuation model goes here

protected:

  double L, R;
  double dL, dR, dN;
  int nL, nR;

  Antenna fTX;
  Antenna fRX;
  Cascade fCS;
  std::vector<double> fTCS; // currently given in cm^2

  double cD, sD;              // cosine and sine of the angle


  // Always done at construction, no reason to be changed.
  // void SetAtDirection(Antenna &at);
  // void SetAtDirCenter(Antenna &at);

  /* Segment the cascade in N segments and compute the properties for each segment:
    Positions, distances, times and E fields */
  void SetSegments(const int& nSeg);

  /* Computes the interference of the segements over the time values */
  void RunScatter();


private:

  std::vector<double> fPhase;
  std::vector<double> fArrivalTime;
  std::vector<double> fAttenuation;
  std::vector<double> fPolarization;
  std::vector<double> fDirectivity;
  std::vector<std::vector<double>> fSegmentCoord;
  // Length, xpos, ypos, zpos, R_TX, R_RX

  std::vector<double> fDuration; // [ns]
  std::vector<double> fWaveform; // [V/m]
  std::vector<double> fPower;    // [W]
  std::vector<double> fRCS;      // [m^2]


  std::vector<std::vector<double>> fTCSTime;
  std::vector<std::vector<double>> fRCSTime;
  std::vector<std::vector<double>> fPhaseTime;
  std::vector<std::vector<double>> fWaveformTime;


};


class Cascade1D: public Scatter {
public:

  Cascade1D(Antenna& tx, Antenna& rx, Cascade& cs);


  // Accesors

  std::vector<std::vector<double>> Radius();
  std::vector<std::vector<double>> Density();
  double Rwaist();
  std::vector<double> Rcrit();
  std::vector<std::vector<double>> PlasmaFreq();
  std::vector<std::vector<double>> Absorption();
  std::vector<std::vector<double>> SkinDepth();
  std::vector<std::vector<double>> Reflectance();
  std::vector<std::vector<double>> Opacity();
  std::vector<double> TCS();

  // std::vector<double> radar_cs();


private:

  double fDamping;

  std::vector<std::vector<double>> fCSLength;
  std::vector<std::vector<double>> fCSRadius;

  std::vector<std::vector<double>> fDensity;
  double fRwaist;
  std::vector<double> fRcrit; // [cm]
  std::vector<std::vector<double>> fPlasmaFrequency;
  std::vector<std::vector<double>> fAbsorption;
  std::vector<std::vector<double>> fSkinDepth;
  std::vector<std::vector<double>> fReflectance;
  std::vector<std::vector<double>> fTransparency;
  std::vector<std::vector<double>> fOpacity;

  void SetCascadeCoodinates();
  void Density(const std::vector<std::vector<double>> &fCSLength,
                  const std::vector<std::vector<double>> &fCSRadius );
  void Rcrit(const std::vector<std::vector<double>> &density, const double & freq);
  void PlasmaFreq(const std::vector<std::vector<double>> &density);
  void Absorption(const std::vector<std::vector<double>> &density, const double & freq);
  void SkinDepth(const std::vector<std::vector<double>> &density, const double & freq);
  void Reflectance(const std::vector<std::vector<double>> &absorption);
  void Opacity(const std::vector<std::vector<double>> &absorption);
  void TCS(const std::vector<std::vector<double>> &reflectance);

};

class Cylinder1D: public Scatter{
public:

  Cylinder1D(Antenna& tx, Antenna& rx, Cascade& cs);
};

class Line1D: public Scatter{
public:

  Line1D(Antenna& tx, Antenna& rx, Cascade& cs);
};

// std::vector<Scatter> run_scatter_events(Detector det, std::vector<Cascade> cascade_list){
//   std::vector<Scatter> event_list;
//   // For every cascade
//   for (auto& cs: cascade_list){
//     // For every transmitter
//     for (auto& tx : det.get_transmitters()){
//       // For every RX
//       for (auto& rx : det.get_receivers()){
//         // Make the bistatic event.
//         event_list.push_back(Scatter(tx,rx,cs));
//       }
//     }
//   }
// }


#endif
