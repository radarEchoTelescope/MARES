// Calculates the return power for a bi-static radar setup
// Branched of power.C at 23/1/2020
// Enrique Huesca Santiago, 10-2019
// Original FORTRAN code by Krijn D. de Vries 20-10-2019

#ifndef CASCADE1D
#define CASCADE1D

#include "cascade.hh"
#include "scatter1D.hh"
// #define NDEBUG     // Turn off debug.

class Cascade1D: public Scatter1D, public Cascade {
public:

  Cascade1D(Antenna& tx, Antenna& rx, Cascade& cs,
        const double deltaL = 1, const double deltaR = 0.1, const double deltaN = 1,
        const double sampling = 100);

  void SetInDirection(double rand_seed = 42);

  std::vector<double> TCS();
  std::vector<std::vector<double>> Radius();
  std::vector<std::vector<double>> Length();

private:

// Dimensions of the TX frame axis.
  double L, R;
/*
  If delta = pi/2: Perpendicular incidence
    - L (TX frame) = L (CS frame)
    - R (TX frame) = 2R (CS frame)

  If delta = 0:
    - L (TX frame) = 2R (CS fame)
    - R (TX frame) = L (CS frame)
*/
// The resolution in the frame that integrates the cascade values.
// The free parameter is the step size.
  double dL, dR, dN;
  // nbins  = size / step_size.
  int nL, nR;

// cD = cosine delta, sD, sine delta, the angle between TX_dir and CS_dir.
  double cD, sD;

// The 2D matrices that hold the length and radial values of the cascade in 
// the transmitter frame. 
  std::vector<std::vector<double>> fCSLength;
  std::vector<std::vector<double>> fCSRadius;

// Find the frame points and the respective density for each point.
// This also fills the density matrix fDensity (from Cascade) to save time. 
  void SetTXFrame();

/* The density is computed not in  the lab frame, but in the plane between the
    tx.direction vector and the cs direction vector. Each one defines a frame
    in this plane their own perpendicular direction. The two frame share the
    same normal direction.

    The two frames are separated by an angle delta. This angle delta behaves like
    the declination, is only defined between 0 and pi. Due to the cascade radial
    symmetry, and the choice of plane (that cuts the cascade longitudinally)
    the solutions for delta and minus delta in the plane frame are equivalent.
*/

  std::vector<double> fTCS;
  void SetTCS();

};
#endif

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
