// Calculates the return power for a bi-static radar setup
// Branched of power.C at 23/1/2020
// Enrique Huesca Santiago, 10-2019
// Original FORTRAN code by Krijn D. de Vries 20-10-2019

#ifndef CASCADE1D
#define CASCADE1D

// #include "macro_settings.hh"
#include "cascade.hh"
#include "scatter.hh"
// #define NDEBUG     // Turn off debug.

// TO_DO Missing the code for multi-receiver setups: How to compute TCS only once.

class Cascade1D: public Scatter, public Cascade {
public:

  Cascade1D(Antenna& tx, Antenna& rx, Cascade& cs,
        const double deltaL = 1, const double deltaR = 0.1, const double deltaN = 1,
        const double lifetime = 1E-8, const double sampling = 100);

  std::vector<std::vector<double>> Radius();
  std::vector<std::vector<double>> Length();

  std::vector<double> TCS();

  // // Inner product between CS direction and TX-CS (line of sight).
  // double  Dot();
  // // Angle of the inner product (between 0 and pi).
  // double  Delta();
  // double Rwaist();
  // std::vector<double> Rcrit();


private:
// The resolution in the frame that integrates the cascade values.
// The free parameter is the step size.
  double dL, dR, dN;

// These are the dimensions of the TX frame axis. [cm]
// In cm since everything in the cascade should be in cm.
  double L, R;
/*
  If delta = pi/2: Perpendicular incidence
    - L (TX frame) = L (CS frame)
    - R (TX frame) = 2R (CS frame)

  If delta = 0:
    - L (TX frame) = 2R (CS fame)
    - R (TX frame) = L (CS frame)
*/

// nbins  = size / step_size.
  int nL, nR;

// The three values that you need to find the TCS of a layer are
// Radius (from the core), Density and Transparency of that point.

  std::vector<std::vector<double>> fCSLength;
  std::vector<std::vector<double>> fCSRadius;

  std::vector<double> fTCS;

/* The density is computed not in  the lab frame, but in the plane between the
    tx.direction vector and the cs direction vector. Each one defines a frame
    in this plane their own perpendicular direction. The two frame share the
    same normal direction.

    The two frames are separated by an angle delta. This angle delta behaves like
    the declination, is only defined between 0 and pi. Due to the cascade radial
    symmetry, and the choice of plane (that cuts the cascade longitudinally)
    the solutions for delta and minus delta in the plane frame are equivalent.
*/
// Find the frame points and the respective density for each point.
  void SetTXFrame();

  void TCS(const std::vector<std::vector<double>> &radius,
           const std::vector<std::vector<double>> &density);

  // Before running the scatter, you need to place your scatterers:
// This can be done before passing to Scatter or with SetInLine()

// [m] distance from the shower head (starting point)
  void SetInLine( std::vector<double> vertex,
                std::vector<double> direction,
                double length, double rand_seed = 42);
// Giving a non-uniform position in the cascade gets rid of artifacts in the FFT.
// RN_uniform rand_line(-0.5, 0.5, 42);          // Same seed for debugging.
// RN_uniform rand_line(-0.5, 0.5, time(0));  // Different seed for random, independent runs.

  void SetWRTCSDirection();

  std::vector<double> fTXMaxLength;
  std::vector<double> fTXMaxRadius;
  void SegmentMaxCoords(const std::vector<std::vector<double>> &reflectance);

  void SetWRTCSPosition(std::vector<double> & CSLength,
    const std::vector<double>& CSRadius);


  // void Rcrit(const std::vector<std::vector<double>> &density, const double & freq);
  // double fDamping;
  //
  // std::vector<ScatterPoint> fSegments;
  // double fRwaist;
  // std::vector<double> fRcrit; // [cm]
};

// class Line1D: public Scatter{
  // public:
  //
  //   Line1D(Antenna& tx, Antenna& rx, Cascade& cs);
  // };

// class Cylinder1D: public Scatter{
// public:
//
//   Cylinder1D(Antenna& tx, Antenna& rx, Cascade& cs);
// };

#endif

//
// /* Computes crital radius, the radial values of the line where the plasma frequency
// equals the detection frequency (overdense line).
// The density associated with the plasma frequency is the critical density [(#e-) cm^-3]
// If the local density is higher than that then the medium is overdense.
//
// The shower waist is by definition the maximum radius for the same plasma density.
// */
// void Cascade::SetRcrit(std::vector<std::vector<double>> &density, const double & freq){
// 	double r, dens_crit = memp*pow(freq/8980,2);
// 	fRcrit = std::vector<double>(density.size(), 0);
//
// 	// Loop over density matrix rows and then columms.
// 	for (int i = 0; i < density.size(); i++){
// 		r = 0;
// 		for (int j = 0; j < density[0].size(); j++){
// 			// If overdense, remember the layer.
// 			if (density[i][j] > dens_crit) { r = j;}   // Layers
// 		}
// 		fRcrit[i] = r*fRdiv; // [cm]
// 	}
//
// 	fRwaist = *max_element(fRcrit.begin(), fRcrit.end());
// }

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
