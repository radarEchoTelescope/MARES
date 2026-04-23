// Calculates the return power for a bi-static radar setup
// Branched of power.C at 23/1/2020
// Enrique Huesca Santiago, 10-2019
// Original FORTRAN code by Krijn D. de Vries 20-10-2019

#ifndef SCATTER
#define SCATTER

// #include "macro_settings.hh"
// #include "cascade.hh"
#include "antenna.hh"
#include "IceRayTracing.hh"
#include "cascade.hh"
// #include "cascade1D.hh"
// #define NDEBUG     // Turn off debug.

struct ScatterPoint{

  // These values should be set at the creation of the point
  int ID;                                        // An identifier
  double TCS;
  std::vector<double> Size{1,1,1};

  // These values are set first, according to number of dimensions
  double L;                                // Distance to the interaction vertex
  std::vector<double> Position{0,0,0};

  // These values are set later, according to the mode of propatation of light. 
  std::vector<double> TXDir{0, 0, 0};      // Vector direction to cs point.
  std::vector<double> RXDir{0, 0, 0};      // Vector direction to cs point.
  double RTX;                              // Module of distance to cs point.
  double RRX ;                             // Module of distance to cs point.
  double Phase;
  double StartTime;
  double ArrivalTime;
  double TXArrivalTime;
  double Attenuation;
  double Directivity;
  double PolEff;
  std::vector<double> Polarization{0,0,0};
  std::vector<double> EFieldAtRX{0,0,0};

// These are needed for raytracing purposes (IceRayTracing)
  double TXRayTime[2];
  double TXRayDistance[2];
  double TXRayStartAngle[2];
  double TXRayEndAngle[2];
  double TXRayAttenuation[2];

  double RXRayTime[2];
  double RXRayDistance[2];
  double RXRayStartAngle[2];
  double RXRayEndAngle[2];
  double RXRayAttenuation[2];
};

// Possible adittions to scatterpoint that are not needed now
  // std::vector<double> fSphericalAngles{0,0};      // Theta and phi in lab frame.

  // double fDot = 0;             // Dot (inner) product with reference (cascade) direction.
  // // This is the same as the cosine of the angle between them.
  // double fDelta = 0;          // The angle.
  // // fDelta is the projection angle, defined between 0 and pi only.



/* Simple scatter event with bistatic configuration*/
class Scatter {
friend class Cascade1D;
public:

  void UpdateTCS(double time);

  Scatter(Antenna& tx, Antenna& rx, std::vector<ScatterPoint> points,
          const double& sampling = _sampling);

  Scatter(Antenna& tx, Antenna& rx, const double& sampling = _sampling);

  /* The RX re-setter*/
  void UpdateRX(const Antenna& new_RX);

  void AddPoint(const ScatterPoint& p);
  void AddPoints(const std::vector<ScatterPoint> new_points);

  /* Positioning functions */ 

  void SetInDirection(std::vector<double> vertex,
                      std::vector<double> direction,
                      double length, double rand_seed = 42);
  
  /* SetInDirection places the ScatterPoints along a direction,
   starting from a vertex position, equally spaced over a given length.

  To get rid of artifacts in the FFT, the positions can be given a small
    non-uniform randomisation along the direction of propagation.
    [Default seed = 42]
    [No randomisation = 0]
  */

  void SetInPlane(std::vector<double> vertex,
                    std::vector<double> l_direction,
                    std::vector<double> r_direction,
                    const std::vector<double> & l_vals,
                    const std::vector<double>& r_vals);
 /* SetInPlane places the ScatterPoints in the lab frame
 following an internal frame (the cascade frame usually) given by the 
 directions of the frame and the positions within the frame {l_vals, r_vals}.
 */

 /* Propagation functions */ 

/* The choice of media for propagation affects the effective directions,
 distances and times(?) of proapgation of the radio waves.
 The choice of propagation _might_ change the other properties of the ScatterPoint. 
  Spatial phase, Attenuation, Polarization, Directivity

Careful, because the medium change and refraction means that the relative
  position of the antennas does not correlate with the distance that the light travels.
*/ 

// Simplest model, rays are straight lines, constant n.
  void SetInConstIce();
// Uzair's IRT assumes that the air-ice interface is the plane z = 0. ({0,0,1,0})
  void SetWithIRT();
//The plane should be a vector with coordinates A,B,C,D: Ax+By+Cz+D=0.
  // void SetInBoundary(const std::vector<double> plane, const double& na, const double& nb);
  //The plane should be a vector with coordinates A,B,C,D: Ax+By+Cz+D=0.
  void SetInBeam(const std::vector<double> plane, const double& na, const double& nb);
  
/* Computes the time integral of the interference of the points at the RX*/
  void RunScatter(const bool save2Dmatrices = false);
  void RunScatterFMCW(const bool save2Dmatrices=false);
  void RunEvent(const bool save2Dmatrices = false);
  
  // Accessors
  Antenna TX();
  Antenna RX();
  std::vector<ScatterPoint> Points();

    // For retrieving a particular variable over the full point collection.

  // TODO: ADD COORDINATES ACCESOR
  std::vector<std::vector<double>> Position();
  std::vector<double> Phase();
  std::vector<double> ArrivalTime();
  std::vector<double> Attenuation();
  // TODO GET RID OF DIRECTIVITY
  std::vector<double> Directivity();
  std::vector<double> Polarization();
  std::vector<double> TCS();

  // The time integral results after RunScatter()
  std::vector<double> Duration();
  std::vector<double> Voltage();
  std::vector<double> Power();
  std::vector<double> RCS();
  std::vector<double> TransmitSignal();
  
  std::vector<std::vector<double>> Phase_time();
  std::vector<std::vector<double>> RCS_time();
  std::vector<std::vector<double>> E_time();
  std::vector<std::vector<double>> AngularFreq_time();


  void save_output_files(const std::string& output_path, const std::array<bool, 17>& flags = {0});
  /*In order to save the member functions, you need to 
  pass an array of flags choosing what members to save:

    True = Save this field
    False = Don't save this field

  0 = Duration
  1 = Voltage
  2 = Power
  3 = TCS
  4 = RCS
  5 = Position / Coordinates
  6 = Phase
  7 = Arrival_time 
  8 = Attenuation
  9 = Polarization
  10 = Phase_time
  11 = RCS_time
  12 = E/V_time
  */

protected:

  Antenna fTX;
  Antenna fRX;
  std::vector<ScatterPoint> fPoints;
  double nP = fPoints.size();

  // double tau;
  double sampling_ratio;


private:
  void SetInConstIce(ScatterPoint& p);
// TO-DO Finish adding IRT. 
  void SetWithIRT(ScatterPoint& p);
  // void SetInBoundary(ScatterPoint& p, const std::vector<double> interface_plane,
  //                const double& na, const double& nb);
  void SetInBeam(ScatterPoint& p, const std::vector<double> interface_plane,
                 const double& na, const double& nb);
  
// This function estimates the point in the interface where the reflection happens.
// This does a bisection seach on the lb value that provides an accurate dt.
  double estimate_interface_point( const double& ra, const double& rb,
    const double& lt, const double& n,
    const double& tolerance = 1E-4, const int& n_tries = 1E4  );

  /* TO PLACE THE SCATTER EVENT IN A BEAM TEST

Coordiantes assigned to match RS/GEANT convention. 
WARNING: This does not match IRT coordinates???
  z as the direction of the beam.
  x as the perpendicular to the beam in the plane where RX lives.
  y as the vertical == 0.

    (P) - - - - - - - - - - - (lt) - - - - - > beam (l direction = z beam direction)
    |   \                      |
    |(rb,lb)  \                |     (B == "Beam" medium)
    |------------I------------------ (interface)
    |              \   (ra, la)|     (A == "Air" medium)
    |                 \        |
    |                    \     |
    (rt) -------------------(RX)
    |
    |
    (r direction == x beam direction)

  rb is the distance between the beam line and the interface.
  (The size of the beam target medium)
  ra is the distance between the interface and the receiver.
  rt is the total distance between the point and the antenna
    in the "vertical" direction of reflection.
  
  We need to figure out the la, lb values to know where the reflection happens,
  and what are the total path lengths of the rays.
*/


  std::vector<double> fDuration;  // [ns]
  std::vector<double> fVoltage;   // [V/m]
  std::vector<double> fPower;     // [W]
  std::vector<double> fRCS;       // [m^2]
  std::vector<double> fTransmitSignal; // a reference signal so unit is not important as it will be normalised when doing the FMCW analysis

  std::vector<std::vector<double>> fRCSTime;
  std::vector<std::vector<double>> fPhaseTime;
  std::vector<std::vector<double>> fAngularFreqTime;
  std::vector<std::vector<double>> fVoltageTime;
};

#endif
