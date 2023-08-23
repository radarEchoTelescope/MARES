#ifndef SCATTER1D
#define SCATTER1D

#include "scatter.hh"

// TO-DO: Finish cleanup this class

class Scatter1D: public Scatter {
    friend class Cascade1D;
public:
// The constructor(s) makes a Scatter object without changes
    Scatter1D(Antenna& tx, Antenna& rx, std::vector<ScatterPoint> points,
                const double sampling );

    Scatter1D(Antenna& tx, Antenna& rx, const double sampling);

protected:

// TO-DO: This needs to know np somehow
    void SetInDirection( std::vector<double> vertex,
                std::vector<double> direction,
                double length, double rand_seed = 42);
// [m] distance from the shower head (starting point)
// Giving a non-uniform position in the cascade gets rid of artifacts in the FFT.
// RN_uniform rand_line(-0.5, 0.5, 42);          // Same seed for debugging.
// RN_uniform rand_line(-0.5, 0.5, time(0));  // Different seed for random, independent runs.

    void SetInPosition(std::vector<double> vertex,
                    std::vector<double> direction,
                    std::vector<double> & l_vals,
                    const std::vector<double>& r_vals);

    void SetInMax(std::vector<double> vertex,
                    std::vector<double> direction, 
                    const double& dL, const double& dR,
                    const std::vector<std::vector<double>> &reflectance);


};

#endif