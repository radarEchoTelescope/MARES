#ifndef SCATTER1D
#define SCATTER1D

#include "scatter.hh"

class Scatter1D: public Scatter {
    friend class Cascade1D;
public:
// The constructor(s) makes a Scatter object without changes
    Scatter1D(Antenna& tx, Antenna& rx, std::vector<ScatterPoint> points,
                const double sampling );

    Scatter1D(Antenna& tx, Antenna& rx, const double sampling);

protected:

/*
Set the points along a direction, starting from a vertex position, 
equally spaced over a given length.

// To get rid of artifacts in the FFT, the positions can be given a small
    non-uniform randomisation along the direction of propagation.
    [Default seed = 42]
    [No randomisation = 0]
*/
    void SetInDirection( std::vector<double> vertex,
                std::vector<double> direction,
                double length, double rand_seed = 42);


/* 
Set the points in the lab frame using a set of coordinates given in the cascade frame
{l_vals, r_vals}. 

So far, the plane R_TX, L_TX was created to cover the cascade at an angle and
the values within were rotated to find the values in the cascade frame (density, etc).
For geometry, a point P with coordinates (r,l) in the TX frame is placed in the
3D lab frame at:

P = CS_vertex + r*e(R_TX) +l*e(L_TX) (capitals == vector)

The choice of L_TX is the direction normal to R_TX that is consistent with the
rotation matrix used in SetTX(), which is the clockwise rotation
of the TX frame with positive angle, or the cascade rotating anti-clockwise from
the frame with positive delta.
*/
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