#include "scatter1D.hh"

Scatter1D::Scatter1D(Antenna& tx, Antenna& rx, std::vector<ScatterPoint> points,
                const double sampling ):
  Scatter(tx, rx, points, sampling){}  

Scatter1D::Scatter1D(Antenna& tx, Antenna& rx, const double sampling):
  Scatter(tx,rx,sampling){}

// Set the position and arrival time of your segments

void Scatter1D::SetInDirection(std::vector<double> vertex,
                        std::vector<double> direction,
                        double length, double rand_seed){

  double dSeg = length/nP;
  // This does not work
  // if(rand_seed){ RN_uniform rand_line(-0.5, 0.5, rand_seed);}
  RN_uniform rand_line(-0.5, 0.5, rand_seed);

  for(int i = 0 ; i < nP; i++){
    ScatterPoint& p = fPoints[i];
    
    if(rand_seed){ 
      p.L = (i + rand_line.get()) * dSeg;
    } else{ 
      p.L = i * dSeg; 
    }

    p.StartTime = p.L/c_vac;

    // Set position
    p.Position[0] = p.L*direction[0] + vertex[0];
    p.Position[1] = p.L*direction[1] + vertex[1];
    p.Position[2] = p.L*direction[2] + vertex[2];

  }
}


/* So far, the plane R_TX, L_TX was created to cover the cascade at an angle and
the values within were rotated to find the values in the cascade frame (density, etc).
For geometry, a point P with coordinates (r,l) in the TX frame is placed in the
3D lab frame at:

P = CS_vertex + r*e(R_TX) +l*e(L_TX) (capitals == vector)

The choice of L_TX is the direction normal to R_TX that is consistent with the
rotation matrix used in SetTX(), which is the clockwise rotation
of the TX frame with positive angle, or the cascade rotating anti-clockwise from
the frame with positive delta.
*/
void Scatter1D::SetInPosition(std::vector<double> vertex,
                    std::vector<double> direction,
                    std::vector<double> & l_vals,
                    const std::vector<double>& r_vals){


  std::vector<double> R_TX_dir{0,0,0}, L_TX_dir{0,0,0},
                      seg_pos{0,0,0}, tx_seg{0,0,0}, rx_seg{0,0,0};
  R_TX_dir = normalize(fTX.Dir()) ;

// This only breaks for R_TX_dir == this->fDirection, exactly, we can guard against it easily
// A deviation of 10^-20 in the direction is well within numerical errors, and avoids
// making the cross product exactly zero.

  if( ( R_TX_dir == normalize( direction ) ) ) {
    R_TX_dir[0] += 1E-20;
    R_TX_dir[1] += 1E-20;
    R_TX_dir[2] += 1E-20;
  }

  L_TX_dir = cross_product( normalize( cross_product(R_TX_dir, normalize(direction) ) ), R_TX_dir);
  L_TX_dir = normalize(L_TX_dir);

  // To separate between the cascade pointing towards TX and away from TX
  double s = sgn( projection(R_TX_dir, direction) ); // cos(delta)
  if(s < 0){  std::reverse( l_vals.begin(), l_vals.end() ); }
  
  // Segment loop
  // The length of l_vals and r_vals must be nL (nP), as that is the number of
  // segments to be placed.
  ScatterPoint p;
  for(int i = 0 ; i < l_vals.size(); i++){
    p = fPoints[i];

    // Set distance from the shower head (starting point)
    p.Position[0]  = r_vals[i]*R_TX_dir[0]*s + l_vals[i]*L_TX_dir[0];
    p.Position[1]  = r_vals[i]*R_TX_dir[1]*s + l_vals[i]*L_TX_dir[1];
    p.Position[2]  = r_vals[i]*R_TX_dir[2]*s + l_vals[i]*L_TX_dir[2];
    p.L = norm(p.Position);

// TO_DO: insert randomisation too.

  // // Giving a non-uniform position in the cascade gets rid of artifacts in the FFT.
  // RN_uniform rand_line(-0.5, 0.5, 42);          // Same seed for debugging.
  // // RN_uniform rand_line(-0.5, 0.5, time(0));  // Different seed for random, independent runs.

    // Set Position in lab frame (TX)
    p.Position[0] += vertex[0];
    p.Position[1] += vertex[1];
    p.Position[2] += vertex[2];
    }
}


/*  For each segment in the L_TX direction, we can find the coordinates (r_max,l_max)
  of the point where the reflectivity is the highest (it should correlate with the
  peak density for the segment, too), and place the segment at that point.
*/
void Scatter1D::SetInMax(std::vector<double> vertex, std::vector<double> direction, 
                          const double& dL, const double& dR,
                          const std::vector<std::vector<double>> &reflectance){
// i lives in L_TX, j lives in R_TX
// (The integral of reflectance happens along j)
    int i, j;
    double l,r;
    std::vector<double> fTXMaxLength;
    std::vector<double> fTXMaxRadius;
    double nL = reflectance.size();
    fTXMaxLength = std::vector<double> (nL);
    fTXMaxRadius = std::vector<double> (nL);
  for (i = 0; i < nL; i++){
    // We want the index, not the value, therefore distance()
    j = std::distance(reflectance[i].begin(),
            max_element(reflectance[i].begin(), reflectance[i].end() ) );
    fTXMaxLength[i] = i*dL;   // This is "b" before the rotation
    fTXMaxRadius[i] = j*dR;   // This is "a" before the rotation
  }
  
  SetInPosition(vertex, direction, fTXMaxLength, fTXMaxRadius);
}
