#include "scatter.hh"

void Scatter::SetInConstIce(){ for(auto& p: fPoints) {SetInConstIce(p);} }

void Scatter::SetInConstIce(ScatterPoint& p){
// Set distances
  p.TXDir = direction( fTX.Pos(), p.Position) ;
  p.RTX = norm(p.TXDir);

  p.RXDir = direction(p.Position, fRX.Pos());
  p.RRX = norm(p.RXDir);

/* Set times
    // (Retarded) time where the scattered radio signal by the segment is produced.
    // production = birth - fRTX/c_ice;

    // (Advanced) time where the scattered signal by the segments arrives in the receiver.
    // arrival = birth + Rr/c_ice;

    // T0 = 0 by definition when the cascade begins (first element = head).
    // p.ArrivalTime = p.L/c_vac + p.RRX/c_ice;
    // time where the segment starts scattering.
*/
  p.ArrivalTime = p.StartTime + p.RRX/c_ice;
  // std::cout<<"c_ice in propagation calc: "<<c_ice<<std::endl;
  // std::cout<<"_c_ice in propagation calc: "<<_c_ice<<std::endl;
  // std::cout<<"lifetime: "<<_lifetime<<std::endl;

  // Set phase
  p.Phase = fTX.Wavenumber()*(p.RTX + p.RRX) - pi/2;

  /*
  tan^-1(f_coll/f_TX) is the phase shift in the oscillation caused by collisions.
  This would a term: 
   + atan2( f_coll, - fTX.Freq() );
  in our regime (f_Coll >>> f_TX), it's almost perfeclty -pi/2.
  */ 

  // E field attenuation at reciever from position dependant factors:
  // e^-r/Latt from medium attenuation
  p.Attenuation =  pow(e, -(p.RTX + p.RRX)/(att_length) ) ;

  /* Attenuation model example
  // Parametrized attenuation length for the Ross Ice Shelf, South Pole.
  double att(double freq_obs){
    double a1=469;                  // [m] Attenuation length parameter
    double a2=-0.205;               // Attenuation length parameter
    double a3=4.87E-5;              // Attenuation length parameter
    return a1+a2*freq_obs/1E6+a3*pow(freq_obs/1E6,2);
  }
  */

  // Segment's polarization direction == E field at segment.
  p.Polarization = cross_product(normalize(p.TXDir), cross_product(normalize(p.TXDir), fTX.Pol() ) );

  // Direction of electric field at the receiver.
  p.EFieldAtRX   = cross_product(normalize(p.RXDir), cross_product(normalize(p.RXDir), p.Polarization) );
  
  // Directivity is hardcoded as a small Herztian dipole.
  p.Directivity = 1.5*norm(p.EFieldAtRX);
   // 3/2*sin(theta_R)*sin(theta_T)

  // Set polarization efficiency
  p.PolEff = abs(projection_pol((p.EFieldAtRX), fRX.Pol()));
   // |sin(theta_R)*sin(theta_T)*cos(theta_R)|
  //In the thin-wire theory, only the  component of  the electric- field vector parallel to the wire
  // axis can interact to form a scattered  wave. That is not our case, our layers will scatter as a free charge

}

void Scatter::SetWithIRT(){ for(auto& p: fPoints) {SetWithIRT(p);} }

void Scatter::SetWithIRT(ScatterPoint& p) {
  /* code call to IRT::GetRayTracing solutions */

  int IgnoreCh[2]={0,0};
  double A0=1;
  double frequency=fTX.Freq();
  // IRT default value for a non-valid ray is -1000
  double RayTime[2]             ={-1E3,-1E3};
  double RayPath[2]             ={-1E3,-1E3};
  double LaunchAngle[2]         ={-1E3,-1E3};
  double RecieveAngle[2]        ={-1E3,-1E3};
  double IncidenceAngleInIce[2] ={-1E3,-1E3};
  double AttRay[2]              ={-1E3,-1E3};

  // double x0=0;/////always has to be zero
  // double x1=sqrt(pow(TxCor[0]-RxCor[0],2)+pow(TxCor[1]-RxCor[1],2));
  // x1 is the distance in the xy plane. 

  //**  IRT expects distances in m (angles in deg) - needs conversion from MARES [mm] units
  double x1=distance(fTX.Pos()[0], fTX.Pos()[1], p.Position[0], p.Position[1])/m;
  double z0=fTX.Pos()[2]/m;
  double z1=p.Position[2]/m;

  // cout<<"dL, TCS "<<p.L<<" "<<p.TCS<<endl;
  // cout<<"CScoords "<<p.Position[0]/m<<" "<<p.Position[1]/m<<" "<<p.Position[2]/m<<endl;
  // cout<<"TXcoords "<<fTX.Pos()[0]/m<<" "<<fTX.Pos()[1]/m<<" "<<fTX.Pos()[2]/m<<endl;
  // cout<<"RXcoords "<<fRX.Pos()[0]/m<<" "<<fRX.Pos()[1]/m<<" "<<fRX.Pos()[2]/m<<endl;

  // cout<<"TX-CS IRT input, z0: "<<z0<<"  x1: "<<x1<<"  z1: "<<z1<<endl;
  
  // This returns the two optimal solutions: Any of the direct, the reflected or two refracted rays. 
  // It will return the shortest ray first, and a second ray if it can find it. 
  IceRayTracing::GetRayTracingSolutions(z1, x1, z0, RayTime, RayPath, 
          LaunchAngle, RecieveAngle, IgnoreCh, IncidenceAngleInIce, A0, frequency, AttRay);

   // void IceRayTracing::GetRayTracingSolutions(double RxDepth, double Distance, double TxDepth, double TimeRay[2], double PathRay[2], 
   // double LaunchAngle[2], double RecieveAngle[2], int IgnoreCh[2], double IncidenceAngleInIce[2], double A0, double frequency, double AttRay[2]){

  // Ignore channel, 0 = direct ray, 1 is refracted ray.
  if(IgnoreCh[0] != 0){
    p.TXRayTime[0]        = RayTime[0];
    p.TXRayDistance[0]    = RayPath[0];
    p.TXRayStartAngle[0]  = LaunchAngle[0];
    p.TXRayEndAngle[0]    = RecieveAngle[0];
    p.TXRayAttenuation[0] = AttRay[0];
  } else {

    cout<<"CScoords "<<p.Position[0]/m<<" "<<p.Position[1]/m<<" "<<p.Position[2]/m<<endl;
    cout<<"TXcoords "<<fTX.Pos()[0]/m<<" "<<fTX.Pos()[1]/m<<" "<<fTX.Pos()[2]/m<<endl;
    cout<<"MARES IRT input [m] "<<z0<<" "<<x1<<" "<<z1<<endl;
    cout<<"IRTout: "<<RayTime[0]*s<<" "<<RayPath[0]<<" "<<LaunchAngle[0]<<" "<<RecieveAngle[0]<<" "<<AttRay[0]<<endl;
    cout<<"ignoreCh[1]: "<<IgnoreCh[1]<<endl;

    std::cerr << "IceRayTracing could not find a ray for TX" << std::endl;
    exit(EXIT_FAILURE); 
  }
// The exit condition can be safely removed if necessary (Multi-RX setups). 

  // 0 = refracted or 1 = reflected. 
  if(IgnoreCh[1] != 0) {
    p.TXRayTime[1]        = RayTime[1];
    p.TXRayDistance[1]    = RayPath[1];
    p.TXRayStartAngle[1]  = LaunchAngle[1];
    p.TXRayEndAngle[1]    = RecieveAngle[1];
    p.TXRayAttenuation[1] = AttRay[1];
  }

  // cout<<"IRT TXraytime: "<<p.TXRayTime[0]<<"  RTX: "<<p.TXRayDistance[0]<<"  TXangle0: "<<p.TXRayStartAngle[0]<<endl;
  // cout<<"TXangle1: "<<p.TXRayEndAngle[0]<<"  TXatten: "<<p.TXRayAttenuation[0]<<endl;
    
// And now, we do the same for the RX

// Reset the default values
  std::fill(IgnoreCh, IgnoreCh+1, 0);
  std::fill(RayTime, RayTime+1, -1E3);
  std::fill(RayPath, RayPath+1, -1E3);
  std::fill(LaunchAngle, LaunchAngle+1, -1E3);
  std::fill(RecieveAngle, RecieveAngle+1, -1E3);
  std::fill(IncidenceAngleInIce, IncidenceAngleInIce+1, -1E3);
  std::fill(AttRay, AttRay+1, -1E3);
  
  x1 = distance(p.Position[0], p.Position[1], fRX.Pos()[0], fRX.Pos()[1])/m;
  z0 = p.Position[2]/m;
  z1 = fRX.Pos()[2]/m;

  // cout<<"CScoords "<<p.Position[0]<<" "<<p.Position[1]<<" "<<p.Position[2]<<endl;
  // // cout<<"TXcoords "<<fTX.Pos()[0]<<" "<<fTX.Pos()[1]<<" "<<fTX.Pos()[2]<<endl;
  // cout<<"RXcoords "<<fRX.Pos()[0]<<" "<<fRX.Pos()[1]<<" "<<fRX.Pos()[2]<<endl;

  // cout<<"CS-RX IRT input, z0: "<<z0<<"  x1: "<<x1<<"  z1: "<<z1<<endl;
  
  IceRayTracing::GetRayTracingSolutions(z1, x1, z0, RayTime, RayPath, 
          LaunchAngle, RecieveAngle, IgnoreCh, IncidenceAngleInIce, A0, frequency, AttRay);

  if(IgnoreCh[0] != 0){
    p.RXRayTime[0]        = RayTime[0];
    p.RXRayDistance[0]    = RayPath[0];
    p.RXRayStartAngle[0]  = LaunchAngle[0];
    p.RXRayEndAngle[0]    = RecieveAngle[0];
    p.RXRayAttenuation[0] = AttRay[0];
  } else {
    cout<<"CScoords "<<p.Position[0]/m<<" "<<p.Position[1]/m<<" "<<p.Position[2]/m<<endl;
    cout<<"RXcoords "<<fRX.Pos()[0]/m<<" "<<fRX.Pos()[1]/m<<" "<<fRX.Pos()[2]/m<<endl;
    // cout<<"MARES IRT input [m] "<<z0<<" "<<x1<<" "<<z1<<endl;

    cout<<"IRTout: "<<RayTime[0]<<" "<<RayPath[0]<<" "<<LaunchAngle[0]<<" "<<RecieveAngle[0]<<" "<<AttRay[0]<<endl;

    std::cerr << "IceRayTracing could not find a ray for RX" << std::endl;
    exit(EXIT_FAILURE); 
  }
// The exit condition can be safely removed if necessary (Multi-RX setups). 

  if(IgnoreCh[1] != 0) {
    p.RXRayTime[1]        = RayTime[1];
    p.RXRayDistance[1]    = RayPath[1];
    p.RXRayStartAngle[1]  = LaunchAngle[1];
    p.RXRayEndAngle[1]    = RecieveAngle[1];
    p.RXRayAttenuation[1] = AttRay[1];
  }

  // cout<<"IRT RXraytime: "<<p.RXRayTime[0]<<"  RRX: "<<p.RXRayDistance[0]<<"  RXangle0: "<<p.RXRayStartAngle[0]<<endl;
  // cout<<"RXangle1: "<<p.RXRayEndAngle[0]<<"  RXatten: "<<p.RXRayAttenuation[0]<<endl;

  // TRIPLE CHECK IRT UNITS, IT SHOULD BE METERS AND DEGREES. 

  p.TXDir = direction( fTX.Pos(), p.Position) ;
  p.RTX = norm(p.TXDir);
  // p.RTX = p.TXRayDistance[0]*m;

  p.RXDir = direction(p.Position, fRX.Pos());
  p.RRX = norm(p.RXDir);
  // p.RRX = p.RXRayDistance[0]*m;

  p.ArrivalTime = p.RXRayTime[0]*s;
  // p.ArrivalTime = p.StartTime + p.RRX/c_ice;

  cout<<"RTX: "<<p.RTX<<"  RRX: "<<p.RRX<<endl;
  // // TXDir: TX -> CS (x1 - x0), (y1 - y0), (z1 - z0)
  cout<<"TXDir: "<<p.TXDir<<endl;
  cout<<"RXDir: "<<p.RXDir<<endl;
  cout<<"TXDir n: "<<normalize(p.TXDir)<<endl;
  cout<<"RXDir n: "<<normalize(p.RXDir)<<endl;
  // cout<<"IRT RX arrival time: "<<p.ArrivalTime<<endl;

  // TODO TRIPLE CHECK THAT THIS IS STILL VALID
  // Set phase
  // According to uzair,
  // Phase = k_vac*Optical path length = k_vac*propagation time*c_vac
  // K_vac*RayPath
  // Waiting for confirmation from Krijn 
  
  // p.Phase = fTX.Wavenumber()*(p.RTX + p.RRX)
  //               + atan2( f_coll, - fTX.Freq() );

  p.Phase = fTX.Wavenumber()*(p.RTX + p.RRX) - pi/2;
  // tan^-1(f_coll/f_TX) is the phase shift in the oscillation caused by collisions.
  // in our regime, it could be almost fixed to -Pi/2.

  // DOUBLE CHECK MINUS SIGN

// DOES IRT INCLUDE ICE ATTENUATION?
  // E field attenuation at reciever from position dependant factors:
  // e^-r/Latt from medium attenuation
  // p.Attenuation =  pow(e, -(p.RTX + p.RRX)/(2*att_length) ) ;
  // cout<<"IRT atten: "<<p.Attenuation<<endl;

  p.Attenuation = p.TXRayAttenuation[0]*p.RXRayAttenuation[0];
  // cout<<"IRT atten: "<<p.Attenuation<<endl;

  // Segment's polarization direction == E field at segment.
  p.Polarization = cross_product(normalize(p.TXDir), cross_product(normalize(p.TXDir), fTX.Pol() ) );

  // Direction of electric field at the receiver.
  p.EFieldAtRX   = cross_product(normalize(p.RXDir), cross_product(normalize(p.RXDir), p.Polarization) );

  // Set polarization efficiency
  p.PolEff = abs(projection(p.EFieldAtRX, fRX.Pol()));
  //In the thin-wire theory, only the  component of  the electric- field vector parallel to the wire
  // axis can interact to form a scattered  wave. That is not our case, our layers will scatter as a free charge
  cout<<"test: "<<p.ArrivalTime<<endl;


  // Directivity is hardcoded as a small Herztian dipole.
  // p.Directivity = 1.5*norm(p.EFieldAtRX);
   // 3/2*sin(theta_R)*sin(theta_T)
}

// void Scatter::SetInBoundary( const std::vector<double> plane, const double& na, const double& nb){
//   // If you wish to make a lookup table, it should be made here.
//    for(auto& p: fPoints) {SetInBoundary(p, plane, na, nb);}
//  }


// void Scatter::SetInBoundary(ScatterPoint& p, Antenna& at, const std::vector<double> interface_plane,
//                  const double& na, const double& nb){
//   // Temporary variables: a is the antenna's medium and b is the point's medium. 
//   double ka, kb, ra, rb, la, lb, lt;
//   double n  = nb/na;

//   // Distance vector from the antenna to the point
//   std::vector<double> at_dir = direction(at.Pos(), p.Position);

//   // Store the direction properly ?? 
//   (at.Power() != 0) ? p.TXDir = at_dir : p.RXDir = -1.* at_dir ;

//     // Direction normal to interface plane.
//   std::vector<double> e_n = {interface_plane[0], interface_plane[1], interface_plane[2]};
//   e_n = normalize(e_n);

//   // Find the distances in the plane of reflection to the interface plane. 
//     // Shortest distance from the TX antenna to the plane.  
//   ra = shortest_distance(interface_plane, fTX.Pos() );
//     // Shortest distance fom the point p to the plane.
//   rb = shortest_distance(interface_plane, p.Position);


//   // lt, the total distance in the plane is basic trigonometry. 
//   lt = sqrt(pow( norm(p.TXDir),2) - (pow(ra + rb,2)) );
//   // Estimate l coordinate for interaction point 
//   la = estimate_interface_point( ra, rb, lt, n);
//   lb = lt - la ;
//   // Distances for paths in both media
//   ka = sqrt(pow(ra,2) + pow(la,2));
//   kb = sqrt(pow(rb,2) + pow(lb,2));

    
//   // Total distance travelled by the TX ray of light
//   p.RTX = ka + kb;

//     // The incidence angle is:
//   double theta_i  = atan2(la,ra);
//     // The refracted angle is:
//   double theta_r  = atan2(lb,rb);
//   p.TXDir = normalize(p.TXDir);

//   /* The normal to the interface plane must be alinged with the direction of indicence.
//   in order to have a consistent geometry frame */

//   if( sgn(projection(p.TXDir, e_n)) < 0 ) { e_n *=-1.; } 
  

//   // Direction perpendicular to the plane of incidence. 
//   std::vector<double> e_s = normalize(cross_product(p.TXDir, e_n));
//   /* If TXDir is perpendicular to the interface incidence, 
//     the plane of incidence is degenerate. 
//     The model is symmetric and we can pick any plane we want. 
//   */
//   if ( norm(e_s) == 0 ){
//     e_s = normalize(cross_product(find_perpendicular(e_n), e_n));
//   }
//   // Direction shared by the plane of incidence and the interface plane
//   std::vector<double> e_l = normalize(cross_product(e_n, e_s));
//   // Vector direction to interface point
//   std::vector<double> e_ka = normalize(la*e_l + ra*e_n);
//   // The direction of polarization in the plane.
//   std::vector<double> e_p1 = normalize(cross_product(e_s, e_ka));
//   // The direction to the point in after the interface is:
//   // e_kb  = rotate(e_ka, e_s, pi + theta_i - theta_r);
//   // Rotation of the plane component due to refraction;
//   std::vector<double> e_p2 = rotate(e_p1, e_s, pi + theta_i - theta_r );

//   // sine of the incidence angle, e_n and k are normalized already. 
//   double s_i = norm(cross_product(e_ka, e_n));
//   // cosine of the incidence angle
//   double c_i = projection(e_n, e_ka);
//   // Check that sin^2 and cos^2 give 1?

//   // Frenel coefficient for the perpendicular component
//   double sqrt_expr = sqrt( pow(n,2) - pow(s_i,2) );
//   double F_s = (2.*c_i)/( c_i + sqrt_expr );
//   // Frenel coefficient for the parallel component; 
//   double F_p = (2.* nb *c_i)/(n*nb*c_i + na * sqrt_expr );


//   // The radiation field BEFORE the interface
//   std::vector<double> E_rad = cross_product( e_ka, cross_product( e_ka, fTX.Pol() ));
//   // E_rad is contains some polarization efficency term.
//   p.PolEff = abs(norm(E_rad)); // == abs( norm(cross_product(e_ka, fTX.Pol())) )


//   // The radiation field AFTER the interface
//   // The component in the plane of the interface (perpendicular)
//   std::vector<double> E_s = projection(E_rad, e_s)*F_s*e_s;
//   // The component in the plane of incidence (parallel)
//   std::vector<double> E_p = projection(E_rad, e_p2)*F_p*e_p2;
// }

void Scatter::SetInBeam( const std::vector<double> plane, const double& na, const double& nb){
  // If you wish to make a lookup table, it should be made here.
   for(auto& p: fPoints) {SetInBeam(p, plane, na, nb);}
 }


// This SetInBeam sets first all TX parameters and then all RX parameters. 
void Scatter::SetInBeam(ScatterPoint& p, const std::vector<double> interface_plane,
                        const double& na, const double& nb){
  double n  = nb/na;
  double ra, rb, la, lb, lt;
  double k1,k2,k3,k4;

    // Direction normal to interface plane.
  std::vector<double> e_n = {interface_plane[0], interface_plane[1], interface_plane[2]};
  e_n = normalize(e_n);

  // Find the distances in the plane of reflection to the interface plane. 
    // Shortest distance from the TX antenna to the plane.  
  ra = shortest_distance(interface_plane, fTX.Pos() );
    // Shortest distance fom the point p to the plane.
  rb = shortest_distance(interface_plane, p.Position);

  // Distance vector from the point to the transmitter
  p.TXDir = (direction(fTX.Pos(), p.Position));

  // lt, the total distance in the plane is basic trigonometry. 
  lt = sqrt(pow( norm(p.TXDir),2) - (pow(ra + rb,2)) );
  // Estimate l coordinate for interaction point 
  la = estimate_interface_point( ra, rb, lt, n);
  lb = lt - la ;
  // Distances for paths in both media
  k1 = sqrt(pow(ra,2) + pow(la,2));
  k2 = sqrt(pow(rb,2) + pow(lb,2));
  // Total distance travelled by the TX ray of light
  p.RTX = k1 + k2;

  /* The normal to the interface plane must be alinged with the direction of indicence.
  in order to have a consistent geometry frame */
  p.TXDir = normalize(p.TXDir);
  if( sgn(projection(p.TXDir, e_n)) < 0 ) { e_n *=-1.; } 
  
  // Direction perpendicular to the plane of incidence. 
  std::vector<double> e_s = normalize(cross_product(p.TXDir, e_n));
  /* If TXDir is perpendicular to the interface incidence, 
    the plane of incidence is degenerate. 
    The model is symmetric and we can pick any plane we want. 
  */
  if ( norm(e_s) == 0 ){
    e_s = normalize(cross_product(find_perpendicular(e_n), e_n));
  }
  // Direction shared by the plane of incidence and the interface plane
  std::vector<double> e_l = normalize(cross_product(e_n, e_s));
  // Vector direction to interface point
  std::vector<double> e_k1 = normalize(la*e_l + ra*e_n);
  // The direction of polarization in the plane.
  std::vector<double> e_p1 = normalize(cross_product(e_s, e_k1));

  // sine of the incidence angle, e_n and k are normalized already. 
  double s_i = norm(cross_product(e_k1, e_n));
  // cosine of the incidence angle
  double c_i = projection(e_n, e_k1);

// Now, we need to rotate the vector once it goes into the dense medium

  // We get the incident and refacted angles from geometry.
    // The incidence angle is:
  double theta_i  = atan2(la,ra);
    // The refracted angle is:
  double theta_r  = atan2(lb,rb);

// Sanity checks time! ------------------------
  //  First, Snell's law must be obeyed
  assert(abs(sin(theta_i) - n*sin(theta_r)) < 1e-5 && "Snell's law is not conserved");
  // Second, sine and cosine from both definitions have to match
  assert(abs(s_i - sin(theta_i)) < 1e-5 && "Sine of incoming angle mismatch");
  assert(abs(c_i - cos(theta_i)) < 1e-5 && "Cosine of incoming angle mismatch");
  // Third, check that sin^2 + cos^2 = 1 (why not)
  assert(abs(pow(s_i,2) + pow(c_i,2) - 1) < 1e-5 && "Math relations not conserved");
// --------------------------------------------

  // The direction to the point in after the interface is:
  // e_k2  = rotate(e_k1, e_s, pi + theta_i - theta_r);
  // Rotation of the plane component due to refraction;
  std::vector<double> e_p2 = rotate(e_p1, e_s, pi + theta_i - theta_r );

  // Frenel coefficient for the perpendicular component
  double sqrt_expr = sqrt( pow(n,2) - pow(s_i,2) );
  double F_s = (2.*c_i)/( c_i + sqrt_expr );
  // Frenel coefficient for the parallel component; 
  double F_p = (2.* nb *c_i)/(n*nb*c_i + na * sqrt_expr );

  // The radiation field BEFORE the interface
  std::vector<double> E_rad = cross_product( e_k1, cross_product( e_k1, fTX.Pol() ));
  // MARES stores the amplitude of the polarization efficiency (mismatch) separately. 
  p.PolEff = abs(norm(E_rad)); // == abs( norm(cross_product(e_k1, fTX.Pol())) )
  // = |sin(theta_t)|

  // We remove it from the field, so we're not double-counting it. 
  E_rad = normalize(E_rad);

  // The radiation field AFTER the interface
  // The component in the plane of the interface (perpendicular)
  std::vector<double> E_s = projection(E_rad, e_s)*F_s*e_s;
  // The component in the plane of incidence (parallel)
  std::vector<double> E_p = projection(E_rad, e_p2)*F_p*e_p2;

  // We store the attenuation factor in its own variable too. 
  // The attenuation in this case is due to the Frenel coefficients,
  // not medium propagation (negligible in small baselines).
  p.Attenuation = sqrt( pow( norm(E_s),2) + pow( norm(E_p),2) );

   // E field at segment == Segment's polarization direction.
  p.Polarization = normalize(E_p + E_s);

  // And now we perform the same steps but from the p to the RX. 
  // CAREFUL! The geometry is inverted w.r.t the TX's case above.
  // We need make sure that the rays go from the
  // from "new" (dense) medium to "original" (air) medium now.

  // Find the distances in the plane of reflection to the interface plane. 
    // Shortest distance from the TX antenna to the plane.  
  ra = shortest_distance(interface_plane, fRX.Pos() );
    // Shortest distance fom the point p to the plane. (already known)
  // rb = shortest_distance(interface_plane, p.Position);
  p.RXDir = direction(p.Position, fRX.Pos() ) ;

  // lt, the total distance in the plane of incidence. 
  lt = sqrt(pow( norm(p.RXDir),2) - (pow(ra + rb,2)) );
  // Estimate l coordinate for interaction point 
  la = estimate_interface_point( ra, rb, lt, n);
  lb = lt - la ;

  // Distances for paths in both media
  k3 = sqrt(pow(rb,2) + pow(lb,2));
  k4 = sqrt(pow(ra,2) + pow(la,2));
  // Total distance travelled by the RX ray of light
  p.RRX = k3 + k4;

  // Distance vector from the point to the receiver
  p.RXDir = normalize( p.RXDir );

  /* The frame geometry needs to be preserved 
  e_n must always point "outwards".
  The half-plane of the viewer might not be the same as the TX's, 
  but the equations are mantained. 
  */

  if( sgn(projection(p.RXDir, e_n)) < 0 ) { e_n *=-1.; } 
  // Direction perpendicular to the plane of incidence. 
  
  /* If RXDir is perpendicular to the interface incidence, 
    we can keep the e_s that we computed before.
  */
  if ( norm(cross_product(p.RXDir, e_n)) != 0. ){
    e_s = normalize(cross_product(p.RXDir, e_n));
  }

  // Direction shared by the plane of incidence and the interface plane
  e_l = cross_product(e_n, e_s);

  // Vector direction to interface point
  std::vector<double> e_k3 = normalize(lb*e_l + rb*e_n); 
  // This is the part in the dense medium

  // Vector normal to propagation in incidence plane
  std::vector<double> e_p3 = normalize(cross_product(e_s, e_k3));
  // This is the part back in the original medium.  

// Once again, we need to rotate this vector once it gets out of the dense medium.

  // sine of the incidence angle, e_n and k are normalised already. 
  s_i = norm(cross_product(e_k3, e_n));
  // cosine of the incidence angle
  c_i = projection(e_k3,e_n);

  // We get the incident and refacted angles from geometry
    // The incidence angle is:
  theta_i  = atan2(lb,rb);
  // The refracted angle is:
  theta_r  = atan2(la,ra);
  // These definitions are flipped w.r.t the previous case above.
  
  // Sanity checks time again! ---------------------------------
  //  First, Snell's law must be obeyed
  assert(abs(n*sin(theta_i) - sin(theta_r)) < 1e-5 && "Snell's law is not conserved");
  // Second, sine and cosine have to match
  assert(abs(s_i - sin(theta_i)) < 1e-5 && "Sine of incoming angle mismatch");
  assert(abs(c_i - cos(theta_i)) < 1e-5 && "Cosine of incoming angle mismatch");
  // Third, check that sin^2 + cos^2 = 1 (why not)
  assert(abs(pow(s_i,2) + pow(c_i,2) - 1) < 1e-5 && "Math relations not conserved");
  //  -----------------------------------------------------------

  // The direction to the point in after the interface is:
  // e_k4  = rotate(e_k3, e_s, pi + theta_i - theta_r);
  // Rotation of the plane component due to refraction;
  std::vector<double> e_p4 = rotate(e_p3, e_s, pi + theta_i - theta_r );

  sqrt_expr = sqrt( pow(n,2) - pow(s_i,2) );
  // Check for total internal reflection,
  // it only happens when the expression above becomes imaginary
  if (std::isnan(sqrt_expr)) {
    F_s = 0.;
    F_p = 0.;
  } else {
    // Frenel coefficient for the perpendicular component
    F_s = (2.*c_i)/( c_i + sqrt_expr );
    // Frenel coefficient for the parallel component; 
    F_p = (2.* na *c_i)/(n*na*c_i + nb *sqrt_expr );
  }
  
  // The radiation field before the interface:
  E_rad = cross_product( e_k3, cross_product( e_k3, p.Polarization ));
  // We store the second component of the polarization,
  // the radiation field of the Herztian dipole.
  p.PolEff *= abs(norm(E_rad)); // == abs( norm(cross_product(e_k3, p.Pol)) )
  //  = |sin(theta_t)*sin(theta_r)|

  // We remove (again) it from the field.
  E_rad = normalize(E_rad);

  // The component in the plane of the interface (perpendicular)
  E_s = projection(E_rad, e_s)*F_s*e_s;
  // The component in the plane of incidence (parallel)
  E_p = projection(E_rad, e_p4)*F_p*e_p4;

  // The total attenuation includes the Frenel coefficients of the two processes.
  p.Attenuation *= sqrt( pow( norm(E_s),2) + pow( norm(E_p),2) );


  // Direction of electric field at the receiver.
  p.EFieldAtRX = normalize(E_p + E_s);

  // Third component of the polarization efficiency
  p.PolEff *= abs(projection(p.EFieldAtRX, fRX.Pol()));
  // |sin(theta_t)*sin(theta_r)*cos(theta_r)|
  
  // Set ArrivalTime
  p.ArrivalTime = p.StartTime + (nb*k3 + na*k4)/c_vac;

  // Set phase
  p.Phase = fTX.Wavenumber()*(na*k1 + nb*k2 + nb*k3 + na*k4)
                + atan2( f_coll, - fTX.Freq() );
  // tan^-1(f_coll/f_TX) is the phase shift in the oscillation caused by collisions.
  // in our regime, it could be almost fixed to -Pi/2.
}

double Scatter::estimate_interface_point( const double& ra, const double& rb,
                             const double& lt, const double& n,
                             const double& tolerance, const int& n_tries ){
  int steps = 0;
  double da, tmp, lt_guess= 0, la_guess,  lb_guess, l_start = 0, l_end = lt;
  do {
    la_guess = (l_start + l_end)/2.;
    // db = distance between interface point guess and p
    da = sqrt(pow(la_guess,2)+pow(ra,2));
    tmp = la_guess/ (da*n);
    lb_guess = rb*tmp / sqrt(1 - pow(tmp,2));
    lt_guess = la_guess + lb_guess;
    lt > lt_guess ? l_start = la_guess : l_end = la_guess;
    steps++;
    if (steps > n_tries){
      throw std::length_error("Did not converge in the set number of tries.");
    }
  }
  while (abs(lt-lt_guess) > tolerance);
  return la_guess;
}