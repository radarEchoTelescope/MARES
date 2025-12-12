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
  in our regime (f_Coll >>> f_TX), it's almost perfectly -pi/2.
  */ 

  // E field attenuation at reciever from position dependant factors:
  // e^-r/Latt from medium attenuation
  p.Attenuation =  pow(e, -(p.RTX + p.RRX)/(2*att_length) ) ;

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
  p.PolEff = abs(projection((p.EFieldAtRX), fRX.Pol()));
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

  // First IRT call, for TX -> CS
  // IRT uses metres, degrees -> requires conversion to and from MARES units
  // double x0=0; // always has to be zero, hardcoded in IRT
  double x1 = distance(fTX.Pos()[0], fTX.Pos()[1], p.Position[0], p.Position[1])/m;  // x1 is the distance in the xy plane b\n the two points ('tx' and 'rx')
  double z0 = fTX.Pos()[2]/m;                                                        // z0 is the 'tx' depth
  double z1 = p.Position[2]/m;                                                       // z1 is the 'rx' depth
  
  // This returns the two optimal solutions: Any of the direct, the reflected or two refracted rays. 
  // It will return the shortest ray first, and a second ray if it can find it. 
  //  cout<<" Running IRT TX -> CS:"<<endl;
  IceRayTracing::GetRayTracingSolutions(z1, x1, z0, RayTime, RayPath, 
          LaunchAngle, RecieveAngle, IgnoreCh, IncidenceAngleInIce, A0, frequency, AttRay);

  // Nan check, in case IRT runs into problems. Can be disabled, but may introduce errors into final output.
  for(int i=0;i<2;i++){
    if((std::isnan(RayTime[i])==true)      || (std::isnan(RayPath[i])==true) || (std::isnan(LaunchAngle[i])==true) ||
       (std::isnan(RecieveAngle[i])==true) || (std::isnan(AttRay[i] )==true)){
      cout<<"Nan found in TX -> CS IRT solutions, for Ch "<<i<<endl;
      cout<<"IRTinput [m], z0: "<<z0<<"   x1: "<<x1<<"  z1: "<<z1<<endl;
      cout<<"IRToutput, time: "<<RayTime[i]*s<<"  path: "<<RayPath[i]<<"   Langle:"<<LaunchAngle[i]<<"   Rangle: "<<RecieveAngle[i]<<"   Attenuation: "<<AttRay[i]<<endl;
      exit(EXIT_FAILURE);
    }
  }

   // void IceRayTracing::GetRayTracingSolutions(double RxDepth, double Distance, double TxDepth, double TimeRay[2], double PathRay[2], 
   // double LaunchAngle[2], double RecieveAngle[2], int IgnoreCh[2], double IncidenceAngleInIce[2], double A0, double frequency, double AttRay[2]){

  // Ignore channel, 0 = direct ray, 1 is refracted ray.
  if(IgnoreCh[0] != 0){
    p.TXRayTime[0]        = RayTime[0]*s;
    p.TXRayDistance[0]    = RayPath[0]*m;
    p.TXRayStartAngle[0]  = LaunchAngle[0]*deg;
    p.TXRayEndAngle[0]    = RecieveAngle[0]*deg;
    p.TXRayAttenuation[0] = AttRay[0];
  } else {
    // std::cerr << "IceRayTracing could not find a ray for TX" << std::endl;
    // exit(EXIT_FAILURE); 

    // (temporary) fix for segments with no solution - we use a the -1E3 value as a flag for 'no solutions found'
    // Note - we no longer stop the simulation if this occurs. 
    p.TXRayTime[0]        = -1E3;
    p.TXRayDistance[0]    = -1E3;
    p.TXRayStartAngle[0]  = -1E3;
    p.TXRayEndAngle[0]    = -1E3;
    p.TXRayAttenuation[0] = -1E3;
  }
// The exit condition can be safely removed if necessary (Multi-RX setups). 

  // 0 = refracted or 1 = reflected. 
  if(IgnoreCh[1] != 0) {
    p.TXRayTime[1]        = RayTime[1]*s;
    p.TXRayDistance[1]    = RayPath[1]*m;
    p.TXRayStartAngle[1]  = LaunchAngle[1]*deg;
    p.TXRayEndAngle[1]    = RecieveAngle[1]*deg;
    p.TXRayAttenuation[1] = AttRay[1];
  } else {// no secondary solution 
    // (temporary) fix for segments with no solution - we use a the -1E3 value as a flag for 'no solutions found'
    p.TXRayTime[1]        = -1E3;
    p.TXRayDistance[1]    = -1E3;
    p.TXRayStartAngle[1]  = -1E3;
    p.TXRayEndAngle[1]    = -1E3;
    p.TXRayAttenuation[1] = -1E3;
  }
    
// Second IRT call, for CS -> RX
// Reset the default values
  std::fill(IgnoreCh, IgnoreCh+2, 0);
  std::fill(RayTime, RayTime+2, -1E3);
  std::fill(RayPath, RayPath+2, -1E3);
  std::fill(LaunchAngle, LaunchAngle+2, -1E3);
  std::fill(RecieveAngle, RecieveAngle+2, -1E3);
  std::fill(IncidenceAngleInIce, IncidenceAngleInIce+2, -1E3);
  std::fill(AttRay, AttRay+2, -1E3);
  
  x1 = distance(p.Position[0], p.Position[1], fRX.Pos()[0], fRX.Pos()[1])/m;
  z0 = p.Position[2]/m;
  z1 = fRX.Pos()[2]/m;
  
  //  cout<<" Running IRT CS -> RX:"<<endl;
  IceRayTracing::GetRayTracingSolutions(z1, x1, z0, RayTime, RayPath, 
          LaunchAngle, RecieveAngle, IgnoreCh, IncidenceAngleInIce, A0, frequency, AttRay);

  // Nan check, in case IRT runs into problems. Can be disabled, but may introduce errors into final output.
  for(int i=0;i<2;i++){
    if((std::isnan(RayTime[i])==true) || (std::isnan(RayPath[i])==true) || (std::isnan(LaunchAngle[i])==true) ||
        (std::isnan(RecieveAngle[i])==true) || (std::isnan(AttRay[i])==true)){
      cout<<"Nan found in CS -> RX IRT solutions, for Ch "<<i<<endl;
      cout<<"IRTinput [m], z0: "<<z0<<"   x1: "<<x1<<"  z1: "<<z1<<endl;
      cout<<"IRToutput, time: "<<RayTime[i]*s<<"  path: "<<RayPath[i]<<"   Langle:"<<LaunchAngle[i]<<"   Rangle: "<<RecieveAngle[i]<<"   Attenuation: "<<AttRay[i]<<endl;
      exit(EXIT_FAILURE);
    }
  }

  if(IgnoreCh[0] != 0){
    p.RXRayTime[0]        = RayTime[0]*s;
    p.RXRayDistance[0]    = RayPath[0]*m;
    p.RXRayStartAngle[0]  = LaunchAngle[0]*deg;
    p.RXRayEndAngle[0]    = RecieveAngle[0]*deg;
    p.RXRayAttenuation[0] = AttRay[0];
  } else {
    // No direct solution found (no solutions at all.) In this case, we want to 'turn off' the cascade segment. 
    // std::cerr << "IceRayTracing could not find a ray for RX" << std::endl;
    // exit(EXIT_FAILURE); 

    // (temporary) fix for segments with no solution - we use a the -1E3 value as a flag for 'no solutions found'
    // Note - we no longer stop the simulation if this occurs. 
    p.RXRayTime[0]        = -1E3;
    p.RXRayDistance[0]    = -1E3;
    p.RXRayStartAngle[0]  = -1E3;
    p.RXRayEndAngle[0]    = -1E3;
    p.RXRayAttenuation[0] = -1E3;
  }
// The exit condition can be safely removed if necessary (Multi-RX setups). 

  if(IgnoreCh[1] != 0) { // fill in the secondary solution (reflected/refracted)
    p.RXRayTime[1]        = RayTime[1]*s;
    p.RXRayDistance[1]    = RayPath[1]*m;
    p.RXRayStartAngle[1]  = LaunchAngle[1]*deg;
    p.RXRayEndAngle[1]    = RecieveAngle[1]*deg;
    p.RXRayAttenuation[1] = AttRay[1];
  } else {  // no secondary solution 
    // (temporary) fix for segments with no solution - we use a the -1E3 value as a flag for 'no solutions found'
    p.RXRayTime[1]        = -1E3; 
    p.RXRayDistance[1]    = -1E3;
    p.RXRayStartAngle[1]  = -1E3;
    p.RXRayEndAngle[1]    = -1E3;
    p.RXRayAttenuation[1] = -1E3;
  }

  if (p.TXRayTime[0] < 0.0 || p.RXRayTime[0] < 0.0){ // turn off segment if no IRT solution found (Ray times are -ve, from -1E3 flag)
    // In this case, we keep the ArrivalTime of the segment as -1E3, to signal to the voltage calculation in scatter.cc to ignore that segment, and set all other variables to 0.0
    p.ArrivalTime = -1E3;
    p.RTX         = 0.0;
    p.RRX         = 0.0;
    p.Phase       = 0.0;
    p.Attenuation = 0.0;
    p.PolEff      = 0.0;

  } else{ // If not the case, we can go ahead and run the scatter propagation calculations.
  
  // Choose whether to use the first (=0) or second (=1) solution:
  int DoR = 0; // 0 may be either a direct or refracted solution, and 1 is normally the reflected solution.
  

  p.TXDir = direction( fTX.Pos(), p.Position) ;   // Straight line direction vector between TX and CS
  p.RTX = p.TXRayDistance[DoR];                   // IRT distance of the TX -> CS ray path 

  p.RXDir = direction(p.Position, fRX.Pos());     // Straight line direction vector between CS and RX
  p.RRX = p.RXRayDistance[DoR];                   // IRT distance of the CS -> RX ray path 


  // TODO TRIPLE CHECK THAT THIS IS STILL VALID
  // Set phase
  // According to uzair,
  // Phase = k_vac*Optical path length = k_vac*propagation time*c_vac
  // K_vac*RayPath
  // Waiting for confirmation from Krijn 
  
  // p.Phase = fTX.Wavenumber()*(p.RTX + p.RRX)+ atan2( f_coll, - fTX.Freq() );

  p.Phase = fTX.Wavenumber()*(p.RTX + p.RRX) - pi/2;  // Set the point phase
  // tan^-1(f_coll/f_TX) is the phase shift in the oscillation caused by collisions.
  // in our regime, it could be almost fixed to -Pi/2.

  // DOUBLE CHECK MINUS SIGN
// DOES IRT INCLUDE ICE ATTENUATION?
  // E field attenuation at reciever from position dependant factors:
  // e^-r/Latt from medium attenuation
  // p.Attenuation =  pow(e, -(p.RTX + p.RRX)/(2*att_length) ) ;

  p.Attenuation = p.TXRayAttenuation[DoR]*p.RXRayAttenuation[DoR];  // Point attenuation calculated from IRT output
  p.ArrivalTime = p.StartTime + p.RXRayTime[DoR];                   // Point arrival time, with starttime calculated as p.L/c_vac

  ///////  Ray bending implementation: 
  // This is effectively four rotations applied to the ray as it travels from TX -> CS -> RX.
  // It is required in order for the polarisation factors to be correctly calculated with the IRT launch/receive angles.

  // Find the rotation axis needed for the launch angle direction vector (cross product between z axis and the straight line TXDir)
  p.TXRotAxis   = cross_product( p.zDir , normalize(p.TXDir) ) ;

  // Find the launch direction by rotating the zdir vector by the launch angle, in the zdir/txdir plane. Produces clockwise rotation.
  p.TXDir_l     = rotate( p.zDir, normalize( p.TXRotAxis ), p.TXRayStartAngle[DoR] );

  // Polarisation at the TX; cross product between the launch direction and the TX polarisation (=sin(theta_TX))
  p.PolarizationAtTX = cross_product(normalize(p.TXDir_l), cross_product(normalize(p.TXDir_l), fTX.Pol() ) );
  
  // Find the rotation angle over the TX -> CS path (difference between the launch and receive angles)
  double RotAngleTX = p.TXRayEndAngle[DoR] - p.TXRayStartAngle[DoR] ;
  if(RotAngleTX < 0.){  //Check whether the receive angle is smaller than the launch angle
    p.TXRotAxis = -1. * p.TXRotAxis; // If true, the rotation direction should be flipped
  } // Might have to catch whether the two angles are equal (no rotation is needed.)

  // Get the CS polarisation by rotating the TX polarisation as needed.
  p.Polarization  = rotate( p.PolarizationAtTX, normalize(p.TXRotAxis), abs(RotAngleTX)  );

  // Repeat the process for the CS -> RX ray path
  p.CSRotAxis     = cross_product( p.zDir , normalize(p.RXDir)) ;
  p.CSDir_l       = rotate( p.zDir, normalize(p.CSRotAxis), p.RXRayStartAngle[DoR] );

  // Cross product between the point polarisation and the RX launch angle ray direction (=sin(theta_TX)*sin(theta_CS))
  p.EFieldatCS    = cross_product(normalize(p.CSDir_l), cross_product(normalize(p.CSDir_l), p.Polarization) );
  double RotAngleCS = p.RXRayEndAngle[DoR] - p.RXRayStartAngle[DoR] ;
  if(RotAngleCS < 0.){ // If true, the rotation direction should be flipped
    p.CSRotAxis = -1. * p.CSRotAxis ; }
  p.EFieldAtRX = rotate( p.EFieldatCS, normalize(p.CSRotAxis), abs(RotAngleCS)  );

  // Get the final RX polarisation efficiency! (=cos(theta_RX)*sin(theta_TX)*sin(theta_CS))
  p.PolEff    = abs(projection_pol(p.EFieldAtRX, fRX.Pol()));

  // Directivity is hardcoded as a small Herztian dipole.
  p.Directivity = 1.5*norm(p.EFieldAtRX);
  //  3/2*sin(theta_R)*sin(theta_T)

  }
  //////// IRT Debugging statements, uncomment as needed (will be removed in future, but useful for now).
  // // TX -> CS Ray path:
  // cout<<"TX Launch Angle: "<<p.TXRayStartAngle[DoR]/deg<<endl;
  // cout<<"TX Recieve Angle: "<<p.TXRayEndAngle[DoR]/deg<<endl;
  // cout<<"TX straight line angle: "<< std::acos(projection(normalize(p.TXDir),p.zDir))/deg <<endl;
  // cout<<"TXRotAxis: "<<normalize(p.TXRotAxis);  //debug cout
  // cout<<"TXDir_l: "<<p.TXDir_l;  //debug cout
  // cout<<"TXDir straight line: "<<normalize(p.TXDir);
  // cout<<"RotAngleTX: "<<RotAngleTX/deg << endl;
  // cout<<"PolarizationAtTX: "<<p.PolarizationAtTX;
  // cout<<"PolarizationCS: "<<p.Polarization;
  // cout<<"Polarization at TX length check: "<<norm(p.PolarizationAtTX)<<endl;
  // cout<<"Polarization at CS length check: "<<norm(p.Polarization)<<endl;
  // //
  // // CS -> RX Ray path:
  // cout<<"RX Launch Angle: "<<p.RXRayStartAngle[DoR]/deg<<endl;
  // cout<<"RX Recieve Angle: "<<p.RXRayEndAngle[DoR]/deg<<endl;
  // cout<<"RX straight line angle: "<< std::acos(projection(p.zDir,normalize(p.RXDir)))/deg <<endl;
  // cout<<"CSRotAxis: "<<normalize(p.CSRotAxis);  //debug cout
  // cout<<"CSDir_l: "<<p.CSDir_l;  //debug cout
  // cout<<"RXDir straight line: "<<normalize(p.RXDir); //debug cout
  // cout<<"RotAngleCS: "<< RotAngleCS/deg << endl;
  // cout<<"EFieldAtCS: "<<p.EFieldatCS;
  // cout<<"EFieldAtRX: "<<p.EFieldAtRX;
  // cout<<"EField at CS length check: "<<norm(p.EFieldatCS)<<endl;
  // cout<<"EField at RX length check: "<<norm(p.EFieldAtRX)<<endl;
  // cout<<"PolEff: "<<p.PolEff<<endl;
  // //
  // // Final output values:
  // cout<<"IRT attenuation: "<<p.Attenuation<<endl;
  // cout<<"IRT RTX: "<<p.RTX<<endl;
  // cout<<"IRT RRX: "<<p.RRX<<endl;
  // cout<<"IRT arrival time: "<<p.ArrivalTime<<endl;
  // cout<<"IRT phase: "<<p.Phase<<endl;
  
  // Final nan checks (should not be triggered, but just in case):
  if(std::isnan(p.PolEff)==true) { 
    cout<<"PolEff nan flag, at depth "<<p.Position[2]/m<<endl;}
  if(std::isnan(p.RRX)==true) {
    cout<<"RRX nan flag, at depth  "<<p.Position[2]/m<<endl;}
  if(std::isnan(p.RTX)==true) {
    cout<<"RTX nan flag, at depth  "<<p.Position[2]/m<<endl;}
  if(std::isnan(p.ArrivalTime)==true) {
      cout<<"Arrival time nan flag, at depth  "<<p.Position[2]/m<<endl;}
  //////////////////
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