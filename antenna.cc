#include "antenna.hh"
// using namespace std;

Antenna::Antenna(double power, double position, double polarization, double G,
  double frequency):
  power(power), pos{position}, polar{polarization}, f_obs(frequency){

  k_obs = w_obs/c_med;
  l_obs = c_med/f_obs;
  w_obs = 2*pi *f_obs;

  if (!G){ power ? gain = 1 : gain = pow(l_obs,2.0); }
  /* If G = 0, default , unphysical value, then set default gains
  Default TX gain (power nonzero) =  1, Isotropic emission.
  Defalut RX gain (power zero): Effective area. */
}

Antenna::Antenna(double power, double position) :
Antenna::Antenna(power, position, (0,0), 0, 1E9){} //default freq_obs

/* Set the antennas directions, module and dot product with cs */
void Antenna::set_direction(const Cascade & cs){

  // Direction to cascade
  cs_dir[0] = (cs.pos[0] - pos[0]);
  cs_dir[1] = (cs.pos[1] - pos[1]);
  cs_dir[2] = (cs.pos[2] - pos[2]);

  // Correct orientation for the Rx case (cs_to_at).
  if (power == 0) {   // Rx
    cs_dir[0] = - cs_dir[0];
    cs_dir[1] = - cs_dir[1];
    cs_dir[2] = - cs_dir[2];
  }

  // Module of distance to cascade
  cs_dist= sqrt(pow(cs_dir[0],2.0)+pow(cs_dir[1],2.0)+pow(cs_dir[2],2.0));

  // Sanity check
  assert(cs_dist!= 0 && "Cascade overlaps antenna");
  /* THIS SHOULD THROW AND EXCEPTION SO YOU CAN CHOOSE HOW TO SOLVE IT
      FOR NOW:
 */
 if(cs_dist== 0) {cs_dist = 10*l_obs;}

  //Determine inner product between point_tc and cascade direction
  cs_dot =  (cs_dir[0]*cs.dir[0] + cs_dir[1]*cs.dir[1] + cs_dir[2]*cs.dir[2])
            / cs_dist; //Normalize
}


/* Set the antennas positions and directions with IceRayTracing corrections */

/*To add, how to use IRT*/

// double * at_to_cs = IceRayTracing::IceRayTracing(x0,z0,x1,z1);

// Get the 2 rays between tx and cascade

// Get the proper distances from the solutions.
// Module of distance to cascade is now the optical path length.
// The cs_angle is the IRT recieved angle for the at_to_cs and the emitted
// angle for cs_to_rx case.

void Antenna::IRT_set_direction(const Cascade & cs){
  double at_r,at_z,cs_r,cs_z;
  double * at_to_cs;

  at_r = sqrt(pow(pos[0],2) + pow(pos[1],2));
  at_z = pos[2];
  cs_r = sqrt(pow(cs.pos[0],2)+pow(cs.pos[1],2));
  cs_z = cs.pos[2];

  if(power) { // TX case
    // The antenna emits, the cascade is the one that recieves the signal.
    at_to_cs = IceRayTracing::IceRayTracing(at_r,at_z,cs_r,cs_z);

  } else { // RX case
    // The "emitter" now is the cascade and the "antenna" the receiver.
    at_to_cs = IceRayTracing::IceRayTracing(cs_r,cs_z,at_r,at_z);

  }

  // Now, let's unpack the values that we recieve from IRT.

  if(at_to_cs[6] == 0){               // No direct ray
    if(at_to_cs[4] < at_to_cs[5]){    // Reflected is shortest
      IRT_cs_dist[0] = at_to_cs[4]*c_med;        // [m]
      IRT_cs_dist[1] = at_to_cs[5]*c_med;   // [m]
      IRT_cs_angle[0] = at_to_cs[7];      // degrees
      IRT_cs_angle[1] = at_to_cs[8];      // degrees

    } else {
      IRT_cs_dist[0] = at_to_cs[5]*c_med;        // [m]
      IRT_cs_dist[1] = at_to_cs[4]*c_med;
      IRT_cs_angle[0] = at_to_cs[8];
      IRT_cs_angle[1] = at_to_cs[7];      // degrees

    }
  }

  if(at_to_cs[7] == 0){               // No reflected ray
    if(at_to_cs[3] < at_to_cs[5]){
      IRT_cs_dist[0] = at_to_cs[3]*c_med;        // [m]
      IRT_cs_dist[1] = at_to_cs[5]*c_med;
      IRT_cs_angle[0] = at_to_cs[6];      // degrees
      IRT_cs_angle[1] = at_to_cs[8];      // degrees

    } else {
      IRT_cs_dist[0] = at_to_cs[5]*c_med;        // [m]
      IRT_cs_dist[1] = at_to_cs[3]*c_med;
      IRT_cs_angle[0] = at_to_cs[8];
      IRT_cs_angle[1] = at_to_cs[6];      // degrees

    }
  }

  if(at_to_cs[8] == 0){               // No refracted ray
    if(at_to_cs[3] < at_to_cs[4]){
      IRT_cs_dist[0] = at_to_cs[3]*c_med;        // [m]
      IRT_cs_dist[1] = at_to_cs[4]*c_med;
      IRT_cs_angle[0] = at_to_cs[6];      // degrees
      IRT_cs_angle[1] = at_to_cs[7];      // degrees

    } else {
      IRT_cs_dist[0] = at_to_cs[4]*c_med;        // [m]
      IRT_cs_dist[1] = at_to_cs[3]*c_med;
      IRT_cs_angle[0] = at_to_cs[7];      // degrees
      IRT_cs_angle[1] = at_to_cs[6];      // degrees

    }
  }

  // Sanity check
  // In theory, only the first value (shortest distance) needs to be evaluated.
  assert(IRT_cs_dist[0] != 0 && "Cascade overlaps antenna");
  /* THIS SHOULD THROW AND EXCEPTION SO YOU CAN CHOOSE HOW TO SOLVE IT
      FOR NOW:
 */

 // But both need to be fixed.
 if(IRT_cs_dist[0] == 0) {IRT_cs_dist[0] = 10*l_obs;}
 if(IRT_cs_dist[1] == 0) {IRT_cs_dist[1] = 10*l_obs;}

}

// Accesors

double Antenna::get_power(){return power;}
double Antenna::get_gain(){return gain;}
double* Antenna::get_position(){return pos;}
double* Antenna::get_polarization(){return polar;}
double Antenna::get_f_obs(){return f_obs;}
double Antenna::get_l_obs(){return l_obs;}
double Antenna::get_w_obs(){return w_obs;}
double Antenna::get_k_obs(){return k_obs;}

double* Antenna::get_cs_direction(){return cs_dir;}
double Antenna::get_cs_distance(){return cs_dist;}
double Antenna::get_cs_projection(){return cs_dot;}

double* Antenna::get_IRT_distances(){return IRT_cs_dist;}
double* Antenna::get_IRT_angles(){return IRT_cs_angle;}

// --------------------------------------------

Detector::Detector(){

  Transmitters.push_back( Antenna(1, (0,0,0)) );
  Receivers.push_back( Antenna(0, (500,0,0)) );
};

//
// Detector::Detector(string det_name){
//   if (det_name == "RET_CR"){
//     Antenna* tx(1, 1E9, (0,0,0));
//     Transmitters.push_back(*tx);
//
//     Antenna rx0(double power = 0, double frequency = 1E9,
//       double position = (200,0,0), double polarization = (0,0));
//     Receivers.push_back(rx0);
//
//     Antenna rx1(double power = 0, double frequency = 1E9,
//         double position = (100,-100,0), double polarization = (0,0));
//     Receivers.push_back(rx1);
//
//     Antenna rx2(double power = 0, double frequency = 1E9,
//         double position = -100,-100,0), double polarization = (0,0));
//     Receivers.push_back(rx2);
//
//     Antenna rx3(double power = 0, double frequency = 1E9,
//         double position = (-100,100,0), double polarization = (0,0));
//     Receivers.push_back(rx3);
// k_obs
//     Antenna rx4(double power = 0, double frequency = 1E9,
//         double position = (200,200,0), double polarization = (0,0));
//     Receivers.push_back(rx4);
//
//     Antenna rx5(double power = 0, double frequency = 1E9,
//         double position = (200,-200,0), double polarization = (0,0));
//     Receivers.push_back(rx5);
//
//     Antenna rx6(double power = 0, double frequency = 1E9,
//         double position = (-200,-200,0), double polarization = (0,0));
//     Receivers.push_back(rx6);
//
//     Antenna rx7(double power = 0, double frequency = 1E9,
//         double position = (-200,200,0), double polarization = (0,0));
//     Receivers.push_back(rx7);
  //
  // } else {
  //   cout << " There is no compatible detector config with " + det_name << endl;
  // }
// };

std::vector<Antenna> Detector::get_transmitters(){return Transmitters;}
std::vector<Antenna> Detector::get_receivers(){return Receivers;}
