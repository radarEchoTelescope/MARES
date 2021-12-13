/* Radar computation - return electric field of a list of events*/
#include "macro_scatter.hh"
int main(int argc, char** argv){

// For now, we stick with the variables being defined either as write-in or
// fed through the command line.

  // I/O Setup

  // As default, the result is stored in the same folder as the executable.
  std::string path_out = "";
  // std::string path_out = cs_filepath.substr(0,cs_filepath.find_last_of("."));

  std::string identifier  = argv[1];

  /* Make the cascade OR parse through load_cascade_file() */
  int evt = 0;
  double cenergy = 1E10;
  // double xpos = 0;
  double xpos = -(0.336854*5)/2.0;
  double ypos = 0 ;
  double zpos = 0;
  double czenith = 90;
  double cazimuth = 0;
  czenith = deg2rad(czenith);
  cazimuth = deg2rad(cazimuth);

  double nenergy = 1E9;
  double nzenith = 0;
  double nazimuth = 0;
  double oneweight = 1;


  Cascade cascade(evt, cenergy, xpos, ypos, zpos, czenith, cazimuth,
                  nenergy, nzenith, nazimuth, oneweight);

  Cascade& cs = cascade;
  // cout << cs.Sph()[0] << '\t' << cs.Sph()[1] << endl;

  // Make your detector
  Detector lab;

  double R = 200; // Antenna - Wire distance [m]

  std::vector<double> TX_angles = {
                                    // 50
                                    10,30,50,70,90
                                    // -10,-30,-50,-70,-90
                                  };

  for(auto& i : TX_angles){
    double ii = deg2rad( i );
    // cout << R*cos(ii) << '\t' << R*sin(ii) << endl;
    lab.add_antenna( Antenna(1E4,  R*cos(ii), R*sin(ii),0 ) );
  }

  std::vector<double> RX_angles = arange<double>(0, 360, 1);
  // std::vector<double> RX_angles = {
  //   // 90
  //                                 // 90, -90
  //                                 // -10,-20,-30,-40,-50,-60,-70,-80,-90,-100,-110,-120,-130,-140,-150,-160,-170,
  //                                 // -180, -170, -160, -150, -140, -130, -120, -110, -100, -90, -80, -70, -60, -50, -40, -30, -20, -10,
  //                                 0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,
  //                                 190, 200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360
  //                                 };
  // for(auto& j : RX_angles){
  //   // cout << j << '\t'<< R*cos( jj ) << '\t' << R*sin( jj ) << endl;
  //   lab.add_antenna( Antenna(0, j, R , 0 ) );
  // }


  // cout << '\n' << endl;
  for(auto& j : RX_angles){
    double jj = deg2rad(j);
    // cout << j << '\t'<< R*cos( jj ) << '\t' << R*sin( jj ) << endl;
    lab.add_antenna( Antenna(0,  R*cos(jj) , R*sin(jj), 0 ) );

  }

  int i = 0;
  int j;
  // For every transmitter-receiver pair
  for (auto& tx : lab.Transmitters()){
    j = 0;
    for (auto& rx : lab.Receivers()){
      Line1D thinwire(tx,rx,cs);

      std::string identifier_l = path_out
                                  + "TX_" + std::to_string((int) TX_angles[i]) + "_deg"
                                  + "RX_" + std::to_string((int) RX_angles[j]) + "_deg";

      /* Choose what to write out by uncommenting the lines. */
      write_2D_array(thinwire.Coordinates(),      identifier_l + "_coords.txt", 1);
      // write_1D_array(thinwire.Attenuation(),   identifier_l + "_attenuation.txt", 1);
      // // write_1D_array(thinwire.arrivals(),   identifier_l + "_t_arrivals.txt", 1);
      write_1D_array(thinwire.TCS(),              identifier_l + "_TCS.txt", 1);
      write_1D_array(thinwire.Duration(),         identifier_l + "_duration.txt", 1);
      write_1D_array(thinwire.Directivity(),      identifier_l + "_directiviy.txt", 1);
      write_1D_array(thinwire.Waveform(),         identifier_l + "_waveform.txt", 1);
      write_1D_array(thinwire.Power(),            identifier_l + "_power.txt", 1);
      write_1D_array(thinwire.RCS(),              identifier_l + "_RCS.txt", 1);
      // write_2D_array(thinwire.phase_time(),    identifier_l + "_phase_time_profile.txt", 1);
      // write_2D_array(thinwire.wave_time(),     identifier_l + "_Er_time_profile.txt", 1);

      j++;
    }
  i++;
  }
}
// End
