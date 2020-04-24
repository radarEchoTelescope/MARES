// Calculates the return power for a bi-static radat setup
// Branched of power.C at 23/1/2020
// Enrique Huesca Santiago, 10-2019
// Original FORTRAN code by Krijn D. de Vries 20-10-2019

#ifndef MACRO_SCATTER
#define MACRO_SCATTER

#include "cascade.hh"
#include "antenna.hh"
#include "macro_settings.hh"

// #define NDEBUG     // Turn off debug.

//---- Free (user) parameters---------------------------------------------------
const string identifier = "return_power_0D";
const bool write_out_arrays = 1;

//---- I/O functions -----------------------------------------------------------
void write_1D_array(std::vector<double> array, string output_path, const bool trigger){
  if(trigger){
    ofstream output_file(output_path);
    ostream_iterator<double> output_iterator(output_file, "\t");
    copy(array.begin(), array.end(), output_iterator);
    output_file.close();
  }
}

void write_2D_array(std::vector<std::vector<double>> array, string output_path, const bool trigger){
  if(trigger){
    ofstream output_file(output_path);
    ostream_iterator<double> output_iterator(output_file, "\t");
    for (int i = 0; i < array.size(); i++){
      copy(array.at(i).begin(), array.at(i).end(), output_iterator);
      output_file << '\n';
    }
    output_file.close();
  }
}

std::vector<Cascade> load_cascade_file(const string& cs_filepath){

  ifstream cs_file(cs_filepath);

  if (!cs_file.is_open()){
    cout << "The input file is not opening" << endl;
    //outfile << "The input file is not opening" << endl;
    //exit (EXIT_FAILURE);
    throw 1;
  }

  int eventnr;
  double nzenith,nazimuth,nenergy,czenith,cazimuth,cenergy,xpos,ypos,zpos,oneweight;

  std::vector<Cascade> Cascades;

  std::vector<double> rejected;
  std::vector<std::vector<double>> rejected_particles;

  while (cs_file  >> eventnr >> nzenith >> nazimuth >> nenergy >> czenith >> cazimuth >>
  cenergy >> xpos >> ypos >> zpos >> oneweight) {
  // Read line from file

    if (cenergy < 1E6) {              // Energy check
       cout << "Skipping cascade " << eventnr <<
       " since it has energy below the 1 PeV threshold" << endl;
       rejected = {eventnr, nzenith, nazimuth, nenergy, czenith, cazimuth,
       cenergy, xpos, ypos, zpos, oneweight};
    }

/*  This is omitted for testing purposes
     else if((zpos - 1389) > 0){   // Position check
      // 1389 = 2778/2, is the middle of the ice shelf [m]
      cout << "The cascade "<< eventnr <<
        " has a positive z-position (out-of-ice)" << endl;
      rejected = {eventnr, nzenith, nazimuth, energy, czenith, cazimuth,
      cenergy, xpos, ypos, zpos, oneweight};
    }
*/

    if (rejected.empty()){
        Cascade cs(eventnr, (xpos,ypos,zpos),
          cenergy, cazimuth, czenith, nenergy, nzenith, nazimuth, oneweight);
        Cascades.push_back(cs);
    } else {
      rejected_particles.push_back(rejected);
      rejected.clear();
    }

  }

  if (!rejected_particles.empty()){
    string path_out = cs_filepath.substr(0,cs_filepath.find_last_of("."))
      + "_rejected.out";
    write_2D_array(rejected_particles, path_out, 1);
  }
  return Cascades;
}
// -------------------------------------------



#endif

//------------------------------------------------------------------------------
/* NO LONGER IN USE */



/*      cout
      << "txpower = " << tx.power << '\t'
      <<"txgain = " << tx.gain << '\t'
      << "rxarea = " << rx.area << '\t'
      << "eta = " << eta << '\t'
      << "exp = " << exp(-2.0*(tx.R+rx.R)/Latt) << '\t'
      << "pow = " << pow(4*pi*tx.R*rx.R,2) << '\t'
      << endl;
/**/

/*      cout <<
      "Event:" << scientific << '\t' <<
      "R_t = " << tx.R + rx.R << '\t' <<
      //"ud_cs = " << ud_cs << '\t' <<
      "od_tot = " << od_tot << '\t' <<
      "rad_cs = " << rad_cs << '\t' <<
      "P_r = " << P_r << defaultfloat << endl;
/**/

/* 1-D reflectance
    std::vector<double> reflectance1D, reflectivity1D;

    // Loop over the layers, from the outside in, with increasing skin depth.
    for (int k = 0; k <nbin ; k++){
      k_mid = nbin - k + 0.5;
      // This takes care that the last layer is smaller, since we look at midpoints.
      k == nbin ? dr = (0.5)*r_bin : dr = r_bin;              // [cm]

      wplasma=8980*sqrt(dens(X_max,k_mid*r_bin,cs.energy))*sqrt(1/mme);  // [Hz]
      skin=cmed_cm/(2*wplasma);

      reflectance_tmp = (1-reflectivity_tmp)*(1-exp(-1*dr/skin));
      reflectivity_tmp += reflectance_tmp;
      assert(reflectivity_tmp < 1 && "Reflectivity larger than 1!");

      reflectance1D.push_back(reflectance_tmp);
      reflectivity1D.push_back(reflectivity_tmp);

      //reflectance.insert(reflectance.begin(), reflectance_tmp);
      //reflectivity.insert(reflectivity.begin(), reflectivity_tmp);
    }

    write_1D_array(reflectance1D, "1D_reflectance.txt");
    write_1D_array(reflectivity1D, "1D_reflectivity.txt");
/* */

/*
// Parametrized attenuation length for the Ross Ice Shelf, South Pole.
double att(double freq_obs){
double a1=469;                  // [m] Attanuation length parameter
double a2=-0.205;               // Attanuation length parameter
double a3=4.87E-5;              // Attanuation length parameter
return a1+a2*freq_obs/1E6+a3*pow(freq_obs/1E6,2);
}

double intensity(double gam, double L, double lambda){
  return pow(sinc((pi*L/lambda)*sin(gam)),2);
}

double av_intensity(double gam, double L, double lambda ){
  double av_int = 0, ang, bins = 200000;
  ang = -pi/2 +pi*(i+1.0)/bins;
  for (int i = 0; i < bins; i++){
    av_int += intensity(ang, L, lambda);
  }
  av_int /= bins;
  return av_int;
}
*/

/*
// Old FORTRAN defintion
double av_intensity(double gam, double L, double lambda ){
  double av_int = 0, xx=0, var, dvar;
  for (int i = 0; i < 200000; i++){
    xx = -1+2.0*(i+1)/200000;     // Make the diivison a double with the 2.0!
    var = xx * pi/2;
    dvar = pi/200000;
    av_int += intensity(var,L,lambda)*dvar/pi;
  }
  return av_int;
}

*/
