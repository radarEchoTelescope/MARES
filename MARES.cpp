/* Radar computation - return electric field of a list of events*/

#include "cascade1D.hh"
using namespace libconfig;

int main(int argc, char** argv){

// For now, we stick with the variables being defined either as write-in or
// fed through the command line.

  // I/O Setup

  // As default, the result is stored in the same folder as the executable.
  // std::string path_out = "";
  // std::string path_out = cs_filepath.substr(0,cs_filepath.find_last_of("."));

  const char* config_file;
  Config cfg;
  // If no argument with a config file has been passed, expect a "MARES.cfg file". 
  argc < 2 ? config_file = "MARES_example.cfg": config_file = argv[1];

  // Read the file. If there is an error, report it and exit.
  try{
    cfg.readFile(config_file);
  } catch(const FileIOException &fioex) {
      std::cerr << "I/O error while reading file." << std::endl;
      return(EXIT_FAILURE);
  } catch(const ParseException &pex) {
      std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
      return(EXIT_FAILURE);
  }

  std::string path_out = "./";
  if(cfg.lookupValue("path_out", path_out)){
    std::cout << "The output will be saved in: " << path_out << std::endl;
  } else {
    std::cout << "The output will be saved locally." << std::endl;
  }

  std::string identifier = "test";
  if(cfg.lookupValue("name", identifier)){
    std::cout << "Starting run: " << identifier << std::endl;
  }

  Setting &root = cfg.getRoot();

  const Setting& cs_list = root["cascade"];
  const Setting& tx_list = root["detector"]["transmitter"];
  const Setting& rx_list = root["detector"]["receiver"];
  // try{
  // } catch(const SettingNotFoundException &nfex) {
  //     std::cerr << "Your config file is missing the basic elements" << std::endl;
  //     return(EXIT_FAILURE);
  // }

  // Make your detector
  Detector lab;
  // Antenna placeholder variables.
  double xpos, ypos, zpos, xpol, ypol, zpol, power, freq, gaindB;

  // Add the transmitters
  for (int i = 0; i < tx_list.getLength(); i++){
    // Grab the next antenna
    const Setting &tx = tx_list[i];

    // Try to load the values from the antenna
    if( !(
          tx["position"].lookupValue("x", xpos)   &&
          tx["position"].lookupValue("y", ypos)   &&
          tx["position"].lookupValue("z", zpos)   &&
          tx["polarization"].lookupValue("x", xpol)      &&
          tx["polarization"].lookupValue("y", ypol)      &&
          tx["polarization"].lookupValue("z", zpol)      &&
          tx.lookupValue("power", power)               &&
          tx.lookupValue("frequency", freq)           &&
          tx.lookupValue("gaindB", gaindB)
        )
      ){
      std::cerr << "Transmitter " << i << " was skipped" << std::endl;
      // Let's skip this antenna
      continue;
    } 

    // You place here the sanity checks on your antenna

    // Order matters here!
    // We add the units here!  
    lab.add_antenna( 
      Antenna(
        xpos *m, ypos *m, zpos *m,
        xpol, ypol, zpol,
        power *W, freq *Hz, gaindB
      )
    );
  }

  for (int i = 0; i < rx_list.getLength(); i++){
    // Grab the next antenna
    const Setting &rx = rx_list[i];

    // Try to load the values from the antenna
    if( !(rx["position"].lookupValue("x", xpos)   &&
          rx["position"].lookupValue("y", ypos)  &&
          rx["position"].lookupValue("z", zpos)  &&
          rx["polarization"].lookupValue("x", xpol)      &&
          rx["polarization"].lookupValue("y", ypol)      &&
          rx["polarization"].lookupValue("z", zpol)      &&
          rx.lookupValue("power", power)               &&
          rx.lookupValue("frequency", freq)           &&
          rx.lookupValue("gaindB", gaindB)
        )
      ){
      std::cerr << "Receiver " << i << " was skipped" << std::endl;
      // Let's skip this antenna
      continue;
    } 

    // You place here the sanity checks on your antenna

    // Extra check for receivers
    assert(power == 0.0 && "This antenna is not defined as a receiver" );

    // Order matters here!
    // We add the units here!   
    lab.add_antenna( 
      Antenna(
        xpos *m, ypos *m, zpos *m,
        xpol, ypol, zpol,
        power *W, freq *Hz, gaindB
      )
    );
  }


  std::vector<Cascade> cascade_list;
  /* Make the cascade OR parse through load_cascade_file() */
  int evt = 0, cprimaries;
  double cenergy, czenith, cazimuth;

  for (int i = 0; i < cs_list.getLength(); i++){
    // Grab the next cascade
    const Setting &cs = cs_list[i];
    
    if( !(
          cs.lookupValue("nu_id", evt)      &&
          cs.lookupValue("energy", cenergy)               &&
          cs["position"].lookupValue("x", xpos)   &&
          cs["position"].lookupValue("y", ypos)  &&
          cs["position"].lookupValue("z", zpos)  &&
          cs["direction"].lookupValue("zenith", czenith)      &&
          cs["direction"].lookupValue("azimuth", cazimuth)      &&
          cs.lookupValue("primaries", cprimaries)          
        )
    ){
      // print warning or error!
      std::cerr << "Cascade " << i << " was skipped" << std::endl;
      // Let's skip this antenna, or save it in a text file?
      continue;
    } 
    std::cout << "Cascade " << i << " was parsed correctly" << std::endl;
    // You place here the sanity checks on your cascade
    // For example, z < 0 if you want to work with an IRT-like frame.
    
    // Order matters here!
    // We add the units here!   
    cascade_list.push_back( 
      Cascade( 
        evt, cenergy * GeV,
        xpos *m, ypos *m, zpos *m,
        czenith *deg, cazimuth *deg,
        cprimaries
      )
    );
  }

  // We check for the saving flags, aka what properties are going to be stored after the event is run.

  std::array<bool,8>  cascade_flags{0};
  const Setting& cascade_flag_list = root["save"]["cascade"];
  cascade_flag_list.lookupValue("electron_number",    cascade_flags[0]);
  cascade_flag_list.lookupValue("electron_density",   cascade_flags[1]);
  cascade_flag_list.lookupValue("plasma_frequency",   cascade_flags[2]);
  cascade_flag_list.lookupValue("plasma_absorption",  cascade_flags[3]);
  cascade_flag_list.lookupValue("skin_depth",         cascade_flags[4]);
  cascade_flag_list.lookupValue("reflectance",        cascade_flags[5]);
  cascade_flag_list.lookupValue("opacity",            cascade_flags[6]);
  cascade_flag_list.lookupValue("transparency",       cascade_flags[7]);
  int total_cascade_flags = std::count(cascade_flags.begin(), cascade_flags.end(), true);

  std::array<bool,13> scatter_flags{0};
  const Setting& scatter_flag_list = root["save"]["scatter"];
  scatter_flag_list.lookupValue("duration",             scatter_flags[0]);
  scatter_flag_list.lookupValue("voltage",              scatter_flags[1]);
  scatter_flag_list.lookupValue("power",                scatter_flags[2]);
  scatter_flag_list.lookupValue("target_cs",            scatter_flags[3]);
  scatter_flag_list.lookupValue("radar_cs",             scatter_flags[4]);
  scatter_flag_list.lookupValue("point_position",       scatter_flags[5]);
  scatter_flag_list.lookupValue("points_phase",         scatter_flags[6]);
  scatter_flag_list.lookupValue("points_arrival_time",  scatter_flags[7]);
  scatter_flag_list.lookupValue("points_attenuation",   scatter_flags[8]);
  scatter_flag_list.lookupValue("points_polatization",  scatter_flags[9]);
  scatter_flag_list.lookupValue("phase_vs_time",        scatter_flags[10]);
  scatter_flag_list.lookupValue("radar_vs_time",        scatter_flags[11]);
  scatter_flag_list.lookupValue("voltage_vs_time",      scatter_flags[12]);
  int total_scatter_flags = std::count(scatter_flags.begin(), scatter_flags.end(), true);

  const bool save_time_profiles = {scatter_flags[10] && scatter_flags[11] && scatter_flags[12]};

  if( !total_cascade_flags && !total_scatter_flags){
    std::cerr << "There will be no output from this event! Please enable at least a saving flag." << std::endl;
    exit(1);
  }
  if( !total_cascade_flags ){ std::cout << "There will be no cascade properties being saved" << std::endl; }
  if( !total_scatter_flags){ std::cout << "There will be no scatter properties being saved" << std::endl; }

  // std::cout << scatter_flags << std::endl;
  // std::cout << cascade_flags << std::endl;
// We prepare the event

  int tx_number = lab.Transmitters().size();
  int rx_number = lab.Receivers().size();
  int cs_number = cascade_list.size();

  Cascade1D* event;
  std::string unique_id;

  for (int k = 0; k < cs_number; k++){
    Cascade& cs = cascade_list[k];
    for (int i = 0; i < tx_number; i++){
      auto tx = lab.Transmitters()[i];
      for (int j = 0; j < rx_number; j++){
        auto rx = lab.Receivers()[j];
        
        if(j == 0){
          event = new Cascade1D(tx,rx,cs);

          // If flag
          event -> SetInDirection();

        } else {
          // Update RX,
          // make sure it updates the direction to the antennas in the appropiate function
          event -> UpdateRX(rx);
        }
        // If flag 
        event -> SetInConstIce();
        event -> RunScatter(save_time_profiles);

        unique_id = path_out + "Cascade1D_" + identifier +
                  "_TX_" + std::to_string((int) i) +
                  "_RX_" + std::to_string((int) j) + 
                  "_CS_" + std::to_string((int) cs.Evt() );

        event -> Cascade::save_output_files(unique_id, cascade_flags);
        event -> Scatter::save_output_files(unique_id, scatter_flags);

      }
    }
  }
  // Housekeeping!
  // delete event;

}
// End

  /* arange() works as in python: Start,stop, steps
    For a single value, start = stop and step > 0.
    E.g: 90 = arange<double>(90, 90, 1)
  */