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

  // If no argument with a config file has been passed, expect a "MARES.cfg file". 
  const char* config_file;
  argc < 2 ? config_file = "MARES_example.cfg": config_file = argv[1];

  Config cfg;
  load_config_file(cfg, config_file);

  std::string path_out = "./";
  if(cfg.lookupValue("path_out", path_out)){
    std::cout << "The output will be saved in: " << path_out << std::endl;
  } else {
    std::cout << "The output will be saved locally." << std::endl;
  }

  std::string identifier = "example";
  if(cfg.lookupValue("name", identifier)){
    std::cout << "Starting run: " << identifier << std::endl;
  }

  Setting &root = cfg.getRoot();

  Detector lab = load_detector_config(cfg);
    // std::string path_out = cs_filepath.substr(0,cs_filepath.find_last_of("."))
  //   + "_rejected.out";

  std::vector<Cascade> cascade_list = load_cascade_config(cfg);

// Remeber to check that your cascade list fits your requirements.
// Position, min energy, etc.

  // We check now what is the type of run that we want
  const Setting& options_list = root["options"];

  std::string position = "direction";
  options_list.lookupValue("position",  position);

  std::string propagation = "const";
  options_list.lookupValue("propagation",  propagation);

  // We check for the saving flags, aka what properties are going to be stored after the event is run.

  int total_cascade_flags;
  std::array<bool,8>  cascade_flags{0};
  try{
    const Setting& cascade_flag_list = root["save"]["cascade"];
    cascade_flag_list.lookupValue("electron_number",    cascade_flags[0]);
    cascade_flag_list.lookupValue("electron_density",   cascade_flags[1]);
    cascade_flag_list.lookupValue("plasma_frequency",   cascade_flags[2]);
    cascade_flag_list.lookupValue("plasma_absorption",  cascade_flags[3]);
    cascade_flag_list.lookupValue("skin_depth",         cascade_flags[4]);
    cascade_flag_list.lookupValue("reflectance",        cascade_flags[5]);
    cascade_flag_list.lookupValue("opacity",            cascade_flags[6]);
    cascade_flag_list.lookupValue("transparency",       cascade_flags[7]);
    total_cascade_flags = std::count(cascade_flags.begin(), cascade_flags.end(), true);
  } catch(const libconfig::SettingNotFoundException &nfex) {
    std::cout << "No cascade flags found. There will be no cascade properties being saved" << std::endl;
  }
  if( !total_cascade_flags ){ std::cout << "There will be no cascade properties being saved" << std::endl; }
  
  int total_scatter_flags;
  std::array<bool,13> scatter_flags{0};
  try{
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
    total_scatter_flags = std::count(scatter_flags.begin(), scatter_flags.end(), true);

  } catch(const libconfig::SettingNotFoundException &nfex) {
    std::cout << "No scatter flags found. There will be no scatter properties being saved" << std::endl;
  }
  if( !total_scatter_flags){ std::cout << "There will be no scatter properties being saved" << std::endl; }
  const bool save_time_profiles = {scatter_flags[10] && scatter_flags[11] && scatter_flags[12]};

  if( !total_cascade_flags && !total_scatter_flags){
    std::cerr << "There will be no output from this event! Please enable at least a saving flag." << std::endl;
    exit(1);
  }
  

  // std::cout << scatter_flags << std::endl;
  // std::cout << cascade_flags << std::endl;
  // We prepare the event

  Cascade1D* event;
  std::string unique_id;
  int tx_number = lab.Transmitters().size();
  int rx_number = lab.Receivers().size();
  int cs_number = cascade_list.size();
  
  for (int k = 0; k < cs_number; k++){
    Cascade& cs = cascade_list[k];
    for (int i = 0; i < tx_number; i++){
      auto tx = lab.Transmitters()[i];
      for (int j = 0; j < rx_number; j++){
        auto rx = lab.Receivers()[j];
        
        if(j == 0){
          event = new Cascade1D(tx,rx,cs);

          if(position == "direction"){
            event -> SetInDirection();
            std::cout << "Check" << std::endl;
          } else if (position == "reflectivity"){
            // event -> Transparency( event -> Density ());
            // event -> Se
          }
          

        } else {
          // Update RX,
          // make sure it updates the direction to the antennas in the appropiate function
          event -> UpdateRX(rx);
        }
        // If flag 
        event -> SetInConstIce();
        event -> RunScatter(save_time_profiles);

        unique_id = path_out + "Cascade1D_" + identifier +
                  "_TX_" + std::to_string(i) +
                  "_RX_" + std::to_string(j) + 
                  "_CS_" + std::to_string(k);
                  // You can use the cascade identifiers, too. 

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