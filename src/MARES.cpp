/* Radar computation - return electric field of a list of events*/

#include "cascade1D.hh"

using namespace libconfig;

int main(int argc, char** argv){
  
  // I/O Setup

  // If no argument with a specific config file has been passed, 
  // expect a "MARES_example.cfg file". 
  const char* config_file;
  assert(argc == 2 && "Missing config file, default can be found in examples/example.cfg");
//  argc < 2 ? config_file = "../examples/example.cfg": config_file = argv[1];
  config_file = argv[1];

  Config cfg;
  load_config_file(cfg, config_file);

  // By default, the result is stored in the same folder as the executable.
  std::string path_out = "../output/";
  if(cfg.lookupValue("path_out", path_out)){
    std::cout << "The output will be saved in: " << path_out << std::endl;
  } else {
    std::cout << "The output will be saved locally." << std::endl;
  }

  std::string identifier = "example";
  if(cfg.lookupValue("name", identifier)){
    std::cout << "Starting run: " << identifier << std::endl;
  }
// We keep loading stuff from our config file.
  Setting &root = cfg.getRoot();

  // We check what are the run options
  const Setting& options_list = root["options"];

    // Default position mode
  std::string position = "direction";
  options_list.lookupValue("position",  position);

    // Default propagation mode
  std::string propagation = "const";
  options_list.lookupValue("propagation", propagation);

  // We check the saving flags,
  //  aka what properties are going to be stored after the event is run.

    // Cascade flags
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
  
    // Scatter flags
  int total_scatter_flags;
  std::array<bool,17> scatter_flags{0};
  try{
    const Setting& scatter_flag_list = root["save"]["scatter"];
    scatter_flag_list.lookupValue("duration",             scatter_flags[0]);
    scatter_flag_list.lookupValue("voltage",              scatter_flags[1]);
    scatter_flag_list.lookupValue("power",                scatter_flags[2]);
    scatter_flag_list.lookupValue("target_cs",            scatter_flags[3]);
    scatter_flag_list.lookupValue("radar_cs",             scatter_flags[4]);
    scatter_flag_list.lookupValue("points_position",      scatter_flags[5]);
    scatter_flag_list.lookupValue("points_phase",         scatter_flags[6]);
    scatter_flag_list.lookupValue("points_arrival_time",  scatter_flags[7]);
    scatter_flag_list.lookupValue("points_attenuation",   scatter_flags[8]);
    scatter_flag_list.lookupValue("points_polarization",  scatter_flags[9]);
    scatter_flag_list.lookupValue("phase_vs_time",        scatter_flags[10]);
    scatter_flag_list.lookupValue("radar_vs_time",        scatter_flags[11]);
    scatter_flag_list.lookupValue("voltage_vs_time",      scatter_flags[12]);
    scatter_flag_list.lookupValue("e_field_vs_time",      scatter_flags[13]);
    scatter_flag_list.lookupValue("radar_angles",         scatter_flags[14]);
    scatter_flag_list.lookupValue("angularfreq_vs_time",  scatter_flags[15]);
    scatter_flag_list.lookupValue("transmit_signal",      scatter_flags[16]);

    total_scatter_flags = std::count(scatter_flags.begin(), scatter_flags.end(), true);

  } catch(const libconfig::SettingNotFoundException &nfex) {
    std::cout << "No scatter flags found. There will be no scatter properties being saved" << std::endl;
  }
  if( !total_scatter_flags){ std::cout << "There will be no scatter properties being saved" << std::endl; }
  const bool save_time_profiles = {scatter_flags[10] || scatter_flags[11] || scatter_flags[12] || scatter_flags[13]};

  if( !total_cascade_flags && !total_scatter_flags){
    std::cerr << "There will be no output from this event! Please enable at least a saving flag." << std::endl;
    exit(1);
  }

  // We also need to check if any of our default parameters should be updated:
    try{
    const libconfig::Setting& resolution_list = root["properties"]["resolution"];
    if( resolution_list.lookupValue("dL", _dL)){ _dL *= cm;}
    if( resolution_list.lookupValue("dR", _dR)){ _dR *= mm;}
    resolution_list.lookupValue("sampling_ratio", _sampling);
    resolution_list.lookupValue("Ltot_factor", _Ltot_factor);
    resolution_list.lookupValue("Rtot_factor", _Rtot_factor);
    
  } catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "No user-defined simulation parameters found. Running with the default resolution" << std::endl;
  }


  try{
    const libconfig::Setting& physics_list  = root["properties"];
    if( physics_list["plasma"].lookupValue("lifetime", _lifetime)){_lifetime *= ns;}
    if( physics_list["plasma"].lookupValue("collision_freq", _f_coll)){ _f_coll *= THz;}
    physics_list["plasma"].lookupValue("mass_ratio", _memp);
    
    
    if( physics_list["medium"].lookupValue("refractive_index", _refindex)){
      _c_ice = c_vac/_refindex;
      _Z_ice = Z_0/_refindex;

    }
    if( physics_list["medium"].lookupValue("density", _rho_ice)){ 
      _rho_ice *= g/pow(cm,3);
      _L_0 = _X_0/_rho_ice;
    }
    if( physics_list["medium"].lookupValue("attenuation_length", _att_length)){ _att_length *= m;}

  } catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "No user-defined physics constants found. Running with the default physics constants." << std::endl;
  }

  try{
    const libconfig::Setting& physics_list  = root["properties"];
    if( physics_list["cascade"].lookupValue("ionization_energy", _E_ionization)){ _E_ionization *= eV;}
    if( physics_list["cascade"].lookupValue("moliere_radius", _r_moliere)){ _r_moliere *= cm;}
    if( physics_list["cascade"].lookupValue("deposition_energy", _E_deposition)){ _E_deposition *= MeV/(g/pow(cm,2));}
    if( physics_list["cascade"].lookupValue("critical_energy", _E_c)){ _E_c *= MeV;}
    if( physics_list["cascade"].lookupValue("radiation_column_density", _X_0)){_X_0 *= g/pow(cm,2);_L_0 = _X_0/_rho_ice;
      }    
  } catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "No user-defined cascade constants found. Running with the default physics constants." << std::endl;
  }
  // Load your detector from your config file.
  Detector lab = load_detector_config(cfg);

  // Load your cascades from your config file.
  std::vector<Cascade> cascade_list = load_cascade_config(cfg);
/* Feel free to sanitize your cascades now, if not done already
  Do your cascades fit your requirements: Position, min energy, etc?
  You can always make a file of rejected cascades!
*/

  // Fianlly, we prepare the and run all the TX-CS-RX events.
  Cascade1D* event;
  std::string unique_id;
  int tx_number = lab.Transmitters().size();
  int rx_number = lab.Receivers().size();
  std::cout<< "Number of Transmitters: "<<tx_number <<std::endl;
  std::cout<< "Number of receivers: "<<rx_number <<std::endl;
  int cs_number = cascade_list.size();
  for (int k = 0; k < cs_number; k++){
    Cascade& cs = cascade_list[k];
    for (int i = 0; i < tx_number; i++){
      auto tx = lab.Transmitters()[i];
      for (int j = 0; j < rx_number; j++){
        auto rx = lab.Receivers()[j];

        // For every new TX-CS combination, make and position your segments. 
        if(j == 0){
          // std::cout << _f_coll << std::endl;
          // std::cout << _E_ionization << std::endl;
          event = new Cascade1D(tx,rx,cs);
          
          // Here you need to run all the other cascade functions
          // event -> Absorption(event -> Density(), event -> TX().Freq());

          if(position == "direction"){
            event -> SetInDirection();
          } else if (position == "reflectivity"){
            event -> SetInMaxReflectivity();
          } else {
            std::cerr << "Unknown position mode, stopping the run now" << std::endl;
            exit(1);
          }
          
        // After the first time, you only need to update the RX.  
        } else {
          event -> Cascade1D::UpdateRX(rx);
        }

        // Run the segment propagation
        if(propagation == "const"){
          event -> SetInConstIce();
        } else if(propagation == "iceraytracing"){
          event -> SetWithIRT();
        } else if(propagation == "beam") {
          const Setting& plane = options_list["interface"];
          double n_out= options_list.lookup("n_out");
          double n_beam= options_list.lookup("n_beam");
          std::vector<double> interface = {(double) plane[0], plane[1], plane[2], plane[3]};
          // The first 3 values of the interface are unit vectors of the normal of the plane, 
          // but the last one should have units!
          interface[3] *= m;
          event -> SetInBeam(interface, n_out, n_beam);
        } else {
          std::cerr << "Unknown propagation mode, stopping the run now" << std::endl;
          exit(1);
        }
        
        // Run the scatter proper
        event -> RunEvent(save_time_profiles);

        // write output 
        unique_id = path_out + "Cascade1D_" + identifier +
                  "_TX_" + std::to_string(i) +
                  "_RX_" + std::to_string(j) + 
                  "_CS_" + std::to_string(k);
                  // You can use the cascade identifiers, too. 

        event -> Cascade::save_output_files(unique_id, cascade_flags);
        event -> Scatter::save_output_files(unique_id, scatter_flags);

      }
      // Housekeeping!
      delete event;
    }
  }

}
// End
