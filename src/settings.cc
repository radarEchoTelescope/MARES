#include "settings.hh"

void write_1D_array(std::vector<double> array, std::string output_path, const bool trigger){
  if(trigger){
    std::ofstream output_file(output_path);
    std::ostream_iterator<double> output_iterator(output_file, "\t");
    copy(array.begin(), array.end(), output_iterator);
    output_file.close();
  }
}


void write_2D_array(std::vector<std::vector<double>> array, std::string output_path, const bool trigger){
  if(trigger){
    std::ofstream output_file(output_path);
    std::ostream_iterator<double> output_iterator(output_file, "\t");
    for (int i = 0; i < array.size(); i++){
      copy(array.at(i).begin(), array.at(i).end(), output_iterator);
      output_file << '\n';
    }
    output_file.close();
  }
}

void load_config_file(libconfig::Config& cfg, const char* config_file){
  try{
    cfg.readFile(config_file);
  } catch(const libconfig::FileIOException &fioex) {
      std::cerr << "I/O error while reading file." << std::endl;
      exit(EXIT_FAILURE);
  } catch(const libconfig::ParseException &pex) {
      std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
      exit(EXIT_FAILURE);
  }
}