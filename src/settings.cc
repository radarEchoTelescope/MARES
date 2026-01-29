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


// Function to read and save txt file input as a vector. Adapted from https://stackoverflow.com/questions/40539385/reading-data-from-a-file-and-storing-each-line-in-an-array
std::vector<double> read_1D_array (std::istream & in){

  std::vector<double> values;
  std::string line; 

  //Iterate through each line of the txt file
  for ( int line_number = 1; std::getline(in, line ); line_number++){

    try { //Convert the line to a double, and append to the global array
      double value = std::stod( line );
      values.push_back(value);
    }
    catch (std::bad_alloc &e) { //Handle allocation memory errors
      std::cerr << "Memory error at line "<<line_number<< std::endl;
      throw e ;
    }
    catch (std::exception &e) { //Handle general errors
      std::cerr <<"Error at line "<< line_number << std::endl;
    }
  }
  return values;
}
