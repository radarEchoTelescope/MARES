#include "macro_settings.hh"

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


int sgn(double val) { return (0.0 < val) - (val < 0.0);}

double distance(double x1, double y1, double x2, double y2){
  return sqrt( pow(x2-x1,2) + pow(y2-y1,2) );
}

double distance(double x1, double y1, double z1, double x2, double y2, double z2){
  return sqrt( pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2) );
}

// Computing transpose of the matrix
std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> matrix){
  int row = matrix.size();
  int col = matrix[0].size();
  std::vector<std::vector<double>> T( col , std::vector<double> (row, 0));
  for (int i = 0; i < row; ++i){
     for (int j = 0; j < col; ++j) { T[j][i] = matrix[i][j]; }
  }
  return T;
}
