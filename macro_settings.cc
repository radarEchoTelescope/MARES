#include "macro_settings.hh"

// Math tools ------------------------------------------------------------------
int sgn(double val) { return (0.0 < val) - (val < 0.0);}

double acos2(double y, double x){
  // If you want to check your angle's range, you can do it here.
  double delta = y - x;
  delta < pi ? delta = acos(delta) : delta = - acos(delta);
  return delta;
}


// template <typename T>
//T norm(){
//  T sum{0};
double norm(std::vector<double> vec_a){
  double n{0};
  for (auto a : vec_a){ n += pow(a,2); }
  return std::sqrt(n);
}

// double norm(std::vector<double> a){
//   double n;
//   if (a.size() == 3){ n = sqrt( pow(a[0], 2) + pow(a[1], 2) + pow(a[2],2)  );  }
//   else if (a.size() == 2){ n = sqrt( pow(a[0], 2) + pow(a[1], 2) );  }
//   return n;
// }

std::vector<double> direction(std::vector<double> a, std::vector<double> b){
  std::vector<double> out;
  std::transform(a.begin(), a.end(),  //input1
                 b.begin(),               //input2
                 std::back_inserter(out),                //output
                 [](double x, double y){ return y-x; });
  // std::vector<double> v = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  return out;
}

// std::vector<double> direction(double[] a, double[] b){
//
// }

std::vector<double> direction(double x1, double y1, double z1, double x2, double y2, double z2){
  std::vector<double> a = {x1, y1, z1};
  std::vector<double> b = {x2, y2, z2};
  std::vector<double> out = direction(a,b);
  return out;
}


// std::vector<double> w = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
double distance(std::vector<double> a, std::vector<double> b){
  return norm(direction(a,b));
}

double distance(double x1, double y1, double x2, double y2){
  std::vector<double> a = {x1, y1};
  std::vector<double> b = {x2, y2};

  return norm(direction(a,b));
}

double distance(double x1, double y1, double z1, double x2, double y2, double z2){
  std::vector<double> a = {x1, y1, z1};
  std::vector<double> b = {x2, y2, z2};
  return norm(direction(a,b));
}

double projection(std::vector<double> a, std::vector<double> b){
  double n_a = norm(a) , n_b = norm(b), p;
  // Check for norm 0 and output error!!
  p = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])/ (n_a * n_b); //Normalize
  return p;
}

std::vector<double> normalize(std::vector<double> a){
  double n = norm(a);
  a[0] = a[0]/n;
  a[1] = a[1]/n;
  a[2] = a[2]/n;
  return a;
}

std::vector<double> cross_product(std::vector<double> a, std::vector<double> b){
    // Check for the size of u and v;
    std::vector<double> ab(3);

    ab[0] = a[1] * b[2] - b[1] * a[2];
    ab[1] = b[0] * a[2] - a[0] * b[2];
    ab[2] = a[0] * b[1] - b[0] * a[1];
    return ab;
}

std::vector<double> cross_product(double a1, double a2, double a3, double b1, double b2, double b3){
  std::vector<double> a = {a1, a2, a3};
  std::vector<double> b = {b1, b2, b3};
  std::vector<double> ab = cross_product(a,b);
  return ab;
}

//---- I/O functions -----------------------------------------------------------
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

// Unused ----------------------------------------------------------------------
//
// // Computing transpose of the matrix
// std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> matrix){
//   int row = matrix.size();
//   int col = matrix[0].size();
//   std::vector<std::vector<double>> T( col , std::vector<double> (row, 0));
//   for (int i = 0; i < row; ++i){
//      for (int j = 0; j < col; ++j) { T[j][i] = matrix[i][j]; }
//   }
//   return T;
// }
