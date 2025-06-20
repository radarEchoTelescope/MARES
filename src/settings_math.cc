#include "settings.hh"

// Math tools ------------------------------------------------------------------

int sgn(double val) { return (0.0 < val) - (val < 0.0);}

double acos2(double x1, double x2, double y1, double y2){
  double angle;
  angle = atan2( x1*y2 - x2*y1, x1*y1 + x2*y2 );
  // angle = atan2( det(a,b), dot(a,b) );
  return angle;
}

double sinc(double x){
  double o = 1;
  if(x != 0){ o = sin(x)/x; }
  return o;
}

double norm(std::vector<double> vec_a){
  double n{0};
  for (auto a : vec_a){ n += pow(a,2); }
  return std::sqrt(n);
}

std::vector<double> direction(std::vector<double> a, std::vector<double> b){
  std::vector<double> out;
  std::transform(a.begin(), a.end(),      //input1
                 b.begin(),               //input2
                 std::back_inserter(out),                //output
                 [](double x, double y){ return y-x; });
  return out;
}

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
  double n_a = norm(a) , n_b = norm(b), p = 0;
  if ( (n_a && n_b) != 0){p = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])/ (n_a * n_b);}
  return p;
}

double shortest_distance(std::vector<double> plane, std::vector<double> point){
  double out;
  out  = abs(plane[0]*point[0] + plane[1]*point[1] + plane[2]*point[2] + plane[3])
        / norm( {plane[0], plane[1], plane[2] });
  return out;
}

std::vector<double> normalize(std::vector<double> a){
  double n = norm(a);
  if(n != 0){
    a[0] = a[0]/n;
    a[1] = a[1]/n;
    a[2] = a[2]/n;
  } else {
    a[0] = 0.;
    a[1] = 0.;
    a[2] = 0.;
  }
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

double dot_product(std::vector<double> a, std::vector<double> b){
    // Check for the size of u and v;
    double dot;

    dot=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    return dot;
}

std::vector<double> cross_product(double a1, double a2, double a3, double b1, double b2, double b3){
  std::vector<double> a = {a1, a2, a3};
  std::vector<double> b = {b1, b2, b3};
  std::vector<double> ab = cross_product(a,b);
  return ab;
}

// Rotates the vector counter clockwise over given angle around rotAxis.
// Counter clockwise is defined wrt observer which has rotAxis pointing towards it.
std::vector<double> rotate(std::vector<double> vector_in, std::vector<double> rotAxis, double theta){
    // The rotation angle must be given in rad
    
    // The incoming vector
    double _x = vector_in[0];
    double _y = vector_in[1];
    double _z = vector_in[2];

    // The axis components
    double ux = rotAxis[0];
    double uy = rotAxis[1];
    double uz = rotAxis[2];

    // The rotated vector components
    double xr = ( cos(theta) + ux*ux*(1.-cos(theta)) ) * _x
                    + ( ux*uy*(1.-cos(theta)) - uz*sin(theta) ) * _y
                    + ( ux*uz*(1.-cos(theta)) + uy*sin(theta) ) * _z;
    double yr = ( uy*ux*(1.-cos(theta)) + uz*sin(theta) ) * _x
                    + ( cos(theta) + uy*uy*(1.-cos(theta)) ) * _y
                    + ( uy*uz*(1.-cos(theta)) - ux*sin(theta) ) * _z;
    double zr = ( uz*ux*(1.-cos(theta)) - uy*sin(theta) ) * _x
                    + ( uz*uy*(1.-cos(theta)) + ux*sin(theta) ) * _y
                    + ( cos(theta) + uz*uz*(1.-cos(theta)) ) * _z;
    std::vector<double> rotVect{xr, yr, zr};
    return rotVect;
}

std::vector<double> find_perpendicular(std::vector<double> vector_in){
      // The outgoing vector
    double x = -vector_in[1];
    double y = vector_in[0];
    double z = 0;
    std::vector<double> vector_out{x, y, z};

    if(projection(vector_in,vector_out) == 0){return vector_out;}
    else{
      double x = -vector_in[2];
      double y = 0;
      double z = vector_in[1];
      std::vector<double> vector_out{x, y, z};
      return vector_out;
    }
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