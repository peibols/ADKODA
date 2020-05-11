#include "Util.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

namespace Util {

  double Delta(FourVector v1, FourVector v2) {
    double dphi = std::acos((v1.x()*v2.x() + v1.y()*v2.y())/v1.p3Tabs()/v2.p3Tabs());
    return std::sqrt(std::pow(v1.rapidity() - v2.rapidity(), 2.) + std::pow(dphi, 2.));
  }

double m(const FourVector& v1, const FourVector& v2) {
    double m2 = std::pow(v1.t() + v2.t(), 2.) - std::pow(v1.x() + v2.x(), 2.)
       - std::pow(v1.y() + v2.y(), 2.) - std::pow(v1.z() + v2.z(), 2.);
    return (m2 > 0.) ? std::sqrt(m2) : 0.;
  }

double m2(const FourVector& v1, const FourVector& v2) { return std::pow(m(v1, v2), 2.); }

FourVector Cross(const FourVector& v1, const FourVector& v2) {
  return FourVector(v1.y()*v2.z()-v1.z()*v2.y(),
                    v1.z()*v2.x()-v1.x()*v2.z(),
                    v1.x()*v2.y()-v1.y()*v2.x(), 0.);
}

FourVector BoostForCS(FourVector p, FourVector q) {
  double rsq = p.m();
  double v0 = (p.t()*q.t()-p.x()*q.x()-p.y()*q.y()-p.z()*q.z())/rsq;
  double c1 = (q.t()+v0)/(p.t()+rsq);
  return FourVector(q.x() - p.x()*c1,
                    q.y() - p.y()*c1,
                    q.z() - p.z()*c1, v0);
}
FourVector BoostBackForCS(FourVector p, FourVector q) {
  double rsq = p.m();
  double v0 = (p.t()*q.t()+p.x()*q.x()+p.y()*q.y()+p.z()*q.z())/rsq;
  double c1 = (q.t()+v0)/(p.t()+rsq);
  return FourVector(q.x() + p.x()*c1,
                    q.y() + p.y()*c1,
                    q.z() + p.z()*c1, v0);
}

FourVector Boost( double b[3], FourVector p) {

  double betamod = std::sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
  double gamma = 1.0 / std::sqrt(1.0-betamod*betamod);
  double pscalb = p.x()*b[0] + p.y()*b[1] + p.z()*b[2];
  double spat = ( pscalb * gamma / (1.0 + gamma) - p.t() ) * gamma;

  double ptt = gamma * ( p.t() - pscalb );
  double ptx = p.x() + b[0] * spat;
  double pty = p.y() + b[1] * spat;
  double ptz = p.z() + b[2] * spat;

  return FourVector(ptx, pty, ptz, ptt);
}

FourVector BoostBack( double b[3], FourVector p) { //FIXME can be simplified BoostBack(v, p) == Boost(-v, p)

  double betamod = std::sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
  double gamma = 1.0 / std::sqrt(1.0-betamod*betamod);
  double pscalb = p.x()*b[0] + p.y()*b[1] + p.z()*b[2];
  double spat = ( -pscalb * gamma / (1.0 + gamma) - p.t() ) * gamma;

  double bt = gamma * ( p.t() + pscalb );
  double bx = p.x() - b[0] * spat;
  double by = p.y() - b[1] * spat;
  double bz = p.z() - b[2] * spat;

  return FourVector(bx, by, bz, bt);
}

void VelCOM(FourVector p1, FourVector p2, double b[3]){
  b[0] = ( p1.x() + p2.x() ) / ( p1.t() + p2.t() );
  b[1] = ( p1.y() + p2.y() ) / ( p1.t() + p2.t() );
  b[2] = ( p1.z() + p2.z() ) / ( p1.t() + p2.t() );
  return;
}

//Rotate a vector around a given axis with an angle (right handed).
void Rotation(FourVector &vector, double angle, double axis[3]) {
  //Normalize the axis if its wasn't.
  double k[3];
  double axis_abs = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
  for (int i = 0; i < 3; i++) k[i] = axis[i]/axis_abs;

  double vec[3], vec_prime[3];
  vec[0]=vector.x(), vec[1]=vector.y(), vec[2]=vector.z();
  vec_prime[0] = vec[0] + sin(angle) * (-k[2]*vec[1]+k[1]*vec[2]) + (1.-cos(angle)) * \
                  ((-k[2]*k[2]-k[1]*k[1])*vec[0] + k[0]*k[1]*vec[1] + k[0]*k[2]*vec[2]);
  vec_prime[1] = vec[1] + sin(angle) * (k[2]*vec[0]-k[0]*vec[2]) + (1.-cos(angle)) * \
                  ((-k[2]*k[2]-k[0]*k[0])*vec[1] + k[1]*k[0]*vec[0] + k[1]*k[2]*vec[2]);
  vec_prime[2] = vec[2] + sin(angle) * (-k[1]*vec[0]+k[0]*vec[1]) + (1.-cos(angle)) * \
                  ((-k[0]*k[0]-k[1]*k[1])*vec[2] + k[2]*k[0]*vec[0] + k[2]*k[1]*vec[1]);
  vector.Set(vec_prime[0], vec_prime[1], vec_prime[2], vector.t());
}

void AlignWithZ(FourVector &vec, double &angle, double k[3]) {
  angle = acos(vec.z()/sqrt(vec.x()*vec.x()+vec.y()*vec.y()+vec.z()*vec.z()));
  k[0]  = vec.y();
  k[1]  = -vec.x();
  k[2]  = 0.;
  Rotation(vec, angle, k);
}

int IsFile(std::string file_name) {
  FILE *temp;
  if( (temp = fopen(file_name.c_str(), "r")) == NULL ) return 0;
  else { fclose(temp); return 1; }
}/* IsFile */

// support comments in the parameters file
// comments need to start with #
std::string StringFind4(std::string file_name, std::string str_in) {
    std::string inputname = file_name;
    std::string str = str_in;

    std::string tmpfilename;
    tmpfilename = "input.default";

    // check whether the input parameter file is exist or not
    if (!IsFile(file_name)) {
        if (file_name == "") {
            fprintf(stderr, "No input file name specified.\n");
            fprintf(stderr, "Creating a default file named input.default\n");
        } else {
            cerr << "The file named " << file_name << " is absent." << endl;
            cout << "Creating " << file_name << "..." << endl;
            tmpfilename = file_name;
        }
        std::ofstream tmp_file(tmpfilename.c_str());
        tmp_file << "EndOfData" << endl;
        tmp_file.close();
        exit(1);
    }/* if isfile */

    // pass checking, now read in the parameter file
    std::string temp_string;
    std::ifstream input(inputname.c_str());
    getline(input, temp_string);  // read in the first entry

    int ind = 0;
    std::string para_name;
    std::string para_val;
    while (temp_string.compare("EndOfData") != 0) {
        // check whether it is the end of the file
        std::string para_string;
        std::stringstream temp_ss(temp_string);
        getline(temp_ss, para_string, '#');  // remove the comments
        if (para_string.compare("") != 0
                && para_string.find_first_not_of(' ') != std::string::npos) {
            // check the read in string is not empty
            std::stringstream para_stream(para_string);
            para_stream >> para_name >> para_val;
            if (para_name.compare(str) == 0) {
                // find the desired parameter
                ind++;
                input.close();
                return(para_val);
            }  /* if right, return */
        }
        getline(input, temp_string);  // read in the next entry
    }/* while */
    input.close(); // finish read in and close the file

    // the desired parameter is not in the parameter file, then return "empty"
    if (ind == 0) {
        return("empty");
    }
    // should not cross here !!!
    cout << "Error in StringFind4 !!!\n";
    return("empty");
}/* StringFind4 */

} // end namespace Util
