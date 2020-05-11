#ifndef FOURVECTOR_H
#define FOURVECTOR_H

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <climits>

using std::cout;
using std::cerr;
using std::endl;

static double rounding_error = 1e-6;
static double maxN = double(pow(2.0,31.0) - 1.0) ;
static double a_very_large_number = maxN;

class FourVector
{
  public:

    FourVector() { tv=xv=yv=zv=0.0; }; //default constructor

    FourVector(const FourVector& srv): xv(srv.xv) , yv(srv.yv), zv(srv.zv), tv(srv.tv) {};  // copy constructor

    FourVector(double a[4])  // constructor with array input
    {
      tv=a[0];
      xv=a[1];
      yv=a[2];
      zv=a[3];
    };

    FourVector(double x_in, double y_in, double z_in, double t_in) { // constructor with component input
      tv=t_in;
      xv=x_in;
      yv=y_in;
      zv=z_in;
    };

    void clear() { tv=xv=yv=zv=0.0; }

    // constructors do all sets
    void Set(double x_in, double y_in, double z_in, double t_in) {
      tv=t_in;
      xv=x_in;
      yv=y_in;
      zv=z_in;
    };

    void Set(double a[4]) {
      tv=a[0];
      xv=a[1];
      yv=a[2];
      zv=a[3];
    };

    void print() { std::cout << xv << " " << yv << " " << zv << " " << tv << std::endl; };

    // all gets are done with name calls e.g., vec.x()
    double x() const { return(xv); };
    double y() const { return(yv); };
    double z() const { return(zv); };
    double t() const { return(tv); };

    const double comp(int i) const {
      switch (i) {
      case 0:
	     return(tv);
	      break;
      case 1:
	     return(xv);
	      break;
      case 2:
	     return(yv);
	      break;
      case 3:
	     return(zv);
	      break;
      default:
	     cout << " component index beyond 0-3! Returning garbage ..." << endl ;
	      return(a_very_large_number);
	       break;
      }
    }

    double plus() { return ( (zv+tv)/sqrt(2.0) ); };
    double minus() { return ( (tv-zv)/sqrt(2.0) ); };

    double rapidity() {
      if (this->minus()>0.0) return ( std::log(this->plus()/this->minus() )/2.0  );
      cout << endl << "ERROR: z component exceeds t component, cannot calculate rapidity" << endl;
      return (0);
    };

    double phi() {
      if ( fabs(x())<rounding_error && fabs(y())<rounding_error ) {
	       return 0;
      }
      double phi=atan2( y(), x() );
      while (phi<0) phi+=2.0*M_PI;
      return phi;
    };

    double m() { return sqrt(tv*tv - xv*xv - yv*yv - zv*zv); };
    double m2() { return this->m()*this->m(); };

    double p3Tabs() { return sqrt(xv*xv + yv*yv); };
    double p3Tabs2() { return this->p3Tabs()*this->p3Tabs(); };

    double p3abs() { return sqrt(xv*xv + yv*yv + zv*zv); };
    double p3abs2() { return this->p3abs()*this->p3abs(); };


    FourVector &operator=(FourVector &c) {
      tv = c.t();
      xv = c.x();
      yv = c.y();
      zv = c.z();
      return (*this);
    };

    FourVector &operator=(const FourVector &c) {
      tv = c.tv;
      xv = c.xv;
      yv = c.yv;
      zv = c.zv;
      return (*this);
    };

    inline FourVector operator-() const {FourVector tmp; tmp.xv = -xv; tmp.yv = -yv;
      tmp.zv = -zv; tmp.tv = -tv; return tmp;}
    inline FourVector& operator+=(const FourVector& v) {xv += v.xv; yv += v.yv; zv += v.zv;
      tv += v.tv; return *this;}
    inline FourVector& operator-=(const FourVector& v) {xv -= v.xv; yv -= v.yv; zv -= v.zv;
      tv -= v.tv; return *this;}
    inline FourVector operator+(const FourVector& v) const {
      FourVector tmp = *this; return tmp += v;}
    inline FourVector operator-(const FourVector& v) const {
      FourVector tmp = *this; return tmp -= v;}
    inline FourVector& operator*=(double f) {xv *= f; yv *= f; zv *= f;
      tv *= f; return *this;}
    inline FourVector& operator/=(double f) {xv /= f; yv /= f; zv /= f;
      tv /= f; return *this;}
    inline FourVector operator*(double f) const {FourVector tmp = *this; return tmp *= f;}
    inline FourVector operator/(double f) const {FourVector tmp = *this; return tmp /= f;}
    inline double operator*(const FourVector& v) const {return tv*v.tv - xv*v.xv - yv*v.yv - zv*v.zv;}

    void rotate_around_z(double theta) {
      double new_xv, new_yv;
      new_xv = xv*cos(theta) - yv*sin(theta);
      new_yv = yv*cos(theta) + xv*sin(theta);
      xv = new_xv;
      yv = new_yv;
    };

  private:
    double xv,yv,zv,tv;

};

#endif // FOURVECTOR_H
