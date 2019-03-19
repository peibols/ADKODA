#include <cstdlib>
#include <cmath>
#include <iostream>
#include <random>

using namespace std;

double Sudakov(double t, double t0);
double generate_t2(double t0, double t1, double R);

int main() {

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.,1.);

  double t_max=1000.;
  double x_max=1.;
  double t0=2.;

  //Parton Shower
  double t1=t_max;
  double x1=x_max;

  double t2=-1000.;
  while (true)
  {
    double R=dis(gen);
    if (R>Sudakov(t1,t0)/Sudakov(2.*t0,t0))
    {
      double t2_temp=generate_t2(2.*t0,t1,R);
      if (t2_temp>2.*t0) t2=t2_temp;
      else break;
      cout << " t2= " << t2 << endl;
      t1=t2;
    }
    else break;
  }

  
  return 0;
}
