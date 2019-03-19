#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include "multi_dim_gsl.h"

using namespace std;

double alpha_s(double t);

//Unintegrated gg splitting function
double P_gg(double z)
{
  double CA=3.;
  return 2.*CA*(z/(1.-z)+1./2.*z*(1.-z));
}

//Integrated gg splitting function
double int0_P_gg(double t0, double t2)
{
  double T=t0/t2;
  double CA=3.;
  return alpha_s(t2)/2./M_PI*2.*CA*(2.*T*T*T/3.-T*T+4.*T+2.*(log(1.-T)-log(T))-11./6.);
}

//Up to x Integrated g splitting function 
double intx_P_gg(double x, double t0, double t2)
{
  double T=t0/t2;
  double CA=3.;
  return alpha_s(t2)/2./M_PI*2.*CA/6.*(T*(12.+T*(2.*T-3.))+x*(x*(3.-2.*x)-12.)-6.*log(T*(x-1.)/(x*(T-1.))));
}

struct find_x2_params
{
  double t0, t2, R, denom;
};

//Find next z fraction function
double find_x2(double x, void *params)
{
  struct find_x2_params *p = (struct find_x2_params*) params;

  double t0=p->t0;
  double t2=p->t2;
  double R=p->R;
  double denom=p->denom;

  int limit=100;
  IntegrationWorkspace wsp(limit);
  double numer_result, numer_abserr;
  double z_lo=t0/t2;
  double z_hi=x;
  auto numer = make_gsl_function( [&](double z) {return alpha_s(t2*z*(1.-z))/2./M_PI * P_gg(z) ;} );
  gsl_integration_qags(numer, z_lo, z_hi, 0, 1e-7, limit,
		       wsp, &numer_result, &numer_abserr);

  return numer_result/denom-R;
}

//Next z fraction iterator
double generate_x2(double t0, double t2, double R)
{
  double x_lo=t0/t2, x_hi=1.-t0/t2; 

  //Integrate denominator (this one time)
  int limit=100;
  double denom_result, denom_abserr;
  IntegrationWorkspace wsp(limit); 
  auto denom = make_gsl_function( [&](double z) {return alpha_s(t2*z*(1.-z))/2./M_PI * P_gg(z) ;} );
  gsl_integration_qags(denom, x_lo, x_hi, 0, 1e-7, limit,
		       wsp, &denom_result, &denom_abserr);

  //Root finder
  int status;
  int iter=0, max_iter=100;

  double r=0.;
 
  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc (T); 

  gsl_function F;
  struct find_x2_params params = {t0,t2,R,denom_result};
  F.function = &find_x2;
  F.params = &params;

  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
 
  do {
    //cout << " at step= " << iter << " LowF= " << GSL_FN_EVAL(&F,x_lo) << " HighF= " << GSL_FN_EVAL(&F,x_hi) << endl;
    
    iter++;
    status=gsl_root_fsolver_iterate(s);
    r=gsl_root_fsolver_root(s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);

    //if (status == GSL_SUCCESS)
      //printf ("Converged:\n");

  } while (status==GSL_CONTINUE && iter<max_iter); 

  gsl_root_fsolver_free (s);

  return r;  
  
}

