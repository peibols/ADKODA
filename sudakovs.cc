#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include "multi_dim_gsl.h"

using namespace std;

double P_gg(double z);

//Running alpha_s(t)
double alpha_s(double t)
{
  double Lam=0.2;
  double Lam2=Lam*Lam;
  double nf=3.;
  double b=(33.-2.*nf)/12./M_PI;
  return 1./(b*log(t/Lam2));
}

//Sudakov integral params
struct sudakov_integral_params
{
  double t0;
};

//Sudakov integral, with analytic z integral
double sudakov_integral(double T, void *params)
{
 struct sudakov_integral_params *p = (struct sudakov_integral_params*) params;

 double t0=p->t0;

 double CA=3.;

 double int_temp=2.*CA*1./T*(-11./6.+4.*T-T*T+2.*T*T*T/3.+2.*(log(1.-T)-log(T)));
 int_temp*=alpha_s(t0/T);
 return int_temp;
}

//Compute log of the Sudakov
double zanalytic_Sudakov_log(double t, double t0)
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  F.function=&sudakov_integral;
  F.params=&t0;

  double result, error;
  gsl_integration_qags (&F, t0/t, 0.5, 0, 1e-7, 1000,
                        w, &result, &error);
  
  return -1./2./M_PI*result;
}

//Compute log of the Sudakov through double integral
double Sudakov_log(double t, double t0)
{
  size_t limit = 100;
  double result, abserr, inner_result, inner_abserr;

  IntegrationWorkspace wsp1(limit);
  IntegrationWorkspace wsp2(limit);

  auto outer = make_gsl_function( [&](double tp) {
    double z_lo=t0/tp;
    double z_hi=1.-t0/tp;
    auto inner = make_gsl_function( [&](double z) {return -1./tp * alpha_s(tp*(z*(1.-z)))/(2.*M_PI) * P_gg(z) ;} );
    //inner.params=&t0;
    gsl_integration_qags(inner, z_lo, z_hi, 0, 1e-7, limit, 
		         wsp1, &inner_result, &inner_abserr);
    return inner_result;
  } );

  double t_lo=2.*t0;
  double t_hi=t;
  gsl_integration_qags(outer, t_lo, t_hi, 0, 1e-7, limit,
		       wsp2, &result, &abserr);  

  return result;
}

//Compute the Sudakov
double Sudakov(double t, double t0)
{
  return exp(Sudakov_log(t,t0));
}

//Compute the zanalytic Sudakov
double zanalytic_Sudakov(double t, double t0)
{
  return exp(zanalytic_Sudakov_log(t,t0));
}

struct find_t2_params
{
  double t0, R, numer;
};

//Find next virtuality function
double find_t2(double t2, void *params)
{
  struct find_t2_params *p = (struct find_t2_params*) params;

  double t0=p->t0;
  double R=p->R;
  double numer=p->numer;

  return numer-Sudakov_log(t2,t0)-log(R);     
}

//Next virtuality iterator
double generate_t2(double t0, double t1, double R, double numer)
{
  //Root finder
  int status;
  int iter=0, max_iter=100;

  double r=0.;
  double x_lo=2.*t0, x_hi=t1;
 
  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc (T); 

  gsl_function F;
  struct find_t2_params params = {t0,R,numer};
  F.function = &find_t2;
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
