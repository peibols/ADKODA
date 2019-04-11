

void EvolveParton(Parton parton, double tmax, bool iter)
{

  double E=parton.p.t();
  double cutoff=1./2.*(1.-sqrt(E*E-m_min*m_min)/E); 

  double t_temp, z_temp; 
  while (true) {
    
    //Pick t
    double Rt=dis(gen);
    double t_temp=inv_prim_envel(prim_envel(tmax,cutoff)+log(Rt));

    //Pick z
    double Rz=dis(gen);
    double z_temp=inv_zprim_envel(zprim_envel(cutoff)+Rz*(zprim_envel(1.-cutoff)-zprim_envel(cutoff)));     

    //Evaluate ratio
    double Rhm=dis(gen);
    double ratio=integrand(t_temp,z_temp)/envel(t_temp,z_temp); 

    if (ratio>Rhm) break;
  }

  //Put on-mass-shell if below m_min
  if (t_temp<m_min*m_min) t_temp=0., z_temp=1., parton.set_status(-1);

  parton.set_virt(t_temp);
  parton.set_z(z_temp);
   
}

//g(t,z)
double envel(double t, double z)
{
  double CA=3.;
  double kernel=CA*(1./z+1./(1.-z));
  double max_alphas=0.5;
  return kernel * max_alphas/2./M_PI / t;
}

//f(t,z)
double integrand(double t, double z)
{
  double kernel=P_gg(z);
  return kernel * alpha_s(t*z(1.-z))/2./M_PI / t;
}

//G(t)
double prim_envel(double t, double cutoff)
{
  double CA=3.;
  double max_alphas=0.5;
  double int_kernel=CA*4.*arctanh(1.-2.*cutoff);
  double result= int_kernel * max_alphas/2./M_PI * log(t);
  return result; 
}

//G^-1
double inv_prim_envel(double y, double cutoff)
{
  double CA=3.;
  double max_alphas=0.5;
  double int_kernel=CA*4.*arctanh(1.-2.*cutoff);
  double result= exp( y / (int_kernel * max_alphas/2./M_PI) ); 
  return result;
}

//H(z)
double zprim_envel(double z);
{
  return log(z)-log(1.-z);
}

//H^-1
double inv_zprim_envel(double y)
{
  return exp(y)/(exp(y)+1.);     
}
