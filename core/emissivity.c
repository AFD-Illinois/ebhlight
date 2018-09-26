/******************************************************************************
 *                                                                            *
 * JNU.C                                                                      *
 *                                                                            *
 * FORMULAS FOR SPECIFIC AND SOLID ANGLE-INTEGRATED EMISSIVITIES              *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
// jnu: dE/(dt dnu dV dOmega)
// Jnu: dE/(dt dnu dV)

double jnu_gray(double nu, const struct of_microphysics *m);
double jnu_brem(double nu, const struct of_microphysics *m);
double jnu_synch(double nu, const struct of_microphysics *m, double theta);
double Jnu_gray(double nu, const struct of_microphysics *m);
double Jnu_brem(double nu, const struct of_microphysics *m);
double Jnu_synch(double nu, const struct of_microphysics *m);
double J_gray(double Ne, double Thetae);
double J_brem(double Ne, double Thetae);
double J_synch(double Ne, double Thetae, double B);

void init_emiss_tables();

void init_emissivity() {
  #if SYNCHROTRON
  init_emiss_tables();
  #endif
}

double jnu(double nu, const struct of_microphysics *m, double theta)
{
  double jnu = 0.;

  #if GRAYABSORPTION
    jnu += jnu_gray(nu, m);
  #endif
  #if BREMSSTRAHLUNG
    jnu += jnu_brem(nu, m);
  #endif
  #if SYNCHROTRON
    jnu += jnu_synch(nu, m, theta);
  #endif

  return jnu;
}

double Jnu(double nu, const struct of_microphysics *m)
{
  double Jnu = 0.;

  #if GRAYABSORPTION
    Jnu += Jnu_gray(nu, m);
  #endif
  #if BREMSSTRAHLUNG
    Jnu += Jnu_brem(nu, m);
  #endif
  #if SYNCHROTRON
    Jnu += Jnu_synch(nu, m);
  #endif

  return Jnu;
}
double get_J(double Ne, double Thetae, double Bmag)
{
  double J  = 0.;

  #if GRAYABSORPTION
    J += J_gray(Ne, Thetae);
  #endif
  #if BREMSSTRAHLUNG
    J += J_brem(Ne, Thetae);
  #endif
  #if SYNCHROTRON
    J += J_synch(Ne, Thetae, Bmag);
  #endif

  return J;
}

// Rybicki & Lightman 1979
double jnu_brem(double nu, const struct of_microphysics *m)
{
  double Ne = m->Ne;
  double Thetae = m->Thetae;
  double Te = Thetae*ME*CL*CL/KBOL;
  double rel = (1. + 4.4e-10*Te);
  double x,efac;
  double gff = 1.2;

  x = HPL*nu/(KBOL*Te);

  if(x < 1.e-3) {
    efac = (24 - 24*x + 12*x*x - 4.*x*x*x + x*x*x*x)/24.;
  } else {
    efac = exp(-x);
  }

  double jv = 1./(4.*M_PI)*pow(2,5)*M_PI*pow(EE,6)/(3.*ME*pow(CL,3));
  jv *= pow(2.*M_PI/(3.*KBOL*ME),1./2.);
  jv *= pow(Te,-1./2.)*Ne*Ne;
  jv *= efac*rel*gff;

  return jv;
}

double Jnu_brem(double nu, const struct of_microphysics *m)
{
  return 4.*M_PI*jnu_brem(nu, m);
}

double J_brem(double Ne, double Thetae)
{
  double Te = Thetae*ME*CL*CL/KBOL;
  double gff = 1.2;
  double rel = (1. + 4.4e-10*Te);
  double Jb = sqrt(2.*M_PI*KBOL*Te/(3.*ME));
  Jb *= pow(2,5)*M_PI*pow(EE,6)/(3*HPL*ME*pow(CL,3));
  Jb *= Ne*Ne*gff*rel;

  return Jb;
}

double jnu_integrand_synch(double th, void *params);
double F_eval_synch(double Thetae, double Bmag, double nu);
double F_synch[NU_BINS_EMISS+1], K2[NU_BINS_EMISS+1];
static double lK_min, dlK;
static double lT_min, dl_T;
#define THETAE_MIN  0.3 // Synchrotron fitting function invalid for low Thetae
#define EPSABS 0.
#define EPSREL 1.e-6
#define KMIN (0.002)
#define KMAX (1.e7)
#define TMIN (THETAE_MIN)
#define TMAX (1.e2)
void init_emiss_tables()
{
  int k;
  double result,err,K,T;
  gsl_function func;
  gsl_integration_workspace *w;

  func.function = &jnu_integrand_synch;
  func.params = &K;

  lK_min = log(KMIN);
  dlK = log(KMAX/KMIN)/(NU_BINS_EMISS);

  lT_min = log(TMIN);
  dl_T = log(TMAX/TMIN)/(NU_BINS_EMISS);

  //  build table for F(K) where F(K) is given by
  // \int_0^\pi ( (K/\sin\theta)^{1/2} + 2^{11/12}(K/\sin\theta)^{1/6})^2 \exp[-(K/\sin\theta)^{1/3}]
  // so that J_{\nu} = const.*F(K)

  w = gsl_integration_workspace_alloc(1000);
  for(k = 0; k <= NU_BINS_EMISS; k++){
    K = exp(k*dlK + lK_min);
    gsl_integration_qag(&func, 0., M_PI/2., EPSABS, EPSREL, 1000,
      GSL_INTEG_GAUSS61, w, &result, &err);
    F_synch[k] = log(4*M_PI*result);
  }
  gsl_integration_workspace_free(w);

  // Build table for quick evaluation of the bessel function K2 for emissivity
  for(k = 0; k <= NU_BINS_EMISS; k++){
    T = exp(k*dl_T + lT_min);
    K2[k] = log(gsl_sf_bessel_Kn(2,1./T));
  }
}

double linear_interp_K2_synch(double Thetae)
{
  int i;
  double di,lT;

  lT = log(Thetae);

  di = (lT - lT_min)/dl_T;
  i = (int)di;
  di = di - i;

  return exp((1.-di)*K2[i] + di*K2[i+1]);
}

double linear_interp_F_synch(double K)
{
  if (K < KMIN || K > KMAX) return 0.;

  int i;
  double di,lK;

  lK = log(K);

  di = (lK - lK_min)/dlK;
  i = (int)di;
  di = di-i;

  return exp((1.-di)*F_synch[i] + di*F_synch[i+1]);
}

double K2_eval_synch(double Thetae)
{
  if(Thetae < THETAE_MIN) return 0.;
  if(Thetae > TMAX) return 2.*Thetae*Thetae;

  return linear_interp_K2_synch(Thetae);
}

#define CST 1.88774862536 // 2^{11/12}
double jnu_synch(double nu, const struct of_microphysics *m, double theta)
{
  double Thetae = m->Thetae;
  double Ne = m->Ne;
  double B = m->B;
  double K2,nuc,nus,x,f,j,sth,xp1, xx ;

  if(Thetae < THETAE_MIN) return 0.;

  K2 = K2_eval_synch(Thetae);

  nuc = EE*B/(2.*M_PI*ME*CL) ;
  sth = sin(theta) ;
  nus = (2./9.)*nuc*Thetae*Thetae*sth;
  if(nu > 1.e12*nus) return(0.) ;
  x = nu/nus ;
  xp1 = pow(x,1./3.);
  xx = sqrt(x) + CST*sqrt(xp1);
  f = xx*xx;
  j = (M_SQRT2*M_PI*EE*EE*Ne*nus/(3.*CL*K2)) * f * exp(-xp1) ;
  
  if (isinf(j) ||isnan(j) || j > 1.e200 || j < -1.e200)
    return 0.;

  return (j);
}
#undef CST

double Jnu_synch(double nu, const struct of_microphysics *m)
{
  double Thetae = m->Thetae;
  double Ne     = m->Ne;
  double Bmag   = m->B;
  if (Thetae < THETAE_MIN)
    return 0.;
  
  double nuC = EE*Bmag/(2*M_PI*ME*CL);

  double F_interp  = F_eval_synch(Thetae, Bmag, nu);
  double K2_interp = K2_eval_synch(Thetae);

  double Jnu = sqrt(2.)*M_PI*EE*EE*Ne*(2./9.*nuC)*Thetae*Thetae;
  Jnu       /= (3.*CL*K2_interp);
  Jnu       *= F_interp;

  if (isnan(Jnu) || isinf(Jnu))
    return 0.;

  return Jnu;
}

double J_synch(double Ne, double Thetae, double Bmag)
{
  double K2 = K2_eval_synch(Thetae);
  
  double J = 32.*sqrt(2.)*(10. + pow(2,5./6.) + 4*pow(2,11./12.));
  J *= Bmag*Bmag*Ne*M_PI*EE*EE*EE*EE*pow(Thetae, 5);
  J /= (81.*CL*CL*CL*K2*ME*ME*M_PI);

  if (isinf(J) || isnan (J) || J > 1.e200 || J < -1.e200)
    return 0.;

  return J;
}

#define CST 1.88774862536 /* 2^{11/12} */
double jnu_integrand_synch(double th, void *params)
{
  double K = *(double *)params;
  double sth = sin(th);
  double x = K/sth;

  if(sth < 1.e-150 || x > 2.e8) return 0.;

  return sth*sth*pow(sqrt(x) + CST*pow(x,1./6.),2.)*exp(-pow(x,1./3.));
}
#undef CST

#define KFAC  (9*M_PI*ME*CL/EE)
double F_eval_synch(double Thetae, double Bmag, double nu)
{
  double K, x;
  double linear_interp_F_synch(double);

  K = KFAC*nu/(Bmag*Thetae*Thetae);

  if(K > KMAX) {
    return 0.;
  } else if(K < KMIN) {
    /* use a good approximation */
    x = pow(K, 0.333333333333333333);
    return (x*(37.67503800178 + 2.240274341836*x));
  } else if (isnan(K)) {
    return 0;
  } else {
    return linear_interp_F_synch(K);
  }
}
#undef KFAC
#undef KMIN
#undef KMAX
#undef EPSABS
#undef EPSREL
#undef THETAE_MIN
#undef TMIN
#undef TMAX

double jnu_gray(double nu, const struct of_microphysics *m)
{
  double opacity = MP*(m->Ne)*kappa;

  return opacity*Bnu_inv(nu, m)*nu*nu*nu;
}

double Jnu_gray(double nu, const struct of_microphysics *m)
{
  return 4.*M_PI*jnu_gray(nu, m);
}

double J_gray(double Ne, double Thetae)
{
  return 0.;
}
#endif // RADIATION

