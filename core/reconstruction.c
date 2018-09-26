/******************************************************************************
 *                                                                            *
 * RECONSTRUCTION.C                                                           *
 *                                                                            *
 * RECONSTRUCTION ALGORITHMS                                                  *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Performs the slope-limiting for the numerical flux calculation
double slope_lim(double y1, double y2, double y3)
{
  double Dqm, Dqp, Dqc, s;

  // Limiter choice now hard-coded
  lim = MC;

  // Woodward, or monotonized central, slope limiter
  if (lim == MC) {
    Dqm = 2. * (y2 - y1);
    Dqp = 2. * (y3 - y2);
    Dqc = 0.5 * (y3 - y1);
    s = Dqm * Dqp;
    if (s <= 0.)
      return 0.;
    else {
      if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
        return (Dqm);
      else if (fabs(Dqp) < fabs(Dqc))
        return (Dqp);
      else
        return (Dqc);
    }
  }

  // van Leer slope limiter
  else if (lim == VANL) {
    Dqm = (y2 - y1);
    Dqp = (y3 - y2);
    s = Dqm * Dqp;
    if (s <= 0.)
      return 0.;
    else
      return (2. * s / (Dqm + Dqp));
  }

  // Minmod slope limiter (crude but robust)
  else if (lim == MINM) {
    Dqm = (y2 - y1);
    Dqp = (y3 - y2);
    s = Dqm * Dqp;
    if (s <= 0.)
      return 0.;
    else if (fabs(Dqm) < fabs(Dqp))
      return Dqm;
    else
      return Dqp;
  }

  fprintf(stderr, "unknown slope limiter\n");
  exit(10);

  return (0.);
}

void linear_mc(double x1, double x2, double x3, double *lout, double *rout)
{
  double Dqm,Dqp,Dqc,s;

  Dqm = 2. * (x2 - x1);
  Dqp = 2. * (x3 - x2);
  Dqc = 0.5 * (x3 - x1);

  s = Dqm * Dqp;

  if (s <= 0.)
    s = 0.;
  else {
    if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
      s = Dqm;
    else if (fabs(Dqp) < fabs(Dqc))
      s = Dqp;
    else
      s = Dqc;
  }

  // Reconstruct left, right
  *lout = x2 - 0.5*s;
  *rout = x2 + 0.5*s;
}

// WENO interpolation. See Tchekhovskoy et al. 2007 (T07), Shu 2011 (S11)
// Implemented by Monika Moscibrodzka
void weno(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout)
{
  // S11 1, 2, 3
  double vr[3], vl[3];
  vr[0] =  (3./8.)*x1 - (5./4.)*x2 + (15./8.)*x3;
  vr[1] = (-1./8.)*x2 + (3./4.)*x3 + (3./8.)*x4;
  vr[2] =  (3./8.)*x3 + (3./4.)*x4 - (1./8.)*x5;

  vl[0] =  (3./8.)*x5 - (5./4.)*x4 + (15./8.)*x3;
  vl[1] = (-1./8.)*x4 + (3./4.)*x3 + (3./8.)*x2;
  vl[2] =  (3./8.)*x3 + (3./4.)*x2 - (1./8.)*x1;

  // Smoothness indicators, T07 A18 or S11 8
  double beta[3];
  beta[0] = (13./12.)*pow(x1 - 2.*x2 + x3, 2) +
            (1./4.)*pow(x1 - 4.*x2 + 3.*x3, 2);
  beta[1] = (13./12.)*pow(x2 - 2.*x3 + x4, 2) +
            (1./4.)*pow(x4 - x2, 2);
  beta[2] = (13./12.)*pow(x3 - 2.*x4 + x5, 2) +
            (1./4.)*pow(x5 - 4.*x4 + 3.*x3, 2);

  // Nonlinear weights S11 9
  double den, wtr[3], Wr, wr[3], wtl[3], Wl, wl[3], eps;
  eps=1.e-26;

  den = eps + beta[0]; den *= den; wtr[0] = (1./16.)/den;
  den = eps + beta[1]; den *= den; wtr[1] = (5./8. )/den;
  den = eps + beta[2]; den *= den; wtr[2] = (5./16.)/den;
  Wr = wtr[0] + wtr[1] + wtr[2];
  wr[0] = wtr[0]/Wr ;
  wr[1] = wtr[1]/Wr ;
  wr[2] = wtr[2]/Wr ;

  den = eps + beta[2]; den *= den; wtl[0] = (1./16.)/den;
  den = eps + beta[1]; den *= den; wtl[1] = (5./8. )/den;
  den = eps + beta[0]; den *= den; wtl[2] = (5./16.)/den;
  Wl = wtl[0] + wtl[1] + wtl[2];
  wl[0] = wtl[0]/Wl;
  wl[1] = wtl[1]/Wl;
  wl[2] = wtl[2]/Wl;

  *lout = vl[0]*wl[0] + vl[1]*wl[1] + vl[2]*wl[2];
  *rout = vr[0]*wr[0] + vr[1]*wr[1] + vr[2]*wr[2];
}

void reconstruct_lr_lin(double Ptmp[NMAX+2*NG][NVAR], int N,
  double P_l[NMAX+2*NG][NVAR], double P_r[NMAX+2*NG][NVAR])
{
  double dqtmp[NMAX + 2*NG][NVAR];

  ISLOOP(-1, N) {
    PLOOP {
      dqtmp[i][ip] = 0.;
    }
  }

  // Calculate slopes
  ISLOOP(-1, N) {
    PLOOP {
      dqtmp[i][ip] = slope_lim(Ptmp[i-1][ip], Ptmp[i][ip], Ptmp[i+1][ip]);
    }
  }

  // Reconstruct left
  ISLOOP(0, N) {
    PLOOP {
      P_l[i][ip] = Ptmp[i][ip] - 0.5*dqtmp[i][ip];
    }
  }

  // Reconstruct right
  ISLOOP(-1, N-1) {
    PLOOP {
      P_r[i][ip] = Ptmp[i][ip] + 0.5*dqtmp[i][ip];
    }
  }
}

void reconstruct_lr_weno(double Ptmp[NMAX+2*NG][NVAR], int N,
  double P_l[NMAX+2*NG][NVAR], double P_r[NMAX+2*NG][NVAR])
{
  ISLOOP(-1,N) {
    PLOOP {
      weno(Ptmp[i-2][ip],
           Ptmp[i-1][ip],
           Ptmp[i][ip],
           Ptmp[i+1][ip],
           Ptmp[i+2][ip],
           &P_l[i][ip],
           &P_r[i][ip]);
    }
  }
}

void reconstruct(double Ptmp[NMAX+2*NG][NVAR], int N,
  double P_l[NMAX+2*NG][NVAR], double P_r[NMAX+2*NG][NVAR])
{
  #if RECONSTRUCTION == LINEAR
    reconstruct_lr_lin(Ptmp, N, P_l, P_r);
  #elif RECONSTRUCTION == WENO
    reconstruct_lr_weno(Ptmp, N, P_l, P_r);
  #else
    fprintf(stderr, "Reconstruction algorithm not specified!\n");
    exit(-1);
  #endif
}

