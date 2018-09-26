/******************************************************************************
 *                                                                            *
 * SCATTERING.C                                                               *
 *                                                                            *
 * RELATIVISTIC SCATTERING KERNEL                                             *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
void sample_electron(double Ktetrad[NDIM], double Pelectron[NDIM],
  double Thetae);
void sample_beta(double Thetae, double *gamma_e, double *beta_e);
double sample_y(double Thetae);
double sample_mu(double beta_e);
void sample_scattered_photon(double k[NDIM], double p[NDIM], double kp[NDIM]);
void boost(double v[NDIM], double u[NDIM], double vp[NDIM]);
double sample_thomson();
double sample_klein_nishina(double k0);
double klein_nishina(double a, double ap);

int scatter_superphoton(grid_prim_type P, struct of_photon *ph, double X[NDIM],
  double Kcov[NDIM], double Kcon[NDIM])
{
  double Pelectron[NDIM], gcov[NDIM][NDIM], gcon[NDIM][NDIM];
  double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];
  
  double Ktetrad[NDIM], Ktetrad_scatt[NDIM];

  int i, j, k;
  Xtoijk(X, &i, &j, &k);

  set_gcov(X, gcov);
  gcon_func(gcov, gcon);

  normalize_null(gcov, Kcon);
  normalize_null_cov(gcon, Kcov);

  // Quality control
  if (Kcon[0] < 0. || isnan(Kcon[0]) || Kcov[0] > 0. || isnan(Kcov[0])) {
    //fprintf(stderr, "Normalization problem in scattering\n");
    //printf("Kcon[] = %e %e %e %e\n", Kcon[0], Kcon[1], Kcon[2], Kcon[3]);
    //printf("Kcov[] = %e %e %e %e\n", Kcov[0], Kcov[1], Kcov[2], Kcov[3]);
    //printf("X[] = %e %e %e %e\n", X[0], X[1], X[2], X[3]);
    
    return 0;
  }

  double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
  struct of_microphysics m;
  get_fluid_zone(i, j, k, P, &m, Ucon, Ucov, Bcon, Bcov);

  if (m.Thetae < 10.*SMALL) {
    return 0;
  }

  make_tetrad(i, j, k, Ucon, Bcon, gcov, Econ, Ecov);

  coord_to_tetrad(Ecov, Kcon, Ktetrad);

  for (int mu = 0; mu < NDIM; mu++) {
    if (isnan(Ktetrad[mu]) || isinf(Ktetrad[mu])) {
      printf("Ktetrad[%i] = %e!\n", mu, Ktetrad[mu]);
      printf("Ucon[%i] = %e\n", mu, Ucon[mu]);
      printf("Ucov[%i] = %e\n", mu, Ucov[mu]);
      printf("Bcon[%i] = %e\n", mu, Bcon[mu]);
      printf("Bcov[%i] = %e\n", mu, Bcov[mu]);
    }
  }

  sample_electron(Ktetrad, Pelectron, m.Thetae);

  for (int mu = 0; mu < NDIM; mu++) {
    if (isnan(Pelectron[mu]) || isinf(Pelectron[mu])) {
      printf("m.Thetae = %e m.Ne = %e m.B = %e\n", m.Thetae, m.Ne, m.B);
      printf("Pelectron[%i] = %e!\n", mu, Pelectron[mu]);
      printf("K = %e %e %e %e\n", Ktetrad[0], Ktetrad[1], Ktetrad[2], Ktetrad[3]);
      printf("k.k = %e\n", -Ktetrad[0]*Ktetrad[0] + Ktetrad[1]*Ktetrad[1] + Ktetrad[2]*Ktetrad[2] + Ktetrad[3]*Ktetrad[3]);
    }
  }

  sample_scattered_photon(Ktetrad, Pelectron, Ktetrad_scatt);

  for (int mu = 0; mu < NDIM; mu++) {
    if (isnan(Ktetrad_scatt[mu]) || isinf(Ktetrad_scatt[mu])) {
      printf("Ktetrad_scatt[%i] = %e!\n", mu, Ktetrad_scatt[mu]);
    }
  }

  // Check NAN after each of these

  tetrad_to_coord(Econ, Ktetrad_scatt, ph->Kcon[2]);

  for (int mu = 0; mu < NDIM; mu++) {
    if (isnan(ph->Kcon[2][mu]) || isinf(ph->Kcon[2][mu])) {
      printf("ph->Kcon[2][%i] = %e!\n", mu, ph->Kcon[2][mu]);
    }
  }
  
  // Ensure scattered superphoton is sane
  normalize_null(gcov, ph->Kcon[2]);

  for (int mu = 0; mu < NDIM; mu++) {
    if (isnan(ph->Kcon[2][mu]) || isinf(ph->Kcon[2][mu])) {
      printf("after norm ph->Kcon[2][%i] = %e!\n", mu, ph->Kcon[2][mu]);
    }
  }

  lower(ph->Kcon[2], gcov, ph->Kcov[2]);

  for (int mu = 0; mu < NDIM; mu++) {
    if (isnan(ph->Kcov[2][mu]) || isinf(ph->Kcov[2][mu])) {
      printf("after norm ph->Kcov[2][%i] = %e!\n", mu, ph->Kcov[2][mu]);
    }
  }

  return 1;
}

// Procedure from Canfield et al. 1987
void sample_electron(double k[NDIM], double p[NDIM], double Thetae)
{
  double beta_e, mu, phi, cphi, sphi, gamma_e, sigma_KN;
  double K, sth, cth, x1, n0dotv0, v0, v1;
  double n0x, n0y, n0z;
  double v0x, v0y, v0z;
  double v1x, v1y, v1z;
  double v2x, v2y, v2z;
  int sample_cnt = 0;

  do {
    sample_beta(Thetae, &gamma_e, &beta_e);
    mu = sample_mu(beta_e);

    // Sometimes |mu| > 1 from roundoff error. Fix it
    if (mu > 1.) mu = 1.;
    else if (mu < -1.) mu = -1;

    // Frequency in electron rest frame
    K = gamma_e*(1. - beta_e*mu)*k[0];

    // Avoid numerical craziness for small K
    if (K < 1.e-3) {
      sigma_KN = 1. - 2.*K + 5.2*K*K - 13.3*K*K*K + 1144*K*K*K*K/35.;
    } else{
      // Klein-Nishina cross-section / Thomson
      sigma_KN = (3./(4.*K*K))*(2. + K*K*(1. + K)/
        ((1. + 2.*K)*(1. + 2.*K)) + (K*K - 2.*K - 2.)/(2.*K)*log(1. + 2.*K));
    }

    x1 = get_rand();

    sample_cnt++;

    if(sample_cnt > 1000000) {
      fprintf(stderr,"in sample_electron mu, gamma_e, K, sigma_KN, x1: %g %g %g %g %g %g\n",
        Thetae, mu, gamma_e, K, sigma_KN, x1);

      // Kluge to prevent stalling for large values of \Theta_e
      Thetae *= 0.5 ;
      sample_cnt = 0 ;
    }
  } while (x1 >= sigma_KN);

  // First unit vector for coordinate system
  v0x = k[1];
  v0y = k[2];
  v0z = k[3];
  v0 = sqrt(v0x*v0x + v0y*v0y + v0z*v0z);
  v0x /= v0;
  v0y /= v0;
  v0z /= v0;

  // Pick zero-angle for coordinate system
  get_ran_dir_3d(&n0x, &n0y, &n0z);
  n0dotv0 = v0x*n0x + v0y*n0y + v0z*n0z;

  // Second unit vector
  v1x = n0x - (n0dotv0)*v0x;
  v1y = n0y - (n0dotv0)*v0y;
  v1z = n0z - (n0dotv0)*v0z;

  // Normalize
  v1 = sqrt(v1x*v1x + v1y*v1y + v1z*v1z);
  v1x /= v1;
  v1y /= v1;
  v1z /= v1;

  // Find one more unit vector using cross product; automatically normalized
  v2x = v0y*v1z - v0z*v1y;
  v2y = v0z*v1x - v0x*v1z;
  v2z = v0x*v1y - v0y*v1x;

  // Resolve new momentum vector along unit vectors and create a four-vector p
  phi = get_rand()*2.*M_PI; // uniform orientation
  sphi = sin(phi);
  cphi = cos(phi);
  cth = mu;
  sth = sqrt(1. - mu*mu);

  p[0] = gamma_e;
  p[1] = gamma_e*beta_e*(cth*v0x + sth*(cphi*v1x + sphi*v2x));
  p[2] = gamma_e*beta_e*(cth*v0y + sth*(cphi*v1y + sphi*v2y));
  p[3] = gamma_e*beta_e*(cth*v0z + sth*(cphi*v1z + sphi*v2z));

  if (beta_e < 0) {
    fprintf(stderr, "betae error: %g %g %g %g\n", p[0], p[1], p[2], p[3]);
  }
}

void sample_beta(double Thetae, double *gamma_e, double *beta_e)
{
  double y = sample_y(Thetae);

  *gamma_e = y*y*Thetae + 1.;
  *beta_e = sqrt(1. - 1./((*gamma_e)*(*gamma_e)));
}

double sample_y(double Thetae)
{
  double S_3, pi_3, pi_4, pi_5, pi_6, y, x1, x2, x, prob, num, den;

  pi_3 = sqrt(M_PI) / 4.;
  pi_4 = sqrt(0.5 * Thetae) / 2.;
  pi_5 = 3. * sqrt(M_PI) * Thetae / 8.;
  pi_6 = Thetae * sqrt(0.5 * Thetae);

  S_3 = pi_3 + pi_4 + pi_5 + pi_6;

  pi_3 /= S_3;
  pi_4 /= S_3;
  pi_5 /= S_3;
  pi_6 /= S_3;

  int max_samp = 100000;
  int n = 0;
  do {
    n++;
    x1 = get_rand();

    if (x1 < pi_3) {
      x = get_chisq(3);
    } else if (x1 < pi_3 + pi_4) {
      x = get_chisq(4);
    } else if (x1 < pi_3 + pi_4 + pi_5) {
      x = get_chisq(5);
    } else {
      x = get_chisq(6);
    }

    // Translate between Canfield et al. and standard chisq distribution
    y = sqrt(x/2);

    x2 = get_rand();
    num = sqrt(1. + 0.5*Thetae*y*y);
    den = 1. + y*sqrt(0.5*Thetae);

    prob = num/den;

  } while (x2 >= prob && n < max_samp);

  if (n >= max_samp) {
    fprintf(stderr, "FAILED TO SAMPLE Y! Thetae = %e\n", Thetae);
    exit(-1);
  }

  return y;
}

double sample_mu(double beta_e)
{
  double mu, x1;

  x1 = get_rand();
  mu = (1. - sqrt(1. + 2.*beta_e + beta_e*beta_e - 4.*beta_e*x1))/beta_e;
  return mu;
}

void sample_scattered_photon(double k[NDIM], double p[NDIM], double kp[NDIM])
{
  double ke[4], kpe[4];
  double k0p;
  double n0x, n0y, n0z, n0dotv0, v0x, v0y, v0z, v1x, v1y, v1z, v2x,
      v2y, v2z, v1, dir1, dir2, dir3;
  double cth, sth, phi, cphi, sphi;

  boost(k, p, ke);

  k0p = sample_klein_nishina(ke[0]);
  cth = 1. - 1./k0p + 1./ke[0];
  sth = sqrt(fabs(1. - cth*cth));

  // Unit vector 1 for scattering coordinate system is oriented along initial
  // photon wavevector

  // Explicitly compute kemag instead of using ke[0] to ensure that photon is
  // created normalized and doesn't inherit the light cone errors from the
  // original photon
  double kemag = sqrt(ke[1]*ke[1] + ke[2]*ke[2] + ke[3]*ke[3]);
  v0x = ke[1]/kemag;
  v0y = ke[2]/kemag;
  v0z = ke[3]/kemag;

  // Randomly pick zero-angle for scattering coordinate system.
  get_ran_dir_3d(&n0x, &n0y, &n0z);
  n0dotv0 = v0x*n0x + v0y*n0y + v0z*n0z;

  // Unit vector 2
  v1x = n0x - (n0dotv0)*v0x;
  v1y = n0y - (n0dotv0)*v0y;
  v1z = n0z - (n0dotv0)*v0z;
  v1 = sqrt(v1x*v1x + v1y*v1y + v1z*v1z);
  v1x /= v1;
  v1y /= v1;
  v1z /= v1;

  // Find one more unit vector using cross product; automatically normalized
  v2x = v0y*v1z - v0z*v1y;
  v2y = v0z*v1x - v0x*v1z;
  v2z = v0x*v1y - v0y*v1x;

  // Resolve new momentum vector along unit vectors
  // Create a four-vector p
  // Solve for orientation of scattered photon

  // Find phi for new photon
  phi = 2.*M_PI*get_rand();
  sphi = sin(phi);
  cphi = cos(phi);

  p[1] *= -1.;
  p[2] *= -1.;
  p[3] *= -1.;

  dir1 = cth*v0x + sth*(cphi*v1x + sphi*v2x);
  dir2 = cth*v0y + sth*(cphi*v1y + sphi*v2y);
  dir3 = cth*v0z + sth*(cphi*v1z + sphi*v2z);

  kpe[0] = k0p;
  kpe[1] = k0p*dir1;
  kpe[2] = k0p*dir2;
  kpe[3] = k0p*dir3;

  // Transform k back to lab frame
  boost(kpe, p, kp);

  if(kp[0] < 0 || isnan(kp[0])) {
    fprintf(stderr,"in sample_scattered_photon:\n") ;
    fprintf(stderr,"kp[0], kpe[0]: %g %g\n",kp[0],kpe[0]) ;
    fprintf(stderr,"kpe: %g %g %g %g\n",kpe[0],kpe[1],kpe[2],kpe[3]) ;
    fprintf(stderr,"k:  %g %g %g %g\n",k[0],k[1],k[2],k[3]) ;
    fprintf(stderr,"p:   %g %g %g %g\n",p[0],p[1],p[2],p[3]) ;
    fprintf(stderr,"kp:  %g %g %g %g\n",kp[0],kp[1],kp[2],kp[3]) ;
  }
}

void boost(double v[NDIM], double u[NDIM], double vp[NDIM])
{
  double g, V, n1, n2, n3, gm1;

  g = u[0];
  V = sqrt(fabs(1. - 1./(g*g)));
  n1 = u[1]/(g*V + SMALL);
  n2 = u[2]/(g*V + SMALL);
  n3 = u[3]/(g*V + SMALL);
  gm1 = g - 1.;

  // Lorentz boost into frame u from lab frame
  vp[0] =  u[0]*v[0] - (          u[1])*v[1] - (          u[2])*v[2] -
    (          u[3])*v[3];
  vp[1] = -u[1]*v[0] + (1. + n1*n1*gm1)*v[1] + (     n1*n2*gm1)*v[2] +
    (     n1*n3*gm1)*v[3];
  vp[2] = -u[2]*v[0] + (     n2*n1*gm1)*v[1] + (1. + n2*n2*gm1)*v[2] +
    (     n2*n3*gm1)*v[3];
  vp[3] = -u[3]*v[0] + (     n3*n1*gm1)*v[1] + (     n3*n2*gm1)*v[2] +
    (1. + n3*n3*gm1)*v[3];
}

double sample_thomson()
{
  double x1, x2;

  do {
    x1 = 2.*get_rand() - 1.;
    x2 = (3./4.)*get_rand();
  } while (x2 >= (3./8.)*(1. + x1*x1)) ;

  return x1;
}

double sample_klein_nishina(double k0)
{
  double k0pmin, k0pmax, k0p_tent, x1;
  int n = 0;

  // A low efficiency sampling algorithm, particularly for large k0. Limiting
  // efficiency is log(2 k0)/(2 k0)
  k0pmin = k0/(1. + 2.*k0); // at theta = Pi
  k0pmax = k0;              // at theta = 0

  do {
    // Tentative value
    k0p_tent = k0pmin + (k0pmax - k0pmin)*get_rand();

    // Rejection sample in box of height = kn(kmin)
    x1 = 2.*(1. + 2.*k0 + 2.*k0*k0)/(k0*k0*(1. + 2.*k0));
    x1 *= get_rand();

    n++;
  } while (x1 >= klein_nishina(k0, k0p_tent));

  return k0p_tent;
}

double klein_nishina(double a, double ap)
{
  double ch = 1. + 1./a - 1./ap;
  double kn = (a/ap + ap/a - 1. + ch*ch)/(a*a);

  return kn;
}

#define NW 220
#define NT 90
#define MINW (1.e-16)
#define MAXW (1.e10)
#define MINT (0.001)
#define MAXT (1.e11)
#define HOTCROSS "hotcross.dat"

double total_compton_cross_num(double w, double Thetae);
double dNdgammae(double thetae, double gammae);
double boostcross(double w, double mue, double gammae);
double hc_klein_nishina(double we);

double table[NW+1][NT+1];
double dlw, dlT, lminw, lmint;

void init_hotcross()
{
  int nread;
  double lw, lT;
  FILE *fp;

  dlw = log10(MAXW/MINW)/NW;
  dlT = log10(MAXT/MINT)/NT;
  lminw = log10(MINW);
  lmint = log10(MINT);

  // Create file if needed using IO proc
  if (mpi_io_proc()) {
    fp = fopen(HOTCROSS, "r");
    if (fp == NULL) {
      fprintf(stdout, "Making lookup table for Compton cross section...\n");
      #pragma omp parallel for collapse(2)
      for (int i = 0; i <= NW; i++) {
        for (int j = 0; j <= NT; j++) {
          lw = lminw + i*dlw;
          lT = lmint + j*dlT;
          table[i][j] = log10(total_compton_cross_num(pow(10.,lw),
            pow(10.,lT)));
          if(isnan(table[i][j])) {
            fprintf(stderr, "NAN for Compton cross section: %d %d %g %g\n",
              i, j, lw, lT);
            exit(0);
          }
        }
      }
      fprintf(stdout, "Lookup table created.\n\n");

      fprintf(stdout, "Writing lookup table to file...\n");
      fp = fopen(HOTCROSS, "w");
      if (fp == NULL) {
        fprintf(stderr, "Couldn't write to file\n");
        exit(0);
      }
      for(int i = 0; i <= NW; i++) {
        for(int j = 0; j <= NT; j++) {
          lw = lminw + i*dlw;
          lT = lmint + j*dlT;
          fprintf(fp, "%d %d %g %g %15.10g\n",i,j,lw,lT,table[i][j]);
        }
      }
      fprintf(stderr, "Lookup table written.\n\n");
    } // fp == NULL
    fclose(fp);
  } // mpi_io_proc()

  mpi_barrier();

  // Read lookup table with every MPI processor
  fp = fopen(HOTCROSS, "r");
  if (fp == NULL) {
    fprintf(stderr, "rank %i: file %s not found.\n", mpi_myrank(), HOTCROSS);
    exit(-1);
  }
  for(int i = 0; i <= NW; i++) {
    for(int j = 0; j <= NT; j++) {
      nread = fscanf(fp, "%*d %*d %*f %*f %lf\n", &table[i][j]);
      if(isnan(table[i][j]) || nread != 1) {
        fprintf(stderr, "Error on table read: %d %d\n", i, j) ;
        exit(0) ;
      }
    }
  }
  fclose(fp);
}

double total_compton_cross_lkup(double w, double thetae)
{
  int i, j;
  double lw, lT, di, dj, lcross;

  // Cold/low-energy: Use Thomson cross section
  if (w*thetae < 1.e-6) return(THOMSON);

  // Cold, but possible high-energy photon: use Klein-Nishina
  if (thetae < MINT) return (hc_klein_nishina(w)*THOMSON);

  // In-bounds for table
  if ((w > MINW && w < MAXW) && (thetae > MINT && thetae < MAXT)) {
    lw = log10(w);
    lT = log10(thetae);
    i = (int)((lw - lminw)/dlw);
    j = (int)((lT - lmint)/dlT);
    di = (lw - lminw)/dlw - i;
    dj = (lT - lmint)/dlT - j;

    lcross = (1. - di)*(1. - dj)*table[i][j] + di*(1. - dj)*table[i+1][j] +
             (1. - di)*dj*table[i][j+1]      + di*dj*table[i+1][j+1];

    if(isnan(lcross)) {
      fprintf(stderr, "NAN cross section: %g %g %d %d %g %g\n",
        lw, lT, i, j, di, dj);
    }

    return pow(10.,lcross);
  }

  fprintf(stderr, "Cross section out of bounds: %g [%g,%g] %g [%g,%g]\n",
    w, MINW, MAXW, thetae, MINT, MAXT);

  return total_compton_cross_num(w, thetae);
}

#define MAXGAMMA (12.)
#define DMUE     (0.05)
#define DGAMMAE  (0.05)

double total_compton_cross_num(double w, double thetae)
{
  double dmue, dgammae, mue, gammae, f, cross;

  if(isnan(w)) {
    fprintf(stderr, "NAN Compton cross section: %g %g\n", w, thetae);
    return 0.;
  }

  // Check for easy limits
  if (thetae < MINT && w < MINW) return THOMSON;
  if (thetae < MINT) return hc_klein_nishina(w)*THOMSON;

  dmue = DMUE;
  dgammae = thetae*DGAMMAE;

  // Integrate over mu_e and gamma_e, where mu_e is the cosine of the angle
  // between K and U_e, and the angle k is assumed to lie, wlog, along the z
  // z axis
  cross = 0.;
  for (mue = -1. + 0.5*dmue; mue < 1.; mue += dmue)
  for (gammae = 1. + 0.5*dgammae;
       gammae < 1. + MAXGAMMA*thetae;
       gammae += dgammae) {
    f = 0.5*dNdgammae(thetae, gammae);

    cross += dmue*dgammae*boostcross(w, mue, gammae)*f;

    if(isnan(cross)) {
      fprintf(stderr, "NAN cross section: %g %g %g %g %g %g\n",
        w, thetae, mue, gammae, dNdgammae(thetae, gammae),
        boostcross(w, mue, gammae));
    }
  }

  return cross*THOMSON;
}

// Normalized (per unit proper electron number density) electron distribution
double dNdgammae(double thetae, double gammae)
{
  double K2f;

  if (thetae > 1.e-2) {
    K2f = gsl_sf_bessel_Kn(2, 1./thetae)*exp(1./thetae);
  } else {
    K2f = sqrt(M_PI*thetae/2.) + 15./8.*sqrt(M_PI/2.)*pow(thetae,1.5) +
          105./128.*sqrt(M_PI/2.)*pow(thetae,2.5) -
          315./1024.*sqrt(M_PI/2.)*pow(thetae,3.5);
  }

  return (gammae*sqrt(gammae*gammae - 1.)/(thetae*K2f))*exp(-(gammae-1.)/thetae);
}

double boostcross(double w, double mue, double gammae)
{
  double we, boostcross, v;

  // Energy in electron rest frame
  v = sqrt(gammae*gammae - 1.)/gammae;
  we = w*gammae*(1. - mue*v);

  boostcross = hc_klein_nishina(we)*(1. - mue*v);

  if(boostcross > 2) {
    fprintf(stderr, "w, mue, gammae: %g %g %g\n", w, mue, gammae);
    fprintf(stderr, "v, we, boostcross: %g %g %g\n", v, we, boostcross);
    fprintf(stderr, "kn: %g %g %g\n", v, we, boostcross);
  }

  if(isnan(boostcross)) {
    fprintf(stderr, "isnan: %g %g %g\n", w, mue, gammae);
    return 0.;
  }

  return boostcross;
}

double hc_klein_nishina(double we)
{
  double sigma;

  if (we < 1.e-3) return(1. - 2.*we);

  sigma = (3./4.)*(
    2./(we*we) +
    (1./(2.*we) -
    (1. + we)/(we*we*we))*log(1. + 2.*we) +
    (1. + we)/((1. + 2.*we)*(1. + 2.*we))
    );

  return sigma;
}
#endif // RADIATION

