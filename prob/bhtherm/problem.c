/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR THERMAL SCHWARZSCHILD ATMOSPHERE                    *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#define INIT_FILE ("init.txt")

double *r_init, *P_cgs_init, *rho_cgs_init;

void set_problem_params() {
}

void init_zone_thermal(int i, int j, int k, double dndlnu[NU_BINS_EMISS+1]);

void init_prob()
{
  // Read in 1D array of fluid data and r from semianalytic solution
  FILE *infil = fopen(INIT_FILE, "r");
  char buf[STRLEN];
  fgets(buf, STRLEN, infil);
  int nsup;
  sscanf(buf, "%d", &nsup);
  r_init = safe_malloc(nsup*sizeof(double));
  P_cgs_init = safe_malloc(nsup*sizeof(double));
  rho_cgs_init = safe_malloc(nsup*sizeof(double));
  int n = 0;
  while (fgets(buf, STRLEN, infil) != NULL) {
    sscanf(buf, "%lf %lf %lf\n", &r_init[n], &P_cgs_init[n], &rho_cgs_init[n]);
    n++;
  }
  fclose(infil);

  double X[NDIM], r, th;
  ZLOOP {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);

    // Get solution
    n = 0;
    while (!(r_init[n+1] > r && r_init[n] < r)) {
      n++;
    }
    double P_cgs = (r - r_init[n])*(P_cgs_init[n+1] - P_cgs_init[n])/(r_init[n+1] - r_init[n]) + P_cgs_init[n];
    double rho_cgs = (r - r_init[n])*(rho_cgs_init[n+1] - rho_cgs_init[n])/(r_init[n+1] - r_init[n]) + rho_cgs_init[n];
  
    // Get primitive 3-velocity
    struct of_geom geom = ggeom[i][j][CENT];
    double alpha = geom.alpha;
    double beta[NDIM];
    for (int l = 1; l < NDIM; l++) {
      beta[l] = geom.gcon[0][l]*alpha*alpha;
    }
    
    P[i][j][k][RHO] = rho_cgs/RHO_unit;
    P[i][j][k][UU] = P_cgs/(gam - 1.)/U_unit;
    P[i][j][k][U1] = beta[1]/(sqrt(alpha*alpha-beta[1]*beta[1]*geom.gcov[1][1]));
    P[i][j][k][U2] = beta[2]/(sqrt(alpha*alpha-beta[2]*beta[2]*geom.gcov[2][2]));
    P[i][j][k][U3] = beta[3]/(sqrt(alpha*alpha-beta[3]*beta[3]*geom.gcov[3][3]));
    //P[i][j][k][U2] = 0.;
    //P[i][j][k][U3] = 0.;
    P[i][j][k][B1] = 1.e-10;
    P[i][j][k][B2] = 0.;
    P[i][j][k][B3] = 0.;
   // printf("T: %e\n", 0.5*(gam-1)*P[i][j][k][UU]/P[i][j][k][RHO]*MP*CL*CL/KBOL);
  }
  //exit(-1);

  fixup(P);
  bound_prim(P);
  bound_prim(Ph);

  int photons_per_zone = 100;

  #pragma omp parallel for collapse(2)
  ZLOOP {
    //double Rcoord[NDIM][NDIM]={0}, Rfluid[NDIM][NDIM]={0}, 
    double Ranalytic[NDIM][NDIM]={0}, Rcoordanalytic[NDIM][NDIM]={0};
    
    double dndlnu[NU_BINS_EMISS+1];
    init_zone_thermal(i, j, k, dndlnu);
    double lnumin = log(numin_emiss);
    double lnumax = log(numax_emiss);
    double dlnu = (lnumax - lnumin)/NU_BINS_EMISS;
    //struct of_photon *ph = photon_lists[omp_get_thread_num()];
    struct of_photon *ph = NULL, *tail = NULL;

    double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
    struct of_microphysics m;
    get_fluid_zone(i, j, k, P, &m, Ucon, Ucov, Bcon, Bcov);
    double Er = AR*pow(m.Thetae*ME*CL*CL/KBOL,4)/U_unit;
    //printf("Er = %e UU = %e Tr = %e\n", Er, P[i][j][k][UU], m.Thetae*ME*CL*CL/KBOL);
    Ranalytic[0][0] = Er;
    Ranalytic[1][1] = Er/3.;
    Ranalytic[2][2] = Er/3.;
    Ranalytic[3][3] = Er/3.;

    double Bhatcon[NDIM], Econ[NDIM][NDIM], Ecov[NDIM][NDIM];
    for (int l = 0; l < NDIM; l++) Bhatcon[l] = Bcon[l]/m.B;
    //Bhatcon[1] = 1.;
    make_tetrad(i, j, k, Ucon, Bhatcon, ggeom[i][j][CENT].gcov, Econ, Ecov);

    //DLOOP2 printf("Econ[%i][%i] = %e\n", mu, nu, Econ[mu][nu]);
    //DLOOP2 printf("Ecov[%i][%i] = %e\n", mu, nu, Ecov[mu][nu]);

    double X[NDIM], r, th;
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);

    for (int mu = 0; mu < NDIM; mu++) {
      for (int nu = 0; nu < NDIM; nu++) {
        for (int lam = 0; lam < NDIM; lam++) {
          for (int kap = 0; kap < NDIM; kap++) {
            Rcoordanalytic[mu][nu] += Econ[mu][lam]*Econ[nu][kap]*Ranalytic[lam][kap];
          }
        }
      }
    }
    //printf("Er_c = %e, rat = %e\n", Rcoordanalytic[0][0], Rcoordanalytic[0][0]/Ranalytic[0][0]);

    double Etot_tet = 0., Etot_coord = 0.;

    struct of_photon *tmp[2];
    double sth[2], cth[2], sphi[2], cphi[2], phi, K_tetrad[NDIM];
    for (int m = 0; m < photons_per_zone; m++) {
      tmp[0] = safe_malloc(sizeof(struct of_photon));
      tmp[1] = safe_malloc(sizeof(struct of_photon));
      tmp[0]->next = tmp[1];
      tmp[1]->next = ph;
      ph = tmp[0];
      if (m == 0) tail = tmp[1];
    
      double nu;
      do {
        nu = exp(get_rand()*(lnumax - lnumin) + lnumin);
      } while (get_rand() > linear_interp_log(nu, dndlnu, lnumin, dlnu));
      
      // Equal E per superphoton
      //nu = 1.e20;
      double w = 1./nu;

      // Equal superphotons per d log nu
      //struct of_microphysics mm;
      //get_fluid_zone(i, j, k, P, &mm, Ucon, Ucov, Bcon, Bcov);
      //double numax = 5.879e10*mm.Thetae*ME*CL*CL/KBOL;
      //double Bnu = Bnu_inv(nu, &mm)*pow(nu,3);
      //double Bnu_max = Bnu_inv(numax, &mm)*pow(numax,3);
      //double w = Bnu/Bnu_max;

      cth[0] = 2.*get_rand() - 1.;
      sth[0] = sqrt(1. - cth[0]*cth[0]);
      phi = 2.*M_PI*get_rand();
      cphi[0] = cos(phi);
      sphi[0] = sin(phi);

      cth[1] = -cth[0];
      sth[1] = sth[0];
      cphi[1] = -cphi[0];
      sphi[1] = -sphi[0];

      double E = nu*HPL/(ME*CL*CL);

      for (int n = 0; n < 2; n++) {
        memset(tmp[n]->X, 0, NSUP*NDIM*sizeof(double));
        memset(tmp[n]->Kcov, 0, NSUP*NDIM*sizeof(double));
        memset(tmp[n]->Kcon, 0, NSUP*NDIM*sizeof(double));

        tmp[n]->X[0][0] = -1.;
        tmp[n]->X[1][0] = -1.;

        tmp[n]->X[2][0] = 0.;
        coord(i, j, k, CENT, tmp[n]->X[2]);

        K_tetrad[0] = -E;
        K_tetrad[1] = E*cth[n];
        K_tetrad[2] = E*cphi[n]*sth[n];
        K_tetrad[3] = E*sphi[n]*sth[n];

        tetrad_to_coord(Ecov, K_tetrad, tmp[n]->Kcov[2]);

        K_tetrad[0] *= -1.;
        tetrad_to_coord(Econ, K_tetrad, tmp[n]->Kcon[2]);

        tmp[n]->w = 0.5*w*tmp[n]->Kcon[2][0]/K_tetrad[0]*ggeom[i][j][CENT].g;
        //tmp[n]->w = w;//0.5*w/(K_tetrad[0]/tmp[n]->Kcon[2][0]);
        //tmp[n]->w = 0.5*w*(K_tetrad[0]/(tmp[n]->Kcon[2][0]*ggeom[i][j][CENT].g));
        //tmp[n]->w = 0.5*w/tmp[n]->Kcon[2][0]*K_tetrad[0]/ggeom[i][j][CENT].g;

    
        double d3V = dx[1]*dx[2]*dx[3];
        //double d3v_fluid = ggeom[i][j][CENT].g*d3V*tmp[n]->Kcon[2][0]/K_tetrad[0];
        //Etot_tet += K_tetrad[0]*w/2*kphys_to_num/d3v_fluid;
        //Etot_tet += K_tetrad[0]*tmp[n]->w*kphys_to_num/d3v_fluid;
        Etot_tet += K_tetrad[0]*K_tetrad[0]*tmp[n]->w*kphys_to_num/(tmp[n]->Kcon[2][0]*ggeom[i][j][CENT].g*dx[1]*dx[2]*dx[3]);
        Etot_coord += tmp[n]->Kcon[2][0]*tmp[n]->w*kphys_to_num/(ggeom[i][j][CENT].g*d3V);

        tmp[n]->nscatt = 0;
        tmp[n]->origin[0] = 0;
        tmp[n]->origin[1] = i;
        tmp[n]->origin[2] = j;
        tmp[n]->origin[3] = k;
        tmp[n]->t0 = 0.;
      }
    }
    //printf("Etot_tet = %e Etot_coord = %e rat = %e\n", Etot_tet, Etot_coord, Etot_coord/Etot_tet);
    
    double Rs_fluid[NDIM][NDIM] = {0};
    double Rs_coord[NDIM][NDIM] = {0};
    double Rs_fluid_transform[NDIM][NDIM] = {0};

    // Normalize superphoton weights
    double Es_final = 0.;
    struct of_photon *phtmp = ph;
    while (phtmp != NULL) {
      //phtmp->w *= Rcoordanalytic[0][0]/Etot_coord;
      phtmp->w *= Ranalytic[0][0]/Etot_tet;
      Es_final += phtmp->Kcon[2][0]*phtmp->w*kphys_to_num/(ggeom[i][j][CENT].g*dx[1]*dx[2]*dx[3]);
      
      double d3V = dx[1]*dx[2]*dx[3];
      double K_tetrad[NDIM] = {0};
      coord_to_tetrad(Ecov, phtmp->Kcon[2], K_tetrad);
      double d3v_coord = ggeom[i][j][CENT].g*d3V;
      double d3v_fluid = ggeom[i][j][CENT].g*d3V*phtmp->Kcon[2][0]/K_tetrad[0];
      //d3v_fluid = dx[1]*dx[2]*dx[3];
      for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
          Rs_coord[mu][nu] += phtmp->Kcon[2][mu]*phtmp->Kcon[2][nu]/phtmp->Kcon[2][0]*phtmp->w*kphys_to_num/d3v_coord;
          Rs_fluid[mu][nu] += K_tetrad[mu]*K_tetrad[nu]/K_tetrad[0]*phtmp->w*kphys_to_num/d3v_fluid;
        }
      } 
      
      phtmp = phtmp->next;
    }
    for (int mu = 0; mu < NDIM; mu++) {
      for (int nu = 0; nu < NDIM; nu++) {
        for (int lam = 0; lam < NDIM; lam++) {
          for (int kap = 0; kap < NDIM; kap++) {
            Rs_fluid_transform[mu][nu] += Ecov[mu][lam]*Ecov[nu][kap]*Rs_coord[lam][kap];
          }
        }
      }
    }
    /*
    printf("Es_final = %e\n", Es_final);
    for (int mu = 0; mu < NDIM; mu++) {
      for (int nu = 0; nu < NDIM; nu++) {
        printf("Rs_fluid[%i][%i] = %e\n", mu, nu, Rs_fluid[mu][nu]);
      }
    } 
    for (int mu = 0; mu < NDIM; mu++) {
      for (int nu = 0; nu < NDIM; nu++) {
        printf("Rs_fluid_transform[%i][%i] = %e\n", mu, nu, Rs_fluid_transform[mu][nu]);
      }
    } */
    //DLOOP2 printf("Ranalytic[%i][%i] = %e\n", mu, nu, Ranalytic[mu][nu]);
    //DLOOP2 printf("Rcoordanalytic[%i][%i] = %e\n", mu, nu, Rcoordanalytic[mu][nu]);
    //DLOOP2 printf("Rscoord[%i][%i] = %e\n", mu, nu, Rs_coord[mu][nu]);
    //double Tg = m.Thetae*ME*CL*CL/KBOL;
    //printf("Tg = %e\n", Tg);
    //double Ert = Rs_coord[0][0];
    //double Tr = pow(Ert*U_unit/AR,1./4.);
    //printf("Tr = %e\n", Tr);
    //double Er_analytic = Rcoordanalytic[0][0];
    //double Tr_analytic = pow(Er_analytic*U_unit/AR,1./4.);
    //printf("Tr_a = %e\n", Tr_analytic);

    //exit(-1);
    //phtmp = ph;
    //int nzone = 0;
    //while (phtmp != NULL) {
    //  nzone++;
    //  phtmp = phtmp->next;
    //}
    //printf("[%i] nzone = %i\n", i, nzone);

    //ph->next = photon_lists[omp_get_thread_num()];
    tail->next = photon_lists[omp_get_thread_num()];
    photon_lists[omp_get_thread_num()] = ph;
  }
}

void bound_gas_prob_x1l(int i, int j, int k, grid_prim_type P)
{
  double X[NDIM], r, th;
  coord(i, j, k, CENT, X);
  bl_coord(X, &r, &th);

  // Get semianalytic solution
  int n = 0;
  while (!(r_init[n+1] > r && r_init[n] < r)) {
    n++;
  }
  double P_cgs = (r - r_init[n])*(P_cgs_init[n+1] - P_cgs_init[n])/(r_init[n+1] - r_init[n]) + P_cgs_init[n];
  double rho_cgs = (r - r_init[n])*(rho_cgs_init[n+1] - rho_cgs_init[n])/(r_init[n+1] - r_init[n]) + rho_cgs_init[n];

  // Get primitive 3-velocity
  struct of_geom geom = ggeom[i][j][CENT];
  double alpha = geom.alpha;
  double beta[NDIM];
  for (int l = 1; l < NDIM; l++) {
    beta[l] = geom.gcon[0][l]*alpha*alpha;
  }

  P[i][j][k][RHO] = rho_cgs/RHO_unit;
  P[i][j][k][UU] = P_cgs/(gam - 1.)/U_unit;
  P[i][j][k][U1] = beta[1]/(sqrt(alpha*alpha-beta[1]*beta[1]*geom.gcov[1][1]));
  P[i][j][k][U2] = beta[2]/(sqrt(alpha*alpha-beta[2]*beta[2]*geom.gcov[2][2]));
  P[i][j][k][U3] = beta[3]/(sqrt(alpha*alpha-beta[3]*beta[3]*geom.gcov[3][3]));
  //P[i][j][k][U2] = 0.;
  //P[i][j][k][U3] = 0.;
  P[i][j][k][B1] = 1.e-10;
  P[i][j][k][B2] = 0.;
  P[i][j][k][B3] = 0.;
}

void bound_gas_prob_x1r(int i, int j, int k, grid_prim_type P)
{
  double X[NDIM], r, th;
  coord(i, j, k, CENT, X);
  bl_coord(X, &r, &th);

  // Get semianalytic solution
  int n = 0;
  while (!(r_init[n+1] > r && r_init[n] < r)) {
    n++;
  }
  double P_cgs = (r - r_init[n])*(P_cgs_init[n+1] - P_cgs_init[n])/(r_init[n+1] - r_init[n]) + P_cgs_init[n];
  double rho_cgs = (r - r_init[n])*(rho_cgs_init[n+1] - rho_cgs_init[n])/(r_init[n+1] - r_init[n]) + rho_cgs_init[n];

  // Get primitive 3-velocity
  struct of_geom geom = ggeom[i][j][CENT];
  double alpha = geom.alpha;
  double beta[NDIM];
  for (int l = 1; l < NDIM; l++) {
    beta[l] = geom.gcon[0][l]*alpha*alpha;
  }

  P[i][j][k][RHO] = rho_cgs/RHO_unit;
  P[i][j][k][UU] = P_cgs/(gam - 1.)/U_unit;
  P[i][j][k][U1] = beta[1]/(sqrt(alpha*alpha-beta[1]*beta[1]*geom.gcov[1][1]));
  P[i][j][k][U2] = beta[2]/(sqrt(alpha*alpha-beta[2]*beta[2]*geom.gcov[2][2]));
  P[i][j][k][U3] = beta[3]/(sqrt(alpha*alpha-beta[3]*beta[3]*geom.gcov[3][3]));
  //P[i][j][k][U2] = 0.;
  //P[i][j][k][U3] = 0.;
  P[i][j][k][B1] = 1.e-10;
  P[i][j][k][B2] = 0.;
  P[i][j][k][B3] = 0.;
}

void init_zone_thermal(int i, int j, int k, double dndlnu[NU_BINS_EMISS+1])
{
  double dlnu = (log(numax_emiss) - log(numin_emiss))/NU_BINS_EMISS;
  double nusamp[NU_BINS_EMISS+1];
  for (int n = 0; n <= NU_BINS_EMISS; n++) {
    nusamp[n] = exp(n*dlnu + log(numin_emiss));
    dndlnu[n] = 0.;
  }

  double dndlnu_max = 0.;
  
  double g;
  double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
  struct of_microphysics m;
  get_fluid_zone(i, j, k, P, &m, Ucon, Ucov, Bcon, Bcov);
  g = ggeom[i][j][CENT].g;
  double dV = dx[1]*dx[2]*dx[3]*pow(L_unit,3)*g;
  

  for (int n = 0; n <= NU_BINS_EMISS; n++) {
    double Bnu = Bnu_inv(nusamp[n], &m)*pow(nusamp[n],3);
    dndlnu[n] = 4.*M_PI*dV*Bnu/(nusamp[n]*CL);
    if (dndlnu[n] > dndlnu_max) {
      dndlnu_max = dndlnu[n];
    }
  }
  
  for (int n = 0; n <= NU_BINS_EMISS; n++) {
    // Equal energy per superphoton
    dndlnu[n] /= dndlnu_max;
    
    // Equal superphotons per d log nu
    //dndlnu[n] = 1.;
  }
}

