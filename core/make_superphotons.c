/******************************************************************************
 *                                                                            *
 * MAKE_SUPERPHOTONS.C                                                        *
 *                                                                            *
 * EMISSION OF MONTE CARLO SAMPLES                                            *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
static double lnu_min, lnu_max, dlnu, nusamp[NU_BINS_EMISS+1], Ns;
double dnzs[N1+2*NG][N2+2*NG][N3+2*NG];

void sample_photon(int i, int j, int k, double t, double dt,
  double dndlnu[NU_BINS_EMISS+1], struct of_photon *tmp, double Econ[NDIM][NDIM],
  double Ecov[NDIM][NDIM], const struct of_microphysics *m,
  double Bcon[NDIM], double Ucon[NDIM]);
void get_dndlnu(int i, int j, int k, double dt, double dndlnu[NU_BINS_EMISS+1],
  const struct of_microphysics *m);

void init_make_statics() {
  double my_tune_emiss = tune_emiss;
  if (my_tune_emiss < 0) my_tune_emiss = 1.;
  Ns = my_tune_emiss/(pow(sim_vol,1./3.)*T_unit/CL);
  lnu_min = log(numin_emiss);
  lnu_max = log(numax_emiss);
  dlnu = (lnu_max - lnu_min)/NU_BINS_EMISS;
  for (int n = 0; n <= NU_BINS_EMISS; n++) {
    nusamp[n] = exp(n*dlnu + lnu_min);
  }
}

double get_wgt(double nu)
{
  return wgtC/nu;
}

void make_superphotons(grid_prim_type Prad, double t, double dt)
{
  #if EMISSION
  timer_start(TIMER_MAKE);
  get_dnz(Prad);

  int step_made_local = 0;

  #pragma omp parallel reduction(+:step_made_local)
  {
    struct of_photon *tmp, *head = photon_lists[omp_get_thread_num()];
    struct of_microphysics m;
    double dndlnu[NU_BINS_EMISS+1];
    double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];
    double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
    double X[NDIM];
    int nz;

    ZLOOP
    {
      nz = (int)dnzs[i][j][k];
      if (dnzs[i][j][k] - nz > get_rand()) nz++;
      //printf("dnzs[%i] = %e\n", i, dnzs[i][j][k]);

      if (nz > 0) {
        // Set up zone
        coord(i, j, k, CENT, X);
        get_fluid_zone(i, j, k, Prad, &m, Ucon, Ucov, Bcon, Bcov);

        make_tetrad(i, j, k, Ucon, Bcon, ggeom[i][j][CENT].gcov, Econ, Ecov);
        get_dndlnu(i, j, k, dt, dndlnu, &m);

        // Create superphotons in pairs
        for (int n = 0; n < nz; n++) {
          tmp = safe_malloc(sizeof(struct of_photon));
          tmp->next = safe_malloc(sizeof(struct of_photon));

          sample_photon(i, j, k, t, dt, dndlnu, tmp, Econ, Ecov, &m, Bcon, Ucon);

          (tmp->next)->next = head;
          head = tmp;
        } // n < nz

        //#pragma omp atomic
        step_made_local += 2*nz;
      } // nz > 0
    } // ZLOOP

    // Prepend created superphotons to each thread's global list
    photon_lists[omp_get_thread_num()] = head;
  } // omp parallel

  step_made += step_made_local;

  timer_stop(TIMER_MAKE);
  #endif // EMISSION
}

void sample_photon(int i, int j, int k, double t, double dt,
  double dndlnu[NU_BINS_EMISS+1], struct of_photon *ph, double Econ[NDIM][NDIM],
  double Ecov[NDIM][NDIM], const struct of_microphysics *m,
  double Bcon[NDIM], double Ucon[NDIM])
{
  double nu, th, cth[2], sth[2], phi, sphi[2], cphi[2];
  double K_tetrad[NDIM];
  struct of_photon *tmp[2];
  tmp[0] = ph;
  tmp[1] = ph->next;
  tmp[1]->next = NULL;

  // Sample emissivity to get frequency
  do {
    nu = exp(get_rand()*(lnu_max - lnu_min) + lnu_min);
  } while (get_rand() > linear_interp_log(nu, dndlnu, lnu_min, dlnu));

  // Get weight from global weight parameter
  double weight = get_wgt(nu);

  // Sample emissivity in solid angle
  double jmax = jnu(nu, m, 0.5*M_PI);
  do {
    cth[0] = 2.*get_rand() - 1.;
    th = acos(cth[0]);
  } while (get_rand() > jnu(nu, m, th)/jmax);

  sth[0] = sqrt(1. - cth[0]*cth[0]);
  phi = 2.*M_PI*get_rand();
  cphi[0] = cos(phi);
  sphi[0] = sin(phi);

  // Second photon antiparallel in fluid frame
  cth[1]  = -cth[0];
  sth[1]  =  sth[0];
  cphi[1] = -cphi[0];
  sphi[1] = -sphi[0];

  double E = nu*HPL/(ME*CL*CL);

  for (int n = 0; n < 2; n++) {
    // Initial zeros
    memset(tmp[n]->X, 0, NSUP*NDIM*sizeof(double));
    memset(tmp[n]->Kcov, 0, NSUP*NDIM*sizeof(double));
    memset(tmp[n]->Kcon, 0, NSUP*NDIM*sizeof(double));

    // Kludge to keep 3rd-order interpolation happy
    tmp[n]->X[0][0] = -1.;
    tmp[n]->X[1][0] = -1.;

    // Set position
    tmp[n]->X[2][0] = t + dt/2.;
    coord(i, j, k, CENT, tmp[n]->X[2]);

    // Randomize phi for visualization if in axisymmetry and BH spacetime
    if (N3TOT == 1 && (METRIC == MKS || METRIC == MMKS))
      tmp[n]->X[2][3] = 2.*M_PI*get_rand();

    // Get coordinate frame wavevector
    K_tetrad[0] = -E;
    K_tetrad[1] = E*cth[n];
    K_tetrad[2] = E*cphi[n]*sth[n];
    K_tetrad[3] = E*sphi[n]*sth[n];

    tetrad_to_coord(Ecov, K_tetrad, tmp[n]->Kcov[2]);

    K_tetrad[0] *= -1.;
    tetrad_to_coord(Econ, K_tetrad, tmp[n]->Kcon[2]);

    // Re-do this to ensure k.k == 0?

    // Set superphoton weight
    tmp[n]->w = 0.5*weight;

    // Diagnostics
    tmp[n]->nscatt = 0;
    tmp[n]->origin[0] = nstep;
    tmp[n]->origin[1] = i;
    tmp[n]->origin[2] = j;
    tmp[n]->origin[3] = k;

    tmp[n]->t0 = t + dt/2.;

    if (!is_null(tmp[n]->Kcov[2], tmp[n]->Kcon[2], tmp[n]->Kcov[2][0], 0.,
          &(tmp[n]->KdotKprev)))
    {
      fprintf(stderr, "Error! K.K != 0 initially!\n");
      fprintf(stderr, "K.K make err [%i %i %i] nu = %e w = %e n = %i K.K = %e\n", i,j,k,nu,weight,n, tmp[n]->KdotKprev);
      double gamma;
      mhd_gamma_calc(P[i][j][k], &(ggeom[i][j][CENT]), &gamma);
      fprintf(stderr, "K_0 = %e gamma = %e\n", tmp[n]->Kcov[2][0], gamma);
    
    }

    if (tmp[n]->Kcov[2][0] > 0.) {
      tmp[n]->w = 0.;
    }

    // Record radiation four-force
    double Gcov[NDIM];
    for (int mu = 0; mu < NDIM; mu++) {
      Gcov[mu] = -1/(ggeom[i][j][CENT].g*dt*dx[1]*dx[2]*dx[3])*kphys_to_num*tmp[n]->w*tmp[n]->Kcov[2][mu];
      #pragma omp atomic
      radG[i][j][k][mu] += Gcov[mu]*ggeom[i][j][CENT].g;
    }

    // du_e / dtau
    #pragma omp atomic
    Jrad[0][i][j][k] -= dot(Ucon, Gcov);

    #pragma omp atomic
    Nem[i][j][k] += 1;

    if (get_rand() < ((double)Nph_to_track)/(nph_per_proc*mpi_nprocs())) {
      tmp[n]->is_tracked = 1;
    } else {
      tmp[n]->is_tracked = 0;
    }
  }
}

#define TINY (1.e-200)
void get_dndlnu(int i, int j, int k, double dt, double dndlnu[NU_BINS_EMISS+1],
  const struct of_microphysics *m)
{
  for (int n = 0; n < NU_BINS_EMISS; n++) {
    dndlnu[n] = 0.;
  }

  double dndlnu_max = -1.e100;
  for (int n = 0; n <= NU_BINS_EMISS; n++) {
    double Jsamp = Jnu(nusamp[n], m);
    Jsamp *= dx[1]*dx[2]*dx[3]*pow(L_unit,3.)*ggeom[i][j][CENT].g;

    double wgt = get_wgt(nusamp[n]);

    //dndlnu[n] = Jsamp/(wgtC/nusamp[n]*HPL + TINY);
    dndlnu[n] = Jsamp/(wgt*HPL + TINY);

    if (dndlnu[n] > dndlnu_max) {
      dndlnu_max = dndlnu[n];
    }
  }

  for (int n = 0; n <= NU_BINS_EMISS; n++) {
    dndlnu[n] /= dndlnu_max;
  }
}
#undef TINY

void set_weight(grid_prim_type Prad)
{
  double Jtot;
  double zoneVol = dV*L_unit*L_unit*L_unit;

  // Set static variables
  init_make_statics();
  //Ns = tune_emiss/(pow(sim_vol,1./3.)*T_unit/CL);
  //lnu_min = log(numin_emiss);
  //lnu_max = log(numax_emiss);
  //dlnu = (lnu_max - lnu_min)/NU_BINS_EMISS;
  //for (int n = 0; n <= NU_BINS_EMISS; n++) {
  //  nusamp[n] = exp(n*dlnu + lnu_min);
  //}
  Jtot = 0.;

  #pragma omp parallel
  {
    #pragma omp for collapse(3) reduction(+:Jtot)
    ZLOOP
    {
      struct of_microphysics m;
      double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
      get_fluid_zone(i, j, k, Prad, &m, Ucon, Ucov, Bcon, Bcov);

      for (int n = 0; n <= NU_BINS_EMISS; n++) {
        Jtot += Jnu(nusamp[n], &m)*zoneVol*ggeom[i][j][CENT].g;
      }
    } // ZLOOP
  } // omp parallel

  Jtot = mpi_reduce(Jtot);

  wgtC = Jtot/(HPL*Ns)*nusamp[0];
}

double f (double x, void *microphysics) {
  struct of_microphysics *m = (struct of_microphysics *) microphysics;
  double nu     = exp(x);
  double Jsamp  = Jnu(nu, m)*nu;
  double wgt = get_wgt(nu);
  if (isinf(Jsamp) || wgt < SMALL) {
    return 0.;
  } else {
    return Jsamp/(nu*wgt);
  }
}

// Calculate number of superphotons to produce per thread in each zone
void get_dnz(grid_prim_type Prad)
{
  #pragma omp parallel
  {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function F;
    F.function = &f;

    double zoneVol = dV*L_unit*L_unit*L_unit;
    struct of_microphysics m;
    double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];

    #pragma omp for collapse(3) schedule(dynamic)
    ZLOOP {
      // Ignore emission outside region of interest
      double X[NDIM];
      coord(i, j, k, CENT, X);
      if (X[1] < startx_rad[1] || X[1] > stopx_rad[1]) {
        dnzs[i][j][k] = 0.;
        continue;
      }

      get_fluid_zone(i, j, k, Prad, &m, Ucon, Ucov, Bcon, Bcov);

      // Get number of superphotons to be emitted
      F.params = &m;
      gsl_integration_qags(&F, lnu_min, lnu_max, 1.e100, 1.e-4, 1000, w,
        &result, &error);
      result /= HPL;
      result *= zoneVol;
      result *= ggeom[i][j][CENT].g;
      result *= dt*T_unit;

      if (isnan(result/nthreads)) {
        dnzs[i][j][k] = 0.;
      } else {
        dnzs[i][j][k] = result/nthreads;
      }
    } // ZLOOP
    gsl_integration_workspace_free(w);
  } // pragma omp parallel
}

#endif // RADIATION

