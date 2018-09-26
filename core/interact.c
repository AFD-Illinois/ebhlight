/******************************************************************************
 *                                                                            *
 * INTERACT.C                                                                 *
 *                                                                            *
 * PROCESS ABSORPTION AND SCATTERING                                          *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#if RADIATION

grid_fourvector_type Ucon_grid, Ucov_grid, Bcon_grid, Bcov_grid;
struct of_microphysics m_grid[N1+2*NG][N2+2*NG][N3+2*NG];

void get_fluid_zone_grid(int i, int j, int k, struct of_microphysics *m, 
  double Ucon[NDIM], double Ucov[NDIM], double Bcon[NDIM], double Bcov[NDIM])
{
  m->Thetae = m_grid[i][j][k].Thetae;
  m->Ne = m_grid[i][j][k].Ne;
  m->B = m_grid[i][j][k].B;
  DLOOP1 {
    Ucon[mu] = Ucon_grid[i][j][k][mu];
    Ucov[mu] = Ucov_grid[i][j][k][mu];
    Bcon[mu] = Bcon_grid[i][j][k][mu];
    Bcov[mu] = Bcov_grid[i][j][k][mu];
  }
}

double get_scatt_bias(double nu, const struct of_microphysics *m, double uph)
{
  double Thetae = m->Thetae;
  double amp = 1. + 4.*Thetae - 2.*pow(Thetae,3./2.) + 16.*pow(Thetae,2.);
  double bias = tune_scatt*amp;

  // Ensure bias is in reasonable physical bounds
  if (bias < 1.) {
    bias = 1.;
  }

  // Another BOUND_BIAS method. Assumes large hotspots, may work fine instead
  // assuming ~GM/c^2 length scale for hot spots.
  double dl = Rout_rad*L_unit;
  double dtau = (alpha_inv_scatt(nu, m)/nu)*dl;
  bias = MY_MIN(bias, 1./dtau);
  bias = MY_MAX(bias, 1.);

  return bias;
}
#undef BOUND_BIAS

#define MAX_INTERACTIONS 100
void interact(grid_prim_type P, double t, double dt)
{
  timer_start(TIMER_INTERACT);  

  // Set Ucon, Ucov, Bcon, Bcov, m for each zone
  #pragma omp parallel for collapse(2)
  ZLOOPALL {
    get_fluid_zone(i, j, k, P, &m_grid[i][j][k], Ucon_grid[i][j][k], 
    Ucov_grid[i][j][k], Bcon_grid[i][j][k], Bcov_grid[i][j][k]);
  }

  #if ABSORPTION || SCATTERING
  double d3x = dx[1]*dx[2]*dx[3];
  int step_tot_local = 0;
  int step_abs_local = 0;
  int step_fail_local = 0;
  int step_scatt_local = 0;
  #pragma omp parallel reduction(+:step_tot_local) reduction(+:step_abs_local) \
                       reduction(+:step_fail_local) reduction(+:step_scatt_local)
  {
    int i, j, k;
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    struct of_photon *prev = NULL;
    struct of_photon *head = ph;
    struct of_microphysics m;
    double X[NDIM], Kcon[NDIM], Kcov[NDIM];
    double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
    double dlam, nu, dtau_abs = 0., dtau_scatt = 0., bias_scatt = 1.;
    double xabs, xscatt;
    double dlfrac;
    int nint;
    while (ph != NULL)
    {
      if (ph->w < SMALL) {
        prev = ph;
        ph = ph->next;
        continue;
      }
      //double tmin = t;
      double tmin = MY_MAX(t, ph->t0);
      // Loop with decreasing d \lambda to account for multiple interactions
      // per superphoton per timestep (scattering only)
      nint = 0;
      while (nint < MAX_INTERACTIONS) {
        nint++;

        double dtint = t+dt-tmin;
        int status = get_X_K_interp(ph, MY_MAX(t+dt/2.,ph->t0), X, Kcov, Kcon);
        if (status == SPH_INTERP_FAIL) {
          prev = ph;
          ph = ph->next;
          break;
        }
        
        dlam = dtint/Kcon[0];
        
        Xtoijk(X, &i, &j, &k);

        if (i < 0 || i > N1+2*NG-1 ||
            j < 0 || j > N2+2*NG-1 ||
            k < 0 || k > N3+2*NG-1)
        {
          fprintf(stderr, "PHOTON IN BAD PLACE! %i %i %i\n", i,j,k);
          print_ph_diag(ph);
          list_remove(&ph, &head, &prev);
          step_tot_local--;
          break;
        }

        // Get quantities relating radiation and fluid
        //get_fluid_zone(i, j, k, P, &m, Ucon, Ucov, Bcon, Bcov);
        get_fluid_zone_grid(i, j, k, &m, Ucon, Ucov, Bcon, Bcov);

        if (m.Thetae < 10.*SMALL) {
          prev = ph;
          ph = ph->next;
          break;
        }
        
        double sigma = pow(m.B/B_unit,2.)/(m.Ne/Ne_unit);
        if (sigma > sigma_max) {
          prev = ph;
          ph = ph->next;
          break;
        }

        nu = get_fluid_nu(X, Kcov, Ucon);

        if (nu <= 0 || isnan(nu)) {
          fprintf(stderr, "Bad NU in interact [%i %i %i]\n", i, j, k);
          fprintf(stderr, "X[] = %e %e %e %e\n", X[0], X[1], X[2], X[3]);
          fprintf(stderr, "Kcov[] = %e %e %e %e\n", 
            Kcov[0], Kcov[1], Kcov[2], Kcov[3]);
          fprintf(stderr, "Ucon[] = %e %e %e %e\n", 
            Ucon[0], Ucon[1], Ucon[2], Ucon[3]);
          double gamma;
          mhd_gamma_calc(P[i][j][k], &(ggeom[i][j][CENT]), &gamma);
          list_remove(&ph, &head, &prev);
          step_tot_local--;
          break;
        }

        // Calculate and sample optical depths along step
        #if ABSORPTION
        double theta = get_bk_angle(X, Kcon, Ucov, Bcov, m.B);
        dtau_abs = (HPL*L_unit/(ME*CL*CL))*dlam*alpha_inv_abs(nu, &m, theta);
        #endif

        #if SCATTERING
        dtau_scatt = (HPL*L_unit/(ME*CL*CL))*dlam*alpha_inv_scatt(nu, &m);
        bias_scatt = get_scatt_bias(nu, &m, 
          HPL*nu*ph->w/(ggeom[i][j][CENT].g*d3x*pow(L_unit,3)));
        #endif

        xabs = -log(get_rand());
        xscatt = -log(get_rand())/bias_scatt;

        // No interaction
        if (xabs > dtau_abs && xscatt > dtau_scatt) {
          prev = ph;
          ph = ph->next;
          break;
        }

        // Absorption
        #if ABSORPTION
        else if (xabs/(dtau_abs + SMALL) < xscatt/(dtau_scatt + SMALL)) {
          dlfrac = xabs/dtau_abs;
          double tabs = tmin+dlfrac*dtint;
          int status = push_to_X_K(tabs, ph, X, Kcov, Kcon);
          if (status == PUSH_FAIL) {
            prev = ph;
            ph = ph->next;
            break;
          }

          Xtoijk(X, &i, &j, &k);
        if (i < 0 || i > N1+2*NG-1 ||
            j < 0 || j > N2+2*NG-1 ||
            k < 0 || k > N3+2*NG-1)
        {
          printf("BAD!\n");
          printf("i = %i j = %i k = %i\n", i, j, k);
          exit(-1);
          }

          if (!strcmp(PROBLEM_NAME, "bhtherm")) {
            if (i < NG && global_start[1] == 0) i = NG;
            if (i >= NG + N1 && global_stop[1] == N1TOT) i = NG + N1 - 1;
            if (j < NG && global_start[2] == 0) j = NG;
            if (j >= NG + N2 && global_stop[2] == N2TOT) j = NG + N2 - 1;
          }

          // Boundary transport cannot use MPI with one zone
          if (N1 == 1) i = NG;
          if (N2 == 1) j = NG;
          if (N3 == 1) k = NG;

          double Gcov[NDIM];
          for (int mu = 0; mu < NDIM; mu++) {
            Gcov[mu] = 1./(ggeom[i][j][CENT].g*dt*d3x)*ph->w*kphys_to_num*Kcov[mu];
            #pragma omp atomic
            radG[i][j][k][mu] += Gcov[mu]*ggeom[i][j][CENT].g;
          }
          
          // du_e / dtau
          #pragma omp atomic
          Jrad[1][i][j][k] += dot(Ucon_grid[i][j][k], Gcov);

          step_abs_local++;

          #pragma omp atomic
          Nabs[i][j][k]++;

          ph->w = 0.;

          tmin = tabs;

          prev = ph;
          ph = ph->next;
          break;
        }
        #endif

        else {
          dlfrac = xscatt/dtau_scatt;
          double tscatt = tmin + dlfrac*dtint;
          int status = push_to_X_K(tscatt, ph, X, Kcov, Kcon);
          if (status == PUSH_FAIL) {
            prev = ph;
            ph = ph->next;
            break;
          }

          struct of_photon *phscatt = safe_malloc(sizeof(struct of_photon));

          if (get_rand() < Nph_to_track/(nph_per_proc*mpi_nprocs())) {
            phscatt->is_tracked = 1;
          } else {
            phscatt->is_tracked = 0;
          }

          Xtoijk(X, &i, &j, &k);

          // Initialize scattered photon at position of scattering event
          for (int mu = 0; mu < NDIM; mu++) {
            phscatt->X[2][mu] = X[mu];
            phscatt->Kcov[2][mu] = Kcov[mu];
            phscatt->Kcon[2][mu] = Kcon[mu];
            for (int n = 0; n < 2; n++) {
              phscatt->X[n][mu] = 0.;
              phscatt->Kcov[n][mu] = 0.;
              phscatt->Kcon[n][mu] = 0.;
            }
          } 
          phscatt->X[0][0] = -1.;
          phscatt->X[1][0] = -1.;
          phscatt->t0 = tscatt;
          phscatt->w = ph->w/bias_scatt;
          phscatt->nscatt = ph->nscatt + 1;
          phscatt->origin[0] = nstep;
          phscatt->origin[1] = i;
          phscatt->origin[2] = j;
          phscatt->origin[3] = k;
          ph->w = (1. - 1./bias_scatt)*ph->w;

          int success = scatter_superphoton(P, phscatt, X, Kcov, Kcon);

          if (!success) {
            step_fail_local++;
            free(phscatt);
            prev = ph;
            ph = ph->next;
            break;
          }

          // Need to reset K.K
          phscatt->KdotKprev = 0.;
          for (int mu = 0; mu < NDIM; mu++) {
            phscatt->KdotKprev += phscatt->Kcov[2][mu]*phscatt->Kcon[2][mu];
          }

          // Boundary transport cannot use MPI with one zone
          if (N1 == 1) i = NG;
          if (N2 == 1) j = NG;
          if (N3 == 1) k = NG;

          // Apply four-force at interaction site
          double Gcov[NDIM];
          for (int mu = 0; mu < NDIM; mu++) {
            Gcov[mu] = 1./(ggeom[i][j][CENT].g*dt*d3x)*phscatt->w*kphys_to_num*(Kcov[mu] - phscatt->Kcov[2][mu]);
            #pragma omp atomic
            radG[i][j][k][mu] += Gcov[mu]*ggeom[i][j][CENT].g;
          }
          
          // du_e / dtau
          int nscatt = MY_MIN(ph->nscatt, MAXNSCATT - 1);
          #pragma omp atomic
          Jrad[nscatt+2][i][j][k] += dot(Ucon_grid[i][j][k], Gcov);
          
          //push phscatt as far as possible
          status = push_superphoton(phscatt, cour*dt_light[i][j]);
          if (status == PUSH_FAIL) {
            free(phscatt);
            prev = ph;
            ph = ph->next;
            step_fail_local++;
            break;
          }

          step_scatt_local++;

          step_tot_local++;

          #pragma omp atomic
          Nsc[i][j][k]++;

          // Add scattered superphoton to list
          phscatt->next = ph->next;
          ph->next = phscatt;

          // Allow for additional scatterings
          tmin = tscatt;
          continue;
        }
      }
    } // ph != NULL
    photon_lists[omp_get_thread_num()] = head;

  } // omp parallel
  step_tot += step_tot_local;
  step_abs += step_abs_local;
  step_fail += step_fail_local;
  step_scatt += step_scatt_local;
  #endif // ABSORPTION || SCATTERING

  timer_stop(TIMER_INTERACT);
}
#endif // RADIATION

