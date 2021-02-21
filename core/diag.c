/******************************************************************************
 *                                                                            *
 * DIAG.C                                                                     *
 *                                                                            *
 * DIAGNOSTIC OUTPUT                                                          *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

void reset_log_variables()
{
  #if RADIATION
  step_tot = step_lost = step_made = step_abs = step_scatt = step_rec = 0;
  step_sent = step_rcvd = step_fail = 0;
  #endif
}

void reset_dump_variables()
{
  #if RADIATION
  memset(nuLnu, 0, (MAXNSCATT+1)*NTH*NPHI*NU_BINS_SPEC*sizeof(double));
  memset(Jrad, 0,
    (MAXNSCATT+2)*(N1+2*NG)*(N2+2*NG)*(N3+2*NG)*sizeof(double));
  memset(Nem, 0, (N1+2*NG)*(N2+2*NG)*(N3+2*NG)*sizeof(int));
  memset(Nabs, 0, (N1+2*NG)*(N2+2*NG)*(N3+2*NG)*sizeof(int));
  memset(Nsc, 0, (N1+2*NG)*(N2+2*NG)*(N3+2*NG)*sizeof(int));
  memset(Esuper, 0, (N1+2*NG)*(N2+2*NG)*(N3+2*NG)*sizeof(double));
  memset(Nsuper, 0, (N1+2*NG)*(N2+2*NG)*(N3+2*NG)*sizeof(int));
  #endif
}

void diag(int call_code)
{
  double U[NVAR], divb;
  double pp = 0. ;
  double divbmax = 0.;
  double rmed = 0.;
  double e = 0.;
  struct of_geom *geom;
  struct of_state q;
  static FILE *ener_file;
  
  // Write diagnostics to dump directory
  char log_fnam[STRLEN];
  strcpy(log_fnam, dumpdir);
  strcat(log_fnam, "diag.out");

  if (call_code == DIAG_INIT) {
    // Set things up
    if(mpi_io_proc()) {
      if (is_restart) {
        ener_file = fopen(log_fnam, "a");
      } else {
        ener_file = fopen(log_fnam, "w");
      }
      if (ener_file == NULL) {
        fprintf(stderr,
          "error opening energy output file\n");
        exit(1);
      }
    }
  }

  // Calculate conserved quantities
  if ((call_code == DIAG_INIT || call_code == DIAG_LOG ||
       call_code == DIAG_FINAL) && !failed) {
    pp = 0.;
    e = 0.;
    rmed = 0.;
    divbmax = 0.;
    ZSLOOP(0, N1 - 1, 0, N2 - 1, 0, N3 - 1) {
      geom = get_geometry(i, j, k, CENT) ;
      get_state(P[i][j][k], geom, &q);
      primtoflux(P[i][j][k], &q, 0, geom, U);

      rmed += U[RHO] * dV;
      pp += U[U3] * dV;
      e += U[UU] * dV;
      
      divb = flux_ct_divb(i, j, k);

      // Don't count divbs at simulation volume boundaries, unless periodic
      if (global_start[1] == 0 && i <= NG && N1TOT > 1) continue;
      if (global_stop[1] == N1TOT && i >= N1+NG-1 && N1TOT > 1) continue;
      if (global_start[2] == 0 && j <= NG && N2TOT > 1) continue;
      if (global_stop[2] == N2TOT && j >= N2+NG-1 && N2TOT > 1) continue;
      #if X3L_GAS_BOUND != PERIODIC
      if (global_start[3] == 0 && j <= NG && N3TOT > 1) continue;
      if (global_stop[3] == N3TOT && j >= N3+NG-1 && N3TOT > 1) continue;
      #endif

      if (divb > divbmax) {
        divbmax = divb;
      }
    }
  }

  rmed = mpi_io_reduce(rmed);
  pp = mpi_io_reduce(pp);
  e = mpi_io_reduce(e);
  divbmax = mpi_io_max(divbmax);

  #if RADIATION
  set_Rmunu();
  #endif

  // Get total mass and energy
  double mass_proc = 0.;
  double egas_proc = 0.;
  #if RADIATION
  double erad_proc = 0.;
  #endif
  double Phi_proc = 0.;
  double jet_EM_flux_proc = 0.;
  double lum_eht_proc = 0.; 
  ZLOOP {
    struct of_state q;
    double U[NVAR];
    get_state(P[i][j][k], &(ggeom[i][j][CENT]), &q);
    primtoflux(P[i][j][k], &q, 0, &(ggeom[i][j][CENT]), U);
    mass_proc += U[0]*dx[1]*dx[2]*dx[3];
    egas_proc += U[1]*dx[1]*dx[2]*dx[3];
    double rho = P[i][j][k][RHO];
    double Pg = (gam - 1.)*P[i][j][k][UU];
    double bsq = dot(q.bcon, q.bcov);
    double Bmag = sqrt(bsq);
    double C_eht = 0.2;
    double j_eht = pow(rho,3.)*pow(Pg,-2.)*exp(-C_eht*pow(rho*rho/(Bmag*Pg*Pg),1./3.));
    lum_eht_proc += j_eht*dx[1]*dx[2]*dx[3]*ggeom[i][j][CENT].g;
    if (global_start[1] == 0 && i == 5+NG) {
      Phi_proc += 0.5*fabs(P[i][j][k][B1])*dx[2]*dx[3]*ggeom[i][j][CENT].g;

      double P_EM[NVAR];
      PLOOP P_EM[ip] = P[i][j][k][ip];
      P_EM[RHO] = 0.;
      P_EM[UU] = 0.;
      get_state(P_EM, &(ggeom[i][j][CENT]), &q);
      double sig = dot(q.bcon, q.bcov)/P[i][j][k][RHO];
      if (sig > 1.) {
        primtoflux(P_EM, &q, 1, &(ggeom[i][j][CENT]), U);
        jet_EM_flux_proc += -U[U1]*dx[2]*dx[3];
      }
    }

    #if RADIATION
    erad_proc += Rmunu[i][j][k][0][0]*dx[1]*dx[2]*dx[3]*ggeom[i][j][CENT].g;
    #endif
  }
  double mass = mpi_reduce(mass_proc);
  double egas = mpi_reduce(egas_proc);
  double Phi = sqrt(4.*M_PI)*mpi_reduce(Phi_proc); // cgs unit convention
  double phi = Phi/sqrt(-mdot + SMALL);
  double jet_EM_flux = mpi_reduce(jet_EM_flux_proc);
  double lum_eht = mpi_reduce(lum_eht_proc);
  #if RADIATION
  double erad = mpi_reduce(erad_proc);
  double lum_proc = 0.;

  // Get last shell entirely enclosed by stopx_rad[1] (assume r(X^1) independent
  // of X^2, X^3)
  int stopi_rad = -1;
  for (int i = NG+1; i <= N1 + NG; i++) {
    double X[NDIM], Xprev[NDIM];
    coord(i-1, NG, NG, FACE1, Xprev);
    coord(i, NG, NG, FACE1, X);
    if (X[1] > stopx_rad[1] && Xprev[1] < stopx_rad[1]) {
      stopi_rad = i - 2;
    }
  }
  if (stopi_rad == -1) stopi_rad = N1 + NG - 1;

  if (stopi_rad >= 0) {
    JSLOOP(0, N2-1) {
      KSLOOP(0, N3-1) {
        lum_proc -= Rmunu[stopi_rad][j][k][1][0]*dx[2]*dx[3]*ggeom[stopi_rad][j][CENT].g;
			}
    }
  }
  double lum = mpi_reduce(lum_proc);
  double eff = lum/(mdot+SMALL);
  
  #if ELECTRONS
  int num_super = 0.;
  double lum_super = 0.;
  ZLOOP {
    num_super += Nsuper[i][j][k];
    lum_super += Esuper[i][j][k];
  }
  num_super = mpi_reduce_int(num_super);
  lum_super = mpi_reduce(lum_super)/DTd;
  #endif
  #endif // RADIATION

  if (call_code == DIAG_DUMP) {
    dump();
    dump_cnt++;
  }

  if (call_code == DIAG_FINAL) {
    dump();
  }

  if (call_code == DIAG_INIT ||
      call_code == DIAG_LOG || call_code == DIAG_FINAL) {
    double mdot_all = mpi_io_reduce(mdot);
    double edot_all = mpi_io_reduce(edot);
    double ldot_all = mpi_io_reduce(ldot);
    double mdot_eh_all = mpi_io_reduce(mdot_eh);
    double edot_eh_all = mpi_io_reduce(edot_eh);
    double ldot_eh_all = mpi_io_reduce(ldot_eh);

    if(mpi_io_proc()) {
      fprintf(stdout, "LOG      t=%g \t divbmax: %g\n",
        t,divbmax);
      fprintf(ener_file, "%10.5g %10.5g %10.5g %10.5g %15.8g %15.8g ",
        t, rmed, pp, e,
        P[N1/2][N2/2][N3/2][UU]*pow(P[N1/2][N2/2][N3/2][RHO], -gam),
        P[N1/2][N2/2][N3/2][UU]);
      fprintf(ener_file, "%15.8g %15.8g %15.8g ", mdot_all, edot_all, ldot_all);
      fprintf(ener_file, "%15.8g %15.8g ", mass, egas);
      fprintf(ener_file, "%15.8g %15.8g %15.8g ", Phi, phi, jet_EM_flux);
      fprintf(ener_file, "%15.8g ", divbmax);
      #if RADIATION
      printf("step_abs = %i step_abs_tot = %i\n", step_abs, step_abs_all);
      fprintf(ener_file, "%i %i %i %i %i %i %i %i ", step_made, step_abs,
        step_scatt, step_lost, step_rec, step_tot, step_sent, step_rcvd);
      fprintf(ener_file, "%i %i %i %i %i %i %i %i %i ", step_made_all, 
        step_abs_all, step_scatt_all, step_lost_all, step_rec_all, step_tot_all, 
        step_sent_all, step_rcvd_all, step_fail_all);
      fprintf(ener_file, "%15.8g %15.8g ", tune_emiss, tune_scatt);
      fprintf(ener_file, "%15.8g %15.8g %15.8g ", erad, lum, eff);
      #if ELECTRONS
      fprintf(ener_file, "%i %15.8g ", num_super, lum_super);
      #endif
      #endif
      fprintf(ener_file, "%15.8g ", lum_eht);
      fprintf(ener_file, "%15.8g %15.8g %15.8g ", mdot_eh_all, edot_eh_all, ldot_eh_all);
      fprintf(ener_file, "%15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g ", 
        get_time_per_step(TIMER_UPDATE), get_time_per_step(TIMER_FLUXCALC),
        get_time_per_step(TIMER_FIXUP), get_time_per_step(TIMER_BOUND),
        get_time_per_step(TIMER_DIAG), get_time_per_step(TIMER_OUT),
        get_time_per_step(TIMER_MAKE), get_time_per_step(TIMER_PUSH), 
        get_time_per_step(TIMER_INTERACT), get_time_per_step(TIMER_ALL));
      #if ELECTRONS
      fprintf(ener_file, "%15.8g ", get_time_per_step(TIMER_ELECTRON));
      #endif
      fprintf(ener_file, "\n");
      fflush(ener_file);
    }

    #if RADIATION
    track_ph();
    #endif
  }
}

// Diagnostic routines
void fail(int fail_type)
{
  failed = 1;

  fprintf(stderr, "\n\nFAIL: [%d %d %d] %d\n", icurr, jcurr, kcurr, fail_type);

  area_map(icurr, jcurr, kcurr, P);

  diag(DIAG_FINAL);

  exit(0);
}

// Map out region around failure point
void area_map(int i, int j, int k, grid_prim_type prim)
{
  fprintf(stderr, "*** AREA MAP ***\n");

  PLOOP {
    fprintf(stderr, "variable %d \n", ip);
    fprintf(stderr, "i = \t %12d %12d %12d\n", i - 1, i,
      i + 1);
    fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j + 1,
      prim[i - 1][j + 1][k][ip], prim[i][j + 1][k][ip],
      prim[i + 1][j + 1][k][ip]);
    fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j,
      prim[i - 1][j][k][ip], prim[i][j][k][ip],
      prim[i + 1][j][k][ip]);
    fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j - 1,
      prim[i - 1][j - 1][k][ip], prim[i][j - 1][k][ip],
      prim[i + 1][j - 1][k][ip]);
  }

  fprintf(stderr, "****************\n");
}

// Evaluate flux based diagnostics; put results in global variables
void diag_flux(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3)
{
  mdot = edot = ldot = 0.;
  mdot_eh = edot_eh = ldot_eh = 0.;
  int iEH = NG + 5;
  if (global_start[1] == 0) {
    #pragma omp parallel for \
      reduction(+:mdot) reduction(+:edot) reduction(+:ldot) \
      reduction(+:mdot_eh) reduction(+:edot_eh) reduction(+:ldot_eh) \
      collapse(2)
    JSLOOP(0, N2 - 1) {
      KSLOOP(0, N3 - 1) {
        mdot += F1[NG][j][k][RHO]*dx[2]*dx[3];
        edot += (F1[NG][j][k][UU] - F1[NG][j][k][RHO])*dx[2]*dx[3];
        ldot += F1[NG][j][k][U3]*dx[2]*dx[3];
        mdot_eh += F1[iEH][j][k][RHO]*dx[2]*dx[3];
        edot_eh += (F1[iEH][j][k][UU] - F1[iEH][j][k][RHO])*dx[2]*dx[3];
        ldot_eh += F1[iEH][j][k][U3]*dx[2]*dx[3];
      }
    }
  }
}

double flux_ct_divb(int i, int j, int k)
{
  return fabs(0.25*(
    P[i][j][k][B1]*ggeom[i][j][CENT].g
    + P[i][j-1][k][B1]*ggeom[i][j-1][CENT].g
    + P[i][j][k-1][B1]*ggeom[i][j][CENT].g
    + P[i][j-1][k-1][B1]*ggeom[i][j-1][CENT].g
    - P[i-1][j][k][B1]*ggeom[i-1][j][CENT].g
    - P[i-1][j-1][k][B1]*ggeom[i-1][j-1][CENT].g
    - P[i-1][j][k-1][B1]*ggeom[i-1][j][CENT].g
    - P[i-1][j-1][k-1][B1]*ggeom[i-1][j-1][CENT].g
    )/dx[1] +
    0.25*(
    P[i][j][k][B2]*ggeom[i][j][CENT].g
    + P[i-1][j][k][B2]*ggeom[i-1][j][CENT].g
    + P[i][j][k-1][B2]*ggeom[i][j][CENT].g
    + P[i-1][j][k-1][B2]*ggeom[i-1][j][CENT].g
    - P[i][j-1][k][B2]*ggeom[i][j-1][CENT].g
    - P[i-1][j-1][k][B2]*ggeom[i-1][j-1][CENT].g
    - P[i][j-1][k-1][B2]*ggeom[i][j-1][CENT].g
    - P[i-1][j-1][k-1][B2]*ggeom[i-1][j-1][CENT].g
    )/dx[2] +
    0.25*(
    P[i][j][k][B3]*ggeom[i][j][CENT].g
    + P[i][j-1][k][B3]*ggeom[i][j-1][CENT].g
    + P[i-1][j][k][B3]*ggeom[i-1][j][CENT].g
    + P[i-1][j-1][k][B3]*ggeom[i-1][j-1][CENT].g
    - P[i][j][k-1][B3]*ggeom[i][j][CENT].g
    - P[i][j-1][k-1][B3]*ggeom[i][j-1][CENT].g
    - P[i-1][j][k-1][B3]*ggeom[i-1][j][CENT].g
    - P[i-1][j-1][k-1][B3]*ggeom[i-1][j-1][CENT].g
    )/dx[3]);
}

#if RADIATION
void record_superphoton(double X[NDIM], struct of_photon *ph)
{
  double lnumin = log(numin_spec);
  double lnumax = log(numax_spec);
  double dlnu = (lnumax - lnumin)/NU_BINS_SPEC;
  int i, j, k;
  Xtoijk(X, &i, &j, &k);

  int nscatt = ph->nscatt;
  nscatt = MY_MIN(nscatt, MAXNSCATT);

  // Preserve index sanity
  if (j < NG)       j = N2 + j;
  if (j >= N2 + NG) j = j - N2;
  if (k < NG) {
    if (N3CPU == 1)
      k = N3 + k;
    else
      k = NG;
  }
  if (k >= N3 + NG) {
    if (N3CPU == 1)
      k = k - N3;
    else
      k = N3+NG-1;
  }

  // Assume X0 symmetry in metric
  double nu = -ph->Kcov[2][0]*ME*CL*CL/HPL;

  int thbin, phibin, nubin = (log(nu) - lnumin)/dlnu;
  double dOmega = get_nuLnu_bin(X, &thbin, &phibin);

  // Store dE / dlognu dOmega dt
  if (nubin >= 0 && nubin < NU_BINS_SPEC) {
    #pragma omp atomic
    nuLnu[nscatt][thbin][phibin][nubin] -= ph->w*ph->Kcov[2][0]*ME*CL*CL/
                                  (dlnu*DTd*T_unit) * 4.*M_PI / dOmega;
    #pragma omp atomic
    step_rec++;
  }
}
#endif // RADIATION

