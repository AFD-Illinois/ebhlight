/******************************************************************************
 *                                                                            *
 * RAD_UTILS.C                                                                *
 *                                                                            *
 * HELPER FUNCTIONS FOR RADIATION INFRASTRUCTURE                              *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
void init_rad(grid_prim_type Prad)
{
  set_units();

  ZLOOP {
    sim_vol += ggeom[i][j][CENT].g*dx[1]*dx[2]*dx[3]*L_unit*L_unit*L_unit;
  }
  sim_vol = mpi_reduce(sim_vol);

  init_emissivity();

  #if SCATTERING
  init_hotcross();
  #endif
}

void init_superphoton_resolution()
{
  init_make_statics();
  
  made_tune_proc = abs_tune_proc = scatt_tune_proc = 0;

  if (tune_emiss < 0) {
    fprintf(stdout, "tune_emiss < 0 initially: guessing what tune_emiss should be\n");
    // Estimate dnz
    tune_emiss = 1.;
    printf("tune_emiss = %e\n", tune_emiss);
    get_dnz(P);
    double dnz_tot = 0;
    ZLOOP {
      dnz_tot += dnzs[i][j][k]*nthreads;
    }
    dnz_tot = mpi_reduce(dnz_tot);
    printf("dnz_tot: %e\n", dnz_tot);
    tune_emiss = nph_per_proc*dt/Rout_rad/dnz_tot;
    init_make_statics();
    set_weight(P);
    printf("tune_emiss = %e\n", tune_emiss);
    get_dnz(P);
    dnz_tot = 0;
    ZLOOP {
      dnz_tot += dnzs[i][j][k]*nthreads;
    }
    dnz_tot = mpi_reduce(dnz_tot);
    printf("dnz_tot: %e\n", dnz_tot);
  }

  #if METRIC == MKS || METRIC == MMKS
  dt_tune_emiss = 0.5;
  dt_tune_scatt = Rout_rad;
  // Reasonable (low) initial guess
  tune_scatt = 1./((Rout_rad*L_unit*THOMSON*Ne_unit)*(16.*pow(10,2)));
  #else
  dt_tune_emiss = tf;
  dt_tune_scatt = MY_MAX(MY_MAX(N1TOT*dx[1], N2TOT*dx[2]), N3TOT*dx[3]);
  dt_tune_scatt = MY_MAX(dt_tune_scatt, tf);
  #endif
  if (t_tune_emiss <= 0.0) t_tune_emiss = dt_tune_emiss;
  if (t_tune_scatt <= 0.0) t_tune_scatt = 2.*dt_tune_scatt;
  #if ELECTRONS
  if (strcmp(init_from_grmhd, "No") == 0 && fel0 < 0.1) {
    t_tune_emiss = 500.;
    t_tune_scatt = 500.;
  }
  #endif
}

void update_superphoton_resolution(grid_prim_type Prad)
{
  double L, real, ideal, correction;
  int made_tune, abs_tune, scatt_tune;

  made_tune_proc += step_made;
  abs_tune_proc += step_abs;
  scatt_tune_proc += step_scatt;

  #if METRIC == MKS || METRIC == MMKS
  L = Rout_rad;
  #else
  L = MY_MAX(MY_MAX(N1TOT*dx[1], N2TOT*dx[2]), N3TOT*dx[3]);
  #endif

  if (t >= t_tune_emiss) {
    made_tune = mpi_reduce_int(made_tune_proc)/mpi_nprocs();
    abs_tune = mpi_reduce_int(abs_tune_proc)/mpi_nprocs();

    real = made_tune - abs_tune;
    ideal = dt_tune_emiss*nph_per_proc/L;
    correction = ideal/real;

    // Limit strength of correction
    if (correction < 0.) {
      correction = 4./3.;
    } else {
      correction = MY_MIN(correction, 4./3.);
      correction = MY_MAX(correction, 3./4.);
    }

    // If no superphotons are being emitted (yet) don't modify emission strength
    //if (real < SMALL) correction = 1.;

    tune_emiss *= correction;

    if (mpi_io_proc()) {
      fprintf(stdout, "Emission correction! tune = %e correction = %e\n",
        tune_emiss, correction);
    }

    t_tune_emiss += dt_tune_emiss;
    set_weight(Prad);
    made_tune_proc = abs_tune_proc = 0;
  }

  scatt_tune = mpi_reduce_int(scatt_tune_proc)/mpi_nprocs();
  real = scatt_tune;
  ideal = dt_tune_scatt*nph_per_proc/L;
  correction = ideal/real;
  if (t >= t_tune_scatt || correction < 0.25) {
    real = scatt_tune;
    ideal = dt_tune_scatt*nph_per_proc/L;
    correction = ideal/real;

    // Limit strength of correction
    correction = MY_MIN(correction, 2.);
    correction = MY_MAX(correction, 0.5);

    // If no superphotons are being emitted (yet) don't modify emission strength
    if (real < SMALL) correction = 1.;

    tune_scatt *= correction;

    if (mpi_io_proc()) {
      fprintf(stdout, "Scattering correction! tune = %e correction = %e\n",
        tune_scatt, correction);
    }

    t_tune_scatt += dt_tune_scatt;
    scatt_tune_proc = 0;
  }
}

double linear_interp_log(double x, double *table, double lx_min, double dlx)
{
  double lx = log(x);
  double dn = (lx - lx_min)/dlx;
  int n = (int)dn;
  dn = dn - n;

  return (1. - dn)*table[n] + dn*table[n+1];
}

void set_units()
{
  #if METRIC == MKS || METRIC == MMKS
  L_unit = GNEWT*Mbh/(CL*CL);
  #endif
  T_unit = L_unit/CL;
  RHO_unit = M_unit*pow(L_unit,-3.);
  U_unit = RHO_unit*CL*CL;
  B_unit = CL*sqrt(4.*M_PI*RHO_unit);
  Ne_unit = RHO_unit/(MP + ME);

  #if ELECTRONS
  Thetae_unit = MP/ME;
  #else
  Thetae_unit = (gam-1.)*MP/ME/(1. + tp_over_te);
  #endif
  kphys_to_num = ME/M_unit;
}

// Remove superphoton from list and release memory
void list_remove(struct of_photon **ph, struct of_photon **ph_head,
  struct of_photon **ph_prev)
{
  if(*ph_prev != NULL)
  {
    (*ph_prev)->next = (*ph)->next;
    free(*ph);
    *ph = (*ph_prev)->next;
  } else {
    *ph_head = (*ph)->next;
    free(*ph);
    *ph = *ph_head;
  }
}

double get_Thetae(double Prad[NVAR])
{
  double Thetae;
  #if ELECTRONS
  Thetae = Prad[KEL]*pow(Prad[RHO],game-1.)*Thetae_unit;
  #else
  Thetae = Prad[UU]/Prad[RHO]*Thetae_unit;
  #endif

  return MY_MIN(Thetae, thetae_max);
}

void get_fluid_zone(int i, int j, int k, grid_prim_type Prad,
  struct of_microphysics *m,
  double Ucon[NDIM], double Ucov[NDIM],
  double Bcon[NDIM], double Bcov[NDIM])
{
  double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;

  m->Ne = Prad[i][j][k][RHO]*Ne_unit;
  m->Thetae = get_Thetae(Prad[i][j][k]);

  Bp[1] = Prad[i][j][k][B1]*B_unit;
  Bp[2] = Prad[i][j][k][B2]*B_unit;
  Bp[3] = Prad[i][j][k][B3]*B_unit;

  Vcon[1] = Prad[i][j][k][U1];
  Vcon[2] = Prad[i][j][k][U2];
  Vcon[3] = Prad[i][j][k][U3];

  // Get Ucov
  VdotV = 0.;
  for(int l = 1; l < NDIM; l++) {
    for(int m = 1; m < NDIM; m++) {
      VdotV += ggeom[i][j][CENT].gcov[l][m]*Vcon[l]*Vcon[m];
    }
  }
  Vfac = sqrt(-1./ggeom[i][j][CENT].gcon[0][0]*(1. + fabs(VdotV)));
  Ucon[0] = -Vfac*ggeom[i][j][CENT].gcon[0][0];
  for(int l = 1; l < NDIM; l++)
    Ucon[l] = Vcon[l] - Vfac*ggeom[i][j][CENT].gcon[0][l];
  lower(Ucon, ggeom[i][j][CENT].gcov, Ucov);

  // Get Bcon, Bcov, and B
  UdotBp = 0.;
  for(int l = 1; l < NDIM; l++)
    UdotBp += Ucov[l]*Bp[l];
  Bcon[0] = UdotBp;
  for(int l = 1; l < NDIM; l++)
    Bcon[l] = (Bp[l] + Ucon[l]*UdotBp)/Ucon[0];
  lower(Bcon, ggeom[i][j][CENT].gcov, Bcov);
  m->B = sqrt(Bcon[0]*Bcov[0] + Bcon[1]*Bcov[1] + Bcon[2]*Bcov[2]
	      + Bcon[3]*Bcov[3]);

  // Prevent highly magnetized regions from emitting due to bad internal energy
  double sigma = pow(m->B/B_unit,2.)/(m->Ne/Ne_unit);
  if (sigma > sigma_max) {
    m->Thetae = SMALL;
  }
}

int is_null(double Kcov[NDIM], double Kcon[NDIM], double K0, double KdotKprev,
  double *KdotK)
{
  *KdotK = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    *KdotK += Kcov[mu]*Kcon[mu];
  }

  double K0sqr = pow(K0, 2.);

  if (fabs(*KdotK - KdotKprev)/K0sqr < kdotk_tol) {
    return 1;
  } else {
    return 0;
  }
}

void Xtoijk(double X[NDIM], int *i, int *j, int *k)
{
  *i = (X[1] - startx[1])/dx[1] + NG - global_start[1];
  *j = (X[2] - startx[2])/dx[2] + NG - global_start[2];
  *k = (X[3] - startx[3])/dx[3] + NG - global_start[3];
}

void copy_photon(struct of_photon *ph, struct of_photon *phc)
{
  for (int mu = 0; mu < NDIM; mu++) {
    for (int n = 0; n < NSUP; n++) {
      phc->X[n][mu] = ph->X[n][mu];
      phc->Kcov[n][mu] = ph->Kcov[n][mu];
      phc->Kcon[n][mu] = ph->Kcon[n][mu];
    }
    phc->origin[mu] = ph->origin[mu];
  }
  phc->w = ph->w;
  phc->KdotKprev = ph->KdotKprev;
  phc->nscatt = ph->nscatt;
  phc->t0 = ph->t0;
}

void print_ph_diag(struct of_photon *ph)
{
  printf(" --- PHOTON STATE --- \n");
  for (int n = 0; n < 3; n++) {
    for (int mu = 0; mu < NDIM; mu++) {
      printf("[%i][%i] X = %e Kcov = %e Kcon = %e\n", n, mu, ph->X[n][mu],
        ph->Kcov[n][mu], ph->Kcon[n][mu]);
    }
    double r, th;
    bl_coord(ph->X[n], &r, &th);
    printf("r, th = %e %e\n", r, th);
  }
  printf("origin = %i %i %i %i\n", ph->origin[0], ph->origin[1], ph->origin[2],
    ph->origin[3]);
  printf("w = %e\n", ph->w);
  printf("K.Kprev = %e\n", ph->KdotKprev);
  printf("nscatt = %i\n", ph->nscatt);
  printf("t0 = %e\n", ph->t0);
}

// Use second-order interpolation to get X^{\mu}, K_{\mu} at time t_interp
int get_X_K_interp(struct of_photon *ph, double t_interp, double X[NDIM],
  double Kcov[NDIM], double Kcon[NDIM])
{
  if (t_interp == ph->X[2][0]) {
    for (int mu = 0; mu < NDIM; mu++) {
      X[mu] = ph->X[2][mu];
      Kcov[mu] = ph->Kcov[2][mu];
      Kcon[mu] = ph->Kcon[2][mu];
    }
    return SPH_INTERP_SUCCESS;
  }

  double KdotKprev;
  KdotKprev = ph->KdotKprev;
  double c1 = (t_interp - ph->X[2][0])*(t_interp - ph->X[1][0])/
             ((ph->X[0][0] - ph->X[2][0])*(ph->X[0][0] - ph->X[1][0]));
  double c2 = (t_interp - ph->X[0][0])*(t_interp - ph->X[2][0])/
             ((ph->X[1][0] - ph->X[0][0])*(ph->X[1][0] - ph->X[2][0]));
  double c3 = (t_interp - ph->X[0][0])*(t_interp - ph->X[1][0])/
             ((ph->X[2][0] - ph->X[0][0])*(ph->X[2][0] - ph->X[1][0]));
  DLOOP1 {
    X[mu] = c1*ph->X[0][mu] + c2*ph->X[1][mu] + c3*ph->X[2][mu];
    Kcov[mu] = c1*ph->Kcov[0][mu] + c2*ph->Kcov[1][mu] + c3*ph->Kcov[2][mu];
    Kcon[mu] = c1*ph->Kcon[0][mu] + c2*ph->Kcon[1][mu] + c3*ph->Kcon[2][mu];
  }

  double kdotk = dot(Kcon, Kcov);
  int status;
  
  // If interpolation is terrible, we need to replace with an honest push
  // Also do this if there aren't yet three support points
  if (fabs((kdotk-KdotKprev)/(Kcov[0]*Kcov[0])) > 1e-3 || ph->X[0][0] < 0) {
    if (t_interp < (1. - 1.e-50)*ph->X[1][0] && ph->X[0][0] > 0) {
      for (int mu = 0; mu < NDIM; mu++) {
        X[mu] = ph->X[0][mu];
        Kcov[mu] = ph->Kcov[0][mu];
        Kcon[mu] = ph->Kcon[0][mu];
      }
      double KdotKprev = dot(Kcov, Kcon);
      status = push_X_K(X, Kcov, Kcon, KdotKprev, 
        t_interp-ph->X[0][0]); 
    } else if (ph->X[1][0] > 0) {
      for (int mu = 0; mu < NDIM; mu++) {
        X[mu] = ph->X[1][mu];
        Kcov[mu] = ph->Kcov[1][mu];
        Kcon[mu] = ph->Kcon[1][mu];
      }
      double KdotKprev = dot(Kcov, Kcon);
      status = push_X_K(X, Kcov, Kcon, KdotKprev, 
        t_interp-ph->X[1][0]);
    } else {
      for (int mu = 0; mu < NDIM; mu++) {
        X[mu] = ph->X[2][mu];
        Kcov[mu] = ph->Kcov[2][mu];
        Kcon[mu] = ph->Kcon[2][mu];
      }
      double KdotKprev = dot(Kcov, Kcon);
      status = push_X_K(X, Kcov, Kcon, KdotKprev, 
        t_interp-ph->X[2][0]);
    }
    if (status == PUSH_FAIL) {
      fprintf(stderr, "get_X_K_interp failed!\n");
      for (int n = 0; n < 3; n++) {
        printf("ph->X[%i][] = %e %e %e %e\n", n, ph->X[n][0], ph->X[n][1], ph->X[n][2], ph->X[n][3]);
      }
      fprintf(stderr, "X[] = %e %e %e %e\n", X[0], X[1], X[2], X[3]);
      return SPH_INTERP_FAIL;
    }
  }

  return SPH_INTERP_SUCCESS;
}

int push_to_X_K(double t, struct of_photon *ph, double X[NDIM], 
  double Kcov[NDIM], double Kcon[NDIM])
{
  int status;
  if (t < (1. - 1.e-50)*ph->X[1][0]) {
    for (int mu = 0; mu < NDIM; mu++) {
      X[mu] = ph->X[0][mu];
      Kcov[mu] = ph->Kcov[0][mu];
      Kcon[mu] = ph->Kcon[0][mu];
    }
    double KdotKprev = dot(Kcov, Kcon);
    status = push_X_K(X, Kcov, Kcon, KdotKprev, t - ph->X[0][0]); 
  } else {
    for (int mu = 0; mu < NDIM; mu++) {
      X[mu] = ph->X[1][mu];
      Kcov[mu] = ph->Kcov[1][mu];
      Kcon[mu] = ph->Kcon[1][mu];
    }
    double KdotKprev = dot(Kcov, Kcon);
    status = push_X_K(X, Kcov, Kcon, KdotKprev, t - ph->X[1][0]);
  }

  return status;
}

// Does superphoton need to be pushed from step t to t + dt?
int to_be_pushed(double t, double dt, struct of_photon *ph)
{
  int i, j, k;
  Xtoijk(ph->X[2], &i, &j, &k);

  if (ph->X[2][0] < t + dt) {
    return 1;
  } else {
    return 0;
  }
}

// Move photon from donor list to head of recipient list; advance donor list
void swap_ph(struct of_photon **donor, struct of_photon **recipient)
{
  struct of_photon *tmp;
  if (*recipient == NULL) {
    *recipient = *donor;
    *donor = (*donor)->next;
    (*recipient)->next = NULL;
  } else {
    tmp = *donor;
    *donor = (*donor)->next;
    tmp->next = *recipient;
    *recipient = tmp;
  }
}

// These atomic calls can be very slow on certain problems
void set_Rmunu() 
{
  memset((void*)Rmunu, 0,
    (N1+2*NG)*(N2+2*NG)*(N3+2*NG)*NDIM*NDIM*sizeof(double));
  memset((void*)Nsph, 0, (N1+2*NG)*(N2+2*NG)*(N3+2*NG)*sizeof(int));
  memset((void*)nph, 0, (N1+2*NG)*(N2+2*NG)*(N3+2*NG)*sizeof(double));
  #pragma omp parallel
	{
		struct of_photon *ph = photon_lists[omp_get_thread_num()];
		double X[NDIM], Kcov[NDIM], Kcon[NDIM];
		while (ph != NULL) {
			int i, j, k;
			get_X_K_interp(ph, t, X, Kcov, Kcon);
			Xtoijk(X, &i, &j, &k);

			double volume = ggeom[i][j][CENT].g*dx[1]*dx[2]*dx[3];

			for (int mu = 0; mu < NDIM; mu++) {
				for (int nu = 0; nu < NDIM; nu++) {
					#pragma omp atomic
					Rmunu[i][j][k][mu][nu] += kphys_to_num*Kcon[mu]*Kcov[nu]*
																		ph->w/(Kcon[0]*volume);
				}
			}
      #pragma omp atomic
      Nsph[i][j][k] += 1;

      #pragma omp atomic
      nph[i][j][k] += ph->w/(volume*pow(L_unit,3));

			ph = ph->next;
		}
	} // omp parallel
}

void get_nuLnu_bin(double X[NDIM], int *thbin, int *phibin)
{
  double r, th, phi;
  bl_coord(X, &r, &th); 
  phi = fmod(X[3], 2.*M_PI);
  //phi = X[3] % (2.*M_PI);

  double dth = M_PI/NTH;
  double dphi = 2.*M_PI/NPHI;

  *thbin = (int)(phi/dphi);
  *phibin = (int)(th/dth);
}
#endif // RADIATION

