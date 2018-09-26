/******************************************************************************
 *                                                                            *
 * DECS.H                                                                     *
 *                                                                            *
 * GLOBAL MACROS, FUNCTION DEFINITIONS, INCLUDES, AND DECLARATIONS            *
 *                                                                            *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

#include "constants.h"
#include "params.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923132169164
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.4142135623730950488016887242
#endif

/*******************************************************************************
      COMPILE-TIME PARAMETERS :
*******************************************************************************/

// Number of active zones on each MPI process
#define N1       (N1TOT/N1CPU)
#define N2       (N2TOT/N2CPU)
#define N3       (N3TOT/N3CPU)

// Max size for 1D slice is NMAX
#define N12      (N1 > N2 ? N1 : N2)
#define NMAX     (N12 > N3 ? N12 : N3)

#define NDIM       (4)    // Number of total dimensions
#define NPG        (5)    // Number of positions on grid for grid functions
#define NG         (3)    // Number of ghost zones

// Fixup parameters
#define RHOMINLIMIT (1.e-20)
#define UUMINLIMIT  (1.e-20)
#define RHOMIN  (1.e-5)
#define UUMIN (1.e-7)

// Numerical convenience to represent a small (<< 1) non-zero quantity
#define SMALL (1.e-20)

// Maximum value of gamma, the Lorentz factor
#define GAMMAMAX (50.)

// Maximum fractional increase in timestep per timestep
#define SAFE  (1.3)

// Superphoton diagnostics
#define MAXNSCATT (3)

// Whether to move polar axis slightly off of coordinate singularity
#define COORDSINGFIX 1
#define SINGSMALL (1.E-20)

// Whether to record a representative sample of superphoton positions
#define TRACK_PH (0)

// I/O format strings
#define FMT_DBL_OUT "%28.18e"
#define FMT_INT_OUT "%10d"
#define STRLEN (2048)

// Reconstruction algorithms
#define LINEAR (0)
#define WENO   (2)

// Primitive and conserved variables
#define RHO (0)
#define UU  (1)
#define U1  (2)
#define U2  (3)
#define U3  (4)
#define B1  (5)
#define B2  (6)
#define B3  (7)
#define NVAR_BASE (B3 + 1)
#if ELECTRONS
#define KEL  (NVAR_BASE - 1 + 1)
#define KTOT (NVAR_BASE - 1 + 2)
#define NVAR_EL (2)
#else
#define NVAR_EL (0)
#endif
#define NVAR (NVAR_BASE + NVAR_EL)

// Centering of grid functions
#define FACE1 (0)
#define FACE2 (1)
#define FACE3 (2)
#define CORN  (3)
#define CENT  (4)
#define FACESTART (FACE1)
#define FACEEND (FACE3 + 1)
#define PGSTART (FACE1)
#define PGEND (NPG)

// Slope limiter
#define MC   (0)
#define VANL (1)
#define MINM (2)

// Fluid and radiation boundaries
#define BC_OUTFLOW  (0)
#define BC_PERIODIC (1)
#define BC_POLAR    (2)
#define BC_PROB     (3)
#define BC_ESCAPE   (4)
#define BC_CAMERA   (5)
#define BC_EQUILIB  (6)
#define BC_REFLECT  (7)

// Metric
#define MINKOWSKI (0) // Historial reasons
#define CARTESIAN (0)
#define MKS       (1)
#define MMKS      (2)

// Diagnostic calls
#define DIAG_INIT  (0)
#define DIAG_DUMP  (1)
#define DIAG_LOG   (2)
#define DIAG_FINAL (3)

// Types of restarts
#define RESTART_TEMP (0)
#define RESTART_PERM (1)

// Failure modes
#define FAIL_UTOPRIM     (0)
#define FAIL_VCHAR_DISCR (1)
#define FAIL_COEFF_NEG   (2)
#define FAIL_COEFF_SUP   (3)
#define FAIL_GAMMA       (4)
#define FAIL_METRIC      (5)

// Geodesic integration and interpolation
#define PUSH_FAIL    (0)
#define PUSH_SUCCESS (1)
#define SPH_INTERP_FAIL    (0)
#define SPH_INTERP_SUCCESS (1)

// Root finding
#define ROOT_SUCCESS (1)
#define ROOT_FAIL    (0)
#define FCOUNT_NBINS (6)
#define FCOUNT_MORE  (FCOUNT_NBINS-1)

// Timers
#define TIMER_UPDATE   (0)
#define TIMER_FLUXCALC (1)
#define TIMER_FIXUP    (2)
#define TIMER_BOUND    (3)
#define TIMER_DIAG     (4)
#define TIMER_OUT      (5)
#define TIMER_ELECTRON (6)
#define TIMER_MAKE     (7)
#define TIMER_PUSH     (8)
#define TIMER_INTERACT (9)
#define TIMER_ALL      (10)
#define NUM_TIMERS     (11)

/*******************************************************************************
    GLOBAL ARRAYS
*******************************************************************************/
typedef double grid_prim_type[N1 + 2*NG][N2 + 2*NG][N3 + 2*NG][NVAR];
typedef double grid_double_type[N1 + 2*NG][N2 + 2*NG][N3 + 2*NG];
typedef double grid_fourvector_type[N1 + 2*NG][N2 + 2*NG][N3 + 2*NG][NDIM];
typedef double grid_tensor_type[N1 + 2*NG][N2 + 2*NG][N3 + 2*NG][NDIM][NDIM];
typedef int grid_int_type[N1 + 2*NG][N2 + 2*NG][N3 + 2*NG];

extern grid_prim_type P;        // Primitive variables
extern grid_prim_type F1;       // X1 fluxes
extern grid_prim_type F2;       // X2 fluxes
extern grid_prim_type F3;       // X3 fluxes
extern grid_prim_type Ph;       // Half-step primitives
extern grid_prim_type Psave;       // Half-step primitives
extern grid_int_type pflag;     // Failure points
extern grid_int_type fail_save;
extern grid_fourvector_type jcon;
#if RADIATION
extern grid_fourvector_type radG; // Radiation four-force
extern grid_fourvector_type radG_prev; // Radiation four-force
extern grid_fourvector_type radG_buf;
extern grid_tensor_type Rmunu;    // Radiation stress-energy tensor
extern grid_int_type Nsph;
extern grid_double_type nph;
extern struct of_photon **photon_lists;
extern struct of_photon **photon_mpi_lists;
extern double nuLnu[MAXNSCATT+1][NTH][NPHI][NU_BINS_SPEC];
//extern double dOmega[N2+2*NG][N3+2*NG];
extern double Jrad[MAXNSCATT+2][N1+2*NG][N2+2*NG][N3+2*NG];
extern double Jrad_buf[MAXNSCATT+2][N1+2*NG][N2+2*NG][N3+2*NG];
extern grid_int_type Nem, Nabs, Nsc;
extern grid_int_type Nsuper;
extern grid_double_type Esuper;
extern grid_prim_type psupersave;
#endif
#if ELECTRONS
extern grid_double_type Qvisc_e, Qvisc_p, Qcoul;
#endif // ELECTRONS

/*******************************************************************************
    GLOBAL VARIABLES SECTION
*******************************************************************************/
// Command line arguments
extern char outputdir[STRLEN], dumpdir[STRLEN], restartdir[STRLEN];
extern char xmfdir[STRLEN];
extern char init_from_grmhd[STRLEN];
extern char metric[STRLEN], reconstruction[STRLEN];

// Physics parameters
extern double a;
extern double gam;
extern double M_unit;
extern double Reh;
extern double Risco;
#if RADIATION
extern double mbh, Mbh, L_unit, T_unit, M_unit, RHO_unit, U_unit, B_unit;
extern double Ne_unit, Thetae_unit, kphys_to_num;
extern double tp_over_te, thetae_max, sigma_max, kdotk_tol;
#endif

// Numerical parameters
extern double Rin, Rout, Rout_vis, hslope;
extern double poly_norm, poly_xt, poly_alpha, mks_smooth;
#if RADIATION
extern double Rout_rad;
extern double nph_per_proc;
extern double tune_emiss, t_tune_emiss, dt_tune_emiss;
extern double tune_scatt, t_tune_scatt, dt_tune_scatt;
extern int made_tune_proc, abs_tune_proc, scatt_tune_proc;
extern double numin_emiss, numax_emiss, numin_spec, numax_spec;
extern double kappa;
extern double startx_rad[NDIM], stopx_rad[NDIM];
extern double wgtC;
extern int step_made, step_abs, step_scatt, step_lost, step_rec, step_tot;
extern int step_sent, step_rcvd, step_fail;
extern int step_made_all, step_abs_all, step_scatt_all, step_lost_all;
extern int step_rec_all, step_sent_all, step_rcvd_all, step_tot_all;
extern int step_fail_all;
//extern int thbin, phibin;
extern double Nph_to_track;
extern double sim_vol;
#endif
#if ELECTRONS
extern double tptemin, tptemax;
#endif
extern double cour;
extern double dV, dx[NDIM], startx[NDIM], stopx[NDIM], startx_proc[NDIM], stopx_proc[NDIM];
extern double x1Min, x1Max, x2Min, x2Max, x3Min, x3Max;
extern double dt, dtsave;
extern double t, tf;
extern int nstep;
extern int is_restart;
extern double dnzs[N1+2*NG][N2+2*NG][N3+2*NG];

// Output parameters
extern double DTd;
extern double DTl;
extern double DTr;
extern int DNr;
extern int DTp;
extern int DTf;
extern double DTw;
extern int dump_cnt;
extern int rdump_cnt;
extern double tdump, trestart, tlog;
extern int root_fcount[FCOUNT_NBINS];

// Global flags
extern int failed;
extern int lim;

// Diagnostics
extern double mdot, mdot_eh;
extern double edot, edot_eh;
extern double ldot, ldot_eh;
extern int icurr, jcurr, kcurr;

// Parallelism
extern int nthreads;

// Electrons
#if ELECTRONS
extern double game, gamp;
extern double fel0;
#endif

// Set global variables that indicate current local metric, etc.
struct of_geom {
  double gcon[NDIM][NDIM];
  double gcov[NDIM][NDIM];
  double g;
  double alpha;
};

struct of_state {
  double ucon[NDIM];
  double ucov[NDIM];
  double bcon[NDIM];
  double bcov[NDIM];
};

#if RADIATION
// WARNING: if you change struct_of_photon, be sure to change
// the the photon hdf5 types in io.c.
#define NSUP 3
struct of_photon {
  // NSUP >=3 X^{\mu}, K^{\mu}, K_{\mu} so photon data always available anywhere
  // between n and n+1
  double X[NSUP][NDIM];
  double Kcov[NSUP][NDIM];
  double Kcon[NSUP][NDIM];
  double w;
  double KdotKprev;
  int nscatt;
  int origin[NDIM];
  double t0;
  int is_tracked;
  struct of_photon *next;
};

#define PH_ELEM (10)
struct of_track_photon {
  double X1;
  double X2;
  double X3;
  int nscatt;
};
struct of_microphysics {
  double Thetae;
  double Ne;
  double B;
};
#endif // RADIATION

// More grid functions. Axisymmetry assumed.
extern double conn[N1 + 2*NG][N2 + 2*NG][NDIM][NDIM][NDIM];
extern struct of_geom ggeom[N1+2*NG][N2+2*NG][NPG];
#if RADIATION
//extern double dt_light, dt_light_min;
extern double dt_light[N1+2*NG][N2+2*NG], dt_light_min;
#endif

// MPI-specific stuff
extern int global_start[NDIM];
extern int global_stop[NDIM];

/*******************************************************************************
    MACROS
*******************************************************************************/
#define ILOOP \
  for (int i = 0 + NG; i < N1 + NG; i++)
#define ILOOPALL \
  for (int i = 0; i < N1 + 2*NG; i++)
#define JLOOP \
  for (int j = 0 + NG; j < N2 + NG; j++)
#define JLOOPALL \
  for (int j = 0; j < N2 + 2*NG; j++)
#define KLOOP \
  for (int k = 0 + NG; k < N3 + NG; k++)
#define KLOOPALL \
  for (int k = 0; k < N3 + 2*NG; k++)
#define ZLOOP \
  for (int i = 0 + NG; i < N1 + NG; i++) \
  for (int j = 0 + NG; j < N2 + NG; j++) \
  for (int k = 0 + NG; k < N3 + NG; k++)
#define ZLOOPALL \
  ILOOPALL JLOOPALL KLOOPALL
#define ISLOOP(istart,istop) \
  for (int i = istart + NG; i <= istop + NG; i++)
#define JSLOOP(jstart,jstop) \
  for (int j = jstart + NG; j <= jstop + NG; j++)
#define KSLOOP(kstart,kstop) \
  for (int k = kstart + NG; k <= kstop + NG; k++)
#define ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) \
  for (int i = istart + NG; i <= istop + NG; i++) \
  for (int j = jstart + NG; j <= jstop + NG; j++) \
  for (int k = kstart + NG; k <= kstop + NG; k++)
// Loop over faces
#define FACELOOP for (int face = FACESTART; face < FACEEND; face++)
// Loop over all locations
#define LOCLOOP for (int loc = PGSTART; loc < PGEND; loc++)
// Loop over primitive variables
#define PLOOP for(int ip = 0; ip < NVAR; ip++)
#define BASELOOP for (int ip = 0; ip < NVAR_BASE; ip++)

// Loop over spacetime indices
#define DLOOP1 for (int mu = 0; mu < NDIM; mu++)
#define DLOOP2 for (int mu = 0; mu < NDIM; mu++) \
               for (int nu = 0; nu < NDIM; nu++)

#define MY_MIN(fval1,fval2) ( ((fval1) < (fval2)) ? (fval1) : (fval2))
#define MY_MAX(fval1,fval2) ( ((fval1) > (fval2)) ? (fval1) : (fval2))
#define MY_SIGN(fval) ( ((fval) <0.) ? -1. : 1. )

#define delta(i,j) ((i == j) ? 1. : 0.)
#define dot(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3])

/*******************************************************************************
    FUNCTION DECLARATIONS
*******************************************************************************/

// bounds.c
void bound_prim(grid_prim_type prim);
void fix_flux(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3);
#if RADIATION
void bound_superphotons(double t, double dt);
int rad_error_check(struct of_photon **ph);
void bound_rad_transport(double X[NDIM], struct of_photon *ph, int is_transported);
int bound_rad_isactive(double X[NDIM], struct of_photon *ph);
int rad_mpi_transport(struct of_photon **ph, struct of_photon **prev, struct
  of_photon **head, double X[NDIM], int active);
#endif

// coord.c
void coord(int i, int j, int k, int loc, double X[NDIM]);
double r_of_X(const double X[NDIM]);
double th_of_X(const double X[NDIM]);
void jac_harm_to_bl(const double X[NDIM],
  double Jcov[NDIM][NDIM], double Jcon[NDIM][NDIM]);
void jac_bl_to_cart(const double X[NDIM],
  double Jcov[NDIM][NDIM], double Jcon[NDIM][NDIM]);
void jac_harm_to_cart(const double X[NDIM],
  double Jcov[NDIM][NDIM], double Jcon[NDIM][NDIM]);
void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM]);
void bl_coord(const double X[NDIM], double *r, double *th);
int bl_i_of_r(double r);
void cart_coord(const double X[NDIM], double Xcart[NDIM]);
void set_gcov(double X[NDIM], double gcov[NDIM][NDIM]);
void set_points();
void zero_arrays(void);
void set_grid(void);

// current.c
void current_calc();

// diag.c
void reset_log_variables();
void reset_dump_variables();
void diag(int call_code);
void fail(int fail_type);
void area_map(int i, int j, int k, grid_prim_type prim);
void diag_flux(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3);
double flux_ct_divb(int i, int j, int k);
#if RADIATION
void record_superphoton(double X[NDIM], struct of_photon *ph);
#endif

// electrons.c
#if ELECTRONS
void init_electrons();
void heat_electrons(grid_prim_type Pi, grid_prim_type Ps, grid_prim_type Pf, 
  double Dt);
double get_fel(int i, int j, int k, double p[NVAR]);
#if RADIATION
void coulomb(grid_prim_type Pi, grid_prim_type Ps, grid_prim_type Pf, 
  double Dt);
void apply_rad_force_e(grid_prim_type Prh, grid_prim_type Pr,
  grid_fourvector_type radG, double Dt);
#endif // RADIATION
void fixup_electrons(grid_prim_type p);
#endif

// emissivity.c
#if RADIATION
double jnu(double nu, const struct of_microphysics *m, double theta);
double Jnu(double nu, const struct of_microphysics *m);
double get_J(double Ne, double Thetae, double Bmag);
void init_emissivity();
#endif

// fixup.c
void fixup(grid_prim_type pv);
void fixup1zone(int i, int j, int k, double pv[NVAR]);
void fixup_utoprim(grid_prim_type pv);

// interact.c
#if RADIATION
void interact(grid_prim_type P, double t, double dt);
#endif

// input.c
void init_params(char *pfname);
void set_param(char *key, void *data);

// io.c
//void set_core_params();
//void set_param(char *key, void *data);
//void read_params(char *pfname);
void init_io();
void init_fluid_restart();
#if RADIATION
void track_ph();
#endif
void dump_grid();
void dump();
void restart_write(int restart_type);
void restart_read(char *fname);
int restart_init();

// make_superphotons.c
#if RADIATION
void init_make_statics();
void make_superphotons(grid_prim_type Prad, double t, double dt);
void set_weight(grid_prim_type Prad);
void get_dnz(grid_prim_type Prad);
#endif

// metric.c
double gcon_func(double lgcov[][NDIM], double lgcon[][NDIM]);
void conn_func(double *X, struct of_geom *geom, double conn[][NDIM][NDIM]);
void lower(double ucon[NDIM], double gcov[NDIM][NDIM], double ucov[NDIM]);
void raise(double ucov[NDIM], double gcon[NDIM][NDIM], double ucon[NDIM]);
struct of_geom *get_geometry(int ii, int jj, int kk, int loc);
void blgset(int i, int j, struct of_geom *geom);
double bl_gdet_func(double r, double th);
void bl_gcov_func(double r, double th, double gcov[][NDIM]);
void bl_gcon_func(double r, double th, double gcon[][NDIM]);
double MINOR(double m[16], int r0, int r1, int r2, int c0, int c1, int c2);
void adjoint(double m[16], double adjOut[16]);
double determinant(double m[16]);
double invert(double *m, double *invOut);

// mpi.c
void init_mpi();
void sync_mpi_boundaries_X1L(grid_prim_type Pr);
void sync_mpi_boundaries_X1R(grid_prim_type Pr);
void sync_mpi_boundaries_X2L(grid_prim_type Pr);
void sync_mpi_boundaries_X2R(grid_prim_type Pr);
void sync_mpi_boundaries_X3L(grid_prim_type Pr);
void sync_mpi_boundaries_X3R(grid_prim_type Pr);
#if RADIATION
void sync_radG();
void sync_Jrad();
void sync_mpi_photons(struct of_photon **ph_mpi, double t, double dt);
void mpi_reduce_nuLnu();
#endif
int mpi_nprocs();
double mpi_max(double f);
double mpi_min(double f);
double mpi_reduce(double f);
int mpi_reduce_int(int f);
int mpi_io_proc();
void mpi_int_broadcast(int *val);
void mpi_int_broadcast_array(int *val, int size);
void mpi_int_broadcast_proc(int *val, int root);
void mpi_dbl_broadcast(double *val);
double mpi_io_reduce(double val);
double mpi_io_max(double val);
int mpi_myrank();
void mpi_sync_output();
void mpi_barrier();
int mpi_is_periodic(int dir);

// phys.c
void primtoflux(double *pr, struct of_state *q, int dir, struct of_geom *geom,
  double *flux);
void bcon_calc(double *pr, double *ucon, double *ucov, double *bcon);
void mhd_calc(double *pr, int dir, struct of_state *q, double *mhd);
void source(double *ph, struct of_geom *geom, int ii, int jj, double *dU,
  double Dt);
double bsq_calc(double *pr, struct of_geom *geom);
void get_state(double *pr, struct of_geom *geom, struct of_state *q);
void ucon_calc(double *pr, struct of_geom *geom, double *ucon);
int mhd_gamma_calc(double *pr, struct of_geom *geom, double *gamma);
void mhd_vchar(double *pr, struct of_state *q, struct of_geom *geom, int js,
  double *vmax, double *vmin);

// problem.c
void set_problem_params();
void init_prob();
void bound_gas_prob_x1l(int i, int j, int k, grid_prim_type P);
void bound_gas_prob_x1r(int i, int j, int k, grid_prim_type P);
void bound_gas_prob_x2l(int i, int j, int k, grid_prim_type P);
void bound_gas_prob_x2r(int i, int j, int k, grid_prim_type P);
void bound_gas_prob_x3l(int i, int j, int k, grid_prim_type P);
void bound_gas_prob_x4r(int i, int j, int k, grid_prim_type P);

// push_superphotons.c
#if RADIATION
int push_X_K(double X[NDIM], double Kcov[NDIM], double Kcon[NDIM],
  double KdotKprev, double dtpush);
int push_superphoton(struct of_photon *ph, double dtpush);
void push_superphotons(double dt);
#endif

// rad_utils.c
#if RADIATION
void init_rad(grid_prim_type Prad);
void init_superphoton_resolution();
void update_superphoton_resolution(grid_prim_type Prad);
double linear_interp_log(double x, double *table, double lmin, double dl);
void set_units();
void list_remove(struct of_photon **ph, struct of_photon **ph_head,
  struct of_photon **ph_prev);
double get_Thetae(double P[NVAR]);
void get_fluid_zone(int i, int j, int k, grid_prim_type Prad,
  struct of_microphysics *m,
  double Ucon[NDIM], double Ucov[NDIM],
  double Bcon[NDIM], double Bcov[NDIM]);
int is_null(double Kcov[NDIM], double Kcon[NDIM], double K0, double KdotKprev,
  double *KdotK);
void set_Rmunu();
//int is_null(struct of_photon *ph, double *KdotK);
void Xtoijk(double X[NDIM], int *i, int *j, int *k);
void copy_photon(struct of_photon *ph, struct of_photon *phc);
void print_ph_diag(struct of_photon *ph);
int get_X_K_interp(struct of_photon *ph, double t_interp, double X[NDIM],
  double Kcov[NDIM], double Kcon[NDIM]);
int push_to_X_K(double t, struct of_photon *ph, double X[NDIM],                  
  double Kcov[NDIM], double Kcon[NDIM]);
int to_be_pushed(double t, double dt, struct of_photon *ph);
void swap_ph(struct of_photon **donor, struct of_photon **recipient);
void get_nuLnu_bin(double X[NDIM], int *thbin, int *phibin);
#endif

// radiation.c
#if RADIATION
double Bnu_inv(double nu, const struct of_microphysics *m); // TODO?
double jnu_inv(double nu, const struct of_microphysics *m, double theta);
double alpha_inv_scatt(double nu, const struct of_microphysics *m);
double alpha_inv_abs(double nu, const struct of_microphysics *m, double theta);
double get_fluid_nu(double X[NDIM], double Kcov[NDIM], double Ucon[NDIM]);
double get_bk_angle(double X[NDIM], double K[NDIM], double Ucov[NDIM],
  double Bcov[NDIM], double B);
#endif

// random.c
void init_random(int seed);
double get_rand();
double get_chisq(double nu);
void get_ran_dir_3d(double *nx, double *ny, double *nz);

// reconstruction.c
void reconstruct(double ptmp[NMAX+2*NG][NVAR], int N,
  double p_l[NMAX+2*NG][NVAR], double p_r[NMAX+2*NG][NVAR]);

// scatter_superphoton.c
#if RADIATION
int scatter_superphoton(grid_prim_type P, struct of_photon *ph, double X[NDIM],
  double Kcov[NDIM], double Kcon[NDIM]);
void init_hotcross();
double total_compton_cross_lkup(double w, double thetae);
#endif

// step.c
void step();

// tetrads.c
#if RADIATION
void coord_to_tetrad(double Ecov[NDIM][NDIM], double Kcoord[NDIM],
  double Ktetrad[NDIM]);
void tetrad_to_coord(double Econ[NDIM][NDIM], double Ktetrad[NDIM],
  double Kcoord[NDIM]);
void make_tetrad(int i, int j, int k, double Ucon[NDIM], double trial[NDIM],
  double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM], double Ecov[NDIM][NDIM]);
void normalize_null(double Gcov[NDIM][NDIM], double K[NDIM]);
void normalize_null_cov(double Gcon[NDIM][NDIM], double K[NDIM]);
#endif

// timing.c
void time_set(int n, double val);
double time_read(int n);
double get_time_per_step(int timerCode);
void time_init();
void timer_start(int timerCode);
void timer_stop(int timerCode);
void report_performance();

void timers_reset();

// utils.c
void *safe_malloc(int size);
void safe_system(const char *command);
void safe_fscanf(FILE *stream, const char *format, ...);

// utop.c
int Utoprim(double U[NVAR], struct of_geom *geom, double prim[NVAR]);

// xdmf_output.c
void write_xml_file(int dump_id, double t, const char vnams[NVAR][STRLEN]);
