/******************************************************************************
 *                                                                            *
 * DEFS.H                                                                     *
 *                                                                            *
 * GLOBAL VARIABLE DEFINITIONS                                                *
 *                                                                            *
 ******************************************************************************/

/*******************************************************************************
    GLOBAL ARRAYS
*******************************************************************************/
grid_prim_type P;
grid_prim_type F1;
grid_prim_type F2;
grid_prim_type F3;
grid_prim_type Ph;
grid_int_type pflag;
grid_int_type fail_save;
grid_prim_type Psave;
grid_fourvector_type jcon;
#if RADIATION
grid_fourvector_type radG; // Radiation four-force
grid_fourvector_type radG_prev; // Radiation four-force
grid_fourvector_type radG_buf;
grid_tensor_type Rmunu;
grid_int_type Nsph;
grid_double_type nph;
struct of_photon **photon_lists;
struct of_photon **photon_mpi_lists;
double nuLnu[MAXNSCATT+1][NTH][NPHI][NU_BINS_SPEC];
double Jrad[MAXNSCATT+2][N1+2*NG][N2+2*NG][N3+2*NG];
double Jrad_buf[MAXNSCATT+2][N1+2*NG][N2+2*NG][N3+2*NG];
grid_int_type Nem, Nabs, Nsc;
grid_int_type Nsuper;
grid_double_type Esuper;
grid_prim_type psupersave;
#endif
#if ELECTRONS
grid_double_type Qvisc_e, Qvisc_p, Qcoul;
#endif // ELECTRONS

double conn[N1 + 2*NG][N2 + 2*NG][NDIM][NDIM][NDIM];
struct of_geom ggeom[N1+2*NG][N2+2*NG][NPG] ;
#if RADIATION
double dt_light[N1+2*NG][N2+2*NG], dt_light_min;
#endif

/*******************************************************************************
    GLOBAL VARIABLES
*******************************************************************************/
char outputdir[STRLEN], dumpdir[STRLEN], restartdir[STRLEN];
char xmfdir[STRLEN];
char init_from_grmhd[STRLEN];
char metric[STRLEN], reconstruction[STRLEN];

double a;
double gam;
double M_unit;
double Reh;
double Risco;
#if RADIATION
double mbh, Mbh, L_unit, T_unit, M_unit, RHO_unit, U_unit, B_unit;
double Ne_unit, Thetae_unit, kphys_to_num;
double tp_over_te, thetae_max, sigma_max, kdotk_tol;
#endif

double Rin, Rout, Rout_vis, hslope;
double poly_norm, poly_xt, poly_alpha, mks_smooth;
#if RADIATION
double Rout_rad;
double nph_per_proc;
double tune_emiss, t_tune_emiss, dt_tune_emiss;
double tune_scatt, t_tune_scatt, dt_tune_scatt;
int made_tune_proc, abs_tune_proc, scatt_tune_proc;
double numin_emiss, numax_emiss, numin_spec, numax_spec;
double kappa;
double startx_rad[NDIM], stopx_rad[NDIM];
double wgtC;
int step_made, step_abs, step_scatt, step_lost, step_rec, step_tot;
int step_sent, step_rcvd, step_fail;
int step_made_all, step_abs_all, step_scatt_all, step_lost_all;
int step_rec_all, step_sent_all, step_rcvd_all, step_tot_all, step_fail_all;
//int thbin, phibin;
double Nph_to_track;
double sim_vol;
#endif
#if ELECTRONS
double tptemin, tptemax;
#endif
double cour;
double dV, dx[NDIM], startx[NDIM], stopx[NDIM], startx_proc[NDIM], stopx_proc[NDIM];
double x1Min, x1Max, x2Min, x2Max, x3Min, x3Max;
double dt, dtsave;
double t, tf;
double rcurr, hcurr;
int istart, istop, jstart, jstop;
int nstep;
int is_restart;

double DTd;
double DTl;
double DTr;
int DNr;
int DTp;
int DTf;
double DTw;
int dump_cnt;
int rdump_cnt;
double tdump, trestart, tlog;
int root_fcount[FCOUNT_NBINS];

int failed;
int lim;

double mdot = 0., mdot_eh = 0.;
double edot = 0., edot_eh = 0.;
double ldot = 0., ldot_eh = 0.;
int icurr, jcurr, kcurr;

int nthreads;

#if ELECTRONS
double game, gamp;
double fel0;
#endif

int global_start[NDIM];
int global_stop[NDIM];

