/******************************************************************************
 *                                                                            *
 * FIXUP.C                                                                    *
 *                                                                            *
 * REPAIR INTEGRATION FAILURES                                                *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Apply floors to density, internal energy
void fixup(grid_prim_type Pv)
{
  timer_start(TIMER_FIXUP);
  #pragma omp parallel for collapse(3) //schedule(dynamic)
  ZLOOP fixup1zone(i, j, k, Pv[i][j][k]);
  timer_stop(TIMER_FIXUP);
}

void ucon_to_utcon(double ucon[NDIM], struct of_geom *geom, double utcon[NDIM])
{
  double alpha, beta[NDIM], gamma;

  alpha = 1./sqrt(-geom->gcon[0][0]);
  for (int i = 1; i < NDIM; i++) {
    beta[i] = geom->gcon[0][i]*alpha*alpha;
  }
  gamma = alpha*ucon[0];

  utcon[0] = 0.;
  for (int i = 1; i < NDIM; i++) {
    utcon[i] = ucon[i] + gamma*beta[i]/alpha;
  }
}

void ut_calc_3vel(double vcon[NDIM], struct of_geom *geom, double *ut)
{
  double AA, BB, CC, DD, one_over_alpha_sq;

  AA = geom->gcov[0][0];
  BB = 2.*(geom->gcov[0][1]*vcon[1] +
           geom->gcov[0][2]*vcon[2] +
           geom->gcov[0][3]*vcon[3]);
  CC = geom->gcov[1][1]*vcon[1]*vcon[1] +
       geom->gcov[2][2]*vcon[2]*vcon[2] +
       geom->gcov[3][3]*vcon[3]*vcon[3] +
       2.*(geom->gcov[1][2]*vcon[1]*vcon[2] +
           geom->gcov[1][3]*vcon[1]*vcon[3] +
           geom->gcov[2][3]*vcon[2]*vcon[3]);
  DD = 1./(AA + BB + CC);

  one_over_alpha_sq = -geom->gcon[0][0];

  if (DD < one_over_alpha_sq) {
    DD = one_over_alpha_sq;
  }

  *ut = sqrt(DD);
}

#define BSQORHOMAX (50.)
#define BSQOUMAX (2500.)
#define UORHOMAX (50.)
void fixup1zone(int i, int j, int k, double pv[NVAR])
{
  double rhoflr, uflr, f, gamma;
  double rhoscal, uscal;
  struct of_geom *geom;
  struct of_state q;

  #if METRIC == MKS || METRIC == MMKS
  double r, th, X[NDIM];
  coord(i, j, k, CENT, X);
  bl_coord(X, &r, &th);
  rhoscal = pow(r, -2.);
  uscal = pow(rhoscal, gam);
  rhoflr = RHOMIN*rhoscal;
  uflr = UUMIN*uscal;
  #elif METRIC == CARTESIAN
  rhoscal = 1.e-2;
  uscal = 1.e-2;
  rhoflr = RHOMIN*rhoscal;
  uflr = UUMIN*uscal;
  if (rhoflr > pv[RHO]) {
    double drho = rhoflr - pv[RHO];
    pv[RHO] += drho;
  }
  if (uflr > pv[UU]) {
    double du = uflr - pv[UU];
    pv[UU] += du;
  }
  geom = get_geometry(i, j, k, CENT);
  if (mhd_gamma_calc(pv, geom, &gamma)) {
    pflag[i][j][k] = -333;
  } else {
    if (gamma > GAMMAMAX) {
      f = sqrt((GAMMAMAX*GAMMAMAX - 1.)/(gamma*gamma - 1.));
      pv[U1] *= f;
      pv[U2] *= f;
      pv[U3] *= f;
    }
  }
  return;
  #endif

  rhoflr = MY_MAX(rhoflr, RHOMINLIMIT);
  uflr = MY_MAX(uflr, UUMINLIMIT);

  geom = get_geometry(i, j, k, CENT);
  get_state(pv, geom, &q);
  double bsq = dot(q.bcon, q.bcov);
  rhoflr = MY_MAX(rhoflr, bsq/BSQORHOMAX);
  uflr = MY_MAX(uflr, bsq/BSQOUMAX);
  rhoflr = MY_MAX(rhoflr, pv[UU]/UORHOMAX);
  
  if (rhoflr > pv[RHO] || uflr > pv[UU]) {
    double Padd[NVAR], Uadd[NVAR];
    PLOOP Padd[ip] = 0.;
    PLOOP Uadd[ip] = 0.;
    Padd[RHO] = MY_MAX(0., rhoflr - pv[RHO]); 
    Padd[UU] = MY_MAX(0., uflr - pv[UU]);
    
    get_state(Padd, &ggeom[i][j][CENT], &q);

    primtoflux(Padd, &q, 0, &ggeom[i][j][CENT], Uadd);

    double Utot[NVAR];
    get_state(pv, &ggeom[i][j][CENT], &q);
    primtoflux(pv, &q, 0, &ggeom[i][j][CENT], Utot);

    PLOOP Utot[ip] += Uadd[ip];

    PLOOP pv[ip] += Padd[ip];

    // Record fails here?
    Utoprim(Utot, &ggeom[i][j][CENT], pv);
  }

  #if ELECTRONS
  // Reset entropy after floors
  pv[KTOT] = (gam - 1.)*pv[UU]/pow(pv[RHO],gam);
  #endif // ELECTRONS

  // Limit gamma with respect to normal observer
  if (mhd_gamma_calc(pv, geom, &gamma)) {
    pflag[i][j][k] = -333;
  } else {
    if (gamma > GAMMAMAX) {
      f = sqrt((GAMMAMAX*GAMMAMAX - 1.)/(gamma*gamma - 1.));
      pv[U1] *= f;
      pv[U2] *= f;
      pv[U3] *= f;
    }
  }
}

grid_prim_type Pv_tmp;
grid_int_type pflag_tmp;
// Replace bad points with values interpolated from neighbors
#define FLOOP for(int ip=0;ip<B1;ip++)
#if ELECTRONS
double uel[N1+2*NG][N2+2*NG][N3+2*NG];
#endif
void fixup_utoprim(grid_prim_type Pv)
{
  timer_start(TIMER_FIXUP);

  // Flip the logic of the pflag[] so that it now indicates which cells are good
  #pragma omp parallel for collapse(3)
  ZSLOOP(-NG, (N1-1 + NG), -NG, (N2-1 + NG), -NG, (N3-1 + NG)) {
    pflag[i][j][k] = !pflag[i][j][k];
  }
  #pragma omp parallel for collapse(3)
  ZSLOOP(-NG,N1+NG-1,-NG,N2+NG-1,-NG,N3+NG-1) {
    FLOOP Pv_tmp[i][j][k][ip] = Pv[i][j][k][ip];
  }
  #if ELECTRONS
  // Average electron internal energy rather than electron entropy
  ZLOOPALL {
    uel[i][j][k] = 1./(game-1.)*Pv[i][j][k][KEL]*pow(Pv[i][j][k][RHO],game);
  }
  #endif

  // Make sure we are not using ill defined corner regions, while respecting MPI structure
  for(int i = 0; i < NG; i++) {
    for(int j = 0; j < NG; j++) {
      for(int k = 0; k < NG; k++) {
        if (global_start[1] == 0 && global_start[2] == 0 && global_start[3] == 0) pflag[i][j][k] = 0;
        if (global_stop[1] == N1TOT && global_start[2] == 0 && global_start[3] == 0) pflag[i+N1+NG][j][k] = 0;
        if (global_start[1] == 0 && global_stop[2] == N2TOT && global_start[3] == 0) pflag[i][j+N2+NG][k] = 0;
        if (global_start[1] == 0 && global_start[2] == 0 && global_stop[3] == N3TOT) pflag[i][j][k+N3+NG] = 0;
        if (global_stop[1] == N1TOT && global_stop[2] == N2TOT && global_start[3] == 0) pflag[i+N1+NG][j+N2+NG][k] = 0;
        if (global_stop[1] == N1TOT && global_start[2] == 0 && global_stop[3] == N3TOT) pflag[i+N1+NG][j][k+N3+NG] = 0;
        if (global_start[1] == 0 && global_stop[2] == N2TOT && global_stop[3] == N3TOT) pflag[i][j+N2+NG][k+N3+NG] = 0;
        if (global_stop[1] == N1TOT && global_stop[2] == N2TOT && global_stop[3] == N3TOT) pflag[i+N1+NG][j+N2+NG][k+N3+NG] = 0;
      }
    }
  }
  
  #pragma omp parallel for collapse(2)
  ZLOOP {
    if (pflag[i][j][k] == 0) {
      double wsum = 0.;
      double sum[B1] = {0};
      #if ELECTRONS
      double sume = 0.;
      #endif
      
      // Average primitive variables over good neighbors, weighted by distance
      for(int l = -1; l < 2; l++) {
        for(int m = -1; m < 2; m++) {
          for(int n= -1; n < 2; n++) {
            double w = 1./(abs(l)+abs(m)+abs(n) + 1)*pflag[i+l][j+m][k+n];
            wsum += w;
            FLOOP sum[ip] += w*Pv_tmp[i+l][j+m][k+n][ip];
            #if ELECTRONS
            sume += w*uel[i+l][j+m][k+n];
            #endif
          }
        }
      }
      
      // No usable neighbors? Cruder average, and set velocity to zero in normal observer frame
      if (wsum < 1.e-10) {
        wsum = 0.;
        fprintf(stderr,"[%i][istart=%i][%i %i %i] fixup_utoprim problem: No usable neighbors!\n",
            mpi_myrank(), global_start[1], i, j, k);
        FLOOP sum[ip] = 0.;
        #if ELECTRONS
        sume = 0.;
        #endif
        for(int l = -1; l < 2; l++) {
          for(int m = -1; m < 2; m++) {
            for(int n= -1; n < 2; n++) {
              double w = 1./(abs(l)+abs(m)+abs(n) + 1);
              wsum += w;
              FLOOP sum[ip] += w*Pv_tmp[i+l][j+m][k+n][ip];
              #if ELECTRONS
              sume += w*uel[i+l][j+m][k+n];
              #endif
            }
          }
        }
        Pv_tmp[i][j][k][U1] = 0.;
        Pv_tmp[i][j][k][U2] = 0.;
        Pv_tmp[i][j][k][U3] = 0.;
      }
      
      FLOOP Pv[i][j][k][ip] = sum[ip]/wsum;

      #if ELECTRONS
      Pv[i][j][k][KEL] = (gam-1.)*uel[i][j][k]/pow(Pv[i][j][k][RHO],game);
      Pv[i][j][k][KTOT] = (gam-1.)*(Pv[i][j][k][UU])/pow(Pv[i][j][k][RHO],gam);
      #endif

      fixup1zone(i, j, k, Pv[i][j][k]);
    }
  }
  
  timer_stop(TIMER_FIXUP);
}
#undef FLOOP

