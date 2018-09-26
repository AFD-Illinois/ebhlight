/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR KOMISSAROV SHOCKTUBES                               *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

static int shock;
void set_problem_params()
{
  set_param("shock", &shock);
}

void init_prob()
{
	double X[NDIM];

  double rhoL, PL, u1L, u2L, u3L, B1L, B2L, B3L;
  double rhoR, PR, u1R, u2R, u3R, B1R, B2R, B3R;

  if (shock == 0) { // Fast shock
    tf = 2.5;
    cour = 0.4;

    rhoL = 1.;
    PL = 1.;
    u1L = 25.;
    u2L = 0.;
    u3L = 0.;
    B1L = 20.;
    B2L = 25.02;
    B3L = 0.;
    rhoR = 25.48;
    PR = 367.5;
    u1R = 1.091;
    u2R = 0.3923;
    u3R = 0.;
    B1R = 20.;
    B2R = 49.;
    B3R = 0.;
  } else if (shock == 1) { // Slow shock
    tf = 2.;
    
    rhoL = 1.;
    PL = 10.;
    u1L = 1.53;
    u2L = 0.;
    u3L = 0.;
    B1L = 10.;
    B2L = 18.28;
    B3L = 0.;
    rhoR = 3.323;
    PR = 55.36;
    u1R = 0.9571;
    u2R = -0.6822;
    u3R = 0.;
    B1R = 10.;
    B2R = 14.49;
    B3R = 0.;
  } else if (shock == 2) { // Switch-on
    tf = 2.;

    rhoL = 1.78e-3;
    PL = 0.1;
    u1L = -0.765;
    u2L = -1.386;
    u3L = 0.;
    B1L = 1.;
    B2L = 1.022;
    B3L = 0.;
    rhoR = 0.01;
    PR = 1.;
    u1R = 0.;
    u2R = 0.;
    u3R = 0.;
    B1R = 1.;
    B2R = 0.;
    B3R = 0.;
  }  else if (shock == 3) { // Switch-off 
    tf = 1.;

    rhoL = 0.1;
    PL = 1.;
    u1L = -2.;
    u2L = 0.;
    u3L = 0.;
    B1L = 2.;
    B2L = 0.;
    B3L = 0.;
    rhoR = 0.562;
    PR = 10.;
    u1R = -0.212;
    u2R = -0.590;
    u3R = 0.;
    B1R = 2.;
    B2R = 4.71;
    B3R = 0.;
  } else if (shock == 4) { // Shocktube 1
    tf = 1.;

    rhoL = 1.;
    PL = 1000.;
    u1L = 0.;
    u2L = 0.;
    u3L = 0.;
    B1L = 1.;
    B2L = 0.;
    B3L = 0.;
    rhoR = 0.1;
    PR = 1.;
    u1R = 0.;
    u2R = 0.;
    u3R = 0.;
    B1R = 1.;
    B2R = 0.;
    B3R = 0.;
  } else if (shock == 5) { // Shocktube 2
    tf = 1.;

    rhoL = 1.;
    PL = 30.;
    u1L = 0.;
    u2L = 0.;
    u3L = 0.;
    B1L = 0.;
    B2L = 20.;
    B3L = 0.;
    rhoR = 0.1;
    PR = 1.;
    u1R = 0.;
    u2R = 0.;
    u3R = 0.;
    B1R = 0.;
    B2R = 0.;
    B3R = 0.;
  } else if (shock == 6) { // Collision
    tf = 1.22;

    rhoL = 1.;
    PL = 1.;
    u1L = 5.;
    u2L = 0.;
    u3L = 0.;
    B1L = 10.;
    B2L = 10.;
    B3L = 0.;
    rhoR = 1.;
    PR = 1.;
    u1R = -5.;
    u2R = 0.;
    u3R = 0.;
    B1R = 10.;
    B2R = -10.;
    B3R = 0.;
  } else {
    printf("Shock %i not supported! Exiting\n", shock);
    exit(-1);
  }

  ZLOOP {
    coord(i, j, k, CENT, X);

    P[i][j][k][RHO] = (X[1] < 0.) ? rhoL : rhoR;
    P[i][j][k][UU]  = ((X[1] < 0.) ? PL : PR)/(gam - 1.);
    P[i][j][k][U1]  = (X[1] < 0.) ? u1L : u1R;
    P[i][j][k][U2]  = (X[1] < 0.) ? u2L : u2R;
    P[i][j][k][U3]  = (X[1] < 0.) ? u3L : u3R;
    P[i][j][k][B1]  = (X[1] < 0.) ? B1L : B1R;
    P[i][j][k][B2]  = (X[1] < 0.) ? B2L : B2R;
    P[i][j][k][B3]  = (X[1] < 0.) ? B3L : B3R;
  } // ZLOOP
}

