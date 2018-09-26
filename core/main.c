/******************************************************************************
 *                                                                            *
 * MAIN.C                                                                     *
 *                                                                            *
 * ESTABLISH MPI COMMUNICATION, LOOP OVER TIME, COMPLETE                      *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include "defs.h"
#include <time.h>

void init_core();
void init_final();

int main(int argc, char *argv[])
{
  #if !OPENMP
  omp_set_num_threads(1);
  #endif
  
  // Check for minimal required MPI thread support
  int threadSafety;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &threadSafety);
  if (threadSafety < MPI_THREAD_FUNNELED) {
    fprintf(stderr, "Thread support < MPI_THREAD_FUNNELED. Unsafe.\n");
    exit(1);
  }

  init_mpi();
  
  if (mpi_myrank() == 0) {
    fprintf(stdout, "\n          ************************************************************\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *                          BHLIGHT                         *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *          Ryan, Dolence & Gammie ApJ 807:31, 2015         *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *  B R Ryan                                                *\n");
    fprintf(stdout, "          *  J C Dolence                                             *\n");
    fprintf(stdout, "          *  C F Gammie                                              *\n");
    fprintf(stdout, "          *  S M Ressler                                             *\n");
    fprintf(stdout, "          *  J M Miller                                              *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *                         RESOLUTION                       *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *  N1TOT : %04d (N1CPU: %04d)                              *\n", N1TOT, N1CPU);
    fprintf(stdout, "          *  N2TOT : %04d (N2CPU: %04d)                              *\n", N2TOT, N2CPU);
    fprintf(stdout, "          *  N3TOT : %04d (N3CPU: %04d)                              *\n", N3TOT, N3CPU);
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *                          OPTIONS                         *\n");
    fprintf(stdout, "          *                                                          *\n");
    #if OPENMP
    fprintf(stdout, "          *  OPENMP                                                  *\n");
    #endif
    #if RADIATION
    fprintf(stdout, "          *  RADIATION                                               *\n");
    #endif
    #if ELECTRONS
    fprintf(stdout, "          *  ELECTRONS      (Ressler et al. MNRAS 454:1848, 2015)    *\n");
    #endif
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *                          SYNTAX                          *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *  -p /path/to/param_file                                  *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          ************************************************************\n\n");
  }


  // Read command line arguments
  strcpy(outputdir, ""); // Default value
  char pfname[STRLEN];
  for (int n = 0; n < argc; n++) {
    // Check for argv[n] of the form '-*'
    if (*argv[n] == '-' && *(argv[n]+1) != '\0' && *(argv[n]+2) == '\0' &&
        n < argc-1) {
      if (*(argv[n]+1) == 'p') { // Set parameter file path
        strcpy(pfname, argv[++n]);
      } 
    }
  }

  #pragma omp parallel
  {
    #pragma omp master
    {
      nthreads = omp_get_num_threads();
    } // omp master
  } // omp parallel

  // Initialize some global variables
  reset_log_variables();
  nstep = 0;

  // Perform initializations, either directly, via checkpoint, or initialize
  // with GRMHD data
  init_params(pfname);
  init_io();
  //init_random(300*mpi_myrank());
  init_random(1 + omp_get_wtime() + 300.*mpi_myrank());
  is_restart = restart_init();
  time_init();
  if (!is_restart) {
    init_core();
   
    // Set primitive variables
    printf("Init from GRMHD? %i\n", strcmp(init_from_grmhd, "No") != 0);
    if (strcmp(init_from_grmhd, "No") == 0) { // Initialize from problem.c
      init_prob();
      #if ELECTRONS
      init_electrons();
      #endif
    } else { // Initialize fluid from GRMHD restart file
      init_fluid_restart();
    }
    
    init_final();
    if (mpi_myrank() == 0) 
      fprintf(stdout, "Initial conditions generated\n\n");
  }

  // Initial diagnostics
  diag(DIAG_INIT);

  if (mpi_io_proc())
    fprintf(stdout, "t = %e tf = %e\n", t, tf);
  if (!is_restart)
    diag(DIAG_DUMP);

  if (mpi_io_proc())
    fprintf(stdout, "\nEntering main loop\n");
  int dumpThisStep = 0;
  
  // Loop over time
  while (t < tf) {
    dumpThisStep = 0;
    timer_start(TIMER_ALL);
    
    // Step variables forward in time
    step();
    nstep++;

    // Per-step output
    #if RADIATION
    step_made_all = mpi_reduce_int(step_made);
    step_abs_all = mpi_reduce_int(step_abs);
    step_scatt_all = mpi_reduce_int(step_scatt);
    step_lost_all = mpi_reduce_int(step_lost);
    step_rec_all = mpi_reduce_int(step_rec);
    step_sent_all = mpi_reduce_int(step_sent);
    step_rcvd_all = mpi_reduce_int(step_rcvd);
    step_fail_all = mpi_reduce_int(step_fail);
    step_tot_all = mpi_reduce_int(step_tot);
    double load_frac = 0, min_load = 0, max_load = 0;
    if (nstep % DTp == 0) {
      load_frac = ((double)mpi_nprocs()*step_tot)/(step_tot_all + SMALL);
      min_load = mpi_min(load_frac);
      max_load = mpi_max(load_frac);
    }
    #endif // RADIATION
    mpi_sync_output();
    if (mpi_io_proc()) {
      fprintf(stdout, "t = %10.5g dt = %10.5g n = %8d\n", t, dt, nstep);
      #if RADIATION
      fprintf(stdout, "[%i] made = %d abs = %d scatt = %d lost = %d rec = %d sent = %d rcvd = %d fail = %d tot = %d\n",
        mpi_myrank(), step_made, step_abs, step_scatt, step_lost, step_rec, 
        step_sent, step_rcvd, step_fail, step_tot);
      fprintf(stdout, "ALL made = %d abs = %d scatt = %d lost = %d rec = %d sent = %d rcvd = %d fail = %d tot = %d\n", 
        step_made_all, step_abs_all, step_scatt_all, step_lost_all, step_rec_all,
        step_sent_all, step_rcvd_all, step_fail_all, step_tot_all);
      
      if (step_sent_all != step_rcvd_all) {
        printf("PHOTON MPI ANOMALY!\n");
        printf("sent = %i rcvd = %i\n", step_sent_all, step_rcvd_all);
      }
      if (step_tot > 1.e8) {
        printf("Too many superphotons on this node!\n");
        exit(-1);
      }
      #endif
    }

    #if RADIATION
    update_superphoton_resolution(P);
    #endif

    // File I/O with set frequencies
    if (t < tf) {
      if (t >= tdump) {
        dumpThisStep = 1;
        diag(DIAG_DUMP);
        tdump += DTd;
      }
      if (t >= tlog) {
        diag(DIAG_LOG);
        tlog += DTl;
      }
      if (nstep % DNr == 0) {
        restart_write(RESTART_TEMP);
      }
      if (t >= trestart) {
        restart_write(RESTART_PERM);
        trestart += DTr;
      }
    }
  
    reset_log_variables();

    timer_stop(TIMER_ALL);

    if (nstep == 1) timers_reset();

    if (nstep % DTp == 0) {
      #if RADIATION
      // Load balancing report
      if (mpi_io_proc())
      fprintf(stdout, "  Superphoton load MIN: %e MAX: %e\n", min_load, max_load);
      #endif
      report_performance();
    }
  }
  
  if (dumpThisStep == 0) diag(DIAG_FINAL);

  MPI_Finalize();

  return 0;
}

void init_core()
{
  zero_arrays();
  reset_dump_variables();
  set_grid();
    
  failed = 0;
  t = 0.;
  dump_cnt = 0;
  rdump_cnt = 0;

  #if RADIATION
  set_units();
  dt = cour*dt_light_min;

  photon_lists = malloc(nthreads*sizeof(struct of_photon*));
  photon_mpi_lists = malloc(nthreads*sizeof(struct of_photon*));
  #pragma omp parallel
  {
    photon_lists[omp_get_thread_num()] = NULL;
    photon_mpi_lists[omp_get_thread_num()] = NULL;
  }
  #endif
}

void init_final()
{
  #if RADIATION
  init_rad(P);
  set_weight(P); 
  init_superphoton_resolution();
  #endif

  fixup(P);
  bound_prim(P);

  tdump = DTd;
  trestart = DTr;
  tlog  = DTl;
  dtsave = dt;

  dump_grid();
}

