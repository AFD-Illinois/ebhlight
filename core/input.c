/******************************************************************************
 *                                                                            *
 * INPUT.C                                                                    *
 *                                                                            *
 * READ IN PARAMETER FILES                                                    *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include <ctype.h>

struct param {
  char *key;
  void *data;
};

#define MAXPARAMS (1000)
static struct param table[MAXPARAMS];
static int nparam = 0;
static int nparamset = 0;

void set_param(char *key, void *data)
{
  table[nparam].key = key;
  table[nparam].data = data;
  nparam++;
}

int get_param(char *key, void **data)
{
  int n = 0;
  while (strncmp(key, table[n].key, strlen(key)) != 0) {
    n++;
    if (n >= nparam) {
      *data = NULL;
      return 0;
    }
  }
  *data = table[n].data;

  return 1;
}

char metric[STRLEN], reconstuction[STRLEN];
void set_core_params()
{
  #if RECONSTRUCTION == LINEAR
  sprintf(reconstruction, "linear");
  #elif RECONSTRUCTION == WENO
  sprintf(reconstruction, "weno");
  #endif

  set_param("tf", &tf);
  set_param("dt", &dt);
  #if METRIC == CARTESIAN
  sprintf(metric, "CARTESIAN");
  set_param("x1Min", &x1Min);
  set_param("x1Max", &x1Max);
  set_param("x2Min", &x2Min);
  set_param("x2Max", &x2Max);
  set_param("x3Min", &x3Min);
  set_param("x3Max", &x3Max);
  #elif METRIC == MKS
  sprintf(metric, "MKS");
  set_param("a", &a);
  set_param("hslope", &hslope);
  if (N2TOT < NG) hslope = 1.; // Surprising errors if you don't do this
  set_param("Rout", &Rout);
  set_param("Rout_vis", &Rout_vis);
  #if RADIATION
  set_param("Rout_rad", &Rout_rad);
  #endif
  #elif METRIC == MMKS
  sprintf(metric, "MMKS");
  set_param("a", &a);
  set_param("hslope", &hslope);
  set_param("poly_xt", &poly_xt);
  set_param("poly_alpha", &poly_alpha);
  set_param("mks_smooth", &mks_smooth);
  if (N2TOT < NG) hslope = 1.;
  set_param("Rout", &Rout);
  set_param("Rout_vis", &Rout_vis);
  #if RADIATION
  set_param("Rout_rad", &Rout_rad);
  #endif
  #endif

  set_param("cour", &cour);
  set_param("gam", &gam);

  #if RADIATION
  #if METRIC == CARTESIAN
  set_param("L_unit", &L_unit);
  set_param("M_unit", &M_unit);
  #elif METRIC == MKS || METRIC == MMKS
  set_param("mbh", &mbh);
  set_param("M_unit", &M_unit);
  #endif
  #endif

  #if ELECTRONS
  set_param("game", &game);
  set_param("gamp", &gamp);
  set_param("fel0", &fel0);
  set_param("tptemin", &tptemin);
  set_param("tptemax", &tptemax);
  #endif

  #if RADIATION
  #if !ELECTRONS
  set_param("tp_over_te", &tp_over_te);
  #endif
  set_param("nph_per_proc", &nph_per_proc);
  set_param("numin_emiss", &numin_emiss);
  set_param("numax_emiss", &numax_emiss);
  set_param("numin_spec", &numin_spec);
  set_param("numax_spec", &numax_spec);
  set_param("tune_emiss", &tune_emiss);
  set_param("tune_scatt", &tune_scatt);
  set_param("t0_tune_emiss", &t_tune_emiss);
  set_param("t0_tune_scatt", &t_tune_scatt);
  set_param("thetae_max", &thetae_max);
  set_param("sigma_max", &sigma_max);
  set_param("kdotk_tol", &kdotk_tol);
  set_param("Nph_to_track", &Nph_to_track);
  #if GRAYABSORPTION
  set_param("kappa", &kappa);
  #endif
  //set_param("thbin", &thbin);
  //set_param("phibin", &phibin);
  #endif

  set_param("init_from_grmhd", &init_from_grmhd);

  set_param("DTd", &DTd);
  set_param("DTl", &DTl);
  set_param("DTr", &DTr);
  set_param("DNr", &DNr);
  set_param("DTp", &DTp);
  set_param("DTf", &DTf);
  set_param("outputdir", &outputdir);
}

void read_params(char *pfname)
{
  void *ptr;

  FILE *fp = fopen(pfname, "r");                                                 
  if (fp == NULL) {                                                              
    fprintf(stderr, "Cannot open parameter file: %s\n", pfname);                 
    exit(-1);                                                                    
  }                                                                              
                                                                                 
  char line[STRLEN];                                                             
  while (fgets(line, STRLEN, fp)) {                                              
    // Ignore comments, newlines, and leading whitespace                         
    if (line[0] == '#' || line[0] == '\n' || isspace(line[0]))                   
      continue;                                                                  
                                                                                 
    // Is key in dictionary, and is variable empty?                              
    char test[STRLEN], key[STRLEN];                                              
    test[0] = '\0';                                                              
    sscanf(line, "%*s %s %*s %s", key, test);                                    
    char *word = test;                                                           
    while(isspace(*word)) {                                                      
      word++;                                                                    
    }                                                                            
    if (word[0] == '\0') {                                                       
      continue;                                                                  
    }                                                                            
                                                                                 
    // Read in parameter depending on datatype                                   
    char type[6];                                                                
    strncpy(type, line, 5);                                                      
    type[5] = 0;                                                                 
    if (get_param(key, &ptr)) {                                                  
      if (!strncmp(type, "[int]", 5)) {                                          
        int buf;                                                                 
        sscanf(line, "%*s %s %*s %d", key, &buf);                                
        *((int*)ptr) = buf;                                                      
        nparamset++;                                                             
      } else if (!strncmp(type, "[dbl]", 5)) {                                   
        double buf;                                                              
        sscanf(line, "%*s %s %*s %lf", key, &buf);                               
        *((double*)ptr) = buf;                                                   
        nparamset++;                                                             
      } else if (!strncmp(type, "[str]", 5)) {                                   
        char buf[STRLEN];                                                        
        sscanf(line, "%*s %s %*s %s", key, buf);                                 
        strcpy((char*)ptr, buf);                                                 
        nparamset++;                                                             
      }                                                                          
    }                                                                            
  }
                                                                                 
  #if (METRIC == MKS || METRIC == MMKS) && RADIATION                                               
  Mbh = mbh*MSUN;                                                                
  #endif                                                                         
                                                                                 
  if (nparamset != nparam && mpi_io_proc()) {                                    
    fprintf(stderr, "Set %i parameters, needed %i!\n", nparamset, nparam);       
    exit(-1);                                                                    
  }                                                                              
                                                                                 
  fclose(fp);                                                                    
                                                                                 
  if (mpi_io_proc()) fprintf(stdout, "Parameter file read\n\n"); 
}

void init_params(char *pfname)
{
  set_core_params();
  set_problem_params();
  read_params(pfname);
}

