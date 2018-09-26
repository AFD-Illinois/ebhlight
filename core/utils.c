/******************************************************************************
 *                                                                            *
 * UTILS.C                                                                    *
 *                                                                            *
 * GENERIC HELPER FUNCTIONS                                                   *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

void *safe_malloc(int size)
{
  // malloc(0) may or may not return NULL, depending on compiler.
  if (size == 0) return NULL;

  void *A = malloc(size);
  if (A == NULL) {
    fprintf(stderr, "Failed to malloc\n");
    exit(-1);
  }
  return A;
}

// Error-handling wrappers for standard C functions                              
void safe_system(const char *command)                                            
{                                                                                
  int systemReturn = system(command);                                            
  if (systemReturn == -1) {                                                      
    fprintf(stderr, "system() call %s failed! Exiting!\n", command);             
    exit(-1);                                                                    
  }                                                                              
}                                                                                
                                                                                 
void safe_fscanf(FILE *stream, const char *format, ...)                          
{                                                                                
  va_list args;                                                                  
  va_start(args, format);                                                        
  int vfscanfReturn = vfscanf(stream, format, args);                             
  va_end(args);                                                                  
  if (vfscanfReturn == -1) {                                                     
    fprintf(stderr, "fscanf() call failed! Exiting!\n");                         
    exit(-1);                                                                    
  }                                                                              
}  

