#ifndef IO_H
#define IO_H

#include <hdf5.h>

// global file pointer, defined in io.c
extern hid_t file_id;

// utilities for writing to an open file_id
void hdf5_make_directory(const char *name);
void hdf5_set_directory(const char *path);
void write_scalar(void *data, const char *name, hsize_t type);
void read_scalar(void *data, const char *name, hsize_t type);
void write_array(void *data, const char *name, hsize_t rank,
  hsize_t *fdims, hsize_t *fstart, hsize_t *fcount, hsize_t *mdims,
  hsize_t *mstart, hsize_t type);
void read_array(void *data, const char *name, hsize_t rank,
  hsize_t *fdims, hsize_t *fstart, hsize_t *fcount, hsize_t *mdims,
  hsize_t *mstart, hsize_t type);
hsize_t product_hsize_t(hsize_t a[], int size);

// some macro tricks to reduce the number of lines of code
#define TYPE_FLOAT H5T_NATIVE_FLOAT
#define TYPE_DBL H5T_NATIVE_DOUBLE
#define TYPE_INT H5T_NATIVE_INT
#define TYPE_STR H5T_NATIVE_CHAR

#define WRITE_HDR(x, type) write_scalar((void*)&x, #x, type)
#define READ_HDR(x, type) read_scalar((void*)&x, #x, type)

#define WRITE_ARRAY(x, rank, fdims, fstart, fcount, mdims, mstart, type) \
  write_array((void*)x, #x, rank, fdims, fstart, fcount, mdims, mstart, type)
#define WRITE_GRID(x, name, type) write_array((void*)x, name, 3, fdims_grid, \
  fstart_grid, fcount_grid, mdims_grid, mstart_grid, type)
#define WRITE_GRID_NO_GHOSTS(x, type) write_array((void*)x, #x, 3, fdims_grid, \
  fstart_grid, fcount_grid, mdims_grid_noghost, mstart_grid_noghost, type)
#define WRITE_PRIM(x, name, type) write_array((void*)x, name, 4, fdims_prim, \
  fstart_prim, fcount_prim, mdims_prim, mstart_prim, type)
#define WRITE_VEC(x, type) write_array((void*)x, #x, 4, fdims_vec, \
  fstart_vec, fcount_vec, mdims_vec, mstart_vec, type)
#define WRITE_VEC_NO_GHOSTS(x, type) \
  write_array((void*)x, #x, 4, fdims_vec, \
  fstart_vec, fcount_vec, mdims_vec_noghost, mstart_vec_noghost, type)
#define WRITE_TENSOR(x, type) write_array((void*)x, #x, 5, fdims_tens, \
  fstart_tens, fcount_tens, mdims_tens, mstart_tens, type)
#define WRITE_TENSOR_NO_GHOSTS(x, type) \
  write_array((void*)x, #x, 5, fdims_tens,        \
  fstart_tens, fcount_tens, mdims_tens_noghost, mstart_tens_noghost, type)
#define READ_ARRAY(x, rank, fdims, fstart, fcount, mdims, mstart, type) \
  read_array((void*)x, #x, rank, fdims, fstart, fcount, mdims, mstart, type)
#define READ_GRID(x, type) read_array((void*)x, #x, 3, fdims_grid, \
  fstart_grid, fcount_grid, mdims_grid, mstart_grid, type)
#define READ_PRIM(x, type) read_array((void*)x, #x, 4, fdims_prim, \
  fstart_prim, fcount_prim, mdims_prim, mstart_prim, type)
#define READ_VEC(x, type) read_array((void*)x, #x, 4, fdims_vec, \
  fstart_vec, fcount_vec, mdims_vec, mstart_vec, type)



#endif // IO_H
