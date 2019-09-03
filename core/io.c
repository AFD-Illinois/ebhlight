/*****************************************************************************
 *                                                                            *
 * IO.C                                                                       *
 *                                                                            *
 * HDF5 OUTPUT AND RESTART                                                    *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#include <hdf5.h>

#define MAX_GRID_DIM (5)

// Lifted from Ben Prather and George Wong's hdf5_utils.c
static char hdf5_cur_dir[STRLEN] = "/";
// Make a directory (in the current directory) with given name
// This doesn't take a full path, just a name
void hdf5_make_directory(const char *name, hid_t file_id)
{
  // Add current directory to group name
  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strlen(path));

  hid_t group_id = H5Gcreate2(file_id, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group_id);
}
void hdf5_set_directory(const char *path)
{
  strncpy(hdf5_cur_dir, path, STRLEN-1);
}
void hdf5_add_att(const void *att, const char *att_name, const char *data_name, hid_t file_id, hsize_t hdf5_type)
{                                                                                
  char path[STRLEN];                                                             
  strncpy(path, hdf5_cur_dir, STRLEN);                                           
  strncat(path, data_name, STRLEN - strlen(path));                               
                                                                                 
  hid_t attribute_id = H5Acreate_by_name(file_id, path, att_name, hdf5_type, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, hdf5_type, att);                                        
  H5Aclose(attribute_id);                                                        
}   
// Add an attribute named "units"
void hdf5_add_units(const char *name, const char *unit, hid_t file_id)
{
  hid_t string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, strlen(unit)+1);
  hdf5_add_att(unit, "units", name, file_id, string_type);
  H5Tclose(string_type);
}
void hdf5_write_str_list (const void *data, const char *name, hid_t file_id,
                     size_t str_len, size_t len)
{
  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strlen(path));

  // Adapted (stolen) from https://support.hdfgroup.org/ftp/HDF5/examples/C/
  hsize_t dims_of_char_dataspace[] = {len};

  hid_t vlstr_h5t = H5Tcopy(H5T_C_S1);
  H5Tset_size(vlstr_h5t, str_len);

  hid_t dataspace = H5Screate_simple(1, dims_of_char_dataspace, NULL);
  hid_t dataset = H5Dcreate(file_id, path, vlstr_h5t, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  herr_t err = H5Dwrite(dataset, vlstr_h5t, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if (err < 0) fprintf(stderr, "hdf5_write_str_list failed writing %s: %d", path, err); // If anything above fails, the write should too

  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Tclose(vlstr_h5t);
}  

void write_scalar(void *data, const char *name, hsize_t type);
void read_scalar(void *data, const char *name, hsize_t type);
void write_array(void *data, const char *name, hsize_t rank,
  hsize_t *fdims, hsize_t *fstart, hsize_t *fcount, hsize_t *mdims,
  hsize_t *mstart, hsize_t type);
void read_array(void *data, const char *name, hsize_t rank,
  hsize_t *fdims, hsize_t *fstart, hsize_t *fcount, hsize_t *mdims,
  hsize_t *mstart, hsize_t type);
hsize_t product_hsize_t(hsize_t a[], int size);

// Some macro tricks to reduce the number of lines of code
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
  write_array((void*)x, #x, 5, fdims_tens,				\
  fstart_tens, fcount_tens, mdims_tens_noghost, mstart_tens_noghost, type)
#define READ_ARRAY(x, rank, fdims, fstart, fcount, mdims, mstart, type) \
  read_array((void*)x, #x, rank, fdims, fstart, fcount, mdims, mstart, type)
#define READ_GRID(x, type) read_array((void*)x, #x, 3, fdims_grid, \
  fstart_grid, fcount_grid, mdims_grid, mstart_grid, type)
#define READ_PRIM(x, type) read_array((void*)x, #x, 4, fdims_prim, \
  fstart_prim, fcount_prim, mdims_prim, mstart_prim, type)
#define READ_VEC(x, type) read_array((void*)x, #x, 4, fdims_vec, \
  fstart_vec, fcount_vec, mdims_vec, mstart_vec, type)

hid_t file_id, plist_id;
hid_t filespace, memspace;
hid_t filespace_hdr, memspace_hdr;
hsize_t hdr_dims[1] = {1}, zero = 0, one = 1;
hsize_t mem_start[MAX_GRID_DIM], file_grid_start[MAX_GRID_DIM];
hsize_t file_grid_count[MAX_GRID_DIM], file_grid_dims[MAX_GRID_DIM];
hsize_t file_hdr_start[1] = {0}, file_hdr_count[1], file_hdr_dims[1] = {1};
hsize_t hdr_rank = 1, grid_rank;

hsize_t fdims_grid[3] = {N1TOT, N2TOT, N3TOT};
hsize_t fstart_grid[3];
hsize_t fcount_grid[3] = {N1, N2, N3};
hsize_t mdims_grid[3] = {N1+2*NG, N2+2*NG, N3+2*NG};
hsize_t mstart_grid[3] = {NG, NG, NG};

hsize_t mdims_grid_noghost[3] = {N1, N2, N3};
hsize_t mstart_grid_noghost[3] = {0, 0, 0};

hsize_t fdims_prim[4] = {N1TOT, N2TOT, N3TOT, NVAR};
hsize_t fstart_prim[4];
hsize_t fcount_prim[4] = {N1, N2, N3, NVAR};
hsize_t mdims_prim[4] = {N1+2*NG, N2+2*NG, N3+2*NG, NVAR};
hsize_t mstart_prim[4] = {NG, NG, NG, 0};

hsize_t fdims_vec[4] = {N1TOT, N2TOT, N3TOT, NDIM};
hsize_t fstart_vec[4];
hsize_t fcount_vec[4] = {N1, N2, N3, NDIM};
hsize_t mdims_vec[4] = {N1+2*NG, N2+2*NG, N3+2*NG, NDIM};
hsize_t mstart_vec[4] = {NG, NG, NG, 0};

hsize_t mdims_vec_noghost[4] = {N1, N2, N3, NDIM};
hsize_t mstart_vec_noghost[4] = {0, 0, 0, 0};

hsize_t fdims_tens[5] = {N1TOT, N2TOT, N3TOT, NDIM, NDIM};
hsize_t fstart_tens[5];
hsize_t fcount_tens[5] = {N1, N2, N3, NDIM, NDIM};
hsize_t mdims_tens[5] = {N1+2*NG, N2+2*NG, N3+2*NG, NDIM, NDIM};
hsize_t mstart_tens[5] = {NG, NG, NG, 0, 0};

hsize_t mdims_tens_noghost[5] = {N1, N2, N3, NDIM, NDIM};
hsize_t mstart_tens_noghost[5] = {0, 0, 0, 0, 0};

hsize_t fdims_hdr[1] = {1};
hsize_t fstart_hdr[1] = {0};
hsize_t fcount_hdr[1] = {1};
hsize_t mdims_hdr[1] = {1};
hsize_t mstart_hdr[1] = {0};

hid_t phfiletype, phmemtype, trackphfiletype, trackphmemtype;


#if ELECTRONS
const char vnams[NVAR][STRLEN] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3", "KEL", "KTOT"};
#else
const char vnams[NVAR][STRLEN] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3"};
#endif
static char version[STRLEN];
static int dump_id = 0, restart_id = 0, restart_perm_id = 0, fdump_id = 0;
static grid_double_type divb;

void init_io()
{
  strcpy(dumpdir, "dumps/");
  strcpy(restartdir, "restarts/");
  strcpy(xmfdir, "xmf/");
  strcpy(version, VERSION);
  int len = strlen(outputdir);
  memmove(dumpdir+len, dumpdir, strlen(dumpdir)+1);
  memmove(restartdir+len, restartdir, strlen(restartdir)+1);
  memmove(xmfdir+len, xmfdir, strlen(xmfdir)+1);
  for (int n = 0; n < len; ++n) {
    dumpdir[n] = outputdir[n];
    restartdir[n] = outputdir[n];
    xmfdir[n] = outputdir[n];
  }

  if (mpi_io_proc()) {
    char mkdircall[STRLEN];
    strcpy(mkdircall, "mkdir -p ");
    strcat(mkdircall, dumpdir);
    strcat(mkdircall, " ");
    strcat(mkdircall, restartdir);
    strcat(mkdircall, " ");
    strcat(mkdircall, xmfdir);
    safe_system(mkdircall);
  }

  fstart_grid[0] = global_start[1];
  fstart_grid[1] = global_start[2];
  fstart_grid[2] = global_start[3];
  fstart_prim[0] = global_start[1];
  fstart_prim[1] = global_start[2];
  fstart_prim[2] = global_start[3];
  fstart_prim[3] = 0;
  fstart_vec[0] = global_start[1];
  fstart_vec[1] = global_start[2];
  fstart_vec[2] = global_start[3];
  fstart_vec[3] = 0;
  fstart_tens[0] = global_start[1];
  fstart_tens[1] = global_start[2];
  fstart_tens[2] = global_start[3];
  fstart_tens[3] = 0;
  fstart_tens[4] = 0;
  if (!mpi_io_proc()) fcount_hdr[0] = 0;

  filespace_hdr = H5Screate_simple(1, file_hdr_dims, NULL);
  file_hdr_count[0] = (mpi_io_proc() ? one : zero);
  H5Sselect_hyperslab(filespace_hdr, H5S_SELECT_SET, file_hdr_start, NULL,
    file_hdr_count, NULL);
  memspace_hdr = H5Screate_simple(hdr_rank, file_hdr_count, NULL);

  #if RADIATION
    // Create custom datatypes for arrays inside struct
  hsize_t array_dim[2] = {NSUP, NDIM};
  hid_t array_tid = H5Tarray_create(H5T_NATIVE_DOUBLE, 2, array_dim);
  hsize_t int_vec_dim[1] = {NDIM};
  hid_t int_vec_tid = H5Tarray_create(H5T_NATIVE_INT, 1, int_vec_dim);

  // Memory contiguous in file
  phfiletype = H5Tcreate(H5T_COMPOUND, sizeof(struct of_photon));
  int offset = 0;
  // Function that accepts phfiletype, name, offset, size
  H5Tinsert(phfiletype, "X", offset, array_tid);
  offset += NSUP*NDIM*sizeof(double);
  H5Tinsert(phfiletype, "Kcov", offset, array_tid);
  offset += NSUP*NDIM*sizeof(double);
  H5Tinsert(phfiletype, "Kcon", offset, array_tid);
  offset += NSUP*NDIM*sizeof(double);
  H5Tinsert(phfiletype, "w", offset, H5T_NATIVE_DOUBLE);
  offset += sizeof(double);
  H5Tinsert(phfiletype, "KdotKprev", offset, H5T_NATIVE_DOUBLE);
  offset += sizeof(double);
  H5Tinsert(phfiletype, "nscatt", offset, H5T_NATIVE_INT);
  offset += sizeof(int);
  H5Tinsert(phfiletype, "origin", offset, int_vec_tid);
  offset += NDIM*sizeof(int);
  H5Tinsert(phfiletype, "t0", offset, H5T_NATIVE_DOUBLE);
  offset += sizeof(double);
  H5Tinsert(phfiletype, "is_tracked", offset, H5T_NATIVE_INT);

  // Use HOFFSET to account for struct padding in memory
  phmemtype = H5Tcreate(H5T_COMPOUND, sizeof(struct of_photon));
  H5Tinsert(phmemtype, "X", HOFFSET(struct of_photon, X), array_tid);
  H5Tinsert(phmemtype, "Kcov", HOFFSET(struct of_photon, Kcov),
    array_tid);
  H5Tinsert(phmemtype, "Kcon", HOFFSET(struct of_photon, Kcon),
    array_tid);
  H5Tinsert(phmemtype, "w", HOFFSET(struct of_photon, w),
    H5T_NATIVE_DOUBLE);
  H5Tinsert(phmemtype, "KdotKprev",
    HOFFSET(struct of_photon, KdotKprev), H5T_NATIVE_DOUBLE);
  H5Tinsert(phmemtype, "nscatt", HOFFSET(struct of_photon, nscatt),
    H5T_NATIVE_INT);
  H5Tinsert(phmemtype, "origin", HOFFSET(struct of_photon, origin),
    int_vec_tid);
  H5Tinsert(phmemtype, "t0", HOFFSET(struct of_photon, t0),
    H5T_NATIVE_DOUBLE);
  H5Tinsert(phmemtype, "is_tracked", HOFFSET(struct of_photon, is_tracked),
    H5T_NATIVE_INT);

  trackphfiletype = H5Tcreate(H5T_COMPOUND, sizeof(struct of_photon));
  offset = 0;
  H5Tinsert(trackphfiletype, "X1", offset, H5T_NATIVE_DOUBLE);
  offset += sizeof(double);
  H5Tinsert(trackphfiletype, "X2", offset, H5T_NATIVE_DOUBLE);
  offset += sizeof(double);
  H5Tinsert(trackphfiletype, "X3", offset, H5T_NATIVE_DOUBLE);
  offset += sizeof(double);
  H5Tinsert(trackphfiletype, "nscatt", offset, H5T_NATIVE_INT);

  trackphmemtype = H5Tcreate(H5T_COMPOUND, sizeof(struct of_track_photon));
  H5Tinsert(trackphmemtype, "X1", HOFFSET(struct of_track_photon, X1),
    H5T_NATIVE_DOUBLE);
  H5Tinsert(trackphmemtype, "X2", HOFFSET(struct of_track_photon, X2),
    H5T_NATIVE_DOUBLE);
  H5Tinsert(trackphmemtype, "X3", HOFFSET(struct of_track_photon, X3),
    H5T_NATIVE_DOUBLE);
  H5Tinsert(trackphmemtype, "nscatt", HOFFSET(struct of_track_photon, nscatt),
    H5T_NATIVE_INT);
  #endif // RADIATION
}

#if RADIATION
void track_ph()
{
  #if TRACK_PH
  static int track_id = 0;
  // Count number of tracked photons on this processor
  int num_tracked = 0;
  #pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      if (ph->is_tracked == 1) {
        #pragma omp atomic
        num_tracked++;
      }
      ph = ph->next;
    }
  }

  // Create buffer for output
  struct of_track_photon *io_track = safe_malloc(
    num_tracked*sizeof(struct of_track_photon));
  int nio = 0;
  double X[NDIM], Kcov[NDIM], Kcon[NDIM];
  for (int n = 0; n < nthreads; n++) {
    struct of_photon *ph = photon_lists[n];
    while (ph != NULL) {
      if (ph->is_tracked == 1) {
        get_X_K_interp(ph, t, X, Kcov, Kcon);

        io_track[nio].X1 = X[1];
        io_track[nio].X2 = X[2];
        io_track[nio].X3 = X[3];
        io_track[nio].nscatt = ph->nscatt;
        nio++;
      }
      ph = ph->next;
    }
  }

  if (Nph_to_track > 0.) {

    char name[STRLEN];
    char fname[STRLEN];
    sprintf(fname, "trackph_%08d.h5", track_id);
    strcpy(name, dumpdir);
    strcat(name, fname);

    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    // One dataset per processor; each processor must create all datasets
    hid_t *ph_dsets = safe_malloc(mpi_nprocs()*sizeof(hid_t));
    for (int n = 0; n < mpi_nprocs(); n++) {
      char dsetnam[STRLEN];
      sprintf(dsetnam, "trackph_%08d", n);
      int num_tracked_buf = num_tracked;
      mpi_int_broadcast_proc(&num_tracked_buf, n);
      hsize_t dims[1] = {num_tracked_buf};
      hid_t space = H5Screate_simple(1, dims, NULL);
      ph_dsets[n] = H5Dcreate(file_id, dsetnam, trackphfiletype, space,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(space);
    }
    if (num_tracked > 0)
      H5Dwrite(ph_dsets[mpi_myrank()], trackphmemtype, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, io_track);

    // Close and release resources
    free(io_track);
    for (int n = 0; n < mpi_nprocs(); n++) {
      H5Dclose(ph_dsets[n]);
    }
    free(ph_dsets);
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    track_id++;
  }
  #endif // TRACK_PH
}
#endif // RADIATION

void dump_grid()
{
  char name[STRLEN], fname[STRLEN];
  sprintf(fname, "grid.h5");
  strcpy(name, dumpdir);
  strcat(name, fname);
  if (mpi_io_proc()) {
    fprintf(stdout, "WRITING GEOMETRY TO %s\n", name);
  }

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  if (file_id < 0) {
    fprintf(stderr, "Could not create grid file! Exiting...\n");
    exit(-1);
  }
  H5Pclose(plist_id);

  {
    double *Xharm = safe_malloc(N1*N2*N3*NDIM*sizeof(double));
    int n = 0;
    ZLOOP
    {
      coord(i, j, k, CENT, &Xharm[n]);
      n += 4;
    }
    WRITE_VEC_NO_GHOSTS(Xharm, TYPE_DBL);
    free(Xharm);
  }

  { // If metric == Minkowski, then Xcart == Xharm
    double *Xcart = safe_malloc(N1*N2*N3*NDIM*sizeof(double));
    double X[NDIM],Xcart_loc[NDIM];
    int n = 0;
    ZLOOP {
      coord(i, j, k, CENT, X);
      cart_coord(X,Xcart_loc);
      for (int l = 0; l < NDIM; l++) Xcart[n+l] = Xcart_loc[l];
      n += 4;
    }
    WRITE_VEC_NO_GHOSTS(Xcart, TYPE_DBL);
    free(Xcart);
  }

  { // Face locations, in HARM and Cartesian coordinates
    #define RANK (4)
    hsize_t fdims[RANK] = {N1TOT+1, N2TOT+1, N3TOT+1, NDIM};
    hsize_t fstart[RANK] = {global_start[1],global_start[2],global_start[3],0};
    hsize_t fcount[RANK] = {N1, N2, N3, NDIM};
    hsize_t mdims[RANK] = {N1, N2, N3, NDIM};
    for (int d = 0; d < 3; d++) {
      if (global_stop[d+1] == fdims[d]-1) {
	fcount[d]++;
	mdims[d]++;
      }
    }
    hsize_t mstart[RANK] = {0, 0, 0, 0};
    hsize_t memsize_tot = product_hsize_t(mdims, RANK);
    double *XFharm = safe_malloc(memsize_tot*sizeof(double));
    double *XFcart = safe_malloc(memsize_tot*sizeof(double));
    double X[NDIM], Xcart_loc[NDIM];
    int n = 0;
    ZSLOOP(0, mdims[0]-1, 0, mdims[1]-1, 0, mdims[2]-1) {
      coord(i, j, k, CORN, X);
      cart_coord(X, Xcart_loc);
      for (int l = 0; l < NDIM; l++) {
	XFharm[n+l] = X[l];
	XFcart[n+l] = Xcart_loc[l];
      }
      n += 4;
    }
    WRITE_ARRAY(XFharm, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
    WRITE_ARRAY(XFcart, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
    free(XFharm);
    free(XFcart);
    #undef RANK
  }

  #if METRIC == MKS || METRIC == MMKS
  {
    double *Xbl = safe_malloc(N1*N2*N3*NDIM*sizeof(double));
    int n = 0;
    ZLOOP
    {
      double X[NDIM], r, th;
      coord(i, j, k, CENT, X);
      bl_coord(X, &r, &th);
      Xbl[n+1] = r;
      Xbl[n+2] = th;
      Xbl[n+3] = X[3];
      n += 4;
    }
    WRITE_VEC_NO_GHOSTS(Xbl, TYPE_DBL);
    free(Xbl);
  }
  #endif

  {
    double *gcov = safe_malloc(N1*N2*N3*NDIM*NDIM*sizeof(double));
    int n = 0;
    ZLOOP {
      DLOOP2 {
	gcov[n] = ggeom[i][j][CENT].gcov[mu][nu];
	n++;
      }
    }
    WRITE_TENSOR_NO_GHOSTS(gcov, TYPE_DBL);
    free(gcov);
  }

  {
    double *gcon = safe_malloc(N1*N2*N3*NDIM*NDIM*sizeof(double));
    int n = 0;
    ZLOOP {
      DLOOP2 {
	gcon[n] = ggeom[i][j][CENT].gcon[mu][nu];
	n++;
      }
    }
    WRITE_TENSOR_NO_GHOSTS(gcon, TYPE_DBL);
    free(gcon);
  }

  {
    double *gdet = safe_malloc(N1*N2*N3*sizeof(double));
    int n = 0;
    ZLOOP {
      gdet[n] = ggeom[i][j][CENT].g;
      n++;
    }
    WRITE_GRID_NO_GHOSTS(gdet, TYPE_DBL);
    free(gdet);
  }

  {
    double *alpha = safe_malloc(N1*N2*N3*sizeof(double));
    int n = 0;
    ZLOOP {
      alpha[n] = ggeom[i][j][CENT].alpha;
      n++;
    }
    WRITE_GRID_NO_GHOSTS(alpha, TYPE_DBL);
    free(alpha);
  }

  #if METRIC == MKS || METRIC == MMKS
  {
    double X[NDIM];
    double Lambda_con_local[NDIM][NDIM];
    double Lambda_cov_local[NDIM][NDIM];
    double *Lambda_h2bl_con = safe_malloc(N1*N2*N3*NDIM*NDIM*sizeof(double));
    double *Lambda_h2bl_cov = safe_malloc(N1*N2*N3*NDIM*NDIM*sizeof(double));
    int n = 0;
    ZLOOP {
      coord(i, j, k, CENT, X);
      jac_harm_to_bl(X, Lambda_cov_local, Lambda_con_local);
      DLOOP2 {
	Lambda_h2bl_con[n] = Lambda_con_local[mu][nu];
	Lambda_h2bl_cov[n] = Lambda_cov_local[mu][nu];
	n++;
      }
    }
    WRITE_TENSOR_NO_GHOSTS(Lambda_h2bl_con, TYPE_DBL);
    WRITE_TENSOR_NO_GHOSTS(Lambda_h2bl_cov, TYPE_DBL);
    free(Lambda_h2bl_con);
    free(Lambda_h2bl_cov);
  }

  {
    double X[NDIM];
    double Lambda_con_local[NDIM][NDIM];
    double Lambda_cov_local[NDIM][NDIM];
    double *Lambda_bl2cart_con = safe_malloc(N1*N2*N3*NDIM*NDIM*sizeof(double));
    double *Lambda_bl2cart_cov = safe_malloc(N1*N2*N3*NDIM*NDIM*sizeof(double));
    int n = 0;
    ZLOOP {
      coord(i, j, k, CENT, X);
      jac_bl_to_cart(X, Lambda_cov_local, Lambda_con_local);
      DLOOP2 {
	Lambda_bl2cart_con[n] = Lambda_con_local[mu][nu];
	Lambda_bl2cart_cov[n] = Lambda_cov_local[mu][nu];
	n++;
      }
    }
    WRITE_TENSOR_NO_GHOSTS(Lambda_bl2cart_con, TYPE_DBL);
    WRITE_TENSOR_NO_GHOSTS(Lambda_bl2cart_cov, TYPE_DBL);
    free(Lambda_bl2cart_con);
    free(Lambda_bl2cart_cov);
  }
  #endif // METRIC == MKS

  {
    double X[NDIM];
    double Lambda_con_local[NDIM][NDIM];
    double Lambda_cov_local[NDIM][NDIM];
    double *Lambda_h2cart_con = safe_malloc(N1*N2*N3*NDIM*NDIM*sizeof(double));
    double *Lambda_h2cart_cov = safe_malloc(N1*N2*N3*NDIM*NDIM*sizeof(double));
    int n = 0;
    ZLOOP {
      coord(i, j, k, CENT, X);
      jac_harm_to_cart(X, Lambda_cov_local, Lambda_con_local);
      DLOOP2 {
	Lambda_h2cart_con[n] = Lambda_con_local[mu][nu];
	Lambda_h2cart_cov[n] = Lambda_cov_local[mu][nu];
	n++;
      }
    }
    WRITE_TENSOR_NO_GHOSTS(Lambda_h2cart_con, TYPE_DBL);
    WRITE_TENSOR_NO_GHOSTS(Lambda_h2cart_cov, TYPE_DBL);
    free(Lambda_h2cart_con);
    free(Lambda_h2cart_cov);
  }

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);
}

void dump()
{
  timer_start(TIMER_OUT);

  char name[STRLEN], fname[STRLEN];
  sprintf(fname, "dump_%08d.h5", dump_id);
  strcpy(name, dumpdir);
  strcat(name, fname);
  if (mpi_io_proc()) {
    fprintf(stdout, "DUMP %s\n", name);
  }
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  if (file_id < 0) {
    fprintf(stderr, "Could not create dump file! Exiting...\n");
    exit(-1);
  }
  H5Pclose(plist_id);

  hdf5_set_directory("/");
  hdf5_make_directory("header", file_id);
  hdf5_set_directory("/header/");
 
  WRITE_HDR(version, TYPE_STR);

  #if RADIATION
  int has_radiation = RADIATION; WRITE_HDR(has_radiation, TYPE_INT);
  #endif
  #if ELECTRONS
  int has_electrons = ELECTRONS; WRITE_HDR(has_electrons, TYPE_INT);
  #endif

  WRITE_HDR(metric, TYPE_STR);
  char gridfile[STRLEN] = "grid.h5";
  WRITE_HDR(gridfile, TYPE_STR);
  WRITE_HDR(reconstruction, TYPE_STR);
  int n1 = N1TOT; WRITE_HDR(n1, TYPE_INT);
  int n2 = N2TOT; WRITE_HDR(n2, TYPE_INT);
  int n3 = N3TOT; WRITE_HDR(n3, TYPE_INT);
  int n_prim = NVAR; WRITE_HDR(n_prim, TYPE_INT);
  hdf5_write_str_list(vnams, "prim_names", file_id, STRLEN, NVAR);  

  WRITE_HDR(gam, TYPE_DBL);
  WRITE_HDR(cour, TYPE_DBL);
  WRITE_HDR(tf, TYPE_DBL);
  hdf5_add_units("tf", "code", file_id);
  #if ELECTRONS
  double gam_e = game; WRITE_HDR(gam_e, TYPE_DBL);
  double gam_p = gamp; WRITE_HDR(gam_p, TYPE_DBL);
  WRITE_HDR(tptemin, TYPE_DBL);
  WRITE_HDR(tptemax, TYPE_DBL);
  WRITE_HDR(fel0, TYPE_DBL);
  #endif
  #if RADIATION
  WRITE_HDR(tp_over_te, TYPE_DBL);
  WRITE_HDR(Mbh, TYPE_DBL);
  hdf5_add_units("Mbh", "g", file_id);
  hdf5_make_directory("units", file_id);
  int maxnscatt = MAXNSCATT; WRITE_HDR(maxnscatt, TYPE_INT);
  int nubins_emiss = NU_BINS_EMISS; WRITE_HDR(nubins_emiss, TYPE_INT);
  WRITE_HDR(numin_emiss, TYPE_DBL);
  hdf5_add_units("numin_emiss", "Hz", file_id);
  WRITE_HDR(numax_emiss, TYPE_DBL);
  hdf5_add_units("numax_emiss", "Hz", file_id);
  int nubins_spec = NU_BINS_SPEC; WRITE_HDR(nubins_spec, TYPE_INT);
  WRITE_HDR(numin_spec, TYPE_DBL);
  hdf5_add_units("numin_spec", "Hz", file_id);
  WRITE_HDR(numax_spec, TYPE_DBL);
  hdf5_add_units("numax_spec", "Hz", file_id);
  int nth = NTH; WRITE_HDR(nth, TYPE_INT);
  int nphi = NPHI; WRITE_HDR(nphi, TYPE_INT);
  WRITE_HDR(thetae_max, TYPE_DBL);
  WRITE_HDR(sigma_max, TYPE_DBL);
  WRITE_HDR(kdotk_tol, TYPE_DBL);

  hdf5_set_directory("/header/units/");
  WRITE_HDR(L_unit, TYPE_DBL);
  hdf5_add_units("L_unit", "cm", file_id);
  WRITE_HDR(M_unit, TYPE_DBL);
  hdf5_add_units("M_unit", "g", file_id);
  WRITE_HDR(T_unit, TYPE_DBL);
  hdf5_add_units("T_unit", "s", file_id);
  WRITE_HDR(RHO_unit, TYPE_DBL);
  hdf5_add_units("RHO_unit", "g cm^-3", file_id);
  WRITE_HDR(U_unit, TYPE_DBL);
  hdf5_add_units("U_unit", "erg", file_id);
  WRITE_HDR(B_unit, TYPE_DBL);
  hdf5_add_units("B_unit", "Gauss", file_id);
  WRITE_HDR(Ne_unit, TYPE_DBL);
  hdf5_add_units("Ne_unit", "cm^-3", file_id);
  WRITE_HDR(Thetae_unit, TYPE_DBL);
  hdf5_set_directory("/header/");
  #endif

  hdf5_make_directory("geom", file_id);
  hdf5_set_directory("/header/geom/");
  double startx1 = startx[1]; WRITE_HDR(startx1, TYPE_DBL);
  double startx2 = startx[2]; WRITE_HDR(startx2, TYPE_DBL);
  double startx3 = startx[3]; WRITE_HDR(startx3, TYPE_DBL);
  double dx1 = dx[1]; WRITE_HDR(dx1, TYPE_DBL);
  double dx2 = dx[2]; WRITE_HDR(dx2, TYPE_DBL);
  double dx3 = dx[3]; WRITE_HDR(dx3, TYPE_DBL);
  int n_dim = NDIM; WRITE_HDR(n_dim, TYPE_INT);
  #if METRIC == MKS
  hdf5_make_directory("mks", file_id);
  hdf5_set_directory("/header/geom/mks/");
  double r_eh = Reh; WRITE_HDR(r_eh, TYPE_DBL);
  hdf5_add_units("r_eh", "code", file_id);
  double r_in = Rin; WRITE_HDR(r_in, TYPE_DBL);
  hdf5_add_units("r_in", "code", file_id);
  double r_isco = Risco; WRITE_HDR(r_isco, TYPE_DBL);
  hdf5_add_units("r_isco", "code", file_id);
  double r_out = Rout; WRITE_HDR(r_out, TYPE_DBL);
  hdf5_add_units("r_out", "code", file_id);
  WRITE_HDR(a, TYPE_DBL);
  WRITE_HDR(hslope, TYPE_DBL);
    #if RADIATION
    double r_out_rad = Rout_rad; WRITE_HDR(r_out_rad, TYPE_DBL);
    hdf5_add_units("r_out_rad", "code", file_id);
    #endif
  #endif

  #if METRIC == MMKS
  hdf5_make_directory("mmks", file_id);
  hdf5_set_directory("/header/geom/mmks/");
  double r_eh = Reh; WRITE_HDR(r_eh, TYPE_DBL);
  hdf5_add_units("r_eh", "code", file_id);
  double r_in = Rin; WRITE_HDR(r_in, TYPE_DBL);
  hdf5_add_units("r_in", "code", file_id);
  double r_isco = Risco; WRITE_HDR(r_isco, TYPE_DBL);
  hdf5_add_units("r_isco", "code", file_id);
  double r_out = Rout; WRITE_HDR(r_out, TYPE_DBL);
  hdf5_add_units("r_out", "code", file_id);
  WRITE_HDR(a, TYPE_DBL);
  WRITE_HDR(hslope, TYPE_DBL);
  WRITE_HDR(poly_alpha, TYPE_DBL);
  WRITE_HDR(poly_xt, TYPE_DBL);
  WRITE_HDR(mks_smooth, TYPE_DBL);
    #if RADIATION
    double r_out_rad = Rout_rad; WRITE_HDR(r_out_rad, TYPE_DBL);
    hdf5_add_units("r_out_rad", "code", file_id);
    #endif
  #endif
 
  hdf5_set_directory("/");
  int is_full_dump = (dump_id % DTf == 0); WRITE_HDR(is_full_dump, TYPE_INT);
  WRITE_HDR(t, TYPE_DBL);
  hdf5_add_units("t", "code", file_id);
  WRITE_HDR(dt, TYPE_DBL);
  hdf5_add_units("dt", "code", file_id);
  int n_step = nstep; WRITE_HDR(n_step, TYPE_INT);
  int n_dump = dump_id; WRITE_HDR(n_dump, TYPE_INT);
  double dump_cadence = DTd; WRITE_HDR(dump_cadence, TYPE_DBL);
  hdf5_add_units("dump_cadence", "code", file_id);
  double full_dump_cadence = DTd*DTf; WRITE_HDR(full_dump_cadence, TYPE_DBL);
  hdf5_add_units("full_dump_cadence", "code", file_id);

  hdf5_set_directory("/");
  hdf5_make_directory("extras", file_id);
  hdf5_set_directory("/extras/");
  double r_out_vis = Rout_vis; WRITE_HDR(r_out_vis, TYPE_DBL);
  hdf5_add_units("r_out_vis", "code", file_id);
  
  hdf5_set_directory("/");
  WRITE_PRIM(P, "prims", TYPE_FLOAT);
  hdf5_add_units("prims", "code", file_id);

  current_calc();
  WRITE_VEC(jcon, TYPE_FLOAT);
  hdf5_add_units("jcon", "code", file_id);

  if (is_full_dump) {
    #pragma omp parallel for collapse(3)
    ZLOOP divb[i][j][k] = flux_ct_divb(i, j, k);
    WRITE_GRID(divb, "divb", TYPE_FLOAT);
    hdf5_add_units("divb", "code", file_id);
    WRITE_GRID(fail_save, "fail", TYPE_INT); 
   
    #if ELECTRONS
    WRITE_GRID(Qvisc_e, "Qvisc_e", TYPE_FLOAT);
    hdf5_add_units("Qvisc_e", "code", file_id);
    WRITE_GRID(Qvisc_p, "Qvisc_p", TYPE_FLOAT);
    hdf5_add_units("Qvisc_p", "code", file_id);
    #endif

    #if RADIATION
      #if ELECTRONS
      WRITE_GRID(Qcoul, "Qcoul", TYPE_FLOAT);
      hdf5_add_units("Qcoul", "code", file_id);
      #endif
      WRITE_GRID(Nsph, "Nsph", TYPE_INT);
      WRITE_GRID(nph, "nph", TYPE_FLOAT);
      WRITE_GRID(Nem, "Nem", TYPE_INT);
      WRITE_GRID(Nabs, "Nabs", TYPE_INT);
      WRITE_GRID(Nsc, "Nsc", TYPE_INT);
      WRITE_GRID(Esuper, "Esuper", TYPE_FLOAT);
      WRITE_GRID(Nsuper, "Nsuper", TYPE_INT);

      {
        #define RANK (4)
        hsize_t fdims[RANK] = {MAXNSCATT+2, N1TOT, N2TOT, N3TOT};
        hsize_t fstart[RANK] = {0, global_start[1], global_start[2],
          global_start[3]};
        hsize_t fcount[RANK] = {MAXNSCATT+2, N1, N2, N3};
        hsize_t mdims[RANK] = {MAXNSCATT+2, N1+2*NG, N2+2*NG, N3+2*NG};
        hsize_t mstart[RANK] = {0, NG, NG, NG};
        WRITE_ARRAY(Jrad, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_FLOAT);
        #undef RANK
      }
      hdf5_add_units("Jrad", "code", file_id);

      WRITE_TENSOR(Rmunu, TYPE_FLOAT);
      hdf5_add_units("Rmunu", "code", file_id);

      {
        mpi_reduce_nuLnu();
        #define RANK (4)
        hsize_t fdims[RANK] = {MAXNSCATT+1, NTH, NPHI, NU_BINS_SPEC};
        hsize_t fstart[RANK] = {0, 0, 0, 0};
        hsize_t fcount[RANK] = {MAXNSCATT+1, NTH, NPHI, NU_BINS_SPEC};
        hsize_t mdims[RANK] = {MAXNSCATT+1, NTH, NPHI, NU_BINS_SPEC};
        hsize_t mstart[RANK] = {0, 0, 0, 0};
        if (!mpi_io_proc()) {
          fcount[0] = 0;
          fcount[1] = 0;
          fcount[2] = 0;
          fcount[3] = 0;
        }
        WRITE_ARRAY(nuLnu, RANK, fdims, fstart, fcount, mdims, mstart,
          TYPE_DBL);
        #undef RANK
      }
      hdf5_add_units("nuLnu", "erg s^-1", file_id);
    #endif
  }

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  //if (mpi_io_proc()) {
  //  write_xml_file(dump_id, t, vnams);
  //}

  dump_id++;
  fdump_id++;
  reset_dump_variables();

  timer_stop(TIMER_OUT);
}

void restart_write(int restart_type)
{
  timer_start(TIMER_OUT);

  char name[STRLEN], fname[STRLEN];
  if (restart_type == RESTART_TEMP) {
    sprintf(fname, "restart_%08d.h5", restart_id);
  } else if (restart_type == RESTART_PERM) {
    sprintf(fname, "restart_perm_%08d.h5", restart_perm_id);
    restart_perm_id++;
  }
  strcpy(name, restartdir);
  strcat(name, fname);

  if (restart_type == RESTART_TEMP) {
    char lastname[STRLEN];
    sprintf(fname, "restart.last");
    strcpy(lastname, restartdir);
    strcat(lastname, fname);

    restart_id++;
    FILE *fp = fopen(lastname,"w");
    fprintf(fp,"%s\n", name);
    fclose(fp);
  }

  if(mpi_io_proc()) {
    fprintf(stdout, "RESTART %s\n", name);
  }

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  if (file_id < 0) {
    printf("Could not create restart file! Exiting...\n");
    exit(-1);
  }
  H5Pclose(plist_id);

  WRITE_HDR(version, TYPE_STR);
  WRITE_HDR(metric, TYPE_STR);
  int electrons = ELECTRONS; WRITE_HDR(electrons, TYPE_INT);
  int radiation = RADIATION; WRITE_HDR(radiation, TYPE_INT);
  int nvar = NVAR; WRITE_HDR(nvar, TYPE_INT);

  WRITE_HDR(t, TYPE_DBL);
  WRITE_HDR(tf, TYPE_DBL);
  WRITE_HDR(nstep, TYPE_INT);
  WRITE_HDR(dump_id, TYPE_DBL);
  WRITE_HDR(fdump_id, TYPE_DBL);
  WRITE_HDR(restart_id, TYPE_INT);
  WRITE_HDR(restart_perm_id, TYPE_INT);
  WRITE_HDR(dt, TYPE_DBL);
  WRITE_HDR(lim, TYPE_INT);
  WRITE_HDR(failed, TYPE_INT);

  WRITE_HDR(DTd, TYPE_DBL);
  WRITE_HDR(DTl, TYPE_DBL);
  WRITE_HDR(DTr, TYPE_DBL);
  WRITE_HDR(DNr, TYPE_INT);
  WRITE_HDR(DTp, TYPE_INT);
  WRITE_HDR(DTf, TYPE_INT);
  WRITE_HDR(tdump, TYPE_DBL);
  WRITE_HDR(trestart, TYPE_DBL);
  WRITE_HDR(tlog, TYPE_DBL);

  #if METRIC == CARTESIAN
    WRITE_HDR(x1Min, TYPE_DBL);
    WRITE_HDR(x1Max, TYPE_DBL);
    WRITE_HDR(x2Min, TYPE_DBL);
    WRITE_HDR(x2Max, TYPE_DBL);
    WRITE_HDR(x3Min, TYPE_DBL);
    WRITE_HDR(x3Max, TYPE_DBL);
  #elif METRIC == MKS || METRIC == MMKS
    WRITE_HDR(a, TYPE_DBL);
    WRITE_HDR(Rin, TYPE_DBL);
    WRITE_HDR(Rout, TYPE_DBL);
    WRITE_HDR(Rout_vis, TYPE_DBL);
    WRITE_HDR(hslope, TYPE_DBL);
    WRITE_HDR(Reh, TYPE_DBL);
    WRITE_HDR(Risco, TYPE_DBL);
    WRITE_HDR(poly_xt, TYPE_DBL);
    WRITE_HDR(poly_alpha, TYPE_DBL);
    WRITE_HDR(mks_smooth, TYPE_DBL);
    #if RADIATION
      WRITE_HDR(mbh, TYPE_DBL);
      WRITE_HDR(Mbh, TYPE_DBL);
    #endif // RADIATION
  #endif // METRIC

  #if RADIATION
    WRITE_HDR(tp_over_te, TYPE_DBL);
    WRITE_HDR(L_unit, TYPE_DBL);
    WRITE_HDR(T_unit, TYPE_DBL);
    WRITE_HDR(M_unit, TYPE_DBL);
    WRITE_HDR(RHO_unit, TYPE_DBL);
    WRITE_HDR(U_unit, TYPE_DBL);
    WRITE_HDR(B_unit, TYPE_DBL);
    WRITE_HDR(Ne_unit, TYPE_DBL);
    WRITE_HDR(Thetae_unit, TYPE_DBL);
  #endif

  WRITE_HDR(cour, TYPE_DBL);
  WRITE_HDR(gam, TYPE_DBL);

  #if ELECTRONS
    WRITE_HDR(game, TYPE_DBL);
    WRITE_HDR(gamp, TYPE_DBL);
    WRITE_HDR(fel0, TYPE_DBL);
    WRITE_HDR(tptemin, TYPE_DBL);
    WRITE_HDR(tptemax, TYPE_DBL);
  #endif // ELECTRONS

  #if RADIATION
    WRITE_HDR(numin_emiss, TYPE_DBL);
    WRITE_HDR(numax_emiss, TYPE_DBL);
    WRITE_HDR(numin_spec, TYPE_DBL);
    WRITE_HDR(numax_spec, TYPE_DBL);
    WRITE_HDR(tune_emiss, TYPE_DBL);
    WRITE_HDR(tune_scatt, TYPE_DBL);
    WRITE_HDR(t_tune_emiss, TYPE_DBL);
    WRITE_HDR(t_tune_scatt, TYPE_DBL);
    WRITE_HDR(dt_tune_emiss, TYPE_DBL);
    WRITE_HDR(dt_tune_scatt, TYPE_DBL);
    WRITE_HDR(made_tune_proc, TYPE_INT);
    WRITE_HDR(abs_tune_proc, TYPE_INT);
    WRITE_HDR(scatt_tune_proc, TYPE_INT);
    WRITE_HDR(thetae_max, TYPE_DBL);
    WRITE_HDR(sigma_max, TYPE_DBL);
    WRITE_HDR(kdotk_tol, TYPE_DBL);
    WRITE_HDR(Nph_to_track, TYPE_DBL);
    #if METRIC == MKS || METRIC == MMKS
      WRITE_HDR(Rout_rad, TYPE_DBL);
    #endif
  #endif

  WRITE_PRIM(P, "P", TYPE_DBL);

  #if RADIATION
    WRITE_GRID(Nsph, "Nsph", TYPE_INT);

    WRITE_GRID(nph, "nph", TYPE_DBL);

    WRITE_GRID(Nem, "Nem", TYPE_INT);

    WRITE_GRID(Nabs, "Nabs", TYPE_INT);

    WRITE_GRID(Nsc, "Nsc", TYPE_INT);
      
    WRITE_GRID(Esuper, "Esuper", TYPE_DBL);
    
    WRITE_GRID(Nsuper, "Nsuper", TYPE_INT);

    {
      #define RANK (4)
      hsize_t fdims[RANK] = {MAXNSCATT+2, N1TOT, N2TOT, N3TOT};
      hsize_t fstart[RANK] = {0, global_start[1], global_start[2],
        global_start[3]};
      hsize_t fcount[RANK] = {MAXNSCATT+2, N1, N2, N3};
      hsize_t mdims[RANK] = {MAXNSCATT+2, N1+2*NG, N2+2*NG, N3+2*NG};
      hsize_t mstart[RANK] = {0, NG, NG, NG};
      WRITE_ARRAY(Jrad, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
      #undef RANK
    }

    {
      mpi_reduce_nuLnu();
      #define RANK (4)
      hsize_t fdims[RANK] = {MAXNSCATT+1, NTH, NPHI, NU_BINS_SPEC};
      hsize_t fstart[RANK] = {0, 0, 0, 0};
      hsize_t fcount[RANK] = {MAXNSCATT+1, NTH, NPHI, NU_BINS_SPEC};
      hsize_t mdims[RANK] = {MAXNSCATT+1, NTH, NPHI, NU_BINS_SPEC};
      hsize_t mstart[RANK] = {0, 0, 0, 0};
      if (!mpi_io_proc()) {
        fcount[0] = 0;
        fcount[1] = 0;
        fcount[2] = 0;
        fcount[3] = 0;
      }
      WRITE_ARRAY(nuLnu, RANK, fdims, fstart, fcount, mdims, mstart,
        TYPE_FLOAT);
      #undef RANK
    }

    // Superphoton data
    hid_t space;
    char dsetnam[STRLEN];
    struct of_photon *wdata;
    wdata = safe_malloc(step_tot*sizeof(struct of_photon));

    // Copy superphotons into buffer
    int nph = 0;
    for (int n = 0; n < nthreads; n++) {
      struct of_photon *ph = photon_lists[n];
      while (ph != NULL) {
        copy_photon(ph, &(wdata[nph]));
        ph = ph->next;
        nph++;
      }
    }

    // Each processor must create all datasets
    hid_t *ph_dsets = safe_malloc(mpi_nprocs()*sizeof(hid_t));
    for (int n = 0; n < mpi_nprocs(); n++) {
      sprintf(dsetnam, "photons_%08d", n);
      //int step_tot_buf = step_tot;
      int step_tot_buf = nph;
      mpi_int_broadcast_proc(&step_tot_buf, n);
      hsize_t dims[1] = {step_tot_buf};
      space = H5Screate_simple (1, dims, NULL);
      ph_dsets[n] = H5Dcreate(file_id, dsetnam, phfiletype, space, H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose (space);
    }
    if (step_tot > 0) {
      H5Dwrite(ph_dsets[mpi_myrank()], phmemtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
        wdata);
    }

    // Close and release resources
    free(wdata);
    for (int n = 0; n < mpi_nprocs(); n++) {
      H5Dclose(ph_dsets[n]);
    }
    free(ph_dsets);
  #endif // RADIATION

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  // Keep only two most recent restart files
  if (mpi_io_proc() && restart_id >= 3 && restart_type == RESTART_TEMP) {
    char nametodel[STRLEN];
    sprintf(fname, "restart_%08d.h5", restart_id-3);
    strcpy(nametodel, restartdir);
    strcat(nametodel, fname);
    remove(nametodel);
  }

  timer_stop(TIMER_OUT);
}

void restart_read(char *fname)
{
  timer_start(TIMER_OUT);

  if(mpi_io_proc()) {
    fprintf(stderr, "Restarting from %s\n\n", fname);
  }
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);
  if (file_id < 0) {
    fprintf(stderr, "ERROR %s is not a valid restart file\n", fname);
    exit(-1);
  }
  H5Pclose(plist_id);

  READ_HDR(t, TYPE_DBL);
  READ_HDR(tf, TYPE_DBL);
  READ_HDR(nstep, TYPE_INT);
  READ_HDR(dump_id, TYPE_DBL);
  READ_HDR(fdump_id, TYPE_DBL);
  READ_HDR(restart_id, TYPE_INT);
  READ_HDR(restart_perm_id, TYPE_INT);
  READ_HDR(dt, TYPE_DBL);
  READ_HDR(lim, TYPE_INT);
  READ_HDR(failed, TYPE_INT);

  READ_HDR(DTd, TYPE_DBL);
  READ_HDR(DTl, TYPE_DBL);
  READ_HDR(DTr, TYPE_DBL);
  READ_HDR(DNr, TYPE_INT);
  READ_HDR(DTp, TYPE_INT);
  READ_HDR(DTf, TYPE_INT);
  READ_HDR(tdump, TYPE_DBL);
  READ_HDR(trestart, TYPE_DBL);
  READ_HDR(tlog, TYPE_DBL);

  #if METRIC == CARTESIAN
    READ_HDR(x1Min, TYPE_DBL);
    READ_HDR(x1Max, TYPE_DBL);
    READ_HDR(x2Min, TYPE_DBL);
    READ_HDR(x2Max, TYPE_DBL);
    READ_HDR(x3Min, TYPE_DBL);
    READ_HDR(x3Max, TYPE_DBL);
  #elif METRIC == MKS || METRIC == MMKS
    READ_HDR(a, TYPE_DBL);
    READ_HDR(Rin, TYPE_DBL);
    READ_HDR(Rout, TYPE_DBL);
    READ_HDR(Rout_vis, TYPE_DBL);
    READ_HDR(hslope, TYPE_DBL);
    READ_HDR(Reh, TYPE_DBL);
    READ_HDR(Risco, TYPE_DBL);
    READ_HDR(poly_xt, TYPE_DBL);
    READ_HDR(poly_alpha, TYPE_DBL);
    READ_HDR(mks_smooth, TYPE_DBL);
    #if RADIATION
      READ_HDR(mbh, TYPE_DBL);
      READ_HDR(Mbh, TYPE_DBL);
    #endif // RADIATION
  #endif // METRIC

  #if RADIATION
    READ_HDR(tp_over_te, TYPE_DBL);
    READ_HDR(L_unit, TYPE_DBL);
    READ_HDR(T_unit, TYPE_DBL);
    READ_HDR(M_unit, TYPE_DBL);
    READ_HDR(RHO_unit, TYPE_DBL);
    READ_HDR(U_unit, TYPE_DBL);
    READ_HDR(B_unit, TYPE_DBL);
    READ_HDR(Ne_unit, TYPE_DBL);
    READ_HDR(Thetae_unit, TYPE_DBL);
  #endif

  READ_HDR(cour, TYPE_DBL);
  READ_HDR(gam, TYPE_DBL);

  #if ELECTRONS
    READ_HDR(game, TYPE_DBL);
    READ_HDR(gamp, TYPE_DBL);
    READ_HDR(fel0, TYPE_DBL);
    READ_HDR(tptemin, TYPE_DBL);
    READ_HDR(tptemax, TYPE_DBL);
  #endif // ELECTRONS

  #if RADIATION
    READ_HDR(numin_emiss, TYPE_DBL);
    READ_HDR(numax_emiss, TYPE_DBL);
    READ_HDR(numin_spec, TYPE_DBL);
    READ_HDR(numax_spec, TYPE_DBL);
    READ_HDR(tune_emiss, TYPE_DBL);
    READ_HDR(tune_scatt, TYPE_DBL);
    READ_HDR(t_tune_emiss, TYPE_DBL);
    READ_HDR(t_tune_scatt, TYPE_DBL);
    READ_HDR(dt_tune_emiss, TYPE_DBL);
    READ_HDR(dt_tune_scatt, TYPE_DBL);
    READ_HDR(made_tune_proc, TYPE_INT);
    READ_HDR(abs_tune_proc, TYPE_INT);
    READ_HDR(scatt_tune_proc, TYPE_INT);
    READ_HDR(thetae_max, TYPE_DBL);
    READ_HDR(sigma_max, TYPE_DBL);
    READ_HDR(kdotk_tol, TYPE_DBL);
    READ_HDR(Nph_to_track, TYPE_DBL);
    #if METRIC == MKS || METRIC == MMKS
      READ_HDR(Rout_rad, TYPE_DBL);
    #endif
  #endif

  READ_PRIM(P, TYPE_DBL);

  #if RADIATION
    READ_GRID(Nsph, TYPE_INT);

    READ_GRID(nph, TYPE_DBL);

    READ_GRID(Nem, TYPE_INT);

    READ_GRID(Nabs, TYPE_INT);

    READ_GRID(Nsc, TYPE_INT);
      
    READ_GRID(Esuper, TYPE_DBL);
    
    READ_GRID(Nsuper, TYPE_INT);

    {
      #define RANK (4)
      hsize_t fdims[RANK] = {MAXNSCATT+2, N1TOT, N2TOT, N3TOT};
      hsize_t fstart[RANK] = {0, global_start[1], global_start[2],
        global_start[3]};
      hsize_t fcount[RANK] = {MAXNSCATT+2, N1, N2, N3};
      hsize_t mdims[RANK] = {MAXNSCATT+2, N1+2*NG, N2+2*NG, N3+2*NG};
      hsize_t mstart[RANK] = {0, NG, NG, NG};
      READ_ARRAY(Jrad, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
      #undef RANK
    }

    {
      #define RANK (4)
      hsize_t fdims[RANK] = {MAXNSCATT+1, NTH, NPHI, NU_BINS_SPEC};
      hsize_t fstart[RANK] = {0, 0, 0, 0};
      hsize_t fcount[RANK] = {MAXNSCATT+1, NTH, NPHI, NU_BINS_SPEC};
      hsize_t mdims[RANK] = {MAXNSCATT+1, NTH, NPHI, NU_BINS_SPEC};
      hsize_t mstart[RANK] = {0, 0, 0, 0};
      if (!mpi_io_proc()) {
        fcount[0] = 0;
        fcount[1] = 0;
        fcount[2] = 0;
        fcount[3] = 0;
      }
      READ_ARRAY(nuLnu, RANK, fdims, fstart, fcount, mdims, mstart,
        TYPE_DBL);
      #undef RANK
    }
  // Superphoton datasets -- each processor must create all datasets
  hsize_t dims[1];
  int nph_in = 0;
  hid_t *ph_dsets = safe_malloc(mpi_nprocs()*sizeof(hid_t));
  char dsetnam[STRLEN];
  for (int n = 0; n < mpi_nprocs(); n++) {
    sprintf(dsetnam, "photons_%08d", n);
    ph_dsets[n] = H5Dopen(file_id, dsetnam, H5P_DEFAULT);
    hid_t space = H5Dget_space(ph_dsets[n]);
    H5Sget_simple_extent_dims(space, dims, NULL);
    if (n == mpi_myrank()) nph_in = dims[0];
    H5Sclose(space);
  }
  
  struct of_photon *rdata = safe_malloc(nph_in*sizeof(struct of_photon));
  if (nph_in > 0) {
    H5Dread(ph_dsets[mpi_myrank()], phmemtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      rdata);

    struct of_photon *tmp = NULL, *ph_in = NULL;
    for (int n = 0; n < nph_in; n++) {
      tmp = malloc(sizeof(struct of_photon));
      copy_photon(&(rdata[n]), tmp);
      if (tmp->Kcov[2][0] > 0.) {
        printf("BAD PHOTON BEING READ IN!\n");
        printf("Kcov[0] = %e\n", tmp->Kcov[2][0]);
      }

      tmp->next = ph_in;
      ph_in = tmp;
    }

    int len_list = 0;
    struct of_photon *ph_in_tmp = ph_in;
    while (ph_in_tmp != NULL) {
      len_list++;
      ph_in_tmp = ph_in_tmp->next;
    }

    // Divide linked list among openmp processes equitably (fragmented?)
    int thread_start = (int)(get_rand()*nthreads);
    for (int n = thread_start; n < thread_start + nph_in; n++) {
      int index = n % nthreads;
      swap_ph(&ph_in, &(photon_lists[index]));
    }
  }

  free(rdata);
  for (int n = 0; n < mpi_nprocs(); n++) {
    H5Dclose(ph_dsets[n]);
  }
  free(ph_dsets);
  #endif // RADIATION

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  timer_stop(TIMER_OUT);
}

// Use fluid restart data as initial conditions, usually for a GRRMHD simulation
void init_fluid_restart()
{
  hsize_t file_grid_dims[4],file_start[4],file_count[4];
  hsize_t mem_grid_dims[4],mem_start[4];
  int fdims[3] = {N1TOT, N2TOT, N3TOT};
  int mdims[3] = {N1, N2, N3};

  // MAKE SURE N1TOT, N2TOT, N3TOT are right!

  if(mpi_io_proc()) {
    fprintf(stderr, "Initializing fluid from %s\n\n", init_from_grmhd);
  }

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fopen(init_from_grmhd, H5F_ACC_RDONLY, plist_id);
  H5Pclose(plist_id);

  if (file_id < 0) {
    fprintf(stderr, "ERROR GRMHD restart file %s not found\n", init_from_grmhd);
    exit(-1);
  }

  // Read some header data? gam? game? gamp?

  // Read fluid primitive data
  {
    char *name = "P";
    for (int d = 0; d < 3; d++)
      file_grid_dims[d] = fdims[d];
    file_grid_dims[3] = NVAR;  // For vectors
    hid_t filespace = H5Screate_simple(4, file_grid_dims, NULL);
    for (int d = 0; d < 3; d++) {
      file_start[d] = global_start[d+1];
      file_count[d] = mdims[d];
    }
    file_start[3] = 0;
    file_count[3] = NVAR;
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count,
      NULL);

    for (int d = 0; d < 3; d++) {
      mem_grid_dims[d] = mdims[d] + 2*NG;
    }
    mem_grid_dims[3] = NVAR;
    hid_t memspace = H5Screate_simple(4, mem_grid_dims, NULL);
    for (int d = 0; d < 3; d++)
      mem_start[d] = NG;
    mem_start[3] = 0;
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count,
      NULL);

    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    hid_t dset_id = H5Dopen(file_id, name, plist_id);
    H5Pclose(plist_id);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id,
      &P[0][0][0][0]);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
//    READ_PRIM(P, TYPE_FLOAT);
  }

  H5Fflush(file_id,H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  MPI_Barrier(MPI_COMM_WORLD);
}

int restart_init()
{
  char fname[STRLEN], lastname[STRLEN];
  sprintf(fname, "restart.last");
  strcpy(lastname, restartdir);
  strcat(lastname, fname);

  #if RADIATION
    photon_lists = malloc(nthreads*sizeof(struct of_photon*));
    #pragma omp parallel
    {
      photon_lists[omp_get_thread_num()] = NULL;
    }
  #endif // RADIATION

  FILE *fp = fopen(lastname, "r");
  if (fp == NULL) {
    if (mpi_io_proc()) {
      fprintf(stdout, "No restart file\n\n");
    }
    return 0;
  }

  if (mpi_io_proc())
    fprintf(stdout, "Loading restart file\n\n");
  zero_arrays();

  safe_fscanf(fp, "%s\n", fname);
  restart_read(fname);

  ZSLOOP(-NG, N1-1+NG, -NG, N2-1+NG, -NG, N3-1+NG) {
    PLOOP {
      Ph[i][j][k][ip] = P[i][j][k][ip];
    }
  }

  set_grid();

  bound_prim(P);

  #if RADIATION
		set_units();

		ZLOOP {
		 sim_vol += ggeom[i][j][CENT].g*dx[1]*dx[2]*dx[3]*L_unit*L_unit*L_unit;
		}
		sim_vol = mpi_reduce(sim_vol);

		photon_mpi_lists = malloc(nthreads*sizeof(struct of_photon*));
		#pragma omp parallel
		{
			photon_mpi_lists[omp_get_thread_num()] = NULL;
		}

		init_emissivity();

		set_weight(P);

		#if SCATTERING
			init_hotcross();
		#endif
  #endif // RADIATION

  return 1;
}

void write_scalar(void *data, const char *name, hsize_t type)
{
  hid_t memspace, plist_id, dset_id;
  char path[STRLEN];
  strcpy(path, hdf5_cur_dir);
  strcat(path, name);

  if (type == TYPE_STR) {
    type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, strlen(data));
  }

  memspace = H5Screate(H5S_SCALAR);
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, path, type, memspace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, type, memspace, memspace, plist_id, data);

  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(memspace);
}

void read_scalar(void *data, const char *name, hsize_t type)
{
  if (type == TYPE_STR) {
    printf("ERROR TYPE_STR is not supported by read_scalar()!\n");
    exit(-1);
  }
  hid_t memspace, plist_id, dset_id;
  char path[STRLEN];
  strcpy(path, hdf5_cur_dir);
  strcat(path, name);

  memspace = H5Screate(H5S_SCALAR);
  dset_id = H5Dopen(file_id, path, H5P_DEFAULT);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dread(dset_id, type, memspace, memspace, plist_id, data);
  
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(memspace);
}

void write_array(void *data, const char *name, hsize_t rank,
  hsize_t *fdims, hsize_t *fstart, hsize_t *fcount, hsize_t *mdims,
  hsize_t *mstart, hsize_t type)
{
  hid_t filespace, memspace, plist_id, dset_id;
  filespace = H5Screate_simple(rank, fdims, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, fstart, NULL, fcount, NULL);
  memspace = H5Screate_simple(rank, mdims, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mstart, NULL, fcount, NULL);
  plist_id = H5Pcreate(H5P_DATASET_CREATE);

  char path[STRLEN];
  strcpy(path, hdf5_cur_dir);
  strcat(path, name);

  if (type == TYPE_FLOAT) {
    int ntot = 1;
    for (int n = 0; n < rank; n++) {
      ntot *= mdims[n];
    }
    float *buf = safe_malloc(ntot*sizeof(float));
    for (int n = 0; n < ntot; n++) {
      buf[n] = (float)(((double*)data)[n]);
    }
    dset_id = H5Dcreate(file_id, path, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
      plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, buf);
    free(buf);
  } else if (type == TYPE_DBL) {
    dset_id = H5Dcreate(file_id, path, H5T_NATIVE_DOUBLE, filespace,
      H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  } else if (type == TYPE_INT) {
    dset_id = H5Dcreate(file_id, path, H5T_NATIVE_INT, filespace,
      H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
  } else if (type == TYPE_STR) {
    hid_t string_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(string_type, strlen(data));
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, path, string_type, filespace, H5P_DEFAULT,
      plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite(dset_id, string_type, memspace, filespace, plist_id, data);
  } else {
    fprintf(stderr, "type %llu not supported by write_array()! Exiting...\n",
      type);
    exit(-1);
  }
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(memspace);
  H5Sclose(filespace);
}

void read_array(void *data, const char *name, hsize_t rank,
  hsize_t *fdims, hsize_t *fstart, hsize_t *fcount, hsize_t *mdims,
  hsize_t *mstart, hsize_t type)
{
  hid_t filespace, memspace, plist_id, dset_id;
  filespace = H5Screate_simple(rank, fdims, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, fstart, NULL, fcount, NULL);
  memspace = H5Screate_simple(rank, mdims, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mstart, NULL, fcount, NULL);
  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dopen(file_id, name, plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);

  if (type == TYPE_DBL) {
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
    if (rank == 1 && fdims[0] == 1) {
      mpi_dbl_broadcast(data);
    }
  } else if (type == TYPE_INT) {
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
    if (rank == 1 && fdims[0] == 1) {
      mpi_int_broadcast(data);
    }
  } else {
    fprintf(stderr, "type %llu not supported by read_array()! Exiting...\n",
      type);
    exit(-1);
  }
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(memspace);
  H5Sclose(filespace);
}

hsize_t product_hsize_t(hsize_t a[], int size)
{
  hsize_t out = 1;
  for (int i = 0; i < size; i++) {
    out *= a[i];
  }
  return out;
}
