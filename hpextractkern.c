#include<stdio.h>
#include<string.h>
#include<strings.h>
#include<math.h>
#include<stdlib.h>
#include<fitsio.h>

#define max(x,y) x>y?x:y
#define min(x,y) x<y?x:y


typedef struct
{
  int xsize; /* image size in x */
  int ysize; /* image size in y */

  int fw; /* fwKernel */
  int ngauss;
  int ncompker; /* nCompKer */
  int ncomp; /* nComp */
  int nbgvec; /* nBGVectors */
  int ncomptot; /* nCompTotal */
  int order; /* kerOrder */
  int bgorder; /* bgOrder */

  int *deg_fixe;
  float *sigma_gauss;
  double *filter_x;
  double *filter_y;
  double **kernel_vec;
  double *kernel_coeffs;
  double *kernel;

  int region; /* currently selected region */
  double *solution; /* kernelSol ; points to the solution for current region*/
  double *ksumim; /* kSumIm */

  /* region by region info */
  int nreg; /* total number of regions */
  int *rxmins;
  int *rxmaxs;
  int *rymins;
  int *rymaxs;
  double *ksumims;
  double **solutions;

} Kernel;


/* forward declarations */
void get_kernel_vec(Kernel *);
double make_kernel(Kernel *, int, int, int, int);
/*void spatial_convolve(float *, int, int, float *, Kernel *, float, float);*/
void release_kernel(Kernel *ck);


void fitserr(int status, Kernel *ck)
{
  /**
   * Print out cfitsio error messages and exit program 
   */
  if (status) {
    fits_report_error(stderr, status);
    if (ck != NULL)
      release_kernel(ck);
    exit(status);
  }
  return;
}


void prgexit(int status, Kernel *ck)
{
  if (status > 100) {
    /* error code > 100 is from CFITSIO */
    fits_report_error(stderr, status);
  }
  if (ck != NULL)
    release_kernel(ck);
  exit(status);
  return;
}


double *kernel_vector(Kernel *ck, int n, int deg_x, int deg_y, int ig, int *ren)
{
  /**
   * Creates kernel sized entry for kernel_vec for each kernel degree 
   *   Mask of filter_x * filter_y, filter = exp(-x**2 sig) * x^deg 
   *   Subtract off kernel_vec[0] if n > 0
   * NOTE: this does not use any image
   */
  double *vector = NULL, *kernel0 = NULL;
  int i, j, k, dx, dy, ix;
  double sum_x, sum_y, x, qe;
  
  vector = ck->kernel_vec[n];
  dx = (deg_x / 2) * 2 - deg_x;
  dy = (deg_y / 2) * 2 - deg_y;
  sum_x = sum_y = 0.0;
  *ren = 0;
  
  for (ix = 0; ix < ck->fw; ix++) {
    x = (double)(ix - ck->fw / 2);
    k = ix + n * ck->fw;
    qe = exp(-x * x * ck->sigma_gauss[ig]);
    ck->filter_x[k] = qe * pow(x, deg_x);
    ck->filter_y[k] = qe * pow(x, deg_y);
    sum_x += ck->filter_x[k];
    sum_y += ck->filter_y[k];
  }
  
  if (n > 0)
    kernel0 = ck->kernel_vec[0];
  
  sum_x = 1. / sum_x;
  sum_y = 1. / sum_y;
  
  if (dx == 0 && dy == 0) {
    for (ix = 0; ix < ck->fw; ix++) {
      ck->filter_x[ix + n * ck->fw] *= sum_x;
      ck->filter_y[ix + n * ck->fw] *= sum_y;
    }
    for (i = 0; i < ck->fw; i++) {
      for (j = 0; j < ck->fw; j++) {
        vector[i + ck->fw * j] = ck->filter_x[i + n * ck->fw] * ck->filter_y[j + n * ck->fw];
      }
    }
    if (n > 0) {
      for (i = 0; i < ck->fw * ck->fw; i++) {
        vector[i] -= kernel0[i];
      }
      *ren = 1;
    }
  }
  else {
    for (i = 0; i < ck->fw; i++)
      for (j = 0; j < ck->fw; j++)
        vector[i + ck->fw * j] = ck->filter_x[i + n * ck->fw] * ck->filter_y[j + n * ck->fw];
  }
  return vector;
}


int init_kernel(char *kimage, Kernel *ck, int *ostatus)
{
  /**
   * Get all 1-time info from kernel fits header, overriding defaults
   * and command line options.
   */
  fitsfile *pf;
  char hkw[1024];
  int i, dummy, status = 0, info = 0, rxmin, rxmax, rymin, rymax;
  char **regions = NULL;

  /* set default */
  ck->ngauss = 3;
  ck->filter_x = NULL;
  ck->filter_y = NULL;
  ck->kernel_vec = NULL;
  ck->kernel_coeffs = NULL;
  ck->kernel = NULL;
  ck->deg_fixe = NULL;
  ck->sigma_gauss = NULL;
  ck->solution = NULL;
  ck->rxmins = NULL;
  ck->rxmaxs = NULL;
  ck->rymins = NULL;
  ck->rymaxs = NULL;
  ck->ksumims = NULL;
  ck->solutions = NULL;
  ck->nreg = 1;

  if (! fits_open_file(&pf, kimage, 0, &status)) {
    /* required keyword in primary HDU */
    fits_read_key_log(pf, "KERINFO", &info, NULL, &status);
    if (! info) {
      fprintf(stderr, "This image does not appear to contain a kernel table\n");
      status = 1;
      goto error_exit;
    }

    if (fits_movabs_hdu(pf, 1, NULL, &status))
      goto error_exit;

    /* grab all regions in the primary image header, up to 10 in number */
    regions = (char**)malloc(10 * sizeof(char*));
    for (i = 0; i < 10; i++)
      regions[i] = (char*)malloc(80 * sizeof(char));
    if (fits_read_keys_str(pf, "REGION", 0, 9, regions, &(ck->nreg), &status))
      goto error_exit;
    
    /* move to binary kernel table... */
    if (fits_movnam_hdu(pf, BINARY_TBL, "CONVOLUTION KERNEL INFORMATION", 0, &status))
      goto error_exit;

    if (fits_read_key(pf, TINT, "NGAUSS", &(ck->ngauss), NULL, &status) ||
        fits_read_key(pf, TINT, "FWKERN", &(ck->fw), NULL, &status) ||
        fits_read_key(pf, TINT, "CKORDER", &(ck->order), NULL, &status) ||
        fits_read_key(pf, TINT, "BGORDER", &(ck->bgorder), NULL, &status))
      goto error_exit;
  
    ck->deg_fixe = (int *)calloc(ck->ngauss, sizeof(int));
    ck->sigma_gauss = (float *)calloc(ck->ngauss, sizeof(float));
  
    /* read kernel gaussian info */
    for (i = 0; i < ck->ngauss; i++) {
      sprintf(hkw, "DGAUSS%d", i + 1);
      if (fits_read_key(pf, TINT, hkw, &(ck->deg_fixe[i]), NULL, &status))
        goto error_exit;
      sprintf(hkw, "SGAUSS%d", i + 1);
      if (fits_read_key(pf, TFLOAT, hkw, &(ck->sigma_gauss[i]), NULL, &status))
        goto error_exit;
      
      /* important! */
      ck->sigma_gauss[i] = (1.0 / (2.0 * ck->sigma_gauss[i] * ck->sigma_gauss[i]));
    }
  
    /* set array and comp sizes */
    ck->ncompker = 0;
    for (i = 0; i < ck->ngauss; i++)
      ck->ncompker += ((ck->deg_fixe[i] + 1) * (ck->deg_fixe[i] + 2)) / 2;
    ck->ncomp = ((ck->order + 1) * (ck->order + 2)) / 2;
    ck->nbgvec = ((ck->bgorder + 1) * (ck->bgorder + 2)) / 2;
    ck->ncomptot = ck->ncompker * ck->ncomp + ck->nbgvec;

    ck->solution = (double *)realloc(ck->solution, (ck->ncomptot + 1) * sizeof(double));
    ck->filter_x = (double *)malloc(ck->ncompker * ck->fw * sizeof(double));
    ck->filter_y = (double *)malloc(ck->ncompker * ck->fw * sizeof(double));
    ck->kernel = (double *)malloc(ck->fw * ck->fw * sizeof(double));
    ck->kernel_vec = (double **)malloc(ck->ncompker * sizeof(double *));
    ck->kernel_coeffs = (double *)malloc(ck->ncompker * sizeof(double));

    if (ck->solution == NULL || ck->filter_x == NULL || ck->filter_y == NULL
        || ck->kernel == NULL || ck->kernel_vec == NULL
        || ck->kernel_coeffs == NULL) {
      fprintf(stderr, "Out of memory\n");
      status = 1;
      goto error_exit;
    }

    /*get kernel vector; called only once*/
    int ig, idegx, idegy, nvec, ren;
    nvec = 0;
    for (ig = 0; ig < ck->ngauss; ig++) {
      for (idegx = 0; idegx <= ck->deg_fixe[ig]; idegx++) {
        for (idegy = 0; idegy <= ck->deg_fixe[ig] - idegx; idegy++) {
          /* stores kernel weight mask for each order */
          ck->kernel_vec[nvec] = (double *)malloc(ck->fw * ck->fw * sizeof(double));
          kernel_vector(ck, nvec, idegx, idegy, ig, &ren);
          nvec++;
        }
      }
    }

    ck->rxmins = (int*)malloc(ck->nreg * sizeof(int));
    ck->rxmaxs = (int*)malloc(ck->nreg * sizeof(int));
    ck->rymins = (int*)malloc(ck->nreg * sizeof(int));
    ck->rymaxs = (int*)malloc(ck->nreg * sizeof(int));
    ck->ksumims = (double*)malloc(ck->nreg * sizeof(double));
    ck->solutions = (double**)malloc(ck->nreg * sizeof(double*));

    for (i = 0; i < ck->nreg; ++i) {
      if (sscanf(regions[i], "[%d:%d,%d:%d]", &rxmin, &rxmax, &rymin, &rymax) != 4) {
        fprintf(stderr, "Problem with region %d (%s), exiting...\n", i, regions[i]);
        goto error_exit;
      }

      ck->rxmins[i] = rxmin;
      ck->rxmaxs[i] = rxmax;
      ck->rymins[i] = rymin;
      ck->rymaxs[i] = rymax;

      if (fits_movabs_hdu(pf, 1, NULL, &status))
        goto error_exit;

      sprintf(hkw, "KSUM%02d", i);
      if (fits_read_key(pf, TDOUBLE, hkw, &(ck->ksumims[i]), NULL, &status))
        goto error_exit;

      ck->solutions[i] = (double*)malloc((ck->ncomptot + 1) * sizeof(double));

      if (fits_movnam_hdu(pf, BINARY_TBL, "CONVOLUTION KERNEL INFORMATION", 0, &status))
        goto error_exit;

      memset(ck->solutions[i], 0, (ck->ncomptot + 1) * sizeof(double));

      if (fits_read_col(pf, TDOUBLE, i + 1, 1, 1, ck->ncomptot + 1, 0, ck->solutions[i], 0, &status))
        goto error_exit;
    }

    for (i = 0; i < 10; i++)
      if (regions[i]) free(regions[i]);
    if (regions) free(regions);
  
  }
  else {
    goto error_exit;
  }

  return 0;

 error_exit:
  fits_close_file(pf, &dummy);
  *ostatus = status;
  return *ostatus;
}


void release_kernel(Kernel *ck)
{
  /**
   * Free memory allocated in init_kernel() 
   */
  int i, ig, idegx, idegy, nvec = 0;
  if (ck->deg_fixe) free(ck->deg_fixe);
  if (ck->sigma_gauss) free(ck->sigma_gauss);
  if (ck->filter_x) free(ck->filter_x);
  if (ck->filter_y) free(ck->filter_y);
  if (ck->kernel) free(ck->kernel);
  if (ck->kernel_coeffs) free(ck->kernel_coeffs);
  if (ck->solution) free(ck->solution);

  for (ig = 0; ig < ck->ngauss; ig++) {
    for (idegx = 0; idegx <= ck->deg_fixe[ig]; idegx++) {
      for (idegy = 0; idegy <= ck->deg_fixe[ig] - idegx; idegy++) {
        /* stores kernel weight mask for each order */
        if (ck->kernel_vec[nvec]) free(ck->kernel_vec[nvec]);
        nvec++;
      }
    }
  }
  if (ck->kernel_vec) free(ck->kernel_vec);

  if (ck->rxmins) free(ck->rxmins);
  if (ck->rxmaxs) free(ck->rxmaxs);
  if (ck->rymins) free(ck->rymins);
  if (ck->rymaxs) free(ck->rymaxs);
  if (ck->ksumims) free(ck->ksumims);
  for (i = 0; i < ck->nreg; ++i)
    if (ck->solutions[i]) free(ck->solutions[i]);
  if (ck->solutions) free(ck->solutions);
  return;
}


int select_kernel_solution(Kernel *ck, int region, int *status)
{
  /**
   * Set kernel solution to the one for specified region
   */
  if (! (region >= 0 && region < ck->nreg)) {
    fprintf(stderr, "Region %d not defined", region);
    *status = 1;
  }
  else {
    ck->region = region;
    ck->solution = (ck->solutions[ck->region]); 
    ck->ksumim = &(ck->ksumims[ck->region]);
    *status = 0;
  }
  return *status;
}


int get_image_size(char *image, long *axes, int *status)
{
  /**
   * Get the size of an image
   */
  fitsfile *pf;
  if (! fits_open_file(&pf, image, 0, status)) {
    if (! fits_get_img_param(pf, 2, NULL, NULL, axes, status)) {
      fits_close_file(pf, status);
    }
  }
  return *status;
}


int find_region_in_kernel(Kernel *ck, float x, float y)
{
  /**
   * Return the region index for given coordinates
   */
  int i, rxmin, rxmax, rymin, rymax;
  x += 1;
  y += 1;
  for (i = 0; i < ck->nreg; ++i) {
    rxmin = ck->rxmins[i];
    rxmax = ck->rxmaxs[i];
    rymin = ck->rymins[i];
    rymax = ck->rymaxs[i];
    if (x >= rxmin && x <= rxmax && y >= rymin && y <= rymax)
      return i;
  }
  return -1;
}


/**
 * from alard.c 
 */

double make_kernel(Kernel *ck, int xi, int yi, int xsize, int ysize)
{
  /**
   * Create the appropriate kernel (accessible as ck->kernel) at xi,
   * yi within xsize by ysize image size in which kernel solution is
   * defined.  The kernel sum is returned.
   */
  int ix, iy, i, j;
  double ax, ay;
  int k = 2;
  /* solution is mapped to range from -1 to 1 (covers the entire
     image, not individual region) */
  double xf = (xi - 0.5 * xsize) / (0.5 * xsize);
  double yf = (yi - 0.5 * ysize) / (0.5 * ysize);
  
  for (j = 1; j < ck->ncompker; ++j) {
    ck->kernel_coeffs[j] = 0.0;
    ax = 1.0;
    for (ix = 0; ix <= ck->order; ++ix) {
      ay = 1.0;
      for (iy = 0; iy <= ck->order - ix; ++iy) {
        ck->kernel_coeffs[j] += ck->solution[k++] * ax * ay;
        ay *= yf;
      }
      ax *= xf;
    }
  }
  ck->kernel_coeffs[0] = ck->solution[1]; 
  
  for (i = 0; i < ck->fw * ck->fw; ++i)
    ck->kernel[i] = 0.0;

  double ksum = 0.0;
  for (i = 0; i < ck->fw * ck->fw; ++i) {
    for (j = 0; j < ck->ncompker; ++j)
      ck->kernel[i] += ck->kernel_coeffs[j] * ck->kernel_vec[j][i];
    ksum += ck->kernel[i];
  }
  return ksum;
}


void spatial_convolve_stamp(float *image, int xsize, int ysize, float *outim, Kernel *ck,
                            float xconv, float yconv)
{
  /**
   * Take image of xsize by ysize and convolve it using the kernelSol
   * every kernel width
   */
  int hw = (ck->fw - 1) / 2;
  int kcstep = ck->fw; /* kcstep (-kcs option) is the size of step for spatial convolution */
  int nsteps_x = ceil((double)xsize / (double)kcstep);
  int nsteps_y = ceil((double)ysize / (double)kcstep);
  
  /* these are x & y offsets if a single kernel stamp is generated. */
  int x0 = max(0, (int)(xconv - xsize / 2));
  int y0 = max(0, (int)(yconv - ysize / 2));
  
  fprintf(stderr, "%d %d %d %d\n", xsize, ysize, x0, y0);

  int i, j, i0, j0, i1, j1, i2, j2, ic, jc, ik, jk;

  /* j1 & i1 are indices that iterate over y and x by step of kcstep
     size. Note that +hw adds padding to the edges. j0 and i0 are
     lower edges of each step. */
  for (j1 = 0; j1 < nsteps_y; ++j1) {
    j0 = j1 * kcstep + hw;
    for(i1 = 0; i1 < nsteps_x; ++i1) {
      i0 = i1 * kcstep + hw;

      /* want to evaluate kernel at the center of step.  originally,
         coordinates were shifted by hw to do this, but that can
         become way off if kcstep != 2*hw + 1. So here it's shifted by
         kcstep / 2 which is more generally correct.  */
      make_kernel(ck, i0 + kcstep/2 + x0, j0 + kcstep/2 + y0,
                  xsize, ysize);

      /* j2 & i2 iterates over pixels within each step */
      /* j & i points to pixels in unconvolved image */
      for (j2 = 0; j2 < kcstep; ++j2) {
        j = j0 + j2;
        if (j >= ysize - hw)
          break;
        for (i2 = 0; i2 < kcstep; ++i2) {
          i = i0 + i2;
          if (i >= xsize - hw)
            break;

          /* jc and ic points to pixels in unconvolved image */
          /* jk and ik points to kernel value */
          double q = 0.0;
          for (jc = j - hw; jc <= j + hw; ++jc) {
            jk = j - jc + hw;
            for (ic = i - hw; ic <= i + hw; ++ic) {
              ik = i - ic + hw;
              q += image[ic + xsize * jc] * ck->kernel[ik + jk * ck->fw];
            }
          }
          outim[i + xsize * j] = q;
	    }
      }
    }
  }
  return;
}


void spatial_convolve(float *image, int xsize, int ysize, float *outim, Kernel *ck)
{
  /**
   * Take image of xsize by ysize and convolve it using the kernelSol
   * every kernel width
   */
  int hw = (ck->fw - 1) / 2;
  int kcstep = ck->fw; /* kcstep (-kcs option) is the size of step for spatial convolution */
  int nsteps_x = ceil((double)xsize / (double)kcstep);
  int nsteps_y = ceil((double)ysize / (double)kcstep);
  
  int i, j, i0, j0, i1, j1, i2, j2, ic, jc, ik, jk;

  /* j1 & i1 are indices that iterate over y and x by step of kcstep
     size. Note that +hw adds padding to the edges. j0 and i0 are
     lower edges of each step. */
  for (j1 = 0; j1 < nsteps_y; ++j1) {
    j0 = j1 * kcstep + hw;
    for(i1 = 0; i1 < nsteps_x; ++i1) {
      i0 = i1 * kcstep + hw;

      /* want to evaluate kernel at the center of step.  originally,
         coordinates were shifted by hw to do this, but that can
         become way off if kcstep != 2*hw + 1. So here it's shifted by
         kcstep / 2 which is more generally correct.  */
      make_kernel(ck, i0 + kcstep/2, j0 + kcstep/2, xsize, ysize);

      /* j2 & i2 iterates over pixels within each step */
      /* j & i points to pixels in unconvolved image */
      for (j2 = 0; j2 < kcstep; ++j2) {
        j = j0 + j2;
        if (j >= ysize - hw)
          break;
        for (i2 = 0; i2 < kcstep; ++i2) {
          i = i0 + i2;
          if (i >= xsize - hw)
            break;

          /* jc and ic points to pixels in unconvolved image */
          /* jk and ik points to kernel value */
          double q = 0.0;
          for (jc = j - hw; jc <= j + hw; ++jc) {
            jk = j - jc + hw;
            for (ic = i - hw; ic <= i + hw; ++ic) {
              ik = i - ic + hw;
              q += image[ic + xsize * jc] * ck->kernel[ik + jk * ck->fw];
            }
          }
          outim[i + xsize * j] = q;
	    }
      }
    }
  }
  return;
}


int main(int argc, char **argv)
{
  char *iminput = NULL; /* inConv */
  long onaxes[2];
  int dofullimage = 0;
  float xconv = -1, yconv = -1; /* used to be initialized to 1e10 */

  int iarg, i, j, xsize, ysize, ndelta;
  int normalize = 0;
  int status = 0;
  float numKerFW = 2;
  double kSum;
  char *diffim = NULL, *outConv = NULL;
  float *delta = NULL, *dconv = NULL;
  long cNaxes[2];
  char scrStr[256], help[4096];
  fitsfile *pf;

  Kernel ck;

  sprintf(help, "Usage : extractkern [options] diffimage outimage\n");
  sprintf(help, "%sOptions:\n", help);
  sprintf(help, "%s   [-xy x y]        : convolve kernel with delta function at x,y\n", help);
  sprintf(help, "%s   [-nkw numkwidth] : # kernel widths for outimage size (%.1f)\n", help, numKerFW);
  sprintf(help, "%s   [-a]             : sample entire diffimage size with delta functions\n", help);
  sprintf(help, "%s   [-im image]      : convolve fitsfile instead of delta function\n", help);
  sprintf(help, "%s   [-n]             : divide convolved image by kernel sum\n\n", help);
  sprintf(help, "%s   To be used in conjuntion with the diffimage produced by hotpants\n", help);   
  sprintf(help, "%s      using the -hki option.  [-xy] convolves a delta function at\n", help);
  sprintf(help, "%s      the image position x, y with the spatially varying kernel\n", help);
  sprintf(help, "%s      used in the hotpants convolution.  Provides a visual realization \n", help);
  sprintf(help, "%s      of the kernel at that position, and can be useful for cosmic ray\n", help);
  sprintf(help, "%s      discrimination.  Also, if used with the [-im] option, one may\n", help);
  sprintf(help, "%s      reconstruct the entire convolved image to avoid storing it on disk.\n", help);

  /* read in command options. j counts # of required args given */
  for (iarg = 1, j = 0; iarg < argc; iarg++) {
    if (argv[iarg][0] == '-') {
      if (strcasecmp(argv[iarg] + 1, "xy") == 0) {
        sscanf(argv[++iarg], "%f", &xconv);
        sscanf(argv[++iarg], "%f", &yconv);
      }
      else if (strcasecmp(argv[iarg] + 1, "nkw") == 0) {
        sscanf(argv[++iarg], "%f", &numKerFW);
      }
      else if (strcasecmp(argv[iarg] + 1, "a") == 0) {
        dofullimage = 1;
      }
      else if (strcasecmp(argv[iarg] + 1, "im") == 0) {
        iminput = argv[++iarg];
      }
      else if (strcasecmp(argv[iarg] + 1, "n") == 0) {
        normalize = 1;
      }
      else {
        fprintf(stderr, "Unknown option %s\n", argv[iarg]);
        exit(1);
      }
    }
    else {
      diffim = argv[iarg++];
      outConv = argv[iarg++];
    }
  }
  if (iarg < 2) {
    /* not enough command line images...*/
    fprintf(stderr, "%s\n", help);
    exit(1);
  }

  /* insanity checking */
  if ( (xconv < 0 || yconv < 0) && ! dofullimage && ! iminput ) {
    fprintf(stderr, "Sorry, I do not know what to do, exiting...\n");
    exit(1);
  }

  if (init_kernel(diffim, &ck, &status)) {
    prgexit(status, &ck);
  }

  if (select_kernel_solution(&ck, 0, &status)) {
    prgexit(status, &ck);
  }

  if (get_image_size(diffim, onaxes, &status)) {
    prgexit(status, &ck);
  }

  /* Set output image size and shape */
  if (iminput) {
    /* input image to be convolved is given */
    if (fits_open_file(&pf, iminput, 0, &status)
        || fits_get_img_param(pf, 2, NULL, NULL, cNaxes, &status))
      fitserr(status, &ck);
    
    delta = (float*)malloc(cNaxes[0] * cNaxes[1] * sizeof(float));
    if (fits_read_img_flt(pf, 1, 1, cNaxes[0] * cNaxes[1], 0, delta, 0, &status)
        || fits_close_file(pf, &status))
      fitserr(status, &ck);
    
    xsize = cNaxes[0];
    ysize = cNaxes[1];
    onaxes[0] = cNaxes[0];
    onaxes[1] = cNaxes[1];
  }
  else {
    if (! dofullimage) {
      xsize = numKerFW * ck.fw;
      ysize = numKerFW * ck.fw;
      onaxes[0] = xsize;
      onaxes[1] = ysize;
    }
    else {
      xsize = onaxes[0];
      ysize = onaxes[1];
    }

    delta = (float*)calloc(xsize * ysize, sizeof(float));
    
    ndelta = 0;
    for (j = ck.fw - 1; j < (ysize - 1); j += ck.fw) {
      for (i = ck.fw - 1; i < (xsize - 1); i += ck.fw) {
        /* delta function in middle */
        delta[i + xsize * j] = 1.;
        ndelta += 1;
      }
    }
  }

  /* output convolved image */
  dconv = (float *)calloc(xsize * ysize, sizeof(float));

  /* do the convolution, fills kernel and calls make_kernel */
  spatial_convolve(delta, xsize, ysize, dconv, &ck);

  if (normalize) {
    /* TS: divide the convolved image by the kernel sum to preserve flux if needed */
    for (i = ck.fw / 2; i < ysize - ck.fw / 2; i++) {
      for (j = ck.fw / 2; j < xsize - ck.fw / 2; j++) {
        dconv[j + xsize * i] /= (*ck.ksumim);
      }
    }
  }

  /* clobber output image */
  sprintf(scrStr, "!%s", outConv);
  /* create and open new empty output FITS file, using input image as template.*/
  if (fits_create_file(&pf, scrStr, &status)
      || fits_create_img(pf, FLOAT_IMG, 2, onaxes, &status)
      || fits_write_img_flt(pf, 1, 1, xsize * ysize, dconv, &status)
      || fits_close_file(pf, &status))
    fitserr(status, &ck);
   
  /* sanity check - add up all pixels in the image */
  if (! (iminput)) {
    kSum = 0;
    for (i = 0; i < xsize * ysize; i++) {
      kSum += dconv[i];
    }
    fprintf(stderr, " Actual sum of pixels in convolved image : %.6f\n", kSum);
  }
  fprintf(stderr, " Kernel Sum from input image             : %.6f\n", (*ck.ksumim));

  /* free memories */
  if (delta) free(delta);
  if (dconv) free(dconv);

  release_kernel(&ck);
  return 0;
}
