#include<stdio.h>
#include<string.h>
#include<strings.h>
#include<math.h>
#include<stdlib.h>
#include<fitsio.h>


#define max(x,y) x>y?x:y
#define min(x,y) x<y?x:y


#define DEFAULT_KERNEL_FW 2.0


typedef struct
{
  int fw; /* fwKernel */
  int ngauss;
  int ncompker; /* nCompKer */
  int ncomp; /* nComp */
  int nbgvec; /* nBGVectors */
  int ncomptot; /* nCompTotal */
  int order; /* kerOrder */
  int bgorder; /* bgOrder */
  int kcstep; // kernel convolution step
  int fwstamp; // stamp full width

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
void release_kernel(Kernel *ck);


int imin(int a, int b) { return a < b ? a : b; }  
int imax(int a, int b) { return a < b ? b : a; }  


void prgexit(int status, Kernel *ck)
{
  if (status > 100) {
    /* error code > 100 is from CFITSIO */
    fits_report_error(stderr, status);
  }
  if (ck != NULL) release_kernel(ck);
  exit(status);
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
  int i, dummy, status = 0, info = 0;//, rxmin, rxmax, rymin, rymax;
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
    if (! fits_read_keys_str(pf, "REGION", 0, 9, regions, &(ck->nreg), &status)) {
      if (ck->nreg == 0) {
        fprintf(stderr, "no region defined\n");
        goto error_exit;
      }
      fprintf(stdout, "%d region(s) found\n", ck->nreg);
    }
    else goto error_exit;
    
    /* move to binary kernel table... */
    if (fits_movnam_hdu(pf, BINARY_TBL, "CONVOLUTION KERNEL INFORMATION", 0, &status))
      goto error_exit;

    if (fits_read_key(pf, TINT, "NGAUSS", &(ck->ngauss), NULL, &status) ||
        fits_read_key(pf, TINT, "FWKERN", &(ck->fw), NULL, &status) ||
        fits_read_key(pf, TINT, "CKORDER", &(ck->order), NULL, &status) ||
        fits_read_key(pf, TINT, "BGORDER", &(ck->bgorder), NULL, &status))
      goto error_exit;

    if (fits_read_key(pf, TINT, "HPKCS", &(ck->kcstep), NULL, &status)) {
      if (status == VALUE_UNDEFINED || status == KEY_NO_EXIST) {
        ck->kcstep = ck->fw;
        status = 0;
      }
      else
        goto error_exit;
    }
    if (fits_read_key(pf, TINT, "HPFWSTMP", &(ck->fwstamp), NULL, &status)) {
      if (status == VALUE_UNDEFINED || status == KEY_NO_EXIST) {
        ck->fwstamp = 50 * ck->fw;
        ck->fwstamp -= ck->fw;
        ck->fwstamp -= ck->fwstamp % 2 == 0 ? 1 : 0;
        status = 0;
      }
      else
        goto error_exit;
    }
    
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
      if (sscanf(regions[i], "[%d:%d,%d:%d]",
                 &(ck->rxmins[i]), &(ck->rxmaxs[i]),
                 &(ck->rymins[i]), &(ck->rymaxs[i])) != 4) {
        fprintf(stderr, "Problem with region %d (%s), exiting...\n", i, regions[i]);
        goto error_exit;
      }
      // 1-index to 0-index
      ck->rxmins[i]--; ck->rxmaxs[i]--; ck->rymins[i]--; ck->rymaxs[i]--;

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


double make_kernel(Kernel *ck, int xi, int yi, int xsize, int ysize)
{
  /**
   * Create the appropriate kernel (accessible as ck->kernel) at xi,
   * yi within xsize by ysize image size in which kernel solution is
   * defined.  The kernel sum is returned.
   */

  // solution is mapped to range from -1 to 1
  double xf = (xi - 0.5 * xsize) / (0.5 * xsize);
  double yf = (yi - 0.5 * ysize) / (0.5 * ysize);

  int k = 2;
  for (int j = 1; j < ck->ncompker; ++j) {
    ck->kernel_coeffs[j] = 0.0;
    double ax = 1.0;
    for (int ix = 0; ix <= ck->order; ++ix) {
      double ay = 1.0;
      for (int iy = 0; iy <= ck->order - ix; ++iy) {
        ck->kernel_coeffs[j] += ck->solution[k++] * ax * ay;
        ay *= yf;
      }
      ax *= xf;
    }
  }
  ck->kernel_coeffs[0] = ck->solution[1]; 
  
  for (int i = 0; i < ck->fw * ck->fw; ++i)
    ck->kernel[i] = 0.0;

  double ksum = 0.0;
  for (int i = 0; i < ck->fw * ck->fw; ++i) {
    for (int j = 0; j < ck->ncompker; ++j)
      ck->kernel[i] += ck->kernel_coeffs[j] * ck->kernel_vec[j][i];
    ksum += ck->kernel[i];
  }
  return ksum;
}


/* void spatial_convolve_stamp(float *image, int xsize, int ysize, float *outim, Kernel *ck, */
/*                             float xconv, float yconv) */
/* { */
/*   /\** */
/*    * Take image of xsize by ysize and convolve it using the kernelSol */
/*    * every kernel width */
/*    *\/ */
/*   int hw = (ck->fw - 1) / 2; */
/*   int kcstep = ck->fw; /\* kcstep (-kcs option) is the size of step for spatial convolution *\/ */
/*   int nsteps_x = ceil((double)xsize / (double)kcstep); */
/*   int nsteps_y = ceil((double)ysize / (double)kcstep); */
  
/*   /\* these are x & y offsets if a single kernel stamp is generated. *\/ */
/*   int x0 = max(0, (int)(xconv - xsize / 2)); */
/*   int y0 = max(0, (int)(yconv - ysize / 2)); */
  
/*   fprintf(stderr, "%d %d %d %d\n", xsize, ysize, x0, y0); */

/*   int i, j, i0, j0, i1, j1, i2, j2, ic, jc, ik, jk; */

/*   /\* j1 & i1 are indices that iterate over y and x by step of kcstep */
/*      size. Note that +hw adds padding to the edges. j0 and i0 are */
/*      lower edges of each step. *\/ */
/*   for (j1 = 0; j1 < nsteps_y; ++j1) { */
/*     j0 = j1 * kcstep + hw; */
/*     for(i1 = 0; i1 < nsteps_x; ++i1) { */
/*       i0 = i1 * kcstep + hw; */

/*       /\* want to evaluate kernel at the center of step.  originally, */
/*          coordinates were shifted by hw to do this, but that can */
/*          become way off if kcstep != 2*hw + 1. So here it's shifted by */
/*          kcstep / 2 which is more generally correct.  *\/ */
/*       make_kernel(ck, i0 + kcstep/2 + x0, j0 + kcstep/2 + y0, */
/*                   xsize, ysize); */

/*       /\* j2 & i2 iterates over pixels within each step *\/ */
/*       /\* j & i points to pixels in unconvolved image *\/ */
/*       for (j2 = 0; j2 < kcstep; ++j2) { */
/*         j = j0 + j2; */
/*         if (j >= ysize - hw) */
/*           break; */
/*         for (i2 = 0; i2 < kcstep; ++i2) { */
/*           i = i0 + i2; */
/*           if (i >= xsize - hw) */
/*             break; */

/*           /\* jc and ic points to pixels in unconvolved image *\/ */
/*           /\* jk and ik points to kernel value *\/ */
/*           double q = 0.0; */
/*           for (jc = j - hw; jc <= j + hw; ++jc) { */
/*             jk = j - jc + hw; */
/*             for (ic = i - hw; ic <= i + hw; ++ic) { */
/*               ik = i - ic + hw; */
/*               q += image[ic + xsize * jc] * ck->kernel[ik + jk * ck->fw]; */
/*             } */
/*           } */
/*           outim[i + xsize * j] = q; */
/* 	    } */
/*       } */
/*     } */
/*   } */
/*   return; */
/* } */


void spatial_convolve(float *image, int xsize, int ysize, float *outim, Kernel *ck, int kcstep)
{
  /**
   * Take image of xsize by ysize and convolve it using the kernelSol
   * every kernel width
   */
  int hw = (ck->fw - 1) / 2;
  //int kcstep = ck->fw; /* kcstep (-kcs option) is the size of step for spatial convolution */
  int nsteps_x = ceil((double)xsize / (double)kcstep);
  int nsteps_y = ceil((double)ysize / (double)kcstep);
  
  /* j1 & i1 are indices that iterate over y and x by step of kcstep
     size. Note that +hw adds padding to the edges. j0 and i0 are
     lower edges of each step. */
  for (int j1 = 0; j1 < nsteps_y; ++j1) {
    int j0 = j1 * kcstep + hw;
    for(int i1 = 0; i1 < nsteps_x; ++i1) {
      int i0 = i1 * kcstep + hw;

      /* want to evaluate kernel at the center of step.  originally,
         coordinates were shifted by hw to do this, but that can
         become way off if kcstep != 2*hw + 1. So here it's shifted by
         kcstep / 2 which is more generally correct.  */
      make_kernel(ck, i0 + kcstep/2, j0 + kcstep/2, xsize, ysize);

      /* j2 & i2 iterates over pixels within each step */
      /* j & i points to pixels in unconvolved image */
      for (int j2 = 0; j2 < kcstep; ++j2) {
        int j = j0 + j2;
        if (j >= ysize - hw)
          break;
        for (int i2 = 0; i2 < kcstep; ++i2) {
          int i = i0 + i2;
          if (i >= xsize - hw)
            break;

          /* jc and ic points to pixels in unconvolved image */
          /* jk and ik points to kernel value */
          double q = 0.0;
          for (int jc = j - hw; jc <= j + hw; ++jc) {
            int jk = j - jc + hw;
            for (int ic = i - hw; ic <= i + hw; ++ic) {
              int ik = i - ic + hw;
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


void show_help()
{
  char help[4096];
  sprintf(help, "Usage : extractkern [options] kernelimage outimage\n");
  sprintf(help, "%sOptions:\n", help);
  sprintf(help, "%s  [-xy x y]        : convolve kernel with delta function at x,y\n", help);
  sprintf(help, "%s  [-nkw numkwidth] : # kernel widths for outimage size (%.1f)\n", help, DEFAULT_KERNEL_FW);
  sprintf(help, "%s  [-a xsize ysize] : sample entire image size given with delta functions\n", help);
  sprintf(help, "%s  [-im image]      : convolve fitsfile instead of delta function\n", help);
  sprintf(help, "%s  [-n]             : divide convolved image by kernel sum\n\n", help);
  sprintf(help, "%s  To be used in conjuntion with the image with kernel info (kernelimage) produced by hotpants\n", help);   
  sprintf(help, "%s     using the -hki option.  [-xy] convolves a delta function at\n", help);
  sprintf(help, "%s     the image position x, y with the spatially varying kernel\n", help);
  sprintf(help, "%s     used in the hotpants convolution.  Provides a visual realization \n", help);
  sprintf(help, "%s     of the kernel at that position, and can be useful for cosmic ray\n", help);
  sprintf(help, "%s     discrimination.  Also, if used with the [-im] option, one may\n", help);
  sprintf(help, "%s     reconstruct the entire convolved image to avoid storing it on disk.\n", help);
  fprintf(stderr, "%s\n", help);
  exit(1);
  return;
}


int main(int argc, char **argv)
{
  int i, j, k, status = 0;
  fitsfile *pf;

  Kernel ck;
  char *kernelim = NULL, *iminput = NULL, *imoutput = NULL;
  int dofullimage = 0;
  int normalize = 0;
  float xconv = -1, yconv = -1;
  // output image dimension
  long onaxes[2] = {0, 0};

  float numKerFW = DEFAULT_KERNEL_FW;

  // read in command options. j counts # of required args given
  int iarg;
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
        if (sscanf(argv[++iarg], "%ld", &onaxes[0]) != 1)
          onaxes[0] = 0;
        if (sscanf(argv[++iarg], "%ld", &onaxes[1]) != 1)
          onaxes[1] = 0;
      }
      else if (strcasecmp(argv[iarg] + 1, "im") == 0) {
        iminput = argv[++iarg];
        onaxes[0] = 0;
        onaxes[1] = 0;
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
      kernelim = argv[iarg++];
      imoutput = argv[iarg++];
    }
  }
  // check input integrity
  if (iarg < 2)
    show_help(); // not enough command line args provided
  if ( (xconv < 0 || yconv < 0) && ! dofullimage && ! iminput )
    show_help(); // inconsistent input

  // prepare kernel
  if (init_kernel(kernelim, &ck, &status))
    prgexit(status, &ck);

  // need image dimension over which kernel has been defined!!
  // TODO: somehow we *must* get this info, but if the kernel info are
  // not stored in diff image, the following function may fail.
  if (! (onaxes[0] > 0 && onaxes[1] > 0)) {
    if (get_image_size((iminput ? iminput : kernelim), onaxes, &status))
      prgexit(status, &ck);
    if (! (onaxes[0] > 0 && onaxes[1] > 0)) {
      fprintf(stderr, "Output image dimension undefined\n");
      prgexit(status, &ck);
    }
  }

  // full image boundary indices (zero indexed)
  int xmin = 0, xmax = onaxes[0] - 1, ymin = 0, ymax = onaxes[1] - 1;

  // set output image data
  float *tbconv = NULL;  // to be convolved

  if (iminput) {
    // input image to be convolved is given
    if (! fits_open_file(&pf, iminput, 0, &status)) {
      tbconv = (float*)malloc(onaxes[0] * onaxes[1] * sizeof(float));
      if (fits_read_img_flt(pf, 1, 1, onaxes[0] * onaxes[1], 0, tbconv, 0, &status))
        goto error_exit;
      fits_close_file(pf, &status);
    }
    else {
      fits_close_file(pf, &status);
      goto error_exit;
    }
  }
  else {
    if (! dofullimage) {
      // TODO: this mode is not working yet
      onaxes[0] = numKerFW * ck.fw;
      onaxes[1] = numKerFW * ck.fw;
    }

    tbconv = (float*)calloc(onaxes[0] * onaxes[1], sizeof(float));
    // generate delta function at the center of each "stamp"
    for (j = ck.fw - 1; j < onaxes[1] - 1; j += ck.fw)
      for (i = ck.fw - 1; i < onaxes[0] - 1; i += ck.fw)
        tbconv[i + onaxes[0] * j] = 1.;
  }

  // storage for output convolved image
  float *dconv = NULL;
  dconv = (float *)calloc(onaxes[0] * onaxes[1], sizeof(float));

  // buffer size
  int hwbuf;
  if (ck.nreg > 1) {
    /*int fwstamp = imin(onaxes[0], onaxes[1]) / 23;
    fwstamp -= ck.fw;
    fwstamp -= fwstamp % 2 == 0 ? 1 : 0;*/
    hwbuf = ck.fwstamp / 2;
  }
  else {
    hwbuf = ck.fw / 2;
  }

  // iterate over regions
  for (int r = 0; r < ck.nreg; ++r) {
    if (select_kernel_solution(&ck, r, &status))
      prgexit(status, &ck);
    
    int rxmin = ck.rxmins[r], rxmax = ck.rxmaxs[r];
    int rymin = ck.rymins[r], rymax = ck.rymaxs[r];
    
    // buffered - add an extra half-stamp/-kernel size for merging regions
    int rxbmin, rxbmax, rybmin, rybmax;
    rxbmin = imax(xmin, rxmin - hwbuf);
    rybmin = imax(ymin, rymin - hwbuf);
    rxbmax = imin(xmax, rxmax + hwbuf);
    rybmax = imin(ymax, rymax + hwbuf);
    
    // buffer size used when merging output image sections
    int xbuflo = rxmin - rxbmin;
    int xbufhi = rxbmax - rxmax;
    int ybuflo = rymin - rybmin;
    int ybufhi = rybmax - rymax;

    // size of buffered region; important since the normalization of
    // coordinates for the kernel solution depends on these
    int rpixx = rxbmax - rxbmin + 1;
    int rpixy = rybmax - rybmin + 1;
    
    // limits of output images
    int fpixeloutx = rxbmin + xbuflo;
    int fpixelouty = rybmin + ybuflo;
    int lpixeloutx = fpixeloutx + (rpixx - xbufhi - xbuflo - 1);
    int lpixelouty = fpixelouty + (rpixy - ybufhi - ybuflo - 1);

    fprintf(stdout, "Region %2d : %d:%d, %d:%d\n", r, rxmin + 1, rxmax + 1, rymin + 1, rymax + 1);
    fprintf(stdout, "  buffered: %d:%d, %d:%d\n", rxbmin + 1, rxbmax + 1, rybmin + 1, rybmax + 1);
    fprintf(stdout, "  good pix: %d:%d, %d:%d\n", rxmin + 1, rxmax + 1, rymin + 1, rymax + 1);

    // working storage for this region
    float *wtbconv = (float *)calloc(rpixx * rpixy, sizeof(float));
    float *wdconv = (float *)calloc(rpixx * rpixy, sizeof(float));
    if (wtbconv == NULL || wdconv == NULL) {
      if (wtbconv) free(wtbconv);
      if (wdconv) free(wdconv);
      fprintf(stderr, "Out of memory\n");
      prgexit(status, &ck);
    }
    
    // read image
    for (j = rybmin, k = 0; j <= rybmax; ++j)
      for (i = rxbmin; i <= rxbmax; ++i)
        wtbconv[k++] = tbconv[i + j * onaxes[0]];
    
    // do the convolution
    spatial_convolve(wtbconv, rpixx, rpixy, wdconv, &ck, ck.kcstep);
    
    // copy image back to the parent
    //float ksum = normalize ? make_kernel(&ck, rpixx/2, rpixy/2, rpixx, rpixy) : 1.;
    float ksum = normalize ? (*ck.ksumim) : 1.;
    for (j = rybmin, k = 0; j <= rybmax; ++j) {
      if (j >= fpixelouty && j <= lpixelouty) {
        k += fpixeloutx - rxbmin;
        for (i = fpixeloutx; i <= lpixeloutx; ++i) {
          dconv[i + j * onaxes[0]] = wdconv[k] / ksum;
          k += 1;
        }
        k += rxbmax - lpixeloutx;
      }
      else {
        k += rpixx;
      }
    }

    free(wtbconv);
    free(wdconv);
  }

  // print out information

  // sanity check - add up all pixels in the image
  double ksum = 0;
  for (i = 0; i < onaxes[0] * onaxes[1]; ++i)
    ksum += dconv[i];
  fprintf(stderr, " Actual sum of pixels in convolved image : %.6f\n", ksum);
  fprintf(stderr, " Kernel Sum from input image             : %.6f\n", (*ck.ksumim));

  // generate output image; clobber existing one!!
  char scrstr[256];
  sprintf(scrstr, "!%s", imoutput);
  if (fits_create_file(&pf, scrstr, &status)
      || fits_create_img(pf, FLOAT_IMG, 2, onaxes, &status)
      || fits_write_img_flt(pf, 1, 1, onaxes[0] * onaxes[1], dconv, &status)
      || fits_close_file(pf, &status))
    prgexit(status, &ck);

  // cleaning up...
  if (tbconv) free(tbconv);
  if (dconv) free(dconv);
  release_kernel(&ck);
  return 0;

 error_exit:
  prgexit(status, &ck);
}
