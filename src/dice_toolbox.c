/*-----------------------------------------------------------------------------
  /
  / Filename: dice_init.c
  / Author: Valentin Perret
  / Author's email: perret.valentin@gmail.com
  / Description: DICE creates galaxies.
  /
  /	       DICE uses the GNU Scientific Library (GSL). You can
  /	       download the GSL source code from:
  /
  /		http://www.gnu.org/software/gsl
  /
  /	       or replace it with another math library.
  /
  / Copyright Information:
  /
  / Copyright (c) 2014       Valentin Perret
  /
  / This program is free software; you can redistribute it and/or modify
  / it under the terms of the GNU General Public License as published by
  / the Free Software Foundation; either version 2 of the License, or
  / (at your option) any later version.
  /
  / This program is distributed in the hope that it will be useful,
  / but WITHOUT ANY WARRANTY; without even the implied warranty of
  / MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  / GNU General Public License for more details.
  /
  / You should have received a copy of the GNU General Public License
  / along with this program; if not, write to the Free Software
  / Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  /
  / The license is also available at:
  /
  /		http://www.gnu.org/copyleft/gpl.html .
  /
  / Date: September 2015
  /
  */

#include "dice.h"

// Find the minimum between two values
double min(double a, double b) {
  if (a < b)
    return a;
  else
    return b;
}

// Find the maximum between two values
double max(double a, double b) {
  if (a > b)
    return a;
  else
    return b;
}

// Sum a 1D array
double sum_dbl(double *tab, int length) {
  int i;
  double sum;

  sum = 0.;
  for (i = 0; i < length; i++) {
    sum += tab[i];
  }
  return sum;
}

// Derivate a function using a 4-point central scheme
double deriv_central4(galaxy *gal, double x, double h, function_to_derivate F) {
  double new_x, f1, f2, f3, f4, derivative;

  new_x = x + 2.0 * h;
  f1 = F(new_x, gal);
  new_x = x + h;
  f2 = F(new_x, gal);
  new_x = x - h;
  f3 = F(new_x, gal);
  new_x = x - 2.0 * h;
  f4 = F(new_x, gal);
  derivative = (-f1 + 8.0 * f2 - 8.0 * f3 + f4) / (12.0 * h);

  return derivative;
}

// Derivate a function using a 4-point central scheme
double deriv_central2(galaxy *gal, double x, double h, function_to_derivate F) {
  double new_x, f1, f2, f3, f4, derivative;

  new_x = x + h;
  f1 = F(new_x, gal);
  new_x = x - h;
  f2 = F(new_x, gal);

  derivative = (f1 - f2) / (2.0 * h);

  return derivative;
}

// Derivate a function using a 2-point forward scheme
double deriv_forward(galaxy *gal, double x, double h, function_to_derivate F) {
  double f1, f2, f3, f4, derivative;

  f1 = F(x, gal);
  f2 = F(x + h, gal);

  derivative = (f2 - f1) / (h);

  return derivative;
}

double standard_dev(double *array, unsigned long int start,
                    unsigned long int end) {
  double stddev, mean;
  unsigned long int i;

  mean = 0.;
  for (i = start; i < end; i++) {
    mean += array[i];
  }
  mean /= (double)(end - start);
  stddev = 0.;
  for (i = start; i < end; i++) {
    stddev += pow(array[i], 2);
  }
  stddev = sqrt(stddev / (double)(end - start));

  return stddev;
}

// Linear interpolation between two points
double interpol(double x0, double x1, double x, double y0, double y1) {
  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

// Smooth in to zero function
double smooth_in(double x, double xc, double sigma) {
  return 0.5 * (1 + erf((x / xc - 1.0) / (sigma * sqrt(2))));
}

// Smooth out to zero function
double smooth_out(double x, double xc, double sigma) {
  return 1 - 0.5 * (1 + erf((x / xc - 1.0) / (sigma * sqrt(2))));
}

// This function calculates the potential due to a galactic disk using Cloud-
// In-Cell mass assignment with vacuum (isolated) boundary conditions on a
// Cartesian grid. The method was adapted from the discussion found in Hockney
// and Eastwood, "Computer Simulation Using Particles," 1981.
int set_galaxy_random_field_grid(galaxy *gal, double scale_inj,
                                 double scale_diss, double nspec, long seed) {
  unsigned long int ii, start, end, l;
  // Loop variables
  int i, j, k;
  // Nodes coordinates
  int node_x, node_y, node_z;
  int ngrid_padded[3];
  int sc;
  double dx, dy, dz, tx, ty, tz, n;
  double x, y, z, r, theta;
  double S_diss, S_inj, K_diss, K_inj, stddev, mean;
  // The particle-mesh grid size and the Green's function and potential
  // storage buffers. Global to keep from hitting the stack limit for large grid
  // size.
  double ***kernel_grid;
  // size_t p_x[gal->num_part[1]], p_y[gal->num_part[1]], p_z[gal->num_part[1]];
  fftw_plan fft_kernel, fft_gaussian_field, fftinv_gaussian_field;
  fftw_complex *kernel, *gaussian_field;

  ngrid_padded[0] = gal->ngrid_gauss[0];
  ngrid_padded[1] = gal->ngrid_gauss[1];
  ngrid_padded[2] = gal->ngrid_gauss[2];

  // Security
  if (scale_inj < scale_diss)
    scale_diss = scale_inj;
  S_diss = (scale_diss / gal->dx_gauss);
  S_inj = (scale_inj / gal->dx_gauss);
  K_diss = ngrid_padded[0] / S_diss;
  K_inj = ngrid_padded[0] / S_inj;

  // Setup fftw threads
#if USE_THREADS == 1
  // if (verbose) printf("/////\t\t-Setting up FFT call with %d
  // threads\n",AllVars.Nthreads);
  fflush(stdout);
  fftw_init_threads();
  fftw_plan_with_nthreads(AllVars.Nthreads);
  fflush(stdout);
#endif

  gsl_rng *rng;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, seed);

  // Allocate grid storage variables
  if (!(kernel = calloc(ngrid_padded[0] * ngrid_padded[1] * ngrid_padded[2],
                        sizeof(fftw_complex)))) {
    fprintf(stderr, "Unable to allocate space for kernel's function.\n");
    return -1;
  }
  if (!(gaussian_field =
            calloc(ngrid_padded[0] * ngrid_padded[1] * ngrid_padded[2],
                   sizeof(fftw_complex)))) {
    fprintf(stderr, "Unable to allocate space for gaussian_field buffer.\n");
    return -1;
  }

  if (!(kernel_grid = calloc(ngrid_padded[0], sizeof(double *)))) {
    fprintf(stderr, "Unable to create kernel's function x axis.\n");
    return -1;
  }
  // kernel grid allocation
  // x-axis
  for (i = 0; i < ngrid_padded[1]; ++i) {
    // y-axis
    if (!(kernel_grid[i] = calloc(ngrid_padded[1], sizeof(double *)))) {
      fprintf(stderr, "Unable to create kernel's function y axis.\n");
      return -1;
    }
    // z-axis
    for (j = 0; j < ngrid_padded[2]; ++j) {
      if (!(kernel_grid[i][j] = calloc(ngrid_padded[2], sizeof(double)))) {
        fprintf(stderr, "Unable to create kernel's function z axis.\n");
        return -1;
      }
    }
  }

  // Allocate the fftw complex output value and the fftw dft plan.
  fft_gaussian_field = fftw_plan_dft_3d(
      ngrid_padded[0], ngrid_padded[1], ngrid_padded[2], gaussian_field,
      gaussian_field, FFTW_FORWARD, FFTW_ESTIMATE);
  fft_kernel =
      fftw_plan_dft_3d(ngrid_padded[0], ngrid_padded[1], ngrid_padded[2],
                       kernel, kernel, FFTW_FORWARD, FFTW_ESTIMATE);

  // Normalization constant
  // See FFTW reference guide for more details
  n = (int)(ngrid_padded[0] * ngrid_padded[1] * ngrid_padded[2]);

  // Check for bad grids
  if (gal->ngrid_gauss[0] <= 0) {
    fprintf(stderr,
            "\t\tGrid dimensions must be greater than zero! (ngrid=%d)\n",
            gal->ngrid_gauss[0]);
    return -1;
  }

  // Initialization loop
  for (i = 0; i < ngrid_padded[0]; ++i) {
    for (j = 0; j < ngrid_padded[1]; ++j) {
      for (k = 0; k < ngrid_padded[2]; ++k) {
        gal->gaussian_field[i][j][k] = gsl_ran_gaussian(rng, 1.0);
        sc = gal->selected_comp[0];
        if (gal->comp_turb_gradient[sc] == 1) {
          x = ((double)i - (double)(ngrid_padded[0] / 2) + 0.5) * gal->dx_gauss;
          y = ((double)j - (double)(ngrid_padded[1] / 2) + 0.5) * gal->dx_gauss;
          z = ((double)k - (double)(ngrid_padded[2] / 2) + 0.5) * gal->dx_gauss;
          r = sqrt(x * x + y * y);
          theta = atan2(y, x);
          gal->gaussian_field[i][j][k] *= density_functions_pool(
              gal, r, theta, z, 1, gal->comp_model[sc], sc);
        }
      }
    }
  }

  if (scale_inj < gal->dx_gauss) {
    // Free fftw plan.
    // Kill the storage arrays since they are no longer needed.
    fftw_free(kernel);
    fftw_free(gaussian_field);
    for (i = 0; i < ngrid_padded[1]; ++i) {
      for (j = 0; j < ngrid_padded[2]; ++j) {
        free(kernel_grid[i][j]);
      }
      free(kernel_grid[i]);
    }
    free(kernel_grid);
    gsl_rng_free(rng);

#if USE_THREADS == 1
    fftw_cleanup_threads();
#endif
    // Flag the galaxy structure
    return 0;
  }
  // Define kernel's function.
  // These are the grid points as measured from the center of kernel's function
  // and the local value of the truncation function. The density is also packed
  // into a buffer here.
  //
  // kernel's function is defined eight times here, once per octant, to take
  // care of the isolated (vacuum) boundary conditions. See Hockney and
  // Eastwood, 1980, ch. 6 for a discussion. The octants start in the lower left
  // at (p,q) = (0,0) and progress counter clockwise.
  for (i = 0; i < ngrid_padded[0] / 2; ++i) {
    for (j = 0; j < ngrid_padded[1] / 2; ++j) {
#pragma omp parallel for private(dx, dy, dz) shared(kernel_grid, i, j)
      for (k = 0; k < ngrid_padded[2] / 2; ++k) {
        dx = (double)(i + 0.5);
        dy = (double)(j + 0.5);
        dz = (double)(k + 0.5);
        // Octant 1
        // Use a power spectrum
        // define wave number function
        double K_loc = sqrt(dx * dx + dy * dy + dz * dz);
        kernel_grid[i][j][k] = pow(K_loc, nspec / 2.0);
        // Suppress high & low frequencies
        if (K_loc <= K_inj) {
          kernel_grid[i][j][k] = pow(K_inj, nspec / 2.0);
          kernel_grid[i][j][k] =
              kernel_grid[i][j][k] *
              exp(-pow(K_loc - K_inj, 2.0) / (2.0 * pow(S_inj, 2.0)));
        }
        if (K_loc >= K_diss) {
          kernel_grid[i][j][k] = pow(K_diss, nspec / 2.0);
          kernel_grid[i][j][k] =
              kernel_grid[i][j][k] *
              exp(-pow(K_loc - K_diss, 2.0) / (2.0 * pow(S_diss, 2.0)));
        }
        // Octant 2
        kernel_grid[ngrid_padded[0] - 1 - i][j][k] = kernel_grid[i][j][k];
        // Octant 3
        kernel_grid[ngrid_padded[0] - 1 - i][ngrid_padded[1] - 1 - j][k] =
            kernel_grid[i][j][k];
        // Octant 4
        kernel_grid[i][ngrid_padded[1] - 1 - j][k] = kernel_grid[i][j][k];
        // Octant 5
        kernel_grid[i][j][ngrid_padded[2] - 1 - k] = kernel_grid[i][j][k];
        // Octant 6
        kernel_grid[ngrid_padded[0] - 1 - i][j][ngrid_padded[2] - 1 - k] =
            kernel_grid[i][j][k];
        // Octant 7
        kernel_grid[ngrid_padded[0] - 1 - i][ngrid_padded[1] - 1 - j]
                   [ngrid_padded[2] - 1 - k] = kernel_grid[i][j][k];
        // Octant 8
        kernel_grid[i][ngrid_padded[1] - 1 - j][ngrid_padded[2] - 1 - k] =
            kernel_grid[i][j][k];
      }
    }
  }

  // Pack kernel's function and the density into 1D arrays
  l = 0;
  for (i = 0; i < ngrid_padded[0]; ++i) {
    for (j = 0; j < ngrid_padded[1]; ++j) {
      for (k = 0; k < ngrid_padded[2]; ++k) {
        kernel[l] = kernel_grid[i][j][k];
        gaussian_field[l] = gal->gaussian_field[i][j][k];
        l++;
      }
    }
  }

  // Perform the fourier transforms. Density first, kernel's function second.
  fftw_execute(fft_gaussian_field);
  // Do not run fftw_execute(fft_kernel) because the kernel was built in Fourier
  // space FFT is computed, we can free the memory
  fftw_destroy_plan(fft_gaussian_field);
  fftw_destroy_plan(fft_kernel);
  // Allocating memory for the inverse fourier computation
  fftinv_gaussian_field = fftw_plan_dft_3d(
      ngrid_padded[0], ngrid_padded[1], ngrid_padded[2], gaussian_field,
      gaussian_field, FFTW_BACKWARD, FFTW_ESTIMATE);
  // Multiply the density by kernel's function to find the k-space
  // gaussian_field and Multiply the density by kernel's function to find the
  // k-space gaussian_field and invert for the real potenital. Second, normalize
  // the system and, finally, put the gaussian_field information into the grid.
  for (i = 0; i < n; ++i) {
    // Convolve the gaussian_field
    gaussian_field[i] = kernel[i] * gaussian_field[i];
  }
  // Inversion
  fftw_execute(fftinv_gaussian_field);
  fftw_destroy_plan(fftinv_gaussian_field);
  // Normalization
  mean = 0.;
  for (i = 0; i < n; ++i) {
    // Convolve the gaussian_field
    gaussian_field[i] = gaussian_field[i] / n;
    // Compute mean and standard deviation
    mean += gaussian_field[i];
  }
  mean /= (double)n;
  stddev = 0.;
  for (i = 0; i < n; ++i) {
    gaussian_field[i] -= mean;
    stddev += pow(gaussian_field[i], 2.0);
  }
  stddev = sqrt(stddev / n);

  l = 0;
  for (i = 0; i < ngrid_padded[0]; ++i) {
    for (j = 0; j < ngrid_padded[1]; ++j) {
      for (k = 0; k < ngrid_padded[2]; ++k) {
        // Fix the grid info
        gal->gaussian_field[i][j][k] = gaussian_field[l] / stddev;
        l++;
      }
    }
  }
  // Free fftw plan.
  // Kill the storage arrays since they are no longer needed.
  fftw_free(kernel);
  fftw_free(gaussian_field);
  for (i = 0; i < ngrid_padded[1]; ++i) {
    for (j = 0; j < ngrid_padded[2]; ++j) {
      free(kernel_grid[i][j]);
    }
    free(kernel_grid[i]);
  }
  free(kernel_grid);
  gsl_rng_free(rng);

#if USE_THREADS == 1
  fftw_cleanup_threads();
#endif
  // Flag the galaxy structure
  return 0;
}

// This function calculates the potential due to a galactic disk using Cloud-
// In-Cell mass assignment with vacuum (isolated) boundary conditions on a
// Cartesian grid. The method was adapted from the discussion found in Hockney
// and Eastwood, "Computer Simulation Using Particles," 1981.
int set_stream_random_field_grid(stream *st, double scale_inj,
                                 double scale_diss, double nspec, long seed) {
  unsigned long int ii, start, end, l;
  // Loop variables
  int i, j, k;
  // Nodes coordinates
  int node_x, node_y, node_z;
  int ngrid_padded[3];
  int sc;
  double dx, dy, dz, tx, ty, tz, n;
  double x, y, z, r, theta;
  double S_diss, S_inj, K_diss, K_inj, stddev, mean;
  // The particle-mesh grid size and the Green's function and potential
  // storage buffers. Global to keep from hitting the stack limit for large grid
  // size.
  double ***kernel_grid;
  // size_t p_x[st->num_part[1]], p_y[st->num_part[1]], p_z[st->num_part[1]];
  fftw_plan fft_kernel, fft_gaussian_field, fftinv_gaussian_field;
  fftw_complex *kernel, *gaussian_field;

  ngrid_padded[0] = st->ngrid_gauss[0];
  ngrid_padded[1] = st->ngrid_gauss[1];
  ngrid_padded[2] = st->ngrid_gauss[2];

  // Security
  if (scale_inj < scale_diss)
    scale_diss = scale_inj;
  S_diss = (scale_diss / st->dx_gauss);
  S_inj = (scale_inj / st->dx_gauss);
  K_diss = ngrid_padded[0] / S_diss;
  K_inj = ngrid_padded[0] / S_inj;

  // Setup fftw threads
#if USE_THREADS == 1
  // if (verbose) printf("/////\t\t-Setting up FFT call with %d
  // threads\n",AllVars.Nthreads);
  fflush(stdout);
  fftw_init_threads();
  fftw_plan_with_nthreads(AllVars.Nthreads);
  fflush(stdout);
#endif

  gsl_rng *rng;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, seed);

  // Allocate grid storage variables
  if (!(kernel = calloc(ngrid_padded[0] * ngrid_padded[1] * ngrid_padded[2],
                        sizeof(fftw_complex)))) {
    fprintf(stderr, "Unable to allocate space for kernel's function.\n");
    return -1;
  }
  if (!(gaussian_field =
            calloc(ngrid_padded[0] * ngrid_padded[1] * ngrid_padded[2],
                   sizeof(fftw_complex)))) {
    fprintf(stderr, "Unable to allocate space for gaussian_field buffer.\n");
    return -1;
  }

  if (!(kernel_grid = calloc(ngrid_padded[0], sizeof(double *)))) {
    fprintf(stderr, "Unable to create kernel's function x axis.\n");
    return -1;
  }
  // kernel grid allocation
  // x-axis
  for (i = 0; i < ngrid_padded[1]; ++i) {
    // y-axis
    if (!(kernel_grid[i] = calloc(ngrid_padded[1], sizeof(double *)))) {
      fprintf(stderr, "Unable to create kernel's function y axis.\n");
      return -1;
    }
    // z-axis
    for (j = 0; j < ngrid_padded[2]; ++j) {
      if (!(kernel_grid[i][j] = calloc(ngrid_padded[2], sizeof(double)))) {
        fprintf(stderr, "Unable to create kernel's function z axis.\n");
        return -1;
      }
    }
  }

  // Allocate the fftw complex output value and the fftw dft plan.
  fft_gaussian_field = fftw_plan_dft_3d(
      ngrid_padded[0], ngrid_padded[1], ngrid_padded[2], gaussian_field,
      gaussian_field, FFTW_FORWARD, FFTW_ESTIMATE);
  fft_kernel =
      fftw_plan_dft_3d(ngrid_padded[0], ngrid_padded[1], ngrid_padded[2],
                       kernel, kernel, FFTW_FORWARD, FFTW_ESTIMATE);

  // Normalization constant
  // See FFTW reference guide for more details
  n = (int)(ngrid_padded[0] * ngrid_padded[1] * ngrid_padded[2]);

  // Check for bad grids
  if (st->ngrid_gauss <= 0) {
    fprintf(stderr,
            "\t\tGrid dimensions must be greater than zero! (ngrid=%d)\n",
            st->ngrid_gauss[0]);
    return -1;
  }
  fflush(stdout);

  // Initialization loop
  for (i = 0; i < ngrid_padded[0]; ++i) {
    for (j = 0; j < ngrid_padded[1]; ++j) {
      for (k = 0; k < ngrid_padded[2]; ++k) {
        st->gaussian_field[i][j][k] = gsl_ran_gaussian(rng, 1.0);
        sc = st->selected_comp[0];
        if (st->comp_turb_gradient[sc] == 1) {
          x = ((double)i - (double)(ngrid_padded[0] / 2) + 0.5) * st->dx_gauss;
          y = ((double)j - (double)(ngrid_padded[1] / 2) + 0.5) * st->dx_gauss;
          z = ((double)k - (double)(ngrid_padded[2] / 2) + 0.5) * st->dx_gauss;
          r = sqrt(x * x + y * y);
          theta = atan2(y, x);
          st->gaussian_field[i][j][k] *= density_functions_stream_pool(
              st, r, theta, z, st->comp_model[sc], sc);
        }
      }
    }
  }

  if (scale_inj < st->dx_gauss) {
    // Free fftw plan.
    // Kill the storage arrays since they are no longer needed.
    fftw_free(kernel);
    fftw_free(gaussian_field);
    for (i = 0; i < ngrid_padded[0]; ++i) {
      for (j = 0; j < ngrid_padded[0]; ++j) {
        free(kernel_grid[i][j]);
      }
      free(kernel_grid[i]);
    }
    free(kernel_grid);
    gsl_rng_free(rng);

#if USE_THREADS == 1
    fftw_cleanup_threads();
#endif
    // Flag the galaxy structure
    return 0;
  }

  // Define kernel's function.
  // These are the grid points as measured from the center of kernel's function
  // and the local value of the truncation function. The density is also packed
  // into a buffer here.
  //
  // kernel's function is defined eight times here, once per octant, to take
  // care of the isolated (vacuum) boundary conditions. See Hockney and
  // Eastwood, 1980, ch. 6 for a discussion. The octants start in the lower left
  // at (p,q) = (0,0) and progress counter clockwise.
  for (i = 0; i < ngrid_padded[0] / 2; ++i) {
    for (j = 0; j < ngrid_padded[1] / 2; ++j) {
#pragma omp parallel for private(dx, dy, dz) shared(kernel_grid, i, j)
      for (k = 0; k < ngrid_padded[2] / 2; ++k) {
        dx = (double)(i + 0.5);
        dy = (double)(j + 0.5);
        dz = (double)(k + 0.5);
        // Octant 1
        // Use a power spectrum
        // define wave number function
        double K_loc = sqrt(dx * dx + dy * dy + dz * dz);
        kernel_grid[i][j][k] = pow(K_loc, nspec / 2.0);
        // Suppress high & low frequencies
        if (K_loc <= K_inj) {
          kernel_grid[i][j][k] = pow(K_inj, nspec / 2.0);
          kernel_grid[i][j][k] =
              kernel_grid[i][j][k] *
              exp(-pow(K_loc - K_inj, 2.0) / (2.0 * pow(S_inj, 2.0)));
        }
        if (K_loc >= K_diss) {
          kernel_grid[i][j][k] = pow(K_diss, nspec / 2.0);
          kernel_grid[i][j][k] =
              kernel_grid[i][j][k] *
              exp(-pow(K_loc - K_diss, 2.0) / (2.0 * pow(S_diss, 2.0)));
        }
        // Octant 2
        kernel_grid[ngrid_padded[0] - 1 - i][j][k] = kernel_grid[i][j][k];
        // Octant 3
        kernel_grid[ngrid_padded[0] - 1 - i][ngrid_padded[1] - 1 - j][k] =
            kernel_grid[i][j][k];
        // Octant 4
        kernel_grid[i][ngrid_padded[1] - 1 - j][k] = kernel_grid[i][j][k];
        // Octant 5
        kernel_grid[i][j][ngrid_padded[2] - 1 - k] = kernel_grid[i][j][k];
        // Octant 6
        kernel_grid[ngrid_padded[0] - 1 - i][j][ngrid_padded[2] - 1 - k] =
            kernel_grid[i][j][k];
        // Octant 7
        kernel_grid[ngrid_padded[0] - 1 - i][ngrid_padded[1] - 1 - j]
                   [ngrid_padded[2] - 1 - k] = kernel_grid[i][j][k];
        // Octant 8
        kernel_grid[i][ngrid_padded[1] - 1 - j][ngrid_padded[2] - 1 - k] =
            kernel_grid[i][j][k];
      }
    }
  }

  // Pack kernel's function and the density into 1D arrays
  l = 0;
  for (i = 0; i < ngrid_padded[0]; ++i) {
    for (j = 0; j < ngrid_padded[1]; ++j) {
      for (k = 0; k < ngrid_padded[2]; ++k) {
        kernel[l] = kernel_grid[i][j][k];
        gaussian_field[l] = st->gaussian_field[i][j][k];
        l++;
      }
    }
  }

  // Perform the fourier transforms. Density first, kernel's function second.
  fftw_execute(fft_gaussian_field);
  // Do not run fftw_execute(fft_kernel) because the kernel was built in Fourier
  // space FFT is computed, we can free the memory
  fftw_destroy_plan(fft_gaussian_field);
  fftw_destroy_plan(fft_kernel);
  // Allocating memory for the inverse fourier computation
  fftinv_gaussian_field = fftw_plan_dft_3d(
      ngrid_padded[0], ngrid_padded[0], ngrid_padded[0], gaussian_field,
      gaussian_field, FFTW_BACKWARD, FFTW_ESTIMATE);
  // Multiply the density by kernel's function to find the k-space
  // gaussian_field and invert for the real potenital. Second, normalize the
  // system and, finally, put the gaussian_field information into the grid.
  for (i = 0; i < n; ++i) {
    // Convolve the gaussian_field
    gaussian_field[i] = kernel[i] * gaussian_field[i];
  }
  // Inversion
  fftw_execute(fftinv_gaussian_field);
  fftw_destroy_plan(fftinv_gaussian_field);
  // Normalization
  mean = 0.;
  for (i = 0; i < n; ++i) {
    // Convolve the gaussian_field
    gaussian_field[i] = gaussian_field[i] / n;
    // Compute mean and standard deviation
    mean += gaussian_field[i];
  }
  mean /= (double)n;
  stddev = 0.;
  for (i = 0; i < n; ++i) {
    gaussian_field[i] -= mean;
    stddev += pow(gaussian_field[i], 2.0);
  }
  stddev = sqrt(stddev / n);

  l = 0;
  for (i = 0; i < ngrid_padded[0]; ++i) {
    for (j = 0; j < ngrid_padded[1]; ++j) {
      for (k = 0; k < ngrid_padded[2]; ++k) {
        // Fix the grid info
        st->gaussian_field[i][j][k] = gaussian_field[l] / stddev;
        l++;
      }
    }
  }
  // Free fftw plan.
  // Kill the storage arrays since they are no longer needed.
  fftw_free(kernel);
  fftw_free(gaussian_field);
  for (i = 0; i < ngrid_padded[0]; ++i) {
    for (j = 0; j < ngrid_padded[0]; ++j) {
      free(kernel_grid[i][j]);
    }
    free(kernel_grid[i]);
  }
  free(kernel_grid);
  gsl_rng_free(rng);

#if USE_THREADS == 1
  fftw_cleanup_threads();
#endif
  // Flag the galaxy structure
  return 0;
}

// This function calculates the turbulence at a point x,y,z
// by interpolating between grid points on the particle mesh. The
// interpolation routine uses the CIC kernel, which oddly enough is just
// a bilinear interpolation scheme...
//
// If the point lies off of the particle mesh, it approximates the potential
// as a function of 1/r.
double galaxy_gaussian_field_func(galaxy *gal, double x, double y, double z) {

  int node_x, node_y, node_z, offset;
  int ngrid_padded[3];
  double gauss1, gauss2, gauss3, gauss4, gauss5, gauss6, gauss7, gauss8;
  double a, dx, dy, dz, tx, ty, tz, gaussian_field;
  double xtemp, ytemp, ztemp, theta, phi, xmax, ymax, zmax, rnorm;

  ngrid_padded[0] = gal->ngrid_gauss[0];
  ngrid_padded[1] = gal->ngrid_gauss[1];
  ngrid_padded[2] = gal->ngrid_gauss[2];

  // Scale the coordinates
  xtemp = x / gal->dx_gauss + ((double)(ngrid_padded[0] / 2) - 0.5 - 0.5);
  ytemp = y / gal->dx_gauss + ((double)(ngrid_padded[1] / 2) - 0.5 - 0.5);
  ztemp = z / gal->dx_gauss + ((double)(ngrid_padded[2] / 2) - 0.5 - 0.5);
  offset = 0;
  // Determine the parent node.
  node_x = floor(xtemp);
  node_y = floor(ytemp);
  node_z = floor(ztemp);

  if (node_x < 0 || node_x >= ngrid_padded[0] - 1)
    return 0.;
  if (node_y < 0 || node_y >= ngrid_padded[1] - 1)
    return 0.;
  if (node_z < 0 || node_z >= ngrid_padded[2] - 1)
    return 0.;

  // Check to see if (x,y,z) is a grid point.
  if (xtemp == (double)node_y && ytemp == (double)node_y &&
      ztemp == (double)node_z) {
    // If (x,y,z) is a grid point, return its potential.
    gaussian_field = gal->gaussian_field[node_x][node_y][node_z];
  } else {
    // If (x,y,z) is not a grid point, use the CIC
    // interpolation function to calculate the potential.
    gauss1 = gal->gaussian_field[node_x][node_y][node_z];
    gauss2 = gal->gaussian_field[node_x + 1][node_y][node_z];
    gauss3 = gal->gaussian_field[node_x][node_y + 1][node_z];
    gauss4 = gal->gaussian_field[node_x][node_y][node_z + 1];
    gauss5 = gal->gaussian_field[node_x][node_y + 1][node_z + 1];
    gauss6 = gal->gaussian_field[node_x + 1][node_y + 1][node_z];
    gauss7 = gal->gaussian_field[node_x + 1][node_y][node_z + 1];
    gauss8 = gal->gaussian_field[node_x + 1][node_y + 1][node_z + 1];
    // CIC fractions
    dx = 1.0 - (xtemp - (double)node_x);
    dy = 1.0 - (ytemp - (double)node_y);
    dz = 1.0 - (ztemp - (double)node_z);
    tx = 1.0 - dx;
    ty = 1.0 - dy;
    tz = 1.0 - dz;
    // Return the interpolated field
    gaussian_field = dx * dy * dz * gauss1 + tx * dy * dz * gauss2 +
                     dx * ty * dz * gauss3 + dx * dy * tz * gauss4 +
                     dx * ty * tz * gauss5 + tx * ty * dz * gauss6 +
                     tx * dy * tz * gauss7 + tx * ty * tz * gauss8;
  }
  return gaussian_field;
}

// This function calculates the turbulence at a point x,y,z
// by interpolating between grid points on the particle mesh. The
// interpolation routine uses the CIC kernel, which oddly enough is just
// a bilinear interpolation scheme...
//
// If the point lies off of the particle mesh, it approximates the potential
// as a function of 1/r.
double stream_gaussian_field_func(stream *st, double x, double y, double z) {

  int node_x, node_y, node_z, offset;
  int ngrid_padded[3];
  double gauss1, gauss2, gauss3, gauss4, gauss5, gauss6, gauss7, gauss8;
  double a, dx, dy, dz, tx, ty, tz, gaussian_field;
  double xtemp, ytemp, ztemp, theta, phi, xmax, ymax, zmax, rnorm;

  ngrid_padded[0] = 2 * st->ngrid_gauss[0];
  ngrid_padded[1] = 2 * st->ngrid_gauss[1];
  ngrid_padded[2] = 2 * st->ngrid_gauss[2];

  // Scale the coordinates
  xtemp = x / st->dx_gauss + ((double)(ngrid_padded[0] / 2) - 0.5 - 0.5);
  ytemp = y / st->dx_gauss + ((double)(ngrid_padded[1] / 2) - 0.5 - 0.5);
  ztemp = z / st->dx_gauss + ((double)(ngrid_padded[2] / 2) - 0.5 - 0.5);
  offset = 0;
  // Determine the parent node.
  node_x = floor(xtemp);
  node_y = floor(ytemp);
  node_z = floor(ztemp);

  if (node_x < 0 || node_x >= ngrid_padded[0] - 1)
    return 0.;
  if (node_y < 0 || node_y >= ngrid_padded[1] - 1)
    return 0.;
  if (node_z < 0 || node_z >= ngrid_padded[2] - 1)
    return 0.;

  // Check to see if (x,y,z) is a grid point.
  if (xtemp == (double)node_y && ytemp == (double)node_y &&
      ztemp == (double)node_z) {
    // If (x,y,z) is a grid point, return its potential.
    gaussian_field = st->gaussian_field[node_x][node_y][node_z];
  } else {
    // If (x,y,z) is not a grid point, use the CIC
    // interpolation function to calculate the potential.
    gauss1 = st->gaussian_field[node_x][node_y][node_z];
    gauss2 = st->gaussian_field[node_x + 1][node_y][node_z];
    gauss3 = st->gaussian_field[node_x][node_y + 1][node_z];
    gauss4 = st->gaussian_field[node_x][node_y][node_z + 1];
    gauss6 = st->gaussian_field[node_x + 1][node_y + 1][node_z];
    gauss7 = st->gaussian_field[node_x + 1][node_y][node_z + 1];
    gauss8 = st->gaussian_field[node_x + 1][node_y + 1][node_z + 1];
    // CIC fractions
    dx = 1.0 - (xtemp - (double)node_x);
    dy = 1.0 - (ytemp - (double)node_y);
    dz = 1.0 - (ztemp - (double)node_z);
    tx = 1.0 - dx;
    ty = 1.0 - dy;
    tz = 1.0 - dz;
    // Return the interpolated potential.
    gaussian_field = dx * dy * dz * gauss1 + tx * dy * dz * gauss2 +
                     dx * ty * dz * gauss3 + dx * dy * tz * gauss4 +
                     dx * ty * tz * gauss5 + tx * ty * dz * gauss6 +
                     tx * dy * tz * gauss7 + tx * ty * tz * gauss8;
  }
  return gaussian_field;
}
