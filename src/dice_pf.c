/*-----------------------------------------------------------------------------
  /
  / Filename: dice_pf.c
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

// This function calculates the potential due to a galactic disk using Cloud-
// In-Cell mass assignment with vacuum (isolated) boundary conditions on a
// Cartesian grid. The method was adapted from the discussion found in Hockney
// and Eastwood, "Computer Simulation Using Particles," 1981.
int set_galaxy_potential(galaxy *gal, double ***potential, double dx,
                         int ngrid[3], int verbose, double exclude[3]) {
  unsigned long int ii, l;
  // Loop variables
  int i, j, k;
  // Nodes coordinates
  int ngrid_padded[3];
  double n;
  double boxsize[3];
  // The particle-mesh grid size and the Green's function and potential
  // storage buffers. Global to keep from hitting the stack limit for large grid
  // size.
  double ***green_grid;
  fftw_plan fft_green, fft_rho, fftinv_potential;
  fftw_complex *green, *rho;

  // Setup fftw threads
#if USE_THREADS == 1
  fflush(stdout);
  fftw_init_threads();
  fftw_plan_with_nthreads(AllVars.Nthreads);
  fflush(stdout);
#endif

  ngrid_padded[0] = 2 * ngrid[0];
  ngrid_padded[1] = 2 * ngrid[1];
  ngrid_padded[2] = 2 * ngrid[2];
  boxsize[0] = ngrid[0] * dx;
  boxsize[1] = ngrid[1] * dx;
  boxsize[2] = ngrid[2] * dx;

  // Allocate grid storage variables
  if (!(green = calloc(ngrid_padded[0] * ngrid_padded[1] * ngrid_padded[2],
                       sizeof(fftw_complex)))) {
    fprintf(stderr, "[Error] Unable to allocate space for Green's function\n");
    return -1;
  }
  if (!(rho = calloc(ngrid_padded[0] * ngrid_padded[1] * ngrid_padded[2],
                     sizeof(fftw_complex)))) {
    fprintf(stderr, "[Error] Unable to allocate space for rho function\n");
    return -1;
  }

  if (!(green_grid = calloc(ngrid_padded[0], sizeof(double *)))) {
    fprintf(stderr, "[Error] Unable to create Green's function x axis\n");
    return -1;
  }
  // Green grid allocation
  // x-axis
  for (i = 0; i < ngrid_padded[0]; ++i) {
    // y-axis
    if (!(green_grid[i] = calloc(ngrid_padded[1], sizeof(double *)))) {
      fprintf(stderr, "[Error] Unable to create Green's function y axis\n");
      return -1;
    }
    // z-axis
    for (j = 0; j < ngrid_padded[1]; ++j) {
      if (!(green_grid[i][j] = calloc(ngrid_padded[2], sizeof(double)))) {
        fprintf(stderr, "[Error] Unable to create Green's function z axis\n");
        return -1;
      }
    }
  }
  // Allocate the fftw complex output value and the fftw dft plan.
  fft_rho = fftw_plan_dft_3d(ngrid_padded[0], ngrid_padded[1], ngrid_padded[2],
                             rho, rho, FFTW_FORWARD, FFTW_ESTIMATE);
  fft_green =
      fftw_plan_dft_3d(ngrid_padded[0], ngrid_padded[1], ngrid_padded[2], green,
                       green, FFTW_FORWARD, FFTW_ESTIMATE);

  // Normalization constant
  // See FFTW reference guide for more details
  n = (int)(ngrid_padded[0] * ngrid_padded[1] * ngrid_padded[2]);

  // Check for bad grids
  if (ngrid_padded[0] <= 0) {
    fprintf(stderr,
            "/////\t\t[Error] Grid dimensions must be greater than zero "
            "[ngrid=%d]\n",
            ngrid_padded[0]);
    return -1;
  }

  // Initialization loop
  for (i = 0; i < ngrid_padded[0]; ++i) {
    for (j = 0; j < ngrid_padded[1]; ++j) {
      for (k = 0; k < ngrid_padded[2]; ++k) {
        potential[i][j][k] = 0.;
      }
    }
  }
  // Calculate the density using the CIC routine. The positions are shifted
  // such that the particles are in the +x,+y,+z octant. space_* is added
  // to take care of the vacuum boundary conditions. The density values are
  // stored in the potential structure for now.
  for (ii = 0; ii < gal->ntot_part_pot; ++ii) {
    double x = gal->x[ii] / dx + ((double)(ngrid_padded[0] / 2) - 0.5);
    double y = gal->y[ii] / dx + ((double)(ngrid_padded[1] / 2) - 0.5);
    double z = gal->z[ii] / dx + ((double)(ngrid_padded[2] / 2) - 0.5);
    // Figure out which node owns the particle
    int node_x = (int)x;
    int node_y = (int)y;
    int node_z = (int)z;
    if (gal->x[ii] >= exclude[0] || gal->x[ii] <= -exclude[0] ||
        gal->y[ii] >= exclude[1] || gal->y[ii] <= -exclude[1] ||
        gal->z[ii] >= exclude[2] || gal->z[ii] <= -exclude[2]) {
      if (gal->x[ii] >= -boxsize[0] / 2. && gal->x[ii] <= boxsize[0] / 2. &&
          gal->y[ii] >= -boxsize[1] / 2. && gal->y[ii] <= boxsize[1] / 2. &&
          gal->z[ii] >= -boxsize[2] / 2. && gal->z[ii] <= boxsize[2] / 2.) {
        // Check if particle is not outside the potential grid
        if (node_x >= 0 && node_x < ngrid_padded[0] - 1 && node_y >= 0 &&
            node_y < ngrid_padded[1] - 1 && node_z >= 0 &&
            node_z < ngrid_padded[2] - 1) {
          // Set the CIC size fractions
          double d_x = 1.0 - (x - (double)node_x);
          double d_y = 1.0 - (y - (double)node_y);
          double d_z = 1.0 - (z - (double)node_z);
          double t_x = 1.0 - d_x;
          double t_y = 1.0 - d_y;
          double t_z = 1.0 - d_z;
          // Calculate the CIC densities
          potential[node_x][node_y][node_z] +=
              gal->mass[ii] * (d_x * d_y * d_z) / pow(dx, 3);
          potential[node_x + 1][node_y][node_z] +=
              gal->mass[ii] * (t_x * d_y * d_z) / pow(dx, 3);
          potential[node_x][node_y + 1][node_z] +=
              gal->mass[ii] * (d_x * t_y * d_z) / pow(dx, 3);
          potential[node_x][node_y][node_z + 1] +=
              gal->mass[ii] * (d_x * d_y * t_z) / pow(dx, 3);
          potential[node_x][node_y + 1][node_z + 1] +=
              gal->mass[ii] * (d_x * t_y * t_z) / pow(dx, 3);
          potential[node_x + 1][node_y + 1][node_z] +=
              gal->mass[ii] * (t_x * t_y * d_z) / pow(dx, 3);
          potential[node_x + 1][node_y][node_z + 1] +=
              gal->mass[ii] * (t_x * d_y * t_z) / pow(dx, 3);
          potential[node_x + 1][node_y + 1][node_z + 1] +=
              gal->mass[ii] * (t_x * t_y * t_z) / pow(dx, 3);
          // NGP (Nearest Point Grid)
          // potential[node_x][node_y][node_z] += gal->mass[ii]/(pow(dx,3));
        }
      }
    }
  }

  // Define Green's function.
  // These are the grid points as measured from the center of Green's function
  // and the local value of the truncation function. The density is also packed
  // into a buffer here.
  //
  // Green's function is defined eight times here, once per octant, to take care
  // of the isolated (vacuum) boundary conditions. See Hockney and Eastwood,
  // 1980, ch. 6 for a discussion. The octants start in the lower left at (p,q)
  // = (0,0) and progress counter clockwise.
  for (i = 0; i < ngrid_padded[0] / 2; ++i) {
    for (j = 0; j < ngrid_padded[1] / 2; ++j) {
      for (k = 0; k < ngrid_padded[2] / 2; ++k) {
        double d_x = sqrt(pow((double)(i + 0.5), 2.0)) * dx;
        double d_y = sqrt(pow((double)(j + 0.5), 2.0)) * dx;
        double d_z = sqrt(pow((double)(k + 0.5), 2.0)) * dx;
        // Octant 1
        green_grid[i][j][k] = 1.0 / (4.0 * pi *
                                     (sqrt(d_x * d_x + d_y * d_y + d_z * d_z +
                                           pow(gal->softening, 2))));
        // Octant 2
        green_grid[ngrid_padded[0] - 1 - i][j][k] = green_grid[i][j][k];
        // Octant 3
        green_grid[ngrid_padded[0] - 1 - i][ngrid_padded[1] - 1 - j][k] =
            green_grid[i][j][k];
        // Octant 4
        green_grid[i][ngrid_padded[1] - 1 - j][k] = green_grid[i][j][k];
        // Octant 5
        green_grid[i][j][ngrid_padded[2] - 1 - k] = green_grid[i][j][k];
        // Octant 6
        green_grid[ngrid_padded[0] - 1 - i][j][ngrid_padded[2] - 1 - k] =
            green_grid[i][j][k];
        // Octant 7
        green_grid[ngrid_padded[0] - 1 - i][ngrid_padded[1] - 1 - j]
                  [ngrid_padded[2] - 1 - k] = green_grid[i][j][k];
        // Octant 8
        green_grid[i][ngrid_padded[1] - 1 - j][ngrid_padded[2] - 1 - k] =
            green_grid[i][j][k];
      }
    }
  }

  // Pack Green's function and the density into 1D arrays
  l = 0;
  for (i = 0; i < ngrid_padded[0]; ++i) {
    for (j = 0; j < ngrid_padded[1]; ++j) {
      for (k = 0; k < ngrid_padded[2]; ++k) {
        green[l] = green_grid[i][j][k];
        rho[l] = potential[i][j][k];
        l++;
      }
    }
  }

  // Perform the fourier transforms. Density first, Green's function second.
  fftw_execute(fft_rho);
  fftw_execute(fft_green);
  // FFT is computed, we can free the memory
  fftw_destroy_plan(fft_rho);
  fftw_destroy_plan(fft_green);
  // Allocating memory for the inverse fourier computation
  fftinv_potential =
      fftw_plan_dft_3d(ngrid_padded[0], ngrid_padded[1], ngrid_padded[2], rho,
                       rho, FFTW_BACKWARD, FFTW_ESTIMATE);
  // Multiply the density by Green's function to find the k-space potential and
  // invert for the real potenital. Second, normalize the system and, finally,
  // put the potential information into the grid.
  for (i = 0; i < n; ++i) {
    // Convolve the potential
    rho[i] = green[i] * rho[i];
  }
  // Inversion
  fftw_execute(fftinv_potential);
  fftw_destroy_plan(fftinv_potential);
  // Normalization
  for (i = 0; i < n; ++i) {
    rho[i] = rho[i] * pow(dx, 3) / n;
  }
  l = 0;
  for (i = 0; i < ngrid_padded[0]; ++i) {
    for (j = 0; j < ngrid_padded[1]; ++j) {
      for (k = 0; k < ngrid_padded[2]; ++k) {
        // Fix the grid info
        potential[i][j][k] = -4.0 * pi * rho[l] * G;
        l++;
      }
    }
  }
  // Free fftw plan.
  // Kill the storage arrays since they are no longer needed.
  fftw_free(green);
  green = NULL;
  fftw_free(rho);
  rho = NULL;
  for (i = 0; i < ngrid_padded[1]; ++i) {
    for (j = 0; j < ngrid_padded[2]; ++j) {
      free(green_grid[i][j]);
      green_grid[i][j] = NULL;
    }
    free(green_grid[i]);
    green_grid[i] = NULL;
  }
  free(green_grid);
  green_grid = NULL;

#if USE_THREADS == 1
  fftw_cleanup_threads();
#endif
  // Flag the galaxy structure
  gal->potential_defined = 1;
  return 0;
}

// This function calculates the potential due to the disk at a point x,y,z
// by interpolating between grid points on the particle mesh. The
// interpolation routine uses the CIC kernel, which oddly enough is just
// a bilinear interpolation scheme...
//
// If the point lies off of the particle mesh, it approximates the potential
// as a function of 1/r.
double galaxy_potential_func(galaxy *gal, double ***potential, double dx,
                             int *ngrid, double x, double y, double z,
                             int extrapol) {

  int node_x, node_y, node_z, offset;
  double pot1, pot2, pot3, pot4, pot5, pot6, pot7, pot8;
  int ngrid_padded[3];
  double a, r_p, r_max, d_x, d_y, d_z, t_x, t_y, t_z, pot;
  double xtemp, ytemp, ztemp, theta, phi, xmax, ymax, zmax, rnorm;

  if (gal->potential_defined == 0)
    return 0;

  ngrid_padded[0] = 2 * ngrid[0];
  ngrid_padded[1] = 2 * ngrid[1];
  ngrid_padded[2] = 2 * ngrid[2];

  // Scale the coordinates
  r_p = sqrt(x * x + y * y + z * z);
  xtemp = x / dx + ((double)(ngrid_padded[0] / 2) - 0.5 - 0.5);
  ytemp = y / dx + ((double)(ngrid_padded[1] / 2) - 0.5 - 0.5);
  ztemp = z / dx + ((double)(ngrid_padded[2] / 2) - 0.5 - 0.5);
  offset = 0;
  // Determine the parent node.
  node_x = floor(xtemp);
  node_y = floor(ytemp);
  node_z = floor(ztemp);
  r_max = dx * ((double)(ngrid_padded[0] / 4) - 0.5 - 0.5 - offset);

  // Consider points off the grid or continue
  // The real information lies in the original grid size
  // between 0 and Ng/2 for all the dimensions!
  if (r_p > r_max && extrapol == 1) {

    theta = atan2(y, x);
    phi = acos(z / r_p);
    zmax = cos(phi) * r_max;
    xmax = r_max * sin(phi) * cos(theta);
    ymax = r_max * sin(phi) * sin(theta);
    node_x = round(xmax / dx + ((double)(ngrid_padded[0] / 2) - 0.5 - 0.5));
    node_y = round(ymax / dx + ((double)(ngrid_padded[1] / 2) - 0.5 - 0.5));
    node_z = round(zmax / dx + ((double)(ngrid_padded[2] / 2) - 0.5 - 0.5));

    rnorm = sqrt(
        pow((double)(node_x - (ngrid_padded[0] / 2 - 0.5 - 0.5)) * dx, 2.0) +
        pow((double)(node_y - (ngrid_padded[1] / 2 - 0.5 - 0.5)) * dx, 2.0) +
        pow((double)(node_z - (ngrid_padded[2] / 2 - 0.5 - 0.5)) * dx, 2.0));

    pot = potential[node_x][node_y][node_z] * rnorm / r_p;
  } else {
    if ((node_x < 0) || (node_x >= ngrid_padded[0] - 1) || (node_y < 0) ||
        (node_y >= ngrid_padded[1] - 1) || (node_z < 0) ||
        (node_z >= ngrid_padded[2] - 1)) {
      // if(r_p>r_max){
      pot = 0.;
    } else {
      // Check to see if (x,y,z) is a grid point.
      if (xtemp == (double)node_y && ytemp == (double)node_y &&
          ztemp == (double)node_z) {
        // If (x,y,z) is a grid point, return its potential.
        pot = potential[node_x][node_y][node_z];
      } else {
        // If (x,y,z) is not a grid point, use the CIC
        // interpolation function to calculate the potential.
        pot1 = potential[node_x][node_y][node_z];
        pot2 = potential[node_x + 1][node_y][node_z];
        pot3 = potential[node_x][node_y + 1][node_z];
        pot4 = potential[node_x][node_y][node_z + 1];
        pot5 = potential[node_x][node_y + 1][node_z + 1];
        pot6 = potential[node_x + 1][node_y + 1][node_z];
        pot7 = potential[node_x + 1][node_y][node_z + 1];
        pot8 = potential[node_x + 1][node_y + 1][node_z + 1];
        // CIC fractions
        d_x = 1.0 - (xtemp - (double)node_x);
        d_y = 1.0 - (ytemp - (double)node_y);
        d_z = 1.0 - (ztemp - (double)node_z);
        t_x = 1.0 - d_x;
        t_y = 1.0 - d_y;
        t_z = 1.0 - d_z;
        // Return the interpolated potential.
        pot = d_x * d_y * d_z * pot1 + t_x * d_y * d_z * pot2 +
              d_x * t_y * d_z * pot3 + d_x * d_y * t_z * pot4 +
              d_x * t_y * t_z * pot5 + t_x * t_y * d_z * pot6 +
              t_x * d_y * t_z * pot7 + t_x * t_y * t_z * pot8;
      }
    }
  }
  return pot;
}

// A wrapper for the galaxy potential function using the cylindrical radius.
double galaxyr_potential_wrapper_func(double radius, void *params) {

  double x, y, z, pot;
  int tid;

#if USE_THREADS == 1
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  galaxy *gal = (galaxy *)params;

  x = radius * cos(gal->theta_cyl[gal->index[tid]]);
  y = radius * sin(gal->theta_cyl[gal->index[tid]]);
  z = gal->z[gal->index[tid]];

  pot = galaxy_total_potential(gal, x, y, z, 0, 0);

  return pot;
}

// A wrapper for the galaxy potential function using the spherical radius.
double galaxyrsph_potential_wrapper_func(double r_sph, void *params) {

  double x, y, z, r_cyl, pot;
  int tid;

#if USE_THREADS == 1
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  galaxy *gal = (galaxy *)params;

  z = cos(gal->phi_sph[gal->index[tid]]) * r_sph;
  r_cyl = sqrt(r_sph * r_sph - z * z);
  x = r_cyl * cos(gal->theta_cyl[gal->index[tid]]);
  y = r_cyl * sin(gal->theta_cyl[gal->index[tid]]);

  pot = galaxy_total_potential(gal, x, y, z, 0, 0);

  return pot;
}

// A wrapper for the galaxy potential function using the cartesian coordinate z.
double galaxyz_potential_wrapper_func(double z, void *params) {

  int tid;
  double pot, x, y;

#if USE_THREADS == 1
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  galaxy *gal = (galaxy *)params;

  if (gal->pseudo[tid]) {
    pot = galaxy_total_potential(gal, gal->x[gal->index[tid]],
                                 gal->y[gal->index[tid]], z, 0, 0);
  } else {
    pot = galaxy_total_potential(gal, gal->x[gal->index[tid]],
                                 gal->y[gal->index[tid]], z, 1, 0);
  }

  return pot;
}

// This function obtain the potential from the coarse and zoomed grids and
// combine it together
double galaxy_total_potential(galaxy *gal, double x, double y, double z,
                              int circular, int coarse) {
  int i, j;
  double pot, pot_fine, pot_ext, r_sph;
  double sigma, transition_factor1, transition_factor2;
  double transition_factor1x, transition_factor1y, transition_factor1z;

  r_sph = sqrt(x * x + y * y + z * z);
  // Coarse  potential
  pot = galaxy_potential_func(gal, gal->potential[0], gal->dx[0], gal->ngrid[0],
                              x, y, z, 1);
  // Finer potential grids
  if (coarse != 1) {
    for (i = 1; i < gal->nlevel; i++) {
      sigma = 1.0 * gal->dx[i];
      if (circular == 1) {
        transition_factor1 =
            0.5 *
            (1 + erf((r_sph - (0.45 * gal->boxsize[i])) / (sigma * sqrt(2))));
      } else {
        transition_factor1x =
            0.5 * (1 + erf((sqrt(x * x) - (0.5 * gal->boxsize[i])) /
                           (sigma * sqrt(2))));
        transition_factor1y =
            0.5 * (1 + erf((sqrt(y * y) - (0.5 * gal->boxsize[i])) /
                           (sigma * sqrt(2))));
        transition_factor1z =
            0.5 * (1 + erf((sqrt(z * z) - (0.5 * gal->boxsize[i])) /
                           (sigma * sqrt(2))));
        transition_factor1 = max(transition_factor1x,
                                 max(transition_factor1y, transition_factor1z));
      }
      transition_factor2 = 1 - transition_factor1;

      pot_ext = 0.;
      for (j = 0; j < i; j++) {
        pot_ext += galaxy_potential_func(gal, gal->potential_ext[j], gal->dx[j],
                                         gal->ngrid[j], x, y, z, 0);
      }
      pot_fine = galaxy_potential_func(gal, gal->potential[i], gal->dx[i],
                                       gal->ngrid[i], x, y, z, 0);

      pot =
          transition_factor1 * pot + (pot_ext + pot_fine) * transition_factor2;
    }
  }
  return pot;
}

double get_h_value(galaxy *gal, double x, double y, double z, int circular,
                   int coarse) {
  int i;
  double h, r_sph, scale;
  double sigma, transition_factor1, transition_factor2;
  double transition_factor1x, transition_factor1y, transition_factor1z;

  scale = 1.0;
  r_sph = sqrt(x * x + y * y + z * z);
  h = scale * gal->dx[0];
  // Check finer levels
  for (i = 1; i < gal->nlevel; i++) {
    sigma = 2.0 * gal->dx[i];
    if (circular == 1) {
      transition_factor1 =
          0.5 *
          (1 + erf((r_sph - (0.5 * gal->boxsize[i])) / (sigma * sqrt(2))));
    } else {
      transition_factor1x =
          0.5 * (1 + erf((sqrt(x * x) - (0.5 * gal->boxsize[i])) /
                         (sigma * sqrt(2))));
      transition_factor1y =
          0.5 * (1 + erf((sqrt(y * y) - (0.5 * gal->boxsize[i])) /
                         (sigma * sqrt(2))));
      transition_factor1z =
          0.5 * (1 + erf((sqrt(z * z) - (0.5 * gal->boxsize[i])) /
                         (sigma * sqrt(2))));
      transition_factor1 = max(transition_factor1x,
                               max(transition_factor1y, transition_factor1z));
    }
    transition_factor2 = 1 - transition_factor1;

    h = transition_factor1 * h + scale * gal->dx[i] * transition_factor2;
  }
  return h;
}

// A wrapper for the second order derivative of the potential function.
double potential_deriv_wrapper_func(double radius, void *params) {

  galaxy *gal = (galaxy *)params;

  if (radius == 0.) {
    return 0.;
  } else {
    return (pow(v_c_func(gal, fabs(radius)), 2.0)) / (fabs(radius));
  }
}

// This function copies one galaxy to another.
void copy_potential(galaxy *gal_1, galaxy *gal_2, int info) {

  int n, i, j, k;

  if (info == 1)
    printf("/////\tCopying potential grid \n");
  // Copy all the coordinate information.
  for (n = 0; n < gal_1->nlevel; ++n) {
    for (i = 0; i < gal_1->ngrid[n][0] * 2; ++i) {
      for (j = 0; j < gal_1->ngrid[n][2] * 2; ++j) {
        for (k = 0; k < gal_1->ngrid[n][2] * 2; ++k) {
          gal_2->potential[n][i][j][k] = gal_1->potential[n][i][j][k];
        }
      }
    }
  }
  gal_2->potential_defined = 1;
  if (info == 1)
    printf("/////\tPotential grid copied \n");
  return;
}

// This function sets the gravitational potential on all the levels
double set_galaxy_potential_all(galaxy *gal, int verbose) {
  int n;
  double exclude[3];

  for (n = 0; n < gal->nlevel; n++) {
    if (verbose == 1) {
      printf("> %d ", n + gal->level_coarse);
    }
    fflush(stdout);
    exclude[0] = 0.;
    exclude[1] = 0.;
    exclude[2] = 0.;
    if (set_galaxy_potential(gal, gal->potential[n], gal->dx[n], gal->ngrid[n],
                             0, exclude) != 0) {
      fprintf(stderr, "\n[Error] Unable to set the potential for level %d\n",
              n + gal->level_coarse);
      return -1;
    }
    if (n > 0) {
      exclude[0] = gal->boxsize[n] * gal->boxsize_flatx[n] / 2.;
      exclude[1] = gal->boxsize[n] * gal->boxsize_flaty[n] / 2.;
      exclude[2] = gal->boxsize[n] * gal->boxsize_flatz[n] / 2.;
      if (set_galaxy_potential(gal, gal->potential_ext[n - 1], gal->dx[n - 1],
                               gal->ngrid[n - 1], 0, exclude) != 0) {
        fprintf(stderr,
                "\n[Error] Unable to set the external potential for level %d\n",
                n + gal->level_coarse);
        return -1;
      }
    }
  }
  return 0;
}
