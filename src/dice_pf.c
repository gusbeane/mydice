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
/ Date: October 2014
/
*///---------------------------------------------------------------------------

#include "dice.h"

// This function calculates the potential due to a galactic disk using Cloud-
// In-Cell mass assignment with vacuum (isolated) boundary conditions on a
// Cartesian grid. The method was adapted from the discussion found in Hockney
// and Eastwood, "Computer Simulation Using Particles," 1981.
int set_galaxy_potential(galaxy *gal, int verbose) {
	unsigned long int ii, start, end, l;
	// Loop variables
	int i, j, k;
	// Nodes coordinates
	int node_x, node_y, node_z;
	double space_x, space_y, space_z;
	double dx, dy, dz, tx, ty, tz, n;
	double x, y, z;
	// The particle-mesh grid size and the Green's function and potential
	// storage buffers. Global to keep from hitting the stack limit for large grid
	// size.
	double ***green_grid;
	//size_t p_x[gal->num_part[1]], p_y[gal->num_part[1]], p_z[gal->num_part[1]];
	fftw_plan fft_green, fft_potential, fftinv_potential;
	fftw_complex *green,*potential;
	
	if (verbose) printf("/////\tComputing potential\n");
	fflush(stdout);
	// Setup fftw threads
	#if USE_THREADS == 1
		if (verbose) printf("/////\t\t-Setting up FFT call with %d threads\n",AllVars.Nthreads);
		fflush(stdout);
		fftw_init_threads();
		fftw_plan_with_nthreads(AllVars.Nthreads);
		fflush(stdout);
	#endif
	
	// Allocate grid storage variables
	if (!(green = calloc(pow(gal->ngrid_padded,3),sizeof(fftw_complex)))) {
		fprintf(stderr,"Unable to allocate space for Green's function.\n");
		return -1;
	}
	if (!(potential = calloc(pow(gal->ngrid_padded,3),sizeof(fftw_complex)))) {
		fprintf(stderr,"Unable to allocate space for potential buffer.\n");
		return -1;
	}
	
	if (!(green_grid=calloc(gal->ngrid_padded,sizeof(double *)))) {
		fprintf(stderr,"Unable to create Green's function x axis.\n");
		return -1;
	}
	// Green grid allocation
	// x-axis
	for (i = 0; i < gal->ngrid_padded; ++i) {
	        // y-axis
	        if (!(green_grid[i] = calloc(gal->ngrid_padded,sizeof(double *)))) {
			fprintf(stderr,"Unable to create Green's function y axis.\n");
           		return -1;
        	}
	        // z-axis
        	for (j = 0; j < gal->ngrid_padded; ++j) {
			if (!(green_grid[i][j] = calloc(gal->ngrid_padded,sizeof(double)))) {
				fprintf(stderr,"Unable to create Green's function z axis.\n");
				return -1;
			}
		}
	}
	// Define which particles we should use in the potential computation
	// We use all the particles
	start 	= 0;
	end 	= gal->ntot_part_pot;
	/*
	if (!(x = calloc(end-start,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate space for x array buffer.\n");
		return -1;
	}
	if (!(y = calloc(end-start,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate space for y array buffer.\n");
		return -1;
	}
	if (!(z = calloc(end-start,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate space for z array buffer.\n");
		return -1;
	}
	*/
	// Allocate the fftw complex output value and the fftw dft plan.
	fft_potential	= fftw_plan_dft_3d(gal->ngrid_padded,gal->ngrid_padded,gal->ngrid_padded,potential,potential,FFTW_FORWARD,FFTW_ESTIMATE);
	fft_green		= fftw_plan_dft_3d(gal->ngrid_padded,gal->ngrid_padded,gal->ngrid_padded,green,green,FFTW_FORWARD,FFTW_ESTIMATE);
	
	// Normalization constant
	// See FFTW reference guide for more details
	n = (int)pow(gal->ngrid_padded,3.0);
	
	// Check for bad grids
	if (gal->ngrid_padded <= 0) {
		fprintf(stderr,"\t\tGrid dimensions must be greater than zero! (ngrid=%d)\n",gal->ngrid_padded);
		return -1;
	}
	
	// Sort the position arrays and figure out the spacing between grid points.
	// Subtract 2 for a.) the C offset and b.) the CIC offset. Finally, store
	// the values in the galaxy for later use.
	//gsl_sort_index(p_x,gal->x,1,gal->num_part[0]);
	//gsl_sort_index(p_y,gal->y,1,gal->num_part[0]);
	//gsl_sort_index(p_z,gal->z,1,gal->num_part[0]);
	space_x = gal->space[0];
	space_y = gal->space[1];
	space_z = gal->space[2];
	// Print the potential in the xy-plane for z = 0 if the option is set.
	if (verbose) printf("/////\t\t-Grid cell spacings [kpc]: dx = %.3f dy = %.3f dz = %.3f\n",space_x,space_y,space_z);
	fflush(stdout);
	// Make the coordinate information unitless.
	for (ii = 0; ii < end-start; ++ii) {
	    x = gal->x[ii+start]/(space_x)+((double)(gal->ngrid_padded/2)-0.5);
	    y = gal->y[ii+start]/(space_y)+((double)(gal->ngrid_padded/2)-0.5);
	    z = gal->z[ii+start]/(space_z)+((double)(gal->ngrid_padded/2)-0.5);
		//if(x > 3*gal->ngrid_padded/4 || x < gal->ngrid_padded/4 || y > 3*gal->ngrid_padded/4 || y < gal->ngrid_padded/4 || z > 3*gal->ngrid_padded/4 || z < gal->ngrid_padded/4) {
		//	fprintf(stderr,"\t\tGrid cell spacing is too low. Unable to cover the entire particles ditribution.\n\t\t\tIncrease boxsize in the galaxy parameter file.\n");
		//	return -1;
		//}
	}
	// Initialization loop
	for (i = 0; i < gal->ngrid_padded; ++i) {
		for (j = 0; j < gal->ngrid_padded; ++j) {
			for (k = 0; k < gal->ngrid_padded; ++k) {
				gal->potential[i][j][k] = 0.;
			}
		}
	}
	
	// Calculate the density using the CIC routine. The positions are shifted
	// such that the particles are in the +x,+y,+z octant. space_* is added
	// to take care of the vacuum boundary conditions. The density values are
	// stored in the potential structure for now.
	#pragma omp parallel for private(x, y, z, node_x, node_y, node_z, dx, dy, dz, tx, ty, tz) shared(gal)
	for (ii = 0; ii < end-start; ++ii) {
		x = gal->x[ii+start]/(space_x)+((double)(gal->ngrid_padded/2)-0.5);
		y = gal->y[ii+start]/(space_y)+((double)(gal->ngrid_padded/2)-0.5);
		z = gal->z[ii+start]/(space_z)+((double)(gal->ngrid_padded/2)-0.5);
		if(x <= 3*gal->ngrid_padded/4 || x >= gal->ngrid_padded/4 || y <= 3*gal->ngrid_padded/4 || y >= gal->ngrid_padded/4 || z <= 3*gal->ngrid_padded/4 || z >= gal->ngrid_padded/4) {
			// Figure out which node owns the particle
			node_x = (int) x;
			node_y = (int) y;
			node_z = (int) z;
			// Check if particle is not outside the potential grid
			if(node_x>=0 && node_x<gal->ngrid_padded-1 && node_y>=0 && node_y<gal->ngrid_padded-1 && node_z>=0 && node_z<gal->ngrid_padded-1) {
				// Set the CIC size fractions
				dx = 1.0 - (x - (double) node_x);
				dy = 1.0 - (y - (double) node_y);
				dz = 1.0 - (z - (double) node_z);
				tx = 1.0 - dx;
				ty = 1.0 - dy;
				tz = 1.0 - dz;
				// Calculate the CIC densities
				gal->potential[node_x][node_y][node_z]		+= gal->mass[ii+start]*(dx*dy*dz) / (space_x*space_y*space_z);
				gal->potential[node_x+1][node_y][node_z]	+= gal->mass[ii+start]*(tx*dy*dz) / (space_x*space_y*space_z);
				gal->potential[node_x][node_y+1][node_z]	+= gal->mass[ii+start]*(dx*ty*dz) / (space_x*space_y*space_z);
				gal->potential[node_x][node_y][node_z+1]	+= gal->mass[ii+start]*(dx*dy*tz) / (space_x*space_y*space_z);
				gal->potential[node_x][node_y+1][node_z+1]	+= gal->mass[ii+start]*(dx*ty*tz) / (space_x*space_y*space_z);
				gal->potential[node_x+1][node_y+1][node_z]	+= gal->mass[ii+start]*(tx*ty*dz) / (space_x*space_y*space_z);
				gal->potential[node_x+1][node_y][node_z+1]	+= gal->mass[ii+start]*(tx*dy*tz) / (space_x*space_y*space_z);
				gal->potential[node_x+1][node_y+1][node_z+1]	+= gal->mass[ii+start]*(tx*ty*tz) / (space_x*space_y*space_z);
				// NGP (Nearest Point Grid) if anyone needs it
				//gal->potential[node_x][node_y][node_z] += gal->mass[i]/(space_x*space_y*space_z);
			}
		}
	}
	// Define Green's function.
	// These are the grid points as measured from the center of Green's function
	// and the local value of the truncation function. The density is also packed into a
	// buffer here.
	//
	// Green's function is defined eight times here, once per octant, to take care of the
	// isolated (vacuum) boundary conditions. See Hockney and Eastwood, 1980, ch. 6 for a
	// discussion. The octants start in the lower left at (p,q) = (0,0) and progress
	// counter clockwise.
	for (i = 0; i < gal->ngrid_padded/2; ++i) {
	        for (j = 0; j < gal->ngrid_padded/2; ++j) {
	        	#pragma omp parallel for private(dx, dy, dz) shared(green_grid,i,j)
				for (k = 0; k < gal->ngrid_padded/2; ++k) {
		    		dx = sqrt(pow((double)(i+0.5),2.0))*space_x;
		        	dy = sqrt(pow((double)(j+0.5),2.0))*space_y;
		            dz = sqrt(pow((double)(k+0.5),2.0))*space_z;
					// Octant 1
		            green_grid[i][j][k] = 1.0 / (4.0*pi*sqrt(dx*dx + dy*dy + dz*dz));
		             // Octant 2
                	green_grid[gal->ngrid_padded-1-i][j][k] = green_grid[i][j][k];
                	// Octant 3
                	green_grid[gal->ngrid_padded-1-i][gal->ngrid_padded-1-j][k] = green_grid[i][j][k];
                	// Octant 4
                	green_grid[i][gal->ngrid_padded-1-j][k] = green_grid[i][j][k];
                	// Octant 5
                	green_grid[i][j][gal->ngrid_padded-1-k] = green_grid[i][j][k];
                	// Octant 6
                	green_grid[gal->ngrid_padded-1-i][j][gal->ngrid_padded-1-k] = green_grid[i][j][k];
                	// Octant 7
                	green_grid[gal->ngrid_padded-1-i][gal->ngrid_padded-1-j][gal->ngrid_padded-1-k] = green_grid[i][j][k];
                	// Octant 8
                	green_grid[i][gal->ngrid_padded-1-j][gal->ngrid_padded-1-k] = green_grid[i][j][k];
			}
		}
	}
	
	// Pack Green's function and the density into 1D arrays
	l = 0;
	for (i = 0; i < gal->ngrid_padded; ++i) {
	        for (j = 0; j < gal->ngrid_padded; ++j) {
				for (k = 0; k < gal->ngrid_padded; ++k) {
		    		green[l] = green_grid[i][j][k];
		        	potential[l] = gal->potential[i][j][k];
					l++;
				}
		}
	}
    	
	// Perform the fourier transforms. Density first, Green's function second.
	fftw_execute(fft_potential);
	fftw_execute(fft_green);
	// FFT is computed, we can free the memory
	fftw_destroy_plan(fft_potential);
	fftw_destroy_plan(fft_green);
	// Allocating memory for the inverse fourier computation
	fftinv_potential = fftw_plan_dft_3d(gal->ngrid_padded,gal->ngrid_padded,gal->ngrid_padded,potential,potential,FFTW_BACKWARD,FFTW_ESTIMATE);
	// Multiply the density by Green's function to find the k-space potential and
	// invert for the real potenital. Second, normalize the system and, finally,
	// put the potential information into the grid.
	for (i = 0; i < n; ++i) {
	        // Convolve the potential
	        potential[i] = green[i]*potential[i];
	}
	// Inversion
	fftw_execute(fftinv_potential);
	fftw_destroy_plan(fftinv_potential);
	// Normalization
	for (i = 0; i < n; ++i) {
	        potential[i] = potential[i]*space_x*space_y*space_z/n;
	}
	l = 0;
	for (i = 0; i < gal->ngrid_padded; ++i) {
	    for (j = 0; j < gal->ngrid_padded; ++j) {
			for (k = 0; k < gal->ngrid_padded; ++k) {
				// Fix the grid info
				gal->potential[i][j][k] = -4.0*pi*potential[l]*G*unit_mass/kpc;
				l++;
			}
		}
	}
	// Free fftw plan.
	// Kill the storage arrays since they are no longer needed.
	fftw_free(green);
	fftw_free(potential);
	for (i = 0; i < gal->ngrid_padded; ++i) {
	        for (j = 0; j < gal->ngrid_padded; ++j) {
			free(green_grid[i][j]);
		}
		free(green_grid[i]);
	}
	free(green_grid);
	
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
double galaxy_potential_func(galaxy *gal, double x, double y, double z) {
	
	int node_x, node_y, node_z, offset;
	double pot1, pot2, pot3, pot4, pot5, pot6, pot7, pot8;
	double a, r_p, r_max, dx, dy, dz, tx, ty, tz, potential;
	double xtemp,ytemp,ztemp,theta,phi,xmax,ymax,zmax,rnorm;
	
	if(gal->potential_defined == 0) return 0;
	
	// Scale the coordinates
	r_p = sqrt(x*x + y*y + z*z);
	xtemp = x/gal->space[0] + ((double)(gal->ngrid_padded/2)-0.5-0.5);
	ytemp = y/gal->space[1] + ((double)(gal->ngrid_padded/2)-0.5-0.5);
	ztemp = z/gal->space[2] + ((double)(gal->ngrid_padded/2)-0.5-0.5);
	offset = 0;
	// Determine the parent node.
	node_x = floor(xtemp);
	node_y = floor(ytemp);
	node_z = floor(ztemp);
	r_max = gal->space[0]*((double)(gal->ngrid_padded/4)-0.5-0.5-offset);
	
	// Consider points off the grid or continue
	// The real information lies in the original grid size
	// between 0 and Ng/2 for all the dimensions!
	if (r_p > r_max) {
		
		theta = atan2(y,x);
		phi = acos(z/r_p);
		zmax = cos(phi)*r_max;
		xmax = r_max*sin(phi)*cos(theta);
		ymax = r_max*sin(phi)*sin(theta);
		node_x = round(xmax/gal->space[0] + ((double)(gal->ngrid_padded/2)-0.5-0.5));
		node_y = round(ymax/gal->space[1] + ((double)(gal->ngrid_padded/2)-0.5-0.5));
		node_z = round(zmax/gal->space[2] + ((double)(gal->ngrid_padded/2)-0.5-0.5));
		
		rnorm = sqrt( pow((double)(node_x-(gal->ngrid_padded/2-0.5-0.5))*gal->space[0],2.0)
			     +pow((double)(node_y-(gal->ngrid_padded/2-0.5-0.5))*gal->space[1],2.0)
			     +pow((double)(node_z-(gal->ngrid_padded/2-0.5-0.5))*gal->space[2],2.0));
		
		potential = gal->potential[node_x][node_y][node_z]*rnorm/r_p;
	} else {
		// Check to see if (x,y,z) is a grid point.
		if (xtemp == (double) node_y && ytemp == (double) node_y && ztemp == (double) node_z) {
			// If (x,y,z) is a grid point, return its potential.
			potential = gal->potential[node_x][node_y][node_z];
		} else {
			// If (x,y,z) is not a grid point, use the CIC
			// interpolation function to calculate the potential.
			pot1 = gal->potential[node_x][node_y][node_z];
			pot2 = gal->potential[node_x+1][node_y][node_z];
			pot3 = gal->potential[node_x][node_y+1][node_z];
			pot4 = gal->potential[node_x][node_y][node_z+1];
			pot5 = gal->potential[node_x][node_y+1][node_z+1];
			pot6 = gal->potential[node_x+1][node_y+1][node_z];
			pot7 = gal->potential[node_x+1][node_y][node_z+1];
			pot8 = gal->potential[node_x+1][node_y+1][node_z+1];
			// CIC fractions
			dx = 1.0 - (xtemp - (double) node_x);
			dy = 1.0 - (ytemp - (double) node_y);
			dz = 1.0 - (ztemp - (double) node_z);
			tx = 1.0 - dx;
			ty = 1.0 - dy;
			tz = 1.0 - dz;
			// Return the interpolated potential.
			potential = dx*dy*dz*pot1 + tx*dy*dz*pot2 +
			dx*ty*dz*pot3 + dx*dy*tz*pot4 +
			dx*ty*tz*pot5 + tx*ty*dz*pot6 +
			tx*dy*tz*pot7 + tx*ty*tz*pot8;
		}
	}
	return potential;
}

// A wrapper for the disk potential function.
double galaxyr_potential_wrapper_func(double radius, void *params) {
	
	double x, y;
	int tid;

	#if USE_THREADS == 1
		tid = omp_get_thread_num();
	#else
		tid = 0;
	#endif

	galaxy *gal = (galaxy *) params;
	
	x = radius*cos(gal->theta_cyl[gal->index[tid]]);
	y = radius*sin(gal->theta_cyl[gal->index[tid]]);
	
	return galaxy_potential_func(gal,x,y,gal->z[gal->index[tid]]);
}

// A wrapper for the disk potential function.
double galaxyz_potential_wrapper_func(double z, void *params) {
	
	int tid;
	
	#if USE_THREADS == 1
		tid = omp_get_thread_num();
	#else
		tid = 0;
	#endif

	galaxy *gal = (galaxy *) params;
	
	return galaxy_potential_func(gal,gal->x[gal->index[tid]],gal->y[gal->index[tid]],z);
}

// A wrapper for the second derivative of the potential function.
double potential_deriv_wrapper_func(double radius, void *params) {
	
	galaxy *gal = (galaxy *) params;
		
	return (pow(v_c_func(gal,radius),2.0))/(radius*kpc);
}
