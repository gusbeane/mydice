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
int set_galaxy_potential(galaxy *gal, double ***potential, double dx, int ngrid[3], int verbose) {
	unsigned long int ii, l;
	// Loop variables
	int i, j, k;
	// Nodes coordinates
	int node_x, node_y, node_z;
	int ngrid_padded[3];
	double d_x, d_y, d_z, t_x, t_y, t_z, n;
	double x, y, z;
	// The particle-mesh grid size and the Green's function and potential
	// storage buffers. Global to keep from hitting the stack limit for large grid
	// size.
	double ***green_grid;
	//size_t p_x[gal->num_part[1]], p_y[gal->num_part[1]], p_z[gal->num_part[1]];
	fftw_plan fft_green, fft_rho, fftinv_potential;
	fftw_complex *green,*rho;
	
	if (verbose) printf("/////\tComputing potential [%d FFT threads][dx=%.1lf pc][box=%.1lf kpc]\n",AllVars.Nthreads,dx*1e3,dx*ngrid[0]);
	fflush(stdout);
	// Setup fftw threads
	#if USE_THREADS == 1
		fflush(stdout);
		fftw_init_threads();
		fftw_plan_with_nthreads(AllVars.Nthreads);
		fflush(stdout);
	#endif
	
	ngrid_padded[0] = 2*ngrid[0];
	ngrid_padded[1] = 2*ngrid[1];
	ngrid_padded[2] = 2*ngrid[2];
	
	// Allocate grid storage variables
	if (!(green = calloc(ngrid_padded[0]*ngrid_padded[1]*ngrid_padded[2],sizeof(fftw_complex)))) {
		fprintf(stderr,"[Error] Unable to allocate space for Green's function\n");
		return -1;
	}
	if (!(rho = calloc(ngrid_padded[0]*ngrid_padded[1]*ngrid_padded[2],sizeof(fftw_complex)))) {
		fprintf(stderr,"[Error] Unable to allocate space for rho function\n");
		return -1;
	}
	
	if (!(green_grid=calloc(ngrid_padded[0],sizeof(double *)))) {
		fprintf(stderr,"[Error] Unable to create Green's function x axis\n");
		return -1;
	}
	// Green grid allocation
	// x-axis
	for (i = 0; i < ngrid_padded[1]; ++i) {
	        // y-axis
	        if (!(green_grid[i] = calloc(ngrid_padded[1],sizeof(double *)))) {
			fprintf(stderr,"[Error] Unable to create Green's function y axis\n");
           		return -1;
        	}
	        // z-axis
        	for (j = 0; j < ngrid_padded[2]; ++j) {
			if (!(green_grid[i][j] = calloc(ngrid_padded[2],sizeof(double)))) {
				fprintf(stderr,"[Error] Unable to create Green's function z axis\n");
				return -1;
			}
		}
	}
	// Allocate the fftw complex output value and the fftw dft plan.
	fft_rho		= fftw_plan_dft_3d(ngrid_padded[0],ngrid_padded[1],ngrid_padded[2],rho,rho,FFTW_FORWARD,FFTW_ESTIMATE);
	fft_green	= fftw_plan_dft_3d(ngrid_padded[0],ngrid_padded[1],ngrid_padded[2],green,green,FFTW_FORWARD,FFTW_ESTIMATE);
	
	// Normalization constant
	// See FFTW reference guide for more details
	n = (int)(ngrid_padded[0]*ngrid_padded[1]*ngrid_padded[2]);
	
	// Check for bad grids
	if (ngrid_padded[0] <= 0) {
		fprintf(stderr,"/////\t\t[Error] Grid dimensions must be greater than zero [ngrid=%d]\n",ngrid_padded[0]);
		return -1;
	}
	
	// Sort the position arrays and figure out the spacing between grid points.
	// Subtract 2 for a.) the C offset and b.) the CIC offset. Finally, store
	// the values in the galaxy for later use.
	//gsl_sort_index(p_x,gal->x,1,gal->num_part[0]);
	//gsl_sort_index(p_y,gal->y,1,gal->num_part[0]);
	//gsl_sort_index(p_z,gal->z,1,gal->num_part[0]);
	
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
	#pragma omp parallel for private(x, y, z, node_x, node_y, node_z, d_x, d_y, d_z, t_x, t_y, t_z, ii) shared(gal) 
	for (ii = 0; ii < gal->ntot_part_pot; ++ii) {
		x = gal->x[ii]/dx+((double)(ngrid_padded[0]/2)-0.5);
		y = gal->y[ii]/dx+((double)(ngrid_padded[1]/2)-0.5);
		z = gal->z[ii]/dx+((double)(ngrid_padded[2]/2)-0.5);
		// Figure out which node owns the particle
		node_x = (int) x;
		node_y = (int) y;
		node_z = (int) z;
		// Check if particle is not outside the potential grid
		if(node_x>=0 && node_x<ngrid_padded[0]-1 && node_y>=0 && node_y<ngrid_padded[1]-1 && node_z>=0 && node_z<ngrid_padded[2]-1) {
			// Set the CIC size fractions
			d_x = 1.0 - (x - (double) node_x);
			d_y = 1.0 - (y - (double) node_y);
			d_z = 1.0 - (z - (double) node_z);
			t_x = 1.0 - d_x;
			t_y = 1.0 - d_y;
			t_z = 1.0 - d_z;
			// Calculate the CIC densities
			potential[node_x][node_y][node_z]		+= gal->mass[ii]*(d_x*d_y*d_z) / pow(dx,3);
			potential[node_x+1][node_y][node_z]		+= gal->mass[ii]*(t_x*d_y*d_z) / pow(dx,3);
			potential[node_x][node_y+1][node_z]		+= gal->mass[ii]*(d_x*t_y*d_z) / pow(dx,3);
			potential[node_x][node_y][node_z+1]		+= gal->mass[ii]*(d_x*d_y*t_z) / pow(dx,3);
			potential[node_x][node_y+1][node_z+1]	+= gal->mass[ii]*(d_x*t_y*t_z) / pow(dx,3);
			potential[node_x+1][node_y+1][node_z]	+= gal->mass[ii]*(t_x*t_y*d_z) / pow(dx,3);
			potential[node_x+1][node_y][node_z+1]	+= gal->mass[ii]*(t_x*d_y*t_z) / pow(dx,3);
			potential[node_x+1][node_y+1][node_z+1]	+= gal->mass[ii]*(t_x*t_y*t_z) / pow(dx,3);
			// NGP (Nearest Point Grid)
			//potential[node_x][node_y][node_z] += gal->mass[i]/(space_x*space_y*space_z);
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
	for (i = 0; i < ngrid_padded[0]/2; ++i) {
	        for (j = 0; j < ngrid_padded[1]/2; ++j) {
	        	#pragma omp parallel for private(d_x, d_y, d_z, k) shared(green_grid,i,j)
				for (k = 0; k < ngrid_padded[2]/2; ++k) {
		    		d_x = sqrt(pow((double)(i+0.5),2.0))*dx;
		        	d_y = sqrt(pow((double)(j+0.5),2.0))*dx;
		            d_z = sqrt(pow((double)(k+0.5),2.0))*dx;
					// Octant 1
		            green_grid[i][j][k] 														= 1.0 / (4.0*pi*sqrt(d_x*d_x + d_y*d_y + d_z*d_z));
		             // Octant 2
                	green_grid[ngrid_padded[0]-1-i][j][k] 										= green_grid[i][j][k];
                	// Octant 3
                	green_grid[ngrid_padded[0]-1-i][ngrid_padded[1]-1-j][k] 					= green_grid[i][j][k];
                	// Octant 4
                	green_grid[i][ngrid_padded[1]-1-j][k] 										= green_grid[i][j][k];
                	// Octant 5
                	green_grid[i][j][ngrid_padded[2]-1-k] 										= green_grid[i][j][k];
                	// Octant 6
                	green_grid[ngrid_padded[0]-1-i][j][ngrid_padded[2]-1-k] 					= green_grid[i][j][k];
                	// Octant 7
                	green_grid[ngrid_padded[0]-1-i][ngrid_padded[1]-1-j][ngrid_padded[2]-1-k] 	= green_grid[i][j][k];
                	// Octant 8
                	green_grid[i][ngrid_padded[1]-1-j][ngrid_padded[2]-1-k] 					= green_grid[i][j][k];
			}
		}
	}
	
	// Pack Green's function and the density into 1D arrays
	l = 0;
	for (i = 0; i < ngrid_padded[0]; ++i) {
	        for (j = 0; j < ngrid_padded[1]; ++j) {
				for (k = 0; k < ngrid_padded[2]; ++k) {
		    		green[l] 	= green_grid[i][j][k];
		        	rho[l] 		= potential[i][j][k];
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
	fftinv_potential = fftw_plan_dft_3d(ngrid_padded[0],ngrid_padded[1],ngrid_padded[2],rho,rho,FFTW_BACKWARD,FFTW_ESTIMATE);
	// Multiply the density by Green's function to find the k-space potential and
	// invert for the real potenital. Second, normalize the system and, finally,
	// put the potential information into the grid.
	for (i = 0; i < n; ++i) {
	        // Convolve the potential
	        rho[i] = green[i]*rho[i];
	}
	// Inversion
	fftw_execute(fftinv_potential);
	fftw_destroy_plan(fftinv_potential);
	// Normalization
	for (i = 0; i < n; ++i) {
	        rho[i] = rho[i]*pow(dx,3)/n;
	}
	l = 0;
	for (i = 0; i < ngrid_padded[0]; ++i) {
	    for (j = 0; j < ngrid_padded[1]; ++j) {
			for (k = 0; k < ngrid_padded[2]; ++k) {
				// Fix the grid info
				potential[i][j][k] = -4.0*pi*rho[l]*G*unit_mass/kpc;
				l++;
			}
		}
	}
	// Free fftw plan.
	// Kill the storage arrays since they are no longer needed.
	fftw_free(green);
	fftw_free(rho);
	for (i = 0; i < ngrid_padded[1]; ++i) {
	        for (j = 0; j < ngrid_padded[2]; ++j) {
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
double galaxy_potential_func(galaxy *gal, double ***potential, double dx, int ngrid[3], double x, double y, double z, int interp) {
	
	int node_x, node_y, node_z, offset;
	double pot1, pot2, pot3, pot4, pot5, pot6, pot7, pot8;
	int ngrid_padded[3];
	double a, r_p, r_max, d_x, d_y, d_z, t_x, t_y, t_z, pot;
	double xtemp,ytemp,ztemp,theta,phi,xmax,ymax,zmax,rnorm;
	
	if(gal->potential_defined == 0) return 0;
	
	ngrid_padded[0] = 2*ngrid[0];
	ngrid_padded[1] = 2*ngrid[1];
	ngrid_padded[2] = 2*ngrid[2];
	
	// Scale the coordinates
	r_p = sqrt(x*x + y*y + z*z);
	xtemp = x/dx + ((double)(ngrid_padded[0]/2)-0.5-0.5);
	ytemp = y/dx + ((double)(ngrid_padded[1]/2)-0.5-0.5);
	ztemp = z/dx + ((double)(ngrid_padded[2]/2)-0.5-0.5);
	offset = 0;
	// Determine the parent node.
	node_x = floor(xtemp);
	node_y = floor(ytemp);
	node_z = floor(ztemp);
	r_max = dx*((double)(ngrid_padded[0]/4)-0.5-0.5-offset);
	
	// Consider points off the grid or continue
	// The real information lies in the original grid size
	// between 0 and Ng/2 for all the dimensions!
	if (r_p > r_max && interp==1) {
		
		theta = atan2(y,x);
		phi = acos(z/r_p);
		zmax = cos(phi)*r_max;
		xmax = r_max*sin(phi)*cos(theta);
		ymax = r_max*sin(phi)*sin(theta);
		node_x = round(xmax/dx + ((double)(ngrid_padded[0]/2)-0.5-0.5));
		node_y = round(ymax/dx + ((double)(ngrid_padded[1]/2)-0.5-0.5));
		node_z = round(zmax/dx + ((double)(ngrid_padded[2]/2)-0.5-0.5));
		
		rnorm = sqrt( pow((double)(node_x-(ngrid_padded[0]/2-0.5-0.5))*dx,2.0)
			         +pow((double)(node_y-(ngrid_padded[1]/2-0.5-0.5))*dx,2.0)
			         +pow((double)(node_z-(ngrid_padded[2]/2-0.5-0.5))*dx,2.0));
		
		pot = potential[node_x][node_y][node_z]*rnorm/r_p;
	} else {
		if((node_x<0)||(node_x>=ngrid_padded[0]-1)||(node_y<0)||(node_y>=ngrid_padded[1]-1)||(node_z<0)||(node_z>=ngrid_padded[2]-1)){
		//if(r_p>r_max){
			pot = 0.;
		} else {
			// Check to see if (x,y,z) is a grid point.
			if (xtemp == (double) node_y && ytemp == (double) node_y && ztemp == (double) node_z) {
				// If (x,y,z) is a grid point, return its potential.
				pot = potential[node_x][node_y][node_z];
			} else {
				// If (x,y,z) is not a grid point, use the CIC
				// interpolation function to calculate the potential.
				pot1 = potential[node_x][node_y][node_z];
				pot2 = potential[node_x+1][node_y][node_z];
				pot3 = potential[node_x][node_y+1][node_z];
				pot4 = potential[node_x][node_y][node_z+1];
				pot5 = potential[node_x][node_y+1][node_z+1];
				pot6 = potential[node_x+1][node_y+1][node_z];
				pot7 = potential[node_x+1][node_y][node_z+1];
				pot8 = potential[node_x+1][node_y+1][node_z+1];
				// CIC fractions
				d_x = 1.0 - (xtemp - (double) node_x);
				d_y = 1.0 - (ytemp - (double) node_y);
				d_z = 1.0 - (ztemp - (double) node_z);
				t_x = 1.0 - d_x;
				t_y = 1.0 - d_y;
				t_z = 1.0 - d_z;
				// Return the interpolated potential.
				pot = d_x*d_y*d_z*pot1 + t_x*d_y*d_z*pot2 +
				d_x*t_y*d_z*pot3 + d_x*d_y*t_z*pot4 +
				d_x*t_y*t_z*pot5 + t_x*t_y*d_z*pot6 +
				t_x*d_y*t_z*pot7 + t_x*t_y*t_z*pot8;
			}
		}
	}
	return pot;
}

// A wrapper for the galaxy potential function using the cylindrical radius.
double galaxyr_potential_wrapper_func(double radius, void *params) {
	
	double x, y, pot, r_sph;
	double sigma,transition_factor1,transition_factor2;
	int tid;

	#if USE_THREADS == 1
		tid = omp_get_thread_num();
	#else
		tid = 0;
	#endif

	galaxy *gal = (galaxy *) params;
	
	x = radius*cos(gal->theta_cyl[gal->index[tid]]);
	y = radius*sin(gal->theta_cyl[gal->index[tid]]);
	
	r_sph = sqrt(pow(x,2)+pow(y,2)+pow(gal->z[gal->index[tid]],2));

	pot = galaxy_potential_func(gal,gal->potential,gal->dx,gal->ngrid,x,y,gal->z[gal->index[tid]],1);
	if(gal->level_grid_zoom1>gal->level_grid) {
		sigma = 0.5*gal->dx_zoom1;
		transition_factor1 = 0.5*(1+erf((r_sph-(0.5*gal->boxsize_zoom1-0.5*sigma))/(sigma*sqrt(2))));
		transition_factor2 = 1-transition_factor1;	
		
		pot = transition_factor1*pot
			+ (galaxy_potential_func(gal,gal->potential_zoom1,gal->dx_zoom1,gal->ngrid_zoom1,x,y,gal->z[gal->index[tid]],0))*transition_factor2;

		if(gal->level_grid_zoom2>gal->level_grid_zoom1) {
			sigma = 0.5*gal->dx_zoom1;
			transition_factor1 = 0.5*(1+erf((r_sph-(0.5*gal->boxsize_zoom2-0.5*sigma))/(sigma*sqrt(2))));
			transition_factor2 = 1-transition_factor1;	
		
			pot = transition_factor1*pot
				+ (galaxy_potential_func(gal,gal->potential_zoom2,gal->dx_zoom2,gal->ngrid_zoom2,x,y,gal->z[gal->index[tid]],0))*transition_factor2;
		}
	}
	return pot;
}

// A wrapper for the galaxy potential function using the spherical radius.
double galaxyrsph_potential_wrapper_func(double r_sph, void *params) {
	
	double x, y, z, r_cyl, pot;
	double sigma,transition_factor1,transition_factor2;
	int tid;

	#if USE_THREADS == 1
		tid = omp_get_thread_num();
	#else
		tid = 0;
	#endif

	galaxy *gal = (galaxy *) params;
	
    z 			= cos(gal->phi_sph[gal->index[tid]])*r_sph;
    r_cyl 		= sqrt(r_sph*r_sph-z*z);
		
	x = r_cyl*cos(gal->theta_cyl[gal->index[tid]]);
	y = r_cyl*sin(gal->theta_cyl[gal->index[tid]]);
	
	pot = galaxy_potential_func(gal,gal->potential,gal->dx,gal->ngrid,x,y,z,1);	
	if(gal->level_grid_zoom1>gal->level_grid) {
		sigma = 8*gal->dx;
		transition_factor1 = 0.5*(1+erf((r_sph-(0.5*gal->boxsize_zoom1-0.5*sigma))/(sigma*sqrt(2))));
		transition_factor2 = 1-transition_factor1;	
		
		pot = transition_factor1*pot
			+ (galaxy_potential_func(gal,gal->potential_zoom1,gal->dx_zoom1,gal->ngrid_zoom1,x,y,z,0))*transition_factor2;
			
		if(gal->level_grid_zoom2>gal->level_grid_zoom1) {
			sigma = 8*gal->dx_zoom1;
			transition_factor1 = 0.5*(1+erf((r_sph-(0.5*gal->boxsize_zoom2-0.5*sigma))/(sigma*sqrt(2))));
			transition_factor2 = 1-transition_factor1;	
		
			pot = transition_factor1*pot
				+ (galaxy_potential_func(gal,gal->potential_zoom2,gal->dx_zoom2,gal->ngrid_zoom2,x,y,z,0))*transition_factor2;
		}
	}
	return pot;
}

// A wrapper for the galaxy potential function using the cartesian coordinate z.
double galaxyz_potential_wrapper_func(double z, void *params) {
	
	int tid;
	double pot, r_sph;
	double sigma,transition_factor1,transition_factor2;
	
	#if USE_THREADS == 1
		tid = omp_get_thread_num();
	#else
		tid = 0;
	#endif

	galaxy *gal = (galaxy *) params;
	
	r_sph = sqrt(pow(gal->x[gal->index[tid]],2)+pow(gal->y[gal->index[tid]],2)+pow(z,2));
	
	pot = galaxy_potential_func(gal,gal->potential,gal->dx,gal->ngrid,gal->x[gal->index[tid]],gal->y[gal->index[tid]],z,1);
	if(gal->level_grid_zoom1>gal->level_grid) {
		sigma = 8*gal->dx;
		transition_factor1 = 0.5*(1+erf((r_sph-(0.5*gal->boxsize_zoom1-0.5*sigma))/(sigma*sqrt(2))));
		transition_factor2 = 1-transition_factor1;	
	
		pot = transition_factor1*pot
			+ (galaxy_potential_func(gal,gal->potential_zoom1,gal->dx_zoom1,gal->ngrid_zoom1,gal->x[gal->index[tid]],gal->y[gal->index[tid]],z,0))*transition_factor2;
		if(gal->level_grid_zoom2>gal->level_grid_zoom1) {
			sigma = 8*gal->dx_zoom1;
			transition_factor1 = 0.5*(1+erf((r_sph-(0.5*gal->boxsize_zoom2-0.5*sigma))/(sigma*sqrt(2))));
			transition_factor2 = 1-transition_factor1;	
	
			pot = transition_factor1*pot
				+ (galaxy_potential_func(gal,gal->potential_zoom2,gal->dx_zoom2,gal->ngrid_zoom2,gal->x[gal->index[tid]],gal->y[gal->index[tid]],z,0))*transition_factor2;
		}
	}

	return pot;
}

// A wrapper for the second derivative of the potential function.
double potential_deriv_wrapper_func(double radius, void *params) {
	
	galaxy *gal = (galaxy *) params;
	
	if(radius==0.) {
		return 0.;
	} else {	
		return (pow(v_c_func(gal,radius),2.0))/(radius*kpc);
	}
}

// This function copies one galaxy to another.
void copy_potential(galaxy *gal_1, galaxy *gal_2, int info) {
	
	unsigned long int i, j, k;
	
	if(info == 1) printf("/////\tCopying potential grid \n");
	// Copy all the coordinate information.
	for (i = 0; i < gal_1->ngrid[0]*2; ++i) {
		for (j = 0; j < gal_1->ngrid[2]*2; ++j) {
			for(k = 0; k < gal_1->ngrid[2]*2; ++k){
				gal_2->potential[i][j][k] = gal_1->potential[i][j][k];
			}
		}
	}
	gal_2->potential_defined = 1;
	if(info == 1) printf("/////\tPotential grid copied \n");
	return;
}

// Function that computes the shift between a potential grid and a zoomed potential grid
void compute_potential_shift(galaxy *gal, double ***potential, double ***potential_zoom, double dx, double dx_zoom, int ngrid[3], int ngrid_zoom[3]) {
	double x,y,potential_shift,boxlen_zoom;
	int i,j,k,neval1,neval2;

	neval1 		= 100;
	neval2 		= 10000;
	boxlen_zoom = ngrid_zoom[0]*dx_zoom;
	for(i=0; i<neval1; i++){
		for(j=0; j<neval2; j++){
			x = (0.40+(0.1*i)/neval1)*boxlen_zoom*cos(2.0*pi*j/((double)neval2));
			y = (0.40+(0.1*i)/neval1)*boxlen_zoom*sin(2.0*pi*j/((double)neval2));
			potential_shift += (galaxy_potential_func(gal,potential,dx,ngrid,x,y,0.,1)-galaxy_potential_func(gal,potential_zoom,dx_zoom,ngrid_zoom,x,y,0.,0));
		}
	}
	potential_shift /= (double)(neval1*neval2);

	for (i = 0; i < 2*ngrid_zoom[0]; ++i) {
		for (j = 0; j < 2*ngrid_zoom[1]; ++j) {
			for (k = 0; k < 2*ngrid_zoom[2]; ++k) {
				potential_zoom[i][j][k] += potential_shift;
			}
		}
	}
	return;	
}