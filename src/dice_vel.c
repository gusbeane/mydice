/*-----------------------------------------------------------------------------
 /
 / Filename: dice_vel.c
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


// This function calculates the z-axis velocity moment at a given radius
double v2a_z_func(galaxy *gal, gsl_integration_workspace *w, int component) {
    
	int status, tid;
	double integral, error, res, rho, infinity;
	size_t neval;
    
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
    
	gsl_function F;
    
	F.function 	= &dv2a_z_func;
	F.params 	= gal;
	
	gal->selected_comp[tid] = component;
	infinity = 10.*gal->comp_scale_height[gal->selected_comp[tid]];
	gsl_integration_qag(&F,fabs(gal->z[gal->index[tid]]),fabs(gal->z[gal->index[tid]])+infinity,epsabs,epsrel,GSL_WORKSPACE_SIZE,key,w,&integral,&error);
    
    rho 		= density_functions_pool(gal,gal->r_cyl[gal->index[tid]],gal->theta_cyl[gal->index[tid]],gal->z[gal->index[tid]],0,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
	res 		= integral*kpc/rho;
        
	if(AllVars.AcceptImaginary==1) {
		return fabs(res);
	} else {
		return (res>0.?res:0.);
	}
    
}

// This is the integrand for the previous function. It is setup to work with
// the GSL_qags structures.
static double dv2a_z_func(double z, void *params) {
    
    double integrand, radius, theta, rho;
	int tid;
    galaxy *gal = (galaxy *) params;
    
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
    
	radius 		= gal->r_cyl[gal->index[tid]];
	theta 		= gal->theta_cyl[gal->index[tid]];
    rho 		= density_functions_pool(gal,radius,theta,z,0,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
	
    integrand 	= rho*galaxy_zforce_func(gal,z);
        
    return integrand;
}

// This function calculates PART of the phi axis velocity moment. In
// particular, it calculates the derivative of the radial velocity
// moment, which is equal to the z-axis velocity moment.
//
// This function uses the galaxy storage variables a lot! If you have
// editted the code to ill effect, check that the storage variables
// properly reassigned after this function.
//
// The derivative is calculated in the x,y plane using the
// 5-point stencil method.
double v2a_theta_func(galaxy *gal, double radius, int component) {
    
	double z, theta, h, v2a_theta, res, rho, abserr;
    double derivative;
	int tid;
	
	#if USE_THREADS == 1
    	tid = omp_get_thread_num();
	#else
   		tid = 0;
	#endif
	
	// Set the derivative stepsize.
    h 					= 0.5*gal->space[0];
	theta 				= gal->theta_cyl[gal->index[tid]];
    z 					= gal->z[gal->index[tid]];
	derivative 			= deriv_central(gal,radius,h,rho_v2a_r_func);
	gal->selected_comp[tid] 	= component;
	
	rho 		= density_functions_pool(gal,radius,gal->theta_cyl[gal->index[tid]],z,0,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
	res 		= derivative*radius*kpc/rho;
	if(AllVars.AcceptImaginary==1) {
		return fabs(res);
    } else {
		return (res>0.?res:0.);
    }
}

// Function which should be derivated to obtain the phi axis velocity moment
// Take a look to the part of the documentation explaining the velocity dispersion computation
double rho_v2a_r_func(double radius, void *params) {
	
	galaxy *gal = (galaxy *) params;
	double rho;
	int tid;
	
	#if USE_THREADS == 1
    	tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
	
	rho = density_functions_pool(gal,radius,gal->theta_cyl[gal->index[tid]],gal->z[gal->index[tid]],0,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
	return rho*v2a_z_func(gal,w[tid],gal->selected_comp[tid]);
}


// This function checks that the azimuthal velocity dispersion verifies
// the lower limit of the Toomre stability cirterion imposed by the user
double v2a_z_toomre(galaxy *gal, double radius, double v2a_z, int component) {
    
	double gamma_sqrd, kappa_sqrd, h, force, dforcedr, Q, surface_density;
	double abserr;
	int tid;
    
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
    
	// Set the derivative step
	h 				= gal->space[0];
	// Calculate force and force derivative
	force 			= potential_deriv_wrapper_func(radius,gal);
	dforcedr 		= deriv_forward(gal,radius,h,potential_deriv_wrapper_func);

	kappa_sqrd 		= 3.0*force/(radius*kpc) + dforcedr;
	gamma_sqrd 		= 4.0*force/(kappa_sqrd*radius*kpc);
	// We take into account the contribution of the gas to the stability of the system
	surface_density = surface_density_func(gal,radius,gal->theta_cyl[gal->index[tid]],1,component);
	if(gal->comp_type[component]==0){
		Q = sqrt(gal->comp_cs_init[component]*fabs(kappa_sqrd))/(pi*G*surface_density);
	} else {
		Q = sqrt(v2a_z*fabs(kappa_sqrd))/(3.36*G*surface_density);
	}
	// Check the value of the Toomre parameter
	if(Q < gal->comp_Q_lim[component]){
		v2a_z 		= pow(gal->comp_Q_lim[component]*3.36*G*surface_density,2.0)/(fabs(kappa_sqrd));
		Q 			= gal->comp_Q_lim[component];
	}
	gal->comp_Q_min[component] = (Q < gal->comp_Q_min[component] && Q > 0.0) ? Q : gal->comp_Q_min[component];
    
	return v2a_z;
}

// This function calculates the first velocity moment in the azimuthal
// direction for the a given flat component using the axisymmetric drift approximation
double sigma2_theta_disk_func(galaxy *gal, double radius, double v2a_z) {
    
	double gamma_sqrd, kappa_sqrd, h, force, dforcedr, Q, surface_density;
	double res, abserr;
    
	// Set the derivative step
	h = gal->space[0];
	// Calculate force and force derivative
	force = potential_deriv_wrapper_func(radius,gal);
	dforcedr = deriv_forward(gal,radius,h,potential_deriv_wrapper_func);
	
	kappa_sqrd = 3.0*force/(radius*kpc)+dforcedr;
	gamma_sqrd = 4.0*force/(kappa_sqrd*radius*kpc);
    
	res = (v2a_z/gamma_sqrd);
    
	if(AllVars.AcceptImaginary==1) {
		return fabs(res);
    } else {
		return (res>0.?res:0.);
    }
}


// Gas theta squared velocity component
// This formulation is the same as the one used by Springel et al. 2005
double v2_theta_gas_func(galaxy *gal, double radius, double z, int component) {
	double v2_theta_gas, v_c2, pressure_force;
	double h, abserr, density_derivative;
	int tid;
    
    #if USE_THREADS == 1
    	tid = omp_get_thread_num();
	#else
    	tid = 0;
	#endif
    
	radius 					= fabs(radius);
	gal->selected_comp[tid] = component;
	// Set the derivative step
	h 						= 1.5*gal->space_dens[0];
	density_derivative 		= deriv_central(gal,radius,h,gas_density_wrapper_func);
	v_c2 					= pow(v_c_func(gal,radius),2.0);
	pressure_force 			= radius*kpc*(pow(gal->comp_cs_init[component],2.0)*density_derivative)/gas_density_wrapper_func(radius,gal);//-3.0*pow(gal->cs_init,2.0);
	// If the derivate fails
	if(pressure_force>0) 	pressure_force = 0.;
	v2_theta_gas 			= v_c2 + pressure_force;
	if(v2_theta_gas<0) v2_theta_gas = 0.;
	return v2_theta_gas;
}

// A wrapper for the gas density function.
double gas_density_wrapper_func(double radius, void *params) {
    
	double theta, rho, x, y;
	int tid;
	galaxy *gal = (galaxy *) params;
    
	#if USE_THREADS == 1
    	tid = omp_get_thread_num();
	#else
    	tid = 0;
	#endif
	
	x = radius*cos(gal->theta_cyl[gal->index[tid]]);
	y = radius*sin(gal->theta_cyl[gal->index[tid]]);
    
	if(gal->pseudo[tid]) {
		rho = pseudo_density_gas_func(gal,x,y,gal->z[gal->index[tid]],0,gal->selected_comp[tid]);
	}
	else {
		rho = density_functions_pool(gal,fabs(radius),0.,gal->z[gal->index[tid]],0,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
	}
	return rho;
}

// The circular velocity function. This function returns the velocity
// in units of cm/s.
double v_c_func(galaxy *gal, double radius) {
    
	double v_c, rforce;
	int tid;
	
	#if USE_THREADS == 1
    	tid = omp_get_thread_num();
	#else
    	tid = 0;
	#endif
	
	if(radius <= 0.0) return 0.;
	rforce 		= galaxy_rforce_func(gal,radius);
	v_c 		= sqrt(radius*kpc*rforce);
		
	if(rforce<0) return 0.;
	else return v_c;
}

// This function calculates the force on a test particle due to the disk
// at some point z. This is used to calculate the velocity dispersion.
//
// This function simply performs a five point differentiation around
// x,y,z in the z direction with a small stepsize.
//
// This function is specially formatted for the velocity dispersion
// routines.
// Disk_potential function returns the potential of the stellar disk + gaseous disk.
double galaxy_rforce_func(galaxy *gal, double radius) {
	
	double force, h, abserr;
	
	h = 1.5*gal->space[0];
	
	force = deriv_central(gal,radius,h,galaxyr_potential_wrapper_func);
	
	return force;
}

// This function calculates the force on a test particle due to the disk
// at some point z. This is used to calculate the velocity dispersion.
//
// This function simply performs a five point differentiation around
// x,y,z in the z direction with a small stepsize.
//
// This function is specially formatted for the velocity dispersion
// routines.
// Disk_potential function returns the potential of the stellar disk + gaseous disk.
double galaxy_zforce_func(galaxy *gal, double z) {
	
	double force, h, abserr;
	
	h = 1.5*gal->space[2];
	
	force = deriv_central(gal,z,h,galaxyz_potential_wrapper_func);
	
	return force;
}

// This function calculates the potential due to a galactic disk using Cloud-
// In-Cell mass assignment with vacuum (isolated) boundary conditions on a
// Cartesian grid. The method was adapted from the discussion found in Hockney
// and Eastwood, "Computer Simulation Using Particles," 1981.
int set_turbulent_grid(galaxy *gal, int component) {
	unsigned long int ii, start, end, l;
	// Loop variables
	int i, j, k;
	// Nodes coordinates
	int node_x, node_y, node_z;
	double space_x, space_y, space_z;
	double dx, dy, dz, tx, ty, tz, n;
	double x, y, z;
	double rand_vel,rad;
	// The particle-mesh grid size and the Green's function and potential
	// storage buffers. Global to keep from hitting the stack limit for large grid
	// size.
	double ***kernel_grid;
	//size_t p_x[gal->num_part[1]], p_y[gal->num_part[1]], p_z[gal->num_part[1]];
	fftw_plan fft_kernel, fft_turbulence, fftinv_turbulence;
	fftw_complex *kernel,*turbulence;
	

	// Setup fftw threads
	#if USE_THREADS == 1
		//if (verbose) printf("/////\t\t-Setting up FFT call with %d threads\n",AllVars.Nthreads);
		fflush(stdout);
		fftw_init_threads();
		fftw_plan_with_nthreads(AllVars.Nthreads);
		fflush(stdout);
	#endif
	
	// Allocate grid storage variables
	if (!(kernel = calloc(pow(gal->ngrid_turb_padded,3),sizeof(fftw_complex)))) {
		fprintf(stderr,"Unable to allocate space for kernel's function.\n");
		return -1;
	}
	if (!(turbulence = calloc(pow(gal->ngrid_turb_padded,3),sizeof(fftw_complex)))) {
		fprintf(stderr,"Unable to allocate space for turbulence buffer.\n");
		return -1;
	}
	
	if (!(kernel_grid=calloc(gal->ngrid_turb_padded,sizeof(double *)))) {
		fprintf(stderr,"Unable to create kernel's function x axis.\n");
		return -1;
	}
	// kernel grid allocation
	// x-axis
	for (i = 0; i < gal->ngrid_turb_padded; ++i) {
	        // y-axis
	        if (!(kernel_grid[i] = calloc(gal->ngrid_turb_padded,sizeof(double *)))) {
			fprintf(stderr,"Unable to create kernel's function y axis.\n");
           		return -1;
        	}
	        // z-axis
        	for (j = 0; j < gal->ngrid_turb_padded; ++j) {
			if (!(kernel_grid[i][j] = calloc(gal->ngrid_turb_padded,sizeof(double)))) {
				fprintf(stderr,"Unable to create kernel's function z axis.\n");
				return -1;
			}
		}
	}
	// Define which particles we should use in the turbulence computation
	// We use all the particles
	start 	= gal->comp_start_part[component];
	end 	= gal->comp_start_part[component]+gal->comp_npart_pot[component];

	// Allocate the fftw complex output value and the fftw dft plan.
	fft_turbulence	= fftw_plan_dft_3d(gal->ngrid_turb_padded,gal->ngrid_turb_padded,gal->ngrid_turb_padded,turbulence,turbulence,FFTW_FORWARD,FFTW_ESTIMATE);
	fft_kernel		= fftw_plan_dft_3d(gal->ngrid_turb_padded,gal->ngrid_turb_padded,gal->ngrid_turb_padded,kernel,kernel,FFTW_FORWARD,FFTW_ESTIMATE);
	
	// Normalization constant
	// See FFTW reference guide for more details
	n = (int)pow(gal->ngrid_turb_padded,3.0);
	
	// Check for bad grids
	if (gal->ngrid_turb_padded <= 0) {
		fprintf(stderr,"\t\tGrid dimensions must be greater than zero! (ngrid=%d)\n",gal->ngrid_turb_padded);
		return -1;
	}
	
	// Sort the position arrays and figure out the spacing between grid points.
	// Subtract 2 for a.) the C offset and b.) the CIC offset. Finally, store
	// the values in the galaxy for later use.
	//gsl_sort_index(p_x,gal->x,1,gal->num_part[0]);
	//gsl_sort_index(p_y,gal->y,1,gal->num_part[0]);
	//gsl_sort_index(p_z,gal->z,1,gal->num_part[0]);
	space_x = gal->space_turb[0];
	space_y = gal->space_turb[1];
	space_z = gal->space_turb[2];
	// Print the turbulence in the xy-plane for z = 0 if the option is set.
	//if (verbose) printf("/////\t\t-Grid cell spacings [kpc]: dx = %.3f dy = %.3f dz = %.3f\n",space_x,space_y,space_z);
	fflush(stdout);

	// Initialization loop
	for (i = 0; i < gal->ngrid_turb_padded; ++i) {
		for (j = 0; j < gal->ngrid_turb_padded; ++j) {
			for (k = 0; k < gal->ngrid_turb_padded; ++k) {
				gal->turbulence[i][j][k] = gsl_ran_gaussian(r[0],1.0);
			}
		}
	}
	
	// Define kernel's function.
	// These are the grid points as measured from the center of kernel's function
	// and the local value of the truncation function. The density is also packed into a
	// buffer here.
	//
	// kernel's function is defined eight times here, once per octant, to take care of the
	// isolated (vacuum) boundary conditions. See Hockney and Eastwood, 1980, ch. 6 for a
	// discussion. The octants start in the lower left at (p,q) = (0,0) and progress
	// counter clockwise.
	for (i = 0; i < gal->ngrid_turb_padded/2; ++i) {
	    for (j = 0; j < gal->ngrid_turb_padded/2; ++j) {
	        #pragma omp parallel for private(dx, dy, dz) shared(kernel_grid,i,j)
			for (k = 0; k < gal->ngrid_turb_padded/2; ++k) {
		    	dx = sqrt(pow((double)(i+0.5),2.0));
		        dy = sqrt(pow((double)(j+0.5),2.0));
		    	dz = sqrt(pow((double)(k+0.5),2.0));
		        rad = sqrt(dx*dx+dy*dy+dz*dz);
				// Octant 1
		        kernel_grid[i][j][k] = 1.0/(pow(gal->comp_turb_scale[component]*sqrt(2.0*pi),3.0))*exp(-pow(rad,2.0)/(2.0*pow(gal->comp_turb_scale[component],2.0)));
		        // Octant 2
                kernel_grid[gal->ngrid_turb_padded-1-i][j][k] = kernel_grid[i][j][k];
                // Octant 3
                kernel_grid[gal->ngrid_turb_padded-1-i][gal->ngrid_turb_padded-1-j][k] = kernel_grid[i][j][k];
                // Octant 4
                kernel_grid[i][gal->ngrid_turb_padded-1-j][k] = kernel_grid[i][j][k];
                // Octant 5
                kernel_grid[i][j][gal->ngrid_turb_padded-1-k] = kernel_grid[i][j][k];
                // Octant 6
                kernel_grid[gal->ngrid_turb_padded-1-i][j][gal->ngrid_turb_padded-1-k] = kernel_grid[i][j][k];
                // Octant 7
                kernel_grid[gal->ngrid_turb_padded-1-i][gal->ngrid_turb_padded-1-j][gal->ngrid_turb_padded-1-k] = kernel_grid[i][j][k];
                // Octant 8
                kernel_grid[i][gal->ngrid_turb_padded-1-j][gal->ngrid_turb_padded-1-k] = kernel_grid[i][j][k];
			}
		}
	}
	
	
	// Pack kernel's function and the density into 1D arrays
	l = 0;
	for (i = 0; i < gal->ngrid_turb_padded; ++i) {
	    for (j = 0; j < gal->ngrid_turb_padded; ++j) {
			for (k = 0; k < gal->ngrid_turb_padded; ++k) {
		       	kernel[l] = kernel_grid[i][j][k];
		       	turbulence[l] = gal->turbulence[i][j][k];
				l++;
			}
		}
	}
    	
	// Perform the fourier transforms. Density first, kernel's function second.
	fftw_execute(fft_turbulence);
	fftw_execute(fft_kernel);
	// FFT is computed, we can free the memory
	fftw_destroy_plan(fft_turbulence);
	fftw_destroy_plan(fft_kernel);
	// Allocating memory for the inverse fourier computation
	fftinv_turbulence = fftw_plan_dft_3d(gal->ngrid_turb_padded,gal->ngrid_turb_padded,gal->ngrid_turb_padded,turbulence,turbulence,FFTW_BACKWARD,FFTW_ESTIMATE);
	// Multiply the density by kernel's function to find the k-space turbulence and
	// invert for the real potenital. Second, normalize the system and, finally,
	// put the turbulence information into the grid.
	for (i = 0; i < n; ++i) {
	        // Convolve the turbulence
	        turbulence[i] = kernel[i]*turbulence[i];
	}
	// Inversion
	fftw_execute(fftinv_turbulence);
	fftw_destroy_plan(fftinv_turbulence);
	// Normalization
	double stddev;
	for (i = 0; i < n; ++i) {
	        turbulence[i] = turbulence[i]/n;
	        stddev += pow(turbulence[i],2.0);
	}
	stddev = sqrt(stddev/n);

	l = 0;
	for (i = 0; i < gal->ngrid_turb_padded; ++i) {
	        for (j = 0; j < gal->ngrid_turb_padded; ++j) {
				for (k = 0; k < gal->ngrid_turb_padded; ++k) {
					// Fix the grid info
					gal->turbulence[i][j][k] = turbulence[l]/stddev*gal->comp_turb_sigma[component];
					l++;
				}
		}
	}
	// Free fftw plan.
	// Kill the storage arrays since they are no longer needed.
	fftw_free(kernel);
	fftw_free(turbulence);
	for (i = 0; i < gal->ngrid_turb_padded; ++i) {
	        for (j = 0; j < gal->ngrid_turb_padded; ++j) {
			free(kernel_grid[i][j]);
		}
		free(kernel_grid[i]);
	}
	free(kernel_grid);
	
	#if USE_THREADS == 1
		fftw_cleanup_threads();
	#endif
	// Flag the galaxy structure
	return 0;
}


// This function calculates the potential due to the disk at a point x,y,z
// by interpolating between grid points on the particle mesh. The
// interpolation routine uses the CIC kernel, which oddly enough is just
// a bilinear interpolation scheme...
//
// If the point lies off of the particle mesh, it approximates the potential
// as a function of 1/r.
double galaxy_turbulence_func(galaxy *gal, double x, double y, double z) {
	
	int node_x, node_y, node_z, offset;
	double turb1, turb2, turb3, turb4, turb5, turb6, turb7, turb8;
	double a, r_p, r_max, dx, dy, dz, tx, ty, tz, turbulence;
	double xtemp,ytemp,ztemp,theta,phi,xmax,ymax,zmax,rnorm;
		
	// Scale the coordinates
	r_p = sqrt(x*x + y*y + z*z);
	xtemp = x/gal->space_turb[0] + ((double)(gal->ngrid_turb_padded/2)-0.5-0.5);
	ytemp = y/gal->space_turb[1] + ((double)(gal->ngrid_turb_padded/2)-0.5-0.5);
	ztemp = z/gal->space_turb[2] + ((double)(gal->ngrid_turb_padded/2)-0.5-0.5);
	offset = 0;
	// Determine the parent node.
	node_x = floor(xtemp);
	node_y = floor(ytemp);
	node_z = floor(ztemp);
	r_max = gal->space[0]*((double)(gal->ngrid_turb_padded/4)-0.5-0.5-offset);
	
	// Check to see if (x,y,z) is a grid point.
	if (xtemp == (double) node_y && ytemp == (double) node_y && ztemp == (double) node_z) {
		// If (x,y,z) is a grid point, return its potential.
		turbulence = gal->turbulence[node_x][node_y][node_z];
	} else {
		// If (x,y,z) is not a grid point, use the CIC
		// interpolation function to calculate the potential.
		turb1 = gal->turbulence[node_x][node_y][node_z];
		turb2 = gal->turbulence[node_x+1][node_y][node_z];
		turb3 = gal->turbulence[node_x][node_y+1][node_z];
		turb4 = gal->turbulence[node_x][node_y][node_z+1];
		turb5 = gal->turbulence[node_x][node_y+1][node_z+1];
		turb6 = gal->turbulence[node_x+1][node_y+1][node_z];
		turb7 = gal->turbulence[node_x+1][node_y][node_z+1];
		turb8 = gal->turbulence[node_x+1][node_y+1][node_z+1];
		// CIC fractions
		dx = 1.0 - (xtemp - (double) node_x);
		dy = 1.0 - (ytemp - (double) node_y);
		dz = 1.0 - (ztemp - (double) node_z);
		tx = 1.0 - dx;
		ty = 1.0 - dy;
		tz = 1.0 - dz;
		// Return the interpolated potential.
		turbulence = dx*dy*dz*turb1 + tx*dy*dz*turb2 +
		dx*ty*dz*turb3 + dx*dy*tz*turb4 +
		dx*ty*tz*turb5 + tx*ty*dz*turb6 +
		tx*dy*tz*turb7 + tx*ty*tz*turb8;
	}
	return turbulence;
}