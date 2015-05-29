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
	infinity 	= 10.*gal->comp_scale_height[gal->selected_comp[tid]];
	gsl_integration_qag(&F,fabs(gal->z[gal->index[tid]]),fabs(gal->z[gal->index[tid]])+infinity,epsabs,epsrel,AllVars.GslWorkspaceSize,key,w,&integral,&error);

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

    rho 		= density_functions_pool(gal,radius,theta,z,gal->comp_jeans_mass_cut[gal->selected_comp[tid]],gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
	
    integrand 	= rho*galaxy_zforce_func(gal,z);
        
    return integrand;
}


// This function calculates the 1D velocity moment at a given spherical radius
double v2a_1D_func(galaxy *gal, gsl_integration_workspace *w, int component) {
    
	int status, tid;
	double integral, error, res, rho, infinity;
	size_t neval;
    
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
    
	gsl_function F;
    
	F.function 	= &dv2a_1D_func;
	F.params 	= gal;
	
	gal->selected_comp[tid] = component;
	infinity = 10.*gal->comp_scale_length[gal->selected_comp[tid]];
	gsl_integration_qag(&F,fabs(gal->r_sph[gal->index[tid]]),fabs(gal->r_sph[gal->index[tid]])+infinity,epsabs,epsrel,AllVars.GslWorkspaceSize,key,w,&integral,&error);
    
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
static double dv2a_1D_func(double r_sph, void *params) {
    
    double integrand, r_cyl, z, theta_cyl, rho;
	int tid;
    galaxy *gal = (galaxy *) params;
    
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
    
    theta_cyl 	= gal->theta_cyl[gal->index[tid]];
    z 			= cos(gal->phi_sph[gal->index[tid]])*r_sph;
    r_cyl 		= sqrt(r_sph*r_sph-z*z);
	
    rho 		= density_functions_pool(gal,r_cyl,theta_cyl,z,gal->comp_jeans_mass_cut[gal->selected_comp[tid]],gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
	
    integrand 	= rho*galaxy_rsphforce_func(gal,r_sph);
            
    return integrand;
}

// This function calculates PART of the phi axis velocity moment. In
// particular, it calculates the derivative of the radial velocity
// moment, which is equal to the z-axis velocity moment.
//
// This function uses the galaxy storage variables a lot! If you have
// edited the code to ill effect, check that the storage variables
// properly reassigned after this function.
//
// The derivative is calculated in the x,y plane using the
// 5-point stencil method.
double v2a_theta_func(galaxy *gal, double radius, int component) {
    
	double z, theta, h, v2a_theta, res, rho, abserr;
    double derivative, x, y;
	int tid;
	
	#if USE_THREADS == 1
    	tid = omp_get_thread_num();
	#else
   		tid = 0;
	#endif
	
    x 						= radius*cos(gal->theta_cyl[gal->index[tid]]);
	y 						= radius*sin(gal->theta_cyl[gal->index[tid]]); 
	// Set the derivative stepsize.
	h 						= get_h_value(gal,x,y,gal->z[gal->index[tid]],0,0);
	theta 					= gal->theta_cyl[gal->index[tid]];
    z 						= gal->z[gal->index[tid]];
    gal->selected_comp[tid] = component;
	derivative 				= deriv_forward(gal,radius,h,rho_v2a_r_func);
		
	rho 					= density_functions_pool(gal,radius,gal->theta_cyl[gal->index[tid]],z,0,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
	res 					= derivative*radius*kpc/rho;
	
	return res;
}

// Function which should be derivated to obtain the phi axis velocity moment
// Take a look to the part of the documentation explaining the velocity dispersion computation
double rho_v2a_r_func(double radius, void *params) {
	
	galaxy *gal = (galaxy *) params;
	double rho;
	double res, res2,save1, save2, save3, save4, save5;
	int tid;
	
	#if USE_THREADS == 1
    	tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif

	rho = density_functions_pool(gal,radius,gal->theta_cyl[gal->index[tid]],gal->z[gal->index[tid]],gal->comp_jeans_mass_cut[gal->selected_comp[tid]],gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);

	return rho*v2a_z_func(gal,w[tid],gal->selected_comp[tid]);
}


// This function checks that the azimuthal velocity dispersion verifies
// the lower limit of the Toomre stability cirterion imposed by the user
double v2a_z_toomre(galaxy *gal, double radius, double v2a_z, int component) {
    
	double gamma_sqrd, kappa_sqrd, h, force, dforcedr, Q, surface_density;
	double abserr, x, y;
	int tid;
    
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif

    x 				= radius*cos(gal->theta_cyl[gal->index[tid]]);
	y 				= radius*sin(gal->theta_cyl[gal->index[tid]]);    
	// Set the derivative step
	h 				= get_h_value(gal,x,y,gal->z[gal->index[tid]],0,0);
	// Calculate force and force derivative
	force 			= potential_deriv_wrapper_func(radius,gal);
	dforcedr 		= deriv_central(gal,radius,h,potential_deriv_wrapper_func);

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

// This function computes the Toomre criterion
double toomre(galaxy *gal, double radius, double v2a_z, int component) {
    
	double kappa_sqrd, h, force, dforcedr, Q, surface_density;
	double x, y;
	int tid;
    
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
    
    x 				= radius*cos(gal->theta_cyl[gal->index[tid]]);
	y 				= radius*sin(gal->theta_cyl[gal->index[tid]]);
	// Set the derivative step
	h 				= get_h_value(gal,x,y,gal->z[gal->index[tid]],0,0);
	// Calculate force and force derivative
	force 			= potential_deriv_wrapper_func(radius,gal);
	dforcedr 		= deriv_central(gal,radius,h,potential_deriv_wrapper_func);

	kappa_sqrd 		= 3.0*force/(radius*kpc) + dforcedr;
	// We take into account the contribution of the gas to the stability of the system
	surface_density = surface_density_func(gal,radius,gal->theta_cyl[gal->index[tid]],1,component);
	if(gal->comp_type[component]==0){
		Q = sqrt(gal->comp_cs_init[component]*fabs(kappa_sqrd))/(pi*G*surface_density);
	} else {
		Q = sqrt(v2a_z*fabs(kappa_sqrd))/(3.36*G*surface_density);
	}
    
	return Q;
}

// This function calculates the first velocity moment in the azimuthal
// direction for the a given flat component using the axisymmetric drift approximation
double sigma2_theta_disk_func(galaxy *gal, double radius, double v2a_z) {
    
	double gamma_sqrd, kappa_sqrd, h, force, dforcedr, surface_density;
	double res, x, y;
	int tid;
    
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif

    x 			= radius*cos(gal->theta_cyl[gal->index[tid]]);
	y 			= radius*sin(gal->theta_cyl[gal->index[tid]]);
	// Set the derivative step
	h 			= get_h_value(gal,x,y,gal->z[gal->index[tid]],0,0);
	// Calculate force and force derivative
	force 		= potential_deriv_wrapper_func(radius,gal);
	dforcedr 	= deriv_central(gal,radius,h,potential_deriv_wrapper_func);
	
	kappa_sqrd 	= 3.0*force/(radius*kpc)+dforcedr;
	gamma_sqrd 	= 4.0*force/(kappa_sqrd*radius*kpc);
    
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

	double v2_theta_gas, v_c2, pressure_force, save;
	double h, abserr, density_derivative, x, y;
	int tid;
    
    #if USE_THREADS == 1
    	tid = omp_get_thread_num();
	#else
    	tid = 0;
	#endif
    
	radius 					= fabs(radius);
	x 						= radius*cos(gal->theta_cyl[gal->index[tid]]);
	y 						= radius*sin(gal->theta_cyl[gal->index[tid]]);
	gal->selected_comp[tid] = component;
	// Set the derivative step
	h 						= get_h_value(gal,x,y,gal->z[gal->index[tid]],0,0);
	// Save z coordinate
	save 					= gal->z[gal->index[tid]];
	// Integrations done in the z=0 plane	
	gal->z[gal->index[tid]] = 0.;
	density_derivative 		= deriv_central(gal,radius,h,gas_density_wrapper_func);
	v_c2 					= pow(v_c_func(gal,radius),2.0);
	pressure_force 			= radius*kpc*(pow(gal->comp_cs_init[component],2.0)*density_derivative)/gas_density_wrapper_func(radius,gal);
	v2_theta_gas 			= v_c2 + pressure_force;
	if(v2_theta_gas<0) v2_theta_gas = 0.;
	// Restore z coordinate
	gal->z[gal->index[tid]] = save;
	
	return v2_theta_gas;
}

// A wrapper for the gas density function.
double gas_density_wrapper_func(double radius, void *params) {
    
	double theta, rho, x, y, sigma, smooth_in_factor, smooth_out_factor;
	int tid, component;
	galaxy *gal = (galaxy *) params;
    
	#if USE_THREADS == 1
    	tid = omp_get_thread_num();
	#else
    	tid = 0;
	#endif
	
	component 	= gal->selected_comp[tid];
	x 			= radius*cos(gal->theta_cyl[gal->index[tid]]);
	y 			= radius*sin(gal->theta_cyl[gal->index[tid]]);
    
	if(gal->pseudo[tid]) {
		rho = pseudo_density_gas_func(gal,fabs(radius),gal->theta_cyl[gal->index[tid]],gal->z[gal->index[tid]],0,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
	} else {
		rho = density_functions_pool(gal,fabs(radius),gal->theta_cyl[gal->index[tid]],gal->z[gal->index[tid]],0,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
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
// at a point with a cylindrical radius r. This is used to calculate the velocity dispersion.
//
// This function simply performs a five point differentiation around
// x,y,z in the cylindrical radius direction with a small stepsize.
//
// This function is specially formatted for the velocity dispersion
// routines.
double galaxy_rforce_func(galaxy *gal, double radius) {
	
	int tid;
	double force, h, abserr;
	double x,y,r_sph;
	double sigma,transition_factor1,transition_factor2;
	
	#if USE_THREADS == 1
		tid = omp_get_thread_num();
	#else
		tid = 0;
	#endif
	
	x = radius*cos(gal->theta_cyl[gal->index[tid]]);
	y = radius*sin(gal->theta_cyl[gal->index[tid]]);
		
	h = get_h_value(gal,x,y,gal->z[gal->index[tid]],0,0);
	
	force = deriv_central2(gal,radius,h,galaxyr_potential_wrapper_func);
	
	return force;
}

// This function calculates the force on a test particle due to the disk
// at a point with an height z. This is used to calculate the velocity dispersion.
//
// This function simply performs a five point differentiation around
// x,y,z in the z direction with a small stepsize.
//
// This function is specially formatted for the velocity dispersion
// routines.
double galaxy_zforce_func(galaxy *gal, double z) {
	
	int tid;
	double force, h, abserr;
	
	#if USE_THREADS == 1
		tid = omp_get_thread_num();
	#else
		tid = 0;
	#endif
	
	h = get_h_value(gal,gal->x[gal->index[tid]],gal->y[gal->index[tid]],z,0,0);

	force = deriv_central2(gal,z,h,galaxyz_potential_wrapper_func);
	
	return force;
}

// This function calculates the force on a test particle due to the global potential
// at a spherical radius r_psh. This is used to calculate the velocity dispersion.
//
// This function simply performs a five point differentiation around
// x,y,z in the spherical radius direction with a small stepsize.
//
// This function is specially formatted for the velocity dispersion
// routines.
double galaxy_rsphforce_func(galaxy *gal, double r_sph) {
	
	int tid;
	double force, h, abserr;
	
	#if USE_THREADS == 1
		tid = omp_get_thread_num();
	#else
		tid = 0;
	#endif
	
	h = get_h_value(gal,gal->x[gal->index[tid]],gal->y[gal->index[tid]],gal->z[gal->index[tid]],0,0);
	
	force = deriv_central(gal,r_sph,h,galaxyrsph_potential_wrapper_func);
	
	return force;
}