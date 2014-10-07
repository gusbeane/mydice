/*-----------------------------------------------------------------------------
 /
 / Filename: dice_structure.c
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

// Distribution function of the considered component
double fdistrib(galaxy *gal, double radius, double theta, double z, int model, int component) {
	//Density of the component times the integrand of spherical volume element
	return (fabs(radius))*density_functions_pool(gal,fabs(radius),theta,z,1,model,component);
}

double density_functions_pool(galaxy *gal, double radius, double theta, double z, int cut, int model, int component) {
	double h, z0, cut_dens, density, zpz0, m, mass, rho_crit, c, delta_c, rho0_nfw;
	// We consider only positive values
    if(sqrt(radius*radius+z*z) <= 0.) return 0.;

    z0 			= gal->comp_scale_height[component];
    h 			= gal->comp_scale_length[component];
	mass 		= gal->comp_mass[component];
    cut_dens 	= gal->comp_cut_dens[component];
    c			= gal->comp_concentration[component];
	//printf("%le\n",mass);

	// Variable for flattened distribution
	m 			= sqrt(pow(radius/h,2.0)+pow(z/z0,2.0));
	rho_crit 	= (mass/200.)*(3.0/(4.0*pi*pow(gal->r200*kpc,3.0)));
    delta_c 	= (200./3.)*pow(c,3.0)/(log(1.0+c)-c/(1.0+c));
    rho0_nfw 	= rho_crit*delta_c;
	// Select a disk model
    switch(model) {
		case 1:
			// Exponential disk + sech z-profile
			density = mass/(4.0*pi*z0*h*h*pow(kpc,3.0))*exp(-radius/h)/(cosh(z/z0)*cosh(z/z0));
			break;
		case 2:
			// Myamoto-Nagai profile
			zpz0 		= sqrt(z*z+z0*z0);
			density 	= mass*z0*z0*(h*radius*radius+(h+3.0*zpz0)*pow(h+zpz0,2.0))
            /(4.0*pi*pow(kpc,3.0)*pow(radius*radius+pow((h+zpz0),2.0),2.5)*pow(zpz0,3.0));
			break;
		case 3:
			// Exponential disk + exponential z-profile
			density = mass/(4.0*pi*z0*h*h*pow(kpc,3.0))*exp(-radius/h)*exp(-fabs(z)/z0);
			break;
		case 4:
			// Flattened Hernquist profile
			density = mass/(2.0*pi*h*z0*z0*pow(kpc,3.0))*(1.0/(m*pow((m + 1.0), 3.0)));
			break;
		case 5:
            // Flattened Plummer profile
            density = mass/(4.0*pi*h*z0*z0*pow(kpc,3.0))*(pow(1.0+pow(m,2.0),-2.5));
            break;
        case 6:
            // Jaffe profile
            density = mass/(4.0*pi*h*z0*z0*pow(kpc,3.0))*1.0/(pow(m*(m+1),2.0));
            break;
        case 7:
            // Isothermal profile
            density = 0.1*rho0_nfw/(pow(m,2.0));
            break;
        case 8:
            // NFW profile
            //rho_crit 	= (gal->m200/200.)*(3.0/(4.0*pi*pow(gal->r200*kpc,3.0)));
            //density 	= rho_crit*delta_c/((rsph/rs)*(1+rsph/rs));
			density 	= rho0_nfw/(m*pow(1.0+m,2.0));
            break;
    	case 9:
            // Burkert profile
			density 	= rho0_nfw/((m+1)*(m*m+1));
            break;
		default:
			// Exponential disk + sech z-profile
			density = mass/(4.0*pi*z0*h*h*pow(kpc,3.0))*exp(-radius/h)/(cosh(z/z0)*cosh(z/z0));
	}
	// Cutting the density at cutr & cutz
    if(cut && density<cut_dens) return 0.0;
	//if(cut && radius>gal->disk_cut)	return 0.0;
    //if(cut && fabs(z)>gal->disk_cut*gal->disk_temp)	return 0.0;
	return density;
}

void mcmc_metropolis_hasting(galaxy *gal, int component, int density_model) {
	unsigned long int i,j,start_part,npart;
	double pi_x, pi_y, q_x, q_y, prob, *radius, delta_pot;
	double theta, phi, randval, step, step_proposal, proposal, prop_r, prop_theta, prop_z, step_r, step_z, new_step_r, new_step_z;
	double acceptance;
	
	if(gal->comp_npart[component]>0) {
		// Use the Metropolis algorithm to place the disk particles.
		// We start the Monte Carlo Markov Chain with a realistic particle position
        theta 				= 2.0*pi*gsl_rng_uniform_pos(r[0]);
		i 					= gal->comp_start_part[component];
        gal->x[i] 			= 0.;
        gal->y[i] 			= 0.;
		gal->z[i] 			= 0.;
		gal->r_cyl[i]		= sqrt(gal->x[i]*gal->x[i]+gal->y[i]*gal->y[i]);
		gal->theta_cyl[i]	= atan2(gal->y[i],gal->x[i]);
		gal->r_sph[i]		= sqrt(gal->x[i]*gal->x[i]+gal->y[i]*gal->y[i]+gal->z[i]*gal->z[i]);
		gal->theta_sph[i]	= atan2(gal->y[i],gal->x[i]);
		gal->phi_sph[i]		= acos(gal->z[i]/gal->r_sph[i]);
		step_r 				= gal->comp_mcmc_step[component]*gal->comp_scale_length[component];
		step_z 				= gal->comp_mcmc_step[component]*gal->comp_scale_height[component];

		acceptance 			= 0.;
		// Burning period
        for (j = 0; j<100; ++j) {
			// Generating a proposal
			prop_r 		= gal->r_cyl[i] + gsl_ran_gaussian(r[0],step_r);
			prop_theta 	= 2.0*pi*gsl_rng_uniform_pos(r[0]);
			prop_z 		= gal->z[i] + gsl_ran_gaussian(r[0],step_z);
			// Computing the probability of the proposal
			step_proposal = step;
			pi_x = fdistrib(gal,gal->r_cyl[i],gal->theta_cyl[i],gal->z[i],density_model,component);
			pi_y = fdistrib(gal,prop_r,prop_theta,prop_z,density_model,component);
			q_x  = gsl_ran_gaussian_pdf(prop_r-gal->r_cyl[i],step_r)*gsl_ran_gaussian_pdf(prop_z-gal->z[i],step_z);
			q_y  = gsl_ran_gaussian_pdf(gal->r_cyl[i]-prop_r,step_r)*gsl_ran_gaussian_pdf(gal->z[i]-prop_z,step_z);
			prob = min(1.0,(pi_y/pi_x));//*(q_x/q_y));
			randval = gsl_rng_uniform_pos(r[0]);
			if (randval <= prob) {
				gal->r_cyl[i] 		= prop_r;
				gal->theta_cyl[i] 	= prop_theta;
				gal->z[i] 			= prop_z;
			}
		}
		// Updating the coordinate values
		gal->x[i] = gal->r_cyl[i]*cos(gal->theta_cyl[i]);
		gal->y[i] = gal->r_cyl[i]*sin(gal->theta_cyl[i]);
		// Updating the coordinate values
        gal->theta_cyl[i]	= atan2(gal->y[i],gal->x[i]);
        gal->r_sph[i]		= sqrt(gal->x[i]*gal->x[i]+gal->y[i]*gal->y[i]+gal->z[i]*gal->z[i]);
        gal->theta_sph[i]	= atan2(gal->y[i],gal->x[i]);
        gal->phi_sph[i]		= acos(gal->z[i]/gal->r_sph[i]);
		// Filling the Markov Chain
		start_part 	= gal->comp_start_part[component];
		npart 		= gal->comp_npart_pot[component];
        for (i = start_part + 1; i < start_part + npart; ++i) {
			// Generating a proposal
			prop_r 		= gal->r_cyl[i-1] + gsl_ran_gaussian(r[0],step_r);
			prop_theta 	= 2.0*pi*gsl_rng_uniform_pos(r[0]);
			prop_z 		= gal->z[i-1] + gsl_ran_gaussian(r[0],step_z);
			// Computing the probability of the proposal
			step_proposal = step;
			pi_x = fdistrib(gal,gal->r_cyl[i-1],gal->theta_cyl[i-1],gal->z[i-1],density_model,component);
			pi_y = fdistrib(gal,prop_r,prop_theta,prop_z,density_model,component);
			q_x  = gsl_ran_gaussian_pdf(prop_r-gal->r_cyl[i-1],step_r)*gsl_ran_gaussian_pdf(prop_z-gal->z[i-1],step_z);
			q_y  = gsl_ran_gaussian_pdf(gal->r_cyl[i-1]-prop_r,step_r)*gsl_ran_gaussian_pdf(gal->z[i-1]-prop_z,step_z);
			prob = min(1.0,(pi_y/pi_x));//*(q_x/q_y));
			randval = gsl_rng_uniform_pos(r[0]);
			if (randval <= prob) {
				gal->r_cyl[i] 		= prop_r;
				gal->theta_cyl[i] 	= prop_theta;
				gal->z[i] 			= prop_z;
				acceptance+=1.0;
			} else {
				gal->r_cyl[i] 		= gal->r_cyl[i-1];
				gal->theta_cyl[i] 	= gal->theta_cyl[i-1];
				gal->z[i] 			= gal->z[i-1];
			}
			// Updating the coordinate values
			gal->x[i] = gal->r_cyl[i]*cos(gal->theta_cyl[i]);
			gal->y[i] = gal->r_cyl[i]*sin(gal->theta_cyl[i]);
			// Updating the coordinate values
            gal->theta_cyl[i]	= atan2(gal->y[i],gal->x[i]);
            gal->r_sph[i]		= sqrt(gal->x[i]*gal->x[i]+gal->y[i]*gal->y[i]+gal->z[i]*gal->z[i]);
            gal->theta_sph[i]	= atan2(gal->y[i],gal->x[i]);
            gal->phi_sph[i]		= acos(gal->z[i]/gal->r_sph[i]);
		}
		acceptance /= gal->comp_npart_pot[component];
		printf("\t-> Acceptance = %lf \n",acceptance);
		if(acceptance<0.50) printf("/////\t\t\tWarning: MCMC acceptance is low!\n\t\t\tLower mcmc_step%d in the galaxy parameter file.\n",component);
		if(acceptance>0.90) printf("/////\t\t\tWarning: MCMC acceptance is high!\n\t\t\tIncrease mcmc_step%d in the galaxy parameter file.\n",component);
		if (AllVars.MeanPartDist) printf("/////\t\t\t[Mean inter-particle distance: %lf kpc]\n",mean_interparticle_distance(gal,component));
	}
	return;
}

// Surface density function
// The surface density is computed using a numerical integration
// of the density profile over z
double surface_density_func(galaxy *gal, double r, double theta, int cut, int component) {
	
	int status,tid;
	double surface_density,error,h;
	double flat;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(GSL_WORKSPACE_SIZE);
	gsl_function F;
	//h = gal->disk_scale_length;
	//surface_density = gal->disk_mass * exp(-radius/h) / (2.*pi*h*h);
	
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
	
	F.function = &integrand_density_func;
    F.params = gal;
	
	gal->storage1[0][tid] 	= r;
	gal->storage1[1][tid] 	= theta;
	gal->storage1[2][tid] 	= cut;
	
	gal->selected_comp[tid] = component; 
	cut 					= gal->comp_cut[component];
	flat 					= gal->comp_flat[component];
	gsl_integration_qag(&F,-cut*flat,cut*flat,epsabs,epsrel,GSL_WORKSPACE_SIZE,key,w,&surface_density,&error);
    
	gsl_integration_workspace_free(w);
	return surface_density*kpc;
}

// Integrand of the previous numerical integration
static double integrand_density_func(double z, void *params) {
	
	int tid;
	double r_temp,theta_temp,rho;
	int cut;
    
	galaxy *gal = (galaxy *) params;
	
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
	
	r_temp 		= gal->storage1[0][tid];
	theta_temp 	= gal->storage1[1][tid];
	cut 		= gal->storage1[2][tid];
	
	// This function should not be used with the pseudo-density function
	//if(gal->pseudo[tid] && gal->comp_type[gal->selected_comp[tid]]==0) {
	//	rho = pseudo_density_gas_func(gal,r_temp,0.0,z,cut,gal->selected_comp[tid]);
	//}
	//else {
		rho = density_functions_pool(gal,r_temp,theta_temp,z,cut,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
	//}
	return rho;
}

// This function calculates the force on a test particle due to the halo
// along the r axis. This is used to calculate the velocity dispersion.
double cumulative_mass_func(galaxy *gal, double radius, int component) {
	
	int status;
	double result,integral, error;
	int tid;
	
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
	
	gsl_integration_workspace *wk = gsl_integration_workspace_alloc(GSL_WORKSPACE_SIZE);
	gsl_function F;
	
	F.function 	= &d_cumulative_mass_func;
	F.params 	= gal;
		
	gal->selected_comp[tid] = component; 
	
	gsl_integration_qag(&F,0,radius,epsabs,epsrel,GSL_WORKSPACE_SIZE,key,wk,&integral,&error);
	gsl_integration_workspace_free(wk);
	
	result = 2*pi*integral;
	
	return result*kpc;
}

// This is the integrand for the previous function. It is setup to work with
// the GSL_qags structures.
static double d_cumulative_mass_func(double r, void *params) {
	
	int tid;
	double integrand;
	
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
	
	galaxy *gal = (galaxy *) params;
	
	integrand = surface_density_func(gal,r,0.,1,gal->selected_comp[tid])*r*kpc;
    
	return integrand;
}


// This function computes the density of the gas assuming vertical structure hydrostatic equilibirum.
// Using this function into a iterative algorithm allows to converge towards
// a complete equilibrium state
double pseudo_density_gas_func(galaxy *gal, double x, double y, double z, int cut, int component) {
	
	double h, z0, r, density, delta_pot, rho_0;
	if(gal->potential_defined != 1) {
		if(set_galaxy_potential(gal,0) != 0) printf("/////Unable to set the potential. Abording.\n");
	}
	// Computing the potential gradient between 0 and z
	// in cgs unit
	delta_pot = galaxy_potential_func(gal,x,y,z) - galaxy_potential_func(gal,x,y,0.);
    z0 = gal->comp_scale_height[component];
    h =  gal->comp_scale_length[component];
	r = sqrt(x*x+y*y);
	// Density in the xy-plane in cgs unit
	rho_0 = get_midplane_density(gal,x,y);
	// Hydrostatic equilibrium requires following density
	density = rho_0*exp(-delta_pot/(pow(gal->comp_cs_init[component],2.0)));

	// Return a pseudo density in cgs
	if(cut && density<gal->comp_cut_dens[component]) return 0.0;
	return density;
}

double midplane_density_gas_func(galaxy *gal, gsl_integration_workspace *w, double x, double y, int component) {
    int status,pseudo_save;
    int tid;
    double result, integral, error, initial_surface_density, radius;
    
    gsl_function F;
    
    #if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif

	
	gal->storage1[3][tid] = x;
	gal->storage1[4][tid] = y;
	
	radius 					= sqrt(x*x+y*y);
	pseudo_save 			= gal->pseudo[tid];
	gal->pseudo[tid] 		= 0;
	initial_surface_density = surface_density_func(gal,radius,0.,0,component);
	gal->pseudo[tid] 		= pseudo_save;

	if(initial_surface_density==0.) return 0.;
    
    F.function 	= &dmidplane_density_gas_func;
    F.params 	= gal;
    
    gsl_integration_qag(&F,-10.*gal->comp_scale_height[component],10.*gal->comp_scale_height[component],epsabs,epsrel,GSL_WORKSPACE_SIZE,key,w,&integral,&error);
    
	return initial_surface_density/(integral*kpc);
}

// This is the integrand for the previous function. It is setup to work with
// the GSL_qags structures.
static double dmidplane_density_gas_func(double z, void *params) {
	
	int tid;
	double integrand,delta_pot,x,y;
	galaxy *gal = (galaxy *) params;
	
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif

	x = gal->storage1[3][tid];
	y = gal->storage1[4][tid];
    
	delta_pot = galaxy_potential_func(gal,x,y,z)-galaxy_potential_func(gal,x,y,0.);
    
	integrand = exp(-delta_pot/pow(gal->comp_cs_init[gal->selected_comp[tid]],2.0));
	
	return integrand;
}

// This function fills a 2D grid with the value of the gas density in the midplane
// The use of a 2D grid intends to lower the computation time
void fill_midplane_dens_grid(galaxy *gal, int component) {
	int i,j;
	double x,y;
	int tid;
	
		
	#pragma omp parallel for schedule(guided) shared(gal) private(j,x,y)
	for (i=0;i<gal->ngrid_dens;i++) {
		for (j=0;j<gal->ngrid_dens;j++) {
			#if USE_THREADS == 1
	    		tid = omp_get_thread_num();
			#else
	    		tid = 0;
			#endif
			x = (i-((double)(gal->ngrid_dens/2)-0.5-0.5))*gal->space_dens[0];
			y = (j-((double)(gal->ngrid_dens/2)-0.5-0.5))*gal->space_dens[1];
			gal->midplane_dens[i][j] = midplane_density_gas_func(gal,w[tid],x,y,component);
		}
	}
}

//
double get_midplane_density(galaxy *gal, double x, double y) {
	int node_x,node_y;
	double xtemp,ytemp,dens1,dens2,dens3,dens4,dens_interp;
	double dx,dy,tx,ty;
    
	xtemp = x/gal->space_dens[0] + ((double)(gal->ngrid_dens/2)-0.5-0.5);
	ytemp = y/gal->space_dens[1] + ((double)(gal->ngrid_dens/2)-0.5-0.5);
	// Determine the parent node.
	node_x = floor(xtemp);
	node_y = floor(ytemp);
	    
	// Interpolation function to compute the potential.
	dens1 = gal->midplane_dens[node_x][node_y];
	dens2 = gal->midplane_dens[node_x+1][node_y];
	dens3 = gal->midplane_dens[node_x][node_y+1];
	dens4 = gal->midplane_dens[node_x+1][node_y+1];
	// CIC fractions
	dx = 1.0 - (xtemp - (double) node_x);
	dy = 1.0 - (ytemp - (double) node_y);
	tx = 1.0 - dx;
	ty = 1.0 - dy;
	// Return the interpolated potential.
	dens_interp = dx*dy*dens1 + tx*dy*dens2 + dx*ty*dens3 + tx*ty*dens4;
	
	return dens_interp;
	
}

// ------------------------------------------

/*// This function determines the scale length of the disk using
// the fitting formula provided in Mo et al. 1998, if the user want to use it.
// Otherwise, the user chose himself the value for the disk scale length.
double disk_scale_length_func(galaxy *gal) {
	
	int i;
	double base, power, c, f_conc, f, scale, v_c;
	
	c = gal->halo_concentration;
	f_conc = 2.0/3.0+ pow((c/21.5),0.7);
	base = (gal->j_d*gal->lambda)/(0.1*gal->m_d);
	power = (-0.06+2.71*gal->m_d+0.0047*gal->m_d/(gal->j_d*gal->lambda));
	f = pow(base,power)*(1.0-3.0*gal->m_d+5.2*gal->m_d*gal->m_d)*(1.0-0.019*c+0.00025*c*c+0.52/c);
	scale = (1.0/sqrt(2.0))*(gal->j_d/gal->m_d)*gal->lambda*gal->r200*(f/sqrt(f_conc));
	
	return scale;
}

// The next two functions are for calculating the angular momentum
// of the disk. The first is the integral of the disk and the
// second is the integrand. Please note that this returns the
// specific angular momentum, J_d/M_d, not J_d. The units are cm/s.
double j_d_func(galaxy *gal) {
	
	int status;
	double result, error, j_d;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	gsl_function F;
	
	F.function = &dj_d_func;
	F.params = gal;
	gsl_integration_qags(&F,0.0,100.0*gal->disk_scale_length,0.001,1.0e-7,1000,w,&result,&error);
	
	gsl_integration_workspace_free(w);
	
	return result;
}

// Integrand of the previous numerical integration
double dj_d_func(double radius, void *params) {
	
	galaxy *gal = (galaxy *) params;
	double v_c, scale;
	
	scale = gal->disk_scale_length;
	v_c = v_c_func(gal,radius);
	
	return v_c*(radius/scale)*(radius/scale)*exp(-radius/scale);
}*/

// This function returns the energy fraction function of MMW for the
// angular momentum.
double f_c_func(double c) {
	
	double a, upper, lower;
	
	upper = c*(1.0 - 1.0/pow(1.0+c,2.0) - 2.0*log(1.0+c)/(1.0+c));
	a = log(1.0+c) - c/(1.0+c);
	lower = (2.0*pow(a,2.0));
	
	return upper/lower;
}

// This is the integral g_c for determining the halo azimuthal circular
// velocity fraction f_s.
double g_c_func(double c) {
	
	int status;
	double result, error;
	gsl_integration_workspace *w
	= gsl_integration_workspace_alloc(1000);
	gsl_function F;
	
	F.function = &dg_c_func;
	gsl_integration_qags(&F,0.0,c,0.001,1.0e-7,1000,w,&result,&error);
	
	gsl_integration_workspace_free(w);
	
	return result;
}

// This is the integrand for the integral g_c above.
double dg_c_func(double x, void *params) {
	
	double a, b;
	
	a = sqrt(log(1.0 + x) - x/(1.0 + x));
	b = pow(x,1.5)/pow(1.0+x,2.0);
	
	return a*b;
}

// This is the azimuthal circular velocity fraction, f_s. It is needed
// to determine the azimuthal streaming velocity.
double f_s_func(double c, double lambda) {
	
	double f_s, g_c, f_c;
	
	g_c = g_c_func(c);
	f_c = f_c_func(c);
	f_s = 1.5*lambda*sqrt(2.0*c/f_c)*pow(log(1.0+c)-c/(1.0+c),1.5)/g_c;
	
	return f_s;
}

// This function computes the mean inter-particle distance
// in order to help the user to set an appropriate smoothing length
double mean_interparticle_distance(galaxy *gal, int component) {
	double dist,mean_dist;
	
	// Initialize interparticle distance
	dist = 0.;
	#pragma omp parallel shared(gal,dist)
	{
        int i,j;
        double closest_part;
        // Loop over particles
		#pragma omp for schedule(dynamic,100)
        for (i=gal->comp_start_part[component]; i < gal->comp_start_part[component]+gal->comp_npart[component]; i++) {
            // For a given particle, we compute the distances to all the other particles
            // with the same type
            closest_part = gal->boxsize;
            for (j=i+1; j < gal->comp_start_part[component]+gal->comp_npart[component]; j++) {
                closest_part = min(closest_part,sqrt(pow(gal->x[i]-gal->x[j],2)+pow(gal->y[i]-gal->y[j],2)+pow(gal->z[i]-gal->z[j],2)));
            }
            dist += closest_part;
        }
	}
	mean_dist = dist/gal->comp_npart[component];
	
	return mean_dist;
}

// This function lower the mass resolution of the particles
void lower_resolution(galaxy *gal) {
	unsigned long int i;
	int j;

	printf("///// Lowering particule mass resolution:\n");
	for (j=0; j<AllVars.MaxCompNumber; j++) {
    	if(gal->comp_npart[j]>0) printf("/////\tComponent %d -> particle mass %.2e solar mass\n",j,(gal->comp_cutted_mass[j]*1.0E10)/(gal->comp_npart[j]*unit_mass));
		// Filling the arrays of the &galaxy structure
		for (i = gal->comp_start_part[j]; i < gal->comp_start_part[j] + gal->comp_npart_pot[j]; i++) {
			gal->mass[i] 		= gal->comp_cutted_mass[j]/gal->comp_npart[j];
			gal->id[i] 			= i;
			gal->u[i] 			= gal->comp_u_init[j];
		}
	}
	return;
}
