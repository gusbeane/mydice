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


double density_functions_pool(galaxy *gal, double radius, double theta, double z, int cut, int model, int component) {
	double h, hx, hy, hz, h_cut, hx_cut, hy_cut, hz_cut, density, w, k, l, m, n, o, alpha, smooth_factor, sigma, r_sph;
	double x,y;
	// We consider only positive values
    if(sqrt(radius*radius+z*z) < 0.) return 0.;

    h 			= gal->comp_scale_length[component];
    hx			= gal->comp_scale_length[component]*gal->comp_flatx[component];
    hy			= gal->comp_scale_length[component]*gal->comp_flaty[component];
	hz 			= gal->comp_scale_length[component]*gal->comp_flat[component];
	h_cut 		= gal->comp_cut[component];
	hx_cut 		= gal->comp_cut[component]*gal->comp_flatx[component];
	hy_cut 		= gal->comp_cut[component]*gal->comp_flaty[component];
	hz_cut 		= gal->comp_cut[component]*gal->comp_flat[component];
    alpha		= gal->comp_alpha[component];
	
	x 			= radius*cos(theta);
	y 			= radius*sin(theta);
	r_sph 		= sqrt(x*x+y*y+z*z);

	k			= sqrt(pow(z/hz,2.0));
	l			= sqrt(pow(x/hx,2.0)+pow(y/hy,2.0));
	m 			= sqrt(pow(x/hx,2.0)+pow(y/hy,2.0)+pow(z/hz,2.0));
	n 			= sqrt(pow(x/hx_cut,2.0)+pow(y/hy_cut,2.0)+pow(z/hz_cut,2.0));
	o 			= sqrt(pow(x/hx_cut,2.0)+pow(y/hy_cut,2.0));
	
	// Select a disk model
    switch(model) {
		case 1:
			// Exponential disk + sech z-profile
			if(strcmp(gal->comp_profile_name[component],"")==0)
			strcpy(gal->comp_profile_name[component],"Exponential disk / sech z");
			density = gal->comp_scale_dens[component]*exp(-l)/(pow(cosh(z/hz),2));
			break;
		case 2:
			// Myamoto-Nagai profile
			if(strcmp(gal->comp_profile_name[component],"")==0)
			strcpy(gal->comp_profile_name[component],"      Myamoto-Nagai      ");
			w = sqrt(z*z+hz*hz);
			density = gal->comp_scale_dens[component]*(h*radius*radius+(h+3.0*w)*pow(h+w,2.0))
			/(pow(radius*radius+pow((h+w),2.0),2.5)*pow(w,3.0));
			break;
		case 3:
			// Exponential disk + exponential z-profile
			if(strcmp(gal->comp_profile_name[component],"")==0)
			strcpy(gal->comp_profile_name[component],"Exponential disk / exp z ");
			density = gal->comp_scale_dens[component]*exp(-l)*exp(-fabs(z)/hz);
			break;
		case 4:
			// Hernquist profile
			if(strcmp(gal->comp_profile_name[component],"")==0)
			strcpy(gal->comp_profile_name[component],"        Hernquist        ");
			density = gal->comp_scale_dens[component]*1.0/(m*pow((m+1.0),3.0));
			break;
		case 5:
            // Plummer profile
			if(strcmp(gal->comp_profile_name[component],"")==0)
			strcpy(gal->comp_profile_name[component],"          Plummer        ");
			density = gal->comp_scale_dens[component]*pow(1.0+pow(m,2.0),-2.5);
			break;
        case 6:
            // Jaffe profile
			if(strcmp(gal->comp_profile_name[component],"")==0)
			strcpy(gal->comp_profile_name[component],"          Jaffe          ");
			density = gal->comp_scale_dens[component]*1.0/(pow(m*(m+1),2.0));
            break;
        case 7:
            // Isothermal profile
			if(strcmp(gal->comp_profile_name[component],"")==0)
			strcpy(gal->comp_profile_name[component],"       Isothermal        ");
			density = gal->comp_scale_dens[component]*1.0/(pow(m,2.0));
            break;
        case 8:
            // NFW profile
			if(strcmp(gal->comp_profile_name[component],"")==0)
			strcpy(gal->comp_profile_name[component],"           NFW           ");
			density = gal->comp_scale_dens[component]*1.0/(m*pow(1.0+m,2.0));
            break;
    	case 9:
            // Burkert profile
			if(strcmp(gal->comp_profile_name[component],"")==0)
			strcpy(gal->comp_profile_name[component],"         Burkert         ");
			density = gal->comp_scale_dens[component]*1.0/((m+1)*(pow(m,2.0)+1));
            break;
        case 10:
            // Einasto profile
			if(strcmp(gal->comp_profile_name[component],"")==0)
			strcpy(gal->comp_profile_name[component],"         Einasto         ");
			density = gal->comp_scale_dens[component]*exp(-pow(m,alpha));
            break;
        case 11:
            // Mestel profile
			if(strcmp(gal->comp_profile_name[component],"")==0)
			strcpy(gal->comp_profile_name[component],"         Mestel          ");
			density = gal->comp_scale_dens[component]*(acos(o)/radius)*exp(-fabs(z)/hz);
			if(o>1) density = 0.;
            break;
        case 12:
            // Kalnajs profile
			if(strcmp(gal->comp_profile_name[component],"")==0)
			strcpy(gal->comp_profile_name[component],"         Kalnajs         ");
			density = gal->comp_scale_dens[component]*pow(1-o*o,0.5)*exp(-fabs(z)/hz);
			if(o>1) density = 0.;
            break;
		default:
			fprintf(stderr,"[Error] model%d=%d is not a valid value\n",component+1,model);
			exit(0);
	}
	
	sigma 			= 0.05;
	smooth_factor 	= 1-0.5*(1+erf((n-1.0)/(sigma*sqrt(2))));
	
	if(gal->dens_gauss_sigma>0. && gal->comp_dens_gauss[component]==1 && gal->gaussian_field_defined){
		if(radius>gal->comp_cut[component]) {
			density = 0.0;
		} else {
			density *= (galaxy_gaussian_field_func(gal,x,y,z)*gal->dens_gauss_sigma+1.0);
			if(density<0.) density = 0.0;
		}
	}
	// Cutting the density
    if(cut==1) 	density *= smooth_factor;
    // Old cutting method
    //if(cut==1 && density<gal->comp_cut_dens[component]) 	density = 0.;
	//if(cut==1 && radius<gal->comp_cut_in[component]) 		density = 0.;
	// Unit is 10e10 solar mass / kpc^3
	return density;
}

double density_functions_stream_pool(stream *st, double radius, double theta, double z, int model, int component) {
	double alpha, h, rs, density, x, y;
	// We consider only positive values
    alpha		= pi/180.*st->comp_opening_angle[component];
    h			= st->comp_length[component];
    rs			= st->comp_scale[component]*fabs(z)*atan(pi/180.*st->comp_opening_angle[component]/2.);

	x = radius*cos(theta);
	y = radius*sin(theta);
	
	// Select a disk model
    switch(model) {
		case 1:
			// Uniform density cone
			if(strcmp(st->comp_profile_name[component],"")==0)
			strcpy(st->comp_profile_name[component],"      Uniform cone       ");
			if(fabs(radius)<z*atan(alpha) && z>0 && z<h) {
				density = st->comp_dens[component]/unit_nh;
			} else {
				density = 0.;
			}	
			break;
		case 2:
			// Exponential profile density cone
			if(strcmp(st->comp_profile_name[component],"")==0)
			strcpy(st->comp_profile_name[component],"        Exp cone         ");
			if(fabs(radius)<z*atan(alpha) && z>0 && z<h) {
				density = exp(-fabs(radius)/rs)*st->comp_dens[component]/unit_nh;
			} else {
				density = 0.;
			}	
			break;
		default:
			fprintf(stderr,"/////\t\t\tSpecify a valid model for component %d\n",component);
			exit(0);		
	}
	if(st->dens_gauss_sigma>0. && st->comp_dens_gauss[component]==1 && st->gaussian_field_defined){
		density *= (stream_gaussian_field_func(st,x,y,z)*st->dens_gauss_sigma+1.0);
		if(density<0.) density = 0.0;
	}
	return density;
}

void mcmc_metropolis_hasting(galaxy *gal, int component, int density_model) {
	unsigned long int i,j,start_part,npart;
	double pi_x, pi_y, q_x, q_y, prob, *radius;
	double theta, phi, randval, proposal, prop_r, prop_theta, prop_z, step_r, step_z, new_step_r, new_step_z;
	double acceptance;
	
	if(gal->comp_npart[component]>0) {
		// Use the Metropolis algorithm to place the disk particles.
		// We start the Monte Carlo Markov Chain with a realistic particle position
		step_r 				= gal->comp_mcmc_step[component]*gal->comp_scale_length[component];
		step_z 				= gal->comp_mcmc_step[component]*gal->comp_scale_height[component];
        theta 				= 2.0*pi*gsl_rng_uniform_pos(r[0]);
		i 					= gal->comp_start_part[component];
        gal->x[i] 			= max(gal->comp_cut_in[component],step_r)*cos(theta);
        gal->y[i] 			= max(gal->comp_cut_in[component],step_r)*sin(theta);
		gal->z[i] 			= 0.;
		gal->r_cyl[i]		= sqrt(gal->x[i]*gal->x[i]+gal->y[i]*gal->y[i]);
		gal->theta_cyl[i]	= atan2(gal->y[i],gal->x[i]);
		gal->r_sph[i]		= sqrt(gal->x[i]*gal->x[i]+gal->y[i]*gal->y[i]+gal->z[i]*gal->z[i]);
		gal->theta_sph[i]	= atan2(gal->y[i],gal->x[i]);
		gal->phi_sph[i]		= acos(gal->z[i]/gal->r_sph[i]);

		acceptance 			= 0.;
		// Burning period
        for (j = 1; j<(int)(0.1*gal->comp_npart_pot[component]); ++j) {
			// Generating a proposal
			prop_r 		= gal->r_cyl[i] + gsl_ran_gaussian(r[0],step_r);
			prop_theta 	= 2.0*pi*gsl_rng_uniform_pos(r[0]);
			prop_z 		= gal->z[i] + gsl_ran_gaussian(r[0],step_z);
			// Distribution function of the considered component
			// Density of the component times the integrand of spherical volume element
			pi_x = fabs(gal->r_cyl[i])*density_functions_pool(gal,fabs(gal->r_cyl[i]),gal->theta_cyl[i],gal->z[i],1,density_model,component);
			pi_y = fabs(prop_r)*density_functions_pool(gal,fabs(prop_r),prop_theta,prop_z,1,density_model,component);
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
        for (i = gal->comp_start_part[component]+1; i < gal->comp_start_part[component]+gal->comp_npart_pot[component]; ++i) {
			// Generating a proposal
			prop_r 		= gal->r_cyl[i-1] + gsl_ran_gaussian(r[0],step_r);
			prop_theta 	= 2.0*pi*gsl_rng_uniform_pos(r[0]);
			prop_z 		= gal->z[i-1] + gsl_ran_gaussian(r[0],step_z);
			// Distribution function of the considered component
			// Density of the component times the integrand of spherical volume element
			pi_x = fabs(gal->r_cyl[i-1])*density_functions_pool(gal,fabs(gal->r_cyl[i-1]),gal->theta_cyl[i-1],gal->z[i-1],1,density_model,component);
			pi_y = fabs(prop_r)*density_functions_pool(gal,fabs(prop_r),prop_theta,prop_z,1,density_model,component);
			q_x  = gsl_ran_gaussian_pdf(prop_r-gal->r_cyl[i-1],step_r)*gsl_ran_gaussian_pdf(prop_z-gal->z[i-1],step_z);
			q_y  = gsl_ran_gaussian_pdf(gal->r_cyl[i-1]-prop_r,step_r)*gsl_ran_gaussian_pdf(gal->z[i-1]-prop_z,step_z);
			prob = min(1.0,(pi_y/pi_x));//*(q_x/q_y));
			randval = gsl_rng_uniform_pos(r[0]);
			// Proposal accepted
			if (randval <= prob) {
				gal->r_cyl[i] 		= prop_r;
				gal->theta_cyl[i] 	= prop_theta;
				gal->z[i] 			= prop_z;
				gal->rho[i]			= pi_y/fabs(gal->r_cyl[i]);
				acceptance+=1.0;
			// Proposal rejected, the particle keeps the same postion
			} else {
				// To avoid tree codes to break, add a slight displacement for the particle
				gal->r_cyl[i] 		= gal->r_cyl[i-1];
				gal->theta_cyl[i] 	= gal->theta_cyl[i-1];
				gal->z[i] 			= gal->z[i-1];
            	gal->rho[i]			= pi_x/fabs(gal->r_cyl[i]);
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
		printf("[acceptance=%.2lf]\n",acceptance);
		if(acceptance<0.50) printf("/////\t\t[Warning] MCMC acceptance is low -> Decrease mcmc_step%d\n",component+1);
		if(acceptance>0.90) printf("/////\t\t[Warning] MCMC acceptance is high -> Increase mcmc_step%d\n",component+1);
		if (AllVars.MeanPartDist) printf("/////\t\t\tMean inter-particle distance -> %lf [kpc]\n",mean_interparticle_distance(gal,component));
	}
	return;
}

// This function modify the z component of the gaseous particles according to an iterative
// algorithm to reach the hydrostatic equilibrium of the disk. This function was inspired by
// Springel, Di Matteo er al. 2005 method.
int set_hydro_equilibrium(galaxy *gal, int component, int n_iter) {
	
	unsigned long int i,j,k;
	int density_model, neval;
	double mu, z0, pi_x, pi_y, q_x, q_y, prob, *radius;
	double theta, phi, randval, x, y, z, step, previous_step, prop_r, prop_theta, prop_z, acceptance;
	double step_r, step_z, prev_step_z;
    double cut_dens,theta_max_dens;
    double xtemp,ytemp,rtemp;
    
	printf("/////\tTargeting gas azimuthal hydrostatic equilibrium\n");

	fflush(stdout);
	// Allow to use pseudo density functions
	for(i=0;i<AllVars.Nthreads;i++) gal->pseudo[i] = 1;
	density_model = gal->comp_model[component];
	// Looking for gas particles
	if(gal->comp_type[component]==0 && gal->comp_npart_pot[component]>0) {
		printf("/////\t\t- Component %d -> recomputing gas particles position\n",component+1);
		// Starting equilibrium iterations
		for(j = 0; j < n_iter; ++j) {
			printf("/////\t\t\tIteration %d -> evaluating potential [%d threads]",j+1,AllVars.Nthreads);
			fflush(stdout);
			// Now that we've got gas particles positions, we can compute the full potential.
			if(set_galaxy_potential(gal,gal->potential,gal->dx,gal->ngrid,0) != 0) {
				fprintf(stderr,"\n[Error] Unable to set the potential\n");
				exit(0);
			}
			// Zoom 1
			if(gal->level_grid_zoom1>gal->level_grid){
				if(set_galaxy_potential(gal,gal->potential_zoom1,gal->dx_zoom1,gal->ngrid_zoom1,0) != 0) {
					printf("[Error] Unable to set the zoomed potential\n");
					exit(0);
				}
				compute_potential_shift(gal,gal->potential,gal->potential_zoom1,gal->dx,gal->dx_zoom1,gal->ngrid,gal->ngrid_zoom1);
				// Zoom 2
				if(gal->level_grid_zoom2>gal->level_grid_zoom1){
					if(set_galaxy_potential(gal,gal->potential_zoom2,gal->dx_zoom2,gal->ngrid_zoom2,1) != 0) {
						printf("[Error] Unable to set the zoomed potential\n");
						exit(0);
					}
					compute_potential_shift(gal,gal->potential_zoom1,gal->potential_zoom2,gal->dx_zoom1,gal->dx_zoom2,gal->ngrid_zoom1,gal->ngrid_zoom2);
				}
			}
			
			fflush(stdout);
			printf(" / Computing midplane dens");
			fflush(stdout);
			fill_midplane_dens_grid(gal,component);
			// Use the Metropolis algorithm to place the disk particles.
			// We start the Monte Carlo Markov Chain with a realistic particle position
			step_r 				= gal->comp_mcmc_step[component]*gal->comp_scale_length[component];
			step_z 				= gal->comp_mcmc_step[component]*gal->comp_scale_height[component];
			theta 				= 2.0*pi*gsl_rng_uniform_pos(r[0]);
			i 					= gal->comp_start_part[component];
			gal->x[i] 			= max(gal->comp_cut_in[component],step_r)*cos(theta);
			gal->y[i] 			= max(gal->comp_cut_in[component],step_r)*sin(theta);
			gal->z[i] 			= 0.;
			gal->r_cyl[i]		= sqrt(gal->x[i]*gal->x[i]+gal->y[i]*gal->y[i]);
			gal->theta_cyl[i]	= atan2(gal->y[i],gal->x[i]);
			gal->r_sph[i]		= sqrt(gal->x[i]*gal->x[i]+gal->y[i]*gal->y[i]+gal->z[i]*gal->z[i]);
			gal->theta_sph[i]	= atan2(gal->y[i],gal->x[i]);
			gal->phi_sph[i]		= acos(gal->z[i]/gal->r_sph[i]);

			acceptance 			= 0.;
			// Burning period
			for (k = 1; k<(int)(0.1*gal->comp_npart_pot[component]); ++k) {
				prev_step_z = step_z;
				step_z 		= gal->comp_mcmc_step_hydro[component]*sqrt(pow(gal->comp_cs_init[component],2.0)/(2.0*pi*G*pseudo_density_gas_func(gal,fabs(gal->r_cyl[i-1]),gal->theta_cyl[i-1],gal->z[i-1],0,density_model,component)*unit_dens))/kpc;
				// Generating a proposal
				prop_r 		= gal->r_cyl[i] + gsl_ran_gaussian(r[0],step_r);
				prop_theta 	= 2.0*pi*gsl_rng_uniform_pos(r[0]);
				prop_z 		= gal->z[i] + gsl_ran_gaussian(r[0],step_z);
				// Distribution function of the considered component
				// Density of the component times the integrand of spherical volume element
				pi_x = fabs(gal->r_cyl[i])*pseudo_density_gas_func(gal,fabs(gal->r_cyl[i]),gal->theta_cyl[i],gal->z[i],1,density_model,component);
				pi_y = fabs(prop_r)*pseudo_density_gas_func(gal,fabs(prop_r),prop_theta,prop_z,1,density_model,component);
				q_x  = gsl_ran_gaussian_pdf(prop_r-gal->r_cyl[i],step_r)*gsl_ran_gaussian_pdf(prop_z-gal->z[i],prev_step_z);
				q_y  = gsl_ran_gaussian_pdf(gal->r_cyl[i]-prop_r,step_r)*gsl_ran_gaussian_pdf(gal->z[i]-prop_z,step_z);
				prob = min(1.0,(pi_y/pi_x)*(q_x/q_y));
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
	
			printf(" / Computing coords ");
			fflush(stdout);
			// Filling the Markov Chain
			for(i = gal->comp_start_part[component]+1; i < gal->comp_start_part[component]+gal->comp_npart_pot[component]; ++i) {
				prev_step_z = step_z;
				step_z 		= gal->comp_mcmc_step_hydro[component]*sqrt(pow(gal->comp_cs_init[component],2.0)/(2.0*pi*G*pseudo_density_gas_func(gal,fabs(gal->r_cyl[i-1]),gal->theta_cyl[i-1],gal->z[i-1],0,density_model,component)*unit_dens))/kpc;
				// Security against NaN values
				if(step_z != step_z) step_z = prev_step_z;
				// Generating a proposal
				prop_r 		= gal->r_cyl[i-1] + gsl_ran_gaussian(r[0],step_r);
				prop_theta 	= 2.0*pi*gsl_rng_uniform_pos(r[0]);
				prop_z 		= gal->z[i-1] + gsl_ran_gaussian(r[0],step_z);
				// Distribution function of the considered component
				// Density of the component times the integrand of spherical volume element
				pi_x = fabs(gal->r_cyl[i-1])*pseudo_density_gas_func(gal,fabs(gal->r_cyl[i-1]),gal->theta_cyl[i-1],gal->z[i-1],1,density_model,component);
				pi_y = fabs(prop_r)*pseudo_density_gas_func(gal,fabs(prop_r),prop_theta,prop_z,1,density_model,component);
				q_x  = gsl_ran_gaussian_pdf(prop_r-gal->r_cyl[i-1],step_r)*gsl_ran_gaussian_pdf(prop_z-gal->z[i-1],prev_step_z);
				q_y  = gsl_ran_gaussian_pdf(gal->r_cyl[i-1]-prop_r,step_r)*gsl_ran_gaussian_pdf(gal->z[i-1]-prop_z,step_z);
				if((j==n_iter-1)&&(fabs(prop_z)>gal->comp_scale_height[component])) pi_y = 0.;
				prob = min(1.0,(pi_y/pi_x)*(q_x/q_y));
				randval = gsl_rng_uniform_pos(r[0]);
				// Proposal accepted
				if (randval <= prob) {
					gal->r_cyl[i] 		= prop_r;
					gal->theta_cyl[i] 	= prop_theta;
					gal->z[i] 			= prop_z;
					gal->rho[i]			= pi_y/fabs(prop_r);
					acceptance+=1.0;
				// Proposal rejected, the particle keeps the same postion
				} else {
					gal->r_cyl[i] 		= gal->r_cyl[i-1];
					gal->theta_cyl[i] 	= gal->theta_cyl[i-1];
					gal->z[i] 			= gal->z[i-1];
					gal->rho[i]			= gal->rho[i-1];
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
			printf("[acceptance=%.2lf]\n",acceptance);
			if(acceptance<0.80) {
				gal->comp_mcmc_step_hydro[component]/=2.0;
				printf("/////\t\t\t[Warning] MCMC acceptance is low -> Decreasing mcmc_step_hydro%d to %.3lf\n",component+1,gal->comp_mcmc_step_hydro[component]);
			}
			if(acceptance>0.95) {
				gal->comp_mcmc_step_hydro[component]*=2.0;
				printf("/////\t\t\t[Warning] MCMC acceptance is high -> Increasing mcmc_step_hydro%d to %.3lf\n",component+1,gal->comp_mcmc_step_hydro[component]);
			}
			if((acceptance<0.10)&&((j+1)==n_iter)) {
				printf("[Error] MCMC failed to converge\n");
				exit(0);
			}
		}
	}
	return 0;
}

void mcmc_metropolis_hasting_stream(stream *st, int component, int density_model) {
	unsigned long int i,j;
	double pi_x, pi_y, q_x, q_y, prob, *radius;
	double theta, phi, randval, proposal, prop_r, prop_theta, prop_z, step_r, step_z, new_step_r, new_step_z;
	double acceptance;

	if(st->comp_npart[component]>0) {
		// Use the Metropolis algorithm to place the disk particles.
		// We start the Monte Carlo Markov Chain with a realistic particle position
        theta 				= 2.0*pi*gsl_rng_uniform_pos(r[0]);
		i 					= st->comp_start_part[component];
        st->x[i] 			= 0.;
        st->y[i] 			= 0.;
		st->z[i] 			= 0.;
		st->r_cyl[i]		= sqrt(st->x[i]*st->x[i]+st->y[i]*st->y[i]);
		st->theta_cyl[i]	= atan2(st->y[i],st->x[i]);
		st->r_sph[i]		= sqrt(st->x[i]*st->x[i]+st->y[i]*st->y[i]+st->z[i]*st->z[i]);
		st->theta_sph[i]	= atan2(st->y[i],st->x[i]);
		st->phi_sph[i]		= acos(st->z[i]/st->r_sph[i]);
		step_r 				= st->comp_mcmc_step[component]*st->comp_length[component];
		step_z 				= st->comp_mcmc_step[component]*st->comp_length[component];

		acceptance 			= 0.;
		// Burning period
        for (j = 1; j<(int)(0.1*st->comp_npart[component]); ++j) {
			// Generating a proposal
			prop_r 		= st->r_cyl[i] + gsl_ran_gaussian(r[0],step_r);
			prop_theta 	= 2.0*pi*gsl_rng_uniform_pos(r[0]);
			prop_z 		= st->z[i] + gsl_ran_gaussian(r[0],step_z);
			// Distribution function of the considered component
			// Density of the component times the integrand of spherical volume element
			pi_x = fabs(st->r_cyl[i])*density_functions_stream_pool(st,fabs(st->r_cyl[i]),st->theta_cyl[i],st->z[i],density_model,component);
			pi_y = fabs(prop_r)*density_functions_stream_pool(st,fabs(prop_r),prop_theta,prop_z,density_model,component);
			q_x  = gsl_ran_gaussian_pdf(prop_r-st->r_cyl[i],step_r)*gsl_ran_gaussian_pdf(prop_z-st->z[i],step_z);
			q_y  = gsl_ran_gaussian_pdf(st->r_cyl[i]-prop_r,step_r)*gsl_ran_gaussian_pdf(st->z[i]-prop_z,step_z);
			prob = min(1.0,(pi_y/pi_x));//*(q_x/q_y));
			randval = gsl_rng_uniform_pos(r[0]);
			if (randval <= prob) {
				st->r_cyl[i] 		= prop_r;
				st->theta_cyl[i] 	= prop_theta;
				st->z[i] 			= prop_z;
			}
		}
		// Updating the coordinate values
		st->x[i] = st->r_cyl[i]*cos(st->theta_cyl[i]);
		st->y[i] = st->r_cyl[i]*sin(st->theta_cyl[i]);
		// Updating the coordinate values
        st->theta_cyl[i]	= atan2(st->y[i],st->x[i]);
        st->r_sph[i]		= sqrt(st->x[i]*st->x[i]+st->y[i]*st->y[i]+st->z[i]*st->z[i]);
        st->theta_sph[i]	= atan2(st->y[i],st->x[i]);
        st->phi_sph[i]		= acos(st->z[i]/st->r_sph[i]);
		// Filling the Markov Chain
        for (i = st->comp_start_part[component]+1; i<st->comp_start_part[component]+st->comp_npart[component]; ++i) {
			// Generating a proposal
			prop_r 		= st->r_cyl[i-1] + gsl_ran_gaussian(r[0],step_r);
			prop_theta 	= 2.0*pi*gsl_rng_uniform_pos(r[0]);
			prop_z 		= st->z[i-1] + gsl_ran_gaussian(r[0],step_z);
			// Distribution function of the considered component
			// Density of the component times the integrand of spherical volume element
			pi_x = fabs(st->r_cyl[i-1])*density_functions_stream_pool(st,fabs(st->r_cyl[i-1]),st->theta_cyl[i-1],st->z[i-1],density_model,component);
			pi_y = fabs(prop_r)*density_functions_stream_pool(st,fabs(prop_r),prop_theta,prop_z,density_model,component);
			q_x  = gsl_ran_gaussian_pdf(prop_r-st->r_cyl[i-1],step_r)*gsl_ran_gaussian_pdf(prop_z-st->z[i-1],step_z);
			q_y  = gsl_ran_gaussian_pdf(st->r_cyl[i-1]-prop_r,step_r)*gsl_ran_gaussian_pdf(st->z[i-1]-prop_z,step_z);
			prob = min(1.0,(pi_y/pi_x));//*(q_x/q_y));
			randval = gsl_rng_uniform_pos(r[0]);
			if (randval <= prob) {
				st->r_cyl[i] 		= prop_r;
				st->theta_cyl[i] 	= prop_theta;
				st->z[i] 			= prop_z;
				st->rho[i]			= pi_y/fabs(st->r_cyl[i]);
				acceptance+=1.0;
			} else {
				st->r_cyl[i] 		= st->r_cyl[i-1];
				st->theta_cyl[i] 	= st->theta_cyl[i-1];
				st->z[i] 			= st->z[i-1];
				st->rho[i]			= pi_x/fabs(st->r_cyl[i]);
			}
			// Updating the coordinate values
			st->x[i] = st->r_cyl[i]*cos(st->theta_cyl[i]);
			st->y[i] = st->r_cyl[i]*sin(st->theta_cyl[i]);
			// Updating the coordinate values
            st->theta_cyl[i]	= atan2(st->y[i],st->x[i]);
            st->r_sph[i]		= sqrt(st->x[i]*st->x[i]+st->y[i]*st->y[i]+st->z[i]*st->z[i]);
            st->theta_sph[i]	= atan2(st->y[i],st->x[i]);
            st->phi_sph[i]		= acos(st->z[i]/st->r_sph[i]);
		}
		acceptance /= st->comp_npart[component];
		printf("[acceptance=%.2lf]\n",acceptance);
		if(acceptance<0.50) printf("/////\t\t[Warning] MCMC acceptance is low -> Decrease mcmc_step%d\n",component+1);
		if(acceptance>0.90) printf("/////\t\t[Warning] MCMC acceptance is high -> Increase mcmc_step%d\n",component+1);
	}
	return;
}


// Surface density function
// The surface density is computed using a numerical integration
// of the density profile over z
double surface_density_func(galaxy *gal, double r, double theta, int cut, int component) {
	
	int status,tid;
	double surface_density,error,h;
	
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(AllVars.GslWorkspaceSize);
	gsl_function F;
	
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
	
	F.function = &integrand_density_func;
    F.params = gal;
	gal->storage[0][tid] 	= r;
	gal->storage[1][tid] 	= theta;
	gal->storage[2][tid] 	= cut;
	gal->selected_comp[tid] = component;

	gsl_integration_qag(&F,-gal->comp_cut[component]*gal->comp_flat[component],gal->comp_cut[component]*gal->comp_flat[component],epsabs,epsrel,AllVars.GslWorkspaceSize,key,w,&surface_density,&error);
    
	gsl_integration_workspace_free(w);

	return surface_density;
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
	
	r_temp 		= gal->storage[0][tid];
	theta_temp 	= gal->storage[1][tid];
	cut 		= gal->storage[2][tid];
	
	// This function should not be used with the pseudo-density function
	rho = density_functions_pool(gal,r_temp,theta_temp,z,cut,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);

	return rho;
}

// This function calculates the force on a test particle due to the halo
// along the r axis. This is used to calculate the velocity dispersion.
double cumulative_mass_func(galaxy *gal, double radius, int component) {
	
	int status;
	double result,integral,error;
	int tid;
	size_t neval;
	
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
	
	gsl_function F;
	
	F.function 	= &d_cumulative_mass_func1;
	F.params 	= gal;
		
	gal->selected_comp[tid] = component; 
	
	gsl_integration_qng(&F,0.0,radius,epsabs,epsrel,&integral,&error,&neval);
	
	result = integral;

	return result;
}

// This is the integrand for the previous function. It is setup to work with
// the GSL_qags structures.
static double d_cumulative_mass_func1(double r, void *params) {
	
	int tid;
	double result,integrand,error;
	size_t neval;
	
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
	
	gsl_function F;
	
	galaxy *gal = (galaxy *) params;
	
	F.function 	= &d_cumulative_mass_func2;
	F.params 	= gal;
	
	gal->storage[6][tid] = r;
	
	gsl_integration_qng(&F,0.0,2.0*pi,epsabs,epsrel,&integrand,&error,&neval);

	return integrand;
}

// This is the integrand for the previous function. It is setup to work with
// the GSL_qags structures.
static double d_cumulative_mass_func2(double theta, void *params) {
	int tid;
	double surface_density,r;
	
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif

	galaxy *gal = (galaxy *) params;
	
	r = gal->storage[6][tid];
	
	surface_density = r*surface_density_func(gal,r,theta,1,gal->selected_comp[tid]);

	return surface_density;
}


// Surface density function
// The surface density is computed using a numerical integration
// of the density profile over z
double surface_density_func_stream(stream *st, double r, double theta, int component) {
	
	int status,tid;
	double surface_density,error,h;
	
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(AllVars.GslWorkspaceSize);
	gsl_function F;
	
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
	
	F.function = &integrand_density_func_stream;
    F.params = st;
	
	st->storage[0][tid] 	= r;
	st->storage[1][tid] 	= theta;
	
	st->selected_comp[tid] 	= component; 
	gsl_integration_qag(&F,0.,st->comp_length[component],epsabs,epsrel,AllVars.GslWorkspaceSize,key,w,&surface_density,&error);
    
	gsl_integration_workspace_free(w);
	return surface_density;
}

// Integrand of the previous numerical integration
static double integrand_density_func_stream(double z, void *params) {
	
	int tid;
	double r_temp,theta_temp,rho;
    
	stream *st = (stream *) params;
	
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
	
	r_temp 		= st->storage[0][tid];
	theta_temp 	= st->storage[1][tid];
	
	rho = density_functions_stream_pool(st,r_temp,theta_temp,z,st->comp_model[st->selected_comp[tid]],st->selected_comp[tid]);
	
	return rho;
}

// This function calculates the force on a test particle due to the halo
// along the r axis. This is used to calculate the velocity dispersion.
double cumulative_mass_func_stream(stream *st, double radius, int component) {
	
	int status;
	double result,integral, error;
	int tid;
	
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
	
	gsl_integration_workspace *wk = gsl_integration_workspace_alloc(AllVars.GslWorkspaceSize);
	gsl_function F;
	
	F.function 	= &d_cumulative_mass_func_stream;
	F.params 	= st;
		
	st->selected_comp[tid] = component; 
	
	gsl_integration_qag(&F,0,radius,epsabs,epsrel,AllVars.GslWorkspaceSize,key,wk,&integral,&error);
	gsl_integration_workspace_free(wk);
	
	result = 2*pi*integral;
	
	return result;
}

// This is the integrand for the previous function. It is setup to work with
// the GSL_qags structures.
static double d_cumulative_mass_func_stream(double r, void *params) {
	
	int tid;
	double integrand;
	
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
	
	stream *st = (stream *) params;
	
	integrand = surface_density_func_stream(st,r,0.,st->selected_comp[tid])*r;
    
	return integrand;
}


// This function computes the density of the gas assuming vertical structure hydrostatic equilibirum.
// Using this function into a iterative algorithm allows to converge towards
// a complete equilibrium state
double pseudo_density_gas_func(galaxy *gal, double r, double theta, double z, int cut, int density_model, int component) {
	
	int tid;
	double density, delta_pot, rho_0, save1, save2;

	#if USE_THREADS == 1
		tid = omp_get_thread_num();
	#else
		tid = 0;
	#endif

	// Computing the potential gradient between 0 and z
	// in cgs unit	
	save1 							= gal->x[gal->index[tid]];
	save2 							= gal->y[gal->index[tid]];
	gal->x[gal->index[tid]] 		= r*cos(theta);
	gal->y[gal->index[tid]] 		= r*sin(theta);
	delta_pot 						= galaxyz_potential_wrapper_func(z,gal)-galaxyz_potential_wrapper_func(0.,gal);
	// Density in the xy-plane in 1e10 solar mass / kpc^3
	//rho_0 							= density_functions_pool(gal,r,theta,0.,1,density_model,component);
	rho_0 							= get_midplane_density(gal,gal->x[gal->index[tid]],gal->y[gal->index[tid]]);
	gal->x[gal->index[tid]] 		= save1;
	gal->y[gal->index[tid]] 		= save2;

	// Hydrostatic equilibrium requires following density
	density 						= rho_0*exp(-delta_pot/(pow(gal->comp_cs_init[component],2.0)));
	
	if(cut && r>gal->comp_cut[component]) density = 0.;
	
	// Return a pseudo density 
	return density;
}

double midplane_density_gas_func(galaxy *gal, gsl_integration_workspace *w, double x, double y, int component) {
    int status,pseudo_save;
    int tid;
    double result, integral, error, initial_surface_density, radius, theta;
    
    gsl_function F;
    
    #if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
	
	gal->storage[3][tid] 	= x;
	gal->storage[4][tid] 	= y;
	
	radius 					= sqrt(x*x+y*y);
	theta					= atan2(y,x);
	pseudo_save 			= gal->pseudo[tid];
	gal->pseudo[tid] 		= 0;
	initial_surface_density = surface_density_func(gal,radius,theta,0,component);
	gal->pseudo[tid] 		= pseudo_save;

	if(initial_surface_density==0.) return 0.;
    
    F.function 	= &dmidplane_density_gas_func;
    F.params 	= gal;
    
    gsl_integration_qag(&F,-10.*gal->comp_scale_height[component],10.*gal->comp_scale_height[component],epsabs,epsrel,AllVars.GslWorkspaceSize,key,w,&integral,&error);

	return initial_surface_density/integral;
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

	x = gal->storage[3][tid];
	y = gal->storage[4][tid];
    
	delta_pot = galaxy_potential_func(gal,gal->potential,gal->dx,gal->ngrid,x,y,z,1)-galaxy_potential_func(gal,gal->potential,gal->dx,gal->ngrid,x,y,0.,1);
    
	integrand = exp(-delta_pot/pow(gal->comp_cs_init[gal->selected_comp[tid]],2.0));

	return integrand;
}

// This function fills a 2D grid with the value of the gas density in the midplane
// The use of a 2D grid intends to lower the computation time
void fill_midplane_dens_grid(galaxy *gal, int component) {
	int i,j;
	unsigned long int k;
	double x,y;
	
	gsl_integration_workspace **wk;
	
	wk = (gsl_integration_workspace **) malloc(AllVars.Nthreads * sizeof(gsl_integration_workspace *));
	for(i=0;i<AllVars.Nthreads;i++) {
		wk[i] = gsl_integration_workspace_alloc(AllVars.GslWorkspaceSize);
	}
	
	#pragma omp parallel shared(gal) private(x,y,i,j)
	for (i=0;i<gal->ngrid_dens[0];i++) {	
		for (j=0;j<gal->ngrid_dens[1];j++) {
			// Get thread ID
			#if USE_THREADS == 1
            	int tid = omp_get_thread_num();
			#else
            	int tid = 0;
			#endif
			x = ((double)i-((double)(gal->ngrid_dens[0]/2)-0.5))*gal->dx_dens;
			y = ((double)j-((double)(gal->ngrid_dens[1]/2)-0.5))*gal->dx_dens;
			gal->midplane_dens[i][j] = midplane_density_gas_func(gal,wk[tid],x,y,component);
		}
	}
	for(i=0;i<AllVars.Nthreads;i++) {
		gsl_integration_workspace_free(wk[i]);
	}
	return;
}

//
double get_midplane_density(galaxy *gal, double x, double y) {
	int node_x,node_y;
	double xtemp,ytemp,dens1,dens2,dens3,dens4,dens_interp;
	double dx,dy,tx,ty;
    
	xtemp = x/gal->dx_dens + ((double)(gal->ngrid_dens[0]/2)-0.5);
	ytemp = y/gal->dx_dens + ((double)(gal->ngrid_dens[1]/2)-0.5);
	// Determine the parent node.
	node_x = floor(xtemp);
	node_y = floor(ytemp);
		
	if(node_x>=0 && node_x<gal->ngrid_dens[0]-1 && node_y>=0 && node_y<gal->ngrid_dens[1]-1) {    
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
	} else {
		dens_interp = 0.;
	}
	return dens_interp;
	
}

// ------------------------------------------

// This function determines the scale length of a self gravitating disk in a NFW halo using
// the fitting formula provided in Mo et al. 1998, if the user want to use it.
// Otherwise, the user chose himself the value for the disk scale length.
double disk_scale_length_func(galaxy *gal, double c) {
	
	int i;
	double f_r_base, f_r_power, f_c, f_r, disk_scale;
	
	f_c = 2.0/3.0+ pow((c/21.5),0.7);
	
	f_r_base = (gal->j_d*gal->lambda)/(0.1*gal->m_d);
	f_r_power = (-0.06+2.71*gal->m_d+0.0047*gal->m_d/(gal->j_d*gal->lambda));
	f_r = pow(f_r_base,f_r_power)*(1.0-3.0*gal->m_d+5.2*gal->m_d*gal->m_d)*(1.0-0.019*c+0.00025*c*c+0.52/c);
	
	disk_scale = (1.0/sqrt(2.0))*(gal->j_d/gal->m_d)*gal->lambda*gal->r200*(f_r/sqrt(f_c));

	return disk_scale;
}

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
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(AllVars.GslWorkspaceSize);
	gsl_function F;
	
	F.function = &dg_c_func;
	gsl_integration_qag(&F,0.0,c,epsabs,epsrel,AllVars.GslWorkspaceSize,key,w,&result,&error);
	
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

	printf("/////\tLowering particule mass resolution\n");
	for (j=0; j<AllVars.MaxCompNumber; j++) {
    	if(gal->comp_npart[j]>0) printf("/////\t\t- Component %2d -> m=%.2e Msol\n",j+1,(gal->comp_cutted_mass[j]*1.0E10)/(gal->comp_npart[j]));
		// Filling the arrays of the &galaxy structure
		for (i = gal->comp_start_part[j]; i < gal->comp_start_part[j] + gal->comp_npart_pot[j]; i++) {
			gal->mass[i] 		= gal->comp_cutted_mass[j]/gal->comp_npart[j];
		}
	}
	return;
}
