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
  / Date: September 2015
  /
  */

#include "dice.h"


double density_functions_pool(galaxy *gal, double radius, double theta, double z, int cut, int model, int component) {
    double h, hx, hy, hz, h_cut, hx_cut, hy_cut, hz_cut, hx_cut_in, hy_cut_in, hz_cut_in;
    double density, w, k, l, m, n, o, r, s, alpha, beta, smooth_factor1, smooth_factor2, sigma1, sigma2, r_sph;
    double x, y, flatx, flaty, flatz, x2, y2, z2;

    if(fabs(radius)>0.) {
	radius = (fabs(radius)/radius)*max(gal->comp_rcore[component],fabs(radius));
    } else {
 	radius = max(gal->comp_rcore[component],fabs(radius));
    }
    double theta_shift, theta_out, A, B, CEDF, tanh_func;
    double z_shift;

    // The spiral galaxy model is inspired from  the GALFIT software (Peng et al. 2010)
    theta_shift = 0.;
    if(gal->comp_spiral_theta_out[component]>0.) {
        theta_out = gal->comp_spiral_theta_out[component]*pi/180.;
        CEDF = 0.23;
        A = (2*CEDF)/(fabs(theta_out)+CEDF)-1.00001;
        B = (2-1./tanh(A))*(gal->comp_spiral_r_out[component]/(gal->comp_spiral_r_out[component]-gal->comp_spiral_r_in[component]));
        tanh_func = 0.5*(tanh(B*(radius/gal->comp_spiral_r_out[component]-1)+2)+1);
        theta_shift = theta_out*tanh_func*pow(0.5*(radius/gal->comp_spiral_r_out[component]+1),gal->comp_spiral_alpha[component]);
    }
    // Warping
    if(gal->comp_warp_scale[component]>0.) {
        z_shift = gal->comp_warp_scale[component]*gal->comp_scale_length[component]*gal->comp_flatz[component]*(cos(gal->comp_warp_mode[component]*theta)-0.5)*radius/gal->comp_cut[component];
        z += z_shift;
    }

    x = radius*cos(theta+theta_shift);
    y = radius*sin(theta+theta_shift);
    r_sph = sqrt(x*x+y*y+z*z);

    h_cut = gal->comp_cut[component];
    hx_cut = gal->comp_cut[component]*gal->comp_flatx_cut[component];
    hy_cut = gal->comp_cut[component]*gal->comp_flaty_cut[component];
    hz_cut = gal->comp_cut[component]*gal->comp_flatz_cut[component];
    hx_cut_in = max(gal->comp_cut_in[component],0.01*gal->comp_scale_length[component])*gal->comp_flatx_cut[component];
    hy_cut_in = max(gal->comp_cut_in[component],0.01*gal->comp_scale_length[component])*gal->comp_flaty_cut[component];
    hz_cut_in = max(gal->comp_cut_in[component],0.01*gal->comp_scale_length[component])*gal->comp_flatz_cut[component];

    flatx = fabs(r_sph)*(gal->comp_flatx_out[component]-gal->comp_flatx[component])/hx_cut+gal->comp_flatx[component];
    flaty = fabs(r_sph)*(gal->comp_flaty_out[component]-gal->comp_flaty[component])/hy_cut+gal->comp_flaty[component];
    flatz = fabs(r_sph)*(gal->comp_flatz_out[component]-gal->comp_flatz[component])/hz_cut+gal->comp_flatz[component];

    h = gal->comp_scale_length[component];
    hx = gal->comp_scale_length[component]*flatx;
    hy = gal->comp_scale_length[component]*flaty;
    hz = gal->comp_scale_length[component]*flatz;

    alpha = gal->comp_alpha[component];
    beta = gal->comp_beta[component];

    k = sqrt(pow(z/hz,2.0));
    l = sqrt(pow(x/hx,2.0)+pow(y/hy,2.0));
    m = sqrt(pow(x/hx,2.0)+pow(y/hy,2.0)+pow(z/hz,2.0));
    n = sqrt(pow(x/hx_cut,2.0)+pow(y/hy_cut,2.0)+pow(z/hz_cut,2.0));
    o = sqrt(pow(x/hx_cut,2.0)+pow(y/hy_cut,2.0));
    r = sqrt(pow(x*flatx,2.0)+pow(y*flaty,2.0));
    if(gal->comp_symmetry[component]==1) {
        s = sqrt(pow(x/hx_cut_in,2.0)+pow(y/hy_cut_in,2.0));
    } else {
        s = sqrt(pow(x/hx_cut_in,2.0)+pow(y/hy_cut_in,2.0)+pow(z/hz_cut_in,2.0));
    }

    // Select a disk model
    switch(model) {
        case 1:
            // Exponential disk + sech z-profile
            if(strcmp(gal->comp_profile_name[component],"")==0)
                strcpy(gal->comp_profile_name[component]," Exponential disk/sech z ");
            density = gal->comp_scale_dens[component]*exp(-l)/(pow(cosh(z/hz),2));
            break;
        case 2:
            // Myamoto-Nagai profile
            if(strcmp(gal->comp_profile_name[component],"")==0)
                strcpy(gal->comp_profile_name[component],"      Myamoto-Nagai      ");
            w = sqrt(z*z+hz*hz);
            density = gal->comp_scale_dens[component]*(h*r*r+(h+3.0*w)*pow(h+w,2.0))
                /(pow(r*r+pow((h+w),2.0),2.5)*pow(w,3.0));
            break;
        case 3:
            // Exponential disk + exponential z-profile
            if(strcmp(gal->comp_profile_name[component],"")==0)
                strcpy(gal->comp_profile_name[component]," Exponential disk/exp z  ");
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
        case 13:
            // Sersic profile
            if(strcmp(gal->comp_profile_name[component],"")==0)
                strcpy(gal->comp_profile_name[component],"         Sersic          ");
            density = gal->comp_scale_dens[component]*exp(-beta*pow(m,1/alpha)-1);
            break;
        case 14:
            // Toomre-Kuzmin profile
            if(strcmp(gal->comp_profile_name[component],"")==0)
                strcpy(gal->comp_profile_name[component],"      Toomre-Kuzmin      ");
            density = gal->comp_scale_dens[component]/(pow(l*l+1,3./2.)*pow(cosh(z/hz),2));
            break;
        case 15:
            // Uniform profile
            if(strcmp(gal->comp_profile_name[component],"")==0)
                strcpy(gal->comp_profile_name[component],"         Uniform         ");
            density = gal->comp_scale_dens[component];
            break;
        case 16:
            // Pseudo-isothermal profile
            if(strcmp(gal->comp_profile_name[component],"")==0)
                strcpy(gal->comp_profile_name[component],"     Pseudo-Isothermal   ");
            density = gal->comp_scale_dens[component]*1.0/(1.0+pow(m,2.0));
            break;

        default:
            fprintf(stderr,"[Error] model%d=%d is not a valid value\n",component+1,model);
            exit(0);
    }
    sigma1 = gal->comp_sigma_cut[component];
    sigma2 = gal->comp_sigma_cut_in[component];
    smooth_factor1 = 1-0.5*(1+erf((n-1.0)/(sigma1*sqrt(2))));
    smooth_factor2 = 0.5*(1+erf((s-1.0)/(sigma2*sqrt(2))));

    if(gal->dens_gauss_sigma>0. && gal->comp_dens_gauss[component]==1 && gal->gaussian_field_defined) {
        density *= (galaxy_gaussian_field_func(gal,x,y,z)*gal->dens_gauss_sigma+1.0);
        if(density<0.) density = 0.0;
    }
    if(gal->comp_excavate[component]>0 && gal->comp_npart[gal->comp_excavate[component]]>0) {
        double dens_target = density_functions_pool(gal,r,theta,z,1,gal->comp_model[gal->comp_excavate[component]-1],gal->comp_excavate[component]-1);
	if(dens_target>density) density=0.;
    }
    // Cutting the density
    if(cut==1) {
        density *= smooth_factor1;
        if(density<gal->comp_cut_dens[component]/unit_nh) density = 0.;
        if(gal->comp_cut_in[component]>0.) density *= smooth_factor2;
    }

    return density;
}

double density_functions_stream_pool(stream *st, double radius, double theta, double z, int model, int component) {
    double alpha, h, rs, density, x, y;
    // We consider only positive values
    alpha = pi/180.*st->comp_opening_angle[component];
    h = st->comp_length[component];
    rs = st->comp_scale[component]*fabs(z)*atan(pi/180.*st->comp_opening_angle[component]/2.);

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
    if(st->dens_gauss_sigma>0. && st->comp_dens_gauss[component]==1 && st->gaussian_field_defined) {
        density *= (stream_gaussian_field_func(st,x,y,z)*st->dens_gauss_sigma+1.0);
        if(density<0.) density = 0.0;
    }
    return density;
}

void mcmc_metropolis_hasting_ntry(galaxy *gal, int component, int density_model) {
    unsigned long int i,j,start_part,npart;
    int k, selected, tid;
    double prob, *radius, k_poly, d, rc;
    double theta, phi, randval, smooth_factor;
    double step_r, step_x, step_y, step_z, step_r_sph, hx, hy, hz;
    double new_step_x, new_step_y, new_step_z, new_step_r, new_step_r_sph, min_step_z, max_step_z;
    double acceptance, norm, ratio, rho_ref, rho_0, mean_metal, Tpart, Tmax, cs2;
    double *prop_x, *prop_y, *prop_z, *prop_r, *prop_theta, *prop_r_sph, *prop_phi_sph;
    double *pi_x, *pi_y, *q_x, *q_y, *weights, *w_x, *w_y;
    double *ref_x, *ref_y, *ref_z, *ref_r, *ref_theta, *ref_r_sph, *ref_phi_sph, *lambda_x, *lambda_y, *dv_x, *dv_y;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    if(!(prop_x = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate prop_x array\n");
        exit(0);
    }
    if(!(prop_y = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate prop_y array\n");
        exit(0);
    }
    if(!(prop_z = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate prop_z array\n");
        exit(0);
    }
    if(!(prop_r = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate prop_r array\n");
        exit(0);
    }
    if(!(prop_theta = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate prop_theta array\n");
        exit(0);
    }
    if(!(prop_r_sph = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate prop_r_sph array\n");
        exit(0);
    }
    if(!(prop_phi_sph = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate prop_phi_sph array\n");
        exit(0);
    }
    if(!(pi_x = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate pi_x array\n");
        exit(0);
    }
    if(!(pi_y = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate pi_y array\n");
        exit(0);
    }
    if(!(q_x = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate q_x array\n");
        exit(0);
    }
    if(!(q_y = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate q_y array\n");
        exit(0);
    }
    if(!(weights = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate w array\n");
        exit(0);
    }
    if(!(w_x = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate w_x array\n");
        exit(0);
    }
    if(!(w_y = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate w_y array\n");
        exit(0);
    }
    if(!(ref_x = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate ref_x array\n");
        exit(0);
    }
    if(!(ref_y = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate ref_y array\n");
        exit(0);
    }
    if(!(ref_z = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate ref_z array\n");
        exit(0);
    }
    if(!(ref_r = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate ref_r array\n");
        exit(0);
    }
    if(!(ref_theta = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate ref_theta array\n");
        exit(0);
    }
    if(!(ref_r_sph = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate ref_r_sph array\n");
        exit(0);
    }
    if(!(ref_phi_sph = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate ref_phi_sph array\n");
        exit(0);
    }
    if(!(lambda_x = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate lambda_x array\n");
        exit(0);
    }
    if(!(lambda_y = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate lambda_y array\n");
        exit(0);
    }
    if(!(dv_x = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate dv_x array\n");
        exit(0);
    }
    if(!(dv_y = calloc(gal->mcmc_ntry,sizeof(double)))) {
        fprintf(stderr,"[Error] Unable to allocate dv_y array\n");
        exit(0);
    }

    rc = gal->comp_rc_entropy[component];
    hx = gal->comp_mcmc_step[component]*gal->comp_scale_length[component]*gal->comp_flatx_cut[component];
    hy = gal->comp_mcmc_step[component]*gal->comp_scale_length[component]*gal->comp_flaty_cut[component];
    hz = gal->comp_mcmc_step[component]*gal->comp_scale_length[component]*gal->comp_flatz_cut[component];

    if(gal->pseudo[0]) {
        min_step_z = 1e-3*gal->comp_cut[component]*gal->comp_mcmc_step_hydro[component];
        max_step_z = gal->comp_cut[component]*gal->comp_mcmc_step_hydro[component];
    }

    i = gal->comp_start_part[component];
    if(gal->comp_npart[component]>1) {
        // Use the Metropolis algorithm to place the disk particles.
        // We start the Monte Carlo Markov Chain with a realistic particle position
        step_x = hx;
        step_y = hy;
        step_z = hz;
        if(gal->pseudo[0]) {
            step_z = gal->comp_mcmc_step_hydro[component]*gal->comp_scale_length[component]*gal->comp_flatz_cut[component];
            if(isinf(step_z)||step_z<min_step_z) step_z = min_step_z;
            if(step_z>max_step_z) step_z = max_step_z;
            if(gal->comp_spherical_hydro_eq[component]) {
                step_x = step_z;
                step_y = step_z;
            }
        }
        step_r = sqrt(pow(step_x,2)+pow(step_y,2));
        step_r_sph = sqrt(pow(step_x,2)+pow(step_y,2)+pow(step_z,2));
        theta = 2.0*pi*gsl_rng_uniform_pos(r[0]);
        // Single particle always placed at the center
        if(gal->comp_npart[component]==1) {
            gal->x[i] = 0.0;
            gal->y[i] = 0.0;
        } else {
            gal->x[i] = (1.01*gal->comp_cut_in[component]+0.01*gal->comp_cut[component])*cos(theta);
            gal->y[i] = (1.01*gal->comp_cut_in[component]+0.01*gal->comp_cut[component])*sin(theta);
        }
        gal->z[i] = 0.;
        gal->r_cyl[i] = sqrt(gal->x[i]*gal->x[i]+gal->y[i]*gal->y[i]);
        gal->theta_cyl[i] = atan2(gal->y[i],gal->x[i])+pi;
        gal->r_sph[i] = sqrt(gal->x[i]*gal->x[i]+gal->y[i]*gal->y[i]+gal->z[i]*gal->z[i]);
	if(gal->r_sph[i]>0.) {
            gal->phi_sph[i] = acos(gal->z[i]/gal->r_sph[i]);
	} else {
            gal->phi_sph[i] = 0.;
	}
	gal->index[tid] = i;
	// Initialize min max densities
	if(gal->pseudo[0]) {
	    gal->comp_dens_min[component] = pseudo_density_gas_func(gal,fabs(gal->x[i]),0.,0.,1,density_model,component,gal->comp_spherical_hydro_eq[component]);
	} else {
	    gal->comp_dens_min[component] = density_functions_pool(gal,fabs(gal->x[i]),0.,0.,1,density_model,component);
	}
	gal->comp_dens_max[component] = 0.0;
	Tmax = 0.;
	
	if(gal->comp_npart[component]==1) return;

        acceptance = 0.;
        mean_metal = 0.;
        // Reduce the dimensionality of the MCMC if the component is axisymmetric
        gal->comp_symmetry[component] = 0;
        if(gal->comp_flatx[component] == gal->comp_flaty[component]) {
            // Shperical symmetry
            if(gal->comp_flatz[component]==gal->comp_flatx[component] && gal->comp_flatz[component]==gal->comp_flaty[component]) {
                gal->comp_symmetry[component] = 2;
                // Cylindrical symmetry
            } else {
                gal->comp_symmetry[component] = 1;
            }
        }
        // Accept cylindrical symmetry for the gas disk with vertical hydrostatic equilibrium
        // only if all the components display a cylindrical symmetry
        if(gal->pseudo[0] && gal->comp_flatx[component] == gal->comp_flaty[component] && gal->comp_spherical_hydro_eq[component]==0) {
            for(k = 0; k < AllVars.MaxCompNumber; k++) {
                if(gal->comp_flatx[k] == gal->comp_flaty[k]) {
                    gal->comp_symmetry[component] = 1;
                } else {
                    gal->comp_symmetry[component] = 0;
                }
                if(gal->comp_symmetry[component] == 0) break;
            }
        }
        switch(gal->comp_symmetry[component]) {
            case 1:
                printf("[ cylindrical symmetry ]");
                break;
            case 2:
                printf("[  spherical  symmetry ]");
                break;
            default:
                printf("[      no symmetry     ]");
                break;
        }
        fflush(stdout);
        // Filling the Markov Chain
        for(i = gal->comp_start_part[component]+1; i < gal->comp_start_part[component]+gal->comp_npart_pot[component]; ++i) {
            // Generating a proposal
            for(k = 0; k<gal->mcmc_ntry; k++) {
                step_x = hx+gal->comp_mcmc_step_slope[component]*gal->x[i-1]/gal->comp_scale_length[component];
                step_y = hy+gal->comp_mcmc_step_slope[component]*gal->y[i-1]/gal->comp_scale_length[component];
                step_z = hz+gal->comp_mcmc_step_slope[component]*gal->z[i-1]/gal->comp_scale_length[component];
                if(gal->pseudo[0]) {
		    if(gal->comp_spherical_hydro_eq[component]==1) {
                        step_z = gal->comp_mcmc_step_hydro[component]*gal->comp_scale_length[component];
	            } else {
			rho_0 = pseudo_density_gas_func(gal,gal->r_cyl[i-1],gal->theta_cyl[i-1],gal->z[i-1],0,density_model,component,0);
                        cs2 = gal->comp_k_poly[component]*gal->comp_gamma_poly[component]*pow(rho_0,gal->comp_gamma_poly[component]-1.0);
                        step_z = gal->comp_mcmc_step_hydro[component]*sqrt(cs2/(2.0*pi*G*rho_0));
                    }
                    if(isinf(step_z)||step_z<min_step_z) step_z = min_step_z;
                    if(step_z>max_step_z) step_z = max_step_z;
                    if(gal->comp_spherical_hydro_eq[component]) {
                        step_x = step_z;
                        step_y = step_z;
                    }
                }
                step_r = sqrt(pow(step_x,2)+pow(step_y,2));
                step_r_sph = sqrt(pow(step_x,2)+pow(step_y,2)+pow(step_z,2));
                switch(gal->comp_symmetry[component]) {
                    // Cylindrical symmetry
                    case 1:
                        prop_r[k] = gal->r_cyl[i-1] + gsl_ran_gaussian(r[0],step_r);
                        prop_theta[k] = 2.0*pi*gsl_rng_uniform_pos(r[0]);
                        prop_z[k] = gal->z[i-1] + gsl_ran_gaussian(r[0],step_z);
                        q_y[k] = gsl_ran_gaussian_pdf(-gal->r_cyl[i-1]+prop_r[k],step_r);
                        q_y[k] *= gsl_ran_gaussian_pdf(-gal->z[i-1]+prop_z[k],step_z);
                        dv_y[k] = fabs(prop_r[k]);
                        break;
                    // Spherical symmetry
                    case 2:
                        prop_r_sph[k] = gal->r_sph[i-1] + gsl_ran_gaussian(r[0],step_r);
                        prop_phi_sph[k] = acos(2*gsl_rng_uniform_pos(r[0])-1);
                        prop_theta[k] = 2.0*pi*gsl_rng_uniform_pos(r[0]);
                        prop_r[k] = prop_r_sph[k]*sin(prop_phi_sph[k]);
                        prop_z[k] = prop_r_sph[k]*cos(prop_phi_sph[k]);
                        q_y[k] = gsl_ran_gaussian_pdf(-gal->r_sph[i-1]+prop_r_sph[k],step_r_sph);
                        dv_y[k] = fabs(pow(prop_r_sph[k],2));
                        break;
                    // No symmetry
                    default:
                        prop_x[k] = gal->x[i-1] + gsl_ran_gaussian(r[0],step_x);
                        prop_y[k] = gal->y[i-1] + gsl_ran_gaussian(r[0],step_y);
                        prop_z[k] = gal->z[i-1] + gsl_ran_gaussian(r[0],step_z);
                        prop_r[k] = sqrt(pow(prop_x[k],2)+pow(prop_y[k],2));
                        prop_theta[k] = atan2(prop_y[k],prop_x[k]);
                        q_y[k] = gsl_ran_gaussian_pdf(-gal->x[i-1]+prop_x[k],step_x);
                        q_y[k] *= gsl_ran_gaussian_pdf(-gal->y[i-1]+prop_y[k],step_y);
                        q_y[k] *= gsl_ran_gaussian_pdf(-gal->z[i-1]+prop_z[k],step_z);
                        dv_y[k] = 1.0;
                        break;
                }
                // Distribution function of the considered component
                if(gal->pseudo[0]) {
                    pi_y[k] = dv_y[k]*pseudo_density_gas_func(gal,prop_r[k],prop_theta[k],prop_z[k],1,density_model,component,gal->comp_spherical_hydro_eq[component]);
                } else { 
                    pi_y[k] = dv_y[k]*density_functions_pool(gal,prop_r[k],prop_theta[k],prop_z[k],1,density_model,component);
                }
                weights[k] = pi_y[k];
            }
            // Normalize weights
            norm = sum_dbl(weights,gal->mcmc_ntry);
            for(k = 0; k<gal->mcmc_ntry; k++) if(norm != 0.) weights[k] /= norm;
            // Select a proposal according to its probability
            selected = 0;
            randval = gsl_rng_uniform_pos(r[0]);
            for(k = 0; k<gal->mcmc_ntry; k++) {
                if(randval<weights[k]) {
                    selected = k;
                    break;
                }
                randval -= weights[k];
            }

            // Produce a reference set
            for(k = 0; k<gal->mcmc_ntry; k++) {
                new_step_x = hx+gal->comp_mcmc_step_slope[component]*prop_x[k]/gal->comp_scale_length[component];
                new_step_y = hy+gal->comp_mcmc_step_slope[component]*prop_y[k]/gal->comp_scale_length[component];
                new_step_z = hz+gal->comp_mcmc_step_slope[component]*prop_z[k]/gal->comp_scale_length[component];
                if(gal->pseudo[0]) {
		    if(gal->comp_spherical_hydro_eq[component]==1) {
                        new_step_z = gal->comp_mcmc_step_hydro[component]*gal->comp_scale_length[component];
	            } else {
			rho_0 = pseudo_density_gas_func(gal,prop_r[k],prop_theta[k],prop_z[k],0,density_model,component,0);
                        cs2 = gal->comp_k_poly[component]*gal->comp_gamma_poly[component]*pow(rho_0,gal->comp_gamma_poly[component]-1.0);
                        new_step_z = gal->comp_mcmc_step_hydro[component]*sqrt(cs2/(2.0*pi*G*rho_0));
                    }
                    if(isinf(new_step_z)||new_step_z<min_step_z) new_step_z = min_step_z;
                    if(new_step_z>max_step_z) new_step_z = max_step_z;
                    if(gal->comp_spherical_hydro_eq[component]) {
                        new_step_x = new_step_z;
                        new_step_y = new_step_z;
                    }
                }
                new_step_r = sqrt(pow(new_step_x,2)+pow(new_step_y,2));
                new_step_r_sph = sqrt(pow(new_step_x,2)+pow(new_step_y,2)+pow(new_step_z,2));
                if(k == gal->mcmc_ntry-1) {
                    ref_r[k] = gal->r_cyl[i-1];
                    ref_r_sph[k] = gal->r_sph[i-1];
                    ref_theta[k] = gal->theta_cyl[i-1];
                    ref_x[k] = gal->x[i-1];
                    ref_y[k] = gal->y[i-1];
                    ref_z[k] = gal->z[i-1];
                } else {
                    switch(gal->comp_symmetry[component]) {
                        // Cylindrical symmetry
                        case 1:
                            ref_r[k] = prop_r[selected] + gsl_ran_gaussian(r[0],new_step_r);
                            ref_theta[k] = 2.0*pi*gsl_rng_uniform_pos(r[0]);
                            ref_z[k] = prop_z[selected] + gsl_ran_gaussian(r[0],new_step_z);
                            break;
                        // Spherical symmetry
                        case 2:
                            ref_r_sph[k] = prop_r_sph[selected] + gsl_ran_gaussian(r[0],step_r_sph);
                            ref_phi_sph[k] = acos(2*gsl_rng_uniform_pos(r[0])-1);
                            ref_theta[k] = 2.0*pi*gsl_rng_uniform_pos(r[0]);
                            ref_r[k] = ref_r_sph[k]*sin(ref_phi_sph[k]);
                            ref_z[k] = ref_r_sph[k]*cos(ref_phi_sph[k]);
                            break;
                        // No symmetry
                        default:
                            ref_x[k] = prop_x[selected] + gsl_ran_gaussian(r[0],new_step_x);
                            ref_y[k] = prop_y[selected] + gsl_ran_gaussian(r[0],new_step_y);
                            ref_z[k] = prop_z[selected] + gsl_ran_gaussian(r[0],new_step_z);
                            ref_r[k] = sqrt(pow(ref_x[k],2)+pow(ref_y[k],2));
                            ref_theta[k] = atan2(ref_y[k],ref_x[k]);
                            break;
                    }
                }
                switch(gal->comp_symmetry[component]) {
                    // Cylindrical symmetry
                    case 1:
                        q_x[k] = gsl_ran_gaussian_pdf(-prop_r[selected]+ref_r[k],new_step_r);
                        q_x[k] *= gsl_ran_gaussian_pdf(-prop_z[selected]+ref_z[k],new_step_z);
                        dv_x[k] = fabs(ref_r[k]);
                        break;
                    // Spherical symmetry
                    case 2:
                        q_x[k] = gsl_ran_gaussian_pdf(-prop_r_sph[selected]+ref_r_sph[k],new_step_r_sph);
                        dv_x[k] = fabs(pow(ref_r_sph[k],2));
                        break;
                    // No symmetry
                    default:
                        q_x[k] = gsl_ran_gaussian_pdf(-prop_x[selected]+ref_x[k],new_step_x);
                        q_x[k] *= gsl_ran_gaussian_pdf(-prop_y[selected]+ref_y[k],new_step_y);
                        q_x[k] *= gsl_ran_gaussian_pdf(-prop_z[selected]+ref_z[k],new_step_z);
                        dv_x[k] = 1.0;
                        break;
                }
                // Distribution function of the considered component
                if(gal->pseudo[0]) {
                    pi_x[k] = dv_x[k]*pseudo_density_gas_func(gal,ref_r[k],ref_theta[k],ref_z[k],1,density_model,component,gal->comp_spherical_hydro_eq[component]);
                } else {
                    pi_x[k] = dv_x[k]*density_functions_pool(gal,ref_r[k],ref_theta[k],ref_z[k],1,density_model,component);
                }
                lambda_x[k] = 1.0;///q_x[k];
                lambda_y[k] = 1.0;///q_y[k];
                w_x[k] = pi_x[k]*q_x[k]*lambda_x[k];
                w_y[k] = pi_y[k]*q_y[k]*lambda_y[k];
            }
            prob = min(1.0,(sum_dbl(w_y,gal->mcmc_ntry)/sum_dbl(w_x,gal->mcmc_ntry)));
            randval = gsl_rng_uniform_pos(r[0]);
            // Proposal accepted
            if(randval <= prob) {
                gal->r_cyl[i] = prop_r[selected];
                gal->theta_cyl[i] = prop_theta[selected];
                gal->phi_sph[i] = prop_phi_sph[selected];
                gal->z[i] = prop_z[selected];
                gal->rho[i] = pi_y[selected]/dv_y[selected];
		if(gal->comp_type[component]==0) {
		    d = gal->comp_spherical_hydro_eq[component]?sqrt(pow(gal->r_cyl[i],2)+pow(gal->z[i],2)):gal->z[i];
		    k_poly = d>rc?gal->comp_k_poly[component]*pow(d/rc,gal->comp_alpha_entropy[component]):gal->comp_k_poly[component];
		    gal->u[i] = k_poly*pow(gal->rho[i],gal->comp_gamma_poly[component]-1.0)/gamma_minus1;
		}
                acceptance += 1.0;
            // Proposal rejected, the particle keeps the same postion
            } else {
                gal->r_cyl[i] = gal->r_cyl[i-1];
                gal->z[i] = gal->z[i-1];
                gal->theta_cyl[i] = gal->theta_cyl[i-1];
                if(gal->comp_symmetry[component]==1) {
                    gal->theta_cyl[i] = prop_theta[selected];
                }
                if(gal->comp_symmetry[component]==2) {
                    gal->r_sph[i] = gal->r_sph[i-1];
                    gal->theta_cyl[i] = prop_theta[selected];
                    gal->phi_sph[i] = prop_phi_sph[selected];
                    gal->r_cyl[i] = gal->r_sph[i]*sin(gal->phi_sph[i-1]);
                    gal->z[i] = gal->r_sph[i]*cos(gal->phi_sph[i-1]);
                }
                gal->rho[i] = pi_x[gal->mcmc_ntry-1]/dv_x[gal->mcmc_ntry-1];
		// Set gas temperature
		if(gal->comp_type[component]==0) {
		    d = gal->comp_spherical_hydro_eq[component]?gal->r_sph[i]:gal->z[i];
		    k_poly = d>rc?gal->comp_k_poly[component]*pow(d/rc,gal->comp_alpha_entropy[component]):gal->comp_k_poly[component];
		    gal->u[i] = k_poly*pow(gal->rho[i],gal->comp_gamma_poly[component]-1.0)/gamma_minus1;
		    Tpart = gal->u[i]*gamma_minus1*protonmass*mu_mol/boltzmann*unit_energy/unit_mass;
		    if(Tpart>Tmax) Tmax = Tpart;
		}
            }	
	    //if(gal->pseudo[0]) printf("%lf %le\n",gal->r_sph[i],gal->rho[i]);
	    if(gal->rho[i]>gal->comp_dens_max[component]) gal->comp_dens_max[component] = gal->rho[i];
	    if(gal->rho[i]<gal->comp_dens_min[component]) gal->comp_dens_min[component] = gal->rho[i];
            // Temporary assign metallicity to local density value
            if(gal->comp_metal_gradient[component]) {
                gal->metal[i] = gal->rho[i];
                mean_metal = mean_metal+gal->metal[i];
            }

            // Updating the coordinate values
            gal->x[i] = gal->r_cyl[i]*cos(gal->theta_cyl[i]);
            gal->y[i] = gal->r_cyl[i]*sin(gal->theta_cyl[i]);

            // Updating the coordinate values
            gal->theta_cyl[i] = atan2(gal->y[i],gal->x[i])+pi;
            gal->r_sph[i] = sqrt(gal->x[i]*gal->x[i]+gal->y[i]*gal->y[i]+gal->z[i]*gal->z[i]);
            gal->phi_sph[i] = acos(gal->z[i]/gal->r_sph[i]);
	    gal->index[tid] = i;
        }
        // Rescale metallicity to user specified mean value
        if(gal->comp_metal_gradient[component]) {
            mean_metal = mean_metal/gal->comp_npart_pot[component];
            for(i = gal->comp_start_part[component]+1; i < gal->comp_start_part[component]+gal->comp_npart_pot[component]; ++i) {
                gal->metal[i] = gal->metal[i]*gal->comp_metal[component]/mean_metal;
            }
        }
        acceptance /= gal->comp_npart_pot[component];
        printf("[  acceptance=%.2lf  ]",acceptance);
        // Recursive calls if acceptance is outside the range [accept_min,accept_max]
        if(acceptance<gal->comp_accept_min[component]) {
            if(gal->pseudo[0]==1) {
                gal->comp_mcmc_step_hydro[component] /= 2.0;
                printf("\n/////\t\t\t---------------[         Warning         ][ Low acceptance->mcmc_hydro_step%d=%.2le ]\n",
		    component+1,gal->comp_mcmc_step_hydro[component]);
                printf("/////\t\t\t- Component %2d [       recomputing       ]",component+1,gal->comp_profile_name[component]);
            } else {
                gal->comp_mcmc_step[component] /= 2.0;
                printf("\n/////\t\t---------------[         Warning         ][ Low MCMC acceptance->mcmc_step%d=%.2le  ]\n",component+1,gal->comp_mcmc_step[component]);
                printf("/////\t\t- Component %2d [       recomputing       ]",component+1,gal->comp_profile_name[component]);
            }
            fflush(stdout);
            mcmc_metropolis_hasting_ntry(gal,component,gal->comp_model[component]);
        }
        if(acceptance>gal->comp_accept_max[component]) {
            if(gal->pseudo[0]==1) {
                gal->comp_mcmc_step_hydro[component] *= 2.0;
                printf("\n/////\t\t\t---------------[         Warning         ][ High acceptance->mcmc_step_hydro%d=%.2le]\n",
		    component+1,gal->comp_mcmc_step_hydro[component]);
                printf("/////\t\t\t- Component %2d [       recomputing       ]",component+1,gal->comp_profile_name[component]);
            } else {
                gal->comp_mcmc_step[component] *= 2.0;
                printf("\n/////\t\t---------------[         Warning         ][ High MCMC acceptance->mcmc_step%d=%.2le ]\n",component+1,gal->comp_mcmc_step[component]);
                printf("/////\t\t- Component %2d [       recomputing       ]",component+1,gal->comp_profile_name[component]);
            }
            fflush(stdout);
            mcmc_metropolis_hasting_ntry(gal,component,gal->comp_model[component]);
        }
        if(acceptance<gal->comp_accept_max[component] && acceptance>gal->comp_accept_min[component]) {
            printf("[rho_min=%4.2le rho_max=%4.2le H/cc]",gal->comp_dens_min[component]*unit_nh,gal->comp_dens_max[component]*unit_nh);
	    if(gal->comp_type[component]==0 && gal->comp_hydro_eq[component]==1) {
	        printf("[Tmax=%4.2le K]",Tmax);
	    }
	    printf("\n");
	}
        if(AllVars.MeanPartDist) printf("/////\t\t\tMean inter-particle distance -> %lf [kpc]\n",mean_interparticle_distance(gal,component));
    } else {
        gal->x[i] = 0.;
        gal->y[i] = 0.;
        gal->z[i] = 0.;
        gal->r_cyl[i] = 0.;
        gal->theta_cyl[i] = 0.;
        gal->r_sph[i] = 0.;
        gal->phi_sph[i] = 0.;
        printf("\n");
    }
    free(prop_x);
    free(prop_y);
    free(prop_z);
    free(prop_r);
    free(prop_r_sph);
    free(prop_theta);
    free(prop_phi_sph);
    free(ref_x);
    free(ref_y);
    free(ref_z);
    free(ref_r);
    free(ref_r_sph);
    free(ref_theta);
    free(ref_phi_sph);
    free(pi_x);
    free(pi_y);
    free(q_x);
    free(q_y);
    free(lambda_x);
    free(lambda_y);
    free(w_x);
    free(w_y);
    free(weights);
    free(dv_x);
    free(dv_y);

    return;
}

void mcmc_metropolis_hasting_stream(stream *st, int component, int density_model) {
    unsigned long int i,j;
    double pi_x, pi_y, q_x, q_y, prob, *radius;
    double theta, phi, randval, proposal, prop_r, prop_theta, prop_x, prop_y, prop_z;
    double step_r, step_x, step_y, step_z, prev_step_x, prev_step_y, prev_step_z;
    double acceptance, mean_metal;

    if(st->comp_npart[component]>0) {
        // Use the Metropolis algorithm to place the disk particles.
        // We start the Monte Carlo Markov Chain with a realistic particle position
        theta = 2.0*pi*gsl_rng_uniform_pos(r[0]);
        i = st->comp_start_part[component];
        st->x[i] = 0.;
        st->y[i] = 0.;
        st->z[i] = 0.;
        st->r_cyl[i] = sqrt(st->x[i]*st->x[i]+st->y[i]*st->y[i]);
        st->theta_cyl[i] = atan2(st->y[i],st->x[i]);
        st->r_sph[i] = sqrt(st->x[i]*st->x[i]+st->y[i]*st->y[i]+st->z[i]*st->z[i]);
        st->phi_sph[i] = acos(st->z[i]/st->r_sph[i]);
        step_r = st->comp_mcmc_step[component]*st->comp_length[component];
        step_x = st->comp_mcmc_step[component]*st->comp_length[component];
        step_y = st->comp_mcmc_step[component]*st->comp_length[component];
        step_z = st->comp_mcmc_step[component]*st->comp_length[component];

        acceptance = 0.;
        mean_metal = 0.;
        // Burning period
        for(j = 1; j<(int)(0.1*st->comp_npart[component]); ++j) {
            // Generating a proposal
            prop_x = st->x[i] + gsl_ran_gaussian(r[0],step_x);
            prop_y = st->y[i] + gsl_ran_gaussian(r[0],step_y);
            prop_z = st->z[i] + gsl_ran_gaussian(r[0],step_z);
            prop_r = sqrt(pow(prop_x,2)+pow(prop_y,2));
            prop_theta = atan2(prop_y,prop_x);
            // Distribution function of the considered component
            // Density of the component times the integrand of spherical volume element
            pi_x = density_functions_stream_pool(st,st->r_cyl[i],st->theta_cyl[i],st->z[i],density_model,component);
            pi_y = density_functions_stream_pool(st,prop_r,prop_theta,prop_z,density_model,component);
            prob = min(1.0,(pi_y/pi_x)); //*(q_x/q_y));
            randval = gsl_rng_uniform_pos(r[0]);
            if(randval <= prob) {
                st->r_cyl[i] = prop_r;
                st->theta_cyl[i] = prop_theta;
                st->z[i] = prop_z;
            }
        }
        // Updating the coordinate values
        st->x[i] = st->r_cyl[i]*cos(st->theta_cyl[i]);
        st->y[i] = st->r_cyl[i]*sin(st->theta_cyl[i]);
        // Updating the coordinate values
        st->theta_cyl[i] = atan2(st->y[i],st->x[i]);
        st->r_sph[i] = sqrt(st->x[i]*st->x[i]+st->y[i]*st->y[i]+st->z[i]*st->z[i]);
        st->phi_sph[i] = acos(st->z[i]/st->r_sph[i]);
        // Filling the Markov Chain
        for(i = st->comp_start_part[component]+1; i<st->comp_start_part[component]+st->comp_npart[component]; ++i) {
            // Generating a proposal
            prop_x = st->x[i-1] + gsl_ran_gaussian(r[0],step_x);
            prop_y = st->y[i-1] + gsl_ran_gaussian(r[0],step_y);
            prop_z = st->z[i-1] + gsl_ran_gaussian(r[0],step_z);
            prop_r = sqrt(pow(prop_x,2)+pow(prop_y,2));
            prop_theta = atan2(prop_y,prop_x);
            prev_step_x = step_x;
            prev_step_y = step_y;
            prev_step_z = step_z;
            // Distribution function of the considered component
            pi_x = density_functions_stream_pool(st,st->r_cyl[i-1],st->theta_cyl[i-1],st->z[i-1],density_model,component);
            pi_y = density_functions_stream_pool(st,prop_r,prop_theta,prop_z,density_model,component);
            q_x = gsl_ran_gaussian_pdf(-prop_x+st->x[i-1],step_x);
            q_x *= gsl_ran_gaussian_pdf(-prop_y+st->y[i-1],step_y);
            q_x *= gsl_ran_gaussian_pdf(-prop_z+st->z[i-1],step_z);
            q_y = gsl_ran_gaussian_pdf(-st->x[i-1]+prop_x,prev_step_x);
            q_y *= gsl_ran_gaussian_pdf(-st->y[i-1]+prop_y,prev_step_y);
            q_y *= gsl_ran_gaussian_pdf(-st->z[i-1]+prop_z,prev_step_z);
            prob = min(1.0,(pi_y/pi_x)*(q_x/q_y));
            randval = gsl_rng_uniform_pos(r[0]);
            if(randval <= prob) {
                st->r_cyl[i] = prop_r;
                st->theta_cyl[i] = prop_theta;
                st->z[i] = prop_z;
                st->rho[i] = pi_y/fabs(st->r_cyl[i]);
                acceptance += 1.0;
            } else {
                st->r_cyl[i] = st->r_cyl[i-1];
                st->theta_cyl[i] = st->theta_cyl[i-1];
                st->z[i] = st->z[i-1];
                st->rho[i] = pi_x/fabs(st->r_cyl[i]);
            }
            // Temporary assign metallicity to local density value
            if(st->comp_metal_gradient[component]) {
                st->metal[i] = st->rho[i];
                mean_metal = mean_metal+st->metal[i];
            }
            // Updating the coordinate values
            st->x[i] = st->r_cyl[i]*cos(st->theta_cyl[i]);
            st->y[i] = st->r_cyl[i]*sin(st->theta_cyl[i]);
            // Updating the coordinate values
            st->theta_cyl[i] = atan2(st->y[i],st->x[i]);
            st->r_sph[i] = sqrt(st->x[i]*st->x[i]+st->y[i]*st->y[i]+st->z[i]*st->z[i]);
            st->phi_sph[i] = acos(st->z[i]/st->r_sph[i]);
        }
        // Rescale metallicity to user specified mean value
        if(st->comp_metal_gradient[component]) {
            mean_metal = mean_metal/st->comp_npart[component];
            for(i = st->comp_start_part[component]+1; i<st->comp_start_part[component]+st->comp_npart[component]; ++i) {
                st->metal[i] = st->metal[i]*st->comp_metal[component]/mean_metal;
            }
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
    size_t neval;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(AllVars.GslWorkspaceSize);
    gsl_function F;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    F.function = &integrand_density_func;
    F.params = gal;
    gal->storage[0][tid] = r;
    gal->storage[1][tid] = theta;
    gal->storage[2][tid] = cut;
    gal->selected_comp[tid] = component;

    //gsl_integration_qag(&F,gal->comp_cut[component]*gal->comp_flatz[component],gal->comp_cut[component]*gal->comp_flatz[component],epsabs,epsrel,AllVars.GslWorkspaceSize,key,w,&surface_density,&error);
    gsl_integration_qng(&F,-gal->comp_cut[component]*gal->comp_flatz[component],gal->comp_cut[component]*gal->comp_flatz[component],epsabs,epsrel,&surface_density,&error,&neval);

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

    r_temp = gal->storage[0][tid];
    theta_temp = gal->storage[1][tid];
    cut = gal->storage[2][tid];

    // This function should not be used with the pseudo-density function
    if(gal->pseudo[tid]) {
        rho = pseudo_density_gas_func(gal,r_temp,theta_temp,z,cut,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid],gal->comp_spherical_hydro_eq[gal->selected_comp[tid]]);
    } else {
        rho = density_functions_pool(gal,r_temp,theta_temp,z,cut,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
    }

    return rho;
}

// This function calculates the force on a test particle due to the halo
// along the r axis. This is used to calculate the velocity dispersion.
double cumulative_mass_func(galaxy *gal, double radius, int component) {

    int status;
    double integral,error;
    int tid;
    size_t neval;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    gsl_function F;

    F.function = &d_cumulative_mass_func1;
    F.params = gal;

    gal->selected_comp[tid] = component;

    gsl_integration_qng(&F,0.0,radius,epsabs,epsrel,&integral,&error,&neval);
    //gsl_integration_qag(&F,0.0,radius,epsabs,epsrel,AllVars.GslWorkspaceSize,key,w[0],&integral,&error);

    return integral;
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

    F.function = &d_cumulative_mass_func2;
    F.params = gal;

    gal->storage[6][tid] = r;

    gsl_integration_qng(&F,0.0,2.0*pi,epsabs,epsrel,&integrand,&error,&neval);
    //gsl_integration_qag(&F,0.0,2.0*pi,epsabs,epsrel,AllVars.GslWorkspaceSize,key,w[0],&integrand,&error);

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

    st->storage[0][tid] = r;
    st->storage[1][tid] = theta;

    st->selected_comp[tid] = component;
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

    r_temp = st->storage[0][tid];
    theta_temp = st->storage[1][tid];

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

    F.function = &d_cumulative_mass_func_stream;
    F.params = st;

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
double pseudo_density_gas_func(galaxy *gal, double r, double theta, double z, int cut, int density_model, int component, int spherical) {

    int tid;
    double density, delta_phi, delta_phi_rc, delta_phi_core, rho_0, rho_c, save1, save2, rsph, x ,y;
    double sigma1, sigma2, smooth_factor1, smooth_factor2, smooth_factor3, hx, hy, hz, n, m;
    double k0, gamma0, rc, alpha, d;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    hx = gal->comp_cut_hydro_eq[component]*gal->comp_flatx_cut[component];
    hy = gal->comp_cut_hydro_eq[component]*gal->comp_flaty_cut[component];
    hz = gal->comp_cut_hydro_eq[component]*gal->comp_flatz_cut[component];
    k0 = gal->comp_k_poly[component];
    gamma0 = gal->comp_gamma_poly[component];
    rc = gal->comp_rc_entropy[component];
    alpha = gal->comp_alpha_entropy[component];
    if(spherical) {
        rsph = sqrt(r*r+z*z);
        x = r*cos(theta);
        y = r*sin(theta);
        // Computing the potential gradient between 0 and rsph
        save1 = gal->theta_cyl[gal->index[tid]];
        save2 = gal->phi_sph[gal->index[tid]];
        gal->theta_cyl[gal->index[tid]] = theta;
	if(rsph>0.) {
            gal->phi_sph[gal->index[tid]] = acos(z/rsph);
	} else {
            gal->phi_sph[gal->index[tid]] = 0.;
	}
        delta_phi = galaxyrsph_potential_wrapper_func(rsph,gal)-galaxyrsph_potential_wrapper_func(0.,gal);
	if(alpha>0.) {
            delta_phi_rc = galaxyrsph_potential_wrapper_func(rsph,gal)-galaxyrsph_potential_wrapper_func(rc,gal);
            delta_phi_core = galaxyrsph_potential_wrapper_func(rc,gal)-galaxyrsph_potential_wrapper_func(0.,gal);
	}
        // Central density
	rho_0 = gal->comp_dens_init[component]/unit_nh;
        n = sqrt(pow(x/hx,2.0)+pow(y/hy,2.0)+pow(z/hz,2.0));
        gal->theta_cyl[gal->index[tid]] = save1;
        gal->phi_sph[gal->index[tid]] = save2;
    } else {
        // Computing the potential gradient between 0 and z
        save1 = gal->x[gal->index[tid]];
        save2 = gal->y[gal->index[tid]];
        gal->x[gal->index[tid]] = r*cos(theta);
        gal->y[gal->index[tid]] = r*sin(theta);
        delta_phi = galaxyz_potential_wrapper_func(z,gal)-galaxyz_potential_wrapper_func(0.,gal);
	if(alpha>0.) {
           delta_phi_rc = galaxyz_potential_wrapper_func(z,gal)-galaxyz_potential_wrapper_func(rc,gal);
           delta_phi_core = galaxyz_potential_wrapper_func(rc,gal)-galaxyz_potential_wrapper_func(0.,gal);
	}
        // Density in the xy-plane
        if(gamma0<1.00001) {
            rho_0 = get_midplane_density(gal,gal->x[gal->index[tid]],gal->y[gal->index[tid]]);
	} else {
            rho_0 = density_functions_pool(gal,r,theta,0.,1,gal->comp_model[component],component);
	}
        n = sqrt(pow(gal->x[gal->index[tid]]/hx,2.0)+pow(gal->y[gal->index[tid]]/hy,2.0));
        m = sqrt(pow(z/hz,2.0));
        gal->x[gal->index[tid]] = save1;
        gal->y[gal->index[tid]] = save2;
    }

    // Safety value for the polytropic index
    gamma0 = max(gamma0,1.00001);
    // Switch between spherical or cylindrical
    d = spherical?fabs(rsph):fabs(z);
    // Hydrostatic equilibrium requires the following density
    // Entropy core case
    if(alpha>0.) {
        if(d<=rc) {
            density = rho_0*pow(max(1-(gamma0-1)/(k0*gamma0)*delta_phi*pow(rho_0,1-gamma0),0.),1.0/(gamma0-1));
        } else {
            rho_c = rho_0*pow(max(1-(gamma0-1)/(k0*gamma0)*delta_phi_core*pow(rho_0,1-gamma0),0.),1.0/(gamma0-1));
            density = pow(d/rc,-alpha/gamma0)*rho_c*pow(max(1-(gamma0-1)/(k0*gamma0)*pow(rho_c,1.0-gamma0)*delta_phi_rc,0.),1.0/(gamma0-1));
        }
    // No entropy core
    } else {
	density = rho_0*pow(max(1-(gamma0-1)/(k0*gamma0)*delta_phi*pow(rho_0,1-gamma0),0.),1.0/(gamma0-1));
    }

    if(cut==1) {
        sigma1 = gal->comp_sigma_cut[component];
        sigma2 = gal->comp_sigma_cut_in[component];
        smooth_factor1 = 1-0.5*(1+erf((n-1.0)/(sigma1*sqrt(2))));
        smooth_factor2 = 0.5*(1+erf((n-1.0)/(sigma2*sqrt(2))));
        smooth_factor3 = 1-0.5*(1+erf((m-1.0)/(sigma1*sqrt(2))));
        if(spherical){
            density *= smooth_factor1;
            if(fabs(r)<gal->comp_cut_in[component]) density *= smooth_factor2;
        } else {
            if(r>gal->comp_cut[component]+3*gal->comp_sigma_cut[component]) density = 0.;
            density *= smooth_factor3;
            if(r<gal->comp_cut_in[component]) density *= smooth_factor2;
        }
        if(density<gal->comp_cut_dens[component]/unit_nh) density = 0.;
    }

    // Return a pseudo density
    return density;
}

double midplane_density_gas_func(galaxy *gal, gsl_integration_workspace *w, double x, double y) {
    int i,status,pseudo_save;
    int tid;
    double result, integral, error, radius, theta, zmax, rho0, surface_density;
    size_t neval;

    gsl_function F;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    gal->storage[3][tid] = x;
    gal->storage[4][tid] = y;

    radius = sqrt(x*x+y*y);
    theta = atan2(y,x);
    pseudo_save = gal->pseudo[tid];
    gal->pseudo[tid] = 0;

    rho0 = 0.;
    for(i = 0; i<AllVars.MaxCompNumber; i++) {
        if(gal->comp_type[i]==0 && gal->comp_npart[i]>1) {
	    // Planar density is defined by surface_density/integral(exp(phi/cs2)dz)
	    if(gal->comp_gamma_poly[i]<1.00001) {
	        gal->selected_comp[tid] = i;
    	        zmax = gal->comp_scale_height[i];
    	        F.function = &dmidplane_density_gas_func;
    	        F.params = gal;
    	        gsl_integration_qng(&F,-10.*zmax,10.*zmax,epsabs,epsrel,&integral,&error,&neval);
	        surface_density = surface_density_func(gal,radius,theta,1,i);
                rho0 += surface_density/integral;
	    // Planar density is defined by the input density profile
	    } else {
	        rho0 += density_functions_pool(gal,radius,theta,0.,1,gal->comp_model[i],i);
	    }
	}
    }
    gal->pseudo[tid] = pseudo_save;
    
    return rho0;
}

// This is the integrand for the previous function. It is setup to work with
// the GSL_qags structures.
static double dmidplane_density_gas_func(double z, void *params) {

    int tid, component;
    double integrand, delta_pot, x, y, cs2;
    galaxy *gal = (galaxy *) params;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    x = gal->storage[3][tid];
    y = gal->storage[4][tid];

    delta_pot = galaxy_total_potential(gal,x,y,z,1,0)-galaxy_total_potential(gal,x,y,0.,1,0);
    cs2 = pow(gal->comp_cs_init[gal->selected_comp[tid]],2.0);
    integrand = exp(-delta_pot/cs2);
	
    return integrand;
}

// This function fills a 2D grid with the value of the gas density in the midplane
// The use of a 2D grid intends to lower the computation time
void fill_midplane_dens_grid(galaxy *gal) {
    int i,j;
    unsigned long int k;
    double x,y;

    gsl_integration_workspace **wk;

    wk = (gsl_integration_workspace **) malloc(AllVars.Nthreads * sizeof(gsl_integration_workspace *));
    for(i = 0; i<AllVars.Nthreads; i++) {
        wk[i] = gsl_integration_workspace_alloc(AllVars.GslWorkspaceSize);
    }

#pragma omp parallel shared(gal) private(x,y,i,j)
    for(i = 0; i<gal->ngrid_dens[0]; i++) {
        for(j = 0; j<gal->ngrid_dens[1]; j++) {
            // Get thread ID
#if USE_THREADS == 1
            int tid = omp_get_thread_num();
#else
            int tid = 0;
#endif
            x = ((double)i-((double)(gal->ngrid_dens[0]/2)-0.5))*gal->dx_dens;
            y = ((double)j-((double)(gal->ngrid_dens[1]/2)-0.5))*gal->dx_dens;
            gal->midplane_dens[i][j] = midplane_density_gas_func(gal,wk[tid],x,y);
        }
    }
    for(i = 0; i<AllVars.Nthreads; i++) {
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
double disk_scale_length_func(galaxy *gal, double c, int component) {

    int i;
    double f_c, f_r, disk_scale;

    f_c = f_c_func(c);
    f_r = f_r_func(gal,component);

    disk_scale = (1.0/sqrt(2.0))*(gal->comp_angmom_frac[component]/gal->comp_mass_frac[component])*gal->lambda*gal->r200*(f_r/sqrt(f_c));

    return disk_scale;
}

// This function returns the energy fraction function of MMW for the
// angular momentum.
double f_c_func(double c) {

    double a, upper, lower;

    upper = c*(1.0 - 1.0/pow(1.0+c,2.0) - 2.0*log(1.0+c)/(1.0+c));
    a = -log(1.0+c) + c/(1.0+c);
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

// This is the azimuthal circular velocity fraction, f_s. It is needed
// to determine the azimuthal streaming velocity.
double f_r_func(galaxy *gal, int component) {

    double f_r, c, m_d, lambdap, power;

    m_d = gal->comp_mass_frac[component];
    c = gal->comp_concentration[gal->index_halo];
    lambdap = gal->lambda*gal->comp_angmom_frac[component]/m_d;

    power = -0.06+2.71*m_d+0.0047/lambdap;
    f_r = pow(lambdap/0.1,power)*(1.0-3.0*m_d+5.2*m_d*m_d)*(1.0-0.019*c+0.00025*c*c+0.52/c);

    return f_r;
}

// This is the azimuthal circular velocity fraction, f_s. It is needed
// to determine the azimuthal streaming velocity.
double f_v_func(galaxy *gal, int component) {

    double f_v, c, m_d, lambdap, power;

    m_d = gal->comp_mass_frac[component];
    c = gal->comp_concentration[gal->index_halo];
    lambdap = gal->lambda*gal->comp_angmom_frac[component]/gal->comp_mass_frac[component];

    power = -2.67*m_d-0.0038/lambdap+0.2*lambdap;
    f_v = pow(lambdap/0.1,power)*(1.0+4.35*m_d-3.76*m_d*m_d)*(1.0+0.057*c-0.00034*c*c-1.54/c)/sqrt(-c/(1.0+c)+log(1.0+c));

    return f_v;
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
        for(i = gal->comp_start_part[component]; i < gal->comp_start_part[component]+gal->comp_npart[component]; i++) {
            // For a given particle, we compute the distances to all the other particles
            // with the same type
            closest_part = gal->boxsize[0];
            for(j = i+1; j < gal->comp_start_part[component]+gal->comp_npart[component]; j++) {
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
    for(j = 0; j<AllVars.MaxCompNumber; j++) {
        if(gal->comp_npart[j]>0) printf("/////\t\t- Component %2d -> m=%.2e Msol\n",j+1,(gal->comp_mass[j]*unit_mass/solarmass)/(gal->comp_npart[j]));
        // Filling the arrays of the &galaxy structure
        for(i = gal->comp_start_part[j]; i < gal->comp_start_part[j] + gal->comp_npart_pot[j]; i++) {
            gal->mass[i] = gal->comp_mass[j]/gal->comp_npart[j];
        }
    }
    return;
}

// This function computes the halo concentration parameter
// for a given redshift and a given mass according to Dutton et al. 2014
double halo_concentration(double m200, double z) {
    double a, b, log10_c200;

    a = 0.520+(0.905-0.520)*exp(-0.617*pow(z,1.21));
    b = -0.101+0.026*z;
    log10_c200 = a+b*log10(m200*unit_mass/solarmass/(1e12/AllVars.h));

    return pow(10.,log10_c200);
}

// This function computes the stellar mass expected from abundance matching (Behroozi et al. 2013)
double halo_abundance(double m200, double z) {
    double a, nu, log_eps, log_M1, alpha, delta, beta, x, f0, fx, mstar; 

    a = 1./(1.+z);
    nu = exp(-4*a*a);
    log_eps = -1.777+(-0.006*(a-1)+0.000*z)*nu-0.119*(a-1);
    log_M1 = 11.514+(-1.793*(a-1)+(-0.251)*z)*nu;
    alpha = -1.412+(0.731*(a-1))*nu;
    delta = 3.508+(2.608*(a-1)+(-0.043)*z)*nu;
    beta = 0.316+(1.319*(a-1)+0.279*z)*nu;
    x = 0.;
    f0 = -log10(pow(10.,alpha*x)+1)+delta*pow(log10(1.+exp(x)),beta)/(1.+exp(pow(10.,-x)));
    x = log10(m200*unit_mass/solarmass)-log_M1;
    fx = -log10(pow(10.,alpha*x)+1)+delta*pow(log10(1.+exp(x)),beta)/(1.+exp(pow(10.,-x)));
    mstar = pow(10.,(log_M1+log_eps)+fx-f0);

    return mstar*solarmass/unit_mass;
}

// This function computes the metallicity expected for a given stellar mass
// from the Tremonti et al. 2004 Mass-Metallicity relation at z=0
// The relation is shiffted lineary to match the Erb et al. 2006 fit at z=2.0
double galaxy_metallicity(double mstar, double z) {
    double Z,shift;

    // Erb+2006 MZ relation shift
    shift = 0.56/2.0*z;
    // Tremonti+2004 MZ relation
    Z = (-1.492+1.847*log10(mstar*unit_mass/solarmass)-0.08026*pow(log10(mstar*unit_mass/solarmass),2)-shift)/logOH_solar;

    return Z;
}
