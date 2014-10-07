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
 / Date: October 2014
 /
 *///---------------------------------------------------------------------------

// Include the DICE header
#include "dice.h"

// Create a galaxy
int create_galaxy(galaxy *gal, char *fname, int info) {
	
	unsigned long int i;
	int j,nt;
	// Molecular weigth t compute the gas internal energy
	double molecular_weight;
	// Hubble parameter
	double H;
	// Seed for the random number generator
	long seed;
	double halo_mass,disk_mass,gas_mass,bulge_mass;
	double cutted_halo_mass,cutted_disk_mass,cutted_gas_mass,cutted_bulge_mass;
	double halo_npart,disk_npart,gas_npart,bulge_npart;
	double max_gas_radius;
	double baryonic_fraction;
	double BD_fraction,BT_fraction,gas_fraction;
	double effective_mass_factor;
		
	// Allocate component arrays.
	if (!(gal->comp_mass=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mass array.\n");
		return -1;
	}
	if (!(gal->comp_model=calloc(AllVars.MaxCompNumber,sizeof(int)))) {
		fprintf(stderr,"Unable to allocate particle comp_model array.\n");
		return -1;
	}
	if (!(gal->comp_npart=calloc(AllVars.MaxCompNumber,sizeof(unsigned long int)))) {
		fprintf(stderr,"Unable to allocate particle comp_npart array.\n");
		return -1;
	}
	if (!(gal->comp_npart_pot=calloc(AllVars.MaxCompNumber,sizeof(unsigned long int)))) {
		fprintf(stderr,"Unable to allocate particle comp_npart_pot array.\n");
		return -1;
	}
	if (!(gal->comp_start_part=calloc(AllVars.MaxCompNumber,sizeof(unsigned long int)))) {
		fprintf(stderr,"Unable to allocate particle comp_start_part array.\n");
		return -1;
	}
	if (!(gal->comp_cutted_mass=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_cutted_mass array.\n");
		return -1;
	}
	if (!(gal->comp_scale_length=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_scale_length array.\n");
		return -1;
	}
	if (!(gal->comp_scale_height=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_scale_height array.\n");
		return -1;
	}
	if (!(gal->comp_cut=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_cut array.\n");
		return -1;
	}
	if (!(gal->comp_flat=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_flat array.\n");
		return -1;
	}
	if (!(gal->comp_mcmc_step=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mcmc_step array.\n");
		return -1;
	}
	if (!(gal->comp_mcmc_step_hydro=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mcmc_step array.\n");
		return -1;
	}
	if (!(gal->comp_vmax=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_vmax array.\n");
		return -1;
	}
	if (!(gal->comp_mass_frac=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mass_frac array.\n");
		return -1;
	}
	if (!(gal->comp_type=calloc(AllVars.MaxCompNumber,sizeof(int)))) {
		fprintf(stderr,"Unable to allocate particle comp_type array.\n");
		return -1;
	}
	if (!(gal->comp_bool=calloc(AllVars.MaxCompNumber,sizeof(int)))) {
		fprintf(stderr,"Unable to allocate particle comp_bool array.\n");
		return -1;
	}
	if (!(gal->comp_cut_dens=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_cut_dens array.\n");
		return -1;
	}
	if (!(gal->comp_concentration=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_concentration array.\n");
		return -1;
	}
	if (!(gal->comp_streaming_fraction=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_streaming_fraction array.\n");
		return -1;
	}
	if (!(gal->comp_theta_sph=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_theta_sph array.\n");
		return -1;
	}
	if (!(gal->comp_phi_sph=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_phi_sph array.\n");
		return -1;
	}
	if (!(gal->comp_metal=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_metal array.\n");
		return -1;
	}
	if (!(gal->comp_t_init=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_t_init array.\n");
		return -1;
	}
	if (!(gal->comp_u_init=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_u_init array.\n");
		return -1;
	}
	if (!(gal->comp_cs_init=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_cs_init array.\n");
		return -1;
	}
	if (!(gal->comp_mean_age=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mean_age array.\n");
		return -1;
	}
	if (!(gal->comp_min_age=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_min_age array.\n");
		return -1;
	}

	// Parse the galaxy parameters file
	if(parse_galaxy_file(gal,fname) != 0) {
		fprintf(stderr,"Unable to find the galaxy parameters file\n");
		return -1;
	}

	// Allocate pseudo all the threads.
	if (!(gal->pseudo=calloc(AllVars.Nthreads,sizeof(int)))) {
		fprintf(stderr,"Unable to allocate pseudo-density switch.\n");
		return -1;
	}
	
	// Allocate component selector all the threads.
	if (!(gal->selected_comp=calloc(AllVars.Nthreads,sizeof(int)))) {
		fprintf(stderr,"Unable to allocate pseudo-density switch.\n");
		return -1;
	}
    
	// Padding potential grid
	gal->ngrid 			= pow(2,gal->level_grid);
	gal->ngrid_padded 	= 2*gal->ngrid;
	gal->ngrid_dens 	= pow(2,gal->level_grid);
	// Create the random number generator environment, if not already done
	if (random_number_set != 1) {
		gsl_rng_env_setup();
		T = gsl_rng_mt19937;
		r = (gsl_rng **) malloc(AllVars.Nthreads * sizeof(gsl_rng *));
		for(i=0;i<AllVars.Nthreads;i++) {
			r[i] = gsl_rng_alloc(T);
			gsl_rng_set(r[i],gal->seed);
		}
		random_number_set = 1;
	}

	// Set a bunch of constants that defines the galaxy's structure.
	// Set the Hubble parameter
	H = sqrt(H0*H0*(omega_m*pow(1+gal->redshift,3.0)+omega_k*pow(1+gal->redshift,2.0)+omega_l));
	// Get the Virial velocity from the parser
	gal->v200 = gal->v200*1E5;
	// Set the Virial mass
	gal->m200 = (gal->v200*gal->v200*gal->v200)/(10.0*G*H);
	// Set the Virial radius
	gal->r200 = gal->v200/(10.0*H*kpc);
    
	// Reset Q_min to some arbitrary value so it can be calculated later.
	gal->Q_min = 100.0;

	// Set the size of the cell for the Potential-Mesh (PM) computation
	gal->space[0] 		= gal->boxsize/((double)gal->ngrid);
	gal->space[1] 		= gal->boxsize/((double)gal->ngrid);
	gal->space[2] 		= gal->boxsize/((double)gal->ngrid);
	
	gal->boxsize_dens 	= 0.;
	for(i=0; i<AllVars.MaxCompNumber; i++){
		if(gal->comp_type[i]==0 && 2.0*gal->comp_cut[i]>gal->boxsize_dens){
			gal->boxsize_dens = 2.1*gal->comp_cut[i];
		}
	}

	gal->space_dens[0] 	= gal->boxsize_dens/((double)gal->ngrid_dens);
	gal->space_dens[1] 	= gal->boxsize_dens/((double)gal->ngrid_dens);

	allocate_galaxy_storage_variable(gal,10);
	
	// Set the GSL QAG integration workspace
	w = (gsl_integration_workspace **) malloc(AllVars.Nthreads * sizeof(gsl_integration_workspace *));
	for(i=0;i<AllVars.Nthreads;i++) {
		w[i] = gsl_integration_workspace_alloc(GSL_WORKSPACE_SIZE);
	}
	
	// Initialisations
	gal->ntot_part 			= (unsigned long int)0;
	gal->comp_start_part[0] = (unsigned long int)0;
	gal->total_mass 		= 0.;
	max_gas_radius			= 0.0;
	halo_mass				= 0.0;
	disk_mass				= 0.0;
	gas_mass				= 0.0;
	bulge_mass				= 0.0;
	cutted_halo_mass		= 0.0;
	cutted_disk_mass		= 0.0;
	cutted_gas_mass			= 0.0;
	cutted_bulge_mass		= 0.0;
	gal->num_part[0]		= 0;
	gal->num_part[1]		= 0;
	gal->num_part[2]		= 0;
	gal->num_part[3]		= 0;
	gal->num_part_pot[0]	= 0;
	gal->num_part_pot[1]	= 0;
	gal->num_part_pot[2]	= 0;
	gal->num_part_pot[3]	= 0;
	effective_mass_factor	= 0.;
	
	// Do not work with pseudo densities yet
	for(i=0; i<AllVars.Nthreads; i++) gal->pseudo[i] = 0;
	// Set up component properties
	for(i=0; i<AllVars.MaxCompNumber; i++){
		// Total number of particules
		if(gal->comp_npart_pot[i]==0) gal->comp_npart_pot[i] = gal->comp_npart[i];
		if(gal->comp_npart[i]==0) gal->comp_npart_pot[i] = 0;
		gal->ntot_part_pot 				+= gal->comp_npart_pot[i];
		gal->ntot_part 					+= gal->comp_npart[i];
		if(gal->comp_type[i]>1) gal->ntot_part_stars += gal->comp_npart[i];
		effective_mass_factor			+= gal->comp_mass_frac[i];
		// Set the start index for each component
		if(i>0) gal->comp_start_part[i] = gal->comp_start_part[i-1]+gal->comp_npart_pot[i-1];
		// Computing mass component
		gal->comp_mass[i] 			= gal->m200*gal->comp_mass_frac[i];
		// Set the halo scalelength
		if(gal->comp_concentration[i]>0 && gal->comp_type[i]==1) {
			gal->comp_scale_length[i] = gal->r200/gal->comp_concentration[i];
		}
		// Set the component scale height
		gal->comp_scale_height[i] 	= gal->comp_flat[i]*gal->comp_scale_length[i];

		// We check how many components are present in the galaxy
		// in order to compute masses
		// Booleans to test the presence of a given component
		if(gal->comp_npart[i] > 0) {
			gal->comp_bool[i] = 1;
			gal->n_component++;
		}
		else gal->comp_bool[i] 		= 0;
		gal->comp_cut_dens[i] 		= 0.0;
		gal->comp_cut_dens[i] 		= density_functions_pool(gal,gal->comp_cut[i],0.,0.,0,gal->comp_model[i],i);
		gsl_set_error_handler_off();
		gal->comp_cutted_mass[i] 	= cumulative_mass_func(gal,gal->comp_cut[i],i);
		gal->total_mass				+= gal->comp_cutted_mass[i];

		// Check for disk scalelength values lower or equal to 0
		// If it is the case, compute the disk scalelength according to Mo et al. 1998
		//if(gal->comp_scale_length[i] <= 0 && (gal->comp_type[i]==1 || gal->comp_type[i]==2)) gal->comp_scale_length[i] = disk_scale_length_func(gal);
		if(gal->comp_type[i]==0 && gal->comp_scale_length[i]>max_gas_radius) max_gas_radius = gal->comp_scale_length[i];
		if(gal->comp_type[j]==1 && gal->lambda>=0) gal->comp_streaming_fraction[j] = f_s_func(gal->comp_concentration[j],gal->lambda);
		// Checking total mass
		if(gal->comp_type[i]==0) gas_mass 			+= gal->comp_bool[i]*gal->comp_mass[i];
		if(gal->comp_type[i]==1) halo_mass 			+= gal->comp_bool[i]*gal->comp_mass[i];
		if(gal->comp_type[i]==2) disk_mass 			+= gal->comp_bool[i]*gal->comp_mass[i];
		if(gal->comp_type[i]==3) bulge_mass 		+= gal->comp_bool[i]*gal->comp_mass[i];
		// Checking total cutted mass
		if(gal->comp_type[i]==0) cutted_gas_mass 	+= gal->comp_bool[i]*gal->comp_cutted_mass[i];
		if(gal->comp_type[i]==1) cutted_halo_mass 	+= gal->comp_bool[i]*gal->comp_cutted_mass[i];
		if(gal->comp_type[i]==2) cutted_disk_mass 	+= gal->comp_bool[i]*gal->comp_cutted_mass[i];
		if(gal->comp_type[i]==3) cutted_bulge_mass 	+= gal->comp_bool[i]*gal->comp_cutted_mass[i];
		// Checking particle number
		if(gal->comp_type[i]==0) gal->num_part[0] 	+= gal->comp_npart[i];
		if(gal->comp_type[i]==1) gal->num_part[1] 	+= gal->comp_npart[i];
		if(gal->comp_type[i]==2) gal->num_part[2] 	+= gal->comp_npart[i];
		if(gal->comp_type[i]==3) gal->num_part[3] 	+= gal->comp_npart[i];
		if(gal->comp_type[i]==0) gal->num_part_pot[0] 	+= gal->comp_npart_pot[i];
		if(gal->comp_type[i]==1) gal->num_part_pot[1] 	+= gal->comp_npart_pot[i];
		if(gal->comp_type[i]==2) gal->num_part_pot[2] 	+= gal->comp_npart_pot[i];
		if(gal->comp_type[i]==3) gal->num_part_pot[3] 	+= gal->comp_npart_pot[i];
		if(gal->comp_type[i]==0){
			// Set the gas specific internal energy
			// We assume the gas is in an isothermal state
			gal->comp_u_init[i] = (boltzmann / protonmass) * gal->comp_t_init[i];
			gal->comp_u_init[i] *= unit_mass_in_g / unit_energy;
			gal->comp_u_init[i] *= (1.0 / gamma_minus1);
			// Assuming full ionization (HII)
			/*if(gal->t_init > 1.0e4)
    		    molecular_weight = 4./(8.-5.*(1.-hydrogen_massfrac));
			// Assuming neutral gas (HI)
			else
    		    molecular_weight = 4./(1.+3.*hydrogen_massfrac);*/
    		molecular_weight = 4./(1.+3.*hydrogen_massfrac);
			gal->comp_u_init[i] /= molecular_weight;
			// Computing isothermal sound speed
			gal->comp_cs_init[i] = sqrt(gamma_minus1*gal->comp_u_init[i]*unit_energy/unit_mass_in_g);
		}
	}
    
    baryonic_fraction 				= (disk_mass+gas_mass+bulge_mass)/gal->total_mass;
    BT_fraction 				= cutted_bulge_mass/(cutted_bulge_mass+cutted_disk_mass);
    BD_fraction 				= cutted_bulge_mass/cutted_disk_mass;
    gas_fraction				= cutted_gas_mass/(cutted_disk_mass+cutted_gas_mass);
    

    
	// Print some information to screen.
	if (info != 0) {
		printf("/////\t-Redshift = %6.2lf\n",gal->redshift);
		printf("/////\t-V200 = %6.2lf km/s\n",gal->v200/1E5);
		printf("/////\t-R200 = %.3lf kpc\n",gal->r200);
		printf("/////\t-M200 = %6.2lf*1E10 solar mass\n",gal->m200/unit_mass);
		printf("/////\t-M200 [effective] = %6.2lf*1E10 solar mass\n",gal->m200*effective_mass_factor/unit_mass);
		printf("/////\t\t-Disk mass: \t%10.3lf*1E10 solar mass\n",disk_mass/unit_mass);
		printf("/////\t\t-Halo mass: \t%10.3lf*1E10 solar mass\n",halo_mass/unit_mass);
		printf("/////\t\t-Gas mass: \t%10.3lf*1E10 solar mass\n",gas_mass/unit_mass);
		printf("/////\t\t-Bulge mass: \t%10.3lf*1E10 solar mass\n",bulge_mass/unit_mass);
		printf("/////\t\t-Total stellar mass: \t%10.3lf*1E10 solar mass\n",(bulge_mass+disk_mass)/unit_mass);
		printf("/////\t-Total cutted mass = %6.2lf*1E10 solar mass.\n",gal->total_mass/unit_mass);
		printf("/////\t\t-Cutted Disk mass: \t%10.3lf*1E10 solar mass\n",cutted_disk_mass/unit_mass);
		printf("/////\t\t-Cutted Halo mass: \t%10.3lf*1E10 solar mass\n",cutted_halo_mass/unit_mass);
		printf("/////\t\t-Cutted Gas mass: \t%10.3lf*1E10 solar mass\n",cutted_gas_mass/unit_mass);
		printf("/////\t\t-Cutted Bulge mass: \t%10.3lf*1E10 solar mass\n",cutted_bulge_mass/unit_mass);
		printf("/////\t\t-Total cutted stellar mass: \t%10.3lf*1E10 solar mass\n",(cutted_bulge_mass+cutted_disk_mass)/unit_mass);
		printf("/////\t-%ld particles:\n",gal->ntot_part);
		printf("/////\t\t-%ld Gas particles\n",gal->num_part[0]);
		printf("/////\t\t-%ld Dark matter particles\n",gal->num_part[1]);
		printf("/////\t\t-%ld Stellar disk particles\n",gal->num_part[2]);
		printf("/////\t\t-%ld Stellar bulge particles\n",gal->num_part[3]);
		printf("/////\t-%ld particles for potential computation:\n",gal->ntot_part_pot);
		printf("/////\t\t-%ld Gas particles\n",gal->num_part_pot[0]);
		printf("/////\t\t-%ld Dark matter particles\n",gal->num_part_pot[1]);
		printf("/////\t\t-%ld Stellar disk particles\n",gal->num_part_pot[2]);
		printf("/////\t\t-%ld Stellar bulge particles\n",gal->num_part_pot[3]);
		printf("/////\t-Baryonic fraction:\t%10.3lf\n",baryonic_fraction);
		printf("/////\t-BD ratio:\t%10.3lf\n",BD_fraction);
		printf("/////\t-BT ratio:\t%10.3lf\n",BT_fraction);
		printf("/////\t-Gas fraction:\t%10.3lf\n",gas_fraction);
		printf("/////\t-Potential PM-grid has dimensions [nx=%d,ny=%d,nz=%d] \n",gal->ngrid,gal->ngrid,gal->ngrid);
		printf("/////\t-Dimensions are zero-padded to [nx=%d,ny=%d,nz=%d] \n",gal->ngrid_padded,gal->ngrid_padded,gal->ngrid_padded);
	}
	fflush(stdout);

	// Allocate particle id numbers array.
	if (!(gal->id=calloc(gal->ntot_part_pot,sizeof(unsigned long int)))) {
		fprintf(stderr,"Unable to allocate particle ID numbers.\n");
		return -1;
	}

	// Allocate x coordinates for all the particles.
	if (!(gal->x=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle x coordinates.\n");
		return -1;
	}
    
	// Allocate y coordinates for all the particles.
	if (!(gal->y=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle y coordinates.\n");
		return -1;
	}

	// Allocate z coordinates for all the particles.
	if (!(gal->z=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle z coordinates.\n");
		return -1;
	}
    
	// Allocate cylindrical radius for all the particles.
	if (!(gal->r_cyl=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle z coordinates.\n");
		return -1;
	}
    
	// Allocate cylindrical azimuthal angle for all the particles.
	if (!(gal->theta_cyl=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle z coordinates.\n");
		return -1;
	}
    
	// Allocate spherical radius for all the particles.
	if (!(gal->r_sph=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle z coordinates.\n");
		return -1;
	}
    
	// Allocate spherical azimuthal angle for all the particles.
	if (!(gal->theta_sph=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle z coordinates.\n");
		return -1;
	}
    
	// Allocate spherical polar angle for all the particles.
	if (!(gal->phi_sph=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle z coordinates.\n");
		return -1;
	}
    
	// Allocate x velocities for all the particles.
	if (!(gal->vel_x=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle x coordinates.\n");
		return -1;
	}
    
	// Allocate y velocities for all the particles.
	if (!(gal->vel_y=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle y coordinates.\n");
		return -1;
	}
    
	// Allocate z velocities for all the particles.
	if (!(gal->vel_z=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle z coordinates.\n");
		return -1;
	}
    
	// Allocate masses for all the particles.
	if (!(gal->mass=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle masses.\n");
		return -1;
	}
    
	// Allocate internal energies for all the particles.
	if (!(gal->u=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle internal energy.\n");
		return -1;
	}
	
	// Allocate metallicity array for particles.
	if (!(gal->metal=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle metal.\n");
		return -1;
	}
	
	// Allocate age array for all particles.
	if (!(gal->age=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle age.\n");
		return -1;
	}
	
	// Allocate density array for all particles.
	if (!(gal->rho=calloc(gal->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle age.\n");
		return -1;
	}
    
	// Turn on and allocate the potential grid, starting with x-axis
	gal->potential_defined = 1;
	if (!(gal->potential=calloc(gal->ngrid_padded,sizeof(double *)))) {
		fprintf(stderr,"Unable to create potential x axis.\n");
		return -1;
	}
    
	for (i = 0; i < gal->ngrid_padded; ++i) {
		// y-axis
		if (!(gal->potential[i] = calloc(gal->ngrid_padded,sizeof(double *)))) {
			fprintf(stderr,"Unable to create potential y axis.\n");
			return -1;
		}
		// z-axis
		for (j = 0; j < gal->ngrid_padded; ++j) {
			if (!(gal->potential[i][j] = calloc(gal->ngrid_padded,sizeof(double)))) {
				fprintf(stderr,"Unable to create potential z axis.\n");
				return -1;
			}
		}
	}
    
	// Allocate the midplane density grid, starting with x-axis
	if (!(gal->midplane_dens=calloc(gal->ngrid_dens,sizeof(double *)))) {
		fprintf(stderr,"Unable to create midplane_dens x axis.\n");
		return -1;
	}
    
	for (i = 0; i < gal->ngrid_dens; ++i) {
		// y-axis
		if (!(gal->midplane_dens[i] = calloc(gal->ngrid_dens,sizeof(double *)))) {
			fprintf(stderr,"Unable to create midplane_dens y axis.\n");
			return -1;
		}
	}
    
	// Allocate index all the threads.
	if (!(gal->index=calloc(AllVars.Nthreads,sizeof(unsigned long int)))) {
		fprintf(stderr,"Unable to allocate particle index.\n");
		return -1;
	}
	

    
    for (j=0; j<AllVars.MaxCompNumber; j++) {
    	if(gal->comp_npart[j]>0) printf("/////\tComponent %d -> particle mass %.2e solar mass\n",j+1,(gal->comp_cutted_mass[j]*1.0E10)/(gal->comp_npart_pot[j]*unit_mass));
		// Filling the arrays of the &galaxy structure
		for (i = gal->comp_start_part[j]; i < gal->comp_start_part[j] + gal->comp_npart_pot[j]; i++) {
			gal->mass[i] 		= gal->comp_cutted_mass[j]/gal->comp_npart_pot[j];
			gal->id[i] 			= i;
			gal->u[i] 			= gal->comp_u_init[j];
			gal->rho[i]			= 0.;
			gal->metal[i]		= gal->comp_metal[j];
			gal->age[i]			= (2*gal->comp_mean_age[j]-gal->comp_min_age[j])*gsl_rng_uniform_pos(r[0])+gal->comp_min_age[j];
		}
	}
	// No problems detected
	return 0;
}


// Use this to allocate the galaxy storage variable on the fly.
void allocate_galaxy_storage_variable(galaxy *gal, int size) {
	int i;
	
	// Allocate the storage array.
	if (!(gal->storage1=calloc(size,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate storage array.\n");
		return;
	}
	for (i = 0; i < size; ++i) {
		// y-axis
		if (!(gal->storage1[i] = calloc(AllVars.Nthreads,sizeof(double *)))) {
			fprintf(stderr,"Unable to allocate storage array.\n");
			return ;
		}
	}
	return;
}

// Set the positions of the particles in the galaxy to be consistent with the
// various particle distributions. This function uses the Metropolis algorithm
// to allocate particle positions from the distribution functions of each
// galaxy component.
int set_galaxy_coords(galaxy *gal) {
	int i,j;
	
	printf("/////\tSetting particles position\n");
	fflush(stdout);

	for(i=0; i<AllVars.MaxCompNumber; i++) {
		if(gal->comp_npart[i]>0) {
			printf("/////\t\t-Setting component %d particles positions ",i+1);
			fflush(stdout);
			mcmc_metropolis_hasting(gal,i,gal->comp_model[i]);
			rotate_component(gal,gal->comp_theta_sph[i],gal->comp_phi_sph[i],i);
		}
	}
	
	if (AllVars.GasHydrostaticEq) {
		if ((i = set_hydro_equilibrium(gal,AllVars.GasHydrostaticEqIter)) != 0) {
			fprintf(stderr,"Unable to compute azimuthal profile to reach hydro equilibrium. Aborting.\n");
			exit(0);
		}
	}
	// Be nice to the memory
	return 0;
}

// This function modify the z component of the gaseous particles according to an iterative
// algorithm to reach the hydrostatic equilibrium of the disk. This function was inspired by
// Springel, Di Matteo er al. 2005 method.
int set_hydro_equilibrium(galaxy *gal, int n_iter) {
	
	unsigned long int i,j,k;
	double mu, z0, pi_x, pi_y, q_x, q_y, prob, *radius, delta_pot;
	double theta, phi, randval, x, y, z, step, proposal, acceptance;
    
	printf("/////\tTargeting gas azimuthal hydrostatic equilibrium:\n");

	fflush(stdout);
	// Allow to use pseudo density functions
	for(i=0;i<AllVars.Nthreads;i++) gal->pseudo[i] = 1;
	for(k=0; k<AllVars.MaxCompNumber; k++) {
		// Looking for gas particles
		if(gal->comp_type[k]==0 && gal->comp_npart_pot[k]>0) {
			printf("/////\t\tRecomputing gas particles position for component %d\n",k);
			// Starting equilibrium iterations
			for(j = 0; j < n_iter; ++j) {
				printf("/////\t\t\tIteration %ld",j+1);
				fflush(stdout);
				// Now that we've got gas particles positions, we can compute the full potential.
				if(set_galaxy_potential(gal,0) != 0) {
					fprintf(stderr,"\nUnable to set the potential. Aborting.\n");
					return -1;
				}
				// Now, we have to set the gaseous disk density at (0.,0.,0.).
				// We have to make a guess: we use the disk density function to compute a realistic rho_gas(r=0,z=0).
				z0 = gal->comp_scale_height[k];
				// The MCMC algorithm is repeated 10 times to ensure convergence
				// Filling the Markov Chain
				acceptance = 0;
				fill_midplane_dens_grid(gal,k);
				gal->index[0] = gal->comp_start_part[k];
				gal->z[gal->comp_start_part[k]] = 0.0;
				// Recomputing density cut
				//gal->comp_cut_dens[k] = pseudo_density_gas_func(gal,0.99*gal->comp_cut[k],0.0,0.0,0,k);
				for(i = gal->comp_start_part[k]+1; i < gal->comp_start_part[k] + gal->comp_npart_pot[k]; ++i) {
					// Generating a proposal
					// Step is chosen as if the azimuthal distribution is sech^2(z/h)
					// Then there is an analytical solution to set hydro equilibrium
					// h = sqrt(cs^2/2*pi*G*rho0)
					// The Gaussian step is scaled by gal->comp_mcmc_step_hydro[k]
					x = gal->x[i];
					y = gal->y[i];
					z = gal->z[i-1];
					step = gal->comp_mcmc_step_hydro[k]*sqrt(pow(gal->comp_cs_init[k],2.0)/(2.0*pi*G*pseudo_density_gas_func(gal,x,y,0.0,0,k)))/kpc;
					// Computing the probability of the proposal
					// according to the updated potential value
					proposal 	= z + gsl_ran_gaussian(r[0],step);
					pi_x 		= pseudo_density_gas_func(gal,x,y,z,0,k);
					pi_y 		= pseudo_density_gas_func(gal,x,y,proposal,0,k);
					q_x  		= gsl_ran_gaussian_pdf(proposal-z,step);
					q_y  		= gsl_ran_gaussian_pdf(z-proposal,step);
					prob 		= min(1.0,(pi_y/pi_x)*(q_x/q_y));
					// Draw a uniform random value to test the acceptance of the proposal
					randval = gsl_rng_uniform_pos(r[0]);
					if (randval <= prob) {
						gal->z[i] = proposal;
						acceptance += 1.0;
						// Updating coordinates values
						gal->r_sph[i]		= sqrt(gal->x[i]*gal->x[i]+gal->y[i]*gal->y[i]+gal->z[i]*gal->z[i]);
						gal->phi_sph[i]		= acos(gal->z[i]/gal->r_sph[i]);
					} else {
						gal->z[i] 			= z;
						gal->r_sph[i]       = sqrt(gal->x[i]*gal->x[i]+gal->y[i]*gal->y[i]+gal->z[i]*gal->z[i]);
						gal->phi_sph[i]     = acos(gal->z[i]/gal->r_sph[i]);
					}
				}
				acceptance /= gal->comp_npart_pot[k];
				printf(" -> Acceptance = %lf\n",acceptance);
				if(acceptance<0.50) printf("/////\t\t\tWarning: MCMC acceptance is low!\n\t\t\tLower mcmc_step in the galaxy parameter file.\n");
				if(acceptance>0.90) printf("/////\t\t\tWarning: MCMC acceptance is high!\n\t\t\tIncrease mcmc_step in the galaxy parameter file.\n");
			}
			// Final gas midplane density
			fill_midplane_dens_grid(gal,k);
		}
	}
	return 0;
}

// Set the velocities of the particles for each component in the galaxy.
// User can choose between cold/hot particles, ie. without velocity dispersion/with velocity dispersion
// This computation is inspired from the work of Springel et al. 2005 & Hernquist 1993
int set_galaxy_velocity(galaxy *gal, int vel_type) {
    
	unsigned long int i,j;
	int status,warning,k;
	double v_c, v_theta, v_r, v_z;
	double f_s, sigma_theta;
	double v2a_r, v2a_theta, v2_theta, va_theta, vel_x, vel_y, vel_z, streaming_fraction;
	double vesc[4*gal->ngrid], mass_in_shell;
    
	// Warning init
	warning = 0;
    
	// Here we compute the velocity of the particles
	// assuming a Gaussian shaped velocity distribution function.
	// This method is much more realistic, and ensures disk stability
	// against axisymetric perturbations.
	printf("/////\tSetting velocities\n");
	fflush(stdout);
	// We fill an array with the value of the escape velocity sampled at radii up to boxsize    
	for(i=0;i<4*gal->ngrid;i++) {
		mass_in_shell = 0.;
		for (j=0;j<gal->ntot_part;++j) if(gal->r_sph[j]<(i+1)*gal->space[0]/4.) mass_in_shell+=gal->mass[j];
		vesc[i] = sqrt(2.0*G*mass_in_shell/((i+1)*gal->space[0]/4.*kpc));
	}

    // Loop over components
	for(j=0; j<AllVars.MaxCompNumber; j++) {
		// Particle velocities.
		if(gal->comp_npart[j]>0) {
			printf("/////\t\t-Setting component %d particles velocities\n",j+1);
			fflush(stdout);
			#pragma omp parallel for private(v2a_r,v2a_theta,sigma_theta,va_theta,v_c,v_r,v_theta,v_z,vel_x,vel_y,vel_z)
			for (i = gal->comp_start_part[j]; i < gal->comp_start_part[j]+gal->comp_npart[j]; ++i) {
				// Get thread ID
				#if USE_THREADS == 1
            		int tid = omp_get_thread_num();
				#else
            		int tid = 0;
				#endif
				double v2a_z;
				gal->index[tid] = i;

				
				// Specific case of gas particles
				if(gal->comp_type[j]==0) {
					if(AllVars.GasHydrostaticEq) gal->pseudo[tid] = 1; 
					v2_theta = v2_theta_gas_func(gal,gal->r_cyl[i],gal->z[i],j);
					// No velocity dispersion is needed because we are dealing with a fluid
					// ie. collisional particles
					v_r 	= 0.;
					v_theta = gal->comp_streaming_fraction[j]*sqrt(v2_theta);
					v_z 	= 0.;

					//Make sure to divide by 1.0E5 to put the velocities in km/s.
					gal->vel_x[i] = (v_r*cos(gal->theta_cyl[i])-v_theta*sin(gal->theta_cyl[i]))/1.0E5;
					gal->vel_y[i] = (v_r*sin(gal->theta_cyl[i])+v_theta*cos(gal->theta_cyl[i]))/1.0E5;
					gal->vel_z[i] = v_z/1.0E5;
					
				} else {
					v_c 		= v_c_func(gal,fabs(gal->r_cyl[i]));
					v2a_z 		= v2a_z_func(gal,w[tid],j);
					// Enforce Q>Q_min
					if(gal->Q_min>0 && gal->comp_type[j]==1) v2a_z = v2a_z_toomre(gal,fabs(gal->r_cyl[i]),v2a_z,j);
					v2a_r = v2a_z;
            	
					if(gal->comp_type[j]==2 && gal->DispExtCoeff>0. && gal->DispExtCoeff<1.0) {
						double disp_softening = (1.0-exp(-fabs(gal->r_cyl[i])/(-(1.0/log(1.0-gal->DispExtCoeff))*gal->comp_scale_length[j])));
						v2a_r *= disp_softening;
						//v2a_z = v2a_r;
					}
            
					if(gal->comp_type[j]==2 && gal->axisymmetric_drift==1) {
						va_theta = gal->comp_streaming_fraction[j]*v_c;
						// The angular velocity dispersion is attenuated radially
						// Because there is a singularity at r=0
						// Using the function: new_disp(scalelength) = 0.80*disp at the disk scalelenght
						sigma_theta = sqrt(sigma2_theta_disk_func(gal,fabs(gal->r_cyl[i]),v2a_z));
					} else {
						va_theta 		= gal->comp_streaming_fraction[j]*v_c;
						v2a_theta 		= v2a_r + v2a_theta_func(gal,gal->r_cyl[i],j) + v_c*v_c;
						// Check if the dispersion is a real or complex number
						if(AllVars.AcceptImaginary==1) {
							sigma_theta = sqrt(fabs(v2a_theta - va_theta*va_theta));
                		} else {
							sigma_theta = (v2a_theta>va_theta*va_theta)?sqrt(v2a_theta-va_theta*va_theta):0.;
                		}
					}
					v_r 	= gsl_ran_gaussian(r[tid],sqrt(v2a_r));
					v_theta = gsl_ran_gaussian(r[tid],sigma_theta)+va_theta;
					v_z 	= gsl_ran_gaussian(r[tid],sqrt(v2a_z));

					// Make sure to divide by 1.0E5 to put the velocities in km/s.
					vel_x = (v_r*cos(gal->theta_cyl[i])-v_theta*sin(gal->theta_cyl[i]))/1.0E5;
					vel_y = (v_r*sin(gal->theta_cyl[i])+v_theta*cos(gal->theta_cyl[i]))/1.0E5;
					vel_z = v_z/1.0E5;
				
				
					// Escape velocity
					int vesc_index = (int)(gal->r_sph[i]/(gal->space[0]/4.));
            	
					// Let's ensure that the particle velocity is lower than gal->disk_vmax times the escape velocity
					int ct = 0;
					//printf("%le %lf %lf \n",gal->r_sph[i],sqrt(vel_x*vel_x+vel_y*vel_y+vel_z*vel_z),gal->disk_vmax*vesc[vesc_index]);
					/*while(sqrt(vel_x*vel_x+vel_y*vel_y+vel_z*vel_z) >= gal->disk_vmax*vesc[vesc_index]) {
            	 	if(ct >= 100) break;
            	 	v_r 	= gsl_ran_gaussian(r[tid],sqrt(v2a_r));
             		v_theta = gsl_ran_gaussian(r[tid],sigma_theta)+va_theta;
             		v_z 	= gsl_ran_gaussian(r[tid],sqrt(v2a_z));
             		vel_x 	= (v_r*cos(gal->theta_cyl[i])-v_theta*sin(gal->theta_cyl[i]))/1.0E5;
             		vel_y 	= (v_r*sin(gal->theta_cyl[i])+v_theta*cos(gal->theta_cyl[i]))/1.0E5;
             		vel_z 	= v_z/1.0E5;
             		ct++;
					}*/
					ct = 0;
					while(fabs(v_r) > gal->comp_vmax[j]*v_c/*vesc[vesc_index]*/ && gal->r_cyl[i]>gal->comp_scale_length[j]) {
						if(ct >= 10000) {
							v_r = 2.0*gal->comp_vmax[j]*v_c/*vesc[vesc_index]*/*(gsl_rng_uniform_pos(r[tid])-0.5);
							break;
						}
						v_r 	= gsl_ran_gaussian(r[tid],sqrt(v2a_r));
						ct++;
					}
					ct = 0;
					while(fabs(v_theta-va_theta) > gal->comp_vmax[j]*v_c/*vesc[vesc_index]*/ && gal->r_cyl[i]>gal->comp_scale_length[j]) {
						if(ct >= 10000) {
							v_theta = 2.0*gal->comp_vmax[j]*v_c/*vesc[vesc_index]*/*(gsl_rng_uniform_pos(r[tid])-0.5)+va_theta;
							break;
						}
						v_theta = gsl_ran_gaussian(r[tid],sigma_theta)+va_theta;
						ct++;
					}
					ct = 0;
					while(fabs(v_z) > gal->comp_vmax[j]*v_c/*vesc[vesc_index]*/ && gal->r_cyl[i]>gal->comp_scale_length[j]) {
						if(ct >= 10000) {
							v_z = 2.0*gal->comp_vmax[j]*v_c/*vesc[vesc_index]*/*(gsl_rng_uniform_pos(r[tid])-0.5);
							break;
						}
						v_z = gsl_ran_gaussian(r[tid],sqrt(v2a_z));
						ct++;
					}
            
					vel_x 	= (v_r*cos(gal->theta_cyl[i])-v_theta*sin(gal->theta_cyl[i]))/1.0E5;
					vel_y 	= (v_r*sin(gal->theta_cyl[i])+v_theta*cos(gal->theta_cyl[i]))/1.0E5;
					vel_z 	= v_z/1.0E5;
            	
					gal->vel_x[i] = vel_x;
					gal->vel_y[i] = vel_y;
					gal->vel_z[i] = vel_z;
				}
			}
			//if(gal->comp_type[j]==2) printf(" --> Toomre parameter minimum: Q = %f\n",gal->Q_min);
		}
	}

	if(warning) printf("/////\t\tWarning: The potential derivative seems unstable. Increase particule number.\n");
	// Be nice to memory
	free(gal->storage1);
    
	if(AllVars.OutputRc==1){
		double maxrad;
		maxrad = 0.;
		for(j=0; j<AllVars.MaxCompNumber; j++) if(gal->comp_cut[j]>maxrad) maxrad = gal->comp_cut[j];
		printf("/////\t\t Writing rotation curve\n");
		write_galaxy_rotation_curve(gal,maxrad);
	}
	return 0;
}

// Be kind to the memory, destroy a galaxy...
void destroy_galaxy(galaxy *gal, int info) {
	
	int i,j;
	
	if (info != 0) printf("/////Destroying galaxy...");
	// Deallocate the potential grid to be really nice to the memory.
	if(gal->potential_defined == 1){
		for (i = 0; i < gal->ngrid_padded; i++){
			for (j = 0; j < gal->ngrid_padded; j++){
				free(gal->potential[i][j]);
			}
			free(gal->potential[i]);
		}
		free(gal->potential);
	}
	// Deallocate all the small parts of the galaxy to be nice to the memory.
	free(gal->id);
	free(gal->x);
	free(gal->y);
	free(gal->z);
	free(gal->mass);
	free(gal->vel_x);
	free(gal->vel_y);
	free(gal->vel_z);
	if (info != 0) printf("///// All memory unallocated.\n");
	return;
}


// Destroy a system of galaxies.
void destroy_galaxy_system(galaxy *gal, int info) {
	
	if (info != 0) printf("/////Destroying galaxy system...");
	free(gal->id);
	free(gal->x);
	free(gal->y);
	free(gal->z);
	free(gal->mass);
	free(gal->vel_x);
	free(gal->vel_y);
	free(gal->vel_z);
	if (info != 0) printf("///// All memory unallocated.\n");
    return;
}


// This function intend to perform rotation onto a galaxy, in order to test
// different combination of ICs for galaxy colision simulations.
// The orientation angles should be specified in degrees.
int rotate_galaxy(galaxy *gal, double alpha, double delta) {
	
	unsigned long int i;
	double x_temp, y_temp, z_temp;
	double vx_temp, vy_temp, vz_temp;
	// Alpha is the spin angle
	// Delta is the inclination of the disk compare to the XY plane
	alpha = alpha*pi/180.;
	delta = delta*pi/180.;
	for(i = 0; i<gal->ntot_part_pot; ++i) {
		//Rotation around Y axis
        x_temp			= cos(delta)*gal->x[i]+sin(delta)*gal->z[i];
        z_temp			= cos(delta)*gal->z[i]-sin(delta)*gal->x[i];
        vx_temp			= cos(delta)*gal->vel_x[i]+sin(delta)*gal->vel_z[i];
        vz_temp			= cos(delta)*gal->vel_z[i]-sin(delta)*gal->vel_x[i];
		gal->x[i] 		= x_temp;
        gal->z[i] 		= z_temp;
		gal->vel_x[i] 	= vx_temp;
		gal->vel_z[i] 	= vz_temp;
        //Rotation around Z axis
        x_temp			= cos(alpha)*gal->x[i]+sin(alpha)*gal->y[i];
        y_temp			= cos(alpha)*gal->y[i]-sin(alpha)*gal->x[i];
        vx_temp			= cos(alpha)*gal->vel_x[i]+sin(alpha)*gal->vel_y[i];
        vy_temp			= cos(alpha)*gal->vel_y[i]-sin(alpha)*gal->vel_x[i];
        gal->x[i] 		= x_temp;
        gal->y[i] 		= y_temp;
		gal->vel_x[i] 	= vx_temp;
		gal->vel_y[i] 	= vy_temp;
	}
	return 0;
}


// This function intend to perform rotation onto a galaxy, in order to test
// different combination of ICs for galaxy colision simulations.
// The orientation angles should be specified in degrees.
int rotate_component(galaxy *gal, double alpha, double delta, int component) {
	
	unsigned long int i;
	double x_temp, y_temp, z_temp;
	double vx_temp, vy_temp, vz_temp;
	// Alpha is the spin angle
	// Delta is the inclination of the disk compare to the XY plane
	alpha = alpha*pi/180.;
	delta = delta*pi/180.;
	for(i = gal->comp_start_part[component]; i<gal->comp_start_part[component]+gal->comp_npart_pot[component]; ++i) {
		//Rotation around Y axis
        x_temp			= cos(delta)*gal->x[i]+sin(delta)*gal->z[i];
        z_temp			= cos(delta)*gal->z[i]-sin(delta)*gal->x[i];
        vx_temp			= cos(delta)*gal->vel_x[i]+sin(delta)*gal->vel_z[i];
        vz_temp			= cos(delta)*gal->vel_z[i]-sin(delta)*gal->vel_x[i];
		gal->x[i] 		= x_temp;
        gal->z[i] 		= z_temp;
		gal->vel_x[i] 	= vx_temp;
		gal->vel_z[i] 	= vz_temp;
        //Rotation around Z axis
        x_temp			= cos(alpha)*gal->x[i]+sin(alpha)*gal->y[i];
        y_temp			= cos(alpha)*gal->y[i]-sin(alpha)*gal->x[i];
        vx_temp			= cos(alpha)*gal->vel_x[i]+sin(alpha)*gal->vel_y[i];
        vy_temp			= cos(alpha)*gal->vel_y[i]-sin(alpha)*gal->vel_x[i];
        gal->x[i] 		= x_temp;
        gal->y[i] 		= y_temp;
		gal->vel_x[i] 	= vx_temp;
		gal->vel_y[i] 	= vy_temp;
	}
	return 0;
}


// This function intend to setup the gloabal postion velocity vectors of a galaxy, in order to test
// different combination of ICs for galaxy colision simulations.
int set_galaxy_trajectory(galaxy *gal) {
	unsigned long int i;
    
	for(i = 0; i<gal->ntot_part_pot; ++i) {
		gal->x[i] 		+= gal->xc;
		gal->y[i] 		+= gal->yc;
        gal->z[i] 		+= gal->zc;
		gal->vel_x[i] 	+= gal->vel_xc;
		gal->vel_y[i] 	+= gal->vel_yc;
		gal->vel_z[i] 	+= gal->vel_zc;
	}
	return 0;
}

// In order to perform a galaxy collision, it is necessary to combine
// two galaxies into one set of orbiting galaxies. The orbits should
// be defined by a user-defined orbit or the set_orbit_parabolic()
// function and the following function should be called to combine
// the galaxies into one object.
//
// It is necessary to allocate the different parts of the memory here
// because this galactic system does not have quantities like halo
// or disk scale lengths. These are quantities in the parent galaxy,
// not the galaxy-orbit-galaxy combination.
int create_galaxy_system(galaxy *gal_1, galaxy *gal_2, galaxy *gal_3) {
	unsigned long int i,a;
	int j,k;

	// Allocate component arrays.
	if (!(gal_3->comp_mass=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mass array.\n");
		return -1;
	}
	if (!(gal_3->comp_model=calloc(AllVars.MaxCompNumber,sizeof(int)))) {
		fprintf(stderr,"Unable to allocate particle comp_model array.\n");
		return -1;
	}
	if (!(gal_3->comp_npart=calloc(AllVars.MaxCompNumber,sizeof(unsigned long int)))) {
		fprintf(stderr,"Unable to allocate particle comp_model array.\n");
		return -1;
	}
	if (!(gal_3->comp_npart_pot=calloc(AllVars.MaxCompNumber,sizeof(unsigned long int)))) {
		fprintf(stderr,"Unable to allocate particle comp_model array.\n");
		return -1;
	}
	if (!(gal_3->comp_start_part=calloc(AllVars.MaxCompNumber,sizeof(unsigned long int)))) {
		fprintf(stderr,"Unable to allocate particle comp_model array.\n");
		return -1;
	}
	if (!(gal_3->comp_cutted_mass=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle cutted_comp_mass array.\n");
		return -1;
	}
	if (!(gal_3->comp_scale_length=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_scale_length array.\n");
		return -1;
	}
	if (!(gal_3->comp_scale_height=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_scale_height array.\n");
		return -1;
	}
	if (!(gal_3->comp_cut=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_cut array.\n");
		return -1;
	}
	if (!(gal_3->comp_flat=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_flat array.\n");
		return -1;
	}
	if (!(gal_3->comp_mcmc_step=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mcmc_step array.\n");
		return -1;
	}
	if (!(gal_3->comp_mcmc_step_hydro=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mcmc_step array.\n");
		return -1;
	}
	if (!(gal_3->comp_vmax=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_vmax array.\n");
		return -1;
	}
	if (!(gal_3->comp_mass_frac=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mass_frac array.\n");
		return -1;
	}
	if (!(gal_3->comp_type=calloc(AllVars.MaxCompNumber,sizeof(int)))) {
		fprintf(stderr,"Unable to allocate particle comp_type array.\n");
		return -1;
	}
	if (!(gal_3->comp_bool=calloc(AllVars.MaxCompNumber,sizeof(int)))) {
		fprintf(stderr,"Unable to allocate particle comp_bool array.\n");
		return -1;
	}
	if (!(gal_3->comp_cut_dens=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_cut_dens array.\n");
		return -1;
	}
	if (!(gal_3->comp_concentration=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_concentration array.\n");
		return -1;
	}
	if (!(gal_3->comp_streaming_fraction=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_streaming_fraction array.\n");
		return -1;
	}
	if (!(gal_3->comp_theta_sph=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_theta_sph array.\n");
		return -1;
	}
	if (!(gal_3->comp_phi_sph=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_phi_sph array.\n");
		return -1;
	}
	if (!(gal_3->comp_metal=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_metal array.\n");
		return -1;
	}
	if (!(gal_3->comp_t_init=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_t_init array.\n");
		return -1;
	}
	if (!(gal_3->comp_u_init=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_u_init array.\n");
		return -1;
	}
	if (!(gal_3->comp_cs_init=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_cs_init array.\n");
		return -1;
	}
	if (!(gal_3->comp_mean_age=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mean_age array.\n");
		return -1;
	}
	if (!(gal_3->comp_min_age=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_min_age array.\n");
		return -1;
	}
	// Create the galaxy system first!
	j=0;
	for(k=0;k<AllVars.MaxCompNumber;k++) {
		if(gal_1->comp_npart[k]>0) {
			gal_3->comp_npart[j]			= gal_1->comp_npart[k];
			gal_3->comp_npart_pot[j]		= gal_1->comp_npart_pot[k];
			gal_3->comp_start_part[j]		= gal_1->comp_start_part[k];
			gal_3->comp_mass_frac[j]		= gal_1->comp_mass_frac[k];
			gal_3->comp_mass[j] 			= gal_1->comp_mass[k];
			gal_3->comp_model[j]			= gal_1->comp_model[k];
			gal_3->comp_cutted_mass[j]		= gal_1->comp_cutted_mass[k];
			gal_3->comp_scale_length[j]		= gal_1->comp_scale_length[k];
			gal_3->comp_concentration[j]	= gal_1->comp_concentration[k];
			gal_3->comp_scale_height[j] 	= gal_1->comp_scale_height[k];
			gal_3->comp_cut[j]             	= gal_1->comp_cut[k];
			gal_3->comp_flat[j]				= gal_1->comp_flat[k];
			gal_3->comp_mcmc_step[j]		= gal_1->comp_mcmc_step[k];
			gal_3->comp_mcmc_step_hydro[j]	= gal_1->comp_mcmc_step_hydro[k];
			gal_3->comp_vmax[j] 			= gal_1->comp_vmax[k];
			gal_3->comp_type[j] 			= gal_1->comp_type[k];
			gal_3->comp_bool[j] 			= gal_1->comp_bool[k];
			gal_3->comp_streaming_fraction[j] 	= gal_1->comp_streaming_fraction[k];
			gal_3->comp_cut_dens[j]				= gal_1->comp_cut_dens[k];
			gal_3->comp_theta_sph[j]			= gal_1->comp_theta_sph[k];
			gal_3->comp_phi_sph[j]				= gal_1->comp_phi_sph[k];
			gal_3->comp_metal[j]				= gal_1->comp_metal[k];
			gal_3->comp_t_init[j]				= gal_1->comp_t_init[k];
			gal_3->comp_u_init[j]				= gal_1->comp_u_init[k];
			gal_3->comp_cs_init[j]				= gal_1->comp_cs_init[k];
			gal_3->comp_mean_age[j]				= gal_1->comp_mean_age[k];
			gal_3->comp_min_age[j]				= gal_1->comp_min_age[k];
			j++;
		}
	}
	for(k=0;k<AllVars.MaxCompNumber;k++) {
		if(gal_2->comp_npart[k]>0) {
			gal_3->comp_npart[j]			= gal_2->comp_npart[k];
			gal_3->comp_npart_pot[j]		= gal_2->comp_npart_pot[k];
			gal_3->comp_start_part[j]		= gal_1->ntot_part_pot+gal_2->comp_start_part[k];
			gal_3->comp_mass_frac[j]		= gal_2->comp_mass_frac[k];
			gal_3->comp_mass[j] 			= gal_2->comp_mass[k];
			gal_3->comp_model[j]			= gal_2->comp_model[k];
			gal_3->comp_cutted_mass[j]		= gal_2->comp_cutted_mass[k];
			gal_3->comp_scale_length[j]		= gal_2->comp_scale_length[k];
			gal_3->comp_concentration[j]	= gal_2->comp_concentration[k];
			gal_3->comp_scale_height[j] 	= gal_2->comp_scale_height[k];
			gal_3->comp_cut[j]             	= gal_2->comp_cut[k];
			gal_3->comp_flat[j]				= gal_2->comp_flat[k];
			gal_3->comp_mcmc_step[j]		= gal_2->comp_mcmc_step[k];
			gal_3->comp_mcmc_step_hydro[j]	= gal_2->comp_mcmc_step_hydro[k];
			gal_3->comp_vmax[j] 			= gal_2->comp_vmax[k];
			gal_3->comp_type[j] 			= gal_2->comp_type[k];
			gal_3->comp_bool[j] 			= gal_2->comp_bool[k];
			gal_3->comp_streaming_fraction[j] 	= gal_2->comp_streaming_fraction[k];
			gal_3->comp_cut_dens[j]				= gal_2->comp_cut_dens[k];
			gal_3->comp_theta_sph[j]			= gal_2->comp_theta_sph[k];
			gal_3->comp_phi_sph[j]				= gal_2->comp_phi_sph[k];
			gal_3->comp_metal[j]				= gal_2->comp_metal[k];
			gal_3->comp_t_init[j]				= gal_2->comp_t_init[k];
			gal_3->comp_u_init[j]				= gal_2->comp_u_init[k];
			gal_3->comp_cs_init[j]				= gal_2->comp_cs_init[k];
			gal_3->comp_mean_age[j]				= gal_2->comp_mean_age[k];
			gal_3->comp_min_age[j]				= gal_2->comp_min_age[k];
			j++;
		}
	}
	
    gal_3->ntot_part 		= gal_1->ntot_part+gal_2->ntot_part;
    gal_3->num_part[0] 		= gal_1->num_part[0]+gal_2->num_part[0];
    gal_3->num_part[1] 		= gal_1->num_part[1]+gal_2->num_part[1];
    gal_3->num_part[2] 		= gal_1->num_part[2]+gal_2->num_part[2];
    gal_3->num_part[3] 		= gal_1->num_part[3]+gal_2->num_part[3]; 
    
    gal_3->ntot_part_pot 	= gal_1->ntot_part_pot+gal_2->ntot_part_pot;
    gal_3->num_part_pot[0] 	= gal_1->num_part_pot[0]+gal_2->num_part_pot[0];
    gal_3->num_part_pot[1] 	= gal_1->num_part_pot[1]+gal_2->num_part_pot[1];
    gal_3->num_part_pot[2] 	= gal_1->num_part_pot[2]+gal_2->num_part_pot[2];
    gal_3->num_part_pot[3] 	= gal_1->num_part_pot[3]+gal_2->num_part_pot[3]; 
    
    // Allocate particle id numbers array.
    if (!(gal_3->id=calloc(gal_3->ntot_part_pot,sizeof(unsigned long int)))) {
        fprintf(stderr,"Unable to allocate particle ID numbers.\n");
        return -1;
	}
	// Allocate x coordinates for all the particles.
	if (!(gal_3->x=calloc(gal_3->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle x coordinates.\n");
        return -1;
	}
	// Allocate y coordinates for all the particles.
	if (!(gal_3->y=calloc(gal_3->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle y coordinates.\n");
        return -1;
	}
	// Allocate z coordinates for all the particles.
	if (!(gal_3->z=calloc(gal_3->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle z coordinates.\n");
        return -1;
	}
	// Allocate x velocities for all the particles.
	if (!(gal_3->vel_x=calloc(gal_3->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle x coordinates.\n");
        return -1;
	}
	// Allocate y velocities for all the particles.
    if (!(gal_3->vel_y=calloc(gal_3->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle y coordinates.\n");
        return -1;
    }
    // Allocate z velocities for all the particles.
    if (!(gal_3->vel_z=calloc(gal_3->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle z coordinates.\n");
        return -1;
    }
    // Allocate masses for all the particles.
    if (!(gal_3->mass=calloc(gal_3->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle masses.\n");
        return -1;
    }
    // Allocate internal energy for all the particles.
    if (!(gal_3->u=calloc(gal_3->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle internal energy.\n");
        return -1;
    }
	// Allocate metallicity array for particles.
	if (!(gal_3->metal=calloc(gal_3->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle metal.\n");
		return -1;
	}
	// Allocate age array for all particles.
	if (!(gal_3->age=calloc(gal_3->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle age.\n");
		return -1;
	}
	// Allocate density array for all particles.
	if (!(gal_3->rho=calloc(gal_3->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle age.\n");
		return -1;
	}
    // Turn off the galaxy potential.
    gal_3->potential_defined = 0;
		
    // Copy all of the galaxy information.
	a = 0;
    for (i = 0; i < gal_1->ntot_part_pot; ++i) {
        gal_3->mass[a] 		= gal_1->mass[i];
        gal_3->id[a] 		= a;
	    gal_3->x[a] 		= gal_1->x[i];
	    gal_3->y[a] 		= gal_1->y[i];
	    gal_3->z[a] 		= gal_1->z[i];
	    gal_3->vel_x[a] 	= gal_1->vel_x[i];
	    gal_3->vel_y[a] 	= gal_1->vel_y[i];
	    gal_3->vel_z[a] 	= gal_1->vel_z[i];
	    gal_3->u[a]			= gal_1->u[i];
	    gal_3->rho[a]		= gal_1->rho[i];
	    gal_3->metal[a]		= gal_1->metal[i];
	    gal_3->age[a]		= gal_1->age[i];
		a++;
    }
    for (i = 0; i < gal_2->ntot_part_pot; ++i) {
        gal_3->mass[a] 		= gal_2->mass[i];
	    gal_3->id[a] 		= a;
	    gal_3->x[a] 		= gal_2->x[i];
	    gal_3->y[a] 		= gal_2->y[i];
	    gal_3->z[a] 		= gal_2->z[i];
	    gal_3->vel_x[a] 	= gal_2->vel_x[i];
   	    gal_3->vel_y[a] 	= gal_2->vel_y[i];
   	    gal_3->vel_z[a] 	= gal_2->vel_z[i];
   	    gal_3->u[a]			= gal_1->u[i];
	    gal_3->rho[a]		= gal_1->rho[i];
	    gal_3->metal[a]		= gal_1->metal[i];
	    gal_3->age[a]		= gal_1->age[i];
		a++;
   	}
    return 0;
}

// Copy a galaxy strucuture into a new one
int copy_galaxy(galaxy *gal_1, galaxy *gal_2, int info) {
    
	unsigned long int i;
	int j;
	
	destroy_galaxy(gal_2,0);
	// Copy all the non-pointers structures variables
	*gal_2 = *gal_1;
		
	// Allocate component arrays.
	if (!(gal_2->comp_mass=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mass array.\n");
		return -1;
	}
	if (!(gal_2->comp_model=calloc(AllVars.MaxCompNumber,sizeof(int)))) {
		fprintf(stderr,"Unable to allocate particle comp_model array.\n");
		return -1;
	}
	if (!(gal_2->comp_npart=calloc(AllVars.MaxCompNumber,sizeof(unsigned long int)))) {
		fprintf(stderr,"Unable to allocate particle comp_model array.\n");
		return -1;
	}
	if (!(gal_2->comp_npart_pot=calloc(AllVars.MaxCompNumber,sizeof(unsigned long int)))) {
		fprintf(stderr,"Unable to allocate particle comp_model array.\n");
		return -1;
	}
	if (!(gal_2->comp_start_part=calloc(AllVars.MaxCompNumber,sizeof(unsigned long int)))) {
		fprintf(stderr,"Unable to allocate particle comp_model array.\n");
		return -1;
	}
	if (!(gal_2->comp_cutted_mass=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle cutted_comp_mass array.\n");
		return -1;
	}
	if (!(gal_2->comp_scale_length=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_scale_length array.\n");
		return -1;
	}
	if (!(gal_2->comp_scale_height=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_scale_height array.\n");
		return -1;
	}
	if (!(gal_2->comp_cut=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_cut array.\n");
		return -1;
	}
	if (!(gal_2->comp_flat=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_flat array.\n");
		return -1;
	}
	if (!(gal_2->comp_mcmc_step=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mcmc_step array.\n");
		return -1;
	}
	if (!(gal_2->comp_mcmc_step_hydro=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mcmc_step array.\n");
		return -1;
	}
	if (!(gal_2->comp_vmax=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_vmax array.\n");
		return -1;
	}
	if (!(gal_2->comp_mass_frac=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mass_frac array.\n");
		return -1;
	}
	if (!(gal_2->comp_type=calloc(AllVars.MaxCompNumber,sizeof(int)))) {
		fprintf(stderr,"Unable to allocate particle comp_type array.\n");
		return -1;
	}
	if (!(gal_2->comp_bool=calloc(AllVars.MaxCompNumber,sizeof(int)))) {
		fprintf(stderr,"Unable to allocate particle comp_bool array.\n");
		return -1;
	}
	if (!(gal_2->comp_cut_dens=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_cut_dens array.\n");
		return -1;
	}
	if (!(gal_2->comp_concentration=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_concentration array.\n");
		return -1;
	}
	if (!(gal_2->comp_streaming_fraction=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_streaming_fraction array.\n");
		return -1;
	}
	if (!(gal_2->comp_theta_sph=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_theta_sph array.\n");
		return -1;
	}
	if (!(gal_2->comp_phi_sph=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_phi_sph array.\n");
		return -1;
	}
	if (!(gal_2->comp_metal=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_metal array.\n");
		return -1;
	}
	if (!(gal_2->comp_t_init=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_t_init array.\n");
		return -1;
	}
	if (!(gal_2->comp_u_init=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_u_init array.\n");
		return -1;
	}
	if (!(gal_2->comp_cs_init=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_cs_init array.\n");
		return -1;
	}
	if (!(gal_2->comp_mean_age=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_mean_age array.\n");
		return -1;
	}
	if (!(gal_2->comp_min_age=calloc(AllVars.MaxCompNumber,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle comp_min_age array.\n");
		return -1;
	}
	
	for (j = 0; j < AllVars.MaxCompNumber; ++j) {
		gal_2->comp_npart[j]				= gal_1->comp_npart[j];
		gal_2->comp_npart_pot[j]			= gal_1->comp_npart_pot[j];
		gal_2->comp_start_part[j]			= gal_1->comp_start_part[j];
		gal_2->comp_mass_frac[j]			= gal_1->comp_mass_frac[j];
		gal_2->comp_mass[j] 				= gal_1->comp_mass[j];
		gal_2->comp_model[j]				= gal_1->comp_model[j];
		gal_2->comp_cutted_mass[j]			= gal_1->comp_cutted_mass[j];
		gal_2->comp_scale_length[j]			= gal_1->comp_scale_length[j];
		gal_2->comp_concentration[j]		= gal_1->comp_concentration[j];
		gal_2->comp_scale_height[j] 		= gal_1->comp_scale_height[j];
		gal_2->comp_cut[j]             		= gal_1->comp_cut[j];
		gal_2->comp_flat[j]					= gal_1->comp_flat[j];
		gal_2->comp_mcmc_step[j]			= gal_1->comp_mcmc_step[j];
		gal_2->comp_mcmc_step_hydro[j]		= gal_1->comp_mcmc_step_hydro[j];
		gal_2->comp_vmax[j] 				= gal_1->comp_vmax[j];
		gal_2->comp_type[j] 				= gal_1->comp_type[j];
		gal_2->comp_bool[j] 				= gal_1->comp_bool[j];
		gal_2->comp_streaming_fraction[j] 	= gal_1->comp_streaming_fraction[j];
		gal_2->comp_cut_dens[j]				= gal_1->comp_cut_dens[j];
		gal_2->comp_theta_sph[j]			= gal_1->comp_theta_sph[j];
		gal_2->comp_phi_sph[j]				= gal_1->comp_phi_sph[j];
		gal_2->comp_metal[j]				= gal_1->comp_metal[j];
		gal_2->comp_t_init[j]				= gal_1->comp_t_init[j];
		gal_2->comp_u_init[j]				= gal_1->comp_u_init[j];
		gal_2->comp_cs_init[j]				= gal_1->comp_cs_init[j];
		gal_2->comp_mean_age[j]				= gal_1->comp_mean_age[j];
		gal_2->comp_min_age[j]				= gal_1->comp_min_age[j];
	}
	// Allocate particle id numbers array.
    if (!(gal_2->id=calloc(gal_2->ntot_part_pot,sizeof(unsigned long int)))) {
        fprintf(stderr,"Unable to allocate particle ID numbers.\n");
        return -1;
	}
	// Allocate x coordinates for all the particles.
	if (!(gal_2->x=calloc(gal_2->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle x coordinates.\n");
        return -1;
	}
	// Allocate y coordinates for all the particles.
	if (!(gal_2->y=calloc(gal_2->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle y coordinates.\n");
        return -1;
	}
	// Allocate z coordinates for all the particles.
	if (!(gal_2->z=calloc(gal_2->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle z coordinates.\n");
        return -1;
	}
	// Allocate x velocities for all the particles.
	if (!(gal_2->vel_x=calloc(gal_2->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle x coordinates.\n");
        return -1;
	}
	// Allocate y velocities for all the particles.
    if (!(gal_2->vel_y=calloc(gal_2->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle y coordinates.\n");
        return -1;
    }
    // Allocate z velocities for all the particles.
    if (!(gal_2->vel_z=calloc(gal_2->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle z coordinates.\n");
        return -1;
    }
    // Allocate masses for all the particles.
    if (!(gal_2->mass=calloc(gal_2->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle masses.\n");
        return -1;
    }
    // Allocate internal energy for all the particles.
    if (!(gal_2->u=calloc(gal_2->ntot_part_pot,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle internal energy.\n");
        return -1;
    }
	// Allocate metallicity array for particles.
	if (!(gal_2->metal=calloc(gal_2->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle metal.\n");
		return -1;
	}
	// Allocate age array for all particles.
	if (!(gal_2->age=calloc(gal_2->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle age.\n");
		return -1;
	}
	// Allocate density array for all particles.
	if (!(gal_2->rho=calloc(gal_2->ntot_part_pot,sizeof(double)))) {
		fprintf(stderr,"Unable to allocate particle age.\n");
		return -1;
	}
    
	gal_2->potential_defined = 0;
	
    // Copy all the coordinate information.
    for (i = 0; i < gal_1->ntot_part_pot; ++i) {
		gal_2->x[i] 	= gal_1->x[i];
		gal_2->y[i] 	= gal_1->y[i];
		gal_2->z[i] 	= gal_1->z[i];
		gal_2->vel_x[i] = gal_1->vel_x[i];
		gal_2->vel_y[i] = gal_1->vel_y[i];
		gal_2->vel_z[i] = gal_1->vel_z[i];
		gal_2->mass[i] 	= gal_1->mass[i];
		gal_2->u[i] 	= gal_1->u[i];
	}
	return 0;
}


// This function places two galaxies on a Keplerian orbit in the x,y
// plane. The orbit is fully parametrizable going from hyperbolic orbits to elliptic orbits,
// including the special case of parabolic orbits.
void set_orbit_keplerian(galaxy *gal_1, galaxy *gal_2, double sep, double per, double e, double phi, double theta) {
    
	unsigned long int i;
	double nu1, nu2, x_1, y_1, z_1, x_2, y_2, z_2, vx_1, vy_1, vz_1, vx_2, vy_2, vz_2, mu, angmom;
	double radius1,radius2,a1,a2,l1,l2,k1,k2,v_ini,r_ini,specific_orbital_energy,kappa,l,theta_diago, degtorad;
	double xt_1,yt_1,zt_1,vxt_1,vyt_1,vzt_1,xt_2,yt_2,zt_2,vxt_2,vyt_2,vzt_2;
   	double m1,m2;
	// Degrees to radian conversion factor
	degtorad = pi/180.;
	// Mu is the standard gravitational parameter, calculated here in
	// terms of the reduced mass.
	m1 = gal_1->total_mass/unit_mass;
	m2 = gal_2->total_mass/unit_mass;
	mu = (m1*m2)/(m1+m2);
	if(per > sep) {
		fprintf(stderr,"Warning: The pericentral distance must be lower than the initial distance.");
		fprintf(stderr,"Warning: Setting initial distance value equal to pericentral distance value.");
		sep = per;
	}
	// Mass ratios
	k1 = mu/m1;
	k2 = mu/m2;
	// Keeping the center of mass at [0,0,0]
	radius1 = k1*sep;
	radius2 = k2*sep;
	// Semi-major axis
	if( e == 1 ) a1 = k1*per; else a1 = k1*per/fabs(1.0-e);
	if( e == 1 ) a2 = k2*per; else a2 = k2*per/fabs(1.0-e);
	// Semi-latus rectum
	if( e == 1 ) l1 = 2.0*a1; else l1 = a1*fabs(1-e*e);
	if( e == 1 ) l2 = 2.0*a2; else l2 = a2*fabs(1-e*e);
	// True anomaly
	nu1 = acos((((double)l1/(radius1)) - 1.0)/e);
	nu2 = acos((((double)l2/(radius2)) - 1.0)/e);
	// Case where the galaxies are on elliptical orbits
	if(e < 1.0) {
		if(nu1!=nu1)	{
			nu1=-pi;
			radius1 = l1/(1.0-e);
		}
		if(nu2!=nu2) {
			nu2=-pi;
			radius2 = l2/(1.0-e);
		}
	}
	
	kappa = G*(gal_1->total_mass+gal_2->total_mass);
    
	x_1 =  (radius1)*cos(nu1);
	y_1 =  (radius1)*sin(nu1);
	z_1 = 0.;
	x_2 = -(radius2)*cos(nu2);
	y_2 = -(radius2)*sin(nu2);
	z_2 = 0.;
	
	l = per*(1.0+e);
    
	vx_1 =  k1*sqrt(kappa/(l*kpc))*sin(nu1)/1.0E5;
	vy_1 = -k1*sqrt(kappa/(l*kpc))*(e+cos(nu1))/1.0E5;
	vz_1 = 0.;
	vx_2 = -k2*sqrt(kappa/(l*kpc))*sin(nu2)/1.0E5;
	vy_2 =  k2*sqrt(kappa/(l*kpc))*(e+cos(nu2))/1.0E5;
	vz_2 = 0.;
    
	// Putting galaxies on the X-axis
	if(x_1 < 0) theta_diago = asin(y_1/radius1);
	else theta_diago = asin(-y_1/radius1);
	xt_1 = x_1*cos(theta_diago)-y_1*sin(theta_diago);
	yt_1 = x_1*sin(theta_diago)+y_1*cos(theta_diago);
	zt_1 = z_1;
	xt_2 = x_2*cos(theta_diago)-y_2*sin(theta_diago);
	yt_2 = x_2*sin(theta_diago)+y_2*cos(theta_diago);
	zt_2 = z_2;
	vxt_1 = vx_1*cos(theta_diago)-vy_1*sin(theta_diago);
	vyt_1 = vx_1*sin(theta_diago)+vy_1*cos(theta_diago);
	vzt_1 = vz_1;
	vxt_2 = vx_2*cos(theta_diago)-vy_2*sin(theta_diago);
	vyt_2 = vx_2*sin(theta_diago)+vy_2*cos(theta_diago);
	vzt_2 = vz_2;
	// Inclining the orbital plane with repsect to Y-axis
	x_1 = xt_1*cos(phi*degtorad);
	y_1 = yt_1;
	z_1 = -xt_1*sin(phi*degtorad);
	x_2 = xt_2*cos(phi*degtorad);
	y_2 = yt_2;
	z_2 = -xt_2*sin(phi*degtorad);
	vx_1 = vxt_1*cos(phi*degtorad);
	vy_1 = vyt_1;
	vz_1 = -vxt_1*sin(phi*degtorad);
	vx_2 = vxt_2*cos(phi*degtorad);
	vy_2 = vyt_2;
	vz_2 = -vxt_2*sin(phi*degtorad);
	// Rotating around Z-axis
	xt_1 = x_1*cos(theta*degtorad)-y_1*sin(theta*degtorad);
	yt_1 = x_1*sin(theta*degtorad)+y_1*cos(theta*degtorad);
	zt_1 = z_1;
	xt_2 = x_2*cos(theta*degtorad)-y_2*sin(theta*degtorad);
	yt_2 = x_2*sin(theta*degtorad)+y_2*cos(theta*degtorad);
	zt_2 = z_2;
	vxt_1 = vx_1*cos(theta*degtorad)-vy_1*sin(theta*degtorad);
	vyt_1 = vx_1*sin(theta*degtorad)+vy_1*cos(theta*degtorad);
	vzt_1 = vz_1;
	vxt_2 = vx_2*cos(theta*degtorad)-vy_2*sin(theta*degtorad);
	vyt_2 = vx_2*sin(theta*degtorad)+vy_2*cos(theta*degtorad);
	vzt_2 = vz_2;
	// Put back the new coordinates in the orginal variables
	x_1 = xt_1;
	y_1 = yt_1;
	z_1 = zt_1;
	x_2 = xt_2;
	y_2 = yt_2;
	z_2 = zt_2;
	vx_1 = vxt_1;
	vy_1 = vyt_1;
	vz_1 = vzt_1;
	vx_2 = vxt_2;
	vy_2 = vyt_2;
	vz_2 = vzt_2;
	// Computing initial distance
	r_ini = sqrt(pow(x_1-x_2,2)+pow(y_1-y_2,2)+pow(z_1-z_2,2));
	// Computing initial relative velocity
	v_ini = sqrt(pow(vx_1-vx_2,2)+pow(vy_1-vy_2,2)+pow(vz_1-vz_2,2));
	// Computing specific orbital energy in 1E4 km^2.s^-2
	specific_orbital_energy = (0.5*v_ini*v_ini-((kappa/(1E15))/(sep*kpc/1E5)))/(1.0E4);
	// Setting the trajectory informations of the two galaxies
	gal_1->xc = x_1;
	gal_1->yc = y_1;
	gal_1->zc = z_1;
	gal_1->vel_xc = vx_1;
	gal_1->vel_yc = vy_1;
	gal_1->vel_zc = vz_1;
	gal_2->xc = x_2;
	gal_2->yc = y_2;
	gal_2->zc = z_2;
	gal_2->vel_xc = vx_2;
	gal_2->vel_yc = vy_2;
	gal_2->vel_zc = vz_2;
	// Printing some informations about the trajectory to the screen
	printf("/////\n////// Setting Keplerian orbit with:\n");
	printf("/////\t\t-Total mass of the galaxy 1: %8.3lf E10 solar mass\n",m1);
	printf("/////\t\t-Total mass of the galaxy 2: %8.3lf E10 solar mass\n",m2);
	printf("/////\t\t-Initial distance between galaxies: %6.1lf kpc\n",r_ini);
	printf("/////\t\t-Pericentral distance of the trajectory: %6.2lf kpc\n",per);
	printf("/////\t\t-Eccentricity of the trajectory: %6.3lf\n",e);
	printf("/////\t\t-Center of galaxy 1 x=%5.1lf kpc y=%5.1lf kpc z=%5.1lf kpc\n",x_1,y_1,z_1);
	printf("/////\t\t-Center of galaxy 2 x=%5.1lf kpc y=%5.1lf kpc z=%5.1lf kpc\n",x_2,y_2,z_2);
	printf("/////\t\t-Velocity of galaxy 1 vx=%6.1lf km.s^-1 vy=%6.1lf km.s^-1 vz=%6.1lf km.s^-1\n",vx_1,vy_1,vz_1);
	printf("/////\t\t-Velocity of galaxy 2 vx=%6.1lf km.s^-1 vy=%6.1lf km.s^-1 vz=%6.1lf km.s^-1\n",vx_2,vy_2,vz_2);
	printf("/////\t\t-Initial relative velocity: %6.2lf km.s^-1\n",v_ini);
	printf("/////\t\t-Specific Orbital energy: %6.2lf 1E4 km^2.s^-2\n",specific_orbital_energy);
	printf("/////\t\t-Orbital plane's normal vector azimuthal angle: %6.1lf degrees\n",theta);
	printf("/////\t\t-Orbital plane's normal vector polar angle: %6.1lf degrees\n",phi);
	return;
}
