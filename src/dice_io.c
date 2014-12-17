/*-----------------------------------------------------------------------------
 /
 / Filename: dice_io.c
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

// This function prints the release date of the current version of DICE.
void write_dice_version() {
	printf("///// Version %d.%d \n",DICE_VERSION_MAJOR,DICE_VERSION_MINOR);
	return;
}

// This function parses the dice configuration file
// which specifies the location of the galaxy parameter files
// and the initial trajectories of the galaxies
int parse_config_file(char *fname) {
    
#define DOUBLE	1
#define STRING	2
#define INT	3
#define MAXTAGS	300
    
	FILE *fd;
	int i,j;
	char buf[200], buf1[200], buf2[200], junk[200];
	int nt;
	int id[MAXTAGS];
	int read[MAXTAGS];
	int mandatory[MAXTAGS];
	void *addr[MAXTAGS];
	char tag[MAXTAGS][50];
	int  errorFlag = 0;
	
	// By default, we will use the computation of Keplerian orbits
	// If one parameter is missing to compute this trajectory,
	// or if there is too much or not enough galaxies
	// then the galaxies will be initialized unsing [xc,yc,zc,vel_xc,vel_yc,vel_zc] keywords
	//AllVars.SetKeplerian = 1;
	
	if(sizeof(long long) != 8) {
		fprintf(stderr,"\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
		return -1;
	}
	
	if(sizeof(int) != 4) {
		fprintf(stderr,"\nType `int' is not 32 bit on this platform. Stopping.\n\n");
		return -1;
	}
	
	if(sizeof(float) != 4) {
		fprintf(stderr,"\nType `float' is not 32 bit on this platform. Stopping.\n\n");
		return -1;
	}
	
	if(sizeof(double) != 8) {
		fprintf(stderr,"\nType `double' is not 64 bit on this platform. Stopping.\n\n");
		return -1;
	}
    
	
	if((fd = fopen(fname,"r"))) {
		j = 0;
		while(!feof(fd)) {
			*buf = 0;
			fgets(buf, 200, fd);
			if(sscanf(buf, "%s%s", buf1, buf2) < 2) continue;
			if(buf1[0] == '%') continue;
			if(strcmp(buf1,"Galaxy") == 0) {
				strcpy(AllVars.GalaxyFiles[j],buf2);
				j++;
                
			}
            
		}
		fclose(fd);
	}

	AllVars.Ngal = j;
	
	if((fd = fopen(fname,"r"))) {
		j = 0;
		while(!feof(fd)) {
			*buf = 0;
			fgets(buf, 200, fd);
			if(sscanf(buf, "%s%s", buf1, buf2) < 2) continue;
			if(buf1[0] == '%') continue;
			if(strcmp(buf1,"Stream") == 0) {
				strcpy(AllVars.StreamFiles[j],buf2);
				j++;
                
			}
            
		}
		fclose(fd);
	}
	AllVars.Nstream = j;
	
	if(AllVars.Ngal+AllVars.Nstream == 0) {
		fprintf(stderr,"////// Error - No galaxy/stream parameters files specified.\n");
		return -1;
	}

	nt = 0;
	
	strcpy(tag[nt], "SetKeplerian");
	addr[nt] = &AllVars.SetKeplerian;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = INT;

	strcpy(tag[nt], "Eccentricity");
	addr[nt] = &AllVars.Eccentricity;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "Rinit");
	addr[nt] = &AllVars.Rinit;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "Rperi");
	addr[nt] = &AllVars.Rperi;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = DOUBLE;
    
	strcpy(tag[nt], "OrbitPlanePhi");
	addr[nt] = &AllVars.OrbitPlanePhi;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = DOUBLE;
    
	strcpy(tag[nt], "OrbitPlaneTheta");
	addr[nt] = &AllVars.OrbitPlaneTheta;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "Nthreads");
	AllVars.Nthreads = 4;
	addr[nt] = &AllVars.Nthreads;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = INT;
	
	strcpy(tag[nt], "GasHydrostaticEq");
	AllVars.GasHydrostaticEq = 0;
	addr[nt] = &AllVars.GasHydrostaticEq;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = INT;
	
	strcpy(tag[nt], "GasHydrostaticEqIter");
	AllVars.GasHydrostaticEqIter = 10;
	addr[nt] = &AllVars.GasHydrostaticEqIter;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = INT;
	
	strcpy(tag[nt], "MeanPartDist");
	addr[nt] = &AllVars.MeanPartDist;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = INT;
    
	strcpy(tag[nt], "AcceptImaginary");
	AllVars.AcceptImaginary = 0;
	addr[nt] = &AllVars.AcceptImaginary;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = INT;
    
	strcpy(tag[nt], "OutputRc");
	AllVars.OutputRc = 0;
	addr[nt] = &AllVars.OutputRc;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = INT;
	
	strcpy(tag[nt], "MaxCompNumber");
	AllVars.MaxCompNumber = 10;
	addr[nt] = &AllVars.MaxCompNumber;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = INT;
    
	strcpy(tag[nt], "Filename");
	addr[nt] = &AllVars.Filename;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = STRING;
	
	// Default values for cosmological parameters is Planck cosmology
	strcpy(tag[nt], "H0");
	addr[nt] = &AllVars.H0;
	AllVars.H0 = 67.77;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "OmegaM");
	addr[nt] = &AllVars.Omega_m;
	AllVars.Omega_m = 0.30712;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "OmegaL");
	addr[nt] = &AllVars.Omega_l;
	AllVars.Omega_l = 0.691391;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "OmegaK");
	addr[nt] = &AllVars.Omega_k;
	AllVars.Omega_k = 0.00;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = DOUBLE;
    
	printf("///// Reading DICE config file: %s\n",fname);
	if((fd = fopen(fname, "r"))) {
		while(!feof(fd)) {
			*buf = 0;
			fgets(buf, 200, fd);
			if(sscanf(buf, "%s%s", buf1, buf2) < 2) continue;
			if(buf1[0] == '%' || buf1[0] == '#' || strcmp(buf1,"Galaxy") == 0 || strcmp(buf1,"Stream") == 0) continue;
			for(i = 0, j = -1; i < nt; i++)
				if(strcmp(buf1, tag[i]) == 0) {
					j = i;
					read[i] = 1;
					break;
				}
			
			if(j >= 0) {
				switch (id[j]) {
					case DOUBLE:
						*((double *) addr[j]) = atof(buf2);
						break;
					case STRING:
						strcpy(addr[j], buf2);
						break;
					case INT:
						*((int *) addr[j]) = atoi(buf2);
						break;
				}
				
			} else {
				fprintf(stdout, "////// Error in file %s - Tag '%s' not allowed or multiple defined.\n",
                        fname, buf1);
				return -1;
			}
		}
		fclose(fd);
	} else {
		fprintf(stderr,"\tDICE config file %s not found.\n", fname);
		return -2;
	}
	
	for(i = 0; i < nt; i++) {
		if(read[i] == 0 && mandatory[i] == 1) {
			fprintf(stderr,"////// Error - I miss a value for '%s' \n",tag[i]);
			AllVars.SetKeplerian = 0;
			return -1;
		}
	}
    
#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
	
	return 0;
}

// This function parses a galaxy parameter file, i.e. a file containing all the needed physical informations
// in order to build the galaxy
int parse_galaxy_file(galaxy *gal, char *fname) {
	#define DOUBLE	1
	#define STRING	2
	#define INT	3
	#define MAXTAGS	26*AllVars.MaxCompNumber+19
	
	FILE *fd;
	int i,j,n;
	char buf[400], buf1[400], buf2[400];
	int nt;
	int id[MAXTAGS];
	void *addr[MAXTAGS];
	char tag[MAXTAGS][400];
	char temp_tag[400];
	int read[MAXTAGS];
	int mandatory[MAXTAGS];
	int errorFlag = 0;
	
	if(sizeof(long long) != 8) {
		fprintf(stderr,"\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
		return -1;
	}
	
	if(sizeof(int) != 4) {
		fprintf(stderr,"\nType `int' is not 32 bit on this platform. Stopping.\n\n");
		return -1;
	}
	
	if(sizeof(float) != 4) {
		fprintf(stderr,"\nType `float' is not 32 bit on this platform. Stopping.\n\n");
		return -1;
	}
	
	if(sizeof(double) != 8) {
		fprintf(stderr,"\nType `double' is not 64 bit on this platform. Stopping.\n\n");
		return -1;
	}
	
	nt = 0;
	
	strcpy(tag[nt], "v200");
	addr[nt] = &gal->v200;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "m200");
	gal->m200 = 0.;
	addr[nt] = &gal->m200;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "redshift");
	addr[nt] = &gal->redshift;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "j_d");
	addr[nt] = &gal->j_d;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "lambda");
	addr[nt] = &gal->lambda;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "level_grid");
	addr[nt] = &gal->level_grid;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = INT;
	
	strcpy(tag[nt], "level_grid_dens");
	addr[nt] = &gal->level_grid_dens;
	gal->level_grid_dens=7;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = INT;
	
	strcpy(tag[nt], "level_grid_turb");
	addr[nt] = &gal->level_grid_turb;
	gal->level_grid_turb=7;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = INT;
	
	strcpy(tag[nt], "boxsize");
	addr[nt] = &gal->boxsize;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "axisymmetric_drift");
	addr[nt] = &gal->axisymmetric_drift;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = INT;
	
	strcpy(tag[nt], "xc");
	addr[nt] = &gal->xc;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "yc");
	addr[nt] = &gal->yc;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "zc");
	addr[nt] = &gal->zc;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "vel_xc");
	addr[nt] = &gal->vel_xc;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "vel_yc");
	addr[nt] = &gal->vel_yc;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "vel_zc");
	addr[nt] = &gal->vel_zc;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "spin");
	addr[nt] = &gal->spin;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "incl");
	addr[nt] = &gal->incl;
	read[nt] = 0;
	mandatory[nt] = 1;
	id[nt++] = DOUBLE;
	
	strcpy(tag[nt], "seed");
	gal->seed = time(NULL);
	addr[nt] = &gal->seed;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = INT;
	
	for(j=0; j<AllVars.MaxCompNumber; j++) {
	
		n = sprintf(temp_tag,"mass_frac%d",j+1);
		strcpy(tag[nt], temp_tag);
		addr[nt] = &gal->comp_mass_frac[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = DOUBLE;
	
		n = sprintf(temp_tag,"model%d",j+1);
		strcpy(tag[nt], temp_tag);
		addr[nt] = &gal->comp_model[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = INT;

		n = sprintf(temp_tag,"scale_length%d",j+1);
		strcpy(tag[nt], temp_tag);	
		addr[nt] = &gal->comp_scale_length[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = DOUBLE;
	
		n = sprintf(temp_tag,"cut%d",j+1);
		strcpy(tag[nt], temp_tag);	
		addr[nt] = &gal->comp_cut[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = DOUBLE;
	
		n = sprintf(temp_tag,"flat%d",j+1);
		strcpy(tag[nt], temp_tag);	
		addr[nt] = &gal->comp_flat[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = DOUBLE;
	
		n = sprintf(temp_tag,"mcmc_step%d",j+1);
		strcpy(tag[nt], temp_tag);	
		addr[nt] = &gal->comp_mcmc_step[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"mcmc_step_hydro%d",j+1);
		gal->comp_mcmc_step_hydro[j] = 1.0;
		strcpy(tag[nt], temp_tag);
		addr[nt] = &gal->comp_mcmc_step_hydro[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;

		n = sprintf(temp_tag,"vmax%d",j+1);
		strcpy(tag[nt], temp_tag);	
		addr[nt] = &gal->comp_vmax[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = DOUBLE;
	
		n = sprintf(temp_tag,"npart%d",j+1);
		strcpy(tag[nt], temp_tag);
		addr[nt] = &gal->comp_npart[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = INT;
		
		n = sprintf(temp_tag,"npart_pot%d",j+1);
		strcpy(tag[nt], temp_tag);
		addr[nt] = &gal->comp_npart_pot[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = INT;
	
		n = sprintf(temp_tag,"type%d",j+1);
		strcpy(tag[nt], temp_tag);
		gal->comp_type[j] = -1;	
		addr[nt] = &gal->comp_type[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = INT;
		
		n = sprintf(temp_tag,"concentration%d",j+1);
		strcpy(tag[nt], temp_tag);		
		gal->comp_concentration[j] = 0.;
		addr[nt] = &gal->comp_concentration[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
						
		n = sprintf(temp_tag,"streaming_fraction%d",j+1);
		strcpy(tag[nt], temp_tag);		
		gal->comp_streaming_fraction[j] = 0.;
		addr[nt] = &gal->comp_streaming_fraction[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"theta_sph%d",j+1);
		strcpy(tag[nt], temp_tag);		
		gal->comp_theta_sph[j] = 0.;
		addr[nt] = &gal->comp_theta_sph[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"phi_sph%d",j+1);
		strcpy(tag[nt], temp_tag);		
		gal->comp_phi_sph[j] = 0.;
		addr[nt] = &gal->comp_phi_sph[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"metal%d",j+1);
		strcpy(tag[nt], temp_tag);		
		gal->comp_metal[j] = 0.;
		addr[nt] = &gal->comp_metal[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"t_init%d",j+1);
		strcpy(tag[nt], temp_tag);		
		gal->comp_t_init[j] = 1e4;
		addr[nt] = &gal->comp_t_init[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"mean_age%d",j+1);
		strcpy(tag[nt], temp_tag);		
		gal->comp_mean_age[j] = 0.;
		addr[nt] = &gal->comp_mean_age[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"min_age%d",j+1);
		strcpy(tag[nt], temp_tag);		
		gal->comp_min_age[j] = 0.;
		addr[nt] = &gal->comp_min_age[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"alpha%d",j+1);
		strcpy(tag[nt], temp_tag);		
		gal->comp_alpha[j] = 1.0;
		addr[nt] = &gal->comp_alpha[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;

		n = sprintf(temp_tag,"disp_ext%d",j+1);
		strcpy(tag[nt], temp_tag);
		gal->comp_disp_ext[j] = 1.0;
		addr[nt] = &gal->comp_disp_ext[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"radius_nfw%d",j+1);
		strcpy(tag[nt], temp_tag);
		gal->comp_radius_nfw[j] = -1.0;
		addr[nt] = &gal->comp_radius_nfw[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"Q_lim%d",j+1);
		strcpy(tag[nt], temp_tag);
		gal->comp_Q_lim[j] = 0.;
		addr[nt] = &gal->comp_Q_lim[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"turb_sigma%d",j+1);
		strcpy(tag[nt], temp_tag);
		gal->comp_turb_sigma[j] = 0.;
		addr[nt] = &gal->comp_turb_sigma[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"turb_scale%d",j+1);
		strcpy(tag[nt], temp_tag);
		gal->comp_turb_scale[j] = 0.;
		addr[nt] = &gal->comp_turb_scale[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"compute_vel%d",j+1);
		strcpy(tag[nt], temp_tag);
		gal->comp_compute_vel[j] = 1;
		addr[nt] = &gal->comp_compute_vel[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = INT;
	}
	
	printf("/////\n///// Reading galaxy params file: %s\n",fname);
	if((fd = fopen(fname, "r"))) {
		while(!feof(fd)) {
			*buf = 0;
			fgets(buf, 400, fd);
			if(sscanf(buf, "%s%s", buf1, buf2) < 2) continue;
			if(buf1[0] == '%' || buf1[0] == '#' ) continue;
			for(i = 0, j = -1; i < nt; i++)
				if(strcmp(buf1,tag[i]) == 0) {
					j = i;
					read[i] = 1;
					break;
				}
			if(j >= 0) {
				switch (id[j]) {
					case DOUBLE:
						*((double *) addr[j]) = atof(buf2);
						break;
					case STRING:
						strcpy(addr[j], buf2);
						break;
					case INT:
						*((int *) addr[j]) = atoi(buf2);
						break;
				}
			} else {
				fprintf(stderr,"////// Error in file %s - Tag '%s' not allowed or multiple defined.\n",fname, buf1);
				return -1;
			}
		}
		fclose(fd);
	} else {
		fprintf(stderr,"\nParameter file %s not found.\n\n", fname);
		return -2;
	}
	
	for(i = 0; i < nt; i++) {
		if(read[i] == 0 && mandatory[i] == 1) {
			fprintf(stderr,"////// Error - I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
			return -1;
		}
	}
    
	#undef DOUBLE
	#undef STRING
	#undef INT
	#undef MAXTAGS

	return 0;
}


// This function parses a galaxy parameter file, i.e. a file containing all the needed physical informations
// in order to build the galaxy
int parse_stream_file(stream *st, char *fname) {
	#define DOUBLE	1
	#define STRING	2
	#define INT	3
	#define MAXTAGS	500
	
	FILE 	*fd;
	int 	i,j,n;
	char 	buf[400], buf1[400], buf2[400];
	int 	nt;
	int 	id[MAXTAGS];
	void 	*addr[MAXTAGS];
	char 	tag[MAXTAGS][200];
	char	temp_tag[200];
	int 	read[MAXTAGS];
	int 	mandatory[MAXTAGS];
	int 	errorFlag = 0;
	
	if(sizeof(long long) != 8) {
		fprintf(stderr,"\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
		return -1;
	}
	
	if(sizeof(int) != 4) {
		fprintf(stderr,"\nType `int' is not 32 bit on this platform. Stopping.\n\n");
		return -1;
	}
	
	if(sizeof(float) != 4) {
		fprintf(stderr,"\nType `float' is not 32 bit on this platform. Stopping.\n\n");
		return -1;
	}
	
	if(sizeof(double) != 8) {
		fprintf(stderr,"\nType `double' is not 64 bit on this platform. Stopping.\n\n");
		return -1;
	}
	
	nt = 0;
		
	for(j=0; j<AllVars.MaxCompNumber; j++) {
	
		n = sprintf(temp_tag,"dens%d",j+1);
		strcpy(tag[nt], temp_tag);
		addr[nt] = &st->comp_dens[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = DOUBLE;
	
		n = sprintf(temp_tag,"model%d",j+1);
		strcpy(tag[nt], temp_tag);
		addr[nt] = &st->comp_model[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = INT;

		n = sprintf(temp_tag,"length%d",j+1);
		strcpy(tag[nt], temp_tag);	
		addr[nt] = &st->comp_length[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"scale%d",j+1);
		strcpy(tag[nt], temp_tag);	
		addr[nt] = &st->comp_scale[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = DOUBLE;
	
		n = sprintf(temp_tag,"opening_angle%d",j+1);
		strcpy(tag[nt], temp_tag);	
		addr[nt] = &st->comp_opening_angle[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = DOUBLE;
	
		n = sprintf(temp_tag,"mcmc_step%d",j+1);
		strcpy(tag[nt], temp_tag);	
		addr[nt] = &st->comp_mcmc_step[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = DOUBLE;

		n = sprintf(temp_tag,"turb_sigma%d",j+1);
		strcpy(tag[nt], temp_tag);	
		addr[nt] = &st->comp_turb_sigma[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"turb_scale%d",j+1);
		strcpy(tag[nt], temp_tag);	
		addr[nt] = &st->comp_turb_scale[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = DOUBLE;
	
		n = sprintf(temp_tag,"npart%d",j+1);
		strcpy(tag[nt], temp_tag);
		addr[nt] = &st->comp_npart[j];
		read[nt] = 0;
		if(j==0) mandatory[nt] = 1;
		else mandatory[nt] = 0;
		id[nt++] = INT;
		
		n = sprintf(temp_tag,"theta_sph%d",j+1);
		strcpy(tag[nt], temp_tag);		
		st->comp_theta_sph[j] = 0.;
		addr[nt] = &st->comp_theta_sph[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"phi_sph%d",j+1);
		strcpy(tag[nt], temp_tag);		
		st->comp_phi_sph[j] = 0.;
		addr[nt] = &st->comp_phi_sph[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"metal%d",j+1);
		strcpy(tag[nt], temp_tag);		
		st->comp_metal[j] = 0.;
		addr[nt] = &st->comp_metal[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;		
		
		n = sprintf(temp_tag,"t_init%d",j+1);
		strcpy(tag[nt], temp_tag);		
		st->comp_t_init[j] = 1e4;
		addr[nt] = &st->comp_t_init[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;

		n = sprintf(temp_tag,"xc%d",j+1);
		strcpy(tag[nt], temp_tag);		
		st->comp_xc[j] = 0.;
		addr[nt] = &st->comp_xc[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"yc%d",j+1);
		strcpy(tag[nt], temp_tag);		
		st->comp_yc[j] = 0.;
		addr[nt] = &st->comp_yc[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
		
		n = sprintf(temp_tag,"zc%d",j+1);
		strcpy(tag[nt], temp_tag);		
		st->comp_zc[j] = 0.;
		addr[nt] = &st->comp_zc[j];
		read[nt] = 0;
		mandatory[nt] = 0;
		id[nt++] = DOUBLE;
	}

	strcpy(tag[nt], "seed");
	st->seed = time(NULL);
	addr[nt] = &st->seed;
	read[nt] = 0;
	mandatory[nt] = 0;
	id[nt++] = INT;
	
	printf("/////\n///// Reading galaxy params file: %s\n",fname);
	if((fd = fopen(fname, "r"))) {
		while(!feof(fd)) {
			*buf = 0;
			fgets(buf, 400, fd);
			if(sscanf(buf, "%s%s", buf1, buf2) < 2) continue;
			if(buf1[0] == '%' || buf1[0] == '#' ) continue;
			for(i = 0, j = -1; i < nt; i++)
				if(strcmp(buf1,tag[i]) == 0) {
					j = i;
					read[i] = 1;
					break;
				}
			if(j >= 0) {
				switch (id[j]) {
					case DOUBLE:
						*((double *) addr[j]) = atof(buf2);
						break;
					case STRING:
						strcpy(addr[j], buf2);
						break;
					case INT:
						*((int *) addr[j]) = atoi(buf2);
						break;
				}
			} else {
				fprintf(stderr,"////// Error in file %s - Tag '%s' not allowed or multiple defined.\n",fname, buf1);
				return -1;
			}
		}
		fclose(fd);
	} else {
		fprintf(stderr,"\nParameter file %s not found.\n\n", fname);
		return -2;
	}
	
	for(i = 0; i < nt; i++) {
		if(read[i] == 0 && mandatory[i] == 1) {
			fprintf(stderr,"////// Error - I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
			return -1;
		}
	}
    
	#undef DOUBLE
	#undef STRING
	#undef INT
	#undef MAXTAGS

	return 0;
}


// The next three functions write the position coordinates of a galaxy to the
// screen in xyz format.
void write_galaxy_position(galaxy *gal) {
	
	int i,j;
	char filename[200];
	FILE *fp1;
	
	// Total rotation curve
	sprintf(filename, "positions.dat");
	fp1 = fopen(filename, "w");
	
	
    fprintf(stderr,"Writing galaxy positions to standard out...\n");
	
	
	//Print particles positions
	for(j=0; j<AllVars.MaxCompNumber; j++) {
		if (gal->comp_npart[j] > 0) {
			for (i = gal->comp_npart[j]; i < gal->comp_npart[j] + gal->comp_npart[j]; ++i) {
				fprintf(fp1,"%d %le %le %le \n",gal->comp_type[j],gal->x[i],gal->y[i],gal->z[i]);
			}
		}
	}
	
	fclose(fp1);
	
	return;
}

// The next three functions write the velocity coordinates of a galaxy to
// the screen in xyz format.
void write_galaxy_velocity(galaxy *gal, int info) {
	
	int i;
	
	//Print disk positions
	if (info != 0) {
		fprintf(stderr,"Writing galaxy velocities to standard out...\n");
		printf("/////#Disk: \n");
	}
	if (gal->num_part[1] > 0) {
		for (i = gal->num_part[0]; i < gal->num_part[0] + gal->num_part[1]; ++i) {
			printf("/////%le %le %le \n",gal->vel_x[i],gal->vel_y[i],gal->vel_z[i]);
		}
	}
	
	//Print halo positions
	if (info != 0) {
		printf("/////#Halo: \n");
	}
	if (gal->num_part[0] != 0) {
		for (i = 0; i < gal->num_part[0] ; ++i) {
			printf("/////%lf %lf %lf \n",gal->vel_x[i],gal->vel_y[i],gal->vel_z[i]);
		}
	}
	
	return;
}

void write_galaxy_velocity_disk(galaxy *gal) {
	
	int i;
	
	fprintf(stderr,"Writing disk velocities to standard out...\n");
	printf("/////#Disk: \n");
	if (gal->num_part[1] > 0) {
		for (i = gal->num_part[0]; i < gal->num_part[0] + gal->num_part[1]; ++i) {
			printf("/////%le %le %le \n",gal->vel_x[i],gal->vel_y[i],gal->vel_z[i]);
		}
	}
	
	return;
}

void write_galaxy_velocity_halo(galaxy *gal) {
	
	int i;
	
	fprintf(stderr,"Writing halo velocities to standard out...\n");
	printf("/////#Halo: \n");
	if (gal->num_part[0] != 0) {
		for (i = 0; i < gal->num_part[0] ; ++i) {
			printf("/////%lf %lf %lf \n",gal->vel_x[i],gal->vel_y[i],gal->vel_z[i]);
		}
	}
	
	return;
}




// Function to write initial conditions to file in default Gadget2 format.
// The code was originally an input routine, read_snapshot.c, provided by
// Volker Springel with the Gadget2 source code. It has been hacked into a
// write routine.
int write_gadget_ics(galaxy *gal, char *fname) {
    
	FILE *fp1, *fp2;
	int dummy, ntot_withmasses, NumPart, ptype;
	unsigned long int i, j, k;
	int t,n,off,pc,pc_new,pc_sph;
	int files = 1;
	int *Ids;
	char buf[200];
#define SKIP2 fwrite(&dummy, sizeof(dummy), 1, fp1);
    
	
	// Set everything to zero and overwrite it later if needed.
	header1.npart[0] = 0;
	header1.npart[1] = 0;
	header1.npart[2] = 0;
	header1.npart[3] = 0;
	header1.npart[4] = 0;
	header1.npart[5] = 0;
	header1.npartTotal[0] = 0;
	header1.npartTotal[1] = 0;
	header1.npartTotal[2] = 0;
	header1.npartTotal[3] = 0;
	header1.npartTotal[4] = 0;
	header1.npartTotal[5] = 0;
	header1.mass[0] = 0.0;
	header1.mass[1] = 0.0;
	header1.mass[2] = 0.0;
	header1.mass[3] = 0.0;
	header1.mass[4] = 0.0;
	header1.mass[5] = 0.0;
    
	// Set the header values to some defaults.
	header1.npart[0] = gal->num_part[0];
	header1.npart[1] = gal->num_part[1];
	header1.npart[2] = gal->num_part[2];
	header1.npart[3] = gal->num_part[3];
	header1.npartTotal[0] = gal->num_part[0];
	header1.npartTotal[1] = gal->num_part[1];
	header1.npartTotal[2] = gal->num_part[2];
	header1.npartTotal[3] = gal->num_part[3];
	header1.time = 0.0;
	header1.redshift = 0.0;
	header1.flag_sfr = 0.0;
	header1.flag_feedback = 0.0;
	header1.flag_cooling = 0.0;
	header1.num_files = 1;
	header1.BoxSize = 0.0;
	header1.Omega0 = 0.0;
	header1.OmegaLambda = 0.0;
	header1.HubbleParam = 1.0;
  	
	if (!(Ids=malloc(gal->ntot_part*sizeof(int)))) {
		fprintf(stderr,"Unable to create particle Ids structure in memory.");
		exit(0);
	}
    
	if (!(P=malloc(gal->ntot_part*sizeof(struct particle_data)))) {
		fprintf(stderr,"Unable to create particle data structure in memory.");
		exit(0);
	}
	P--;
    
	// Transfer the particle data from the DICE galaxy data structure to the Gadget
	// position_data structure.
	j = 1;
	//We need to transfer in the order of particle type as defined in GADGET2. 0->Gas 1->Disk 2->Halo etc.
	for (ptype=0; ptype<10; ptype++) {
		for (k=0;k<AllVars.MaxCompNumber;k++) {
			if(gal->comp_type[k]==ptype&&gal->comp_npart[k]>0) {
				for (i = gal->comp_start_part[k]; i < gal->comp_start_part[k] + gal->comp_npart[k]; ++i) {
					P[j].Pos[0] = gal->x[i];
					P[j].Pos[1] = gal->y[i];
					P[j].Pos[2] = gal->z[i];
					P[j].Vel[0] = gal->vel_x[i];
					P[j].Vel[1] = gal->vel_y[i];
					P[j].Vel[2] = gal->vel_z[i];
					P[j].U 		= gal->u[i];
					P[j].Rho 	= gal->rho[i];
					P[j].Mass 	= gal->mass[i];
					P[j].Type 	= gal->comp_type[k];
					P[j].Metal  = gal->metal[i];
					P[j].Age	= gal->age[i];
					++j;
				}
			}
		}
	}
    
	for(i=0, pc=1; i<files; i++, pc=pc_new){
		if(files>1) sprintf(buf,"%s.g1.%d",fname,(int)i);
		else sprintf(buf,"%s.g1",fname);
        
		if(!(fp1=fopen(buf,"w"))){
			fprintf(stderr,"can't open file `%s`\n",buf);
			exit(0);
		}
        
		fflush(stdout);
        // Header
		dummy = sizeof(header1);
		fwrite(&dummy, sizeof(dummy), 1, fp1);
		fwrite(&header1, sizeof(header1), 1, fp1);
		fwrite(&dummy, sizeof(dummy), 1, fp1);
        
		for(k=0, ntot_withmasses=0; k<6; k++){
			if(header1.mass[k]==0) ntot_withmasses+= header1.npart[k];
		}
        
		dummy = 3*sizeof(float)*gal->ntot_part;
		SKIP2;
		// Positions
		for(k=0,pc_new=pc;k<6;k++) {
			for(n=0;n<header1.npart[k];n++){
				fwrite(&P[pc_new].Pos[0], sizeof(float), 3, fp1);
				pc_new++;
			}
		}
		SKIP2;
        // Velocities
		SKIP2;
		for(k=0,pc_new=pc;k<6;k++) {
			for(n=0;n<header1.npart[k];n++){
				fwrite(&P[pc_new].Vel[0], sizeof(float), 3, fp1);
				pc_new++;
			}
		}
		SKIP2;
        // Identifiers
		dummy = sizeof(int)*gal->ntot_part;
		SKIP2;
		for(k=0,pc_new=pc;k<6;k++) {
			for(n=0;n<header1.npart[k];n++){
				Ids[pc_new] = pc_new;
				fwrite(&Ids[pc_new], sizeof(int), 1, fp1);
				pc_new++;
			}
		}
		SKIP2;
		// Mass        
		if(ntot_withmasses>0) {
			dummy = sizeof(float)*gal->ntot_part;
			SKIP2;
		}
		for(k=0, pc_new=pc; k<6; k++) {
			for(n=0;n<header1.npart[k];n++) {
				P[pc_new].Type=k;
				if(header1.mass[k]==0) fwrite(&P[pc_new].Mass, sizeof(float), 1, fp1);
				else P[pc_new].Mass = header1.mass[k];
				pc_new++;
			}
		}
		if(ntot_withmasses>0) {
			SKIP2;
		}
        // Gas specific datablocks
		if(header1.npart[0]>0) {
			// Internal energy
			dummy = sizeof(float)*header1.npart[0];
			SKIP2;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++) {
				fwrite(&P[pc_sph].U, sizeof(float), 1, fp1);
				pc_sph++;
			}
			SKIP2;
            // Density
			SKIP2;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++) {
				fwrite(&P[pc_sph].Rho, sizeof(float), 1, fp1);
				pc_sph++;
			}
			SKIP2;
		}
		// Metallicity
		dummy = sizeof(int)*gal->ntot_part;
		SKIP2;
		for(k=0,pc_new=pc;k<6;k++) {
			for(n=0;n<header1.npart[k];n++){
				fwrite(&P[pc_new].Metal, sizeof(float), 1, fp1);
				pc_new++;
			}
		}
		SKIP2;
		// Age for star particles
		if(gal->ntot_part_stars>0){
			dummy = sizeof(int)*(gal->ntot_part_stars);
			SKIP2;
			for(k=2,pc_new=header1.npart[0]+header1.npart[1]+1;k<6;k++) {
				for(n=0;n<header1.npart[k];n++){
					fwrite(&P[pc_new].Age, sizeof(float), 1, fp1);
					pc_new++;
				}
			}
			SKIP2;
		}
	}
    P++;
	free(P);
	fclose(fp1);
	return 0;
}

// This function calculates and writes the rotation curve for a particular galaxy.
// The circular velocity is calculated out to the virial radius.
//
// Unlike some of the other output functions, this function writes the output to
// a set of files -- rcurve.dat, rcurve_disk.dat, rcurve_halo.dat -- These files
// can be plotted in gnuplot with the command
//
//    plot 'rcurve.dat','rcurve_disk.dat','rcurve_halo.dat'
//
// to produce a cumulative rotation curve and its parts on the same graph.
void write_galaxy_rotation_curve(galaxy *gal, double rmax) {
	int i,tid;
	double radius, theta, v_c, save;
	char filename[200];
	FILE *fp1;
    
	// Total rotation curve
	sprintf(filename, "rcurve.dat");
	fp1 = fopen(filename, "w");
	if (fp1 == NULL) {
		printf("/////Could not open rcurve.dat for writing the rotation curve. Aborting.\n");
	}
	#if USE_THREADS == 1
	    tid = omp_get_thread_num();
	#else
	    tid = 0;
	#endif
	gal->index[tid] 		= 0;
	save 					= gal->z[gal->index[tid]];
	gal->z[gal->index[tid]] = 0.;
	for (i = 1; i < 500; ++i) {
		radius = (i/500.)*rmax;
		v_c = v_c_func(gal,radius)/1.0E5;
		// Write the radius and circular velocity to file in kpc and km/s respectively.
		fprintf(fp1,"%lf %lf\n",radius,v_c);
	}
	gal->z[gal->index[tid]] = save;
	fclose(fp1);
	
	return;
}

// This function writes the potential for a particular galaxy. The function
// calculates the potentials of the disk and halo and their sum, (the total
// gravitational potential), and outputs in an nxy format to a file. This
// data can be plotted with xmgrace using a command like
//
//    xmgrace -nxy potential.dat
//
// or with gnuplot using splot and a script like
//
//    splot 'potential.dat' using 1:2:3, \
//         'potential.dat' using 1:2:4, \
//         'potential.dat' using 1:2:5
//
// to produce a detailed picture of the potential structure of the system.
// The user may want to consider adding "w lines" to that script since
// there are a lot of points in the data file.
//
// Please note that the potential of the disk is calculated only in the x,y
// plane. The potentials are calculated to the +/- the virial radius.
void write_galaxy_potential(galaxy *gal, double xmax, double ymax, double z) {
	
	int i,j,k;
	char filename[200];
	FILE *fp1;
	
	if (gal->ntot_part > 0) {
		sprintf(filename, "potential.dat");
		fp1 = fopen(filename, "w");
		if (fp1 == NULL) {
			fprintf(stderr,"Could not open potential.dat for writing the");
			fprintf(stderr," potential. Skipping it...\n");
			return;
		}
		for (i = -(int) xmax; i < (int) xmax; i=i+1) {
			for (j = -(int) ymax; j < (int) ymax; j=j+1) {
				fprintf(fp1,"%d %d %le\n",i,j,galaxy_potential_func(gal,(double) i, (double) j,z));
			}
		}
		fclose(fp1);
	} else {
		fprintf(stderr,"There are no particles in the galaxy!\n");
	}
	return;
}

void write_galaxy_potential_grid(galaxy *gal, int zindex) {
	
	int i,j,k;
	char filename[200];
	FILE *fp1;
	
	if (gal->ntot_part > 0) {
		sprintf(filename, "potential_grid.dat");
		fp1 = fopen(filename, "w");
		if (fp1 == NULL) {
			fprintf(stderr,"Could not open potential.dat for writing the");
			fprintf(stderr," potential. Skipping it...\n");
			return;
		}
		for (i = 0; i < gal->ngrid_padded; i=i+1) {
			for (j = 0; j < gal->ngrid_padded; j=j+1) {
				fprintf(fp1,"%lf %lf %le\n",(double)(i-gal->ngrid_padded/2-1)*gal->space[0],(double)(j-gal->ngrid_padded/2-1)*gal->space[1],gal->potential[i][j][zindex]);
			}
		}
		fclose(fp1);
	} else {
		fprintf(stderr,"There are no particles in the galaxy!\n");
	}
	return;
}

// This function copies one galaxy to another.
void copy_potential(galaxy *gal_1, galaxy *gal_2, int info) {
	
	unsigned long int i, j, k;
	
	if(info == 1) printf("/////\tCopying potential grid \n");
	// Copy all the coordinate information.
	for (i = 0; i < gal_1->ngrid*2; ++i) {
		for (j = 0; j < gal_1->ngrid*2; ++j) {
			for(k = 0; k < gal_1->ngrid*2; ++k){
				gal_2->potential[i][j][k] = gal_1->potential[i][j][k];
			}
		}
	}
	gal_2->potential_defined = 1;
	if(info == 1) printf("/////\tPotential grid copied \n");
	return;
}

/* this routine loads particle data from Gadget's default
 * binary file format (A snapshot may be distributed
 * into multiple files).
 *
 * Author: Volker Springel, Max Planck Institute
 */
int load_snapshot(char *fname, int files) {
	FILE *fd;
	char buf[200];
	int i,j,k,dummy,ntot_withmasses,NumPart,Ngas;
	int t,n,off,pc,pc_new,pc_sph;
	int *Id;
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
	
	for(i=0, pc=1; i<files; i++, pc=pc_new) {
		if(files>1) sprintf(buf,"%s.%d",fname,i);
		else sprintf(buf,"%s",fname);
		
		if(!(fd=fopen(buf,"r")))
		{
			printf("/////can't open file `%s`\n",buf);
			exit(0);
		}
		
		printf("/////reading `%s' ...\n",buf); fflush(stdout);
		
		fread(&dummy, sizeof(dummy), 1, fd);
		fread(&header1, sizeof(header1), 1, fd);
		fread(&dummy, sizeof(dummy), 1, fd);
		
		if(files==1)
		{
			for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
				NumPart+= header1.npart[k];
			Ngas= header1.npart[0];
		}
		else
		{
			for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
				NumPart+= header1.npartTotal[k];
			Ngas= header1.npartTotal[0];
		}
		
		for(k=0, ntot_withmasses=0; k<5; k++)
		{
			if(header1.mass[k]==0)
				ntot_withmasses+= header1.npart[k];
		}
		
		if(i==0) {
			printf("/////allocating memory...\n");
            
			if(!(P=malloc(NumPart*sizeof(struct particle_data)))) {
				fprintf(stderr,"failed to allocate memory.\n");
				exit(0);
			}
			P--;   // start with offset 1
			if(!(Id=malloc(NumPart*sizeof(int)))) {
				fprintf(stderr,"failed to allocate memory.\n");
				exit(0);
			}
			Id--;   // start with offset 1
            
			printf("/////allocating memory...done\n");
		}
		SKIP;
		for(k=0,pc_new=pc;k<6;k++)
		{
			for(n=0;n<header1.npart[k];n++)
			{
				fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
				pc_new++;
			}
		}
		SKIP;
		
		SKIP;
		for(k=0,pc_new=pc;k<6;k++)
		{
			for(n=0;n<header1.npart[k];n++)
			{
				fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
				pc_new++;
			}
		}
		SKIP;
		
		
		SKIP;
		for(k=0,pc_new=pc;k<6;k++)
		{
			for(n=0;n<header1.npart[k];n++)
			{
				fread(&Id[pc_new], sizeof(int), 1, fd);
				pc_new++;
			}
		}
		SKIP;
		
		
		if(ntot_withmasses>0)
			SKIP;
		for(k=0, pc_new=pc; k<6; k++)
		{
			for(n=0;n<header1.npart[k];n++)
			{
				P[pc_new].Type=k;
				
				if(header1.mass[k]==0)
					fread(&P[pc_new].Mass, sizeof(float), 1, fd);
				else
					P[pc_new].Mass= header1.mass[k];
				pc_new++;
			}
		}
		if(ntot_withmasses>0)
			SKIP;
		
		
		if(header1.npart[0]>0)
		{
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
			{
				fread(&P[pc_sph].U, sizeof(float), 1, fd);
				pc_sph++;
			}
			SKIP;
			
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
			{
				fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
				pc_sph++;
			}
			SKIP;
			
			if(header1.flag_cooling)
			{
				SKIP;
				for(n=0, pc_sph=pc; n<header1.npart[0];n++)
				{
					fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
					pc_sph++;
				}
				SKIP;
			}
			else
				for(n=0, pc_sph=pc; n<header1.npart[0];n++)
				{
					P[pc_sph].Ne= 1.0;
					pc_sph++;
				}
		}
		
		fclose(fd);
	}
	return 0;
}

/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 *
 * Author: Volker Springel, Max Planck Institute
 */
int reordering(int NumPart, int *Id) {
	int i,j;
	int idsource, idsave, dest;
	struct particle_data psave, psource;
    
	printf("/////reordering....\n");
    
	for(i=1; i<=NumPart; i++) {
		if(Id[i] != i) {
			psource= P[i];
			idsource=Id[i];
			dest=Id[i];
			do {
				psave= P[dest];
				idsave=Id[dest];
				P[dest]= psource;
				Id[dest]= idsource;
				if(dest == i) break;
				psource= psave;
				idsource=idsave;
				dest=idsource;
			}
			while(1);
		}
	}
	printf("/////done.\n");
	Id++;
	free(Id);
	printf("/////space for particle ID freed\n");
	return 0;
}

// This function frees the memory used when the load_snapshot function is
// called. In general, users will not know that P has to be incremented
// because it was decremented and an attempt to free it will result in
// a segmentation fault. This prevents the problem.
int unload_snapshot() {
    
    P++; free(P);
    
    return 0;
}

// Load a Gadget2 snapshot into a galaxy object.
int load_gadget2_galaxy(char *filename, galaxy *galaxy_1) {
	
	int i,j;
	unsigned long int parts[4];
	double cuts[6];
	
	load_snapshot(filename,1);
	parts[0] = header1.npart[2];
	parts[1] = header1.npart[0];
	parts[2] = header1.npart[1];
	parts[3] = header1.npart[3];
	
	
	// Create the galaxy, but pass false information for the irrelevant parts
	// at the moment.
	create_galaxy(galaxy_1,AllVars.GalaxyFiles[0],0);
	
	i = 1;
	for (j = galaxy_1->num_part[0]+galaxy_1->num_part[1]; j < galaxy_1->num_part[3]; ++j) {
		galaxy_1->x[j] = P[i].Pos[0];
		galaxy_1->y[j] = P[i].Pos[1];
		galaxy_1->z[j] = P[i].Pos[2];
		galaxy_1->vel_x[j] = P[i].Vel[0];
		galaxy_1->vel_y[j] = P[i].Vel[1];
		galaxy_1->vel_z[j] = P[i].Vel[2];
		galaxy_1->mass[j] = P[i].Mass*unit_mass;
		++i;
	}
	for (j = galaxy_1->num_part[0]; j < galaxy_1->num_part[0]+galaxy_1->num_part[1]; ++j) {
		galaxy_1->x[j] = P[i].Pos[0];
		galaxy_1->y[j] = P[i].Pos[1];
		galaxy_1->z[j] = P[i].Pos[2];
		galaxy_1->vel_x[j] = P[i].Vel[0];
		galaxy_1->vel_y[j] = P[i].Vel[1];
		galaxy_1->vel_z[j] = P[i].Vel[2];
		galaxy_1->mass[j] = P[i].Mass*unit_mass;
		++i;
	}
	for (j = 0; j < galaxy_1->num_part[0]; ++j) {
		galaxy_1->x[j] = P[i].Pos[0];
		galaxy_1->y[j] = P[i].Pos[1];
		galaxy_1->z[j] = P[i].Pos[2];
		galaxy_1->vel_x[j] = P[i].Vel[0];
		galaxy_1->vel_y[j] = P[i].Vel[1];
		galaxy_1->vel_z[j] = P[i].Vel[2];
		galaxy_1->mass[j] = P[i].Mass*unit_mass;
		++i;
	}
	// It is important to always unload the snapshot as soon as you use it.
	unload_snapshot();
	
	return 0;
}