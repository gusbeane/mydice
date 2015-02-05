/*-----------------------------------------------------------------------------
/
/ Filename: dice.c
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

// The All-Powerful & All-Mighty Main!
int main (int argc, char **argv) {

	int i, j, k, l, ngal, N;
	// Galaxy pointers
	galaxy *gal, *pgal, *stack1, *stack2;
	stream *st;
	//Starting clock
	clock_start = clock();
	// Set some important variables: The number of disk particles and
	// the number of halo particles, the disk mass fraction, the disk
	// angular momentum fraction, the spin parameter, the halo 
	// concentration, the virial velocity, the potential grid size,
	// and the grid spacing.
	printf("/////  _____     __     ______     ______			\n");
	printf("///// /\\  __-.  /\\ \\   /\\  ___\\   /\\  ___\\		\n");
	printf("///// \\ \\ \\/\\ \\ \\ \\ \\  \\ \\ \\____  \\ \\  __\\	\n");
	printf("/////  \\ \\____-  \\ \\_\\  \\ \\_____\\  \\ \\_____\\	\n");
	printf("/////   \\/____/   \\/_/   \\/_____/   \\/_____/		\n");

	write_dice_version();
	printf("///// Written by Valentin Perret [University of Zurich]\n");
	printf("/////\n");
	printf("/////\t--------------------------------------------------\n");
	// Looking for the galaxy parameter file
	fflush(stdout);
	if(argc != 2) {
		fprintf(stderr,"[Error] No configuration file specified\n");
		exit(0);
	}
	strcpy(AllVars.ParameterFile, argv[1]);
	if(parse_config_file(AllVars.ParameterFile) != 0) {
		fprintf(stderr,"[Error] Unable to read the DICE config file\n");
		exit(0);
	}
	// Since this example is using pointers to galaxies, allocate them.
	// Alternatively, you could declare a non-pointer galaxy and pass its
	// memory address.
	if (!(gal=calloc(1,sizeof(galaxy)))) {
	        fprintf(stderr,"[Error] Unable to allocate galaxy\n");
	        return 0;
	}
	if (!(pgal=calloc(1,sizeof(galaxy)))) {
	        fprintf(stderr,"[Error] Unable to allocate galaxy\n");
	        return 0;
	}
	if (!(stack1=calloc(1,sizeof(galaxy)))) {
	        fprintf(stderr,"[Error] Unable to allocate galaxy\n");
	        return 0;
	}
	if (!(stack2=calloc(1,sizeof(galaxy)))) {
	        fprintf(stderr,"[Error] Unable to allocate galaxy\n");
	        return 0;
	}
	if (!(st=calloc(1,sizeof(stream)))) {
	        fprintf(stderr,"[Error] Unable to allocate stream\n");
	        return 0;
	}
	// Checking threads number
	if (AllVars.Nthreads <= 0) {
		printf("/////\t[Warning] Nthreads <= 0\n");
		printf("/////\tSetting Nthreads to 1\n");
		AllVars.Nthreads = 1;
	}
	// Set the OpenMP thread number
	#if USE_THREADS == 1
		omp_set_num_threads(AllVars.Nthreads);
		printf("/////\t%d OpenMP threads\n",AllVars.Nthreads);
	#else
		AllVars.Nthreads = 1;
	#endif
	// Set the GSL random number generator
	printf("/////\tAllocating GSL random number generators\n");
	gsl_rng_env_setup();
	T = gsl_rng_mt19937;
	r = (gsl_rng **) malloc(AllVars.Nthreads * sizeof(gsl_rng *));
	for(i=0;i<AllVars.Nthreads;i++) r[i] = gsl_rng_alloc(T);
	// Set the GSL QAG integration workspace
	printf("/////\tAllocating GSL integration workspaces\n");
	w = (gsl_integration_workspace **) malloc(AllVars.Nthreads * sizeof(gsl_integration_workspace *));
	for(i=0;i<AllVars.Nthreads;i++) {
		w[i] = gsl_integration_workspace_alloc(GSL_WORKSPACE_SIZE);
	}
	if (AllVars.Ngal>0) {
		printf("/////\t--------------------------------------------------\n");
		if(AllVars.Ngal==1)printf("/////\t%d galaxy to generate\n",AllVars.Ngal);
		else printf("/////\t%d galaxies to generate\n",AllVars.Ngal);
		printf("/////\t--------------------------------------------------\n");
		printf("/////\n");
	}
	fflush(stdout);
	for(k=0; k<AllVars.Ngal; k++) {
		AllVars.CurrentGalaxy = k;
		// Provide values for the free parameters of the galaxy. If the galaxy can
		// not be created, return an error.
		printf("/////\t--------------------------------------------------\n");
		if ((i = create_galaxy(gal,AllVars.GalaxyFiles[k],1))!= 0) {
		        fprintf(stderr,"[Error] Unable to create galaxy\n");
		        exit(0);
		}
		if(k>0) {
			// If the new galaxy is the same as the previous one
			// Then there is no need to recompute everything
			if(strcmp(AllVars.GalaxyFiles[k],AllVars.GalaxyFiles[k-1]) == 0) {
				printf("/////\tSame galaxy -> Using previous computation\n");
				if(copy_galaxy(pgal,gal,0) != 0) {
					printf("[Error] Unable to copy galaxy\n");
					exit(0);
				}
				copy_potential(pgal,gal,1);
			} else {
				// Set up the particles positions, disk potential, and particle velocities of 
				// the particles in the galaxy. The option on set_galaxy_velocity tells the 
				// function to use dispersion.
				if (set_galaxy_coords(gal) != 0) {
				        fprintf(stderr,"[Error] Unable to set coordinates\n");
				        exit(0);
				}
				if(set_galaxy_potential(gal,1) != 0) {
					printf("[Error] Unable to set the potential\n");
					exit(0);
				}
				if(set_galaxy_velocity(gal) != 0) {
					printf("[Error] Unable to set the velocities\n");
					exit(0);
				}
				lower_resolution(gal);
			}
		} else {
			// If the new galaxy is different from the previous one
			// Let's recompute everything
			if (set_galaxy_coords(gal) != 0) {
			        fprintf(stderr,"[Error] Unable to set coordinates\n");
			        exit(0);
			}
			if(set_galaxy_potential(gal,1) != 0) {
				fprintf(stderr,"[Error] Unable to set the potential\n");
				exit(0);
			}
			if(set_galaxy_velocity(gal) != 0) {
				fprintf(stderr,"[Error] Unable to set the velocities\n");
				exit(0);
			}
			lower_resolution(gal);
		}
		// Apply rotation using spherical reference frame
		rotate_galaxy(gal,gal->spin,gal->incl);
		if(k == 1 && AllVars.SetKeplerian == 1) {
			// Case where the user choose to set a Keplerian trajectory
			// between two and only two galaxies
			set_orbit_keplerian(gal,stack1,AllVars.Rinit,AllVars.Rperi,AllVars.Eccentricity,AllVars.OrbitPlanePhi,AllVars.OrbitPlaneTheta);
			set_galaxy_trajectory(gal);
			set_galaxy_trajectory(stack1);
		} else {
			set_galaxy_trajectory(gal);
		}
		printf("/////\tCopying galaxy to the stack\n");
		if(k > 0) {
			// If a galaxy has been created previously
			// let's stack the result of the previous computation in &stack2
			if(create_galaxy_system(gal,stack1,stack2) != 0) {
				fprintf(stderr,"[Error] Unable to build the galaxy system\n");
				exit(0);
			}
			if(copy_galaxy(stack2,stack1,0) != 0) {
				fprintf(stderr,"[Error] Unable to copy galaxy\n");
				exit(0);
			}
			if(copy_galaxy(gal,pgal,0) != 0) {
				fprintf(stderr,"[Error] Unable to copy galaxy\n");
				exit(0);
			}
		} else {
			// If this is the first computation
			if(copy_galaxy(gal,stack1,0) != 0) {
				fprintf(stderr,"[Error] Unable to copy galaxy\n");
				exit(0);
			}
			if(copy_galaxy(gal,pgal,0) != 0) {
				fprintf(stderr,"[Error] Unable to copy galaxy\n");
				exit(0);
			}
		}
		printf("/////\t--------------------------------------------------\n");
		printf("/////\n");
	}
	
	if(AllVars.Nstream > 0) {
		printf("/////\t--------------------------------------------------\n");
		if(AllVars.Nstream==1) printf("/////\t%d stream to generate\n",AllVars.Nstream);
		else printf("/////\t%d streams to generate\n",AllVars.Nstream);
		printf("/////\t--------------------------------------------------\n");
		printf("/////\n");
	}
	for(l=0; l<AllVars.Nstream; l++) {
		// Provide values for the free parameters of the galaxy. If the galaxy can
		// not be created, return an error.
		if ((i = create_stream(st,AllVars.StreamFiles[l],1))!= 0) {
		        fprintf(stderr,"[Error] Unable to create galaxy\n");
		        exit(0);
		}
		// If the new galaxy is different from the previous one
		// Let's recompute everything
		if ((i = set_stream_coords(st)) != 0) {
		        fprintf(stderr,"[Error] Unable to set coordinates\n");
		        exit(0);
		}
		if(set_stream_velocity(st) != 0) {
			fprintf(stderr,"[Error] Unable to set the velocities\n");
			exit(0);
		}
		printf("/////\tCopying stream to the stack\n");
		if(k>0||l>0){
			// If a galaxy has been created previously
			// let's stack the result of the previous computation in &stack2
			if(add_stream_to_system(st,stack1,stack2) != 0) {
				fprintf(stderr,"[Error] Unable to add stream to the galaxy system\n");
				exit(0);
			}
			if(copy_galaxy(stack2,stack1,0) != 0) {
				fprintf(stderr,"[Error] Unable to copy galaxy\n");
				exit(0);
			}
		} else {
			if(stream_to_galaxy(st,stack1,0) != 0) {
				fprintf(stderr,"[Error] Unable to create galaxy from stream\n");
				exit(0);
			}
		}
		printf("/////\t--------------------------------------------------\n");
		printf("/////\n");
	}	
	printf("/////\t--------------------------------------------------\n");
	// Write initial conditions for the massively parallel N-body code Gadget2.
	printf("/////\tWriting %s file\n",AllVars.ICformat);
	if(k+l==1) {
		if(strcmp(AllVars.ICformat,"Gadget1")==0) write_gadget1_ics(stack1,AllVars.Filename);
		if(strcmp(AllVars.ICformat,"Gadget2")==0) write_gadget2_ics(stack1,AllVars.Filename);		
	} else {
		if(strcmp(AllVars.ICformat,"Gadget1")==0) write_gadget1_ics(stack2,AllVars.Filename);
		if(strcmp(AllVars.ICformat,"Gadget2")==0) write_gadget2_ics(stack2,AllVars.Filename);
	}
	// Destroy a galaxy. If the galaxy can not be destroyed, return an error. This
	// function will SEGFAULT if the arrays in the galaxy can not be freed.
	printf("/////\tCleaning memory\n");
	destroy_galaxy(gal,0);
	destroy_galaxy(pgal,0);
	destroy_galaxy(stack1,0);
	destroy_galaxy(stack2,0);
	// Free random number generator & GSL integration workspace
	for(i=0;i<AllVars.Nthreads;i++) {
		gsl_rng_free(r[i]);
		gsl_integration_workspace_free(w[i]);
	}
	free(r);
	free(w);
	clock_end 	= clock();
	cpu_time 	= ((double)(clock_end-clock_start))/CLOCKS_PER_SEC;
	printf("/////\tICs successfully created [%d seconds]\n",(int)ceil(cpu_time));
	return 0;
}

