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
	// Set some important variables: The number of disk particles and
	// the number of halo particles, the disk mass fraction, the disk
	// angular momentum fraction, the spin parameter, the halo 
	// concentration, the virial velocity, the potential grid size,
	// and the grid spacing.
	printf("/////\n");
	printf("/////  _____     __     ______     ______			\n");
	printf("///// /\\  __-.  /\\ \\   /\\  ___\\   /\\  ___\\		\n");
	printf("///// \\ \\ \\/\\ \\ \\ \\ \\  \\ \\ \\____  \\ \\  __\\	\n");
	printf("/////  \\ \\____-  \\ \\_\\  \\ \\_____\\  \\ \\_____\\	\n");
	printf("/////   \\/____/   \\/_/   \\/_____/   \\/_____/		\n");

	printf("/////\n///// D.I.C.E. (Disk Initial Conditions Environment)\n");
	write_dice_version();
	// Looking for the galaxy parameter file
	fflush(stdout);
	if(argc != 2) {
		fprintf(stderr,"You should pass one argument to dice!\nThis argument should be an ASCII dice configuration file.\nTake a look at the documentation for help.\n");
		exit(0);
	}
	strcpy(AllVars.ParameterFile, argv[1]);
	if(parse_config_file(AllVars.ParameterFile) != 0) {
		fprintf(stderr,"Unable to read the DICE config file. Aborting.\n");
		exit(0);
	}
	// Since this example is using pointers to galaxies, allocate them.
	// Alternatively, you could declare a non-pointer galaxy and pass its
	// memory address.
	if (!(gal=calloc(1,sizeof(galaxy)))) {
	        fprintf(stderr,"Unable to allocate galaxy. Aborting.\n");
	        return 0;
	}
	if (!(pgal=calloc(1,sizeof(galaxy)))) {
	        fprintf(stderr,"Unable to allocate galaxy. Aborting.\n");
	        return 0;
	}
	if (!(stack1=calloc(1,sizeof(galaxy)))) {
	        fprintf(stderr,"Unable to allocate galaxy. Aborting.\n");
	        return 0;
	}
	if (!(stack2=calloc(1,sizeof(galaxy)))) {
	        fprintf(stderr,"Unable to allocate galaxy. Aborting.\n");
	        return 0;
	}
	if (!(st=calloc(1,sizeof(stream)))) {
	        fprintf(stderr,"Unable to allocate stream. Aborting.\n");
	        return 0;
	}
	// Checking threads number
	if (AllVars.Nthreads <= 0) {
		printf("///// Nthreads value should be greater than 0\n");
		printf("///// Setting Nthreads to 1\n");
		AllVars.Nthreads = 1;
	}
	// Set the OpenMP thread number
	#if USE_THREADS == 1
		omp_set_num_threads(AllVars.Nthreads);
		printf("///// %d OpenMP threads\n",AllVars.Nthreads);
	#else
		AllVars.Nthreads = 1;
	#endif
	// Set the GSL QAG integration workspace
	w = (gsl_integration_workspace **) malloc(AllVars.Nthreads * sizeof(gsl_integration_workspace *));
	for(i=0;i<AllVars.Nthreads;i++) {
		w[i] = gsl_integration_workspace_alloc(GSL_WORKSPACE_SIZE);
	}
	if (AllVars.Ngal == 1) printf("///// Found %d galaxy to generate\n",AllVars.Ngal);
	else printf("///// Found %d galaxies to generate\n",AllVars.Ngal);
	fflush(stdout);
	for(k=0; k<AllVars.Ngal; k++) {
		// Provide values for the free parameters of the galaxy. If the galaxy can
		// not be created, return an error.
		if ((i = create_galaxy(gal,AllVars.GalaxyFiles[k],1))!= 0) {
		        fprintf(stderr,"Unable to create galaxy. Aborting.\n");
		        exit(0);
		}
		if(k>0) {
			// If the new galaxy is the same as the previous one
			// Then there is no need to recompute everything
			if(strcmp(AllVars.GalaxyFiles[k],AllVars.GalaxyFiles[k-1]) == 0) {
				printf("/////\tSame galaxy! Using previous computation.\n");
				if(copy_galaxy(pgal,gal,0) != 0) {
					printf("/////Unable to copy galaxy. Aborting.\n");
					exit(0);
				}
				copy_potential(pgal,gal,1);
			} else {
				// Set up the particles positions, disk potential, and particle velocities of 
				// the particles in the galaxy. The option on set_galaxy_velocity tells the 
				// function to use dispersion.
				if ((i = set_galaxy_coords(gal)) != 0) {
				        fprintf(stderr,"Unable to set coordinates. Aborting.\n");
				        exit(0);
				}
				if(set_galaxy_potential(gal,1) != 0) {
					printf("/////Unable to set the potential. Aborting.\n");
					exit(0);
				}
				if(set_galaxy_velocity(gal) != 0) {
					printf("/////Unable to set the velocities. Aborting.\n");
					exit(0);
				}
				lower_resolution(gal);
			}
		} else {
			// If the new galaxy is different from the previous one
			// Let's recompute everything
			if ((i = set_galaxy_coords(gal)) != 0) {
			        fprintf(stderr,"Unable to set coordinates. Aborting.\n");
			        exit(0);
			}
			if(set_galaxy_potential(gal,1) != 0) {
				fprintf(stderr,"Unable to set the potential. Aborting.\n");
				exit(0);
			}
			if(set_galaxy_velocity(gal) != 0) {
				fprintf(stderr,"Unable to set the velocities. Aborting.\n");
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
		if(k > 0) {
			// If a galaxy has been created previously
			// let's stack the result of the previous computation in &stack2
			printf("///// Creating galaxy system\n");
			
			if(create_galaxy_system(gal,stack1,stack2) != 0) {
				fprintf(stderr,"Unable to build the galaxy system. Aborting.\n");
				exit(0);
			}
			printf("///// Copying galaxy system to the stack\n");
			if(copy_galaxy(stack2,stack1,0) != 0) {
				fprintf(stderr,"Unable to copy galaxy. Aborting.\n");
				exit(0);
			}
			printf("///// Saving galaxy system\n");
			if(copy_galaxy(gal,pgal,0) != 0) {
				fprintf(stderr,"Unable to copy galaxy. Aborting.\n");
				exit(0);
			}
		} else {
			// If this is the first computation
			printf("///// Copying galaxy to the stack\n");
			if(copy_galaxy(gal,stack1,0) != 0) {
				fprintf(stderr,"Unable to copy galaxy. Aborting.\n");
				exit(0);
			}
			printf("///// Saving galaxy\n");
			if(copy_galaxy(gal,pgal,0) != 0) {
				fprintf(stderr,"Unable to copy galaxy. Aborting.\n");
				exit(0);
			}
		}
	}

	if(AllVars.Nstream>0&&AllVars.Ngal==0) {
		//fprintf(stderr,"Unable to set streams without any galaxy. Aborting.\n");
		//exit(0);
	}
	if (AllVars.Nstream == 1) printf("///// Found %d stream to generate\n",AllVars.Nstream);
	if (AllVars.Nstream > 1) printf("///// Found %d streams to generate\n",AllVars.Nstream);
	for(l=0; l<AllVars.Nstream; l++) {
		// Provide values for the free parameters of the galaxy. If the galaxy can
		// not be created, return an error.
		if ((i = create_stream(st,AllVars.StreamFiles[l],1))!= 0) {
		        fprintf(stderr,"Unable to create galaxy. Aborting.\n");
		        exit(0);
		}
		// If the new galaxy is different from the previous one
		// Let's recompute everything
		if ((i = set_stream_coords(st)) != 0) {
		        fprintf(stderr,"Unable to set coordinates. Aborting.\n");
		        exit(0);
		}
		if(set_stream_velocity(st) != 0) {
			fprintf(stderr,"Unable to set the velocities. Aborting.\n");
			exit(0);
		}
		if(k>0||l>0){
			// If a galaxy has been created previously
			// let's stack the result of the previous computation in &stack2
			printf("///// Adding stream to system\n");
			if(add_stream_to_system(st,stack1,stack2) != 0) {
				fprintf(stderr,"Unable to add stream to the galaxy system. Aborting.\n");
				exit(0);
			}
			printf("///// Copying galaxy system to the stack\n");
			if(copy_galaxy(stack2,stack1,0) != 0) {
				fprintf(stderr,"Unable to copy galaxy. Aborting.\n");
				exit(0);
			}
		} else {
			printf("///// Creating system from streams\n");
			if(stream_to_galaxy(st,stack1,0) != 0) {
				fprintf(stderr,"Unable to create galaxy from stream. Aborting.\n");
				exit(0);
			}
		}
	}	
	
	// Write initial conditions for the massively parallel N-body code Gadget2.
	// If you do not use Gadget, do not use this function -- the output is binary!
	if(k+l==1) {
		write_gadget_ics(stack1,AllVars.Filename);
	} else {
		write_gadget_ics(stack2,AllVars.Filename);
	}
	// Destroy a galaxy. If the galaxy can not be destroyed, return an error. This
	// function will SEGFAULT if the arrays in the galaxy can not be freed.
	printf("///// Cleaning memory\n");
	destroy_galaxy(gal,0);
	destroy_galaxy(pgal,0);
	destroy_galaxy(stack1,0);
	destroy_galaxy(stack2,0);
	// Destroy the random number environment.
	for(i=0;i<AllVars.Nthreads;i++) {
		gsl_rng_free(r[i]);
	}
	for(i=0;i<AllVars.Nthreads;i++) {
		gsl_integration_workspace_free(w[i]);
	}
	free(r);
	printf("///// ICs successfully created\n");
	return 0;
}

