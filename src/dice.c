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
/ Date: September 2015
/
*///---------------------------------------------------------------------------


// Include the DICE header
#include "dice.h"

// The All-Powerful & All-Mighty Main!
int main (int argc, char **argv) {

	int i, j, k, l, ngal, N, neval;
	double x,y;
	// Galaxy pointers
	galaxy *gal, *pgal, *stack;
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
	printf("/////\tReading configuration file [%s]\n",AllVars.ParameterFile);
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
	if (!(stack=calloc(1,sizeof(galaxy)))) {
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
		w[i] = gsl_integration_workspace_alloc(AllVars.GslWorkspaceSize);
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
		printf("/////\t--------------------------------------------------\n");
		printf("/////\tReading galaxy %d params file [%s]\n",k+1,AllVars.GalaxyFiles[k]);
		// First iteration
		if(k==0) {
			if((i = create_galaxy(gal,AllVars.GalaxyFiles[k],1)) != 0) {
			 	fprintf(stderr,"[Error] Unable to create galaxy\n");
				exit(0);
			}
			if(copy_galaxy(gal,stack,0) != 0) {
				fprintf(stderr,"[Error] Unable to copy galaxy\n");
				exit(0);
			}
		// Next iterations
		} else {
			// Same galaxy parameter file
			if(strcmp(AllVars.GalaxyFiles[k],AllVars.GalaxyFiles[k-1]) == 0) {
				printf("/////\tSame galaxy -> Using previous computation\n");
				if(copy_galaxy(pgal,gal,0) != 0) {
					printf("[Error] Unable to copy galaxy\n");
					exit(0);
				}
			// New galaxy parameter file
			} else {
				if((i = create_galaxy(gal,AllVars.GalaxyFiles[k],1)) != 0) {
			        fprintf(stderr,"[Error] Unable to create galaxy\n");
			        exit(0);
				}	
			}
			printf("/////\tCopying galaxy to the stack\n");
			if(add_galaxy_to_system(gal,stack) != 0) {
				fprintf(stderr,"[Error] Unable to build the galaxy system\n");
				exit(0);
			}

		}
		// Store galaxy for next iteration
		if(k>0) trash_galaxy(pgal,0);
		if(copy_galaxy(gal,pgal,0) != 0) {
			fprintf(stderr,"[Error] Unable to copy galaxy\n");
			exit(0);
		}
		// Store galaxy properties
		AllVars.GalMass[k] 		= gal->total_mass;
		if(k==0) {
			AllVars.GalStart[k]	= 0;
		} else {
			AllVars.GalStart[k]	= AllVars.GalStart[k-1]+AllVars.GalNpart[k-1];
		}
		AllVars.GalNpart[k]		= gal->ntot_part;
		trash_galaxy(gal,0);
		printf("/////\t--------------------------------------------------\n");
		printf("/////\n");
	}
	
	// Computing Keplerian trajectories
	for(k=0; k<AllVars.Ngal; k++) {
		if(AllVars.Kepler_Gal1[k]>0 && AllVars.Kepler_Gal2[k]>0) {
			if(AllVars.GalMass[AllVars.Kepler_Gal1[k]-1]>0 && AllVars.GalMass[AllVars.Kepler_Gal2[k]-1]>0) {
				set_orbit_keplerian(AllVars.Kepler_Gal1[k]-1,AllVars.Kepler_Gal2[k]-1,
					AllVars.Kepler_Rinit[k],AllVars.Kepler_Rperi[k],AllVars.Kepler_Ecc[k],
					AllVars.Kepler_OrbitPlanePhi[k],AllVars.Kepler_OrbitPlaneTheta[k],AllVars.Kepler_GalCenter[k]);
			}
		}
	}
	// Position galaxies
	printf("/////\t--------------------------------------------------\n");
	printf("/////\tGalaxies trajectories\n");
	for(k=0; k<AllVars.Ngal; k++) {
    	printf("/////\t\tGalaxy %d -> [x=%5.1lf y=%5.1lf z=%5.1lf][kpc] [vx=%5.1lf vy=%5.1lf vz=%5.1lf][km/s] [spin=%5.1lf incl=%5.1lf][deg]\n",k,
    		AllVars.GalPos[k][0],AllVars.GalPos[k][1],AllVars.GalPos[k][2],
    		AllVars.GalVel[k][0],AllVars.GalVel[k][1],AllVars.GalVel[k][2],
    		AllVars.GalSpin[k],AllVars.GalIncl[k]);
		// Apply rotation using spherical reference frame
		rotate_galaxy(stack,k);
		set_galaxy_trajectory(stack,k);
	}
	printf("/////\t--------------------------------------------------\n");

	if(AllVars.Nstream > 0) {
		printf("/////\t--------------------------------------------------\n");
		if(AllVars.Nstream==1) printf("/////\t%d stream to generate\n",AllVars.Nstream);
		else printf("/////\t%d streams to generate\n",AllVars.Nstream);
		printf("/////\t--------------------------------------------------\n");
		printf("/////\n");
	}
	for(l=0; l<AllVars.Nstream; l++) {
		if(l==0) {
			AllVars.StreamStart[l] = stack->ntot_part;
		} else {
			AllVars.StreamStart[l] = AllVars.StreamStart[l-1]+AllVars.StreamNpart[l-1];
		}
		// Provide values for the free parameters of the galaxy. If the galaxy can
		// not be created, return an error.
		printf("/////\t--------------------------------------------------\n");
		if ((i = create_stream(st,AllVars.StreamFiles[l],1))!= 0) {
		        fprintf(stderr,"[Error] Unable to create galaxy\n");
		        exit(0);
		}
		printf("/////\tCopying stream to the stack\n");
		if(k>0||l>0){
			// If a galaxy has been created previously
			// let's stack the result of the previous computation in &stack2
			if(add_stream_to_system(st,stack) != 0) {
				fprintf(stderr,"[Error] Unable to add stream to the galaxy system\n");
				exit(0);
			}
		} else {
			if(stream_to_galaxy(st,stack,0) != 0) {
				fprintf(stderr,"[Error] Unable to create galaxy from stream\n");
				exit(0);
			}
		}

		AllVars.StreamNpart[l] = st->ntot_part;
		printf("/////\t--------------------------------------------------\n");
		printf("/////\n");
	}
	// Position streams
	if(AllVars.Nstream > 0) {
		printf("/////\t--------------------------------------------------\n");
		printf("/////\tStreams positions\n");
		for(k=0; k<AllVars.Nstream; k++) {
			if(AllVars.StreamNpart[k]>0) {
				printf("/////\t\tStream %d -> [x=%5.1lf y=%5.1lf z=%5.1lf][kpc] [spin=%5.1lf incl=%5.1lf][deg]\n",k,
    				AllVars.StreamPos[k][0],AllVars.StreamPos[k][1],AllVars.StreamPos[k][2],
    				AllVars.StreamSpin[k],AllVars.StreamIncl[k]);
				position_stream(stack,k);
				rotate_stream(stack,k);
			}
		}
	}
	printf("/////\t--------------------------------------------------\n");
	// Write initial conditions for the massively parallel N-body code Gadget2.
	printf("/////\tWriting %s file\n",AllVars.ICformat);
	if(strcmp(AllVars.ICformat,"Gadget1")==0) write_gadget1_ics(stack,AllVars.Filename);
	if(strcmp(AllVars.ICformat,"Gadget2")==0) write_gadget2_ics(stack,AllVars.Filename);		
	printf("/////\tCleaning memory\n");
	trash_galaxy(pgal,0);
	trash_galaxy(stack,0);
	free(gal);
	free(pgal);
	free(stack);
	free(st);
	// Free random number generator & GSL integration workspace
	for(i=0;i<AllVars.Nthreads;i++) {
		gsl_rng_free(r[i]);
		gsl_integration_workspace_free(w[i]);
	}
	free(r);
	free(w);
	// Stop the clock and print execution time
	clock_end 	= clock();
	cpu_time 	= ((double)(clock_end-clock_start))/CLOCKS_PER_SEC;
	printf("/////\tICs successfully created [%d seconds]\n",(int)ceil(cpu_time));
	printf("/////\t--------------------------------------------------\n");
	return 0;
}

