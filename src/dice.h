/*-----------------------------------------------------------------------------
 /
 / Filename: dice.h
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

// Header files to include
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <err.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_ellint.h>
#include <complex.h>
#include <fftw3.h>
#include "dice_config.h"

#if USE_THREADS == 1
	#include <omp.h>
#endif

#define MAXLEN_FILENAME  	100
#define MAXLEN_FILELINE  	500
#define GSL_WORKSPACE_SIZE	100000



// A tolerance parameter for the GSL QAG integration.
#define epsabs 1.0E-10
#define epsrel 1.0E-7
#define key 6

// Some physical constants needed for the computation
#define gamma					(5.0/3.0)			// adiabatic index of simulated gas (mono-atomic)
#define gamma_minus1			(gamma-1.0)			// adiabatic index of simulated gas minus 1
#define hydrogen_massfrac		0.76				// mass fraction of hydrogen, relevant only for radiative cooling
#define boltzmann				1.3806200e-16		// Boltzmann constant in [erg.K^-1] = [g.cm^2.s^-2.K^-1]
#define protonmass				1.6600000e-24		// proton mass in [g]
#define solarmass				1.989E33			// Mass of the Sun in [g]
#define mu_mol					1.2195e0			// Molecular weight
#define pi						3.14159274101257 	// pi = 4.0*atan(1.0)
#define	G						6.67428E-8			// G = 6.67428E-8 [cm^3 g^-1 s^-2]
#define	kpc						3.085678E21			// 1 kpc = 3.085678E21 [cm]
#define unit_mass				1.989E43			// 1.989E43 = 1E10 [Solar masses]
#define unit_velocity			1e5					// km.s^-1 in [cm.s^-1]
#define unit_length				3.085678e21			// kpc in [cm]
#define unit_time				(unit_length/unit_velocity)
#define unit_energy				(unit_mass*(unit_length*unit_length)/(unit_time*unit_time))
#define unit_dens				(unit_mass/(unit_length*unit_length*unit_length))
#define unit_nh					(hydrogen_massfrac/protonmass*unit_dens)


// Global variables for the GSL random number environment.
const gsl_rng_type *T;
gsl_rng **r;
gsl_integration_workspace **w;

// Clock timer
clock_t clock_start,clock_end;
double cpu_time;

// This is a type definition of a galaxy. The thought here is to create galaxies
// as objects and, hopefully, make the code extremely clean. It is essentially a
// a collection of arrays and constants.
typedef struct {
	// Redshift of the galaxy
	double redshift;
	// Spin parameter of the galaxy
	double lambda;
	// Spin fraction of the disk
	double j_d;
	// Baryonic mass fraction of the disk
	double m_d;
	// Total mass of the galactic model
	double total_mass;
	// Density fluctuation dispersion
	double dens_gauss_sigma;
	// Density fluctuation scale
	double dens_gauss_scale;
	// Total number of components
	int n_component;
	int *selected_comp;
	// Component quantities
	char				**comp_profile_name;
	unsigned long int 	*comp_npart;
	unsigned long int 	*comp_npart_pot;
	unsigned long int 	*comp_start_part;
	double 				*comp_mass_frac;
	double 				*comp_mass;
	int 				*comp_model;
	double 				*comp_cutted_mass;
	double 				*comp_scale_length;
	double 				*comp_concentration;
	double 				*comp_scale_height;
	double 				*comp_cut;
	double 				*comp_flat;
	double 				*comp_mcmc_step;
	double 				*comp_mcmc_step_hydro;
	double 				*comp_vmax;
	int 				*comp_type;
	int 				*comp_bool;
	double 				*comp_streaming_fraction;
	double 				*comp_cut_dens;
	double 				*comp_theta_sph;
	double 				*comp_phi_sph;
	double 				*comp_metal;
	double 				*comp_t_init;
	double 				*comp_u_init;
	double 				*comp_cs_init;
	double 				*comp_turb_sigma;
	double 				*comp_turb_scale;
	double 				*comp_mean_age;
	double 				*comp_min_age;
	double 				*comp_alpha;
	double 				*comp_disp_ext;
	double				*comp_scale_dens;
	double				*comp_radius_nfw;
	double 				*comp_Q_lim;
	double				*comp_Q_min;
	int					*comp_compute_vel;
	int 				*comp_hydro_eq;
	int 				*comp_hydro_eq_niter;
	int 				*comp_dens_gauss;
	double 				*comp_cut_in;
	int 				*comp_thermal_eq;
	double 				*comp_part_mass;
	// Virial quantities
	double v200;
	double r200;
	double m200;
	// Coordinates vectors
	double *x;
	double *y;
	double *z;
	double *r_cyl;
	double *theta_cyl;
	double *r_sph;
	double *theta_sph;
	double *phi_sph;
	// Galaxy's barycenter coordinates
	double xc;
	double yc;
	double zc;
	// Velocities vectors
	double *vel_x;
	double *vel_y;
	double *vel_z;
	// Galaxy's barycenter velocity
	double vel_xc;
	double vel_yc;
	double vel_zc;
	// Galaxy's spin orientation spherical angles
	double spin;
	double incl;
	// Mass vector
	double *mass;
	// Density vector
	double *rho;
	// Internal energy vector
	double *u;
	// Age vector
	double *age;
	// Metallicity vector
	double *metal;
	// Particle Mesh potential grid
	double ***potential;
	// Particle Mesh potential grid
	double ***potential_zoom;
	// Particle Mesh gaussian field grid
	double ***gaussian_field;
	// Gas midplane density grid
	double **midplane_dens;
	// Potential grid cell size vector [kpc]
	double dx;
	// Potential zoom grid cell size vector [kpc]
	double dx_zoom;
	// Gas midplane density grid cell size vector [kpc]
	double dx_dens;
	// Gaussian field grid cell size vector [kpc]
	double dx_gauss;
	// Number of particles per particle type
	unsigned long int num_part[4];
	unsigned long int num_part_pot[4];
	// Total potential grid size [kpc]
	double boxsize;
	// Total potential zoom grid size [kpc]
	double boxsize_zoom;
	// Total gas midplane density grid size [kpc]
	double boxsize_dens;
	// Level of refinement of the potential grid
	int level_grid;
	// Level of refinement of the potential zoom grid
	int level_grid_zoom;
	// Level of refinement of the gas density grid
	int level_grid_dens;
	// Level of refinement of the gas turbulence grid
	int level_grid_turb;
	// Level of refinement of the density gaussian fluctuations
	int level_grid_dens_gauss;
	// Number of cells in the potential grid
	int ngrid[3];
	// Number of cells in the zoomed potential grid
	int ngrid_zoom[3];
	// Number of cells in the midplane density grid
	int ngrid_dens[2];
	// Number of cells in the turbulence grid
	int ngrid_gauss[3];
	// Storage array
	double **storage;
	// Identifier of particle
	unsigned long int *index;
	unsigned long int *id;
	// Total number of particles
	unsigned long int ntot_part;
	unsigned long int ntot_part_stars;
	unsigned long int ntot_part_pot;
	// Boolean variable checking the computation of the potential
	int potential_defined;
	// Boolean variable checking the computation of the gaussian field
	int gaussian_field_defined;
	// Boolean which enable the axisymmetric drift approximation for the stellar disk
	int epicycle;
	// Seed for random number generator
	long seed;
	// Pseudo density boolean
	int *pseudo;
	// Shift term for the gravitational potential in the zoom region
	double potential_shift_zoom;
} galaxy;

// This is a type definition of a stream. The thought here is to create streams
// as objects and, hopefully, make the code extremely clean. It is essentially a
// a collection of arrays and constants.
typedef struct {
	int *selected_comp;
	// Stream's peak coordinates
	double 				*comp_xc;
	double 				*comp_yc;
	double 				*comp_zc;
	double 				*comp_mass;
	double 				*comp_dens;
	double 				*comp_opening_angle;
	double 				*comp_turb_sigma;
	double 				*comp_turb_scale;
	double 				*comp_t_init;
	// Stream's spin orientation spherical angles
	double 				*comp_theta_sph;
	double 				*comp_phi_sph;
	double 				*comp_length;
	double 				*comp_scale;
	unsigned long int 	*comp_npart;
	int					*comp_bool;
	int					*comp_model;
	double 				*comp_mcmc_step;
	double 				*comp_metal;
	double 				*comp_u_init;
	double 				*comp_cs_init;
	unsigned long int 	*comp_start_part;
	char 				**comp_profile_name;
	int					*comp_dens_gauss;
	// Coordinates vectors
	double *x;
	double *y;
	double *z;
	double *r_cyl;
	double *theta_cyl;
	double *r_sph;
	double *theta_sph;
	double *phi_sph;
	// Mass vector
	double *mass;
	// Density vector
	double *rho;
	// Internal energy vector
	double *u;
	// Metallicity vector
	double *metal;
	// Velocities vectors
	double *vel_x;
	double *vel_y;
	double *vel_z;
	unsigned long int *id;
	// Level of refinement of the gas turbulence grid
	int level_grid_turb;
	// Level of refinement of the density gaussian fluctuations
	int level_grid_dens_gauss;
	// Number of cells in the turbulence grid
	unsigned long int ngrid_gauss[3];
	// Gaussian field grid cell size vector [kpc]
	double dx_gauss;
	// Particle Mesh gaussian_field grid
	double ***gaussian_field;
	// Storage array
	double **storage;
	unsigned long int ntot_part;
	double total_mass;
	int n_component;
	long seed;
	// Density fluctuation dispersion
	double dens_gauss_sigma;
	// Density fluctuation scale
	double dens_gauss_scale;
	// Boolean variable checking the computation of the gaussian field
	int gaussian_field_defined;
} stream;

//Gadget2-style header for Gadget2 snapshots.
struct io_header_1 {
	//char fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8];
	int npart[6];                        /*!< number of particles of each type in this file */
	double mass[6];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
                                          stored in the mass-block of the snapshot file, otherwise they are omitted */
	double time;                         /*!< time of snapshot file */
	double redshift;                     /*!< redshift of snapshot file */
	int flag_sfr;                        /*!< flags whether the simulation was including star formation */
	int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
	unsigned int npartTotal[6];          /*!< total number of particles of each type in this snapshot. This can be
                                          different from npart if one is dealing with a multi-file snapshot. */
	int flag_cooling;                    /*!< flags whether cooling was included  */
	int num_files;                       /*!< number of files in multi-file snapshot */
	double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
	double Omega0;                       /*!< matter density in units of critical density */
	double OmegaLambda;                  /*!< cosmological constant parameter */
	double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
	int flag_stellarage;                 /*!< flags whether the file contains formation times of star particles */
	int flag_metals;                     /*!< flags whether the file contains metallicity values for gas and star particles */
	unsigned int npartTotalHighWord[6];  /*!< High word of the total number of particles of each type */
	int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */   
	char fill[60]; 
} header1;

//Gadget2-style particle data structure.
struct particle_data {
	float Pos[3];
	float Vel[3];
	float Mass;
	int Type;
	float Pot;
	float U;
	float Rho;
	float Ne;
	float Metal;
	float Age;
	float Hsml;
} *P;

// Global variables such as the parameters of the Keplerian system, or the number of galaxies to generate.
struct GlobalVars {
	// Variable containing the parameter file's name
	char ParameterFile[MAXLEN_FILENAME];
	char GalaxyFiles[64][MAXLEN_FILENAME];
	char StreamFiles[64][MAXLEN_FILENAME];
	char Filename[MAXLEN_FILENAME];
	// Variable contained in the DICE parameter file
	char ICformat[MAXLEN_FILENAME];
	double Eccentricity;
	double Rinit;
	double Rperi;
	int SetKeplerian;
	double OrbitPlaneTheta;
	double OrbitPlanePhi;
	int Ngal;
	int Nstream;
	// Number of threads to launch when using fftw3_threads library
	int Nthreads;
	int MeanPartDist;
	int AcceptImaginary;
	int OutputRc;
	int OutputGasRc;
	int OutputPot;
	int OutputRho;
	int MaxCompNumber;
	double H0;
	double Omega_m;
	double Omega_l;
	double Omega_k;
	int GaussianRejectIter;
	int CurrentGalaxy;
} AllVars;

// ------------------------------------------
// These are the function prototypes for DICE.
// ------------------------------------------
// Initialization and destruction functions
int allocate_component_arrays(galaxy *);
int allocate_variable_arrays(galaxy *);
int allocate_component_arrays_stream(stream *);
int allocate_variable_arrays_stream(stream *);
int create_galaxy(galaxy *, char *, int);
void allocate_galaxy_storage_variable(galaxy *, int);
void allocate_stream_storage_variable(stream *, int);
void allocate_dispersion_grid();
void destroy_dispersion_grid();
int set_galaxy_coords(galaxy *);
int set_galaxy_velocity(galaxy *);
void destroy_galaxy(galaxy *, int);
void destroy_galaxy_system(galaxy *, int);
int create_galaxy_system(galaxy *, galaxy *, galaxy *);
int rotate_galaxy(galaxy *, double, double);
int rotate_component(galaxy *, double, double, int);
int set_galaxy_trajectory(galaxy *);
int copy_galaxy(galaxy *, galaxy *, int);
void set_orbit_keplerian(galaxy *, galaxy*, double, double, double, double, double);
int add_stream_to_system(stream *, galaxy *, galaxy *);
int set_stream_coords(stream *);
int set_stream_velocity(stream *);
int rotate_stream(stream *, double, double, int);
int position_stream(stream *, double, double, double, int);

// Structure functions
double density_functions_pool(galaxy *, double, double, double, int, int, int);
void mcmc_metropolis_hasting(galaxy *, int, int);
void mcmc_metropolis_hasting_stream(stream *, int , int); 
int set_hydro_equilibrium(galaxy *, int, int);
double density_functions_stream_pool(stream *, double, double, double, int, int);

double surface_density_func(galaxy *, double, double, int, int);
static double integrand_density_func(double, void *);
double cumulative_mass_func(galaxy *, double, int);
static double d_cumulative_mass_func1(double, void *);
static double d_cumulative_mass_func2(double, void *);

double surface_density_func_stream(stream *, double, double, int);
static double integrand_density_func_stream(double, void *);
double cumulative_mass_func_stream(stream *, double, int);
static double d_cumulative_mass_func_stream(double, void *);

static double integrand_density_gas_func(double, void *);
double pseudo_density_gas_func(galaxy *, double, double, double, int, int, int);
double midplane_density_gas_func(galaxy *, gsl_integration_workspace *, double, double, int);
static double dmidplane_density_gas_func(double, void *);
void fill_midplane_dens_grid(galaxy *, int);
double get_midplane_density(galaxy *, double, double);

double disk_scale_length_func(galaxy *, double);
double f_c_func(double);
double g_c_func(double);
static double dg_c_func(double, void *);
double f_s_func(double, double);
double mean_interparticle_distance(galaxy *, int);
void lower_resolution(galaxy *);

// Velocity functions
double v_c_func(galaxy *,double);
double v_c_exp_disk_func(galaxy *,double);
double v2a_z_func(galaxy *, gsl_integration_workspace *, int);
static double dv2a_z_func(double, void *);
double v2a_1D_func(galaxy *, gsl_integration_workspace *, int);
static double dv2a_1D_func(double, void *);
double v2a_theta_func(galaxy *, double, int);
double sigma2_theta_func(galaxy *, double, double);
double v2a_z_toomre(galaxy *, double, double, int);
double rho_v2a_r_func(double, void *);
double v2_theta_gas_func(galaxy *, double, double, int);
double gas_density_wrapper_func(double, void *);
int set_turbulent_grid(galaxy *, int);
double galaxy_turbulence_func(galaxy *, double, double, double, int);

// Potential and force functions
int set_galaxy_potential(galaxy *, double ***, double, int [3], int);
double galaxy_potential_func(galaxy *, double ***, double, int [3], double, double, double, int);
double galaxy_zforce_func(galaxy *, double);
double galaxy_rforce_func(galaxy *, double);
double galaxy_rsphforce_func(galaxy *, double);
double galaxyr_potential_wrapper_func(double, void *);
double galaxyrsph_potential_wrapper_func(double, void *);
double galaxyz_potential_wrapper_func(double, void *);
double potential_deriv_wrapper_func(double, void *);
void copy_potential(galaxy *, galaxy *, int);

// Input, output, and manipulation functions
void write_dice_version();
int parse_config_file(char *);
int parse_galaxy_file(galaxy *, char *);
void write_galaxy_rotation_curve(galaxy *, double, char *, double);
void write_galaxy_potential_curve(galaxy *, double, char *, double);
int write_gadget1_ics(galaxy *, char *);
int write_gadget2_ics(galaxy *, char *);
int load_snapshot(char *, int);
int allocate_memory(int);
int reordering(int, int *);
int unload_snapshot();

// Toolbox functions and definitions
double min(double, double);
double max(double, double);
typedef double (*function_to_derivate)(double, void *);
double deriv_central(galaxy *, double, double, function_to_derivate);
double deriv_forward(galaxy *, double, double, function_to_derivate);
int set_galaxy_gaussian_field_grid(galaxy *, double);
int set_stream_gaussian_field_grid(stream *, double);
double galaxy_gaussian_field_func(galaxy *, double, double, double);
double stream_gaussian_field_func(stream *, double, double, double);
