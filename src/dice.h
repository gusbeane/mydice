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
   / Date: September 2015
   /
 */

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
#include <gsl/gsl_monte_vegas.h>
#include <complex.h>
#include <fftw3.h>
#include "dice_config.h"

#if USE_THREADS == 1
#include <omp.h>
#endif

#define MAXLEN_FILENAME     100
#define MAXLEN_FILELINE     500
#define MAX_GAL             64
#define MAX_STREAM          64
#define GSL_WORKSPACE_SIZE  100000
#define NSTORAGE            8


// A tolerance parameter for the GSL QAG integration.
#define epsabs  1.0E-10
#define epsrel  1.0E-7
#define key     1

// Some physical constants needed for the computation
#define gamma                   (5.0/3.0)           // adiabatic index of simulated gas (mono-atomic)
#define gamma_minus1            (gamma-1.0)
#define hydrogen_massfrac       0.76                // mass fraction of hydrogen
#define boltzmann               1.3806200e-16       // Boltzmann constant in [erg.K^-1] = [g.cm^2.s^-2.K^-1]
#define boltzmann_kev           8.6173324e-8        // Boltzmann constant in [kev.K^-1]
#define protonmass              1.6600000e-24       // proton mass in [g]
#define solarmass               1.989E33            // Mass of the Sun in [g]
#define mu_mol                  1.0/(hydrogen_massfrac/1.0+(1.0-hydrogen_massfrac)/4.0) // Molecular weight
#define mu_e                    2.0/(1.0+hydrogen_massfrac) // Molecular weight
#define logOH_solar             8.66
#define pi                      3.14159274101257    // pi = 4.0*atan(1.0)
#define G_cgs                   6.67428E-8          // G = 6.67428E-8 [cm^3 g^-1 s^-2]
#define kpc                     3.085678E21         // 1 kpc = 3.085678E21 [cm]
#define kev         6.242e8             // kilo ElectronVolt [erg]
#define unit_mass               AllVars.UnitMass    // 1.989E43 = 1E10 [Solar masses]
#define unit_velocity           AllVars.UnitVelocity // km.s^-1 in [cm.s^-1]
#define unit_length             AllVars.UnitLength  // kpc in [cm]
#define unit_time               (unit_length/unit_velocity)
#define unit_energy             (unit_mass*(unit_length*unit_length)/(unit_time*unit_time))
#define unit_dens               (unit_mass/(unit_length*unit_length*unit_length))
#define unit_nh                 ((hydrogen_massfrac/protonmass)*unit_dens)
#define unit_ne                 ((((1+hydrogen_massfrac)/2.0)/protonmass)*unit_dens)
#define G                       G_cgs*unit_mass/(unit_length*unit_velocity*unit_velocity)


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
        // Spin parameter of the galaxy
        double lambda;
        // Baryonic mass fraction of the disk
        double m_d;
        // Total mass of the galactic model
        double total_mass;
        // Total mass of the galactic model within R200
        double total_mass_r200;
        // Density fluctuation dispersion
        double dens_fluct_sigma;
        // Density fluctuation injection scale
        double dens_fluct_scale_inj;
        // Density fluctuation injection scale
        double dens_fluct_scale_diss;
        // Density fluctuations spectral power index
        double dens_fluct_nspec;
        // Density fluctuation seed
        long dens_fluct_seed;
        // Total number of components
        int n_component;
        int                 *selected_comp;
        int hydro_eq_niter;
        int hydro_eq;
        double maxrad;
        double maxrad_gas;
        // Component quantities
        char                **comp_profile_name;
        char                **comp_imf_name;
        unsigned long int   *comp_npart;
        unsigned long int   *comp_npart_pot;
        unsigned long int   *comp_start_part;
        double              *comp_mass_frac;
        double              *comp_mass;
        int                 *comp_model;
        double              *comp_scale_length;
        double              *comp_concentration;
        double              *comp_scale_height;
        double              *comp_cut;
        double              *comp_sigma_cut;
        double              *comp_sigma_cut_in;
        double              *comp_flatx;
        double              *comp_flaty;
        double              *comp_flatz;
        double              *comp_flatx_cut;
        double              *comp_flaty_cut;
        double              *comp_flatz_cut;
        double              *comp_flatx_out;
        double              *comp_flaty_out;
        double              *comp_flatz_out;
        double              *comp_flatx_rt;
        double              *comp_flaty_rt;
        double              *comp_flatz_rt;
        double              *comp_flatx_st;
        double              *comp_flaty_st;
        double              *comp_flatz_st;
        char                *comp_flatx_var;
        char                *comp_flaty_var;
        char                *comp_flatz_var;
        double              *comp_mcmc_step;
        double              *comp_mcmc_step_slope;
        double              *comp_mcmc_step_hydro;
        double              *comp_vmax_esc;
        double              *comp_vmax_circ;
        int                 *comp_type;
        int                 *comp_bool;
        double              *comp_stream_frac;
        double              *comp_cut_dens;
        double              *comp_theta_sph;
        double              *comp_phi_sph;
        double              *comp_metal;
        double              *comp_metal_sigma;
        double              *comp_metal_scale;
        long                *comp_metal_seed;
        double              *comp_t_init;
        double              *comp_u_init;
        double              *comp_cs_init;
        int                 *comp_turb_gradient;
        double              *comp_turb_sigma;
        double              *comp_turb_frac;
        double              *comp_turb_scale;
        double              *comp_turb_scale_inj;
        double              *comp_turb_scale_diss;
        double              *comp_turb_nspec;
        long                *comp_turb_seed;
        double              *comp_age_sigma;
        double              *comp_age_scale;
        long                *comp_age_seed;
        double              *comp_mean_age;
        double              *comp_min_age;
        double              *comp_alpha_struct;
        double              *comp_beta_struct;
        double              *comp_scale_dens;
        double              *comp_radius_nfw;
        double              *comp_Q_lim;
        double              *comp_Q_min;
        double              *comp_Q_fixed;
        double              *comp_Q_boost;
        double              *comp_Q_bar;
        double              *comp_t_min;
        int                 *comp_compute_vel;
        int                 *comp_hydro_eq;
        int                 *comp_hydro_eq_mode;
        int                 *comp_spherical_hydro_eq;
        int                 *comp_dens_fluct;
        double              *comp_cut_in;
        int                 *comp_thermal_eq;
        double              *comp_part_mass;
        double              *comp_part_mass_pot;
        int                 *comp_jeans_mass_cut;
        int                 *comp_jeans_anisotropy_model;
        int                 *comp_epicycle;
        int                 *comp_metal_gradient;
        int                 *comp_excavate;
        double              *comp_sfr;
        double              *comp_spiral_theta_out;
        double              *comp_spiral_r_in;
        double              *comp_spiral_r_out;
        double              *comp_spiral_alpha;
        double              *comp_warp_scale;
        int                 *comp_warp_mode;
        int                 *comp_jeans_dim;
        int                 *comp_sigmar_model;
        int                 *comp_sigmaz_model;
        double              *comp_sigmar;
        double              *comp_sigmaz;
        double              *comp_sigmar_radius;
        double              *comp_sigmaz_radius;
        double              *comp_sigmar_scale;
        double              *comp_sigmaz_scale;
        double              *comp_jeans_f_sigma;
        double              *comp_k_stream;
        int         *comp_delete;
        double              *comp_stream_scale;
        int                 *comp_stream_method;
        double              *comp_angmom_frac;
        double              *comp_dens_min;
        double              *comp_dens_max;
        double              *comp_gamma_poly;
        double              *comp_k_poly;
        double              *comp_dens_init;
        double              *comp_accept_min;
        double              *comp_accept_max;
        double              *comp_rcore;
        double      *comp_ggd_beta;
        double      *comp_softening;
        double      *comp_rc_entropy;
        double      *comp_alpha_entropy;
        double      *comp_cut_hydro_eq;
        int                 *comp_symmetry;
        int                 *comp_imf_model;
        double              *comp_mstar_min;
        double              *comp_mstar_max;
        double              *comp_mcmc_step_mass;
        // Virial quantities
        double v200;
        double r200;
        double m200;
        double s200;
        // Total angular momentum within R200
        double J200;
        // Coordinates vectors
        double              *x;
        double              *y;
        double              *z;
        double              *r_cyl;
        double              *theta_cyl;
        double              *r_sph;
        double              *phi_sph;
        // Velocities vectors
        double              *vel_x;
        double              *vel_y;
        double              *vel_z;
        // Mass vector
        double              *mass;
        // Density vector
        double              *rho;
        // Internal energy vector
        double              *u;
        // Age vector
        double              *age;
        // Metallicity vector
        double              *metal;
        // Particle Mesh potential grid
        double              ****potential;
        double              ****potential_ext;
        double      **vr2_tilted_mixed;
        double      **vz2_tilted_mixed;
        double      **vtheta2_mixed;
        // Particle Mesh gaussian field grid
        double              ***gaussian_field;
        // Gas midplane density grid
        double              **midplane_dens;
        // Potential grid cell size vector [kpc]
        double *dx;
        // Gas midplane density grid cell size vector [kpc]
        double dx_dens;
        double dx_jeans;
        // Gaussian field grid cell size vector [kpc]
        double dx_gauss;
        // Number of particles per particle type
        unsigned long int num_part[4];
        unsigned long int num_part_pot[4];
        // Total potential grid size [kpc]
        double *boxsize;
        double *boxsize_flatx;
        double *boxsize_flaty;
        double *boxsize_flatz;
        // Total gas midplane density grid size [kpc]
        double boxsize_dens;
        double boxsize_jeans;
        // Level of refinement of the potential grid
        int level_coarse;
        int nlevel;
        // Level of refinement of the gas density grid
        int level_grid_dens;
        int level_grid_jeans_3D;
        // Level of refinement of the gas turbulence grid
        int level_grid_turb;
        // Level of refinement of the stars age grid
        int level_grid_age;
        // Level of refinement of the metal grid
        int level_grid_metal;
        // Level of refinement of the density gaussian fluctuations
        int level_grid_dens_fluct;
        // Number of cells in the potential grid
        int **ngrid;
        // Number of cells in the midplane density grid
        int ngrid_dens[2];
        int ngrid_jeans[2];
        // Number of cells in the turbulence grid
        int ngrid_gauss[3];
        // Storage array
        double              **storage;
        // Identifier of particle
        unsigned long int   *index;
        unsigned long int   *id;
        // Total number of particles
        unsigned long int ntot_part;
        unsigned long int ntot_part_stars;
        unsigned long int ntot_part_pot;
        // Boolean variable checking the computation of the potential
        int potential_defined;
        // Boolean variable checking the computation of the midplane gas density
        int midplane_dens_defined;
        // Boolean variable checking the computation of the gaussian field
        int gaussian_field_defined;
        int jeans_3D_defined;
        // Seed for random number generator
        long seed;
        // Pseudo density boolean
        int *pseudo;
        // Gravitational softening
        double softening;
        // MCMC multiple try parameter
        int mcmc_ntry;
        // Main halo component index
        int index_halo;
        // Main disk component index
        int index_disk;
        // Main gas disk component index
        int index_gasdisk;
        // First valid component index
        int index_first;
        // Copy tag
        int copy;
        // Masses per type
        double gas_mass;
        double halo_mass;
        double disk_mass;
        double bulge_mass;
        double stellar_mass;
} galaxy;

// This is a type definition of a stream. The thought here is to create streams
// as objects and, hopefully, make the code extremely clean. It is essentially a
// a collection of arrays and constants.
typedef struct {
        int                 *selected_comp;
        // Stream's peak coordinates
        double              *comp_xc;
        double              *comp_yc;
        double              *comp_zc;
        double              *comp_mass;
        double              *comp_dens;
        double              *comp_opening_angle;
        int                 *comp_turb_gradient;
        double              *comp_turb_sigma;
        double              *comp_turb_scale_inj;
        double              *comp_turb_scale_diss;
        double              *comp_turb_nspec;
        long                *comp_turb_seed;
        double              *comp_t_init;
        // Stream's spin orientation spherical angles
        double              *comp_theta_sph;
        double              *comp_phi_sph;
        double              *comp_length;
        double              *comp_scale;
        unsigned long int   *comp_npart;
        int                 *comp_bool;
        int                 *comp_model;
        double              *comp_mcmc_step;
        double              *comp_metal;
        int                 *comp_metal_gradient;
        double              *comp_u_init;
        double              *comp_cs_init;
        unsigned long int   *comp_start_part;
        char                **comp_profile_name;
        int                 *comp_dens_fluct;
        double              *comp_dens_min;
        double              *comp_dens_max;
        double              *comp_k_poly;
        double              *comp_gamma_poly;
        double              *comp_accept_min;
        double              *comp_accept_max;
        // Coordinates vectors
        double              *x;
        double              *y;
        double              *z;
        double              *r_cyl;
        double              *theta_cyl;
        double              *r_sph;
        double              *phi_sph;
        // Mass vector
        double              *mass;
        // Density vector
        double              *rho;
        // Internal energy vector
        double              *u;
        // Metallicity vector
        double              *metal;
        // Velocities vectors
        double              *vel_x;
        double              *vel_y;
        double              *vel_z;
        unsigned long int   *id;
        // Level of refinement of the gas turbulence grid
        int level_grid_turb;
        // Level of refinement of the density gaussian fluctuations
        int level_grid_dens_fluct;
        // Number of cells in the turbulence grid
        int ngrid_gauss[3];
        // Gaussian field grid cell size vector [kpc]
        double dx_gauss;
        // Particle Mesh gaussian_field grid
        double              ***gaussian_field;
        // Storage array
        double              **storage;
        unsigned long int ntot_part;
        double total_mass;
        int n_component;
        long seed;
        // Density fluctuation dispersion
        double dens_fluct_sigma;
        // Density fluctuation scale
        double dens_fluct_scale_inj;
        double dens_fluct_scale_diss;
        double dens_fluct_nspec;
        // Density fluctuation seed
        long dens_fluct_seed;
        // Boolean variable checking the computation of the gaussian field
        int gaussian_field_defined;
        // MCMC multiple try parameter
        int mcmc_ntry;
} stream;

//Gadget2-style header for Gadget2 snapshots.
struct io_header_1 {
        //char fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8];
        int npart[6];                             /*!< number of particles of each type in this file */
        double mass[6];                           /*!< mass of particles of each type. If 0, then the masses are explicitly
                                                     stored in the mass-block of the snapshot file, otherwise they are omitted */
        double time;                              /*!< time of snapshot file */
        double redshift;                          /*!< redshift of snapshot file */
        int flag_sfr;                             /*!< flags whether the simulation was including star formation */
        int flag_feedback;                        /*!< flags whether feedback was included (obsolete) */
        unsigned int npartTotal[6];               /*!< total number of particles of each type in this snapshot. This can be
                                                     different from npart if one is dealing with a multi-file snapshot. */
        int flag_cooling;                         /*!< flags whether cooling was included  */
        int num_files;                            /*!< number of files in multi-file snapshot */
        double BoxSize;                           /*!< box-size of simulation in case periodic boundaries were used */
        double Omega0;                            /*!< matter density in units of critical density */
        double OmegaLambda;                       /*!< cosmological constant parameter */
        double HubbleParam;                       /*!< Hubble parameter in units of 100 km/sec/Mpc */
        int flag_stellarage;                      /*!< flags whether the file contains formation times of star particles */
        int flag_metals;                          /*!< flags whether the file contains metallicity values for gas and star particles */
        unsigned int npartTotalHighWord[6];       /*!< High word of the total number of particles of each type */
        int flag_entropy_instead_u;               /*!< flags that IC-file contains entropy instead of u */
        char fill[60];
} header1;

//Gadget2-style particle data structure.
struct particle_data {
        int Type;
        int Id;
        float Pos[3];
        float Vel[3];
        float Mass;
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
        char ParameterFile[MAXLEN_FILENAME];
        char GalaxyFiles[MAX_GAL][MAXLEN_FILENAME];
        char StreamFiles[MAX_STREAM][MAXLEN_FILENAME];
        char Filename[MAXLEN_FILENAME];
        char ICformat[MAXLEN_FILENAME];
        char UnitMassName[20];
        char UnitVelocityName[20];
        char UnitLengthName[20];
        int Kepler_Gal1[MAX_GAL];
        int Kepler_Gal2[MAX_GAL];
        int Kepler_GalCenter[MAX_GAL];
        int Circular_Gal1[MAX_GAL];
        int Circular_Gal2[MAX_GAL];
        int Circular_GalCenter[MAX_GAL];
        int GalId[MAX_GAL];
        int StreamId[MAX_STREAM];
        int Ngal;
        int Nstream;
        int Nthreads;
        int MeanPartDist;
        int OutputRc;
        int OutputGasRc;
        int OutputPot;
        int OutputRho;
        int OutputSigma;
        int OutputToomre;
        int MaxCompNumber;
        int MaxNlevel;
        int GaussianRejectIter;
        int CurrentGalaxy;
        int NormMassFact;
        int GslIntegrationScheme;
        unsigned long int GalStart[MAX_GAL];
        unsigned long int GalNpart[MAX_GAL];
        unsigned long int StreamStart[MAX_STREAM];
        unsigned long int StreamNpart[MAX_STREAM];
        unsigned long int GslWorkspaceSize;
        double StreamSpin[MAX_STREAM];
        double StreamIncl[MAX_STREAM];
        double StreamPos[MAX_STREAM][3];
        double GalSpin[MAX_GAL];
        double GalIncl[MAX_GAL];
        double GalMass[MAX_GAL];
        double GalPos[MAX_GAL][3];
        double GalVel[MAX_GAL][3];
        double Kepler_Ecc[MAX_GAL];
        double Kepler_Rinit[MAX_GAL];
        double Kepler_Rperi[MAX_GAL];
        double Kepler_OrbitPlaneTheta[MAX_GAL];
        double Kepler_OrbitPlanePhi[MAX_GAL];
        double Circular_Rinit[MAX_GAL];
        double Circular_OrbitPlaneTheta[MAX_GAL];
        double Circular_OrbitPlanePhi[MAX_GAL];
        double Circular_Vc[MAX_GAL];
        double redshift;
        double H0;
        double H;
        double h;
        double Omega_m;
        double Omega_l;
        double Omega_k;
        double UnitMass;
        double UnitVelocity;
        double UnitLength;
} AllVars;

// ------------------------------------------
// These are the function prototypes for DICE.
// ------------------------------------------
// Initialization and destruction functions
int allocate_component_arrays(galaxy *);
int allocate_variable_arrays(galaxy *);
int allocate_galaxy_potential(galaxy *);
int allocate_galaxy_midplane_dens(galaxy *);
int reallocate_variable_arrays(galaxy *, unsigned long int);
int allocate_component_arrays_stream(stream *);
int allocate_variable_arrays_stream(stream *);
int allocate_galaxy_gaussian_grid(galaxy *);
int deallocate_galaxy_gaussian_grid(galaxy *);
int allocate_stream_gaussian_grid(stream *);
int deallocate_stream_gaussian_grid(stream *);
int create_galaxy(galaxy *, char *, int);
int create_stream(stream *, char *, int);
void allocate_galaxy_storage_variable(galaxy *, int);
void allocate_stream_storage_variable(stream *, int);
void allocate_dispersion_grid();
void destroy_dispersion_grid();
int set_galaxy_coords(galaxy *);
int set_galaxy_velocity(galaxy *);
void trash_galaxy(galaxy *, int);
void trash_stream(stream *, int);
int add_galaxy_to_system(galaxy *, galaxy *);
int rotate_galaxy(galaxy *, int);
int rotate_component(galaxy *, double, double, int);
int rotate_all(galaxy *, double, double, int);
int set_galaxy_trajectory(galaxy *, int);
int copy_galaxy(galaxy *, galaxy *, int);
void set_orbit_keplerian(int, int, double, double, double, double, double, int);
void set_orbit_circular(int, int, double, double, double, double, int);
int add_stream_to_system(stream *, galaxy *);
int stream_to_galaxy(stream *, galaxy *);
int set_stream_coords(stream *);
int set_stream_velocity(stream *);
int rotate_stream(galaxy *, int);
int position_stream(galaxy *, int);

// Structure functions
double density_functions_pool(galaxy *, double, double, double, int, int, int);
double imf_functions_pool(galaxy *, double, int, int);
void mcmc_metropolis_hasting_ntry(galaxy *, int, int);
void mcmc_metropolis_hasting_ntry_stream(stream *, int, int);
void mcmc_metropolis_hasting_ntry_mass(galaxy *, int, int);
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
double pseudo_density_gas_func(galaxy *, double, double, double, int, int, int, int);
double midplane_density_gas_func(galaxy *, gsl_integration_workspace *, double, double);
static double dmidplane_density_gas_func(double, void *);
void fill_midplane_dens_grid(galaxy *);
double get_midplane_density(galaxy *, double, double);

double disk_scale_length_func(galaxy *, double, int);
double disk_scale_length_obs_func(galaxy *, int);
double f_c_func(double);
double g_c_func(double);
static double dg_c_func(double, void *);
double f_s_func(double, double);
double f_r_func(galaxy *, int);
double f_v_func(galaxy *, int);
double mean_interparticle_distance(galaxy *, int);
void lower_resolution(galaxy *);

double halo_concentration(double, double);
double halo_abundance(double, double);
double galaxy_metallicity(double, double);

// Velocity functions
double v_c_func(galaxy *,double);
double v2a_r_1D_func(galaxy *, gsl_integration_workspace *, int);
double v2a_r_2D_func(galaxy *, gsl_integration_workspace *, int);
double v2a_theta_2D_func(galaxy *, double, double, double, int);
double rho_v2a_r_2D_func(double, void *);
double v2a_r_3D_func(galaxy *);
double v2a_z_3D_func(galaxy *);
double v2a_theta_3D_func(galaxy *);
static double dv2a_r_1D_func(double, void *);
static double dv2a_r_2D_func(double, void *);
void fill_jeans_3D_grid (galaxy *, int);
double get_jeans_3D_cic(galaxy *, double, double, double **);
double jeans_3D_h_func(galaxy *, double, double, int);
double sigma2_theta_epicycle_func(galaxy *, double, double);
double v2a_r_toomre(galaxy *, double, double, int);
double toomre(galaxy *, double, double, int);
double v2_theta_gas_func(galaxy *, double, double, int);
double density_wrapper_func(double, void *);
int set_turbulent_grid(galaxy *, int);
double galaxy_turbulence_func(galaxy *, double, double, double, int);
double Jtot_func(galaxy *, int);

// Potential and force functions
int set_galaxy_potential(galaxy *, double ***, double, int [3], int, double[3]);
double galaxy_potential_func(galaxy *, double ***, double, int [3], double, double, double, int);
double galaxy_zforce_func(galaxy *, double);
double galaxy_rforce_func(galaxy *, double);
double galaxy_rsphforce_func(galaxy *, double);
double galaxyr_potential_wrapper_func(double, void *);
double galaxyrsph_potential_wrapper_func(double, void *);
double galaxyz_potential_wrapper_func(double, void *);
double potential_deriv_wrapper_func(double, void *);
void copy_potential(galaxy *, galaxy *, int);
double galaxy_total_potential(galaxy *, double, double, double, int, int);
double get_h_value(galaxy *, double, double, double, int, int);
double set_galaxy_potential_all(galaxy *, int);

// Input, output, and manipulation functions
void write_dice_version();
int parse_config_file(char *);
int parse_galaxy_file(galaxy *, char *);
int parse_stream_file(stream *, char *);
void write_galaxy_rz_quantities(galaxy *, double, char *, double);
int write_gadget1_ics(galaxy *, char *);
int write_gadget2_ics(galaxy *, char *);
int load_snapshot(char *, int);
int allocate_memory(int);
int reordering(int, int *);
int unload_snapshot();

// Toolbox functions and definitions
double min(double, double);
double max(double, double);
double sum_dbl(double *, int);
typedef double (*function_to_derivate)(double, void *);
double deriv_central2(galaxy *, double, double, function_to_derivate);
double deriv_central4(galaxy *, double, double, function_to_derivate);
double deriv_forward(galaxy *, double, double, function_to_derivate);
double standard_dev(double *, unsigned long int, unsigned long int);
double interpol(double, double, double, double, double);
double smooth_in(double, double, double);
double smooth_out(double, double, double);
int set_galaxy_random_field_grid(galaxy *, double, double, double, long);
int set_stream_random_field_grid(stream *, double, double, double, long);
double galaxy_gaussian_field_func(galaxy *, double, double, double);
double stream_gaussian_field_func(stream *, double, double, double);
