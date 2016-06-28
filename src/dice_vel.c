/*-----------------------------------------------------------------------------
  /
  / Filename: dice_vel.c
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

// This function calculates the z-axis velocity moment at a given radius
double v2a_r_func(galaxy *gal, gsl_integration_workspace *w, int component) {

    int status, tid, save;
    double integral, error, res, rho, infinity, v2a_r;
    size_t neval,limit;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    gsl_function F;
    F.function = &dv2a_z_func;
    F.params = gal;
    gal->selected_comp[tid] = component;

    // Dispersion from Jeans equation with non zero mixed moment <vr.vz>
    if(gal->comp_sigmar_model[component]==-2) {
        v2a_r = get_jeans_array_cic(gal,fabs(gal->r_cyl[gal->index[tid]]),fabs(gal->z[gal->index[tid]]),gal->vr2_tilted_mixed);
        // Dispersion from Jeans equation
    } else if(gal->comp_sigmar_model[component]==-1) {

        infinity = gal->comp_cut[gal->selected_comp[tid]];
        switch(AllVars.GslIntegrationScheme){
            case 1:
                gsl_integration_qag(&F,fabs(gal->z[gal->index[tid]]),fabs(gal->z[gal->index[tid]])+infinity,epsabs,epsrel,AllVars.GslWorkspaceSize,key,w,&integral,&error);
                break;
            case 2:
                gsl_integration_qagiu(&F,fabs(gal->z[gal->index[tid]]),epsabs,epsrel,AllVars.GslWorkspaceSize,w,&integral,&error);
                break;
            case 3:
                gsl_integration_qng(&F,fabs(gal->z[gal->index[tid]]),fabs(gal->z[gal->index[tid]])+infinity,epsabs,epsrel,&integral,&error,&neval);
                break;
            default:
                fprintf(stderr,"[Error] gsl_integration_scheme=%d is not a valid value\n",AllVars.GslIntegrationScheme);
                exit(0);
        }
        rho = density_functions_pool(gal,gal->r_cyl[gal->index[tid]],gal->theta_cyl[gal->index[tid]],gal->z[gal->index[tid]],0,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
        v2a_r = integral/rho;
        v2a_r = (v2a_r>0. ? v2a_r : 0.);
        // Dispersion from Isothermal sheet
    } else if(gal->comp_sigmar_model[component]==0) {
        v2a_r = pi*G*surface_density_func(gal,gal->r_cyl[gal->index[tid]],gal->theta_cyl[gal->index[tid]],1,component)*gal->comp_scale_length[component]*gal->comp_flatz[component];
        // Dispersion proportional to surface density
    } else if (gal->comp_sigmar_model[component]>0) {
        save = gal->comp_model[gal->selected_comp[tid]];
        gal->comp_model[gal->selected_comp[tid]] = gal->comp_sigmar_model[component];
        v2a_r = gal->comp_sigmar_scale[component]*surface_density_func(gal,gal->r_cyl[gal->index[tid]],gal->theta_cyl[gal->index[tid]],1,gal->selected_comp[tid]);
        gal->comp_model[gal->selected_comp[tid]] = save;
    }

    return v2a_r;
}

// This function calculates the z-axis velocity moment at a given radius
double v2a_z_func(galaxy *gal, gsl_integration_workspace *w, int component) {

    int status, tid, save;
    double integral, error, res, rho, infinity, v2a_z;
    size_t neval,limit;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    gsl_function F;
    F.function = &dv2a_z_func;
    F.params = gal;
    gal->selected_comp[tid] = component;

    // Dispersion from Jeans equation with non zero mixed moment <vr.vz>
    if(gal->comp_sigmaz_model[component]==-2) {
        v2a_z = get_jeans_array_cic(gal,fabs(gal->r_cyl[gal->index[tid]]),fabs(gal->z[gal->index[tid]]),gal->vz2_tilted_mixed);
        // Dispersion from Jeans equation
    } else if(gal->comp_sigmaz_model[component]==-1) {

        infinity = gal->comp_cut[gal->selected_comp[tid]];
        switch(AllVars.GslIntegrationScheme){
            case 1:
                gsl_integration_qag(&F,fabs(gal->z[gal->index[tid]]),fabs(gal->z[gal->index[tid]])+infinity,epsabs,epsrel,AllVars.GslWorkspaceSize,key,w,&integral,&error);
                break;
            case 2:
                gsl_integration_qagiu(&F,fabs(gal->z[gal->index[tid]]),epsabs,epsrel,AllVars.GslWorkspaceSize,w,&integral,&error);
                break;
            case 3:
                gsl_integration_qng(&F,fabs(gal->z[gal->index[tid]]),fabs(gal->z[gal->index[tid]])+infinity,epsabs,epsrel,&integral,&error,&neval);
                break;
            default:
                fprintf(stderr,"[Error] gsl_integration_scheme=%d is not a valid value\n",AllVars.GslIntegrationScheme);
                exit(0);
        }
        rho = density_functions_pool(gal,gal->r_cyl[gal->index[tid]],gal->theta_cyl[gal->index[tid]],gal->z[gal->index[tid]],0,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
        v2a_z = integral/rho;
        v2a_z = (v2a_z>0. ? v2a_z : 0.);
        // Dispersion from Isothermal sheet
    } else if(gal->comp_sigmaz_model[component]==0) {
        v2a_z = pi*G*surface_density_func(gal,gal->r_cyl[gal->index[tid]],gal->theta_cyl[gal->index[tid]],1,component)*gal->comp_scale_length[component]*gal->comp_flatz[component];
        // Dispersion proportional to surface density
    } else if (gal->comp_sigmaz_model[component]>0) {
        save = gal->comp_model[gal->selected_comp[tid]];
        gal->comp_model[gal->selected_comp[tid]] = gal->comp_sigmaz_model[component];
        v2a_z = gal->comp_sigmaz_scale[component]*surface_density_func(gal,gal->r_cyl[gal->index[tid]],gal->theta_cyl[gal->index[tid]],1,gal->selected_comp[tid]);
        gal->comp_model[gal->selected_comp[tid]] = save;
    }

    return v2a_z;
}


// This is the integrand for the previous function. It is setup to work with
// the GSL_qags structures.
static double dv2a_z_func(double z, void *params) {

    double integrand, radius, theta, rho;
    int tid;

    galaxy *gal = (galaxy *) params;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    radius = gal->r_cyl[gal->index[tid]];
    theta = gal->theta_cyl[gal->index[tid]];

    rho = density_functions_pool(gal,radius,theta,z,gal->comp_jeans_mass_cut[gal->selected_comp[tid]],gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);

    integrand = rho*galaxy_zforce_func(gal,z);

    return integrand;
}

// This function calculates the 1D velocity moment at a given spherical radius
double v2a_1D_func(galaxy *gal, gsl_integration_workspace *w, int component) {

    int status, tid;
    double integral, error, v2a, rho, infinity;
    size_t neval;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    gsl_function F;

    F.function = &dv2a_1D_func;
    F.params = gal;

    gal->selected_comp[tid] = component;
    infinity = 10.*gal->comp_cut[gal->selected_comp[tid]];
    //gsl_integration_qag(&F,fabs(gal->r_sph[gal->index[tid]]),fabs(gal->r_sph[gal->index[tid]])+infinity,epsabs,epsrel,AllVars.GslWorkspaceSize,key,w,&integral,&error);
    gsl_integration_qng(&F,fabs(gal->r_sph[gal->index[tid]]),fabs(gal->r_sph[gal->index[tid]])+infinity,epsabs,epsrel,&integral,&error,&neval);

    rho = density_functions_pool(gal,gal->r_cyl[gal->index[tid]],gal->theta_cyl[gal->index[tid]],gal->z[gal->index[tid]],0,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
    v2a = integral/rho;
    v2a = v2a>0.?v2a:0.;

    return v2a;

}

// This is the integrand for the previous function. It is setup to work with
// the GSL_qags structures.
static double dv2a_1D_func(double r_sph, void *params) {

    double integrand, r_cyl, z, theta_cyl, rho;
    int tid;
    galaxy *gal = (galaxy *) params;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    theta_cyl = gal->theta_cyl[gal->index[tid]];
    z = cos(gal->phi_sph[gal->index[tid]])*r_sph;
    r_cyl = sqrt(r_sph*r_sph-z*z);

    rho = density_functions_pool(gal,r_cyl,theta_cyl,z,gal->comp_jeans_mass_cut[gal->selected_comp[tid]],gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);

    integrand = rho*galaxy_rsphforce_func(gal,r_sph);

    return integrand;
}

// This function calculates PART of the phi axis velocity moment. In
// particular, it calculates the derivative of the radial velocity
// moment, which is equal to the z-axis velocity moment.
//
// This function uses the galaxy storage variables a lot! If you have
// edited the code to ill effect, check that the storage variables
// properly reassigned after this function.
//
// The derivative is calculated in the x,y plane using the
// 5-point stencil method.
double v2a_theta_func(galaxy *gal, double radius, double v2a_r, double v_c, int component) {

    double z, theta, h, v2a_theta, rho, abserr;
    double derivative, x, y;
    int tid;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    z = gal->z[gal->index[tid]];
    if(gal->comp_sigmar_model[component]==0) {
        v2a_theta = get_jeans_array_cic(gal,radius,z,gal->vtheta2_mixed);
    } else {
        x = radius*cos(gal->theta_cyl[gal->index[tid]]);
        y = radius*sin(gal->theta_cyl[gal->index[tid]]);
        // Set the derivative stepsize.
        h = get_h_value(gal,x,y,gal->z[gal->index[tid]],0,0);
        theta = gal->theta_cyl[gal->index[tid]];
        gal->selected_comp[tid] = component;
        derivative = deriv_central2(gal,radius,h,rho_v2a_r_func);
        rho = density_functions_pool(gal,radius,gal->theta_cyl[gal->index[tid]],z,0,gal->comp_model[component],component);
        v2a_theta = derivative*radius/rho+v2a_r+pow(v_c,2.0);
    }
    // Reject complex values values
    v2a_theta = v2a_theta<0.?0.:v2a_theta;

    return v2a_theta;
}

// Function which should be derivated to obtain the phi axis velocity moment
// Take a look to the part of the documentation explaining the velocity dispersion computation
double rho_v2a_r_func(double radius, void *params) {

    galaxy *gal = (galaxy *) params;
    double rho, v2a_r;
    double save1;
    int tid;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    rho = density_functions_pool(gal,radius,gal->theta_cyl[gal->index[tid]],gal->z[gal->index[tid]],gal->comp_jeans_mass_cut[gal->selected_comp[tid]],gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
    save1 = gal->r_cyl[gal->index[tid]];
    gal->r_cyl[gal->index[tid]] = radius;
    v2a_r = v2a_r_func(gal,w[tid],gal->selected_comp[tid]);
    gal->r_cyl[gal->index[tid]] = save1;

    return rho*v2a_r;
}

void fill_jeans_mixed_grid (galaxy *gal, int component) {
    int i, j, k, n, nsteps, tid;
    double *qprev, alpha, h1, h2, r1, r2, dr, z1, z2, z, dz, rho, rho1, rho2, p, f, dqdr;
    double rmin,fac, save1, save2, save3, save4, save5;
    double vrvz, vr2, vr2_tilted, vz2_tilted;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    if (!(qprev = calloc(gal->ngrid_jeans[0],sizeof(double *)))) {
        fprintf(stderr,"[Error] Unable to create vz mixed temporary array.\n");
        return ;
    }

    nsteps = 100;
    gal->index[tid] = 0;
    save1 = gal->z[gal->index[tid]];
    save2 = gal->theta_cyl[gal->index[tid]];
    save3 = gal->r_cyl[gal->index[tid]];
    save4 = gal->x[gal->index[tid]];
    save5 = gal->y[gal->index[tid]];
    gal->z[gal->index[tid]] = 0.;
    gal->theta_cyl[gal->index[tid]] = 0.;


    for (k=gal->ngrid_jeans[0]-1; k>=0; k--) {
        if(k==gal->ngrid_jeans[0]-1) {
            for (j=0; j<gal->ngrid_jeans[1]; j++) {
                gal->vz2_tilted_mixed[k][j] = 0.;
            }
        } else {
            z1 = ((double)(k+1)+0.5)*gal->dx_jeans;
            z2 = ((double)k+0.5)*gal->dx_jeans;
            dz = z2-z1;
            for (j=0; j<gal->ngrid_jeans[1]; j++) {
                qprev[j] = gal->vz2_tilted_mixed[k+1][j];
            }
            for (n=0; n<nsteps; n++) {
                z = z1+(dz/nsteps)*n;
                for (j=0; j<gal->ngrid_jeans[1]; j++) {
                    r1 = ((double)j+0.5)*gal->dx_jeans;
                    gal->r_cyl[gal->index[tid]] = r1;
                    gal->x[gal->index[tid]] = gal->r_cyl[gal->index[tid]]*cos(gal->theta_cyl[gal->index[tid]]);
                    gal->y[gal->index[tid]] = gal->r_cyl[gal->index[tid]]*sin(gal->theta_cyl[gal->index[tid]]);
                    gal->z[gal->index[tid]] = z;

                    h1 = h_function(gal,r1,z,component);
                    rho = density_functions_pool(gal,r1,0.,z,gal->comp_jeans_mass_cut[component],gal->comp_model[component],component);
                    p = rho*galaxy_zforce_func(gal,z);

                    if(j==gal->ngrid_jeans[1]-1){
                        dqdr = 0.;
                    } else {
                        r2 = ((double)(j+1)+0.5)*gal->dx_jeans;
                        h2 = h_function(gal,r2,z,component);
                        dr = r2-r1;
                        dqdr = (h2*qprev[j+1]-h1*qprev[j])/dr+qprev[j]*h1/r1;
                    }
                    gal->vz2_tilted_mixed[k][j] = qprev[j]+(-p-dqdr)*(dz/nsteps);
                }
                for (j=0; j<gal->ngrid_jeans[1]; j++) qprev[j] = gal->vz2_tilted_mixed[k][j];
            }
        }
    }

    for (k=0; k<gal->ngrid_jeans[0]; k++) {
        for (j=0; j<gal->ngrid_jeans[1]; j++) {
            r1 = ((double)j+0.5)*gal->dx_jeans;
            z1 = ((double)k+0.5)*gal->dx_jeans;
            alpha = atan(z1/r1);
            f = gal->comp_f_sigma[component];
            h1 = h_function(gal,r1,z1,component);
            gal->r_cyl[gal->index[tid]] = r1;
            gal->x[gal->index[tid]] = gal->r_cyl[gal->index[tid]]*cos(gal->theta_cyl[gal->index[tid]]);
            gal->y[gal->index[tid]] = gal->r_cyl[gal->index[tid]]*sin(gal->theta_cyl[gal->index[tid]]);
            gal->z[gal->index[tid]] = z1;
            rho = density_functions_pool(gal,r1,0.,z1,gal->comp_jeans_mass_cut[component],gal->comp_model[component],component);
            // sigma_z^2
            gal->vz2_tilted_mixed[k][j] = rho>0. ? gal->vz2_tilted_mixed[k][j]/rho : 0.;
            // <v_z.v_r>
            vrvz = gal->vz2_tilted_mixed[k][j]*((f-1)/2*tan(2*alpha))/
                (pow(cos(alpha),2)-f*pow(sin(alpha),2)+(1.0+f)/2.0*sin(2*alpha)*tan(2*alpha));
            // sigma_r^2
            vr2 = gal->vz2_tilted_mixed[k][j]*(f*pow(cos(alpha),2)-pow(sin(alpha),2)+(1.0+f)/2.0*sin(2*alpha)*tan(2*alpha))/
                (pow(cos(alpha),2)-f*pow(sin(alpha),2)+(1.0+f)/2.0*sin(2*alpha)*tan(2*alpha));
            // <v_r^2>
            vr2_tilted = vr2*pow(cos(alpha),2)+2*vrvz*sin(alpha)*cos(alpha)+gal->vz2_tilted_mixed[k][j]*pow(sin(alpha),2);
            // <v_z^2>
            vz2_tilted = vr2*pow(sin(alpha),2)-2*vrvz*sin(alpha)*cos(alpha)+gal->vz2_tilted_mixed[k][j]*pow(cos(alpha),2);
            gal->vr2_mixed[k][j] = vr2;
            gal->vz2_tilted_mixed[k][j] = vz2_tilted;
            gal->vr2_tilted_mixed[k][j] = vr2_tilted;
        }
    }
    for (k=gal->ngrid_jeans[0]-1; k>=0; k--) {
        for (j=0; j<gal->ngrid_jeans[1]; j++) {
            r1 = ((double)j+0.5)*gal->dx_jeans;
            z1 = ((double)k+0.5)*gal->dx_jeans;
            gal->r_cyl[gal->index[tid]] = r1;
            gal->x[gal->index[tid]] = gal->r_cyl[gal->index[tid]]*cos(gal->theta_cyl[gal->index[tid]]);
            gal->y[gal->index[tid]] = gal->r_cyl[gal->index[tid]]*sin(gal->theta_cyl[gal->index[tid]]);
            gal->z[gal->index[tid]] = z1;

            double fz = galaxy_rforce_func(gal,z);
            rho1 = density_functions_pool(gal,r1,0.,z1,gal->comp_jeans_mass_cut[component],gal->comp_model[component],component);

            int j2;
            if(j<gal->ngrid_jeans[1]-2) {
                r2 = ((double)(j+1)+0.5)*gal->dx_jeans;
                j2 = j+1;
            } else {
                r2 = ((double)(j-1)+0.5)*gal->dx_jeans;
                j2 = j-1;
            }
            rho2 = density_functions_pool(gal,r2,0.,z1,gal->comp_jeans_mass_cut[component],gal->comp_model[component],component);
            double vtheta2 = 0.;
            if(rho1>0.) {
                vtheta2 = gal->vr2_tilted_mixed[k][j]+r1*fz+
                    r1/rho1*(rho2*gal->vr2_tilted_mixed[k][j2]-rho1*gal->vr2_tilted_mixed[k][j])/(r2-r1)+
                    r1/rho1*(rho2*h_function(gal,r2,z,component)*gal->vz2_tilted_mixed[k][j2]-rho1*h_function(gal,r1,z1,component)*gal->vz2_tilted_mixed[k][j])/(r2-r1);
            }
            if(vtheta2>0) {
                gal->vtheta2_mixed[k][j] = vtheta2;
            } else {
                gal->vtheta2_mixed[k][j] = 0.;
            }
        }
    }
    gal->z[gal->index[tid]] = save1;
    gal->theta_cyl[gal->index[tid]] = save2;
    gal->r_cyl[gal->index[tid]] = save3;
    gal->x[gal->index[tid]] = save4;
    gal->y[gal->index[tid]] = save5;

    return;
}

double h_function(galaxy *gal, double r, double z, int component) {
    double f, alpha, h;
    int tid;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    alpha = atan(z/r);
    f = gal->comp_f_sigma[component];
    h = ((f-1)/2*tan (2*alpha))/(pow(cos(alpha),2)-f*pow(sin(alpha),2)+(1.0+f)/2.0*sin(2*alpha)*tan(2*alpha));

    return h;
}

double get_jeans_array_cic(galaxy *gal, double r, double z, double **array) {
    int node_r,node_z;
    double rtemp,ztemp,v1,v2,v3,v4,v_interp;
    double dr,dz,tr,tz;

    rtemp = r/gal->dx_jeans;
    ztemp = z/gal->dx_jeans;
    // Determine the parent node.
    node_r = floor(rtemp);
    node_z = floor(ztemp);

    if(node_r>=0 && node_r<gal->ngrid_jeans[0]-1 && node_z>=0 && node_z<gal->ngrid_jeans[1]-1) {
        // Interpolation function to compute the potential.
        v1 = array[node_r][node_z];
        v2 = array[node_r+1][node_z];
        v3 = array[node_r][node_z+1];
        v4 = array[node_r+1][node_z+1];
        // CIC fractions
        dr = 1.0 - (rtemp - (double) node_r);
        dz = 1.0 - (ztemp - (double) node_z);
        tr = 1.0 - dr;
        tz = 1.0 - dz;
        // Return the interpolated potential.
        v_interp = dr*dz*v1 + tr*dz*v2 + dr*tz*v3 + tr*tz*v4;
    } else {
        v_interp = 0.;
    }
    return v_interp;

}


// This function checks that the azimuthal velocity dispersion verifies
// the lower limit of the Toomre stability cirterion imposed by the user
double v2a_r_toomre(galaxy *gal, double radius, double v2a_r, int component) {

    double gamma_sqrd, kappa_sqrd, h, force, dforcedr, Q, surface_density;
    double abserr, x, y, v2a_r_new;
    int tid;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    if(radius==0.) return 0.;

    v2a_r_new = v2a_r;
    x = radius*cos(gal->theta_cyl[gal->index[tid]]);
    y = radius*sin(gal->theta_cyl[gal->index[tid]]);
    // Set the derivative step
    h = get_h_value(gal,x,y,gal->z[gal->index[tid]],0,0);
    // Calculate force and force derivative
    force = potential_deriv_wrapper_func(radius,gal);
    dforcedr = deriv_central2(gal,radius,h,potential_deriv_wrapper_func);

    kappa_sqrd = 3.0*force/(radius) + dforcedr;
    gamma_sqrd = 4.0*force/(kappa_sqrd*radius);
    // We take into account the contribution of the gas to the stability of the system
    surface_density = surface_density_func(gal,radius,gal->theta_cyl[gal->index[tid]],1,component);
    if(gal->comp_type[component]==0) {
        Q = sqrt(gal->comp_cs_init[component]*fabs(kappa_sqrd))/(pi*G*surface_density);
    } else {
        Q = sqrt(v2a_r*fabs(kappa_sqrd))/(3.36*G*surface_density);
    }
    // Check the value of the Toomre parameter
    if(gal->comp_Q_boost[component]>0.) {
        Q = Q + gal->comp_Q_boost[component];
        v2a_r_new = pow(Q*3.36*G*surface_density,2.0)/(fabs(kappa_sqrd));
    }
    if(Q < gal->comp_Q_lim[component] && gal->comp_Q_lim[component]>0) {
        v2a_r_new = pow(gal->comp_Q_lim[component]*3.36*G*surface_density,2.0)/(fabs(kappa_sqrd));
        Q = gal->comp_Q_lim[component];
    }
    if(gal->comp_Q_fixed[component]>0) {
        v2a_r_new = pow(gal->comp_Q_fixed[component]*3.36*G*surface_density,2.0)/(fabs(kappa_sqrd));
        Q = gal->comp_Q_fixed[component];
    }

    gal->comp_Q_min[component] = (Q < gal->comp_Q_min[component] && Q > 0.0) ? Q : gal->comp_Q_min[component];

    return v2a_r_new;
}

// This function computes the Toomre criterion
double toomre(galaxy *gal, double radius, double v2a_r, int component) {

    double kappa_sqrd, h, force, dforcedr, Q, surface_density;
    double x, y, save;
    int tid;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    if(radius==0.) return 0.;

    x = radius*cos(gal->theta_cyl[gal->index[tid]]);
    y = radius*sin(gal->theta_cyl[gal->index[tid]]);
    save = gal->z[gal->index[tid]];
    gal->z[gal->index[tid]] = 0.;
    // Set the derivative step
    h = get_h_value(gal,x,y,gal->z[gal->index[tid]],0,0);
    // Calculate force and force derivative
    force = potential_deriv_wrapper_func(radius,gal);
    dforcedr = deriv_central2(gal,radius,h,potential_deriv_wrapper_func);

    kappa_sqrd = 3.0*force/radius + dforcedr;
    // We take into account the contribution of the gas to the stability of the system
    surface_density = surface_density_func(gal,radius,gal->theta_cyl[gal->index[tid]],1,component);
    if(gal->comp_type[component]==0) {
        Q = sqrt(gal->comp_cs_init[component]*fabs(kappa_sqrd))/(pi*G*surface_density);
    } else {
        Q = sqrt(v2a_r*fabs(kappa_sqrd))/(3.36*G*surface_density);
    }
    gal->comp_Q_min[component] = (Q < gal->comp_Q_min[component] && Q > 0.0) ? Q : gal->comp_Q_min[component];
    gal->z[gal->index[tid]] = save;

    return Q;
}

// This function calculates the first velocity moment in the azimuthal
// direction for the a given flat component using the axisymmetric drift approximation
double sigma2_theta_epicycle_func(galaxy *gal, double radius, double v2a_r) {

    double gamma_sqrd, kappa_sqrd, h, force, dforcedr, surface_density;
    double sigma2_theta, x, y;
    int tid;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    x = radius*cos(gal->theta_cyl[gal->index[tid]]);
    y = radius*sin(gal->theta_cyl[gal->index[tid]]);
    // Set the derivative step
    h = get_h_value(gal,x,y,gal->z[gal->index[tid]],0,0);
    // Calculate force and force derivative
    force = potential_deriv_wrapper_func(radius,gal);
    dforcedr = deriv_central2(gal,radius,h,potential_deriv_wrapper_func);

    if(force==0.) return 0.0;

    kappa_sqrd = 3.0*force/(radius)+dforcedr;
    gamma_sqrd = 4.0*force/(kappa_sqrd*radius);

    sigma2_theta = (v2a_r/gamma_sqrd);
    sigma2_theta = (sigma2_theta>0.?sigma2_theta:0.);

    return sigma2_theta;
}


// Gas theta squared velocity component
// This formulation is the same as the one used by Springel et al. 2005
double v2_theta_gas_func(galaxy *gal, double radius, double z, int component) {

    double v2_theta_gas, v_c2, pressure_force, save;
    double h, abserr, density_derivative, x, y;
    int tid;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    radius = fabs(radius);
    x = radius*cos(gal->theta_cyl[gal->index[tid]]);
    y = radius*sin(gal->theta_cyl[gal->index[tid]]);
    gal->selected_comp[tid] = component;
    // Set the derivative step
    h = get_h_value(gal,x,y,gal->z[gal->index[tid]],0,0);
    // Save z coordinate
    save = gal->z[gal->index[tid]];
    // Integrations done in the z=0 plane
    gal->z[gal->index[tid]] = 0.;
    density_derivative = deriv_central2(gal,radius,h,gas_density_wrapper_func);
    v_c2 = pow(v_c_func(gal,radius),2.0);
    //pressure_force = radius*kpc*(pow(gal->comp_cs_init[component],2.0)*density_derivative)/gas_density_wrapper_func(radius,gal);
    pressure_force = radius*(pow(gal->comp_cs_init[component],2.0)*density_derivative)/gas_density_wrapper_func(radius,gal);
    v2_theta_gas = v_c2 + pressure_force;
    if(v2_theta_gas<0) v2_theta_gas = 0.;
    // Restore z coordinate
    gal->z[gal->index[tid]] = save;

    return v2_theta_gas;
}

// A wrapper for the gas density function.
double gas_density_wrapper_func(double radius, void *params) {

    int i, tid, component;
    double theta, rho, x, y, sigma, smooth_in_factor, smooth_out_factor;
    galaxy *gal = (galaxy *) params;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    component = gal->selected_comp[tid];
    x = radius*cos(gal->theta_cyl[gal->index[tid]]);
    y = radius*sin(gal->theta_cyl[gal->index[tid]]);

    rho = 0.0;
    for(i=0; i<AllVars.MaxCompNumber; i++) {
	if(gal->comp_type[i]==0) {
            if(gal->pseudo[tid]) {
                rho += pseudo_density_gas_func(gal,fabs(radius),gal->theta_cyl[gal->index[tid]],gal->z[gal->index[tid]],0,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid],gal->comp_spherical_hydro_eq[component]);
            } else {
                rho += density_functions_pool(gal,fabs(radius),gal->theta_cyl[gal->index[tid]],gal->z[gal->index[tid]],0,gal->comp_model[gal->selected_comp[tid]],gal->selected_comp[tid]);
            }
	}
    }

    return rho;
}

// The circular velocity function. This function returns the velocity
double v_c_func(galaxy *gal, double radius) {

    double v_c, rforce, save;
    int tid;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    if(radius <= 0.0) return 0.;
    rforce = galaxy_rforce_func(gal,radius);
    v_c = sqrt(radius*rforce);

    if(rforce<0) return 0.;
    else return v_c;
}

// This function calculates the force on a test particle due to the disk
// at a point with a cylindrical radius r. This is used to calculate the velocity dispersion.
//
// This function simply performs a five point differentiation around
// x,y,z in the cylindrical radius direction with a small stepsize.
//
// This function is specially formatted for the velocity dispersion
// routines.
double galaxy_rforce_func(galaxy *gal, double radius) {

    int tid;
    double force, h, abserr;
    double x,y,r_sph;
    double sigma,transition_factor1,transition_factor2;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    x = radius*cos(gal->theta_cyl[gal->index[tid]]);
    y = radius*sin(gal->theta_cyl[gal->index[tid]]);

    h = get_h_value(gal,x,y,gal->z[gal->index[tid]],0,0);

    force = deriv_central2(gal,radius,h,galaxyr_potential_wrapper_func);

    return force;
}

// This function calculates the force on a test particle due to the disk
// at a point with an height z. This is used to calculate the velocity dispersion.
//
// This function simply performs a five point differentiation around
// x,y,z in the z direction with a small stepsize.
//
// This function is specially formatted for the velocity dispersion
// routines.
double galaxy_zforce_func(galaxy *gal, double z) {

    int tid;
    double force, h, abserr;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    h = get_h_value(gal,gal->x[gal->index[tid]],gal->y[gal->index[tid]],z,0,0);

    force = deriv_central2(gal,z,h,galaxyz_potential_wrapper_func);

    return force;
}

// This function calculates the force on a test particle due to the global potential
// at a spherical radius r_psh. This is used to calculate the velocity dispersion.
//
// This function simply performs a five point differentiation around
// x,y,z in the spherical radius direction with a small stepsize.
//
// This function is specially formatted for the velocity dispersion
// routines.
double galaxy_rsphforce_func(galaxy *gal, double r_sph) {

    int tid;
    double force, h, abserr;

#if USE_THREADS == 1
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    h = get_h_value(gal,gal->x[gal->index[tid]],gal->y[gal->index[tid]],gal->z[gal->index[tid]],0,0);

    force = deriv_central2(gal,r_sph,h,galaxyrsph_potential_wrapper_func);

    return force;
}

// This function computes the total angular momentum of a given component using the particles
double Jtot_func(galaxy *gal, int component) {
    unsigned long int i;
    double Jx, Jy, Jz, J;

    Jx = 0.;
    Jy = 0.;
    Jz = 0.;
    for(i=gal->comp_start_part[component]; i<gal->comp_start_part[component]+gal->comp_npart_pot[component]; i++) {
	if(gal->r_sph[i]<=gal->r200) {
            Jx += gal->mass[i]*(gal->y[i]*gal->vel_z[i]-gal->z[i]*gal->vel_y[i]);
            Jy += gal->mass[i]*(gal->z[i]*gal->vel_x[i]-gal->x[i]*gal->vel_z[i]);
            Jz += gal->mass[i]*(gal->x[i]*gal->vel_y[i]-gal->y[i]*gal->vel_x[i]);
	}
    }
    J = sqrt(Jx*Jx+Jy*Jy+Jz*Jz);

    return J;
}
