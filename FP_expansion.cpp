//
//
//
//  Created by Vladi Skokov on 9/26/19.
//
// t = \tau 4 \pi Nc^2 \alpha_s^2 l
// this is why my 1/\theta is different from the notes

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>

#include <iostream>

#include "FP.hpp"

const int L_max=10; // there are some limitation on L_max
//these limitation are due to GSL's factorial
//N! is only defined in GSL for N<=170
//thus  L_max < 85

const double one_over_theta0 = 5.0/(2.0*pow(M_PI,4));
const double kappa = pow(M_PI,2)/3.0 - 2 * 1.2020569031595942;

const double T_init = 1.0; //Temperature - well I do not know what would be appropriate
const double t_init = .1; //Initial t, see above its relation to tau: t = \tau 4 \pi Nc^2 \alpha_s^2 l
const double t_final = 200; //Final t

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double B(int l)
{
    double out = 2.0*(14.0*l*l + 7.0*l - 2.0)/(4.0*l-1.0)/(4.0*l+3.0);
    return out;
}

double U(int l)
{
    double out = - (2.0*l-1)*(2.0*l+1)*(2.0*l+2)/(4.0*l+3)/(4.0*l+5);
    return out;
}

double C(int l)
{
    double out = (2.0*l-1)*2.0*l*(2.0*l+2)/(4.0*l-3.0)/(4.0*l-1.0);
    return out;
}


double fact(int l)
{
    return gsl_sf_gamma(l+1.0);
}

double alpha(int l)
{
    if(l<0) return 0.0;
    if(l==0) return 1.0;
    
    //std::cout <<  l  << "\n";
    
    double out =   gsl_sf_doublefact(2*l-1)/((double) gsl_sf_fact(l));
    return out;
}

double Omega(int l, int m, int n)
{
    //std::cout << m << " " << n << " " << l  << "\n";
    
    double out  =  alpha(m-n+l) * alpha(m+n-l) * alpha(n-m+l)/alpha(m+n+l) * (4*l+1) / (2.0*(n+m+l)+1) ;
    return out;
}

double Sum_Omega(int l, const double c[])
{
    double sum = 0.0;
    for(int m=1; m<L_max; m++)
    for(int n=1; n<L_max; n++)
    {
        if(abs(m-n)<l+1) sum+= Omega(l,m,n) * c[m]*c[n];
    }
    return sum;
}

double Sum_NL(int l, const double c[])
{
    double sum = 0.0;
    for(int n=1; n<L_max; n++)
    {
        sum+=pow(c[n],2)/(4.0*n+1);
    }
    return sum * (2*l-1)*(l+1)*c[l]/3.0;
}

int RHS_f (double t, const double c[], double RHS[], void *params_ptr)
{
    double T  = c[0];
    
    for(int l=1; l<L_max; l++)
        {
                
            double B_bar = B(l) - 4.0/3.0;
            double LHS_119;
            if(l>1) {LHS_119 = 1.0/t * ( U(l) * c[l+1] + (B_bar - 2.0/15.0 * c[1]) + C(l) * c[l-1]);}
            else {LHS_119 = 1.0/t * ( U(1) * c[2] + (B_bar - 2.0/15.0 * c[1]) + C(1));}
 
            double Sum1 = Sum_Omega(l, c);
            double Sum2 = Sum_NL(l, c);
            
            double RHS_119 = - T*one_over_theta0 *(
                                                    (kappa + M_PI*M_PI*l*(2*l+1)/3.0)*c[l]
                                                   + kappa*Sum1
                                                   + kappa*Sum2
                                                   );
            
            RHS[ l ] = -LHS_119 + RHS_119;
        }
    
    RHS[0] = - T/3.0/t *(1.0+0.1*c[1]);
    
    return GSL_SUCCESS; /* GSL_SUCCESS defined in gsl/errno.h as 0 */
}

void initial(double c[])
//Set initial conditions
{
    for(int l=1; l<L_max; l++)
    {
        c[l] = 0.0;
    }
    
    c[0] = T_init; // Initial temperature
    c[L_max]=0.0; //truncation; do not modify
}



void dump (FILE * out, double c[], double t)
{
    fprintf(out,"%.5e ", t);
    for(int l=0; l<L_max; l++)
    {
       fprintf(out,"%.5e ", c[l]);
    }
    fprintf(out,"\n");
}


int main () {
    
    FILE * outdata = fopen ("output_data.csv","w");
    
    const int dimension = L_max;   /* number of differential equations */
    int status;                /* status of driver function */
    const double eps_abs = 1.e-5;  /* absolute error requested  */
    const double eps_rel = 1.e-5; /* relative error requested  */
    double c[dimension+1];            /* current solution vector */
    double t, t_next;             /* current and next independent variable */
    double tmin, tmax, delta_t;   /* range of t and step size for output */
    double h = 1.0e-10;            /* starting step size for ode solver */
    
    gsl_odeiv2_system ode_system;    /* structure with the dfunc function, etc. */
    
    /* load values into the ode_system structure */
    ode_system.function = RHS_f;       /* the right-hand-side of equation */
    ode_system.dimension = dimension;  /* number of diffeq's */
    ode_system.params = NULL;      /* parameters to pass to dfunc */
    
    tmin = t_init;            /* starting t value */
    tmax = t_final;            /* final t value */
    delta_t = 0.1;
    
    /* initial values of c */
    initial(c);
    
    gsl_odeiv2_driver * drv =
        gsl_odeiv2_driver_alloc_y_new (&ode_system, gsl_odeiv2_step_rkf45, h, eps_abs, eps_rel);
    
    t = tmin;             /* initialize t */
    /* step from tmin to tmax */
    //printf ("%.5e %.5e \n", t, f [  int_to_l(N_P/2,N_Pt/2) ] ); /* print at t=t_next */
    
    dump(outdata,c,t);
    fflush(outdata);
    
    for (t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
        {
            status = gsl_odeiv2_driver_apply (drv, &t, t_next, c);
            if (status != GSL_SUCCESS) {
                    printf("Error: status = %d \n", status);
                    break;
                }
            
            dump(outdata,c,t);
            fflush(outdata);
        } // end for
    
      
    gsl_odeiv2_driver_free (drv);
   
    fclose (outdata);
    
    return 0;
}
