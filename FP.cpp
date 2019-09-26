//
//  FP.cpp
//  
//
//  Created by Vladimir Skokov on 9/24/19.
//


// t = \tau 4 \pi Nc^2 \alpha_s^2 l

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_legendre.h>

#include <iostream>

#include "FP.hpp"

const double MinP = -20.0, MaxP = 20.0, MinPt = 0.0, MaxPt = 20.0; //Momentum limits
const int  N_P = 512, N_Pt = 256;
const double a_P = (MaxP - MinP)/N_P;
const double a_Pt = (MaxPt - MinPt)/N_Pt;

const double a_Pt2 = a_Pt*a_Pt;
const double a_P2 = a_P*a_P;
const double Pi = M_PI;


double IR_cut_2 = 1e-2;


struct Parameters
{
    int I;
    //int j;
    double I1;
    double I2;
    double I3;
    
};

int l_to_int_j(int I)
{
    return int(((double) I)/N_P);
}

int l_to_int_i(int I)
{
    return I-N_P*int(((double) I)/N_P);
}

int int_to_l(int i, int j)
{
    int out = i + j*N_P;
    if (out<0) std::cerr << "bounds \n";
    if (out>N_P*N_Pt-1) std::cerr << "bounds \n";
    
    return out;
}

double int_to_P(int i)
{
    return (MinP + a_P*i);
}

double int_to_Pt(int i)
{
    return (MinPt + a_Pt*i);
}

double integration(double t, const double f[], double (*underInt)(double, double, double, double))
{
  double Integ = 0.0;
  double PreFactor = 0.25*a_P*a_Pt;
  int i,j;
    
  for(i=0; i<N_P-1; i++)
      for(j=0;j<N_Pt-1; j++)
        {
            double P = int_to_P(i);
            double Pt = int_to_Pt(j);
            
            Integ += (*underInt)(P, Pt, f [  int_to_l(i,j) ], t);
            Integ += (*underInt)(P+a_P, Pt, f [  int_to_l(i+1,j) ], t);
            Integ += (*underInt)(P, Pt+a_Pt, f [  int_to_l(i,j+1) ], t);
            Integ += (*underInt)(P+a_P, Pt+a_Pt, f [  int_to_l(i+1,j+1) ], t);
        }

  Integ = Integ*PreFactor/(4.0*M_PI*M_PI);
  return(Integ);
}

double Under_cl(double P, double Pt, double f, double l)
{
    double E = sqrt(P*P + Pt*Pt + IR_cut_2);
    return ( Pt*E*(4*l+1)*gsl_sf_legendre_Pl(2*l,P/E)*f );
}

double Under_Int_J(double P, double Pt, double f, double t)
// Denoted by J in the notes
{
  double out;
  out = Pt*f*(1.0+f);
  return (out);
}


double Under_Int_K(double P, double Pt, double f, double t)
// Denoted by K in the notes MODULO factor of 2!
{
  double out;
  double E;
  
  if (Pt<a_Pt) return (0.0);
  E = sqrt(P*P + Pt*Pt + IR_cut_2);
  out = Pt*f/E;
    
  return(out);
}

double Under_E(double P, double Pt, double f, double t)
{
  double out;
  double E;
  if(Pt<a_Pt) return (0.0);

  E = sqrt(Pt*Pt+P*P + IR_cut_2);
  out = Pt*E*f;
  
 return (out);
}

double drift(int i, int j, double P, double Pt, double t,  const double f[])
{
    if(i+1>N_P-1) return 0.0;
    if(i-1<0) return 0.0;
    
    int I_p = int_to_l(i+1, j);
    int I_m = int_to_l(i-1, j);
    
    double f_p = f[I_p];
    double f_m = f[I_m];
    
    double df_dP = 0.5*(f_p-f_m)/a_P;
    
    double out = - P/t * df_dP ;
    
    return (out);
}

double collision(int i, int j, double P, double Pt, double t,  const double f[], double mueller_v3, double K, double J, double L)
{
  double out;
  double df_dP, df_dPt;
  double d2f_d2P, d2f_d2Pt, d2f_dPdPt;
  double nabla_t;
    double E, PreFactor, v3;
    
  //remove the tails
  if((i<5)||(i>(N_P-5))||(j>(N_Pt-5))) return 0.0;
  
  E = sqrt(P*P+Pt*Pt + IR_cut_2);

  if (Pt<a_Pt)
    {
      df_dPt = 0.0;
      d2f_d2Pt = 2.0 * (f[ int_to_l(i,j+1) ] - f[ int_to_l(i,j) ] ) / a_Pt2;
      nabla_t = 2.0*d2f_d2Pt;
    }
  else
    {
      df_dPt = (f[ int_to_l(i,j+1) ] - f[ int_to_l(i,j-1) ])/(2.0*a_Pt);
      d2f_d2Pt = (f[ int_to_l(i,j+1) ] - 2.0*f[ int_to_l(i,j) ] + f[ int_to_l(i,j-1) ] )/(a_Pt2);
      nabla_t = d2f_d2Pt + df_dPt/Pt;
    }
 
  df_dP = (f[ int_to_l(i+1,j) ] - f[ int_to_l(i-1,j) ])/(2.0*a_P);
  d2f_d2P = (f[ int_to_l(i+1,j) ] - 2.0*f[ int_to_l(i,j) ] + f[ int_to_l(i-1,j) ])/a_P2;
  
  out = J*(nabla_t + d2f_d2P)  +  K*4.0*(1.0+f[ int_to_l(i,j) ])*f[ int_to_l(i,j) ]/E + 2.0*K*(1.0+2.0*f[ int_to_l(i,j) ] ) *(Pt*df_dPt+P*df_dP)/E;

 return(out);
}

int RHS_f (double t, const double f[], double RHS[], void *params_ptr)
{
    /*Parameters Params = (Parameters *) params_ptr;
     get parameter(s) from params_ptr
    int I  = Params.I;
    double omega = lparams[1];
    evaluate the right-hand-side functions at t */
    
    double K = integration(t, f, &Under_Int_K);
    double J = integration(t, f, &Under_Int_J);
    
    for(int i=0; i<N_P; i++)
        for(int j=0; j<N_Pt; j++)
        {
            double P = int_to_P(i);
            double Pt = int_to_Pt(j);
            
            double temp = collision(i, j, P, Pt, t, f, 0 , K, J, 0);
            temp -= drift(i, j, P, Pt, t, f);
            
            RHS[ int_to_l(i,j) ] = temp;
            
        }
    
    return GSL_SUCCESS; /* GSL_SUCCESS defined in gsl/errno.h as 0 */
}

void initial(double f[])
{
    for(int i=0; i<N_P; i++)
        for(int j=0; j<N_Pt; j++)
        {
            double P = int_to_P(i);
            double Pt = int_to_Pt(j);
            double E = sqrt( P*P+Pt*Pt+ IR_cut_2);
            
            f[ int_to_l(i,j) ] =0.1* exp(-P*P*0.1) * 0.5* (1+tanh((-Pt+5)*10));
        }
}



void dump (FILE * out, double f[], double t)
{
    for(int i=0; i<N_P; i++)
    {
        for(int j=0; j<N_Pt; j++)
        {
            double P = int_to_P(i);
            double Pt = int_to_Pt(j);
            fprintf(out,"%.5e %.5e %.5e %.5e\n", t, P, Pt, f[ int_to_l(i,j) ]);
        }
        fprintf(out,"\n");
    }
}


int main () {
    
    FILE * outdata = fopen ("output_data.csv","w");
    FILE * outdata_init = fopen ("output_init.csv","w");
    
    const int dimension = N_P*N_Pt;   /* number of differential equations */
    int status;                /* status of driver function */
    const double eps_abs = 1.e-5;  /* absolute error requested  */
    const double eps_rel = 1.e-5; /* relative error requested  */
    double f[dimension];            /* current solution vector */
    double t, t_next;             /* current and next independent variable */
    double tmin, tmax, delta_t;   /* range of t and step size for output */
    double h = 1.0e-10;            /* starting step size for ode solver */
    
    gsl_odeiv2_system ode_system;    /* structure with the dfunc function, etc. */
    
    /* load values into the ode_system structure */
    ode_system.function = RHS_f;       /* the right-hand-side of equation */
    ode_system.dimension = dimension;  /* number of diffeq's */
    ode_system.params = NULL;      /* parameters to pass to dfunc */
    
    tmin = 1;            /* starting t value */
    tmax = 50;            /* final t value */
    delta_t = 0.1;
    /* initial values of f */
    initial(f);
    
    gsl_odeiv2_driver * drv =
        gsl_odeiv2_driver_alloc_y_new (&ode_system, gsl_odeiv2_step_rkf45, h, eps_abs, eps_rel);
    
    t = tmin;             /* initialize t */
    /* step from tmin to tmax */
    //printf ("%.5e %.5e \n", t, f [  int_to_l(N_P/2,N_Pt/2) ] ); /* print at t=t_next */
    
    dump(outdata_init,f,t);
    fflush(outdata_init);
    
    double Energy = integration(t, f, &Under_E);
    printf ("%.5e %.5e ", t, Energy);
    for(int l=0;l<21;l++)
    {
      printf ("%.5e ", integration(l, f, &Under_cl)/Energy);
    }
    printf ("\n");
    
    
    for (t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
        {
            status = gsl_odeiv2_driver_apply (drv, &t, t_next, f);
            if (status != GSL_SUCCESS) {
                    printf("Error: status = %d \n", status);
                    break;
                }
            //printf ("%.5e %.5e \n", t, f [  int_to_l(N_P/2,N_Pt/2) ] ); /* print at t=t_next */
            
            std::cerr<< "h=" << h << "\n";
            Energy = integration(t, f, &Under_E);
            printf ("%.5e %.5e ", t, Energy);
            for(int l=0;l<21;l++)
            {
              printf ("%.5e ", integration(l, f, &Under_cl)/Energy);
            }
            printf ("\n");
            
            dump(outdata,f,t);
            fflush(outdata);
        } // end for
    
    dump(outdata,f,t);
      
    gsl_odeiv2_driver_free (drv);
   
    
    fclose (outdata);
    
    return 0;
}
