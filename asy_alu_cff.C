#include "Riostream.h"
#include <iostream>
#include <fstream>
#include <math.h>

double asy_alu_cff(const dvcsvar_t &xx, double phi, double Im, double Re) 
 {

  const double PI = 3.1416;
  const double MPROTON= 0.93827;
  const double MALPHA = 3.7274;
  const double EBEAM = 11.0; 


  double Q2 = xx.q2;
  double xB = xx.x;
  double t =  xx.t;

  double xA = xB*MPROTON/MALPHA;
  double xA1= 1 - xA;
  double y = Q2/2/MALPHA/xA/EBEAM;
  double e = 2*xA*MALPHA/sqrt(Q2); // epsilon
  double e2 = e*e;
  double Tmin = -Q2 * (2*xA1*(1-sqrt(1+e2))+e2) / (4*xA*xA1 + e2);
 
  // kinematical factors
  double J = (1-y-y*e2/2) * (1+t/Q2)  - (1-xA)*(2-y)*t/Q2;
  double dt = (t - Tmin)/Q2;
  double K_hat = sqrt(Tmin - t) * sqrt(xA1*sqrt(1+e2) + (Tmin - t)*(e2 + 4*xA1*xA)/(4*Q2) );
  double K = sqrt(1 - y + e2*y*y/4)* (K_hat)/(sqrt(Q2));
  double K2 = K*K; 

  // Helium form factor
  double a=0.316;
  double b=0.681;
  double FF4He = (1-pow(a*a*-t/pow(197.327/1000,2),6))*exp(-b*b*-t/pow(197.327/1000,2)); //* 1e-15;
 
  // BH propagators
  double P1_phi = -1.0*(J + 2*K*cos(PI-phi)) / (y * (1+e2));
  double P2_phi =  1 + t/Q2 + (1/(y*(1+e*e))) * (J + 2*K*cos(PI-phi));

  // BH fourier coefficients
  double c0_BH = ( (pow(2-y,2)+pow(y * (1+e2),2)) * (e2*Q2/t + 4*xA1 + (4*xA+e2)*t/Q2) 
                   + 2*e2*(4*(1-y)*(3+2*e2) + y*y*(2-e2*e2))
                   - 4*xA*xA*pow(2-y,2)*(2+e2)*t/Q2
                   + 8*K2*e2*Q2/t) * pow(FF4He,2);
  double c1_BH = -8*(2-y)* K * (2*xA + e2 - e2*Q2/t) * pow(FF4He,2);
  double c2_BH =  8*K2*e2*Q2*pow(FF4He,2)/t;


 // new parametrization  from 2010 KM paper
  // define the fourier harmonic in leptonic tensor 

  double C_DVCS_0 = 2*((2-2*y+y*y + 0.5*e2*y*y)/(1+e2)) ;

         double con_0 = -4*(2-y)*(1 + sqrt(1+e2))/pow(1+e2 ,2) ;
         double con_1 = pow(K_hat*(2-y) ,2)/(Q2*sqrt(1+e2)) ;
         double con_2 = (t/Q2)*(1-y-y*y*e2/4)*(2-xA) ;
         double con_3 = 2*xA*(t/Q2)*(2-xA + ((sqrt(1+e2) -1)/(2)) + (e2/(2*xA)) ) + e2 ;
         double con_4 = (2-xA)*(1+sqrt(1+e2));
  double C_INT_plus_plus_0 = (con_0 * ( con_1 +  con_2 * ( 1 + con_3/con_4 ) ) ) * FF4He;  

         double con_00 = (-16*K*(1 - y - e2*y*y/4))/pow(1+e2 ,5/2) ;
         double con_01 = ( 1 + (1-xA)*((sqrt(1+e2) -1) /(2*xA)) + e2/(4*xA))* (xA*t/Q2) - 3.0*e2/4.0 ;
         double con_02 = 4.0*K*(2-2*y+y*y + e2*y*y/2) * ( (1+sqrt(1+e2) -e2)/pow(1+e2 ,5/2) );
         double con_03 = (1 - (1-3.0*xA)*t/Q2 + (1-sqrt(1+e2)+3*e2)/(1+sqrt(1+e2) -e2)*(xA*t/Q2)); 
  double C_INT_plus_plus_1 =  (con_00 * con_01 - con_02 * con_03) * FF4He;  

         double  con_30 = 8*K*(2-y)*y /(1+e2);
         double  con_31 = ((1 - xA + 0.5*(sqrt(1+e2)-1))/(1+e2))*dt;
  double S_INT_plus_plus_1 = con_30 * ( 1 + con_31 ) * FF4He;  


  double BH_Con = c0_BH + c1_BH*cos(phi) + c2_BH*cos(2*phi);

  double A0 = xA* pow(1+e2, 2) * S_INT_plus_plus_1/y;
  double A1 = 2*xA*xA*t*((1+e2)/Q2) * P1_phi * P2_phi * C_DVCS_0; 
  double A2 = xA*pow(1+e2,2) * C_INT_plus_plus_0 /y;  
  double A3 = xA* pow(1+e2,2) * C_INT_plus_plus_1/y; 

  double domm = (BH_Con + A1*(Re*Re + Im*Im) + A2*Re + A3*Re * cos(phi));
  double ALU =  A0* Im* sin(phi)/domm;

  return ALU;
}



