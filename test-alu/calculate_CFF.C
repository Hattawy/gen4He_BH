#include "Riostream.h"
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <time.h>
#include <TMultiGraph.h>
#include <TVector3.h>
#include <TMath.h>
#include <TProfile.h>
#include "TCutG.h"
#include <TGraphErrors.h>
#include <TLatex.h>
#include "TGraph.h"
#include <TLine.h>
#include <TLorentzVector.h>
#include "TSystem.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TText.h"


double Find_BSA_err(double N_p, double N_m);
double Find_BSA(double N_p, double N_m);
void calculate_ALU(double Q2, double xB, double t,
                   double Im, double Im_err, double Re, double Re_err, double &ALU, double &ALU_err);

void calculate_CFF(double Q2, double xB, double t, 
                   double &A0, double &A1, double &A2, double &A3, 
                   double &c0_BH, double &c1_BH, double &c2_BH)
  {

  double PI = 3.1416;
  double phi = 90.0*PI/180.0;

  double MPROT = 0.93827;
  double MALPH = 3.7274;
  double EBEAM = 6.064;

  double xA = xB*MPROT/MALPH;
  double xA1= 1 - xA;
  double y = Q2/2/MALPH/xA/EBEAM;
  double e = 2*xA*MALPH/sqrt(Q2); // epsilon
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
  double FF4He = (1-pow(a*a*Q2,6))*exp(-b*b*Q2); //* 1e-15;
 
  // BH propagators
  double P1_phi = -( J + 2*K*cos(phi)) / (y * (1+e2));
  double P2_phi =  1 + t/Q2 + (1/(y*(1+e*e))) * (J + 2*K*cos(phi));


  // BH fourier coefficients
   c0_BH = ( (pow(2-y,2) + pow(y*(1+e2),2)) * (e2*Q2/t+4*(1-xA)+(4*xA+e2)*t/Q2)
                  + 2*e2*(4*(1-y)*(3+2*e2)+y*y*(2-e2*e2))
                  - 4*xA*xA*pow(2-y,2)*(2+e2)*t/Q2
                  + 8*K2*e2*Q2/t ) * pow(FF4He,2);
   c1_BH = -8*(2-y)* K * (2*xA+e*e*(1-Q2/t)) * pow(FF4He,2);
   c2_BH = 8*K2*e2*Q2*pow(FF4He,2)/t;


  // redefine the fourier hamoinic to fit for Re and Im of HA 
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


   A0 = xA* pow(1+e2, 2) * S_INT_plus_plus_1/y;
   A1 = 2*xA*xA*t*((1+e2)/Q2) * P1_phi * P2_phi * C_DVCS_0; 
   A2 = xA*pow(1+e2,2) * C_INT_plus_plus_0 /y;  
   A3 = xA* pow(1+e2,2) * C_INT_plus_plus_1/y; 
}


void calculate_ALU(double Q2, double xB, double t, 
                   double Im, double Im_err, double Re, double Re_err, double &ALU, double &ALU_err)
  {

  double PI = 3.1416;
  double phi = 90.0*PI/180.0;

  double MPROT = 0.93827;
  double MALPH = 3.7274;
  double EBEAM = 6.064;

  double xA = xB*MPROT/MALPH;
  double xA1= 1 - xA;
  double y = Q2/2/MALPH/xA/EBEAM;
  double e = 2*xA*MALPH/sqrt(Q2); // epsilon
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
  double FF4He = (1-pow(a*a*Q2,6))*exp(-b*b*Q2); //* 1e-15;
 
  // BH propagators
  double P1_phi = -( J + 2*K*cos(phi)) / (y * (1+e2));
  double P2_phi =  1 + t/Q2 + (1/(y*(1+e*e))) * (J + 2*K*cos(phi));


  // BH fourier coefficients
   double c0_BH = ( (pow(2-y,2) + pow(y*(1+e2),2)) * (e2*Q2/t+4*(1-xA)+(4*xA+e2)*t/Q2)
                  + 2*e2*(4*(1-y)*(3+2*e2)+y*y*(2-e2*e2))
                  - 4*xA*xA*pow(2-y,2)*(2+e2)*t/Q2
                  + 8*K2*e2*Q2/t ) * pow(FF4He,2);
   double c1_BH = -8*(2-y)* K * (2*xA+e*e*(1-Q2/t)) * pow(FF4He,2);
   double c2_BH = 8*K2*e2*Q2*pow(FF4He,2)/t;


  // redefine the fourier hamoinic to fit for Re and Im of HA 
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


   double A0 = xA* pow(1+e2, 2) * S_INT_plus_plus_1/y;
   double A1 = 2*xA*xA*t*((1+e2)/Q2) * (2-2*y+y*y + e2*y*y/2) * P1_phi * P2_phi * C_DVCS_0; 
   double A2 = xA*pow(1+e2,2) * C_INT_plus_plus_0 /y;  
   double A3 = xA* pow(1+e2,2) * C_INT_plus_plus_1/y; 
  
   double domm = (c0_BH + c1_BH*cos(phi) + c2_BH*cos(2*phi) + A1*(Re*Re + Im*Im) + A2*Re + A3*Re * cos(phi));

   ALU =  A0* Im* sin(phi)/ domm;
   ALU_err = ((A0* Im_err* sin(phi)) *domm - (A0* Im* sin(phi)) *(2*A1*Re*Re_err + 2*A1*Im*Im_err + A2 + A3*cos(phi))) / pow(domm, 2);

}




  double Find_BSA_err(double N_p, double N_m)
  {
      double P = 0.837;
      double BSA_err = sqrt(pow(N_p,2)*N_m + pow(N_m,2)*N_p )/ pow((N_p + N_m),2);
      double BSA = (1/P)*(N_p - N_m)/(N_p + N_m);
      double BSA_err_2 = (1/P)*sqrt((1-pow(P*BSA,2))/(N_p + N_m));
    return BSA_err_2;
   }

  double Find_BSA(double N_p, double N_m)
  {
    double BSA = 0.;
    double P = 0.837;    // plarisation of the beam 83.7% +- 0.4
    BSA = (1/P)*(N_p - N_m)/(N_p + N_m);
    return BSA;
   }

