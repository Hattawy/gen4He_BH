#include "Riostream.h"
#include <iostream>
#include <fstream>
#include <math.h>

struct dvcsvar_t { double q2,t,x,phi; };
struct dvcspar_t { double q2_0,alpha,b,beta,c,d; };

double BH_xsection (const dvcsvar_t &x,const dvcspar_t &p)
 {

  // calculate the cross section by taking only the BH contribution
  const double MALPH = 3.7274;
  const double MPROT = 0.93827;
  const double EBEAM = 11.0;
  
  double PI = 3.1416;
  double Q2 = x.q2 ;
  double t =  -0.1;//x.t;
  double xB = x.x;
  double phi = x.phi;

  double xA = xB*MPROT/MALPH;
  double xA1= 1 - xA;
  double y = Q2/2/MALPH/xA/EBEAM;
  double e = 2*xA*MALPH/sqrt(Q2); // epsilon
  double e2 = e*e;
  double T0 = -Q2 * (2*xA1*(1-sqrt(1+e2))+e2) / (4*xA*xA1 + e2);

  // kinematical factors
  double J = (1-y-y*e2/2) * (1+t/Q2) - (1-xA)*(2-y)*t/Q2;
  double dt = (t - T0)/Q2;
  double K2 = -1.0*dt * xA1 * (1 -y - y*y*e2/4) * (sqrt(1+e2) + ((4*xA*xA1+e2)*dt/(4*xA1)) );
  double K = sqrt(K2);
  

  // BH propagators
  double P1_phi = -1.0*(J + 2*K*cos(PI-phi)) / (y * (1+e2));
  double P2_phi =  1 + t/Q2 + (1/(y*(1+e*e))) * (J + 2*K*cos(PI-phi));

  
  // Helium form factor
  double a=0.316;
  double b=0.681;
  double FF4He = (1-pow(a*a*-t/0.197,6))*exp(-b*b*-t/0.197) * 1e-15;
 
 
  // BH fourier coefficients
  double c0_BH = ( (pow(2-y,2)+pow(y * (1+e2),2)) * (e2*Q2/t + 4*xA1 + (4*xA+e2)*t/Q2) 
                   + 2*e2*(4*(1-y)*(3+2*e2) + y*y*(2-e2*e2))
                   - 4*xA*xA*pow(2-y,2)*(2+e2)*t/Q2
                   + 8*K2*e2*Q2/t) * pow(FF4He,2);
  double c1_BH = -8*(2-y)* K * (2*xA + e2 - e2*Q2/t) * pow(FF4He,2);
  double c2_BH =  8*K2*e2*Q2*pow(FF4He,2)/t;
  
  double constant = (1/(pow(137,3)*pow(4*PI,2))) * (1/(pow(Q2,2)*pow(1+e2,2.5))) * (1/(xA*t*P1_phi*P2_phi));
 
  double calculated_cross_section = constant *(c0_BH + c1_BH*cos(PI-phi) + c2_BH*cos(2*(PI-phi)));


return calculated_cross_section/1e-36; // return in pb

}


