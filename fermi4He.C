#include "TLorentzVector.h"
#ifndef __FASTMC__
#include "TRandom3.h"
TRandom3 RANDOM;
#endif
double pfermi4He(const double k)
{
    // Nucleon momentum distribution in 4He.
    // From C.Ciofi degli Atti and S.Simula, PRC 53 (1996) 1689.
    static const double max=0.0907;
    static const double norm=1.0/(4.0*3.1415926535);
    static const double a0=4.33,b0=1.54,c0=0.419,d0=5.49,e0=4.9,f0=0.0;
    static const double a1=0.665,b1=2.15,b2=0.0,c1=0.0244,d1=0.22;
    const double k2=k*k;
    const double n0 = a0*exp(-b0*k2)/pow(1+c0*k2,2) + 
                      d0*exp(-e0*k2)/pow(1+f0*k2,2) ;
    const double tmp1 = a1*exp(-b1*k2)/pow(1+b2*k2,2);
    const double hard = c1*exp(-d1*k2);
    const double n1 = tmp1 + hard;
    return norm*(n0+n1)*k2/max;
}
double pfermi4He()
{
    static const double pmax=300.0; // truncating high-momentum tail
    double pp,rr;
    while (1)
    {
        pp=RANDOM.Rndm()*pmax;
        rr=RANDOM.Rndm();
        if (pfermi4He(pp/197.32696) > rr)
            return pp/1000; // convert to GeV
    }
}
TLorentzVector vfermi4He()
{
    const double pp=pfermi4He();
    const double ee=sqrt(pp*pp+pow(0.93827,2));
    const double the=acos(RANDOM.Uniform(-1,1));
    const double phi=RANDOM.Uniform(0,2*TMath::Pi());
    const TVector3 dir(sin(the)*cos(phi),sin(the)*sin(phi),cos(the));
    return TLorentzVector(dir.Unit()*pp,ee);
}
double pfermi4He(double *x,double *p) { return p[0]*pfermi4He(x[0]/197.32696*1000); }

