#ifndef __FASTMC__
#define __FASTMC__
#include <iostream>
#include "SpecsFID.hh"
#include "TRandom3.h"
//TRandom3 RANDOM;
SpecsFID FID;
SpecsGEO GEO;

#define R2D 57.29578
//const double R2D=TMath::RadToDeg();

// smearing factors from e1dvcs ///////////
#define FDC_MOM 3.4
#define FDC_THE 2.5
#define FDC_PHI 4.0
#define FIC_POS 1.2
//const double FDC_MOM=3.4;
//const double FDC_THE=2.5;
//const double FDC_PHI=4.0;
//const double FIC_POS=1.2;
///////////////////////////////////////////

// additional IC smearing based upon eg6's pi0 width (~10MeV, 30% larger than eg1dvcs),
// and FX's note that the eqn below translates to 7/8 MeV pi0 width for 1GeV gammas
#define FIC_ENE 1.33
//const double FIC_ENE=1.33;

// calorimeter thresholds:
#define IC_THRESH 0.15
#define EC_THRESH 0.15
//const double IC_THRESH=0.15;
//const double EC_THRESH=0.15;

// electron vertex resolution
#define ELEVZRES 0.3
//const double ELEVZRES=0.3; // measured on eg6 downstream beam window

// RTPC vertex resolution: 
#define RTPCVZRES 0.46
//const double RTPCVZRES=0.46; // assuming CLAS-RTPC width is quadratic sum of individual resolutions

// RTPC resolutions:
#define FTPC_MOM 0.1
#define FTPC_THE 4.0
#define FTPC_PHI 2.0
//const double FTPC_MOM=0.1; // relative // this is unmeasured!!!
//const double FTPC_THE=3.0; // deg (measured from elastics)
//const double FTPC_PHI=2.0; // deg (measured from elastics)

double SmearDCmom(const double p,const double theta,const double beta,const double curr)
{
    const double sigma = FDC_MOM * p * 
                         pow(theta*R2D/35,0.7) * 
                         sqrt(pow(0.0033*p,2)+pow(0.0018/beta,2)) *
                         3375/curr;
    return RANDOM.Gaus(p,sigma);
}
double SmearDCthe(const double p,const double theta,const double beta)
{
    const double sigma = FDC_THE * sqrt(pow(0.55,2)+pow(1.39/p/beta,2)); // mrad
    return RANDOM.Gaus(theta,sigma/1000);
}
double SmearDCphi(const double p,const double phi,const double beta)
{
    const double sigma = FDC_PHI*sqrt(pow(3.73,2)+pow(3.14/p/beta,2)); // mrad
    return RANDOM.Gaus(phi,sigma/1000);
}
double SmearECene(const double e)
{
    const double sigma = 0.116*sqrt(e);
    return RANDOM.Gaus(e,sigma);
}
double SmearECang(const double ang)
{
    const double sigma = 0.004;
    return RANDOM.Gaus(ang,sigma);
}
double SmearICene(const double e)
{
    const double sigma = e*sqrt(pow(0.024,2)+pow(0.033/sqrt(e),2)+pow(0.019/e,2));
    return RANDOM.Gaus(e,sigma*FIC_ENE);
}
double SmearICthe(const double e,const double theta)
{
    const double sigma = FIC_POS*sqrt(pow(0.003/sqrt(e),2)+pow(0.013*theta,2));
    return RANDOM.Gaus(theta,sigma);
}
double SmearICphi(const double e,const double phi)
{
    const double sigma = FIC_POS*0.003/sqrt(e);
    return RANDOM.Gaus(phi,sigma);
}

TLorentzVector SmearDC(const TLorentzVector v,const double curr)
{
    const double m=v.M2();
    const double p=SmearDCmom(v.P(),v.Theta(),v.Beta(),curr);
    const double the=SmearDCthe(v.P(),v.Theta(),v.Beta());
    const double phi=SmearDCphi(v.P(),v.Phi(),v.Beta());
    const TVector3 dir(sin(the)*cos(phi),sin(the)*sin(phi),cos(the));
    return TLorentzVector(dir*p,sqrt(p*p+m));
}
TLorentzVector SmearEC(const TLorentzVector v)//,const double dvz=0)
{
    const double m2=v.M2();
    const double p=SmearECene(v.E());
    const double the=SmearECang(v.Theta());
    const double phi=SmearECang(v.Phi());
    const TVector3 dir(sin(the)*cos(phi),sin(the)*sin(phi),cos(the));
    return TLorentzVector(dir*p,sqrt(p*p+m2));
}
TLorentzVector SmearIC(const TLorentzVector v,const double vz)
{
    const double ric=tan(v.Theta())*vz;
    const double dthedz=-ric/(vz*vz+ric*ric);
    const double sigthe=dthedz*RANDOM.Gaus(0,ELEVZRES);

    const double m2=v.M2();
    const double p=SmearICene(v.E());
    const double the=SmearICthe(v.E(),v.Theta())+sigthe;
    const double phi=SmearICphi(v.E(),v.Phi());
    const TVector3 dir(sin(the)*cos(phi),sin(the)*sin(phi),cos(the));
    return TLorentzVector(dir*p,sqrt(p*p+m2));
}
double SmearRTPCthe(const double theta)
{
    // based upon z-vertex resolution:
    // CLAS @ beam window ~3mm
    // CLAS-RTPC ~5.5 mm
    // ---> RTPC ~4.6 mm
    // +0.8deg to put resolution at 3deg for elastics @ 78deg
    
    static const double R=4.5;
    static const double sigmaZ=0.46; // cm
    const double dthedz=sin(2*theta)/R;
    const double sigmaThe=fabs(sigmaZ*dthedz)+0.8/R2D;
    return RANDOM.Gaus(theta,sigmaThe);
}
TLorentzVector SmearRTPC(const TLorentzVector v)
{
    const double m2=v.M2();
    const double sigma = v.P()*FTPC_MOM;
    const double p=RANDOM.Gaus(v.P(),sigma);
    const double the=RANDOM.Gaus(v.Theta(),FTPC_THE/R2D);
    //const double the=SmearRTPCthe(v.Theta());
    const double phi=RANDOM.Gaus(v.Phi(),FTPC_PHI/R2D);
    const TVector3 dir(sin(the)*cos(phi),sin(the)*sin(phi),cos(the));
    return TLorentzVector(dir*p,sqrt(p*p+m2));
}

TVector3 GammaICpos(const TLorentzVector v,const double vz)
{
    const TVector3 vert(0,0,vz);
    return TVector3(vert + v.Vect().Unit()*(-vz/cos(v.Theta())));
}
bool AcceptIC(const TLorentzVector v,const double vz)
{
    if (v.E() < IC_THRESH) return 0;
    const TVector3 icpos=GammaICpos(v,vz);
//    const TVector3 vert(0,0,vz);
//   const TVector3 icpos(vert + v.Vect().Unit()*(-vz/cos(v.Theta())));
    if (!FID.FxIC(icpos.X(),icpos.Y())) return 0;
    return 1;
}
bool AcceptCLAS(const TLorentzVector v,const double vz,const double charge)
{
    // BEWARE!!!!!!!!!!



    if (charge==0)
    {
        // eyeballed this from eg6 EC photons:
        if (v.Theta()*R2D < 15 + 5./20.*(vz+75)) return 0;

        if (v.E() < EC_THRESH) return 0;
    }
    else
    {
        // from FX's e1dvcs analysis note (2099A):
        // NEEDS to BE UPDATED FOR eg6's 1900A?
        // NEED to implement z-dependence
    
        double thcut,phi0,sigmadphi;
        if (!FID.SOL(vz,v.CosTheta())) return 0;
        //    if (!FID.RTPCendplate(vz,v.CosTheta())) return 0;
        
        if (charge==-1) // FOR ELECTRONS
        {
            if (v.P()<0.5) return 0;
            thcut=16.8+0.93/(v.P()-0.52);
            phi0=charge*(1.012+6.632/v.P());
            sigmadphi=27*pow(sin(v.Theta()-thcut/R2D),0.195);
        }
        else if (charge==1) // FOR PROTONS
        {
            if (v.P()<0.3) return 0;
            thcut=16.8;
            phi0=-(-1.4-10.5/v.P()+sin(v.Theta())/0.66);
            sigmadphi=27;
        }
        
        if (v.Theta()*R2D < thcut) return 0;
        const int phisec=60*GEO.GetSector(v);
        const float phi=GEO.GetPhiCLAS(v)*R2D;
        if (fabs(phisec-phi0-phi) > sigmadphi) return 0;
    }



    const TVector3 icpos=FID.Bosted2IC(v,vz,charge);

    // Remove if hits IC:
    if (charge!=0)
    {
        if (FID.BostedIC(icpos.X(),icpos.Y(),2,0)) return 0;
    }
   

    if (!FID.FxInSector(icpos.X(),icpos.Y())) return 0;


    return 1;
}
bool AcceptRTPC(const double vz,const double theta,const double phi)
{
    static const double rinner=3;
    static const double router=6;
    const TVector3 vert(0,0,vz+64);

    const TVector3 dir(cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta));

    // make sure track intersects cathode and gem inside RTPC:
    const TVector3 r1=vert+dir.Unit()*(rinner/fabs(sin(dir.Theta())));
    const TVector3 r2=vert+dir.Unit()*(router/fabs(sin(dir.Theta())));
    if (fabs(r1.Z())>10) return 0;
    if (fabs(r2.Z())>10) return 0;

    // kill top/bottom support regions:
    const double phi2=GEO.GetPhi(phi,0)*R2D;
    if (fabs(phi2- 90)<30) return 0;
    if (fabs(phi2-270)<30) return 0;

    // kill if track goes through target holder at upstream end:
    if (FID.RTPCnose(vz,cos(theta))) return 0;

    return 1;
}
bool AcceptRTPC(const TLorentzVector v,const double vz)
{
    // input vz is in cm in the CLAS coordinate system

    // fiducial cut:
    if (!AcceptRTPC(vz,v.Theta(),v.Phi())) return 0;
    
    // p/theta thresholds estimated from eloss correction:
    double the = v.Theta()*R2D;
    if (the<90) the=180-the;
    if (fabs(v.M()-3.7274)<1e-2) // 4He
    {
        if (v.P()*1000 < 1./sin(the/R2D)*62+200) return 0;
    }
    else if (fabs(v.M()-2.8084)<1e-2) // 3He
    {
        if (v.P()*1000 < 1./sin(the/R2D)*60+155) return 0;
    }
    else if (fabs(v.M()-0.93827)<1e-2) // proton
    {
        if (v.P()*1000 < 1./sin(the/R2D+8/R2D)*8+60) return 0;
    }
    else std::cerr<<"AcceptRTPC not ready for this recoil type:  mass="<<v.M()<<std::endl;
    return 1;
}
void PrintTLV(const TLorentzVector v)
{
    printf("m=%12.3f p=%12.3f the=%12.3f phi=%12.3f\n",
    v.M(),v.P(),v.Theta()*TMath::RadToDeg(),v.Phi()*TMath::RadToDeg());
}
#endif
