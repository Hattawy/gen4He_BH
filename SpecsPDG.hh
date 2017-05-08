#ifndef SPECSPDG__hh
#define SPECSPDG__hh
#include <iostream>
#include <vector>
#include "TLatex.h"
///////////////////////////////////////////////////////////////
class SpecsPDG {
private:
    struct particle_t
    {
        int pid;
        int pidG3;
        double mass;
        double charge;
        string name;
        string latex;
        string greek;
    };
    particle_t unknown;
    particle_t electron;
    particle_t proton;
    particle_t piplus; 
    particle_t piminus;
    particle_t kplus;
    particle_t kminus;
    particle_t deuteron;
    particle_t triton;
    particle_t helion;
    particle_t alpha;
    particle_t Particle(const char* name);
    particle_t Particle(const int pid);
    particle_t ParticleG3(const int pid);
    vector <particle_t> particles;
    // RTPC drift gas specs (NeDME 4:1 ratio):
    static const double DENSITY=1.03E-3; // density (g/cm^3)
    static const double AEFF=126.787;    // effective atomic number (g/mol)
    static const double ZEFF=66;         // effective atomic charge
    static const double IEFF=99.794E-6;  // effective ionization energy (MeV)
    static const double KAPPA=0.307075;  // PDG constant (MeV/g*cm^2)
public:
    SpecsPDG();
    virtual ~SpecsPDG();
    virtual int Pid(const char* name);
    virtual int PidG3(const char* name);
    virtual double Mass(const int pid);
    virtual double MassG3(const int pid);
    virtual double Charge(const int pid);
    virtual double ChargeG3(const int pid);
    virtual double Mass(const char* name);
    virtual double Charge(const char* name);
    virtual const char* Latex(const char* name);
    virtual const char* LatexG(const char* name);
    virtual double BetheBlochMean(const double beta,const double mass,const double charge);
    virtual double BetheBlochMean(const int pid,const double mom);
    virtual double BetheBlochMeanAlpha(const float mom);
    virtual double BetheBlochMPV(const double beta,const double charge);
    virtual double BetheBlochMPV(const int pid,const double mom);
    virtual double BetheBlochMPVAlpha(const float mom);
    virtual double BetheBlochFWHM(const double beta,const double distance,const double charge);
    inline double C() { return 29.9792458; } // [cm/ns]
};
SpecsPDG::SpecsPDG()
{
    unknown.pid=0;
    unknown.pidG3=0;
    unknown.mass=0;
    unknown.charge=0;
    unknown.name="unknown";
    unknown.latex="?";

    electron.pid=11;
    electron.pidG3=0;;
    electron.mass=0.000511;
    electron.charge=-1;
    electron.name="electron";
    electron.latex="e^{-}";

    proton.pid=2212;
    proton.pidG3=2212;
    proton.mass=0.938272;
    proton.charge=1;
    proton.name="proton";
    proton.latex="p";
    proton.greek="p";
    
    piplus.pid=211;
    piplus.pidG3=0;
    piplus.mass=0.13957018;
    piplus.charge=1;
    piplus.name="piplus";
    piplus.latex="#pi^{+}";
    
    piminus.pid=-211;
    piminus.pidG3=0;
    piminus.mass=0.13957018;
    piminus.charge=-1;
    piminus.name="piminus";
    piminus.latex="#pi^{-}";
    
    kplus.pid=321;
    kplus.pidG3=0;
    kplus.mass=0.493677;
    kplus.charge=1;
    kplus.name="kplus";
    kplus.latex="K^{+}";
    
    kminus.pid=-321;
    kminus.pidG3=0;
    kminus.mass=0.493677;
    kminus.charge=-1;
    kminus.name="kminus";
    kminus.latex="K^{-}";

    deuteron.pid=1000010020;
    deuteron.pidG3=45;
    deuteron.mass=1.875613;
    deuteron.charge=1;
    deuteron.name="deuteron";
    deuteron.latex="^{2}H";
    deuteron.greek="d";

    triton.pid=1000010030;
    triton.pidG3=46;
    triton.mass=2.808921;
    triton.charge=1;
    triton.name="triton";
    triton.latex="^{3}H";
    triton.greek="T";

    helion.pid=1000020030;
    helion.pidG3=49;
    helion.mass=2.808391;
    helion.charge=2;
    helion.name="helion";
    helion.latex="^{3}He";
    helion.greek="h";

    alpha.pid=1000020040;
    alpha.pidG3=47;
    alpha.mass=3.727379;
    alpha.charge=2;
    alpha.name="alpha";
    alpha.latex="^{4}He";
    alpha.greek="#alpha";

    particles.push_back(unknown);
    particles.push_back(electron);
    particles.push_back(proton);
    particles.push_back(piplus);
    particles.push_back(piminus);
    particles.push_back(kplus);
    particles.push_back(kminus);
    particles.push_back(deuteron);
    particles.push_back(triton);
    particles.push_back(helion);
    particles.push_back(alpha);
}
SpecsPDG::~SpecsPDG(){}
///////////////////////////////////////////////////////////////
SpecsPDG::particle_t SpecsPDG::Particle(const char* name)
{
    for (unsigned int ii=0; ii<particles.size(); ii++)
        if (name==particles[ii].name) return particles[ii]; 
    return unknown;
}
SpecsPDG::particle_t SpecsPDG::Particle(const int pid)
{
    for (unsigned int ii=0; ii<particles.size(); ii++)
        if (pid==particles[ii].pid) return particles[ii]; 
    return unknown;
}
SpecsPDG::particle_t SpecsPDG::ParticleG3(const int pid)
{
    for (unsigned int ii=0; ii<particles.size(); ii++)
    {
        if (particles[ii].pidG3==0) continue; 
        if (pid==particles[ii].pidG3) return particles[ii];
    }
    for (unsigned int ii=0; ii<particles.size(); ii++)
    {
        if (pid==particles[ii].pid) return particles[ii];
    }
    return unknown;
}
int SpecsPDG::Pid(const char* name)
{ 
    return Particle(name).pid;
}
int SpecsPDG::PidG3(const char* name)
{ 
    particle_t pp=Particle(name);
    if (pp.pidG3==0) return pp.pid;
    else             return pp.pidG3;
}
double SpecsPDG::Mass(const int pid)
{
    return Particle(pid).mass;
}
double SpecsPDG::MassG3(const int pid)
{
    return ParticleG3(pid).mass;
}
double SpecsPDG::Charge(const int pid)
{
    return Particle(pid).charge;
}
double SpecsPDG::ChargeG3(const int pid)
{
    return ParticleG3(pid).charge;
}
double SpecsPDG::Mass(const char* name)
{
    return Particle(name).mass;
}
double SpecsPDG::Charge(const char* name)
{
    return Particle(name).charge;
}
const char* SpecsPDG::Latex(const char* name)
{
    return Particle(name).latex.c_str();
}
const char* SpecsPDG::LatexG(const char* name)
{
    return Particle(name).greek.c_str();
}
//////////////////////////////////////////////////////////////////////////////////////
double SpecsPDG::BetheBlochMean(const double beta,const double mass,const double charge)
{
    // input mass must be in MeV
    // from PDG's 2012 "Passage of Particles Through Matter", page 4
    static const double melec=0.511;
    if (beta<0.00001 || beta>0.99999) return 0;
    const double gamma=1/sqrt(1-beta*beta);
    const double Tmax=2*melec*pow(beta*gamma,2)/(1+2*gamma*melec/mass+pow(melec/mass,2)); // MeV
    const double coeff=KAPPA*ZEFF/AEFF*pow(charge/beta,2);
    const double bb=DENSITY*coeff*(log(2*melec*pow(beta*gamma,2)*Tmax/pow(IEFF,2))/2-beta*beta); // MeV/cm
    return bb>0 ? bb : 0;
}
double SpecsPDG::BetheBlochMean(const int pid,const double mom)
{
    // input momentum must be in MeV
    const double mass=MassG3(pid)*1000;
    const double beta = mom / sqrt( mom*mom + mass*mass );
    return BetheBlochMean(beta,mass,ChargeG3(pid));
}
double SpecsPDG::BetheBlochMeanAlpha(const float mom)
{
    // input momentum must be in MeV
    return BetheBlochMean(47,mom);
}
double SpecsPDG::BetheBlochMPV(const double beta,const double charge)
{
    // input mass must be in MeV
    // from PDG's 2012 "Passage of Particles Through Matter", page 12
    // ** I put in charge^2 coefficient because it must be missing
    if (beta<0.00001 || beta>0.99999) return 0; 
    static const double melec=0.511;
    const double gamma=1/sqrt(1-beta*beta);
    const double xi=KAPPA/2*(ZEFF/AEFF)/pow(beta,2);
    const double log1=log(2*melec*pow(beta*gamma,2)/IEFF);
    const double log2=log(xi/IEFF);
    return pow(charge,2)*xi*(log1+log2+0.2-pow(beta,2));
}
double SpecsPDG::BetheBlochMPV(const int pid,const double mom)
{
    // input momentum must be in MeV
    const double mass=MassG3(pid)*1000;
    const double beta=mom / sqrt(mom*mom + mass*mass);
    return BetheBlochMPV(beta,ChargeG3(pid));
}
double SpecsPDG::BetheBlochMPVAlpha(const float mom)
{
    // input momentum must be in MeV
    return BetheBlochMPV(47,mom);
}
double SpecsPDG::BetheBlochFWHM(const double beta,const double distance,const double charge)
{
    // input thickness must be in cm
    // from PDG's 2012 "Passage of Particles Through Matter", page 13
    // but PDG is missing charge-squared
    return 2*KAPPA*(ZEFF/AEFF)*DENSITY*distance/pow(beta,2)*pow(charge,2);
}

#endif
