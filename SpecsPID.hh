class SpecsPID {
private:
    int runnumber;
    virtual float ICVTSIGMA(const float e);
    virtual float ICVTMEAN(const float e);
public:
    SpecsPID();
    virtual ~SpecsPID();
//    bool ElectronID(const int iEVNT,const float mean,const float width,const float nphot);
    bool ICVTCUT(const float e,const float dt);
    bool ICVTCUT(const float t0,const float e,const float t,const float z,const float x,const float y);
};
SpecsPID::SpecsPID()
{
    runnumber=0;
}
SpecsPID::~SpecsPID(){}


/*
bool SpecsPID::ElectronID(MyH10 *h10,const int iEVNT,const float mean,const float width,const float nphot)
{
    const int iEC=h10->ec[iEVNT]-1;

    if (iEC<0 || iEC>=h10->ec_part) return 0;

    if (h10->stat[iEVNT] <= 0 || h10->q[iEVNT] >= 0) return 0;

    if (h10->ec_ei[iEC]<0.08) return 0;
    
    const double ece=std::max(h10->etot[iEC],h10->ec_ei[iEC]+h10->ec_eo[iEC]);
    
    if (fabs(ece/h10->p[iEVNT]-mean) > width) return 0;

    if (nphot>0 && h10->nphe[iEC]<nphot) return 0;
    
    return 1;
}
*/
float SpecsPID::ICVTSIGMA(const float e)
{
    static const double par[7]={
         9.37046e+00,
        -1.94524e+01,
         6.27129e-01,
        -1.16642e-01,
        -6.01624e-03,
         1.09265e-02,
        -1.47874e-03
    };
    return par[0]*exp(par[1]*e)+
        par[2]+par[3]*e+par[4]*e*e+par[5]*e*e*e+par[6]*e*e*e*e;
}
float SpecsPID::ICVTMEAN(const float e)
{
    static const double par[7]={
         1.54689e+01,
        -1.05420e+01,
         8.02862e-01,
        -6.46833e-01,
         1.71032e-02,
         1.28834e-02,
        -1.32820e-03
    };
    return par[0]*exp(par[1]*e)+
        par[2]+par[3]*e+par[4]*e*e+par[5]*e*e*e+par[6]*e*e*e*e;
}
bool SpecsPID::ICVTCUT(const float e,const float dt)
{
    // e is gamma energy in GeV (etc)
    // dt is "IC minus CLAS" : (tc - sqrt(xc*xc+yc*yc+vz*vz)/30 - tr_time)
    // returns true if within 3-sigma time window, else false
    const float mean = ICVTMEAN(e);
    float sigma;
    if      (e<0.2) sigma = ICVTSIGMA(0.2);
    else if (e>5.0) sigma = ICVTSIGMA(5);
    else            sigma = ICVTSIGMA(e);
    return fabs(dt-mean)<3*sigma;
}
bool SpecsPID::ICVTCUT(const float t0,const float e,const float t,const float z,const float x,const float y)
{
    const float r=sqrt(x*x+y*y+z*z);
    const float dt=t0-(t-r/30);
    return ICVTCUT(e,dt);
}
bool ICVTCUT(const float t0,const float e,const float t,const float z,const float x,const float y)
{
    static SpecsPID s;
    const float r=sqrt(x*x+y*y+z*z);
    const float dt=t0-(t-r/30);
    return s.ICVTCUT(e,dt);
}

