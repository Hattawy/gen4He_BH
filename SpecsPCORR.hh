#ifndef __SPECSPCORR__
#define __SPECSPCORR__
#include "MyRootAna.C"
#include "TSystem.h"
#include "TString.h"
class SpecsPCORR {
private:
public:
    virtual double epcorr_phi(const double p,const double cx,const double cy);
    virtual double epcorr_z1(const double p,const double vz);
    virtual double epcorr_the1(const double p,const double cz);
    virtual double epcorr_phi_the(const double p,const double cx,const double cy,const double cz);
    virtual double ElectronPCORR(const double p,const double cx,const double cy,const double cz,const double vz);
    virtual double EG6ICecorr(const double ee,const double rf150);
    virtual double EG6ICecorr3(const double ee,const double rf150);
    virtual double EG6ICeAbsCorrRadial(const double ee,const double radius);
    virtual double EG6ICeRelCorrRadial(const double ee,const double radius);
    virtual double EG6ECsampcorr(const int sector,const double rf150);
    virtual double E16ECphocorr(const double etot,const int sector,const double rf150);
    virtual double SimpleRtpcVertexCorr(const double theta,const double zz);
    virtual double SimpleRtpcThetaCorr(const double theta);
//    virtual double EG6ICmpi0Radius(const double r);
    virtual double EG6ICmpi0RadiusC(const double r);
};

double SpecsPCORR::EG6ICeAbsCorrRadial(const double ee,const double radius)
{   // radius in cm
    // energy in GeV
    // *** Measured via Missing Energy in Coherent Channel ***
    //     Input is mohammad's corrected photon energy
    const double eX = 1.113 * exp( -0.3432 * (radius-3.111) );
    return ee+eX;
}
double SpecsPCORR::EG6ICeRelCorrRadial(const double ee,const double radius)
{   // radius in cm
    // energy in GeV
    // *** Measured via Missing Energy in Coherent Channel ***
    //     Input is mohammad's corrected photon energy
    const double eXovere = 0.3044 * exp( -0.3905 * (radius-2.901) );
    return ee*(eXovere+1);
}
double SpecsPCORR::EG6ICmpi0RadiusC(const double r)
{
    // input is cluster radius in cm
    // output is scale factor to multiply by photon energy 
    // should be used with Mohammad's corrected photon energy
    
    static const int nn1=9;
    static const double pp1[nn1]=
    {
         6.90509265413806617e+01,
        -6.19324582243512154e+01,
         2.41565771428874072e+01,
        -5.27375202981722335e+00,
         7.05591853764797605e-01,
        -5.93279091766953787e-02,
         3.06575949389111067e-03,
        -8.91258233223929261e-05,
        1.11719930566835684e-06
    };
    static const int nn2=4;
    static const double pp2[nn2]=
    {
        3.36555662260284283e+00, 
       -1.33525171666376541e+00,
        2.56227164257084528e-01,
       -1.63874655493372913e-02
    };

    double ff=0;
    const int     nn = r<6 ? nn2 : nn1;
    const double *pp = r<6 ? pp2 : pp1;
    for (int ii=0; ii<nn; ii++) ff += pp[ii]*pow(r,ii);
    return ff;
}
/*
double SpecsPCORR::EG6ICmpi0Radius(const double r)
{
    // derived from:
    // 1)  Event Selection Requirement: fabs(r1-r2)<1
    // 2)  Assumption:  Ec = alpha(r) * E 
    // 3)  Then:  alpha(r) = 0.135/mgg

    static const double p[]=
    {
         1.07785e-01,
         1.43396e+02,
         1.16635e-01,
        -5.95037e+01,
         2.52692e+01,
        -9.81404e+00,
         1.75159e+00,
        -1.60679e-01,
         7.41656e-03,
        -1.36718e-04
    };
  const double mpi0=(p[0]-p[1]*exp(-(r-p[3])*p[2])) *
         (p[4]+p[5]*r+p[6]*r*r+p[7]*r*r*r+p[8]*r*r*r*r+p[9]*r*r*r*r*r);

  return 0.135/mpi0;
}
*/
double SpecsPCORR::epcorr_phi(const double p,const double cx,const double cy)
{
    static TF1 *fallsectors=NULL;
    static TF1 *fsectors[6]={NULL};
    if (!fsectors[0] || !fallsectors)
    {
        TFile *f=GetMyFile("epcorr_phi.root");
        if (!f) return 9999;
        for (int ii=0; ii<6; ii++)
            fsectors[ii]=(TF1*)f->Get(Form("fp6_0_%d",ii))
                ->Clone(Form("fsectors_%d",ii));
        fallsectors=(TF1*)f->Get("fitgraph")->Clone("fallsectors");
        f->Close();
    }
    const int sector=GetSector(cx,cy);
    const double phid=GetPhi(cx,cy,-TMath::Pi()/6)*180/TMath::Pi();

    // return +/- to tell wether it's in a valid phi range:
    if (p < 0)
    {
        if (phid<fsectors[sector]->GetXmin()) return -1;
        if (phid>fsectors[sector]->GetXmax()) return -1;
        else                                  return  1;
    }

    // return the corrected momentum:
    else  return p*fallsectors->Eval(phid); 

/*
    static const int polyorder=4;
    static const double pars[6*(polyorder+1)]={
        1.00052,  -0.000277844,  8.16856e-06, -3.26174e-07,  2.28242e-09,
        0.997162, -0.000142453,  3.39238e-05, -2.73832e-06, -1.22078e-07,
        0.985489,  0.00590335,  -0.000844199,  4.75279e-05, -8.98404e-07,
        0.999446,  1.30988e-05, -4.53816e-07, -3.33824e-07,  1.75588e-08,
        1.00259,  -0.000271954,  3.26552e-05,  3.36167e-07, -5.44711e-08,
        1.08412,  -0.0198513,    0.00167425,  -6.03819e-05,  7.99104e-07 };
    const int sector=GetSector(cx,cy);
    const double phi=GetPhi(cx,cy,-TMath::Pi()/6);
    const double localphi=GetLocalPhi(phi);
    const int ipar0=sector*(polyorder+1);
    double f=0;
    for (int ii=0; ii<polyorder+1; ii++)
        f += pars[ipar0+ii]*pow(localphi,ii);
    return p*f;
*/
}








double SpecsPCORR::epcorr_z1(const double p,const double vz)
{
    static TF1 *fit_z=NULL;
    if (!fit_z)
    {
        TFile *f=GetMyFile("epcorr_z1.root");
        if (!f) return 9999;
        fit_z=(TF1*)f->Get("fpoly_p14_0")->Clone("fit_z");
        f->Close();
    }
    return p*fit_z->Eval(10*vz+640); 
}
double SpecsPCORR::epcorr_the1(const double p,const double cz)
{
    static TF1 *fit_the=NULL;
    if (!fit_the)
    {
        TFile *f=GetMyFile("epcorr_the1.root");
        if (!f) return 9999;
        fit_the=(TF1*)f->Get("fpolycossin_p4_0")->Clone("fit_the");
        f->Close();
    }
    return p*fit_the->Eval(acos(cz)*180/TMath::Pi()); 
}

double SpecsPCORR::epcorr_phi_the(const double p,const double cx,const double cy,const double cz)
{
    return epcorr_the1( epcorr_phi(p,cx,cy), cz);
}
double SpecsPCORR::ElectronPCORR(const double p,const double cx,const double cy,const double cz,const double vz)
{
    return epcorr_z1(epcorr_the1( epcorr_phi(p,cx,cy), cz), vz);
}




double SpecsPCORR::EG6ICecorr(const double etc,const double rf150)
{
//    rf150  == Run# + File# / 150
//
    // Mohammad's 1st IC energy correction
    // run-dependent and energy-dependent

//    Returns the sampling fraction parameterized by piecewise functions of the form:
//    E0 + A0 * ( exp(-ALPHA1*(x-X0)) + exp(-ALPHA2*(x-tX0)) )
//    where x = Run# + File#/150
//    Negative return value is an error.
//    Reads parameters from EG6ICecorr.dat
//
    static bool initialized=0;
    static const int npar=5;

    // table columns:
    static const int ncols=1+2+npar;
    static vector <double> tXLO,tXHI;  // <--- Run# + File#/150
    static vector <double> tE0,tX0,tA0,tALPH1,tALPH2;

    // read the data table:
    if (!initialized)
    {
        const char* stub="fit_M_Pi0_Run_1.dat";
        const char* cparm=gSystem->Getenv("CLAS_PARMS");
        const TString filename=cparm?Form("%s/%s",cparm,stub):stub;

        FILE *fin=fopen(filename.Data(),"r");
        if (!fin)
        {
            fprintf(stderr,"EG6ICecorr:  Missing Input File:  %s\n",filename.Data());
            return -9999;
        }

        int nlines=0;
        char buf[1024];
        double tmp[ncols];

        while ((fgets(buf,1024,fin)) != NULL)
        {
            // skip header line:
            if (nlines++==0) continue;

            for (int ii=0; ii<ncols; ii++)
            {
                if (ii==0) tmp[ii] = atof((char*)strtok(buf," "));
                else       tmp[ii] = atof((char*)strtok(NULL," "));
            
                if (isnan(tmp[ii]))
                {
                    fprintf(stderr,"EG6ICecorr:  Error Reading File:  %s\n",filename.Data());
                    return -9999;
                }
            }

            int jj=1;
            tXLO  .push_back(tmp[jj++]);
            tXHI  .push_back(tmp[jj++]+0.999);
            tE0   .push_back(tmp[jj++]);
            tX0   .push_back(tmp[jj++]);
            tA0   .push_back(tmp[jj++]);
            tALPH1.push_back(tmp[jj++]);
            tALPH2.push_back(tmp[jj++]);
        }

        initialized=1;
        fclose(fin);
    }


    // check for invalid input:
    if (rf150 < tXLO[0] || rf150 > tXHI[tXHI.size()-1])
    {
        fprintf(stderr,"EG6ICecorr:  Invalid rf150:  %f\n",rf150);
        return -9999;
    }

    // binary search:
    int imax=tXLO.size()-1;
    int imin=0;
    int jj=-1;

    while (imax >= imin)
    {
        jj = (imax+imin)/2;
        if      (rf150 > tXHI[jj]) imin=jj+1;
        else if (rf150 < tXLO[jj]) imax=jj-1;
        else break;
    }

    if (rf150<tXLO[jj] || rf150>tXHI[jj])
    {
        fprintf(stderr,"EG6ICecorr:  Missing rf150:  %f\n",rf150);
        return -9999;
    }

    //printf("%d %f %f %f %f %f %f %f\n",
    //        jj+1,tXLO[jj],tXHI[jj],
    //        tE0[jj],tX0[jj],tA0[jj],
    //        tALPH1[jj],tALPH2[jj]);


    // the parameterization:
    const double corr=
        tE0[jj] +
        tA0[jj] * exp(-tALPH1[jj]*(rf150-tX0[jj])) -
        tA0[jj] * exp(-tALPH2[jj]*(rf150-tX0[jj])) ;

/*
    static int ncalls=0;
    if (ncalls++==0)
    {
        TF1 *f[100];
        TH2F *h=new TH2F("h","",1000,61510,61930,1000,-10,10);
        h->DrawClone();
        for (unsigned int ii=0; ii<tX0.size(); ii++)
        {
            f[ii]=new TF1(Form("f__%d",ii),EG6ICecorr2,tXLO[ii],tXHI[ii],0);
            f[ii]->SetLineColor(kRed);
            if (ii==0) f[ii]->Draw("RSAME");
            else       f[ii]->Draw("RSAME");
        }
    }
*/

    // 1st: energy-dependent correction:
    static const double alpha=-0.001278;
    static const double beta=0.137054;
    const double etc2 = etc * 0.135 / (2*alpha*etc+beta);
    
    // 2nd: time-dependent correction:
    const double etc3 = etc2*0.135/corr;

    return etc3;
}

double EG6ICecorr39(double* x,double *)
{
//    (void*) p;
    static SpecsPCORR s;
    return s.EG6ICecorr3(1.5,x[0]);
}
double SpecsPCORR::EG6ICecorr3(const double etc,const double rf150)
{
// mohammad's correction #3

    static bool initialized=0;
    static const int npar=5;

    // table columns:
    static const int ncols=1+2+npar;
    static vector <double> aXLO,aXHI;  // <--- Run# + File#/150
    static vector <double> bXLO,bXHI;  // <--- Run# + File#/150
    static vector <double> a0,a1,a2,a3,a4;
    static vector <double> b0,b1,b2,b3,b4;

    // read the data table:
    if (!initialized)
    {
        const char* cparm=gSystem->Getenv("CLAS_PARMS");
        const char* stubA="alpha_Run_final.dat";
        const char* stubB="beta_Run_final.dat";
        const TString filenameA=cparm?Form("%s/%s",cparm,stubA):stubA;
        const TString filenameB=cparm?Form("%s/%s",cparm,stubB):stubB;

        int nlines=0;
        char buf[1024];
        double tmp[ncols];

        FILE *finA=fopen(filenameA.Data(),"r");
        if (!finA)
        {
            fprintf(stderr,"EG6ICecorr:  Missing Input File:  %s\n",filenameA.Data());
            return -9999;
        }

        while ((fgets(buf,1024,finA)) != NULL)
        {
            // skip header line:
            if (nlines++==0) continue;

            for (int ii=0; ii<ncols; ii++)
            {
                if (ii==0) tmp[ii] = atof((char*)strtok(buf," "));
                else       tmp[ii] = atof((char*)strtok(NULL," "));
            
                if (isnan(tmp[ii]))
                {
                    fprintf(stderr,"EG6ICecorr:  Error Reading File:  %s\n",filenameA.Data());
                    return -9999;
                }
            }

            int jj=1;
            aXLO  .push_back(tmp[jj++]);
            aXHI  .push_back(tmp[jj++]+0.999);
            a0    .push_back(tmp[jj++]);
            a1    .push_back(tmp[jj++]);
            a2    .push_back(tmp[jj++]);
            a3    .push_back(tmp[jj++]);
            a4    .push_back(tmp[jj++]);
        }

        fclose(finA);
       

        FILE *finB=fopen(filenameB.Data(),"r");
        if (!finB)
        {
            fprintf(stderr,"EG6ICecorr:  Missing Input File:  %s\n",filenameB.Data());
            return -9999;
        }

        while ((fgets(buf,1024,finB)) != NULL)
        {
            // skip header line:
            if (nlines++==0) continue;

            for (int ii=0; ii<ncols; ii++)
            {
                if (ii==0) tmp[ii] = atof((char*)strtok(buf," "));
                else       tmp[ii] = atof((char*)strtok(NULL," "));
            
                if (isnan(tmp[ii]))
                {
                    fprintf(stderr,"EG6ICecorr:  Error Reading File:  %s\n",filenameB.Data());
                    return -9999;
                }
            }

            int jj=1;
            bXLO  .push_back(tmp[jj++]);
            bXHI  .push_back(tmp[jj++]+0.999);
            b0    .push_back(tmp[jj++]);
            b1    .push_back(tmp[jj++]);
            b2    .push_back(tmp[jj++]);
            b3    .push_back(tmp[jj++]);
            b4    .push_back(tmp[jj++]);
        }

        fclose(finB);
        
        
        initialized=1;
    }


    // check for invalid input:
    if (rf150 < aXLO[0] || rf150 > aXHI[aXHI.size()-1] ||
        rf150 < bXLO[0] || rf150 > bXHI[bXHI.size()-1] )

    {
        fprintf(stderr,"EG6ICecorr:  Invalid rf150:  %f\n",rf150);
        return -9999;
    }

    // find run range for alpha and beta
    int aa=-1;
    int bb=-1;
    
    // binary search:
    int imin=0;
    int imax=aXLO.size()-1;
    while (imax >= imin)
    {
        aa = (imax+imin)/2;
        if      (rf150 > aXHI[aa]) imin=aa+1;
        else if (rf150 < aXLO[aa]) imax=aa-1;
        else break;
    }

    imin=0;
    imax=bXLO.size()-1;
    while (imax >= imin)
    {
        bb = (imax+imin)/2;
        if      (rf150 > bXHI[bb]) imin=bb+1;
        else if (rf150 < bXLO[bb]) imax=bb-1;
        else break;
    }

    if (rf150<aXLO[aa] || rf150>aXHI[aa] ||
        rf150<bXLO[bb] || rf150>bXHI[bb] )

    {
        fprintf(stderr,"EG6ICecorr:  Missing rf150:  %f\n",rf150);
        return -9999;
    }

    const double alpha = a0[aa] + a2[aa] * ( exp(-a3[aa]*(rf150-a1[aa])) - exp(-a4[aa]*(rf150-a1[aa])) );
    const double beta  = b0[bb] + b2[bb] * ( exp(-b3[bb]*(rf150-b1[bb])) - exp(-b4[bb]*(rf150-b1[bb])) );
    const double corr = 0.135 / ( 2*alpha*etc + beta);
   
    /* 
    static int ncalls=0;
    if (ncalls++==0)
    {
        TF1 *f[100];
        TH2F *h=new TH2F("h","",1000,61510,61930,1000,-0.005,0.005);
        //TH2F *h=new TH2F("h","",1000,61510,61930,1000,0.09,0.17);
        h->DrawClone();
        for (unsigned int ii=0; ii<bXLO.size(); ii++)
        {
            f[ii]=new TF1(Form("f__%d",ii),EG6ICecorr39,bXLO[ii],bXHI[ii],0);
            f[ii]->SetLineColor(kRed);
            if (ii==0) f[ii]->Draw("RSAME");
            else       f[ii]->Draw("RSAME");
        }
    }
    return alpha;
    return beta;
    */

    return corr*etc;
}


double SpecsPCORR::EG6ECsampcorr(const int sector,const double rf150)
{
//    sector == Sector (1-6)
//    rf150  == Run# + File# / 150
//
//    Returns the sampling fraction parameterized by piecewise functions of the form:
//    E0 + A0 * ( exp(-ALPHA1*(x-X0)) + exp(-ALPHA2*(x-tX0)) )
//    where x = Run# + File#/150
//    Negative return value is an error.
//    Reads parameters from EG6ECsampcorr.dat
//
    static bool initialized=0;
    static const int npar=5;

    // table columns:
    static const int ncols=2+6*npar;
    static vector <double> tXLO,tXHI;  // <--- Run# + File#/150
    static vector <double> tE0[6],tX0[6],tA0[6],tALPH1[6],tALPH2[6];

    // read the data table:
    if (!initialized)
    {
        const char* stub="EG6ECsampcorr.dat";
        const char* cparm=gSystem->Getenv("CLAS_PARMS");
        const TString filename=cparm?Form("%s/%s",cparm,stub):stub;
        
        FILE *fin=fopen(filename.Data(),"r");
        if (!fin)
        {
            fprintf(stderr,"EG6ECsampcorr:  Missing Input File:  %s\n",filename.Data());
            return -9999;
        }

        char buf[1024];
        double tmp[ncols];

        while ((fgets(buf,1024,fin)) != NULL)
        {
            for (int ii=0; ii<ncols; ii++)
            {
                if (ii==0) tmp[ii] = atof((char*)strtok(buf," "));
                else       tmp[ii] = atof((char*)strtok(NULL," "));
            
                if (isnan(tmp[ii]))
                {
                    fprintf(stderr,"EG6ECsampcorr:  Error Reading File:  %s\n",filename.Data());
                    return -9999;
                }
            }

            tXLO.push_back(tmp[0]);
            tXHI.push_back(tmp[1]);

            for (int sec=0; sec<6; sec++)
            {
                tE0[sec]   .push_back(tmp[2+npar*sec+0]);
                tX0[sec]   .push_back(tmp[2+npar*sec+1]);
                tA0[sec]   .push_back(tmp[2+npar*sec+2]);
                tALPH1[sec].push_back(tmp[2+npar*sec+3]);
                tALPH2[sec].push_back(tmp[2+npar*sec+4]);
            }
        }

        initialized=1;
        fclose(fin);
    }


    // check for invalid input:
    if (sector<1 || sector>6)
    {
        fprintf(stderr,"EG6ECsampcorr:  Invalid sector:  %d\n",sector);
        return -9999;
    }
    if (rf150 < tXLO[0] || rf150 >= tXHI[tXHI.size()-1])
    {
        fprintf(stderr,"EG6ECsampcorr:  Invalid rf150:  %f\n",rf150);
        return -9999;
    }

    // binary search:
    int imax=tXLO.size()-1;
    int imin=0;
    int jj=-1;

    while (imax >= imin)
    {
        jj = (imax+imin)/2;
        if      (rf150 >= tXHI[jj]) imin=jj+1;
        else if (rf150 <  tXLO[jj]) imax=jj-1;
        else break;
    }

    // the parameterization:
    const int ss=sector-1;
    const double corr=
        tE0[ss][jj] +
        tA0[ss][jj] * exp(-tALPH1[ss][jj]*(rf150-tX0[ss][jj])) +
        tA0[ss][jj] * exp(-tALPH2[ss][jj]*(rf150-tX0[ss][jj])) ;

    return corr;

}

double SpecsPCORR::E16ECphocorr(const double etot,const int sector,const double rf150)
{
    // for photons from different sectors::::
    static const double p[4]={1.21299,-0.0391713,0.00306506,0.000098974};
    const double corr=p[0]+p[1]/etot+p[2]/pow(etot,2)+p[3]/pow(etot,3);

    const double sampfrac=EG6ECsampcorr(sector,rf150);

    return etot / sampfrac / corr;
}

double SpecsPCORR::SimpleRtpcVertexCorr(const double theta,const double zz)
{   // theta = track polar angle [rad]
    // zz    = vertex [cm]
    // return value:  corrected vertex [cm]
    return zz - (1.53 - 0.017*theta*180/3.14159);
}
double SpecsPCORR::SimpleRtpcThetaCorr(const double theta)
{   // theta = track polar angle [rad]
    // return value: corrected angle [rad]
    //return theta + SimpleRtpcVertexCorr(theta,0) * sin(2*theta) / 4.5;
    return theta + SimpleRtpcVertexCorr(theta,0) * pow(sin(2*theta),2) / 4.5;
}
#endif
