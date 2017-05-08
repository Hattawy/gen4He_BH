#ifndef SPECSEG6_hh
#define SPECSEG6_hh
#include "TString.h"
#include "TRegexp.h"
#include "SpecsRTPC.hh"
#include "SpecsPDG.hh"
#include "SpecsFID.hh"
#include "SpecsPCORR.hh"
#include "SpecsELOSS.hh"
#include "SpecsPID.hh"
#include "SpecsGEO.hh"
#include <iostream>
#include <cstdio>
#include <vector>
//////////////////////////////////////////////////////////
class SpecsEG6 {
private:
    int runnumber;
    vector <vector <int> > badfilesC;
    vector <vector <int> > badfilesI;
    vector <int> badrunsC;
    vector <int> badrunsI;
public:
    SpecsPDG PDG;
    SpecsRTPC RTPC;
    SpecsELOSS ELOSS;
    SpecsPCORR PCORR;
    SpecsFID FID;
    SpecsPID PID;
    SpecsGEO GEO;
    SpecsEG6();
    virtual ~SpecsEG6();
    virtual inline int GetTargetCenter(){ return -64; }
    virtual inline void SetRunNumber(const int runno){ runnumber=runno; };
    virtual int SetRunNumber(const char* filestub);
    virtual int GetRunNumber(const char* filestub);
    virtual int GetFileNumber(const char* filestub);
    virtual inline int GetRunNumber(){return runnumber;};
    virtual double GetBeamEnergy(const int runno=-1);
    virtual double GetTargetMass(const int runno=-1);
    virtual double GetTargetCharge(const int runno=-1);
    virtual int    GetTargetPID(const int runno=-1);
    virtual int    GetTargetType(const int runno=-1);
    virtual int    GetHelicity(const int h10hel,const int runno);
    virtual int    GetHWP(const int runno);
//    virtual double VertexCorrCLAS_2(const double z,const double theta,const double phi);
//    virtual double VertexCorrCLAS_2(const double z,const double cx,const double cy,const double cz);
//    virtual double ThetaCorrCLAS_2(const double theta,const double phi);
//    virtual double ThetaCorrCLAS_2(const double cx,const double cy,const double cz);
    virtual double VertexCorrCLAS(const double z,const double theta,const double phi);
    virtual inline double VertexCorrCLAS(const double z,const double cx,const double cy,const double cz){return VertexCorrCLAS(z,acos(cz),atan2(cy,cx));};
    virtual void  Print();
    virtual double PhiCorrRTPCoffset(const double phi,const double x0,const double y0);
    virtual double MomCorrRTPCoffset(const double mom,const double theta,const double phi,const double x0,const double y0);

    virtual double ev2file(const int runno,const int evno);
    virtual double nevents(const int runno);
    virtual double EG6ECsampcorr(const int sector,const int runno,const int evno);
    virtual double EG6ICecorr(const double etc,const int runno,const int evno);
    virtual double E16ECphocorr(const double etot,const int sector,const int runno,const int evno);
    virtual bool   LoadVector1D(const char* filename,vector <int> &xx);
    virtual bool   LoadVector2D(const char* filename,vector <vector <int> > &xx);
    virtual bool   LoadBadRuns(); 
    virtual bool   IsBadRun(const int runno,const TString type);
    virtual bool   IsBadFile(const int runno,const int fileno,const TString type);
    virtual void   PrintBadRuns();
    
    virtual double ElectronSampling(const int sector,const int runno,const int evno,
            const double mom);
            //const double ece,
            //const double nsigma);
    virtual void ElectronSampling(const double mom,double& mean,double& sigma);

};
SpecsEG6::SpecsEG6()
{
    runnumber=0;
    LoadBadRuns();
}
SpecsEG6::~SpecsEG6(){}
//////////////////////////////////////////////////////////
int SpecsEG6::GetFileNumber(const char* filestub)
{
    // convert standard EG6 filename into file number

    TString ss=filestub;
    const TRegexp rr3="[0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9]_";
    const Ssiz_t ind3=ss.Index(rr3);
    if ((size_t)ind3==string::npos) 
    {
        const TRegexp rr2="[0-9][0-9][0-9][0-9][0-9]_[0-9][0-9]_";
        const Ssiz_t ind2=ss.Index(rr2);
        if ((size_t)ind2==string::npos) return -1;
        ss.Remove(0,ind2+6);
        ss.Resize(2);
    }
    else
    {
        ss.Remove(0,ind3+6);
        ss.Resize(3);
    }
    return ss.Atoi();
}
int SpecsEG6::GetRunNumber(const char* filestub)
{
    TString ss=filestub;
    TRegexp rr="[0-9][0-9][0-9][0-9][0-9]";  // any five numbers in a row will match!
    Ssiz_t ind=ss.Index(rr);
    if ((size_t)ind==string::npos) return -1;
    ss.Remove(0,ind);
    ss.Resize(5);
    return ss.Atoi();
}
int SpecsEG6::SetRunNumber(const char* filestub)
{
    const int runno=GetRunNumber(filestub);
    runnumber=runno;
    return runnumber;
//    return runno;
//    if (runno<=0) return 0;
//    else          SetRunNumber(runno);
//    return 1;
}
double SpecsEG6::GetBeamEnergy(const int runno)
{
    if (runno<0)
    {
        if (runnumber<0)
        {
            cerr<<"SpecsEG6::GetBeamEnergy --  Run Number Undefined"<<endl;
            return -99999;
        }
        return GetBeamEnergy(runnumber);
    }
    if (runno>=61225 && runno<=61482) return 1.206;
    if (runno>=61001 && runno<=61210) return 5.776;
    if (runno>=61483 && runno<=61496) return 5.776;
    if (runno>=61510 && runno<=61782) return 6.064;
    if (runno==61789 || runno==61790) return 1.269;
    if (runno>=61791 && runno<=61930) return 6.064;
    if (runno>=61932 && runno<=61966) return 1.269;
    cerr<<"Unkown Beam Energy for Run Number: "<<runno<<endl;
    return 0;
}
int SpecsEG6::GetTargetType(const int runno)
{
    if (runno<0) return GetTargetType(runnumber);
    // hydrogen target:
    if (runno>=61963 && runno<=61966) return 0;
    if (runno>=61432 && runno<=61446) return 0;
    // helium-4 target:
    return 1;
}
double SpecsEG6::GetTargetMass(const int runno)
{
    if (GetTargetType(runno)==0) return PDG.Mass("proton");
    else                         return PDG.Mass("alpha"); 
}
double SpecsEG6::GetTargetCharge(const int runno)
{
    if (GetTargetType(runno)==0) return PDG.Charge("proton");
    else                         return PDG.Charge("alpha"); 
}
int SpecsEG6::GetTargetPID(const int runno)
{
    if (GetTargetType(runno)==0) return PDG.PidG3("proton");
    else                         return PDG.PidG3("alpha"); 
}
/*
double SpecsEG6::VertexCorrCLAS_2(const double z,const double theta,const double phi)
{
    static const double up[3]={0.191,0.168,-79.55}; // (x,y,z) cm
    static const double dn[3]={0.177,0.007,-48.56}; // (x,y,z) cm
    static const double dxdz = (dn[0]-up[0])/(dn[2]-up[2]);
    static const double dydz = (dn[1]-up[1])/(dn[2]-up[2]);
    const double xoff = up[0] + (z-up[2]) * dxdz;
    const double yoff = up[1] + (z-up[2]) * dydz;
    const double roffset = sqrt(xoff*xoff+yoff*yoff);
    const double phioffset = atan2(yoff,xoff);
    return z - roffset * cos(phi - phioffset) / tan(theta);
}
double SpecsEG6::VertexCorrCLAS_2(const double z,const double cx,const double cy,const double cz)
{
    return VertexCorrCLAS_2(z,acos(cz),atan2(cy,cx));
}
double SpecsEG6::ThetaCorrCLAS_2(const double theta,const double phi)
{
    static const double up[3]={0.191,0.168,-79.55}; // (x,y,z) cm
    static const double dn[3]={0.177,0.007,-48.56}; // (x,y,z) cm
    static const double dx = dn[0]-up[0];
    static const double dy = dn[1]-up[1];
    static const double dz = dn[2]-up[2];
    static const double beamtheta = atan(sqrt(dx*dx+dy*dy)/dz);
    static const double beamphi = atan2(dy,dx);
    return theta - beamtheta*cos(phi - beamphi);
}
double SpecsEG6::ThetaCorrCLAS_2(const double cx,const double cy,const double cz)
{
    return ThetaCorrCLAS_2(acos(cz),atan2(cy,cx));
}
*/
double SpecsEG6::VertexCorrCLAS(const double z,const double theta,const double phi)
{
    if (fabs(sin(theta)) < 1e-2) return -9999;
    
    static const int nr=4;
    static const int runranges[nr]={61483,61580,61850,99999};
    static const double xx[nr]={0.155, 0.237, 0.27,  0.30};
    static const double yy[nr]={0.029,-0.040, 0.04, -0.04};
    int therange=0;
    for (int ii=0; ii<nr; ii++) {
        if (runnumber < runranges[ii]) {
            therange=ii;
            break;
        }
    }
    const double rr=sqrt(pow(xx[therange],2)+pow(yy[therange],2));
    const double phi0=atan2(yy[therange],xx[therange])+3.141593;
    return z-rr*cos(phi-phi0)/tan(theta);
/*    
    // z in cm, theta/phi in radians
    const double roffset = runnumber<61483 ? 0.158 : 0.24;   //cm
    const double phi0_    = runnumber<61483 ? 190.5 : -189.6; // deg
    return z-roffset*cos(phi-phi0_/57.29578)/tan(theta);
    return -9999;
*/
}
void SpecsEG6::Print()
{
    cerr<<"ebeam  = "<<GetBeamEnergy()<<" GeV"<<endl;
    cerr<<"ztgt   = "<<GetTargetCenter()<<" cm"<<endl;
    cerr<<"mass   = "<<GetTargetMass()<<" GeV"<<endl;
    cerr<<"charge = "<<GetTargetCharge()<<endl;
    cerr<<"pid    = "<<GetTargetPID()<<endl;

}
double SpecsEG6::PhiCorrRTPCoffset(const double phi,const double x0,const double y0)
{
    // input phi in degrees
    // input x0,y0 in mm
    static const double ratio=2.56;
    const double phir=phi*3.14159/180.;
    return phi  -   ratio * (-y0*cos(phir) + x0*sin(phir));
}
double SpecsEG6::MomCorrRTPCoffset(const double mom,const double theta,const double phi,const double x0,const double y0)
{
    // input momentum in MeV
    // input phi,theta in degrees
    // input x0,y0 in mm
    static const double ref_mom=130;
    static const double ref_theta=78;
    static const double ref_momt=ref_mom*sin(ref_theta*3.14159/180.);
    static const double ratio=0.118;

    const double momt=mom*sin(theta*3.14159/180.);
    const double ratio2=ratio*(momt/ref_momt);
    const double phir=phi*3.14159/180.;
    const double dpoverp= - ratio2*(-y0*cos(phir)+x0*sin(phir));
    return (dpoverp+1.0)*mom;
}
int SpecsEG6::GetHWP(const int runno)
{
    // For EG6 Run Numbers, return position of Half-Wave Plate:
    //  1 is IN
    //  0 is OUT
    //  else UNKNOWN

    // before first HWP change, HWP is OUT(0):
    const int initialstate = 0;

    // these are the first runs AFTER a HWP change:
    static const int nchanges=20;
    static const int changes[nchanges]=
    {
        61626,61642,61670,61688,61691,
        61717,61734,61745,61759,61770,
        61791,61793,61802,61819,61836,
        61852,61873,61886,61906,61931
    };
    
    // FIRST check for bad run numbers where HWP is uncertain:
    if (runno < 61510 || runno > 61930) 
    {
        cerr<<"HWP unknown for run "<<runno<<endl;
        return -99999;
    }
    else if (runno == 61733)
    {
        cerr<<"HWP changed in middle of run "<<runno<<endl;
        return -99999;
    }
    else if (runno > 61779 && runno < 61791) 
    {
        cerr<<"HWP unknown for run "<<runno<<endl;
        return -99999;
    }

    // THEN find the HWP position:   (should replace with binary search)
    for (int ii=0; ii<nchanges; ii++)
        if (runno < changes[ii])
            return (initialstate ^ ii) & 1;

    return -99999;
}

int SpecsEG6::GetHelicity(const int h10hel,const int runno)
{
    // returns 0/1 if helicity is certain, else -99999

    // if it's negative, the true helicity could not be determined:
    if (h10hel<0) return -99999;
    
    // if bits 5-24 are zero, this is an incomplete quartet:
    if ( !(h10hel & 0xFFFFF0) ) return -99999;
  
    // HWP position 0/1 for OUT/IN, else unknown:
    const int hwp=GetHWP(runno);
    if (hwp!=0 && hwp!=1) return -99999;
  
    return (hwp ^ h10hel) & 1;
}

/*
int SpecsEG6::GetHelicity(const int runno,const int fileno,const int evno)
{
    // Get the helicity value from the 4th column of the helcor ascii files

    static int prevrunno=-1;
    static int prevfileno=-1;
    
    static vector <double> flag;
    static vector <int> evlo;
    static vector <int> evhi;
    static vector <int> hel;
    
    static TString cfile;

    if (runno!=prevrunno || fileno!=prevfileno)
    {

        prevrunno=-1;
        prevfileno=-1;

        // this doesn't actually free the memory (which could be slow),
        // just set size() to zero
        flag.clear();
        evlo.clear();
        evhi.clear();
        hel.clear();

        const char dir[100]="/music/clas1/clas/eg6/pass1v1/6gev/helcor";
        if (fileno<100) cfile = Form("%s/hel_%.5d_%.2d_pass1v1.txt",dir,runno,fileno);
        else            cfile = Form("%s/hel_%.5d_%.3d_pass1v1.txt",dir,runno,fileno);
        
        cout<<"SpecsEG6::GetHelicity   Reading Helicty: "<<cfile<<endl;
        
        FILE *fin=fopen(cfile,"r");
        if (!fin)
        {
            cerr<<"SpecsEG6::GetHelicity  Missing File: "<<cfile<<endl;
            return -999;
        }

        char line[1000];
        char* str;
        int nrows=0;

        while (fgets(line,1000,fin))
        {
            str = (char *) strtok(line," ");  flag.push_back(atof(str));
            str = (char *) strtok(NULL," ");  evlo.push_back(atoi(str));
            str = (char *) strtok(NULL," ");  evhi.push_back(atoi(str));
            str = (char *) strtok(NULL," ");   hel.push_back(atoi(str));
            nrows++;
        }

        fclose(fin);

        if (flag.size() != evlo.size() ||
            flag.size() != evhi.size() ||
            flag.size() != hel.size())
        {
            cerr<<"SpecsEG6::GetHelicity  Read File Error: "<<cfile<<endl;
            return -999;
        }

        if (nrows>2e4)
        {
            cerr<<"SpecsEG6::GetHelicity  Warning:  >20K lines in "<<cfile<<endl;
        }

        prevrunno=runno;
        prevfileno=fileno;
    }

    if (evno < evlo[0] || evno > evhi[evhi.size()-1])
    {
        cerr<<"SpecsEG6::GetHelicity  Out of Range Event #"<<evno<<" in "<<cfile<<endl;
        return -999;
    }

    // binary search:

    const double fraction = (double)(evno-evlo[0])/(evhi[evhi.size()-1]-evlo[0]);
    
    int jj=fraction*(evlo.size()-1); // initial guess
    int min=0; 
    int max=flag.size()-1;

    cout<<"START:  "<<jj<<" "<<min<<" "<<max<<" "<<fraction<<endl;

    int ifoundit=-1;
    int niter=0;

    while (1)
    {
        if (evno >= evlo[jj] && evno <= evhi[jj]) 
        {
            ifoundit=jj;
            break;
        }
        
        if   (evno < evlo[jj])  max=jj;
        else                    min=jj;

        jj = min+(max-min+1)/2;

        if (++niter>100000) 
        {
            cerr<<"SpecsEG6::GetHelicity  BINARY SEARCH NEEDS FIX (Event #"<<
                evno<<" in File #"<<fileno<<" in Run #"<<runno<<")"<<endl;
            break;
        }

        cout<<jj<<" "<<min<<" "<<max<<endl;
    }

    if (fabs(flag[ifoundit])>1e-8) return hel[ifoundit];
    return -9999;
}
*/

/*
double SpecsEG6::evno2fileno(const int runno,const int evno)
{
    // Convert from h10's evntid to the file#.
    // Return -1 if runno is not in list of runs.
    // Cast the return value to int to get real file#.
    
    static const int nruns=40;
    static const int runs[nruns] =
    {
        61532, 61536, 61546, 61548, 61558, 61581, 61587, 61610,
        61615, 61623, 61629, 61633, 61638, 61640, 61641, 61655,
        61673, 61679, 61685, 61695, 61704, 61711, 61712, 61713,
        61714, 61718, 61725, 61731, 61763, 61792, 61793, 61797,
        61798, 61802, 61813, 61826, 61854, 61856, 61857, 61876
    };
    static const double slopes[nruns]=
    {
       287203,318112,325309,281010,260842,276504,301957,315423,
       298713,307596,305713,298374,274936,289847,297758,301145,
       304949,296217,306388,471614,318827,312822,328162,320280,
       325433,329275,325785,317948,361102,379684,349286,353864,
       336509,346545,345352,349426,355034,397348,390995,374142};

    if (runno<runs[0] || runno>runs[nruns-1]) return -1;

    // binary search:
    int imax=nruns-1;
    int imin=0;
    int itry=-1;
    while (imax>=imin)
    {
        itry=(imax+imin)/2;
        if      (runno>runs[itry]) imin=itry+1;
        else if (runno<runs[itry]) imax=itry-1;
        else                       break;
        itry=-1;
    }
    
    return itry<0 ? -1 : evno/slopes[itry];
}
*/
double SpecsEG6::ev2file(const int runno,const int evno)
{
//    runno = Run Number (61510-61930)
//    evno  = Event Number
//
//    Returns the interpolated file number.
//    Negative return value is an error.  (e.g. run not found)
//    Reads parameters from EG6ev2file.dat
//
    static bool initialized=0;

    // table columns:
    static const int ncols=4;
    static vector <int> tRUN,tFILE,tEVMIN,tEVMAX;

    // read the data table:
    if (!initialized)
    {
        const char* stub="EG6ev2file.dat";
        const char* cparm=gSystem->Getenv("CLAS_PARMS");
        const TString filename=cparm?Form("%s/%s",cparm,stub):stub;
        
        FILE *fin=fopen(filename.Data(),"r");
        if (!fin)
        {
            fprintf(stderr,"EG6ev2file:  Missing Input File:  %s\n",filename.Data());
            return -9999;
        }

        char buf[256];
        int tmp[ncols];

        while ((fgets(buf,256,fin)) != NULL)
        {
            for (int ii=0; ii<ncols; ii++)
            {
                if (ii==0) tmp[ii] = atoi((char*)strtok(buf," "));
                else       tmp[ii] = atoi((char*)strtok(NULL," "));
            
                if (isnan(tmp[ii]))
                {
                    fprintf(stderr,"EG6ev2file:  Error Reading File:  %s\n",filename.Data());
                    return -9999;
                }
            }

            tRUN .push_back(tmp[0]);
            tFILE.push_back(tmp[1]);
            tEVMIN.push_back(tmp[2]);
            tEVMAX.push_back(tmp[3]);
        }

        initialized=1;
        fclose(fin);
    }
 
    // check for invalid input:
    if ( runno<tRUN[0] || runno > tRUN[tRUN.size()-1])
    {
        fprintf(stderr,"EG6ev2file:  Run Out of Range:  %d\n",runno);
        return -9999;
    }
//    if ( evno<0 )
//    {
//        fprintf(stderr,"EG6ev2file:  Invalid Event #:  %d\n",evno);
//        return -9999;
//    }

    // binary search:
    int imax=tRUN.size()-1;
    int imin=0;
    int jj=-1;

    if (evno<0)
    {
        while (imax>=imin)
        {
            jj = (imax+imin)/2;
            if      (runno > tRUN[jj]) imin=jj+1;
            else if (runno < tRUN[jj]) imax=jj-1;
            else break;
        }
        if (runno != tRUN[jj])
        {
            fprintf(stderr,"EG6ev2file:  Run Not Found:  %d\n",runno);
            return -9999;
        }
        while (1)
        {
            if (jj+1 >= (int)tRUN.size()) return tEVMAX[tEVMAX.size()-1];
            if (runno != tRUN[jj+1])      return tEVMAX[jj];
            jj++;
        }
    }

    while (imax >= imin)
    {
        jj = (imax+imin)/2;

        if      (runno > tRUN[jj]) imin=jj+1;
        else if (runno < tRUN[jj]) imax=jj-1;
        else
        {
            if      (evno > tEVMAX[jj]) imin=jj+1;
            else if (evno < tEVMIN[jj]) imax=jj-1;
            else break;
        }
    }

    if (runno != tRUN[jj])
    {
        fprintf(stderr,"EG6ev2file:  Run Not Found:  %d\n",runno);
        return -9999;
    }

    // interpolate:
    const double nev=tEVMAX[jj]-tEVMIN[jj]+1;
    return tFILE[jj] + (evno-tEVMIN[jj])/nev;
}
double SpecsEG6::nevents(const int runno)
{
    // should just make a separate table for this
    return ev2file(runno,-1);
}
double SpecsEG6::EG6ECsampcorr(const int sector,const int runno,const int evno)
{
//    sector == Sector (1-6)
//    runno  == Run Number (61510-61930)
//    evno   == Event Number
//
    const double file=ev2file(runno,evno);
    if (file<0) return -9999;
    const double rf150 = runno + file/150.0;
    return PCORR.EG6ECsampcorr(sector,rf150);
}
double SpecsEG6::EG6ICecorr(const double etc,const int runno,const int evno)
{
//    ee     == IC energy (h10's etc)
//    runno  == Run Number (61510-61930)
//    evno   == Event Number
//
    const double file=ev2file(runno,evno);
    if (file<0) return -9999;
    const double rf150 = runno + file/150.0;
    return PCORR.EG6ICecorr(etc,rf150);
}
double SpecsEG6::E16ECphocorr(const double etot,const int sector,const int runno,const int evno)
{
    const double file=ev2file(runno,evno);
    if (file<0) return -9999;
    const double rf150 = runno + file/150.0;
    return PCORR.E16ECphocorr(etot,sector,rf150);
}

bool SpecsEG6::IsBadRun(const int runno,const TString type)
{
    vector <int> *badruns = 
        (type.Contains("Inc") || type.Contains("inc"))
        ? &badrunsI : &badrunsC;
    int imin=0; 
    int imax=badruns->size();
    while (imax >= imin)
    {
        if (imax<0 || imin>=(int)badruns->size()) break;
        int ii = (imax+imin)/2;
        if      (runno > badruns->at(ii)) imin = ii+1;
        else if (runno < badruns->at(ii)) imax = ii-1;
        else return 1;
    }
    return 0;
}
bool SpecsEG6::IsBadFile(const int runno,const int fileno,const TString type)
{
    vector <vector <int> > *badfiles = 
        (type.Contains("Inc") || type.Contains("inc"))
        ? &badfilesI : &badfilesC;
            
    int imin=0; 
    int imax=badfiles->size();
    while (imax >= imin)
    {
        if (imax<0 || imin>=(int)badfiles->size()) break;
        int ii = (imax+imin)/2;
        if      (runno  > badfiles->at(ii)[0]) imin = ii+1;
        else if (runno  < badfiles->at(ii)[0]) imax = ii-1;
        else if (fileno > badfiles->at(ii)[1]) imin = ii+1; 
        else if (fileno < badfiles->at(ii)[1]) imax = ii-1; 
        else  return 1; 
    }
    return 0;
}
bool SpecsEG6::LoadVector1D(const char* filename,vector <int> &xx)
{
    FILE *fin=fopen(filename,"r");
    if (!fin)
    {
        fprintf(stderr,"Missing File:   %s\n",filename);
        return 0;
    }
    int itmp;
    char buf[1024];
    while ((fgets(buf,1024,fin)) != NULL)
    {
        itmp = atoi((char*)strtok(buf," "));
        if (isnan(itmp)) 
        {
            fprintf(stderr,"Error Reading File:  %s\n",filename);
            return 0;
        }
        xx.push_back(itmp);
    }
    return 1;
}
bool SpecsEG6::LoadVector2D(const char* filename,vector <vector <int > > &xx)
{
    FILE *fin=fopen(filename,"r");
    if (!fin)
    {
        fprintf(stderr,"Missing File:   %s\n",filename);
        return 0;
    }
    static const int ncols=2;
    char buf[1024];
    int itmp[ncols];
    while ((fgets(buf,1024,fin)) != NULL)
    {
        vector <int> vtmp;
        for (int ii=0; ii<ncols; ii++)
        {
            if (ii==0) itmp[ii] = atoi((char*)strtok(buf," "));
            else       itmp[ii] = atoi((char*)strtok(NULL," "));

            if (isnan(itmp[ii]))
            {
                fprintf(stderr,"Error Reading File:  %s\n",filename);
                return 0;
            }
        }
        vtmp.push_back(itmp[0]);
        vtmp.push_back(itmp[1]);
        xx.push_back(vtmp);
    }
    return 1;
}
bool SpecsEG6::LoadBadRuns()
{
    static bool initialized=0;
    if (initialized) return 1;
    const TString cparm=gSystem->Getenv("CLAS_PARMS");
    if (!LoadVector1D(cparm+"/runlists/CohBadRuns.txt",badrunsC))   return 0;
    if (!LoadVector1D(cparm+"/runlists/IncohBadRuns.txt",badrunsI)) return 0;
    if (!LoadVector2D(cparm+"/runlists/CohBadFileList.txt",badfilesC))   return 0;
//    if (!LoadVector2D(cparm+"/runlists/IncohBadFileList.txt",badfilesI)) return 0;
    return 1;
}
void SpecsEG6::PrintBadRuns()
{
//    for (unsigned int ii=0; ii<badrunsI.size(); ii++)
//    {
//        cout<<ii<<" "<<badrunsI[ii]<<endl;
//    }
    for (unsigned int ii=0; ii<badfilesC.size(); ii++)
    {
        cout<<ii<<" "<<badfilesC[ii][0]<<" "<<badfilesC[ii][1]<<endl;
    }
}
double SpecsEG6::ElectronSampling(
        const int sector,
        const int runno,
        const int evno,
        const double mom)
//        const double ece,
//        const double nsigma)
{
    double mean0=0,sigma0=0;
    ElectronSampling(1.65,mean0,sigma0);

    double mean=0,sigma=0;
    ElectronSampling(mom,mean,sigma);

    // momentum-averaged for this sector
    const double sampS=EG6ECsampcorr(sector,runno,evno);

    const double sampC=sampS*mean/mean0;

    return sampC;//fabs( sampC - ece/mom ) > nsigma*sigma;
}

void SpecsEG6::ElectronSampling(const double mom,double& mean,double& sigma)
{
    // based upon run 61776
    static const int npar=5;
    static const double mpar[npar]={0.2639010, 0.0619894,-0.0205942, 0.00392891,-0.000304636};
    static const double spar[npar]={0.0563971,-0.0355622, 0.0185887,-0.00442158, 0.000371359};
    const double mom2 = (mom>4.6) ? 4.6 : mom;
    mean=0;
    sigma=0;
    for (int ii=0; ii<npar; ii++)
    {
        mean  += mpar[ii]*pow(mom2,ii);
        sigma += spar[ii]*pow(mom2,ii);
    }
    return;
}

#endif

