///////////////////////////////////////
// 4He DVCS and DVMP Event Generator //
//                                   // 
//   k0 = beam                       //
//   q0 = virtual photon             //
//   p0 = target nucleus/nucleon     //
//   k1 = scattered electron         //
//   q1 = real photon or pi0         //
//   p1 = recoil nucleus/nucleon     //
//                                   //
//   L = lab frame                   //
//   C = q0+p0 CM frame              //
//   S = smeared by resolution       //
//   X = reconstructed via missing   //
//                                   //
//  Weighting by FX's pDVCS Model    //
//        Modified for 4He           //
///////////////////////////////////////

#include "cheaders.h"
#include "fermi4He.C"
#include "fastmc.C"
#include "dvcstree.C"
#include "init_track.C"
#include "SpecsGEO.hh"

#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>

using namespace std;

dvcstree_t DVCS_;
//SpecsGEO GEO;

TRACK trk;
#define SWAPINTS(a,b) {a ^= b; b ^= a; a ^= b;}

const double EBEAM=6.064;
const double TCURR=1900.0;
const double R2D=57.29578;


const double MALPHA=3.7274;
const double MPROTON=0.93827;
const double MPI0=0.13498;

//Nathan kienamtics ranges
//const double Q2LIMS[2]={1.0,3.0};
//const double XBLIMS[2]={0.1,0.25};  // NOTE: xb is for proton
//const double TLIMS[2]={0.05,0.3};//0.2};
// Coh : Moh
//const double Q2LIMS[2]={1.0,3.0};
//const double XBLIMS[2]={0.1,0.3};
//const double TLIMS[2]={0.05,0.3};
//InCoh : Moh
const double Q2LIMS[2]={1.0,5.0};
const double XBLIMS[2]={0.1,0.55};  // NOTE: xb is for proton
const double TLIMS[2]={0.05,2.5};

const double VZLIMS[2]={-64-10,-64+10};

const double WCUT=2.0;

// from FX's pDVCS analysis note:
struct dvcsvar_t { double q2,t,x,phi; };
struct dvcspar_t { double q2_0,alpha,b,beta,c,d; };
double xsec_e1dvcs(const dvcsvar_t &x,const dvcspar_t &p)
{
    return pow(p.q2_0/x.q2,p.alpha)     *  1./pow(1.+p.b*x.t,p.beta) *
           1./(1.+pow((x.x-0.5)/p.c,2)) * (1.-p.d*(1.-cos(x.phi)));
}
double asy_alu(const dvcsvar_t &x,const double alpha,const double beta)
{
    return alpha * sin(x.phi) / (1+beta*cos(x.phi));
}
                               //  q20  alpha      b     beta    c     d
//const dvcspar_t PAR_COHDVCS   = {  1.0,  2.0,  -8.8,     7.3,   0.3,  0.4}; // b,beta from t-fit to tr
const dvcspar_t PAR_COHDVCS   = {  1.0,  2.0,  -0.667,  35.1,   0.3,  0.4}; // b,beta from t-fit to tc
const dvcspar_t PAR_INCDVCS   = {  1.0,  1.5,  -1.408,   4.0,   0.2,  0.4};
const dvcspar_t PAR_COHPI0    = {  1.0,  3.0,  -8.8,     7.3,   0.3,  0.0};
const dvcspar_t PAR_INCPI0    = {  1.0,  3.0,  -1.408,   2.0,   0.5,  0.0};
                                 



// print the progress of the genrated events
void ProgressMeter(const double total,const double current)
{
    if (!isatty(fileno(stdout))) return;
    static const int maxdots=40;
    const double frac = current/total;
    int ii=0;
    printf("%3.0f%% [",frac*100);
    for ( ; ii < frac*maxdots; ii++) printf("=");
    for ( ; ii < maxdots;      ii++) printf(" ");
    printf("]\r");
    fflush(stdout);
 }

int main(int argc,char** argv)
{
    const char* usage="\ngen4He [options]\n"
                      "\t-n NEVENTS\n"
                      "\t-o root_OUTFILE\n"
                      "\t-k dat_OUTFILE\n"
                      "\t-s SEED\n"
                      "\t   0 == use current time\n"
                      "\t-f disables fermi motion\n"
                      "\t-a A_{LU}(90deg)\n"
                      "\t-r RTPC\n"
                      "\t   0 == Detect Recoil in CLAS\n" 
                      "\t   1 == Detect Recoil in RTPC\n"
                      "\t-t TYPE\n"
                      "\t   0 == Coherent DVCS\n"
                      "\t   1 == Coherent Pi0\n"
                      "\t   2 == Incoherent DVCS\n"
                      "\t   3 == Incoherent Pi0\n"
                      "\t-T TYPE\n"
                      "\t   pi0/dvcs+coh/inc\n";

    int itmp;
    int seed=1234567;
    int nevents=1000000;
    bool inRTPC=1;
    bool doFermi=1;
    int channel=-1;
    string ofile_root="";
    int n_dat_files = 1 + (int) nevents/4000;
    string ofile_dat[n_dat_files];
    string shannel="";
    float asymmetry=0;

    int elc_id = 11;
    int photon_id = 22;
    int proton_id = 2212;
    int alpha_id = 3312;

    // interpret command line:
    while ( (itmp=getopt(argc,argv,"a:n:s:t:T:o:k:r:fh")) != -1)
    {
        switch (itmp)
        {
            case 't':
                {
                    char* ctmp1=optarg;
                    char* ctmp2=optarg;
                    int errno=0;
                    channel=(int)strtoul(ctmp1,&ctmp2,10);
                    if (errno!=0 || ctmp1==ctmp2 || *ctmp2!=0 ||
                            channel<0 || channel>3)
                    {
                        cerr<<"Invalid Argument for -t:  "<<optarg<<endl;
                        cout<<usage<<endl;
                        exit(1);
                    }
                }
                break;
            case 'T':
                shannel=optarg;
                break;
            case 'r':
                inRTPC=atoi(optarg);
                break;
            case 'n':
                nevents=atoi(optarg);
                break;
            case 'f':
                doFermi=0;
                break;
            case 's':
                seed=atoi(optarg);
                break;
            case 'o':
                ofile_root=optarg;
                break;             
            case 'a':
                asymmetry=atof(optarg);
                break;
            default:
                cout<<usage<<endl;
                exit(1);
        }
    }

    if ( (channel<0 && shannel=="") || (channel>=0 && shannel!=""))
    {
        cerr<<endl<<"Must specify channel with EXACTLY ONE of -T or -t."<<endl;
        cout<<usage<<endl;
        exit(1);
    }

    // interpret -T flag:
    if (shannel!="")
    {
        const int idvcs=shannel.find("dvcs")!=string::npos;
        const int ipi0 =shannel.find("pi0") !=string::npos;
        const int icoh =shannel.find("coh") !=string::npos;
        const int iinc =shannel.find("inc") !=string::npos;
        if ( (idvcs+ipi0 != 1) || (icoh+iinc != 1) )
        {
            cerr<<"-T requires EXACTLY ONE of each of pi0/dvcs and coh/inc."<<endl;
            cout<<usage<<endl;
            exit(1);
        }
        channel=2*iinc+ipi0;
    }
   
    if (fabs(asymmetry)>1)
    {
        cerr<<"Bad asymmetry:  "<<asymmetry<<endl;
        exit(1);
    }

    // interpret channel:
    const bool isPi0 = channel%2;
    const bool isInc = channel/2;

    if (!isInc && !inRTPC)
    {
        cerr<<"Cannot Detect 4He in CLAS."<<endl;
        cout<<usage<<endl;
        exit(1);
    }

    // set filename if not set by user:
    //root file
    if (ofile_root=="") 
    {
        ofile_root = "gen4He_";
        ofile_root += isInc  ? "inc" : "coh";
        ofile_root += isPi0  ? "pi0" : "dvcs";
        ofile_root += inRTPC ? "_R"  : "_C";
        if (isInc && !doFermi) ofile_root += "_noF";
        ofile_root += ".root";
    }
    // dat files the fist is named as 0
    for( int ii=0; ii<n_dat_files; ii++)
     {
       if(ofile_dat[ii]=="")
        {
          ofile_dat[ii] = "gen4He_";
          ofile_dat[ii] += isInc  ? "inc" : "coh";
          ofile_dat[ii] += isPi0  ? "pi0" : "dvcs";
          ofile_dat[ii] += inRTPC ? "_R"  : "_C";
          if (isInc && !doFermi) ofile_dat[ii] += "_noF";
          ofile_dat[ii] += Form("_%d",ii);
          ofile_dat[ii] += ".dat";
         }
       }
   
    // print configuration in pretty format:
    string config="";
    config += isInc  ? "Incoherent" : "Coherent";
    config += isPi0  ? " Pi0" : " DVCS";
    config += isInc  ? ", Proton": ", 4He";
    config += inRTPC ? " in RTPC" : " in CLAS";
    cout<<endl<<"Gen4He Configuration:  "<<config<<endl;
    cout<<"A_{LU}(90deg) = "<<asymmetry<<endl;
    cout<<"Output ROOT File:      "<<ofile_root<<endl<<endl;
    cout<<"Output Dat Files are from:      "<<ofile_dat[0]<< "       to:"<<ofile_dat[n_dat_files-1]<<endl<<endl;
    if (isInc && !doFermi) cout<<"FermiMotion Disabled."<<endl;

    if (!gSystem->AccessPathName(ofile_root.data()))
    {
        cerr<<"Output root File Already Exists."<<endl;
        exit(1);
    }
   
    if (nevents<=0) exit(0);

    RANDOM.SetSeed(seed);
    
    TFile *rfile=new TFile(ofile_root.data(),"CREATE");
    TTree *tree=InitDvcsTree("Tgen4He",DVCS_,config.data());

    // ouput .dat files
    ofstream outfile_[n_dat_files];
    for( int ii=0; ii<n_dat_files; ii++)
     {
       outfile_[ii].open(ofile_dat[ii].data());
      }

    // write parts of the compilations in a txt file
    ofstream ffile;
    ffile.open("out.txt");

    TLorentzVector vtmp;
    double dtmp,max=0;
    int iev=0,nICECerr=0;


    ////////////////////////////////////////////////////////////////////////////

    const TLorentzVector k0(0,0,EBEAM,EBEAM);
    
    // setup some targets:
    const TLorentzVector p0_alpha(0,0,0,MALPHA);
    const TLorentzVector p0_proton(0,0,0,MPROTON);
    
    // choose stationary target for reconstruction:
    const TLorentzVector target = isInc ? p0_proton : p0_alpha;
    
    // choose q1 mass (real photon or pi0):
    const double mq1 = isPi0 ? MPI0 : 0.0;

    // choose cross-section parameters:
    const dvcspar_t xsecpar = isPi0 ? (isInc ? PAR_INCPI0  : PAR_COHPI0 )
                                    : (isInc ? PAR_INCDVCS : PAR_COHDVCS);

    // get maximum cross section for sampling
    const dvcsvar_t xmax = { Q2LIMS[0], -TLIMS[0], 0.3, 0.0 };
    const double xsecmax = xsec_e1dvcs(xmax,xsecpar);
    
   
    ////////////////////////////////////////////////////////////////////////////


    TStopwatch timer;
    timer.Start();

    while (iev<nevents)
    {
        if (nevents>9999 && iev%(nevents/100)==0) ProgressMeter(nevents,iev);



        ////////////  GENERATE KINEMATICS  ////////////////////////////////

        // sample the EM vertex kinematics:
        const double Q2_1 =RANDOM.Uniform(Q2LIMS[0],Q2LIMS[1]);
        const double xb =RANDOM.Uniform(XBLIMS[0],XBLIMS[1]);

        // sample hadron vertex kinematics:
        const double ttt=-RANDOM.Uniform(TLIMS[0],TLIMS[1]);
        const double phi= RANDOM.Uniform(0.0,2*TMath::Pi());

         // sample cross-section:
        dvcsvar_t xsecvar={Q2_1,ttt,xb,phi};
        const double xsec=xsec_e1dvcs(xsecvar,xsecpar);
        if (xsec>max) max=xsec;
        if (RANDOM.Uniform(0,xsecmax)>xsec) continue;

        const int helicity=RANDOM.Integer(2);


        if (fabs(asymmetry) > 1e-6)
        {
            double samp = 1e12;
            xsecvar.phi = RANDOM.Uniform(0.0,2*TMath::Pi());
            while ( samp > 1.0+(2*helicity-1)*asy_alu(xsecvar,asymmetry,0.))
            {
                xsecvar.phi = RANDOM.Uniform(0.0,2*TMath::Pi());
                samp = RANDOM.Uniform(0.0,1.0+asymmetry);
            }
        }


        ////////////  COMPUTE 4-VECTORS ///////////////////////////////////

        // virtual photon energy:
        // NOTE: proton mass used here, because xb is calculated with proton
        const double target_mass = isInc ? MPROTON : MALPHA;
        const int    target_id = isInc ? proton_id : alpha_id;
        const int    parent_A = isInc? 1 : 4;
        const int    target_type = isInc ? 1 : 0; 
        const double nu = Q2_1 / xb / 2 / MPROTON;  // in the case of coherent channel, do we need the mass of 4He!!
                                     
        // scattered electon energy:
        const double eene = EBEAM - nu;
        if (eene<0) continue;

        // scattered electron angle:
        const double sinby2 = sqrt(Q2_1/4/eene/EBEAM);
        if (fabs(sinby2)>1) continue;
        const double ethe = 2*asin(sinby2);

        // electron scattered uniformly in phi:
        const double ephi = RANDOM.Uniform(0,2*TMath::Pi());
        
        // scattered electron 4-vector:
        const TVector3 edir( sin(ethe)*cos(ephi), sin(ethe)*sin(ephi), cos(ethe) );
        const TLorentzVector k1( eene * edir, eene );

        // virtual photon 4-vector:
        const TLorentzVector q0 = ( k0 - k1 );

        // cut on CM energy assuming stationary proton: 
        const double W=(q0+p0_proton).M();
        if (W<WCUT) continue;

        // choose LAB target:

        const TLorentzVector p0 = isInc ? (doFermi ? vfermi4He() : p0_proton ) : p0_alpha;

        // get boost to virtualphoton-nucleus CM:
        const TLorentzVector sysL( q0 + p0 ); // Lab syastem : virtual photon+ rest target
        const TVector3 boostC( sysL.BoostVector() ); // 

        // boost virtualphoton-nucleus system to CM: 
        vtmp=sysL; vtmp.Boost(-boostC);
        const TLorentzVector sysC=vtmp;

        // boost virtualphoton to CM:
        vtmp=q0; vtmp.Boost(-boostC);
        const TLorentzVector q0C=vtmp;

        // momentum/energy of q1 in CM:
        dtmp=pow(sysC.E(),2)-mq1*mq1-p0.M2();
        const double p2q1C=(dtmp*dtmp-4*pow(mq1*p0.M(),2))/4/pow(sysC.E(),2);
        if (p2q1C<0) continue;
        const double pq1C=sqrt(p2q1C);
        const double eq1C=sqrt(p2q1C+mq1*mq1);

        // angle between virtualphoton and q1 in CM:
        const double cosq0q1C = ( (Q2_1+ttt-mq1*mq1)/2 + q0C.E()*eq1C ) / pq1C / q0C.P();
        if (fabs(cosq0q1C)>1) continue;

        // get q1 4-vector in CM:       
        // 1: direction of virtualphoton:
        vtmp=TLorentzVector(q0C.Vect().Unit()*pq1C,eq1C);
        // 2: rotate to make q0-q1 scattering angle:
        vtmp.Rotate(acos(cosq0q1C),(q0C.Vect().Cross(k1.Vect())).Unit());
        // 3: rotate to make phi:
        //vtmp.Rotate(180/R2D-xsecvar.phi,q0C.Vect().Unit()); // THE PHI-CONVENTION IS HERE 
        vtmp.Rotate(phi,q0C.Vect().Unit()); // THE PHI-CONVENTION IS HERE 
        const TLorentzVector q1C=vtmp;

        // boost q1 to LAB:
        vtmp=q1C; vtmp.Boost(boostC);
        const TLorentzVector q1=vtmp;

        // recoil in LAB:
        const TLorentzVector p1L = (q0+p0-q1);


        ////////////  APPLY DETECTOR TO ELECTRON & RECOIL  /////////////////////
        
        // sample z-vertex:
        const double v_x = 0.0;
        const double v_y = 0.0;
        const double v_z =RANDOM.Uniform(VZLIMS[0],VZLIMS[1]);
       
         // recoil RTPC acceptance:
         if (inRTPC) { if (!AcceptRTPC(p1L,v_z))   continue; }
         TLorentzVector p1;
            if (inRTPC) { p1 = SmearRTPC(p1L);  }
             else p1 = p1L; 
        // don't do this until necessary (30% faster now):
        ZeroDvcsTree(DVCS_);
        init_track(trk);


        ////////////  DEAL WITH GAMMA OPTIONS  /////////////////////////////////
        
        // (first) detected gamma in LAB:
         TLorentzVector g, g1, g2;
         
        // DVCS gamma:
        if (!isPi0) 
        {
            g = q1;
            DVCS_.ic = 1;
        }

        // pi0 decay:
        // (anything that's only valid with 2 gammas is here)
        else
        {
            DVCS_.ic = 2;
            // generate gammas in pi0 rest frame:
            const double gene=q1C.M()/2;
            const double gphi=RANDOM.Uniform(0,2*TMath::Pi());
            //const double gthe=RANDOM.Uniform(0,TMath::Pi());
            const double kk = RANDOM.Uniform(-1,1);
            const double gthe = acos(kk);
            const TVector3 gdir(sin(gthe)*cos(gphi),sin(gthe)*sin(gphi),cos(gthe));
            const TVector3 unit_gdir_1 = gene*(gdir.Unit());
            const TVector3 unit_gdir_2 = -gene*(gdir.Unit());
             g1.SetPxPyPzE(unit_gdir_1.X(), unit_gdir_1.Y(), unit_gdir_1.Z(), gene);
             g2.SetPxPyPzE(unit_gdir_2.X(), unit_gdir_2.Y(), unit_gdir_2.Z(), gene);

            // boost to lab frame:
            const TVector3 boostq1=q1.BoostVector();
            g1.Boost(boostq1);
            g2.Boost(boostq1);
         
            // and then reconstruct stuff:
             const TLorentzVector pi0 = (g1+g2);

             DVCS_.m_pi0 = pi0.M();
             DVCS_.e_pi0 = pi0.E();
             DVCS_.phi_pi0 = pi0.Phi();
             DVCS_.the_pi0 = pi0.Theta();
             DVCS_.px_pi0 = pi0.Px();
             DVCS_.py_pi0 = pi0.Py();
             DVCS_.pz_pi0 = pi0.Pz();
             DVCS_.the_ggpi0 = g1.Vect().Angle(g2.Vect())*R2D;

             // constrain mass:   (this is insignificant)
             //vtmp=pi0S;
             //vtmp.SetE(sqrt(pow(pi0S.P(),2)+mq1*mq1));
             //const TLorentzVector pi0M=vtmp;

             const TLorentzVector X = k0 + target - k1 - p1 - pi0;
             const TLorentzVector pi0X = k0 + target - k1 - p1;
             const TLorentzVector p1X = k0 + target - k1 - pi0;

             DVCS_.m_kpi0X    = p1X.M();
             
             DVCS_.m2_kppi0X  = X.M2();
             DVCS_.e_kppi0X   = X.E();
             DVCS_.pt_kppi0X  = X.Pt();

             DVCS_.cop_pi0 = GEO.DihedralAngle(q0,p1,-pi0)*R2D;
             DVCS_.the_pi0_kpX = pi0.Vect().Angle(pi0X.Vect())*R2D;

             DVCS_.dphi_kpi0X_p = GEO.GetPhi(p1X.Phi()-p1.Phi(),-180/R2D)*R2D;
             DVCS_.the_kpi0X_p = p1X.Vect().Angle(p1.Vect())*R2D;
            
        }

       

        ///////// FILL NTUPLE //////////////////////////////////////////////////////////
        // Reconstruct all quantities with gamma #1.
        // Anything requiring 2 gammas was done above.

        const TLorentzVector X=(k0+target)-(k1+p1+g);
        const TLorentzVector gammaX = (k0+target)-(k1+p1);
        const TLorentzVector p1X = (k0+target)-(k1+g);

        const double Q2=-q0.M2();
        const double Wp=(q0+p0_proton).M();
        const double Nu=q0.E();
        const double zS=g.E()/q0.E();
        const double yS=q0.E()/k0.E();
        const double xbS=Q2/2/MPROTON/Nu;
        const double the_q0g=g.Vect().Angle(q0.Vect())*R2D;
        const double thegx=g.Vect().Angle(gammaX.Vect())*R2D;
        const double cop=GEO.DihedralAngle(q0,p1,-g)*R2D;
        //const double alpha=Nu-sqrt(Nu*Nu+Q2)*cos(the_q0g/R2D);
        //const double dphi = GEO.GetPhi(p1X.Phi()-p1.Phi(),-180/R2D)*R2D;
        //const double tcS= (Q2+2*Nu*alpha)/(1+alpha/target.M());

        //const double phi_h = 360-GEO.GetPhi(GEO.DihedralAngle(k1,q0,g),0)*R2D; 
        const double phi_h =   GEO.Moh_PHI_H(k0,k1,p1);;
        //const double phi2 =360-GEO.GetPhi(GEO.DihedralAngle(k1,q0,g),0)*R2D;
        //const double phi3 =180-GEO.GetPhi(GEO.DihedralAngle(k1,q0,p1),-180/R2D)*R2D;

        DVCS_.hel=helicity;

        DVCS_.phi_h=phi_h;
        //DVCS_.phi2=phi3;
        DVCS_.wp=Wp;
        DVCS_.Q2=Q2;
        DVCS_.xb=xbS;
        DVCS_.y_q0k0=yS;
        DVCS_.z_gq0=zS;
        DVCS_.nu=q0.E();
        DVCS_.t=-(target-p1).M2();
        //DVCS_.tc=tcS;
        DVCS_.cop=cop;
        DVCS_.the_q0g=the_q0g;
        DVCS_.thegx=thegx;
        DVCS_.therx=p1.Vect().Angle(p1X.Vect())*R2D;

        // missing quantities if we detect all the final state particles
        DVCS_.m2_kpgX =X.M2();
        DVCS_.e_kpgX  =X.E();
        DVCS_.p_kpgX  =X.P();
        DVCS_.pt_kpgX =X.Pt();
        DVCS_.px_kpgX =X.Px();
        DVCS_.py_kpgX =X.Py();
        DVCS_.pz_kpgX =X.Pz();
        
        DVCS_.m2_kpX=gammaX.M2(); // missing mass squared if we miss the real photon in the final state
        DVCS_.m_kgX=p1X.M();      // missing mass if we miss the recoild proton/4He


        DVCS_.v_x=v_x;
        DVCS_.v_y=v_y;
        DVCS_.v_z=v_z;
        
        // beam energy:
        DVCS_.e_k0=k0.E();
        DVCS_.p_p0=p0.P();


        // scattered electron:
        DVCS_.e_k1=k1.E();
        DVCS_.px_k1=k1.Px();
        DVCS_.py_k1=k1.Py();
        DVCS_.pz_k1=k1.Pz();
        DVCS_.phi_k1 = k1.Phi()*R2D;
        DVCS_.the_k1 = k1.Theta()*R2D;

        // virtual gamma:
        DVCS_.e_q0=q0.E();
        DVCS_.p_q0=q0.P();
        DVCS_.px_q0=q0.Px();
        DVCS_.py_q0=q0.Py();
        DVCS_.pz_q0=q0.Pz();         

        // recoil:
        DVCS_.e_p1=p1.P();
        DVCS_.m_p1 = p1.M();
        DVCS_.px_p1=p1.Px();
        DVCS_.py_p1=p1.Py();
        DVCS_.pz_p1=p1.Pz();
        DVCS_.phi_p1 = p1.Phi()*R2D;
        DVCS_.the_p1 = p1.Theta()*R2D;
        if (DVCS_.ic==1)
         {
           DVCS_.e_g[0]=g.E();
           DVCS_.px_g[0]=g.Px();
           DVCS_.py_g[0]=g.Py();
           DVCS_.pz_g[0]=g.Pz();
           DVCS_.the_g[0]=g.Theta();
           DVCS_.phi_g[0]=g.Phi();
          }
        if(DVCS_.ic==2)
         {
            DVCS_.e_g[0]=g1.E();
            DVCS_.px_g[0]=g1.Px();
            DVCS_.py_g[0]=g1.Py();
            DVCS_.pz_g[0]=g1.Pz();
            DVCS_.the_g[0]=g1.Theta();
            DVCS_.phi_g[0]=g1.Phi();

            DVCS_.e_g[1]=g2.E();
            DVCS_.px_g[1]=g2.Px();
            DVCS_.py_g[1]=g2.Py();
            DVCS_.pz_g[1]=g2.Pz();
            DVCS_.the_g[1]=g2.Theta();
            DVCS_.phi_g[1]=g2.Phi();
     
          } 
        tree->Fill();
   

    ////////////////////////////////////////////////////////////// 
    ////          fill the .dat output file               ////////
    //////////////////////////////////////////////////////////////

    if (DVCS_.ic==1)
       { trk.N_particles = 3; }
     else if (DVCS_.ic==2)      
        { trk.N_particles = 4; }
   
     trk.parent_A = parent_A;
     trk.Q2       = DVCS_.Q2;
     trk.nu       = DVCS_.nu;

     // scattered electron
     trk.Type[0]        = 1;
     trk.Particle_id[0] = elc_id;
     trk.Px[0]          = k1.Px();
     trk.Py[0]          = k1.Py();
     trk.Pz[0]          = k1.Pz();
     trk.E[0]           = k1.E();
     trk.M[0]           = 0.000511;
     trk.X_vertex[0]    = v_x;
     trk.Y_vertex[0]    = v_y;
     trk.Z_vertex[0]    = v_z;

     // scattered target
     trk.Type[1]        = target_type;
     trk.Particle_id[1] = target_id;
     trk.Px[1]          = p1.Px();
     trk.Py[1]          = p1.Py();
     trk.Pz[1]          = p1.Pz();
     trk.E[1]           = p1.E();
     trk.M[1]           = sqrt(p1.E()*p1.E()-p1.P()*p1.P());
     trk.X_vertex[1]    = v_x;
     trk.Y_vertex[1]    = v_y;
     trk.Z_vertex[1]    = v_z;

     // photons
     if (DVCS_.ic==1) {
         trk.Type[2]        = 1;
         trk.Particle_id[2] = photon_id;
         trk.Px[2]          = g.Px();
         trk.Py[2]          = g.Py();
         trk.Pz[2]          = g.Pz();
         trk.E[2]           = g.E();
         trk.M[2]           = 0.0;
         trk.X_vertex[2]    = v_x;
         trk.Y_vertex[2]    = v_y;
         trk.Z_vertex[2]    = v_z;
        }
     else if (DVCS_.ic==2) {
         trk.Type[2]        = 1;
         trk.Particle_id[2] = photon_id;
         trk.Px[2]          = g1.Px();
         trk.Py[2]          = g1.Py();
         trk.Pz[2]          = g1.Pz();
         trk.E[2]           = g1.E();
         trk.M[2]           = 0.0;
         trk.X_vertex[2]    = v_x;
         trk.Y_vertex[2]    = v_y;
         trk.Z_vertex[2]    = v_z;

         trk.Type[3]        = 1;
         trk.Particle_id[3] = photon_id;
         trk.Px[3]          = g2.Px();
         trk.Py[3]          = g2.Py();
         trk.Pz[3]          = g2.Pz();
         trk.E[3]           = g2.E();
         trk.M[3]           = 0.0;
         trk.X_vertex[3]    = v_x;
         trk.Y_vertex[3]    = v_y;
         trk.Z_vertex[3]    = v_z;           
        }
       
        int n_dat_out =  (int) iev/4000;
        outfile_[n_dat_out]<<trk.N_particles<<"   "<<trk.parent_A<<"   "<<trk.Q2<<"   "<<trk.nu<<"   "<<1<<"   "<<1<<endl;
        for(int i=0; i<trk.N_particles; i++){
            outfile_[n_dat_out]<<trk.Type[i]<<"   "<<trk.Particle_id[i]<<"   "<<trk.Particle_id[i]<<"   "<<i+1<<"   "<<i+1<<endl; 
            outfile_[n_dat_out]<<trk.Px[i]<<"   "<< trk.Py[i]<<"   "<< trk.Pz[i]<<"   "<<trk.E[i]<<"   "<<trk.M[i]<<endl;
            outfile_[n_dat_out]<<trk.X_vertex[i]<<"   "<<trk.Y_vertex[i]<<"   "<<trk.Z_vertex[i]<<"   "<<0.0<<"   "<<0.0<<endl;;
          }

    iev++;
    }


    timer.Stop();
    tree->AutoSave();
    rfile->Close();
    ffile.close();
    for( int ii=0; ii<n_dat_files; ii++)
        outfile_[ii].close();

    cout<<endl<<endl;
    cout<<"-------------------------------------------------"<<endl;
    cout<<setiosflags(ios::fixed)<<setprecision(0)<<nevents/timer.RealTime()/1000*60<<"K  Events Per Minute"<<endl;
    cout<<resetiosflags(ios::floatfield)<<setprecision(8)<<"Maximum Cross Section Weight:    "<<max<<" / "<<xsecmax<<endl;
    cout<<"# IC/EC Errors:    "<<nICECerr<<endl;
    cout<<"-------------------------------------------------"<<endl<<endl;

    exit(1);
}

