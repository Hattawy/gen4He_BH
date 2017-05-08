#ifndef SPECSRTPC__hh
#define SPECSRTPC__hh
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TMarker.h"
#include "TSystem.h"
#include "SpecsPDG.hh"
#include <iostream>
#include <fstream>
#include <vector>
//////////////////////////////////////////////////////////////////////////////////////
class SpecsRTPC {
private:
    SpecsPDG PDG;
    static const int NROWS=80;
    static const int NCOLS=40;
    static const int NCHANS=3200;
    static const int NCOLSPERBOARD=8;
    static const int NROWSPERBOARD=2;
    // from gem/rtpc.h :
    static const double DEL_A=4.45;
    static const double DEL_Z=5.0;
    static const double PAD_RAD=69.0;
    // solenoid field strength:
    static const double FIELD=6; // Tesla
    // drift region radii:
    static const double DRIFT_R1=3;
    static const double DRIFT_R2=6;
    // same as in gem/read_scale.c
    static const double MINGAIN=0.01;
    static const double MAXGAIN=100.0;
    // parameters to read from input file:
    vector <int> BADBOARDS;
    vector <int> BADPADS;
    float GAINS[NCHANS];
    float PEDS[NCHANS];
    bool READBADPADS;
    bool READGAINS;
    bool READPEDS;
public:
    SpecsRTPC(const char* gainfile=NULL,const char* badpadfile=NULL,const char* pedafile=NULL);
    virtual ~SpecsRTPC();
    
    virtual void   UnityGains();
    virtual void   ZeroPedestals();
    virtual void   ZeroBadPads();

    virtual void   MarkBadRowCol(const int color=1,const float size=0.5,const int style=21);
    virtual void   MarkPad(const int pad,const int color=1,const float size=0.5,const int style=21);

    virtual bool   CheckRowCol(TH1* h);
    virtual bool   CheckRowCol(const int row,const int col);
    virtual bool   Pad2RowColZPhi(const int chan,int* row,int* col,double* z,double* phi);
    virtual int    RowCol2Pad(const int row,const int col);
    virtual int    Pad2Row(const int chan);
    virtual int    Pad2Col(const int chan);
    virtual double Pad2Z(const int chan);
    virtual double Pad2Phi(const int chan);
    virtual double RowCol2Z(const int row,const int col);
    virtual double RowCol2Phi(const int row,const int col);
    
    virtual bool   ReadGains(const char* ifile,const bool ignorebad=0);
    virtual bool   ReadGains(TH2* h);
    virtual void   ReadGains(double *g);
    virtual bool   ReadGains();
    virtual bool   ReadBadPads(const char* ifile);
    virtual bool   ReadBadPads();
    virtual bool   ReadPedestals(const char* ifile);
    
    virtual bool   WriteGains(const char* ofile,const bool printrowcol=0);
    virtual bool   WriteGains(TH2 *h,const char* ofile,const bool printrowcol=0);
    virtual bool   WriteBadPads(const char* ofile);
    virtual bool   WriteBadPads(const char* ofile,vector<int> pads);
    virtual bool   WriteBadPads(const char* ofile,vector<int> rows,vector<int> cols);
    
    virtual float  GetPedestal(const int chan);
    virtual float  GetPedestal(const int row,const int col);
    virtual float  GetGain(const int chan);
    virtual float  GetGain(const int row,const int col);
    virtual bool   IsBadPad(const int chan);
    virtual bool   IsBadPad(const int row,const int col);
    
    virtual void   PrintChannels(const int rmin,const int rmax,const int cmin,const int cmax);
    virtual double PathLength(const double radcur,const double theta);
    virtual double RadOfCur(const double poverq,const double theta);
    virtual double BetheBlochFWHM(const double poverq,const double mass,const double charge,const double theta);
    
    virtual void   CopyGains(float* g);
    virtual TH2*   GetGainHist(const bool fill=1);
    virtual TH2*   MakeGainHist(double *g);
    virtual TH2*   MakeGainHist(const char* ifile,const bool ignorebad=0);
    virtual bool   PrintOutOfRange(const char* ofile,const double lo,const double hi,const bool printrowcol);
    virtual void   DrawBadPads(const char* badfile=NULL);
    
    virtual int    GetBoard(const int row,const int col);
    virtual int    GetBoard(const int pad);
    virtual bool   InBoard(const int board,const int row,const int col);
    virtual bool   InBoard(const int board,const int pad);
    virtual bool   InBadBoard(const int row,const int col);
    virtual bool   InBadBoard(const int pad);
    virtual void   AddBadBoards(const int type);
    virtual void   AddBadBoard(const int board);
    virtual void   WriteBadBoardPads();

    virtual inline int GetNCHANS(){return NCHANS;}
    virtual inline int GetNROWS(){return NROWS;}
    virtual inline int GetNCOLS(){return NCOLS;}

};
//////////////////////////////////////////////////////////////////////////////////////
SpecsRTPC::SpecsRTPC(const char* gainfile,const char* badpadfile,const char* pedafile)
{
    ReadGains(gainfile);
    ReadBadPads(badpadfile);
    ReadPedestals(pedafile);
}
SpecsRTPC::~SpecsRTPC(){}
//////////////////////////////////////////////////////////////////////////////////////
void SpecsRTPC::UnityGains()
{
    READGAINS=0;
    for (int ipad=0; ipad<NCHANS; ipad++) GAINS[ipad]=1;
}
void SpecsRTPC::ZeroPedestals()
{
    READPEDS=0;
    for (int ipad=0; ipad<NCHANS; ipad++) PEDS[ipad]=1;
}
void SpecsRTPC::ZeroBadPads()
{
    READBADPADS=0;
    BADPADS.clear();
}
void SpecsRTPC::ReadGains(double *g)
{
    for (int ii=0; ii<NCHANS; ii++) GAINS[ii]=g[ii];
    READGAINS=1;
}
bool SpecsRTPC::ReadGains(TH2 *h)
{
    UnityGains();
    if (!CheckRowCol((TH1*)h)) {
        cerr<<"ERROR  SpecsRTPC::ReadGains(TH2*):  Invalid Binning"<<endl;
        return READGAINS;
    }
    for (int irow=0; irow<NROWS; irow++)
    {
        for (int icol=0; icol<NCOLS; icol++)
        {
            const int ipad = RowCol2Pad(irow,icol);
            GAINS[ipad] = h->GetBinContent(irow+1,icol+1);
        }
    }
    READGAINS=1;
    return READGAINS;
}
bool SpecsRTPC::ReadPedestals(const char* ifile)
{
    ZeroPedestals();
    if (ifile==NULL) return READPEDS;
    FILE *f=fopen(ifile,"r");
    if (f==NULL) {
        cerr<<"Error:  SpecsRTPC::ReadPedestals:  Missing file: "<<ifile<<endl;
        return READPEDS;
    }
    bool error=0;
    int npads=0;
    int nbad=0;
    float ped;
    while (fscanf(f,"%f",&ped)==1)
    {
        // pedestal file is sorted by column/row, but PEDS array is ordered by channel#
        const int icol = npads % NCOLS;
        const int irow = npads / NCOLS;
        const int ipad = RowCol2Pad(irow,icol);
       
        if (npads >= NCHANS) {
            cerr<<"Error:  SpecsRTPC::ReadPedestals:  Too Many Pads."<<endl;
            error=1;
            break;
        }
        if (ped < 0) {
            ped=0;
            nbad++;
        }

        PEDS[ipad] = ped;

        if (feof(f)) break;
        if (ferror(f)) break;
        npads++;
    }
    if (npads<NCHANS-1) {
        cerr<<"Error:  SpecsRTPC::ReadPedestals:  Missing Pads."<<endl;
        error=1;
    }
    if (error) 
    {
        ZeroPedestals();
        cout<<"--> Setting all Pedestals to 0."<<endl;
    }
    else {
        READPEDS=1;
        cout<<"RTPC Pedestals Read Successfully from:  "<<ifile<<endl;
        if (nbad>0) cout<<"* "<<nbad<<" pads with Negative pedestals, reset to zero"<<endl;
    }
    return READPEDS;
}
bool SpecsRTPC::ReadGains(const char* ifile,const bool ignorebad)
{
    UnityGains();
    if (ifile==NULL) return READGAINS;
    FILE *f=fopen(ifile,"r");
    if (f==NULL) {
        cerr<<"Error:  SpecsRTPC::ReadGains:  Missing file: "<<ifile<<endl;
        return READGAINS;
    }
    int nbig=0;
    int nlittle=0;
    bool error=0;
    int npads=0;
    float gain;
    while (fscanf(f,"%f",&gain)==1)
    {
        // gains file is sorted by column/row, but GAINS array is ordered by channel#
        // this is the same as $CLAS_PACK/gem, and should have been changed.
        const int icol = npads % NCOLS;
        const int irow = npads / NCOLS;
        const int ipad = RowCol2Pad(irow,icol);
       
        if (npads >= NCHANS) {
            cerr<<"Error:  SpecsRTPC::ReadGains:  Too Many Pads."<<endl;
            error=1;
            break;
        }
        if (!ignorebad && gain < MINGAIN) {
            fprintf(stderr,"Error:  SpecsRTPC::ReadGains:  Gain Tooo Small  p%.4d r%.2d c%.2d %f\n",ipad,irow,icol,gain);
            nlittle++;
//            error=1;
//            break;
        }
        if (!ignorebad && gain > MAXGAIN) {
            fprintf(stderr,"Error:  SpecsRTPC::ReadGains:  Gain Tooo Big  p%.4d r%.2d c%.2d %f\n",ipad,irow,icol,gain);
            nbig++;
//            error=1;
//            break;
        }

        GAINS[ipad] = gain;

        if (feof(f)) break;
        if (ferror(f)) break;
        npads++;
    }
    if (npads<NCHANS-1) {
        cerr<<"Error:  SpecsRTPC::ReadGains:  Missing Pads."<<endl;
        error=1;
    }
    if (error) 
    {
        UnityGains();
        cout<<"--> Setting all Gains to 1."<<endl;
    }
    else {
        READGAINS=1;
        cout<<"RTPC Gains Read Successfully from:  "<<ifile<<endl;
        if (nlittle>0) cout<<"* "<<nlittle<<" pads with tiny gains"<<endl;
        if (nbig  >0) cout<<"* "<<nbig  <<" pads with huge gains"<<endl;
    }
    return READGAINS;
}
bool SpecsRTPC::WriteGains(TH2 *h,const char* ofile,const bool printrowcol)
{
    if (!ReadGains(h)) return 0;
    return WriteGains(ofile,printrowcol);
}
bool SpecsRTPC::WriteGains(const char* ofile,const bool printrowcol)
{
    FILE *f=NULL;
    if (ofile)
    {
        if (!gSystem->AccessPathName(ofile))
        {
            cerr<<"ERROR  SpecsRTPC:WriteGains  file alerady exists: "<<ofile<<endl;
            return 0;
        }
        f=fopen(ofile,"w");
    }
    for (int irow=0; irow<NROWS; irow++)
    {
        for (int icol=0; icol<NCOLS; icol++)
        {
            const int ipad = RowCol2Pad(irow,icol);
            if (f) 
            {
                if (printrowcol) fprintf(f,"%8d %8d %8d %.8f\n",irow,icol,ipad,GAINS[ipad]);
                else             fprintf(f,"%.8f\n",GAINS[ipad]);
            }
            else
            {
                if (printrowcol) printf("%8d %8d %8d %.8f\n",irow,icol,ipad,GAINS[ipad]);
                else             printf("%.8f\n",GAINS[ipad]);
            }
        }
    }
    if (f) fclose(f);
    return 1;
}
bool SpecsRTPC::PrintOutOfRange(const char* ofile,const double lo,const double hi,const bool printrowcol)
{
    FILE *f=NULL;
    if (ofile)
    {
        if (!gSystem->AccessPathName(ofile))
        {
            cerr<<"ERROR  SpecsRTPC:WriteGains  file alerady exists: "<<ofile<<endl;
            return 0;
        }
        f=fopen(ofile,"w");
    }
    for (int ipad=0; ipad<NCHANS; ipad++)
    {
        if (GAINS[ipad]>lo && GAINS[ipad]<hi) continue;
        const int irow = Pad2Row(ipad);
        const int icol = Pad2Col(ipad);
        if (f) 
        {
            if (printrowcol) fprintf(f,"%8d %8d %8d %.8f\n",irow,icol,ipad,GAINS[ipad]);
            else             fprintf(f,"%d\n",ipad);
        }
        else
        {
            if (printrowcol) printf("%8d %8d %8d %.8f\n",irow,icol,ipad,GAINS[ipad]);
            else             printf("%8d\n",ipad);
        }
    }
    if (f) { fprintf(f,"-1000\n"); fclose(f); }
    return 1;
}
bool SpecsRTPC::ReadBadPads(const char* badpadfile)
{
    READBADPADS=0;
    BADPADS.clear();
    if (badpadfile==NULL) return 0;
    FILE *f=fopen(badpadfile,"r");
    if (f==NULL) {
        cerr<<"Error:  SpecsRTPC::ReadBadPads:  Missing file: "<<badpadfile<<endl;
        return 0;
    }
    int pad;
    bool error=0;
    while (fscanf(f,"%d",&pad)==1)
    {
        if (pad<0) break;
        if (pad>NCHANS-1) {
            cerr<<"Error:  SpecsRTPC::ReadBadPads:  Invalid Pad# in Data File: "<<pad<<endl;
            error=1;
            break;
        }
        if (BADPADS.size() >= NCHANS-1) {
            cerr<<"Error:  SpecsRTPC::ReadBadPads:  Too Many Bad Pads."<<endl;
            error=1;
            break;
        }
        BADPADS.push_back(pad);
    }
    if (error) BADPADS.clear();
    else {
        READBADPADS=1;
        cout<<"RTPC "<<BADPADS.size()<<" BadPads Read Successfully from:  "<<badpadfile<<endl;
    }

    return READBADPADS;
}
bool SpecsRTPC::WriteBadPads(const char* ofile)
{
    return WriteBadPads(ofile,BADPADS);
}
bool SpecsRTPC::WriteBadPads(const char* ofile,vector<int> pads)
{
    if (!gSystem->AccessPathName(ofile))
    {
        cerr<<"File Already Exists: "<<ofile<<endl;
        return 0;
    }
    ofstream o;
    o.open(ofile);
    for (unsigned int ii=0; ii<pads.size(); ii++) o<<pads[ii]<<endl;
    o<<"-1000"<<endl;
    o.close();
    return 1;
}
bool SpecsRTPC::WriteBadPads(const char* ofile,vector<int> rows,vector<int> cols)
{
    if (rows.size()!=cols.size())
    {
        cerr<<"Error: SpecsRTPC::WriteBadPads:  Vectors of Different lengths."<<endl;
        return 0;
    }
    if (!gSystem->AccessPathName(ofile))
    {
        cerr<<"File Already Exists: "<<ofile<<endl;
        return 0;
    }
    ofstream o;
    o.open(ofile);
    for (unsigned int ii=0; ii<rows.size(); ii++) o<<rows[ii]<<" "<<cols[ii]<<endl;
    o<<"-1000"<<endl;
    o.close();
    return 1;
}
//////////////////////////////////////////////////////////////////////////////////////
bool SpecsRTPC::ReadGains() { return READGAINS; }
bool SpecsRTPC::ReadBadPads() { return READBADPADS; }
//////////////////////////////////////////////////////////////////////////////////////
float SpecsRTPC::GetGain(const int chan)
{
    if (chan<0 || chan>=NCHANS) {
        cerr<<"Error:  SpecsRTPC::GetGain:  Invalid Pad#: "<<chan<<endl;
        return 1;
    }
    return GAINS[chan];
}
float SpecsRTPC::GetGain(const int row,const int col)
{
    if (!CheckRowCol(row,col)) return 1;
    return GetGain(RowCol2Pad(row,col));
}
bool SpecsRTPC::IsBadPad(const int chan)
{
    if (chan<0 || chan>=NCHANS) {
        cerr<<"Error:  SpecsRTPC::IsBadPad:  Invalid Pad#: "<<chan<<endl;
        return 1;
    }
    for (unsigned int ipad=0; ipad<BADPADS.size(); ipad++)
        if (chan==BADPADS[ipad]) return 1;
    return 0;
}
bool SpecsRTPC::IsBadPad(const int row,const int col)
{
    if (!CheckRowCol(row,col)) return 1;
    return IsBadPad(RowCol2Pad(row,col));
}
float SpecsRTPC::GetPedestal(const int chan)
{
    if (chan<0 || chan>=NCHANS) {
        cerr<<"Error:  SpecsRTPC::GetPedestal:  Invalid Pad#: "<<chan<<endl;
        return 1;
    }
    return PEDS[chan];
}
float SpecsRTPC::GetPedestal(const int row,const int col)
{
    if (!CheckRowCol(row,col)) return 1;
    return GetPedestal(RowCol2Pad(row,col));
}
//////////////////////////////////////////////////////////////////////////////////////
bool SpecsRTPC::Pad2RowColZPhi(const int chan,int* row,int* col,double* zz,double* pp)
{   // given pad#, get row, column, z, phi
    if (chan<0 || chan>=NCHANS) {
        cerr<<"Error:  SpecsRTPC::Chan2RowColZPhi:  Invalid chan: "<<chan<<endl;
        return 0;
    }
    // from gem/rtpc.h:
    static const double phi_start[2]={5.02700,1.88542};
    static const double zoff[4]={-97.50,-96.50,-97.50,-99.50};
    static const int seq_col[2][16] = {
        {4, 5, 2, 3, 7, 7, 1, 0, 1, 3, 4, 6, 2, 0, 6, 5},
        {1, 3, 4, 6, 2, 0, 6, 5, 4, 5, 2, 3, 7, 7, 1, 0} };
    static const int seq_row[2][16] = {
        {0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0},
        {0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0} };
    // from gem/generate_pad_locations.c and gem/pad2rowcol.c :
    const int rtpc_side=chan/1600;
    const int mypad=chan%1600;
    const int conn_num=mypad/16;
    const int conn_row=conn_num/5;
    const int conn_col=conn_num%5;
    const int pad_seq=mypad%16;
    const int odd=conn_num%2;
    *col=8*conn_col+seq_col[odd][pad_seq];
    *row=2*conn_row+seq_row[odd][pad_seq]+40*rtpc_side;
    *zz=zoff[(*row-40*rtpc_side)%4]+DEL_Z*(*col);
    *pp=phi_start[rtpc_side]+(DEL_A/PAD_RAD)*(*row-40*rtpc_side);
    return 1;
}
int SpecsRTPC::RowCol2Pad(const int row,const int col)
{   // given row/column, return pad#
    if (row<0 || row>=NROWS || col<0 || col>= NCOLS) {
        cerr<<"Error:  SpecsRTPC::RowCol2Chan:   Invalid row,col: "<<row<<" "<<col<<endl;
        return -9999;
    }
    // from gem/rtpc.h:
    static const int pcb_map[2][2][8]= {
        {{ 7,  8,  12, 3,  0,  15, 11, 4  }, 
         { 13, 6,  2,  9,  10, 1,  14, 5  }}, 
        {{ 15, 0,  4,  11, 8,  7,  3,  12 },
         { 5,  14, 10, 1,  2,  9,  6,  13 }} };
    // from gem/rowcol2pad.c:
    const int myrow=row%40;
    const int rtpc_side=row/40;
    const int conn_row=myrow/2;
    const int conn_side=myrow%2;
    const int conn_col=col%8;
    const int conn_num=5*conn_row+col/8;
    const int even_odd=conn_num%2;
    return 16*conn_num+pcb_map[even_odd][conn_side][conn_col]+1600*rtpc_side;
}
//////////////////////////////////////////////////////////////////////////////////////
int SpecsRTPC::Pad2Row(const int chan)
{   // given pad channel, return row 
    int col,row;
    double zzz,ppp;
    if (Pad2RowColZPhi(chan,&row,&col,&zzz,&ppp)) return row;
    else return -9999;
}
int SpecsRTPC::Pad2Col(const int chan)
{   // given pad channel, return column
    int col,row;
    double zzz,ppp;
    if (Pad2RowColZPhi(chan,&row,&col,&zzz,&ppp)) return col;
    else return -9999;
}
double SpecsRTPC::Pad2Z(const int chan)
{
    int col,row;
    double zzz,ppp;
    if (Pad2RowColZPhi(chan,&row,&col,&zzz,&ppp)) return zzz;
    else return -9999;
}
double SpecsRTPC::Pad2Phi(const int chan)
{
    int col,row;
    double zzz,ppp;
    if (Pad2RowColZPhi(chan,&row,&col,&zzz,&ppp)) return ppp;
    else return -9999;
}
double SpecsRTPC::RowCol2Z(const int row,const int col)
{
    return Pad2Z(RowCol2Pad(row,col));
}
double SpecsRTPC::RowCol2Phi(const int row,const int col)
{
    return Pad2Phi(RowCol2Pad(row,col));
}
//////////////////////////////////////////////////////////////////////////////////////
bool SpecsRTPC::CheckRowCol(TH1 *h)
{
    if (h->GetNbinsX()==NROWS && h->GetNbinsY()==NCOLS) return 1;
    cerr<<"Error CheckRowCol "<<h->GetName()<<endl;
    return 0;
}
bool SpecsRTPC::CheckRowCol(const int row,const int col)
{
    if (row>=0 && row<NROWS && col>=0 && col<NCOLS) return 1;
    cerr<<"Error CheckRowCol: "<<row<<" "<<col<<endl;
    return 0;
}
void SpecsRTPC::MarkBadRowCol(const int color,const float size,const int style)
{
    for (int ix=0; ix<NROWS; ix++) {
        for (int iy=0; iy<NCOLS; iy++) {
            if (IsBadPad(ix,iy)) {
                TMarker m((float)ix,(float)iy,style);
                m.SetMarkerColor(color);
                m.SetMarkerSize(size);
                m.DrawClone();
            }
        }
    }
}
void SpecsRTPC::MarkPad(const int pad,const int color,const float size,const int style)
{
    const int row=Pad2Row(pad);
    const int col=Pad2Col(pad);
    TMarker m((float)row,(float)col,style);
    m.SetMarkerColor(color);
    m.SetMarkerSize(size);
    m.DrawClone();
}
void SpecsRTPC::PrintChannels(const int rmin,const int rmax,const int cmin,const int cmax)
{
    for (int row=rmin; row<=rmax; row++)
        for (int col=cmin; col<=cmax; col++)
            cout<<RowCol2Pad(row,col)<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////
double SpecsRTPC::RadOfCur(const double poverq,const double theta)
{
    // poverq = drift region momentum p/q (MeV)
    // theta  = track polar angle (rad)
    const double poverq_t = poverq * sin(theta);
    return poverq_t / FIELD / 3; // (cm)
}
double SpecsRTPC::PathLength(const double radcur,const double theta)
{
    // radcur = radius of curvature (cm)
    // theta  = track polar angle (rad)
    if (fabs(radcur-(DRIFT_R2-DRIFT_R1)/2) < 1e-2) return 0;
    if (radcur <= DRIFT_R1/2) return 0;
    double phi1,phi2;
    if (radcur <= DRIFT_R2/2) {
        phi1 = acos((radcur*radcur+DRIFT_R2*DRIFT_R2-DRIFT_R1*DRIFT_R1)/2/DRIFT_R2/radcur);
        phi2 = 0;
    } else {
        phi1 = asin( DRIFT_R1/radcur * sqrt(1-pow(DRIFT_R1/radcur/2,2)) );
        phi2 = asin( DRIFT_R2/radcur * sqrt(1-pow(DRIFT_R2/radcur/2,2)) );
        if (radcur<DRIFT_R2/sqrt(2)) {
            phi2=TMath::Pi()-phi2;
        }
    }
    const double path = radcur * fabs(phi2-phi1) / sin(theta);
    return path >= DRIFT_R2-DRIFT_R1 ? path : 0;
}
double SpecsRTPC::BetheBlochFWHM(const double poverq,const double mass,const double charge,const double theta)
{
    // poverq = drift region momentum p/q (MeV)
    // theta  = track polar angle (rad)
    // mass   = (MeV)
    const double mom = poverq*charge;
    const double beta = mom / sqrt( mom*mom + mass*mass );
    const double radcur = RadOfCur(poverq,theta);
    const double pathlength = PathLength(radcur,theta);
    return pathlength>0 ? PDG.BetheBlochFWHM(beta,pathlength,charge) : 0;
}
void SpecsRTPC::CopyGains(float *g) { for (int ii=0; ii<NCHANS; ii++) g[ii]=GAINS[ii]; }

TH2* SpecsRTPC::MakeGainHist(const char* ifile,const bool ignorebad)
{
    ReadGains(ifile,ignorebad);
    return GetGainHist(1);
}
TH2* SpecsRTPC::MakeGainHist(double *g)
{
    ReadGains(g);
    return GetGainHist(1);
}
    
TH2* SpecsRTPC::GetGainHist(const bool fill)
{
    static int ncalls=0;
    TH2D *h=new TH2D(Form("hgains__n%d",ncalls++),"",
            NROWS, -0.5,NROWS-0.5,
            NCOLS, -0.5,NCOLS-0.5);
    h->GetXaxis()->SetTitle("Row");
    h->GetYaxis()->SetTitle("Column");
    h->GetXaxis()->SetTickLength(-0.01);
    h->GetYaxis()->SetTickLength(-0.005);
    h->GetXaxis()->SetNdivisions(1010);
    h->GetYaxis()->SetTitleOffset(0.7);
    if (!fill) return h;
    for (int irow=0; irow<NROWS; irow++)
    {
        for (int icol=0; icol<NCOLS; icol++)
        {
            const int ipad=RowCol2Pad(irow,icol);
            h->SetBinContent(irow+1,icol+1,GAINS[ipad]);
        }
    }
    return h;
}
void SpecsRTPC::DrawBadPads(const char* badfile)
{
    if (badfile) ReadBadPads(badfile);
    TH2 *h=GetGainHist(0);
    h->Draw();
    MarkBadRowCol();
}
int SpecsRTPC::GetBoard(const int row,const int col)
{
    if (!CheckRowCol(row,col)) return -1;
    //return (NCOLS/NCOLSPERBOARD)*(row/NROWSPERBOARD) + col/NCOLSPERBOARD;
    return 5*(row/2) + col/8;
}
int SpecsRTPC::GetBoard(const int pad)
{
    const int row=Pad2Row(pad);
    const int col=Pad2Col(pad);
    return GetBoard(row,col);
}
bool SpecsRTPC::InBoard(const int board,const int row,const int col)
{
    return board==GetBoard(row,col);
}
bool SpecsRTPC::InBoard(const int board,const int pad)
{
    return board==GetBoard(pad);
}
void SpecsRTPC::AddBadBoard(const int board)
{
    if (board<0 || board>=200)
        cerr<<"Error:  SpecsRTPC::AddBadBoard:  Invalid Board#: "<<board<<endl;
    else
    {
        cout<<"SpecsRTPC:  Added Board:  "<<board<<endl;
        BADBOARDS.push_back(board);
    }
}
bool SpecsRTPC::InBadBoard(const int row,const int col)
{
    if (!CheckRowCol(row,col)) return 1;
    for (unsigned int ii=0; ii<BADBOARDS.size(); ii++)
        if (InBoard(BADBOARDS[ii],row,col)) return 1;
    return 0;
}
bool SpecsRTPC::InBadBoard(const int pad)
{
    return InBadBoard(Pad2Row(pad),Pad2Col(pad));
}
void SpecsRTPC::WriteBadBoardPads()
{
    for (int ii=0; ii<NCHANS; ii++) if (InBadBoard(ii)) cout<<ii<<endl;
}
void SpecsRTPC::AddBadBoards(const int type)
{
    if (type >= 0)
    {
        AddBadBoard(181);
    }
    if (type > 0)
    {
        AddBadBoard(157);
        AddBadBoard(174);
    }
    if (type > 1)
    {
        AddBadBoard(25);
        AddBadBoard(173);
        AddBadBoard(73);
        AddBadBoard(158);
        AddBadBoard(180);
        AddBadBoard(166);
    }
    if (type > 2)
    {
        AddBadBoard(49);
        AddBadBoard(108);
        AddBadBoard(132);
        AddBadBoard(141);
        AddBadBoard(148);
        AddBadBoard(150);
        AddBadBoard(191);
    }
    cout<<"SpecsRTPC::AddBadBoards  -- "<<BADBOARDS.size()<<" Total Bad Boards"<<endl;
}




#endif
