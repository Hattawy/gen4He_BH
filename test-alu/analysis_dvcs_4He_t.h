//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 25 11:07:09 2016 by ROOT version 6.06/02
// from TChain ch/analysis_dvcs_4He_t
//////////////////////////////////////////////////////////

#ifndef analysis_dvcs_4He_t_h
#define analysis_dvcs_4He_t_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class analysis_dvcs_4He_t {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        phi_h;
   Int_t           hel;
   Double_t        wp;
   Double_t        Q2;
   Double_t        xb;
   Double_t        y_q0k0;
   Double_t        z_gq0;
   Double_t        nu;
   Double_t        t;
   Double_t        cop;
   Double_t        thegx;
   Double_t        therx;
   Double_t        m2_kpgX;
   Double_t        e_kpgX;
   Double_t        p_kpgX;
   Double_t        px_kpgX;
   Double_t        py_kpgX;
   Double_t        pz_kpgX;
   Double_t        pt_kpgX;
   Double_t        m_kgX;
   Double_t        m2_kpX;
   Double_t        m_kpi0X;
   Double_t        m2_kppi0X;
   Double_t        e_kppi0X;
   Double_t        pt_kppi0X;
   Double_t        cop_pi0;
   Double_t        the_pi0_kpX;
   Double_t        dphi_kpi0X_p;
   Double_t        the_kpi0X_p;
   Double_t        phi_pi0;
   Double_t        the_pi0;
   Double_t        m_pi0;
   Double_t        e_pi0;
   Double_t        px_pi0;
   Double_t        py_pi0;
   Double_t        pz_pi0;
   Double_t        the_ggpi0;
   Double_t        e_k0;
   Double_t        v_x;
   Double_t        v_y;
   Double_t        v_z;
   Double_t        e_k1;
   Double_t        px_k1;
   Double_t        py_k1;
   Double_t        pz_k1;
   Double_t        phi_k1;
   Double_t        the_k1;
   Double_t        e_q0;
   Double_t        p_q0;
   Double_t        px_q0;
   Double_t        py_q0;
   Double_t        pz_q0;
   Double_t        e_p1;
   Double_t        m_p1;
   Double_t        px_p1;
   Double_t        py_p1;
   Double_t        pz_p1;
   Double_t        p_p1;
   Double_t        phi_p1;
   Double_t        the_p1;
   Double_t        p_p1L;
   Double_t        phi_p1L;
   Double_t        the_p1L;
   Double_t        p1_delta_p;
   Double_t        p1_delta_phi;
   Double_t        p1_delta_theta;
   Int_t           ic;
   Double_t        e_g[1];   //[ic]
   Double_t        px_g[1];   //[ic]
   Double_t        py_g[1];   //[ic]
   Double_t        pz_g[1];   //[ic]
   Double_t        the_g[1];   //[ic]
   Double_t        phi_g[1];   //[ic]
   Double_t        p_p0;

   // List of branches
   TBranch        *b_phi_h;   //!
   TBranch        *b_hel;   //!
   TBranch        *b_wp;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_xb;   //!
   TBranch        *b_y_q0k0;   //!
   TBranch        *b_z_gq0;   //!
   TBranch        *b_nu;   //!
   TBranch        *b_t;   //!
   TBranch        *b_cop;   //!
   TBranch        *b_thegx;   //!
   TBranch        *b_therx;   //!
   TBranch        *b_m2_kpgX;   //!
   TBranch        *b_e_kpgX;   //!
   TBranch        *b_p_kpgX;   //!
   TBranch        *b_px_kpgX;   //!
   TBranch        *b_py_kpgX;   //!
   TBranch        *b_pz_kpgX;   //!
   TBranch        *b_pt_kpgX;   //!
   TBranch        *b_m_kgX;   //!
   TBranch        *b_m2_kpX;   //!
   TBranch        *b_m_kpi0X;   //!
   TBranch        *b_m2_kppi0X;   //!
   TBranch        *b_e_kppi0X;   //!
   TBranch        *b_pt_kppi0X;   //!
   TBranch        *b_cop_pi0;   //!
   TBranch        *b_the_pi0_kpX;   //!
   TBranch        *b_dphi_kpi0X_p;   //!
   TBranch        *b_the_kpi0X_p;   //!
   TBranch        *b_phi_pi0;   //!
   TBranch        *b_the_pi0;   //!
   TBranch        *b_m_pi0;   //!
   TBranch        *b_e_pi0;   //!
   TBranch        *b_px_pi0;   //!
   TBranch        *b_py_pi0;   //!
   TBranch        *b_pz_pi0;   //!
   TBranch        *b_the_ggpi0;   //!
   TBranch        *b_e_k0;   //!
   TBranch        *b_v_x;   //!
   TBranch        *b_v_y;   //!
   TBranch        *b_v_z;   //!
   TBranch        *b_e_k1;   //!
   TBranch        *b_px_k1;   //!
   TBranch        *b_py_k1;   //!
   TBranch        *b_pz_k1;   //!
   TBranch        *b_phi_k1;   //!
   TBranch        *b_the_k1;   //!
   TBranch        *b_e_q0;   //!
   TBranch        *b_p_q0;   //!
   TBranch        *b_px_q0;   //!
   TBranch        *b_py_q0;   //!
   TBranch        *b_pz_q0;   //!
   TBranch        *b_e_p1;   //!
   TBranch        *b_m_p1;   //!
   TBranch        *b_px_p1;   //!
   TBranch        *b_py_p1;   //!
   TBranch        *b_pz_p1;   //!
   TBranch        *b_p_p1;   //!
   TBranch        *b_phi_p1;   //!
   TBranch        *b_the_p1;   //!
   TBranch        *b_p_p1L;   //!
   TBranch        *b_phi_p1L;   //!
   TBranch        *b_the_p1L;   //!
   TBranch        *b_p1_delta_p;   //!
   TBranch        *b_p1_delta_phi;   //!
   TBranch        *b_p1_delta_theta;   //!
   TBranch        *b_ic;   //!
   TBranch        *b_e_g;   //!
   TBranch        *b_px_g;   //!
   TBranch        *b_py_g;   //!
   TBranch        *b_pz_g;   //!
   TBranch        *b_the_g;   //!
   TBranch        *b_phi_g;   //!
   TBranch        *b_p_p0;   //!

   analysis_dvcs_4He_t(TTree *tree=0);
   virtual ~analysis_dvcs_4He_t();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef analysis_dvcs_4He_t_cxx
analysis_dvcs_4He_t::analysis_dvcs_4He_t(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("ch",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("ch","analysis_dvcs_4He_t");
      chain->Add("../gen4He_cohdvcs_R.root/Tgen4He");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

analysis_dvcs_4He_t::~analysis_dvcs_4He_t()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t analysis_dvcs_4He_t::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t analysis_dvcs_4He_t::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void analysis_dvcs_4He_t::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("phi_h", &phi_h, &b_phi_h);
   fChain->SetBranchAddress("hel", &hel, &b_hel);
   fChain->SetBranchAddress("wp", &wp, &b_wp);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("xb", &xb, &b_xb);
   fChain->SetBranchAddress("y_q0k0", &y_q0k0, &b_y_q0k0);
   fChain->SetBranchAddress("z_gq0", &z_gq0, &b_z_gq0);
   fChain->SetBranchAddress("nu", &nu, &b_nu);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("cop", &cop, &b_cop);
   fChain->SetBranchAddress("thegx", &thegx, &b_thegx);
   fChain->SetBranchAddress("therx", &therx, &b_therx);
   fChain->SetBranchAddress("m2_kpgX", &m2_kpgX, &b_m2_kpgX);
   fChain->SetBranchAddress("e_kpgX", &e_kpgX, &b_e_kpgX);
   fChain->SetBranchAddress("p_kpgX", &p_kpgX, &b_p_kpgX);
   fChain->SetBranchAddress("px_kpgX", &px_kpgX, &b_px_kpgX);
   fChain->SetBranchAddress("py_kpgX", &py_kpgX, &b_py_kpgX);
   fChain->SetBranchAddress("pz_kpgX", &pz_kpgX, &b_pz_kpgX);
   fChain->SetBranchAddress("pt_kpgX", &pt_kpgX, &b_pt_kpgX);
   fChain->SetBranchAddress("m_kgX", &m_kgX, &b_m_kgX);
   fChain->SetBranchAddress("m2_kpX", &m2_kpX, &b_m2_kpX);
   fChain->SetBranchAddress("m_kpi0X", &m_kpi0X, &b_m_kpi0X);
   fChain->SetBranchAddress("m2_kppi0X", &m2_kppi0X, &b_m2_kppi0X);
   fChain->SetBranchAddress("e_kppi0X", &e_kppi0X, &b_e_kppi0X);
   fChain->SetBranchAddress("pt_kppi0X", &pt_kppi0X, &b_pt_kppi0X);
   fChain->SetBranchAddress("cop_pi0", &cop_pi0, &b_cop_pi0);
   fChain->SetBranchAddress("the_pi0_kpX", &the_pi0_kpX, &b_the_pi0_kpX);
   fChain->SetBranchAddress("dphi_kpi0X_p", &dphi_kpi0X_p, &b_dphi_kpi0X_p);
   fChain->SetBranchAddress("the_kpi0X_p", &the_kpi0X_p, &b_the_kpi0X_p);
   fChain->SetBranchAddress("phi_pi0", &phi_pi0, &b_phi_pi0);
   fChain->SetBranchAddress("the_pi0", &the_pi0, &b_the_pi0);
   fChain->SetBranchAddress("m_pi0", &m_pi0, &b_m_pi0);
   fChain->SetBranchAddress("e_pi0", &e_pi0, &b_e_pi0);
   fChain->SetBranchAddress("px_pi0", &px_pi0, &b_px_pi0);
   fChain->SetBranchAddress("py_pi0", &py_pi0, &b_py_pi0);
   fChain->SetBranchAddress("pz_pi0", &pz_pi0, &b_pz_pi0);
   fChain->SetBranchAddress("the_ggpi0", &the_ggpi0, &b_the_ggpi0);
   fChain->SetBranchAddress("e_k0", &e_k0, &b_e_k0);
   fChain->SetBranchAddress("v_x", &v_x, &b_v_x);
   fChain->SetBranchAddress("v_y", &v_y, &b_v_y);
   fChain->SetBranchAddress("v_z", &v_z, &b_v_z);
   fChain->SetBranchAddress("e_k1", &e_k1, &b_e_k1);
   fChain->SetBranchAddress("px_k1", &px_k1, &b_px_k1);
   fChain->SetBranchAddress("py_k1", &py_k1, &b_py_k1);
   fChain->SetBranchAddress("pz_k1", &pz_k1, &b_pz_k1);
   fChain->SetBranchAddress("phi_k1", &phi_k1, &b_phi_k1);
   fChain->SetBranchAddress("the_k1", &the_k1, &b_the_k1);
   fChain->SetBranchAddress("e_q0", &e_q0, &b_e_q0);
   fChain->SetBranchAddress("p_q0", &p_q0, &b_p_q0);
   fChain->SetBranchAddress("px_q0", &px_q0, &b_px_q0);
   fChain->SetBranchAddress("py_q0", &py_q0, &b_py_q0);
   fChain->SetBranchAddress("pz_q0", &pz_q0, &b_pz_q0);
   fChain->SetBranchAddress("e_p1", &e_p1, &b_e_p1);
   fChain->SetBranchAddress("m_p1", &m_p1, &b_m_p1);
   fChain->SetBranchAddress("px_p1", &px_p1, &b_px_p1);
   fChain->SetBranchAddress("py_p1", &py_p1, &b_py_p1);
   fChain->SetBranchAddress("pz_p1", &pz_p1, &b_pz_p1);
   fChain->SetBranchAddress("p_p1", &p_p1, &b_p_p1);
   fChain->SetBranchAddress("phi_p1", &phi_p1, &b_phi_p1);
   fChain->SetBranchAddress("the_p1", &the_p1, &b_the_p1);
   fChain->SetBranchAddress("p_p1L", &p_p1L, &b_p_p1L);
   fChain->SetBranchAddress("phi_p1L", &phi_p1L, &b_phi_p1L);
   fChain->SetBranchAddress("the_p1L", &the_p1L, &b_the_p1L);
   fChain->SetBranchAddress("p1_delta_p", &p1_delta_p, &b_p1_delta_p);
   fChain->SetBranchAddress("p1_delta_phi", &p1_delta_phi, &b_p1_delta_phi);
   fChain->SetBranchAddress("p1_delta_theta", &p1_delta_theta, &b_p1_delta_theta);
   fChain->SetBranchAddress("ic", &ic, &b_ic);
   fChain->SetBranchAddress("e_g", e_g, &b_e_g);
   fChain->SetBranchAddress("px_g", px_g, &b_px_g);
   fChain->SetBranchAddress("py_g", py_g, &b_py_g);
   fChain->SetBranchAddress("pz_g", pz_g, &b_pz_g);
   fChain->SetBranchAddress("the_g", the_g, &b_the_g);
   fChain->SetBranchAddress("phi_g", phi_g, &b_phi_g);
   fChain->SetBranchAddress("p_p0", &p_p0, &b_p_p0);
   Notify();
}

Bool_t analysis_dvcs_4He_t::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void analysis_dvcs_4He_t::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t analysis_dvcs_4He_t::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef analysis_dvcs_4He_t_cxx
