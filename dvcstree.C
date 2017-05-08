#ifndef __DVCSTREEC__
#define __DVCSTREEC__
#include "TTree.h"
#include "dvcstree.h"
#define INVLD -99999.0

void ZeroDvcsTree(dvcstree_t& t)
{
    t.hel=INVLD;
    t.P0target = INVLD;
    t.wp=INVLD;
    t.Q2=INVLD;
    t.xb=INVLD;
    t.y_q0k0=INVLD;
    t.z_gq0=INVLD;
    t.nu=INVLD;
    t.t=INVLD;

    t.thegx=INVLD;
    t.therx=INVLD;
    t.cop=INVLD;
    t.phi_h=INVLD;
    //t.phi2=INVLD;
  
    t.ic=0;

    t.m2_kpgX     =INVLD;
    t.e_kpgX      =INVLD;
    t.p_kpgX      =INVLD;
    t.px_kpgX     =INVLD;
    t.py_kpgX     =INVLD;
    t.pz_kpgX     =INVLD;
    t.pt_kpgX     =INVLD;
    t.e_kpgX      =INVLD;
    t.m2_kpX      =INVLD;
    t.m_kgX       =INVLD;
    t.the_q0g     =INVLD;
 

    t.phi_pi0     =INVLD; 
    t.the_pi0     =INVLD;
    t.m_pi0       =INVLD; 
    t.e_pi0       =INVLD; 
    t.px_pi0      =INVLD; 
    t.py_pi0      =INVLD; 
    t.pz_pi0      =INVLD; 
    t.the_ggpi0   =INVLD;
    t.m_kpi0X     =INVLD;
    t.e_kppi0X    =INVLD; 
    t.m2_kppi0X   =INVLD; 
    t.the_kpi0X_p =INVLD;
    t.pt_kppi0X   =INVLD; 
    t.cop_pi0     =INVLD;
    t.the_pi0_kpX =INVLD;
    t.dphi_kpi0X_p=INVLD;

    t.e_k0=INVLD;
    t.e_k1=INVLD; 
    t.px_k1=INVLD; 
    t.py_k1=INVLD; 
    t.pz_k1=INVLD;
    t.phi_k1=INVLD;
    t.the_k1=INVLD;
 
    t.e_p1=INVLD; 
    t.m_p1=INVLD;
    t.px_p1=INVLD; 
    t.py_p1=INVLD; 
    t.pz_p1=INVLD;
    t.p_p1=INVLD;
    t.phi_p1=INVLD;
    t.the_p1=INVLD;

    t.p_p1L=INVLD;
    t.phi_p1L=INVLD;
    t.the_p1L=INVLD;

    t.p1_delta_p=INVLD;
    t.p1_delta_phi=INVLD;
    t.p1_delta_theta=INVLD;

 
    t.p_q0=INVLD;
    t.e_q0=INVLD; 
    t.px_q0=INVLD; 
    t.py_q0=INVLD; 
    t.pz_q0=INVLD;

    for( int i=0; i<2; i++)
     {
      t.e_g[i] = INVLD; 
      t.px_g[i] = INVLD;  
      t.py_g[i] = INVLD; 
      t.pz_g[i] = INVLD;
      t.the_g[i]= INVLD;
      t.phi_g[i]= INVLD;
     }

    
    t.v_x = INVLD;
    t.v_y = INVLD;
    t.v_z = INVLD;

    t.p_p0=INVLD;

    t.gammas.clear();
    t.gammasEC.clear();
    t.gammasIC.clear();

}
TTree *InitDvcsTree(const char* otname,dvcstree_t& t_,const char* title="")
{
    TTree *t=new TTree(otname,title);

    t->Branch("phi_h",&t_.phi_h);
    //t->Branch("phi2",&t_.phi2);
    //t->Branch("phig",&t_.phig);
    t->Branch("hel",&t_.hel);
    t->Branch("wp",&t_.wp);
    t->Branch("Q2",&t_.Q2);
    t->Branch("xb",&t_.xb);
    t->Branch("y_q0k0",&t_.y_q0k0);
    t->Branch("z_gq0",&t_.z_gq0);
    t->Branch("nu",&t_.nu);
    t->Branch("t",&t_.t);

    t->Branch("cop",&t_.cop);
    t->Branch("thegx",&t_.thegx);
    t->Branch("therx",&t_.therx);

    
    t->Branch("m2_kpgX",&t_.m2_kpgX);
    t->Branch("e_kpgX",&t_.e_kpgX);
    t->Branch("p_kpgX",&t_.p_kpgX);
    t->Branch("px_kpgX",&t_.px_kpgX);
    t->Branch("py_kpgX",&t_.py_kpgX);
    t->Branch("pz_kpgX",&t_.pz_kpgX);
    t->Branch("pt_kpgX",&t_.pt_kpgX);
  
    t->Branch("m_kgX",&t_.m_kgX);
    t->Branch("m2_kpX",&t_.m2_kpX);


    t->Branch("m_kpi0X",&t_.m_kpi0X);
    t->Branch("m2_kppi0X",&t_.m2_kppi0X);
    t->Branch("e_kppi0X",&t_.e_kppi0X);
    t->Branch("pt_kppi0X",&t_.pt_kppi0X);
    t->Branch("cop_pi0",&t_.cop_pi0);
    t->Branch("the_pi0_kpX",&t_.the_pi0_kpX);
    t->Branch("dphi_kpi0X_p",&t_.dphi_kpi0X_p);
    t->Branch("the_kpi0X_p",&t_.the_kpi0X_p);
  
    t->Branch("phi_pi0",&t_.phi_pi0);
    t->Branch("the_pi0",&t_.the_pi0);
    t->Branch("m_pi0",&t_.m_pi0);
    t->Branch("e_pi0",&t_.e_pi0);
    t->Branch("px_pi0",&t_.px_pi0);
    t->Branch("py_pi0",&t_.py_pi0);
    t->Branch("pz_pi0",&t_.pz_pi0);
    t->Branch("the_ggpi0",&t_.the_ggpi0);
    
    // incident electrons
    t->Branch("e_k0",&t_.e_k0);
    t->Branch("v_x",&t_.v_x);
    t->Branch("v_y",&t_.v_y);
    t->Branch("v_z",&t_.v_z);

    // scattered electron
    t->Branch("e_k1",&t_.e_k1);
    t->Branch("px_k1",&t_.px_k1);
    t->Branch("py_k1",&t_.py_k1);
    t->Branch("pz_k1",&t_.pz_k1);
    t->Branch("phi_k1",&t_.phi_k1);
    t->Branch("the_k1",&t_.the_k1);

    // virtual photon   
    t->Branch("e_q0",&t_.e_q0);
    t->Branch("p_q0",&t_.p_q0);
    t->Branch("px_q0",&t_.px_q0);
    t->Branch("py_q0",&t_.py_q0);
    t->Branch("pz_q0",&t_.pz_q0);
    
    // recoiled proton/helium
    t->Branch("e_p1",&t_.e_p1);
    t->Branch("m_p1",&t_.m_p1);
    t->Branch("px_p1",&t_.px_p1);
    t->Branch("py_p1",&t_.py_p1);
    t->Branch("pz_p1",&t_.pz_p1);
    t->Branch("p_p1",&t_.p_p1);
    t->Branch("phi_p1",&t_.phi_p1);
    t->Branch("the_p1",&t_.the_p1);

    t->Branch("p_p1L",&t_.p_p1L);
    t->Branch("phi_p1L",&t_.phi_p1L);
    t->Branch("the_p1L",&t_.the_p1L);
    t->Branch("p1_delta_p",&t_.p1_delta_p);
    t->Branch("p1_delta_phi",&t_.p1_delta_phi);
    t->Branch("p1_delta_theta",&t_.p1_delta_theta);

    // photons   
    t->Branch("ic",&t_.ic);
    t->Branch("e_g",t_.e_g,"e_g[ic]/D");
    t->Branch("px_g",t_.px_g,"px_g[ic]/D");
    t->Branch("py_g",t_.py_g,"py_g[ic]/D");
    t->Branch("pz_g",t_.pz_g,"pz_g[ic]/D");
    t->Branch("the_g",t_.the_g,"the_g[ic]/D");
    t->Branch("phi_g",t_.phi_g,"phi_g[ic]/D");

    t->Branch("p_p0",&t_.p_p0);
    
    return t;
}
#endif
