#ifndef __DVCSTREE__H
#define __DVCSTREE__H
#include <vector>
#include "TLorentzVector.h"

/*
   NTUPLE VARIABLES DEFINITION

   Naming conventions are almost always true.

   I use this same ntuple format for data and simulation.
   Some variables are only used in one or the other.

   

   For pi0 case:
   First gamma (g1) is always detected.
   If both gammas are detected, then g1
   is the most energetic of the two.

   
   
   k1 =  scattered electron
   q0 =  virtual photon
   g1 =  real photon (may also be called q1)
   g2 =  real photon
   p1 =  4He/proton

   // common prefixes:
   cx,cy,cz = direction cosines
   e = energy (except in case of recoil, then it's momentum!)
   
   // common missing prefixes:
   mX  = missing mass
   m2X = missing mass squared
   eX  = missing energy
   pX  = missing momentum
   ptX = missing transverse momentum


   // common suffixes:
   x   = reconstructed from other particles
   g   = generated value
   pi0 = this quantity requires two gammas
 

   // 4-momentum transfer:
   t0 = generated
   tg = (q0-g1)^2
   tr = (p0-p1)^2
   tc = constrained t



   // DVCS Exclusivity variables:
   dphi  = azimuthal angle between p1 and p1x(k1+g1)
   thegx = angle between g1 and g1x(k1+p1)
   therx = angle between p1 and p1x(k1+g1)
   cop   = q0/p1/g1 coplanarity
   pX    = Px(k1+g1+p1)
   eX    = Ex(k1+g1+p1)
   m2X   = Mx(k1+g1+p1)
   mXeg  = Mx(k1+g1)
   m2Xer = Mx^2(k1+p1)

   // pi0 Exclusivity variables:
   dphi_kpi0X_p  = azimuthal angle between p1 and p1x(k1+g1+g2)
   cop_pi0   = q0/p1/g1+g2 coplanarity
   m_kpi0X   = Mx(k1+g1+g2)
   e_kppi0X  = Ex(k1+g1+g2+p1)
   m2_kppi0X = Mx^2(k1+g1+g2+p1)
   the_kpi0X_p = angle between p1 and p1x(k1+g1+g2)
   the_pi0_kpX  = angle between g1g2 and g1g2x(k1+p1)

   the_ggpi0 = angle between q0 and q1
   
   phi = DVCS phi

   vz = z-vertex position (CLAS coordinate system)

   pf = p0 momentum (0 for coherent, fermi motion for incoherent)

   ic = 1/2 = dvcs/pi0


*/

struct particle_t
{
    TLorentzVector tlv;
    bool ic;
    int ind;
    double dt;
    double z;

    void Print()
    {
        cerr<<ic<<" "<<ind<<" "<<tlv.E()<<endl;
    };
};
struct dvcstree_t
{
    int hel;
    double wp, Q2, xb, y_q0k0, z_gq0, nu, t;
    double P0target; 
    double the_q0g, thegx,cop,phi_h,phi2,phig,therx;

    double m2_kpgX, p_kpgX, px_kpgX, py_kpgX, pz_kpgX, pt_kpgX, e_kpgX;
    double m2_kpX, m_kgX;
    int ic;

    double m_kpi0X, e_kppi0X, m2_kppi0X, the_kpi0X_p;
    double pt_kppi0X, cop_pi0, the_pi0_kpX, dphi_kpi0X_p;
    double m_pi0, e_pi0, px_pi0, py_pi0, pz_pi0, the_ggpi0;
    double phi_pi0, the_pi0;

    double v_x, v_y, v_z;
    double e_k0;
    double e_k1, px_k1, py_k1, pz_k1, phi_k1, the_k1;
    double p_p0;
    double e_p1, px_p1, py_p1, pz_p1, p_p1, phi_p1, the_p1, m_p1;
    double p_p1L, phi_p1L, the_p1L;
    double p1_delta_p, p1_delta_phi, p1_delta_theta;
    double e_q0, p_q0,  px_q0, py_q0, pz_q0;
    double e_g[2], px_g[2], py_g[2], pz_g[2], the_g[2], phi_g[2];


    vector <particle_t> gammas; 
    vector <particle_t> gammasEC; 
    vector <particle_t> gammasIC; 
};
#endif


