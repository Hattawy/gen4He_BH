#define analysis_dvcs_4He_t_cxx
#include "analysis_dvcs_4He_t.h"
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <time.h>
#include <TMultiGraph.h>
#include <TVector3.h>
#include <TMath.h>
#include <TProfile.h>
#include "TCutG.h"
#include <TGraphErrors.h>
#include <TLatex.h>
#include "TGraph.h"
#include<TLine.h>
#include <TLorentzVector.h>
#include "TSystem.h"
#include "TLegend.h"
#include<algorithm>
#include<vector>

 using namespace std;

void analysis_dvcs_4He_t::Loop()
{

   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetLabelSize(0.03,"xyz"); // size of axis value font
   gStyle->SetTitleSize(0.035,"xyz"); // size of axis title font
   gStyle->SetTitleFont(22,"xyz"); // font option
   gStyle->SetLabelFont(22,"xyz");
   //gStyle->SetTitleOffSet(1.2,"y");
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasBorderSize(0);
   gStyle->SetPadBottomMargin(0.16); //margins...
   gStyle->SetPadTopMargin(0.16);
   gStyle->SetPadLeftMargin(0.16);
   gStyle->SetPadRightMargin(0.16);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetPaperSize(20,24);
   gStyle->SetLabelSize(0.05,"xy");
   gStyle->SetTitleSize(0.06,"xy");




      // ---------------------------------------------------
     // bins for the projections
     // as a function of t in xB bins

      const int n_xB = 3;
      const int n_t = 8;
      const int n_phi = 13;
    
      double xB_lims[n_xB+1] = {0.095, 0.18, 0.24, 0.65};
      double t_lims[n_t+1] = {0.058, 0.07, 0.08, 0.09, 0.11,  0.14,  0.22,  0.36, 0.8};
      double phi_lims[n_phi+1]; 
      double PHI[n_phi];
      double PHI_err[n_phi]; 
 
      for(int jj=0; jj<n_phi; jj++){
            PHI[jj] = (360.0/n_phi)*jj+ (360.0/(2*n_phi));
            PHI_err[jj] = 0.0;
          std::cout<<PHI[jj]<<std::endl;
         }

      for(int jj=0; jj<n_phi+1; jj++){
         phi_lims[jj] = (360.0/n_phi)*jj;
         }

      double dvcs_N_p[n_xB][n_t][n_phi];
      double dvcs_N_m[n_xB][n_t][n_phi];
      double bsa_coh[n_xB][n_t][n_phi];
      double bsa_coh_err[n_xB][n_t][n_phi];
    
        for(int jj=0; jj<n_xB; jj++){
           for(int kk=0; kk<n_t; kk++){
              for(int nn=0; nn<n_phi; nn++){
              dvcs_N_p[jj][kk][nn]  = 0.0;
              dvcs_N_m[jj][kk][nn]  = 0.0;
              bsa_coh[jj][kk][nn] = 0.0;
              bsa_coh_err[jj][kk][nn] = 0.0;
           }}}

 // as a function of -t in xB bins
 TH1D *h_t_xB_Coh[n_xB][n_t];
 TH1D *h_t_Q2_Coh[n_xB][n_t];
 TH1D *h_t_t_Coh[n_xB][n_t];

 for( int i=0; i<n_xB; i++){
    for( int j=0; j<n_t; j++){
        h_t_Q2_Coh[i][j] = new TH1D(Form("h_t_Q2_Coh[%d][%d]",i,j),"Q^{2} of e^{4}He#gamma events",150, 0.5, 10);
        h_t_xB_Coh[i][j] = new TH1D(Form("h_t_xB_Coh[%d][%d]",i,j),"x_{B} of e^{4}He#gamma events",150, 0.05, 0.7);
        h_t_t_Coh[i][j] = new TH1D(Form("h_t_t_Coh_2[%d][%d]",i,j),"-t of e^{4}He#gamma events",150, 0.0, 0.7);
          }
       }


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

 Long64_t nbytes = 0, nb = 0;
for (Long64_t jentry=0; jentry<nentries;jentry++)
  {

    if (jentry% 100000 == 0) printf("still running %d \n",(int)jentry);
   //   if (jentry== 30000) break;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
                  int helicity;
                  if(hel == 1) helicity = 1;
                  else if (hel == 0) helicity = -1;

                      int which_Q2 = -1;
                      int which_xB = -1;
                      int which_t  = -1;
                      int which_phi = -1;
       
                      //for(int ii=0; ii<n_Q2; ii++){
                      //   if(Q2_lims[ii]<=Q2 && Q2 <Q2_lims[ii+1]) which_Q2 = ii ; }
       
                      for(int ii=0; ii<n_xB; ii++){
                         if(xB_lims[ii]<=xb && xb <xB_lims[ii+1]) which_xB = ii ; }
       
                      for(int ii=0; ii<n_t; ii++){
                         if(t_lims[ii]<=t && t <t_lims[ii+1]) which_t = ii ; }
       
                      for(int ii=0; ii<n_phi; ii++){
                        if(phi_lims[ii]<=phi_h && phi_h <phi_lims[ii+1]) which_phi = ii ; }
                        
                       // as a function of -t in xB bins
                       if (which_xB != -1 && which_t != -1  && which_phi != -1 )
                        {
                         if (helicity == 1 ) dvcs_N_p[which_xB][which_t][which_phi] += 1.0 ;
                         else if (helicity == -1) dvcs_N_m[which_xB][which_t][which_phi] += 1.0;
                         h_t_Q2_Coh[which_xB][which_t]->Fill(Q2);
                         h_t_xB_Coh[which_xB][which_t]->Fill(xb);
                         h_t_t_Coh[which_xB][which_t]->Fill(t);
                        }

   }


    double alu_t_x[n_xB][n_t];
    double alu_t_x_err[n_xB][n_t];
    double mean_t[n_xB][n_t];
    double mean_x[n_xB][n_t];
    double mean_Q2[n_xB][n_t];
    double mean_t_err[n_xB][n_t];
    double mean_x_err[n_xB][n_t];
    double mean_Q2_err[n_xB][n_t];

    double c0_BH[n_xB][n_t], c1_BH[n_xB][n_t], c2_BH[n_xB][n_t]; 
    double A0[n_xB][n_t], A1[n_xB][n_t], A2[n_xB][n_t], A3[n_xB][n_t];
    double Im[n_xB][n_t], Im_err[n_xB][n_t], Re[n_xB][n_t], Re_err[n_xB][n_t];

    TF1 *ffit[n_xB][n_t];

    for(int jj=0; jj<n_xB; jj++){
       for(int kk=0; kk<n_t; kk++){
           alu_t_x[jj][kk] = 0.1*jj;
           alu_t_x_err[jj][kk] = 0.0;
           mean_Q2[jj][kk] = h_t_Q2_Coh[jj][kk]->GetMean();
           mean_x[jj][kk] = h_t_xB_Coh[jj][kk]->GetMean();
           mean_t[jj][kk] = h_t_t_Coh[jj][kk]->GetMean();

           mean_Q2_err[jj][kk] = 0.0;
           mean_x_err[jj][kk] = 0.0;
           mean_t_err[jj][kk] = 0.0;

           Im[jj][kk] = 10*jj + 5.0;
           Im_err[jj][kk] = 0.0;
           Re[jj][kk] = -10*jj - 5.0;
           Re_err[jj][kk] = 0.0;

           calculate_CFF(mean_Q2[jj][kk], mean_x[jj][kk], -1*mean_t[jj][kk], 
                         A0[jj][kk], A1[jj][kk], A2[jj][kk], A3[jj][kk], 
                         c0_BH[jj][kk], c1_BH[jj][kk], c2_BH[jj][kk]);

           cout<<A0[jj][kk]<<"   "<< A1[jj][kk]<<"  "<<A2[jj][kk]<<"  "<< A3[jj][kk]<<"  "<<
                c0_BH[jj][kk]<<"  "<< c1_BH[jj][kk]<<"  "<< c2_BH[jj][kk]<<endl;

           ffit[jj][kk] = new TF1(Form("ffit[%d][%d]",jj,kk), Form("%0.7f*[0]*sin(x*3.14/180.0) / (%0.7f+ %0.7f*cos(x*3.14/180.0) + %0.7f*cos(2*x*3.14/180.0) + %0.7f*([0]*[0] + [1]*[1]) + %0.7f*[1] + %0.7f*[1]*cos(x*3.14/180.0))", A0[jj][kk], c0_BH[jj][kk], c1_BH[jj][kk],c2_BH[jj][kk], A1[jj][kk], A2[jj][kk], A3[jj][kk]),0.0,360.0);
           ffit[jj][kk]->SetLineColor(kRed);

           ffit[jj][kk]->SetParName(0,"Im(H_{A})");
           ffit[jj][kk]->SetParName(1,"Re(H_{A})");

          }
        }

   // as a function of -t in xB bins
    for(int jj=0; jj<n_xB; jj++){
       for(int kk=0; kk<n_t; kk++){
          for(int nn=0; nn<n_phi; nn++){
             bsa_coh[jj][kk][nn] = Find_BSA(dvcs_N_p[jj][kk][nn], dvcs_N_m[jj][kk][nn]);
             bsa_coh_err[jj][kk][nn] = Find_BSA_err(dvcs_N_p[jj][kk][nn], dvcs_N_m[jj][kk][nn]);;
            }
         TCanvas *c3 = new TCanvas("c3","",750,600 ); c3->cd(); c3->SetGrid();
         TF1 *myfit = new TF1("myfit","[0]*sin(x*3.1416/180.0)/(1 + [1]*cos(x*3.1416/180.0))",0.0,360.0);

         TGraphErrors *BSA_Coherent_Phi = new TGraphErrors(n_phi,PHI,bsa_coh[jj][kk], PHI_err, bsa_coh_err[jj][kk]);
                       BSA_Coherent_Phi->SetTitle("Coherent A_{LU} ");
                       BSA_Coherent_Phi->GetYaxis()->SetTitle("A_{LU}");
                       BSA_Coherent_Phi->GetXaxis()->SetTitle("#phi_{h} [deg.]");
                       BSA_Coherent_Phi->SetMarkerStyle(21);
                       BSA_Coherent_Phi->GetXaxis()->SetLimits(-5.0,365.0);
                       BSA_Coherent_Phi->GetYaxis()->SetRangeUser(-0.8,1.0);
                       BSA_Coherent_Phi->Draw("AP");
                    // BSA_Coherent_Phi->Fit(Form("ffit[%d][%d]",jj,kk),"E Q I"); 
                    // Im[jj][kk] =  ffit[jj][kk]->GetParameter(0);
                    // Re[jj][kk] =  ffit[jj][kk]->GetParameter(1);   
                    // Im_err[jj][kk] = ffit[jj][kk]->GetParError(0);
                    // Re_err[jj][kk] = ffit[jj][kk]->GetParError(1);
                    // BSA_Coherent_Phi->Fit("myfit");
                       //alu_t_x_err[jj][kk] =  myfit->GetParError(0);
                       // calculate_ALU(mean_Q2[jj][kk], mean_x[jj][kk], -1*mean_t[jj][kk], 
                       //              Im[jj][kk], Im_err[jj][kk], Re[jj][kk], Re_err[jj][kk],
                       //              alu_t_x[jj][kk], alu_t_x_err[jj][kk]);
                       //alu_t_x[jj][kk] = 0.1 *(jj+2);
 

                  TLatex *l2= new TLatex(10.0,-0.1,Form("%.2f< x_{B} <%.2f",xB_lims[jj], xB_lims[jj+1]));
                          l2->Draw("same");
                  TLatex *l3= new TLatex(10.0,-0.25,Form("%.2f< -t <%.2f",t_lims[kk], t_lims[kk+1]));
                          l3->Draw("same");
                  c3->Print(Form("fig/BSA_Coherent_Phi_%d_t%d.png",jj,kk));
                 // c3->Print(Form("fig/BSA_Coherent_Phi_%d_t%d.C",jj,kk));
       }}


       TGraphErrors *BSA_Coherent_Phi_t_x[n_xB];
               BSA_Coherent_Phi_t_x[0]  = new TGraphErrors(6, mean_t[0], alu_t_x[0], mean_t_err[0], alu_t_x_err[0]);
               BSA_Coherent_Phi_t_x[0]->SetMarkerStyle(20);
               BSA_Coherent_Phi_t_x[0]->SetMarkerColor(kBlack);

               BSA_Coherent_Phi_t_x[1]  = new TGraphErrors(6, mean_t[1], alu_t_x[1], mean_t_err[1], alu_t_x_err[1]);
               BSA_Coherent_Phi_t_x[1]->SetMarkerStyle(21);
               BSA_Coherent_Phi_t_x[1]->SetMarkerColor(kRed);

               BSA_Coherent_Phi_t_x[2]  = new TGraphErrors(6, mean_t[2], alu_t_x[2], mean_t_err[2], alu_t_x_err[2]);
               BSA_Coherent_Phi_t_x[2]->SetMarkerStyle(22);
               BSA_Coherent_Phi_t_x[2]->SetMarkerColor(kBlue);

          TCanvas *c4 = new TCanvas("c4","",750,600 ); c4->cd();  c4->SetGrid();
          TMultiGraph *mg = new TMultiGraph();
                       mg->Add(BSA_Coherent_Phi_t_x[0]);
                       mg->Add(BSA_Coherent_Phi_t_x[1]);
                       mg->Add(BSA_Coherent_Phi_t_x[2]);

                    mg->Draw("AP");
                    mg->SetTitle("Coherent A_{LU} projections");
                    mg->GetYaxis()->SetTitle("A_{LU}(90^{#circ})");
                    mg->GetXaxis()->SetTitle("-t [GeV^{2}/c^{2}]");
                    mg->GetXaxis()->SetLimits(0,0.4);
                    mg->GetYaxis()->SetRangeUser(0.15,0.5);
                    c4->SetGrid();
                  
           TLegend* leg_t_Coh = new TLegend(0.65,0.75,0.9,0.9);
                    leg_t_Coh-> SetNColumns(1);  
                    leg_t_Coh->AddEntry(BSA_Coherent_Phi_t_x[0],Form("%.2f <x_{B}< %.2f", xB_lims[0], xB_lims[1]),"P");
                    //leg_t_Coh->AddEntry(BSA_Coherent_Phi_t_x[1],Form("%.2f <x_{B}< %.2f", xB_lims[1], xB_lims[2]),"P");
                    //leg_t_Coh->AddEntry(BSA_Coherent_Phi_t_x[2],Form("%.2f <x_{B}< %.2f", xB_lims[2], xB_lims[3]),"P");
                    leg_t_Coh->Draw(); 
                   
             c4->Print("fig/BSA1_projections_Coherent_t.png");
             c4->Print("fig/BSA1_projections_Coherent_t.C");



       TGraphErrors *Im_CFF_Coherent_t[n_xB];
                     Im_CFF_Coherent_t[0]  = new TGraphErrors(6, mean_t[0], Im[0], mean_t_err[0], Im_err[0]);
                     Im_CFF_Coherent_t[0]->SetMarkerStyle(20);
                     Im_CFF_Coherent_t[0]->SetMarkerColor(kBlack);
                     Im_CFF_Coherent_t[1]  = new TGraphErrors(6, mean_t[1], Im[1], mean_t_err[1], Im_err[1]);
                     Im_CFF_Coherent_t[1]->SetMarkerStyle(21);
                     Im_CFF_Coherent_t[1]->SetMarkerColor(kRed);
                     Im_CFF_Coherent_t[2]  = new TGraphErrors(6, mean_t[2], Im[2], mean_t_err[2], Im_err[2]);
                     Im_CFF_Coherent_t[2]->SetMarkerStyle(22);
                     Im_CFF_Coherent_t[2]->SetMarkerColor(kBlue);

       TGraphErrors *Re_CFF_Coherent_t[n_xB];
                     Re_CFF_Coherent_t[0]  = new TGraphErrors(6, mean_t[0], Re[0], mean_t_err[0], Re_err[0]);
                     Re_CFF_Coherent_t[0]->SetMarkerStyle(20);
                     Re_CFF_Coherent_t[0]->SetMarkerColor(kBlack);
                     Re_CFF_Coherent_t[1]  = new TGraphErrors(6, mean_t[1], Re[1], mean_t_err[1], Re_err[1]);
                     Re_CFF_Coherent_t[1]->SetMarkerStyle(21);
                     Re_CFF_Coherent_t[1]->SetMarkerColor(kRed);
                     Re_CFF_Coherent_t[2]  = new TGraphErrors(6, mean_t[2], Re[2], mean_t_err[2], Re_err[2]);
                     Re_CFF_Coherent_t[2]->SetMarkerStyle(22);
                     Re_CFF_Coherent_t[2]->SetMarkerColor(kBlue);


          TCanvas *c5 = new TCanvas("c5","",750,600 ); c5->cd();  c5->SetGrid();
          TMultiGraph *mmg = new TMultiGraph();
                       mmg->Add(Im_CFF_Coherent_t[0]);
                       mmg->Add(Im_CFF_Coherent_t[1]);
                       mmg->Add(Im_CFF_Coherent_t[2]);
                      // mmg->Add(Re_CFF_Coherent_t[0]);
                      // mmg->Add(Re_CFF_Coherent_t[1]);
                      // mmg->Add(Re_CFF_Coherent_t[2]);

                    mmg->Draw("AP");
                    mmg->SetTitle("CFF H_{A} projections");
                    //mmg->GetYaxis()->SetTitle("H_{A}(Re, Im))");
                    mmg->GetYaxis()->SetTitle("Im(H_{A})");
                    mmg->GetXaxis()->SetTitle("-t [GeV^{2}/c^{2}]");
                    mmg->GetXaxis()->SetLimits(0,0.4);
                    mmg->GetYaxis()->SetRangeUser(-60,60);
                    c5->SetGrid();
                  
           TLegend* leg_t_cff = new TLegend(0.65,0.75,0.9,0.9);
                    leg_t_cff-> SetNColumns(1);  
                    leg_t_cff->AddEntry(Im_CFF_Coherent_t[0],Form("%.2f <x_{B}< %.2f", xB_lims[0], xB_lims[1]),"P");
                    leg_t_cff->AddEntry(Im_CFF_Coherent_t[1],Form("%.2f <x_{B}< %.2f", xB_lims[1], xB_lims[2]),"P");
                    leg_t_cff->AddEntry(Im_CFF_Coherent_t[2],Form("%.2f <x_{B}< %.2f", xB_lims[2], xB_lims[3]),"P");
                    leg_t_cff->Draw(); 
                   
             c5->Print("fig/CFF_projections_Coherent_t.png");
             c5->Print("fig/CFF_projections_Coherent_t.C");


  }


