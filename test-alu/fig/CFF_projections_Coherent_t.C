void CFF_projections_Coherent_t()
{
//=========Macro generated from canvas: c5/
//=========  (Mon May  8 11:10:20 2017) by ROOT version6.09/01
   TCanvas *c5 = new TCanvas("c5", "",0,0,750,600);
   gStyle->SetOptStat(0);
   c5->Range(-0.09411764,-88.23529,0.4941176,88.23529);
   c5->SetFillColor(0);
   c5->SetBorderMode(0);
   c5->SetBorderSize(0);
   c5->SetGridx();
   c5->SetGridy();
   c5->SetLeftMargin(0.16);
   c5->SetRightMargin(0.16);
   c5->SetTopMargin(0.16);
   c5->SetBottomMargin(0.16);
   c5->SetFrameBorderMode(0);
   c5->SetFrameBorderMode(0);
   
   TMultiGraph *multigraph = new TMultiGraph();
   multigraph->SetName("");
   multigraph->SetTitle("CFF H_{A} projections");
   
   Double_t Graph_fx1004[6] = {
   0.06660288,
   0.07487206,
   0.08490888,
   0.09888697,
   0.1238392,
   0.181615};
   Double_t Graph_fy1004[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t Graph_fex1004[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph_fey1004[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   TGraphErrors *gre = new TGraphErrors(6,Graph_fx1004,Graph_fy1004,Graph_fex1004,Graph_fey1004);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_Graph1004 = new TH1F("Graph_Graph1004","Graph",100,0.05510166,0.1931162);
   Graph_Graph1004->SetMinimum(4.9);
   Graph_Graph1004->SetMaximum(6.1);
   Graph_Graph1004->SetDirectory(0);
   Graph_Graph1004->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph1004->SetLineColor(ci);
   Graph_Graph1004->GetXaxis()->SetLabelFont(22);
   Graph_Graph1004->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1004->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph1004->GetXaxis()->SetTitleFont(22);
   Graph_Graph1004->GetYaxis()->SetLabelFont(22);
   Graph_Graph1004->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1004->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph1004->GetYaxis()->SetTitleFont(22);
   Graph_Graph1004->GetZaxis()->SetLabelFont(22);
   Graph_Graph1004->GetZaxis()->SetLabelSize(0.03);
   Graph_Graph1004->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1004->GetZaxis()->SetTitleFont(22);
   gre->SetHistogram(Graph_Graph1004);
   
   multigraph->Add(gre,"");
   
   Double_t Graph_fx1005[6] = {
   0.06633049,
   0.07455528,
   0.08434492,
   0.1002237,
   0.1241009,
   0.1785441};
   Double_t Graph_fy1005[6] = {
   15,
   15,
   15,
   15,
   15,
   15};
   Double_t Graph_fex1005[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph_fey1005[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   gre = new TGraphErrors(6,Graph_fx1005,Graph_fy1005,Graph_fex1005,Graph_fey1005);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(21);
   
   TH1F *Graph_Graph1005 = new TH1F("Graph_Graph1005","Graph",100,0.05510913,0.1897654);
   Graph_Graph1005->SetMinimum(14.9);
   Graph_Graph1005->SetMaximum(16.1);
   Graph_Graph1005->SetDirectory(0);
   Graph_Graph1005->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1005->SetLineColor(ci);
   Graph_Graph1005->GetXaxis()->SetLabelFont(22);
   Graph_Graph1005->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1005->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph1005->GetXaxis()->SetTitleFont(22);
   Graph_Graph1005->GetYaxis()->SetLabelFont(22);
   Graph_Graph1005->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1005->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph1005->GetYaxis()->SetTitleFont(22);
   Graph_Graph1005->GetZaxis()->SetLabelFont(22);
   Graph_Graph1005->GetZaxis()->SetLabelSize(0.03);
   Graph_Graph1005->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1005->GetZaxis()->SetTitleFont(22);
   gre->SetHistogram(Graph_Graph1005);
   
   multigraph->Add(gre,"");
   
   Double_t Graph_fx1006[6] = {
   0.06738587,
   0.07566904,
   0.08497524,
   0.1007518,
   0.1266077,
   0.1851955};
   Double_t Graph_fy1006[6] = {
   25,
   25,
   25,
   25,
   25,
   25};
   Double_t Graph_fex1006[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph_fey1006[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   gre = new TGraphErrors(6,Graph_fx1006,Graph_fy1006,Graph_fex1006,Graph_fey1006);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(22);
   
   TH1F *Graph_Graph1006 = new TH1F("Graph_Graph1006","Graph",100,0.05560491,0.1969765);
   Graph_Graph1006->SetMinimum(24.9);
   Graph_Graph1006->SetMaximum(26.1);
   Graph_Graph1006->SetDirectory(0);
   Graph_Graph1006->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1006->SetLineColor(ci);
   Graph_Graph1006->GetXaxis()->SetLabelFont(22);
   Graph_Graph1006->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1006->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph1006->GetXaxis()->SetTitleFont(22);
   Graph_Graph1006->GetYaxis()->SetLabelFont(22);
   Graph_Graph1006->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1006->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph1006->GetYaxis()->SetTitleFont(22);
   Graph_Graph1006->GetZaxis()->SetLabelFont(22);
   Graph_Graph1006->GetZaxis()->SetLabelSize(0.03);
   Graph_Graph1006->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1006->GetZaxis()->SetTitleFont(22);
   gre->SetHistogram(Graph_Graph1006);
   
   multigraph->Add(gre,"");
   multigraph->Draw("AP");
   multigraph->GetXaxis()->SetTitle("-t [GeV^{2}/c^{2}]");
   multigraph->GetXaxis()->SetLabelFont(22);
   multigraph->GetXaxis()->SetLabelSize(0.05);
   multigraph->GetXaxis()->SetTitleSize(0.06);
   multigraph->GetXaxis()->SetTitleFont(22);
   multigraph->GetYaxis()->SetTitle("Im(H_{A})");
   multigraph->GetYaxis()->SetLabelFont(22);
   multigraph->GetYaxis()->SetLabelSize(0.05);
   multigraph->GetYaxis()->SetTitleSize(0.06);
   multigraph->GetYaxis()->SetTitleFont(22);
   
   TPaveText *pt = new TPaveText(0.3143968,0.9278671,0.6856032,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *AText = pt->AddText("CFF H_{A} projections");
   pt->Draw();
   
   TLegend *leg = new TLegend(0.65,0.75,0.9,0.9,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph","0.10 <x_{B}< 0.18","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","0.18 <x_{B}< 0.24","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","0.24 <x_{B}< 0.65","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   c5->Modified();
   c5->cd();
   c5->SetSelected(c5);
}
