void BSA1_projections_Coherent_t()
{
//=========Macro generated from canvas: c4/
//=========  (Mon May  8 11:10:20 2017) by ROOT version6.09/01
   TCanvas *c4 = new TCanvas("c4", "",0,0,750,600);
   gStyle->SetOptStat(0);
   c4->Range(-0.09411764,0.06764707,0.4941176,0.5823529);
   c4->SetFillColor(0);
   c4->SetBorderMode(0);
   c4->SetBorderSize(0);
   c4->SetGridx();
   c4->SetGridy();
   c4->SetLeftMargin(0.16);
   c4->SetRightMargin(0.16);
   c4->SetTopMargin(0.16);
   c4->SetBottomMargin(0.16);
   c4->SetFrameBorderMode(0);
   c4->SetFrameBorderMode(0);
   
   TMultiGraph *multigraph = new TMultiGraph();
   multigraph->SetName("");
   multigraph->SetTitle("Coherent A_{LU} projections");
   
   Double_t Graph_fx1001[6] = {
   0.06660288,
   0.07487206,
   0.08490888,
   0.09888697,
   0.1238392,
   0.181615};
   Double_t Graph_fy1001[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph_fex1001[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph_fey1001[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   TGraphErrors *gre = new TGraphErrors(6,Graph_fx1001,Graph_fy1001,Graph_fex1001,Graph_fey1001);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_Graph1001 = new TH1F("Graph_Graph1001","Graph",100,0.05510166,0.1931162);
   Graph_Graph1001->SetMinimum(0);
   Graph_Graph1001->SetMaximum(1.1);
   Graph_Graph1001->SetDirectory(0);
   Graph_Graph1001->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph1001->SetLineColor(ci);
   Graph_Graph1001->GetXaxis()->SetLabelFont(22);
   Graph_Graph1001->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1001->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph1001->GetXaxis()->SetTitleFont(22);
   Graph_Graph1001->GetYaxis()->SetLabelFont(22);
   Graph_Graph1001->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1001->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph1001->GetYaxis()->SetTitleFont(22);
   Graph_Graph1001->GetZaxis()->SetLabelFont(22);
   Graph_Graph1001->GetZaxis()->SetLabelSize(0.03);
   Graph_Graph1001->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1001->GetZaxis()->SetTitleFont(22);
   gre->SetHistogram(Graph_Graph1001);
   
   multigraph->Add(gre,"");
   
   Double_t Graph_fx1002[6] = {
   0.06633049,
   0.07455528,
   0.08434492,
   0.1002237,
   0.1241009,
   0.1785441};
   Double_t Graph_fy1002[6] = {
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1};
   Double_t Graph_fex1002[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph_fey1002[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   gre = new TGraphErrors(6,Graph_fx1002,Graph_fy1002,Graph_fex1002,Graph_fey1002);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(21);
   
   TH1F *Graph_Graph1002 = new TH1F("Graph_Graph1002","Graph",100,0.05510913,0.1897654);
   Graph_Graph1002->SetMinimum(0);
   Graph_Graph1002->SetMaximum(1.2);
   Graph_Graph1002->SetDirectory(0);
   Graph_Graph1002->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1002->SetLineColor(ci);
   Graph_Graph1002->GetXaxis()->SetLabelFont(22);
   Graph_Graph1002->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1002->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph1002->GetXaxis()->SetTitleFont(22);
   Graph_Graph1002->GetYaxis()->SetLabelFont(22);
   Graph_Graph1002->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1002->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph1002->GetYaxis()->SetTitleFont(22);
   Graph_Graph1002->GetZaxis()->SetLabelFont(22);
   Graph_Graph1002->GetZaxis()->SetLabelSize(0.03);
   Graph_Graph1002->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1002->GetZaxis()->SetTitleFont(22);
   gre->SetHistogram(Graph_Graph1002);
   
   multigraph->Add(gre,"");
   
   Double_t Graph_fx1003[6] = {
   0.06738587,
   0.07566904,
   0.08497524,
   0.1007518,
   0.1266077,
   0.1851955};
   Double_t Graph_fy1003[6] = {
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2};
   Double_t Graph_fex1003[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph_fey1003[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   gre = new TGraphErrors(6,Graph_fx1003,Graph_fy1003,Graph_fex1003,Graph_fey1003);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(22);
   
   TH1F *Graph_Graph1003 = new TH1F("Graph_Graph1003","Graph",100,0.05560491,0.1969765);
   Graph_Graph1003->SetMinimum(0.1);
   Graph_Graph1003->SetMaximum(1.3);
   Graph_Graph1003->SetDirectory(0);
   Graph_Graph1003->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1003->SetLineColor(ci);
   Graph_Graph1003->GetXaxis()->SetLabelFont(22);
   Graph_Graph1003->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1003->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph1003->GetXaxis()->SetTitleFont(22);
   Graph_Graph1003->GetYaxis()->SetLabelFont(22);
   Graph_Graph1003->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1003->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph1003->GetYaxis()->SetTitleFont(22);
   Graph_Graph1003->GetZaxis()->SetLabelFont(22);
   Graph_Graph1003->GetZaxis()->SetLabelSize(0.03);
   Graph_Graph1003->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1003->GetZaxis()->SetTitleFont(22);
   gre->SetHistogram(Graph_Graph1003);
   
   multigraph->Add(gre,"");
   multigraph->Draw("AP");
   multigraph->GetXaxis()->SetTitle("-t [GeV^{2}/c^{2}]");
   multigraph->GetXaxis()->SetLabelFont(22);
   multigraph->GetXaxis()->SetLabelSize(0.05);
   multigraph->GetXaxis()->SetTitleSize(0.06);
   multigraph->GetXaxis()->SetTitleFont(22);
   multigraph->GetYaxis()->SetTitle("A_{LU}(90^{#circ})");
   multigraph->GetYaxis()->SetLabelFont(22);
   multigraph->GetYaxis()->SetLabelSize(0.05);
   multigraph->GetYaxis()->SetTitleSize(0.06);
   multigraph->GetYaxis()->SetTitleFont(22);
   
   TPaveText *pt = new TPaveText(0.2580965,0.9257692,0.7419035,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *AText = pt->AddText("Coherent A_{LU} projections");
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
   leg->Draw();
   c4->Modified();
   c4->cd();
   c4->SetSelected(c4);
}
