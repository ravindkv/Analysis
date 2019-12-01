#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

using namespace std;
//INPUT FILES
TFile* fitDiagOut = TFile::Open("fitDiagnostics_data_cminDS0.root");

//SAVE HISTOS ON DISK
bool isSaveHisto = true;

//--------------------------------------------//
//various functions
//--------------------------------------------//
TPaveText *paveText(double minX, double minY, double maxX, double maxY, int lineColor, int fillColor, int size, int style, int font );
//--------------------------------------------//
//stack histos
//--------------------------------------------//
void example_stack(TString chName, TString chLable, TString charmCat){
  //Pad
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(2);
  TCanvas *c1 = new TCanvas();
  TGaxis::SetMaxDigits(3);
  //TCanvas *c1 = new TCanvas("c1", "Data_MC", 400, 600);
  const float xpad[2] = {0.,1};
  const float ypad[4] = {0.,0.2351916,0.2351916,0.98};
  c1->Divide(1,2); c1->cd(1);
  gPad->SetBottomMargin(0.0);
  gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]);

  //Data
  TH1F* hTTPreFit = (TH1F*)(fitDiagOut->Get("shapes_prefit/"+chName+"/ttbar"))->Clone("ttbar_pre");
  double scale_prefit = (1-0.1)*(1-0.1);
  hTTPreFit->Scale(1/scale_prefit);
  hTTPreFit->SetMarkerStyle(20); 
  hTTPreFit->SetMarkerSize(0.7);
  hTTPreFit->SetLineColor(kGreen);
  hTTPreFit->SetLineWidth(2);
  hTTPreFit->SetTitle("");
  hTTPreFit->GetYaxis()->SetTitle("Events");
  hTTPreFit->GetYaxis()->SetRangeUser(1, 1.2* hTTPreFit->GetMaximum());
  //hTTPreFit->GetXaxis()->SetTitle("m_{jj} [(x5 + 20)GeV]");
  hTTPreFit->GetYaxis()->SetTitleOffset(0.70);
  hTTPreFit->GetXaxis()->SetTitleOffset(1.00);
  hTTPreFit->GetYaxis()->SetTitleSize(0.07);   
  hTTPreFit->GetXaxis()->SetTitleSize(0.07);
  hTTPreFit->GetXaxis()->SetLabelSize(0.07);   
  hTTPreFit->GetYaxis()->SetLabelSize(0.07);   
  hTTPreFit->GetXaxis()->SetTickLength(0.05); 
  hTTPreFit->GetYaxis()->SetTickLength(0.04); 
  hTTPreFit->Draw("E"); 
  
  TH1F* hTTPostFit = (TH1F*)(fitDiagOut->Get("shapes_fit_b/"+chName+"/ttbar"))->Clone("ttbar_post");
  hTTPostFit->SetMarkerStyle(22); 
  hTTPostFit->SetMarkerSize(0.7);
  hTTPostFit->SetLineColor(kBlue);
  hTTPostFit->SetLineWidth(2);
  hTTPostFit->Draw("Esame"); 
  TLegend* leg = new TLegend(0.5618792,0.7061504,0.6712081,0.8798861,NULL,"brNDC");
  leg->SetFillStyle(0); leg->SetBorderSize(0);
  leg->SetFillColor(10); leg->SetTextSize(0.06);
  leg->AddEntry(hTTPreFit,  "SM t#bar{t} + jets, pre-fit","PE"); 
  leg->AddEntry(hTTPostFit, "SM t#bar{t} + jets, post-fit","PE"); 
  leg->Draw();
  
  //-------------------------------------///
  //  Draw Pave Text 
  //-------------------------------------///
  
  //channel
  TPaveText *ch = paveText(0.203,0.7061504,0.3010067,0.8798861, 0, 19, 1, 0, 132);
  ch->SetTextSize(0.07);
  ch->AddText(chLable +", "+charmCat);
  //CMS prili
  TPaveText *pt = paveText(0.30,0.9354,0.82,0.9362, 0, 19, 1, 0, 132);
  pt->SetTextSize(0.07);
  TText *text = pt->AddText("CMS Preliminary, 13 TeV (35.9 fb^{-1})");
  text->SetTextAlign(11);
  pt->Draw();
  ch->Draw();
  gPad->RedrawAxis();
  c1->Update();
  
  //-------------------------------------///
  // Ratio = DATA/Bkg
  //-------------------------------------///
  c1->cd(2);
  gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[1]);
  gPad->SetTopMargin(0); 
  gPad->SetBottomMargin(0.5); //gPad->SetGridy();
  //gPad->SetLeftMargin(0.10);
  //gPad->SetRightMargin(0.05);
  TH1F *hRatio = (TH1F*)hTTPreFit->Clone("hRatio");
  hRatio->Reset();
  hRatio->Add(hTTPreFit);
  hRatio->GetYaxis()->SetNdivisions(5);
  hRatio->Divide(hTTPostFit); 
  hRatio->SetMarkerStyle(20); 
  hRatio->SetMarkerSize(0.7);
  hRatio->SetMarkerColor(kBlack); 
  hRatio->SetLineColor(kBlack); 
  hRatio->GetYaxis()->SetRangeUser(0.7, 1.3);
  hRatio->GetXaxis()->SetTickLength(0.13); 
  hRatio->GetYaxis()->SetTickLength(0.04); 
  hRatio->GetXaxis()->SetTitle("m_{jj} [(bin x 5 + 20) GeV]"); 
  hRatio->GetXaxis()->SetLabelSize(0.20); 
  hRatio->GetYaxis()->SetLabelSize(0.10); 
  hRatio->GetXaxis()->SetTitleSize(0.20); 
  hRatio->GetYaxis()->SetTitleSize(0.15); 
  hRatio->GetYaxis()->SetTitleOffset(0.27);
  hRatio->GetXaxis()->SetTitleOffset(1.0);
  hRatio->GetYaxis()->SetTitle("#frac{Pre-fit}{Post-fit}"); 
  hRatio->GetYaxis()->CenterTitle();
  hRatio->Draw("E"); // use "P" or "AP"

  //base line at 1
  TF1 *baseLine = new TF1("baseLine","1", -100, 2000); 
  baseLine->SetLineColor(kCyan);
  baseLine->Draw("SAME");
  hRatio->Draw("Esame"); // use "P" or "AP"
  c1->Update();
  if(isSaveHisto){
    TString outFile("$PWD/output/ttbarPostFit_"+chName+".pdf");
    c1->SaveAs(outFile);
    c1->Close();
  }
}

TPaveText *paveText(double minX, double minY, double maxX, double maxY, int lineColor, int fillColor, int size, int style, int font ){
  TPaveText *pt = new TPaveText(minX, minY, maxX, maxY, "brNDC"); // good_v1
  pt->SetBorderSize(size);
  pt->SetFillColor(fillColor);
  pt->SetFillStyle(style);
  pt->SetLineColor(lineColor);
  pt->SetTextFont(font);
  return pt;
}

void MyPostFitTTbar(){
  example_stack("ch1", "#mu + jets", "loose");
  example_stack("ch2", "#mu + jets", "medium");
  example_stack("ch3", "#mu + jets", "tight");
  example_stack("ch4", "e + jets", "loose");
  example_stack("ch5", "e + jets", "medium");
  example_stack("ch6", "e + jets", "tight");
} 

