#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>
#include <cstring>

#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TGraphPainter.h"
#include "TMultiGraph.h"
#include "TTree.h"

using namespace std;


//-------------------------------------------//
///INPUT FILES
//-------------------------------------------//
//data
TFile* fData    = TFile::Open("all_muData.root");
//bkg
TFile* fVV      = TFile::Open("all_VV.root");
TFile* fDY      = TFile::Open("all_DY.root");
TFile* fWJ      = TFile::Open("all_WJets.root");
TFile* fQCD     = TFile::Open("all_QCD.root");
TFile* fST      = TFile::Open("all_ST.root");
TFile* fTT      = TFile::Open("all_TTJetsP.root");
//signal
TFile *fWH80      = TFile::Open("all_Hplus80.root");
TFile *fWH90      = TFile::Open("all_Hplus90.root");
TFile *fWH100     = TFile::Open("all_Hplus100.root");
TFile *fWH120     = TFile::Open("all_Hplus120.root");
TFile *fWH140     = TFile::Open("all_Hplus140.root");
TFile *fWH150     = TFile::Open("all_Hplus150.root");
TFile *fWH155     = TFile::Open("all_Hplus155.root");
TFile *fWH160     = TFile::Open("all_Hplus160.root");
//data driven qcd
///TFile* fQCD_dd = TFile::Open("all_QCD_dd.root"); 

//histogram name and range
TString histName = "pt_bjet" ;
Float_t xMin_ = 0.0 ;
Float_t xMax_ = 300.0 ;
bool isSaveHisto = true;

//-------------------------------------------//
//function to get histogram from root file
//-------------------------------------------//
TH1F* getHisto(TFile *file, TString histName, int histColor, Float_t xMin_ = 0.0, Float_t xMax_ = 500){

  if(file->IsZombie() || file->IsZombie()){
    cout << "Cannot open file "<< endl;
  }
  TH1F* h;
  h = (TH1F*)file->Get("base/Iso/KinFit/"+histName);
  h->SetMarkerColor(kRed);
  h->SetMarkerStyle(1);
  h->SetMarkerSize(1.5);
  h->SetLineColor(histColor);
  h->SetLineWidth(3);

  h->GetXaxis()->SetRangeUser(xMin_, xMax_);
  h->GetYaxis()->SetRangeUser(0, 40000);
  h->GetYaxis()->SetTitleOffset(1.33);
  //h->SetMinimum(yMin_);
  //h->SetMaximum(yMax_);
  h->GetXaxis()->SetTitle(histName);
  h->GetYaxis()->SetTitle("Events");
  h->SetTitle("");
  return h;
}

TH1F* getRatio(TH1F* h1, TH1F* h2, TString histName){
  TH1F *hRatio = (TH1F*)h2->Clone("hRatio");
  hRatio->Divide(h1); 
  hRatio->SetMarkerStyle(20); 
  hRatio->SetMarkerSize(0.8);
  hRatio->SetMarkerColor(kBlack); 
  hRatio->SetLineColor(kBlack); 
  hRatio->GetYaxis()->SetRangeUser(0, 2);
  hRatio->GetXaxis()->SetTitle(histName); 
  hRatio->GetYaxis()->SetTitleOffset(1.33);
  hRatio->GetYaxis()->SetTitle("Ratio"); 
  //hRatio->GetYaxis()->CenterTitle();
  //hRatio->GetYaxis()->SetTitleSize(0.02); 
  //hRatio->GetXaxis()->SetTitleSize(0.02);
  //hRatio->GetXaxis()->SetLabelSize(0.12); 
  //hRatio->GetXaxis()->LabelsOption("u"); // extra
  //hRatio->GetYaxis()->SetLabelSize(0.04);
  return hRatio;
}

//-------------------------------------------//
//overlap two histograms
//-------------------------------------------//
void overLap2Histos(TFile *f1, TFile* f2, TString histName1, TString histName2, TString histLeg1, TString histLeg2){
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas();
  //TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->Divide(1,2);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  TLegend* leg = new TLegend(0.59,0.60,0.46,0.85,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  
  //overlay
  c1->cd(1);
  TH1F* h1 = getHisto(f1, histName1, kBlack, xMin_, xMax_);
  TH1F* h2 = getHisto(f2, histName2, kBlue, xMin_, xMax_);
  h1->Draw("HIST");
  h2->Draw("HISTSAME");
  leg->AddEntry(h1, histLeg1,"PL");
  leg->AddEntry(h2, histLeg2,"PL");
  leg->Draw();
  //ratio
  c1->cd(2);
  //leg->SetHeader(Form("#splitline{CMS Preliminary #sqrt{s}=13 TeV}{35.5 fb^{-1}}"));
  TH1F* hRatio = getRatio(h1, h2, histName);
  hRatio->Draw("E"); // use "P" or "AP"
  //base line at 1
  TF1 *baseLine = new TF1("baseLine","1", -100, 2000); 
  baseLine->SetLineColor(kRed);
  baseLine->Draw("SAME");
  
  if(isSaveHisto){
    TString outFile("$PWD/");
    outFile += histName;
    outFile += "_mu"+histLeg1+".pdf";
    c1->SaveAs(outFile);
    c1->Close();
  }
}

void overLapRatios(){
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas();
  TLegend* leg = new TLegend(0.59,0.60,0.46,0.85,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  
  c1->cd();
  TH1F* hRL_TTWH80 = getRatio(getHisto(fTT, "pt_jet", kBlack, xMin_, xMax_), getHisto(fWH80, "pt_jet", kBlue, xMin_, xMax_), "pt_jet");
  hRL_TTWH80->SetTitle("Pt of Had bjet : electron + jets channel");
  hRL_TTWH80->SetLineColor(1);
  hRL_TTWH80->Draw();
  
  TH1F* hRL_TTWH90 = getRatio(getHisto(fTT, "pt_jet", kBlack, xMin_, xMax_), getHisto(fWH90, "pt_jet", kBlue, xMin_, xMax_), "pt_jet");
  hRL_TTWH90->SetLineColor(2);
  hRL_TTWH90->Draw("SAME");
  
  TH1F* hRL_TTWH100 = getRatio(getHisto(fTT, "pt_jet", kBlack, xMin_, xMax_), getHisto(fWH100, "pt_jet", kBlue, xMin_, xMax_), "pt_jet");
  hRL_TTWH100->SetLineColor(3);
  hRL_TTWH100->Draw("SAME");
 
  TH1F* hRL_TTWH120 = getRatio(getHisto(fTT, "pt_jet", kBlack, xMin_, xMax_), getHisto(fWH120, "pt_jet", kBlue, xMin_, xMax_), "pt_jet");
  hRL_TTWH120->SetLineColor(4);
  hRL_TTWH120->Draw("SAME");

  TH1F* hRL_TTWH140 = getRatio(getHisto(fTT, "pt_jet", kBlack, xMin_, xMax_), getHisto(fWH140, "pt_jet", kBlue, xMin_, xMax_), "pt_jet");
  hRL_TTWH140->SetLineColor(5);
  hRL_TTWH140->Draw("SAME");
  
  TH1F* hRL_TTWH150 = getRatio(getHisto(fTT, "pt_jet", kBlack, xMin_, xMax_), getHisto(fWH150, "pt_jet", kBlue, xMin_, xMax_), "pt_jet");
  hRL_TTWH150->SetLineColor(6);
  hRL_TTWH150->Draw("SAME");

  TH1F* hRL_TTWH155 = getRatio(getHisto(fTT, "pt_jet", kBlack, xMin_, xMax_), getHisto(fWH155, "pt_jet", kBlue, xMin_, xMax_), "pt_jet");
  hRL_TTWH155->SetLineColor(7);
  hRL_TTWH155->Draw("SAME");

  TH1F* hRL_TTWH160 = getRatio(getHisto(fTT, "pt_jet", kBlack, xMin_, xMax_), getHisto(fWH160, "pt_jet", kBlue, xMin_, xMax_), "pt_jet");
  hRL_TTWH160->SetLineColor(8);
  hRL_TTWH160->Draw("SAME");

  leg->AddEntry(hRL_TTWH80,  "WH80/TT","PL");
  leg->AddEntry(hRL_TTWH90,  "WH90/TT","PL");
  leg->AddEntry(hRL_TTWH100, "WH100/TT","PL");
  leg->AddEntry(hRL_TTWH120, "WH120/TT","PL");
  leg->AddEntry(hRL_TTWH140, "WH140/TT","PL");
  leg->AddEntry(hRL_TTWH150, "WH150/TT","PL");
  leg->AddEntry(hRL_TTWH155, "WH155/TT","PL");
  leg->AddEntry(hRL_TTWH160, "WH160/TT","PL");
  leg->Draw();
  //base line at 1
  TF1 *baseLine = new TF1("baseLine","1", -100, 2000); 
  baseLine->SetLineColor(kRed);
  baseLine->Draw("SAME");
  
  if(isSaveHisto){
    TString outFile("$PWD/");
    outFile += histName;
    outFile += "_ele_Ratio.png";
    c1->SaveAs(outFile);
    c1->Close();
  }
}

void overLapHistos(){
  overLap2Histos(fTT,  fTT,  "pt_bjet", "pt_bjetL", "TT:Had bjet",  "TT:Lep bjet");
  overLap2Histos(fWH80,  fWH80,  "pt_bjet", "pt_bjetL", "WH80:Had bjet",  "WH80:Lep bjet");
  overLap2Histos(fWH90,  fWH90,  "pt_bjet", "pt_bjetL", "WH90:Had bjet",  "WH90:Lep bjet");
  overLap2Histos(fWH100, fWH100, "pt_bjet", "pt_bjetL", "WH100:Had bjet", "WH100:Lep bjet");
  overLap2Histos(fWH120, fWH120, "pt_bjet", "pt_bjetL", "WH120:Had bjet", "WH120:Lep bjet");
  overLap2Histos(fWH140, fWH140, "pt_bjet", "pt_bjetL", "WH140:Had bjet", "WH140:Lep bjet");
  overLap2Histos(fWH150, fWH150, "pt_bjet", "pt_bjetL", "WH150:Had bjet", "WH150:Lep bjet");
  overLap2Histos(fWH155, fWH155, "pt_bjet", "pt_bjetL", "WH155:Had bjet", "WH155:Lep bjet");
  overLap2Histos(fWH160, fWH160, "pt_bjet", "pt_bjetL", "WH160:Had bjet", "WH160:Lep bjet");
}

void overLapAllHistos(){
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas();
  TLegend* leg = new TLegend(0.59,0.60,0.46,0.85,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.02);
  leg->SetFillColor(0);
  
  //overlay
  c1->cd();
  TH1F* hTT = getHisto(fTT,       "pt_bjet", 1, xMin_, xMax_);
  TH1F* hWH80 = getHisto(fWH80,  "pt_bjet", 2, xMin_, xMax_);
  TH1F* hWH90 = getHisto(fWH90,  "pt_bjet", 3, xMin_, xMax_);
  TH1F* hWH100 = getHisto(fWH100, "pt_bjet", 4, xMin_, xMax_);
  TH1F* hWH120 = getHisto(fWH120, "pt_bjet", 5, xMin_, xMax_);
  TH1F* hWH140 = getHisto(fWH140, "pt_bjet", 6, xMin_, xMax_);
  TH1F* hWH150 = getHisto(fWH150, "pt_bjet", 7, xMin_, xMax_);
  TH1F* hWH155 = getHisto(fWH155, "pt_bjet", 8, xMin_, xMax_);
  TH1F* hWH160 = getHisto(fWH160, "pt_bjet", 9, xMin_, xMax_);
  
  hTT->Draw("HIST");
  hWH80 ->Draw("HISTSAME");
  hWH90 ->Draw("HISTSAME");
  hWH100->Draw("HISTSAME");
  hWH120->Draw("HISTSAME");
  hWH140->Draw("HISTSAME");
  hWH150->Draw("HISTSAME");
  hWH155->Draw("HISTSAME");
  hWH160->Draw("HISTSAME");
  leg->AddEntry(hTT, "ttbar","PL");
  leg->AddEntry(hWH80 , "WH80 ","PL");
  leg->AddEntry(hWH90 , "WH90 ","PL");
  leg->AddEntry(hWH100, "WH100","PL");
  leg->AddEntry(hWH120, "WH120","PL");
  leg->AddEntry(hWH140, "WH140","PL");
  leg->AddEntry(hWH150, "WH150","PL");
  leg->AddEntry(hWH155, "WH155","PL");
  leg->AddEntry(hWH160, "WH160","PL");
  leg->Draw();
  if(isSaveHisto){
    TString outFile("$PWD/");
    outFile += "_mu.pdf";
    c1->SaveAs(outFile);
    c1->Close();
  }
}
