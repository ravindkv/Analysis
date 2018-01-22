// Gouranga Kole (UC San Diego)
// Usage
// .L compareLimit.C
// compareLimit("limit_csbar_.root", "Test1/limit_csbar_.root", "Test1", "Test2", true, false, 0.00, 0.06);

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


TGraphAsymmErrors* getLimit( TString dir, TString file, TString comment1, int histColor, bool exp, bool obs, Float_t yMin_ = 0.0, Float_t yMax_ = 0.10){

  TFile f(dir+"/"+file,"READ");
  if(f.IsZombie() || f.IsZombie()){
    cout << "Cannot open file "<< endl;
    continue;
  }
  TGraphAsymmErrors* g1;
  if(exp) g1 = (TGraphAsymmErrors*)f.Get("expected");
  else if(obs) g1 = (TGraphAsymmErrors*)f.Get("observed");
  g1->SetMarkerColor(kRed);
  g1->SetMarkerStyle(1);
  //g1->SetMarkerSize(1.5);
  g1->SetLineColor(histColor);
  g1->SetLineWidth(3);

  g1->GetXaxis()->SetLimits(85,165);
  g1->GetYaxis()->SetTitleOffset(1.33);
  g1->SetMinimum(yMin_);
  g1->SetMaximum(yMax_);
  g1->GetXaxis()->SetTitle("m_{H^{+}} (GeV)");
  g1->GetYaxis()->SetTitle("95% CL limit for BR(t#rightarrow bH^{#pm})");
  cout<<comment1<<endl;
  return g1;
}

void compareLimit_13TeV( ){
  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  
  TLegend* leg = new TLegend(0.10,0.60,0.40,0.85,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  
  TGraphAsymmErrors* g11 = getLimit("cutSets_limits", "_comb_mjj_kfit_limit.root",       "ele+jets: mjj_kfit_Inc", 1, true, true, 0.0, 0.03);
  TGraphAsymmErrors* g12 = getLimit("cutSets_limits", "_comb_mjj_kfit_pt_bjetH_limit.root", "ele+jets: mjj_kfitL", 2, true, true, 0.0, 0.03);
  TGraphAsymmErrors* g13 = getLimit("cutSets_limits", "_comb_mjj_kfit_CTagL_limit.root", "mu+jets: mjj_kfitT", 3, true, true, 0.0, 0.03);
  TGraphAsymmErrors* g14 = getLimit("cutSets_limits", "_comb_mjj_kfit_CTagM_limit.root", "mu+jets: mjj_kfitT", 4, true, true, 0.0, 0.03);
  TGraphAsymmErrors* g15 = getLimit("cutSets_limits", "_comb_mjj_kfit_CTagT_limit.root", "mu+jets: mjj_kfitT", 6, true, true, 0.0, 0.03);
  
  c1->cd();
  g11->Draw("ALP");
  g12->Draw("LPsame");
  g13->Draw("LPsame");
  g14->Draw("LPsame");
  g15->Draw("LPsame");
  gPad->Modified();
  //leg->SetHeader(Form("#splitline{CMS Preliminary #sqrt{s}=13 TeV}{35.5 fb^{-1}}"));
  
  leg->AddEntry(g11,"lepton+jets: mjj_kfit_Inc","PL");
  leg->AddEntry(g12,"lepton+jets: mjj_kfit_Inc + pt_bjetH","PL");
  leg->AddEntry(g13,"lepton+jets: mjj_kfit_CTagL","PL");
  leg->AddEntry(g14,"lepton+jets: mjj_kfit_CTagM","PL");
  leg->AddEntry(g15,"lepton+jets: mjj_kfit_CTagT","PL");
  leg->Draw();

  TPaveText *pl2 = new TPaveText(0.63,0.75,0.83,0.88, "brNDC");
  pl2->SetTextSize(0.032);
  pl2->SetFillColor(0);
  pl2->SetTextFont(132);
  pl2->SetBorderSize(0);
  pl2->SetTextAlign(11);
  pl2->AddText("t #rightarrow H^{#pm}b, H^{+} #rightarrow c#bar{s}");
  pl2->AddText("BR(H^{+} #rightarrow c#bar{s}) = 1");
  if(exp)pl2->AddText("Expected Limit");
  else if(obs)pl2->AddText("Observed Limit");
  //pl2->Draw();
  gPad->RedrawAxis();
  //c1->SaveAs("comb_overlay.png");
  //c1->Close();
}
