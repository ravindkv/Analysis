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
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.5);
  g1->SetLineColor(histColor);
  g1->SetLineWidth(3);
  //g1->SetTitle(comment1);
  g1->GetXaxis()->SetLimits(85,165);
  g1->GetYaxis()->SetTitleOffset(1.53);
  g1->SetMinimum(yMin_);
  g1->SetMaximum(yMax_);
  g1->GetXaxis()->SetTitle("m_{H^{+}} (GeV)");
  //g1->GetYaxis()->SetTitle("95% CL limit for BR(t#rightarrow bH^{#pm})");
  g1->SetTitle("95% CL limit for BR(t#rightarrow bH^{#pm})");
  cout<<comment1<<endl;
  return g1;
}

void compareLimit_13TeV( ){
  TCanvas *c1 = new TCanvas("c1","",5,30,950,900);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->Divide(2,2);
  
  TLegend* leg_mu = new TLegend(0.30,0.40,0.60,0.85,NULL,"brNDC");
  //leg->SetHeader(Form("#splitline{CMS Preliminary #sqrt{s}=13 TeV}{35.5 fb^{-1}}"));
  leg_mu->SetBorderSize(0);
  leg_mu->SetTextSize(0.03);
  leg_mu->SetFillColor(0);
  
  TLegend* leg_ele = new TLegend(0.30,0.40,0.60,0.85,NULL,"brNDC");
  leg_ele->SetBorderSize(0);
  leg_ele->SetTextSize(0.03);
  leg_ele->SetFillColor(0);
  
  TLegend* leg_lep = new TLegend(0.30,0.40,0.60,0.85,NULL,"brNDC");
  leg_lep->SetBorderSize(0);
  leg_lep->SetTextSize(0.03);
  leg_lep->SetFillColor(0);
  
  TLegend* leg_comb = new TLegend(0.30,0.32,0.60,0.90,NULL,"brNDC");
  leg_comb->SetBorderSize(0);
  leg_comb->SetTextSize(0.03);
  leg_comb->SetFillColor(0);
  
  //muon channel 
  TGraphAsymmErrors* g_mu1 = getLimit("../limit/mu/mjj_kfit","limit_mu_mjj_kfit.root", "mu+jets", 1, true, false, 0.0, 0.02);
  TGraphAsymmErrors* g_mu2 = getLimit("mu/mjj_CTagCat","limit_mu_mjj_CTagCat.root", "mu+jets", 2, true, false, 0.0, 0.02);
  TGraphAsymmErrors* g_mu3 = getLimit("mu/mjj_bJetPtCat","limit_mu_mjj_bJetPtCat.root", "mu+jets", 3, true, false, 0.0, 0.02);
  TGraphAsymmErrors* g_mu4 = getLimit("mu/mjj_CTagCat_bJetPtCat","limit_mu_mjj_CTagCat_bJetPtCat.root", "mu+jets", 4, true, false, 0.0, 0.02);
  
  //electron channel 
  TGraphAsymmErrors* g_ele1 = getLimit("../limit/ele/mjj_kfit","limit_ele_mjj_kfit.root", "ele+jets", 6, true, false, 0.0, 0.02);
  TGraphAsymmErrors* g_ele2 = getLimit("ele/mjj_CTagCat","limit_ele_mjj_CTagCat.root", "ele+jets", 7, true, false, 0.0, 0.02);
  TGraphAsymmErrors* g_ele3 = getLimit("ele/mjj_bJetPtCat","limit_ele_mjj_bJetPtCat.root", "ele+jets", 29, true, false, 0.0, 0.02);
  TGraphAsymmErrors* g_ele4 = getLimit("ele/mjj_CTagCat_bJetPtCat","limit_ele_mjj_CTagCat_bJetPtCat.root", "ele+jets", 9, true, false, 0.0, 0.02);
  
  //lepton channel 
  TGraphAsymmErrors* g_lep1 = getLimit("../limit/mu_ele/mjj_kfit","limit_mu_ele_mjj_kfit.root", "mu_ele+jets", 12, true, false, 0.0, 0.02);
  TGraphAsymmErrors* g_lep2 = getLimit("mu_ele/mjj_CTagCat","limit_mu_ele_mjj_CTagCat.root", "mu_ele+jets", 22, true, false, 0.0, 0.02);
  TGraphAsymmErrors* g_lep3 = getLimit("mu_ele/mjj_bJetPtCat","limit_mu_ele_mjj_bJetPtCat.root", "mu_ele+jets", 32, true, false, 0.0, 0.02);
  TGraphAsymmErrors* g_lep4 = getLimit("mu_ele/mjj_CTagCat_bJetPtCat","limit_mu_ele_mjj_CTagCat_bJetPtCat.root", "mu_ele+jets", 42, true, false, 0.0, 0.02);
  
  
  c1->cd(1);
  g_mu1->GetYaxis()->SetTitleOffset(1.73);
  g_mu1->Draw("ALP");
  g_mu2->Draw("LPsame");
  g_mu3->Draw("LPsame");
  g_mu4->Draw("LPsame");
  gPad->Modified();
  leg_mu->AddEntry(g_mu1,"mu+jets: M_{jj}","PL");
  leg_mu->AddEntry(g_mu2,"mu+jets: M_{jj}(Pt_{bjet}^{Had} bins)","PL");
  leg_mu->AddEntry(g_mu3,"mu+jets: M_{jj}^{CatL} + M_{jj}^{CatM}+ M_{jj}^{CatT}","PL");
  leg_mu->AddEntry(g_mu4,"mu+jets: M_{jj}^{CatL}(Pt_{bjet}^{Had} bins) + M_{jj}^{CatM}(Pt_{bjet}^{Had} bins)+ M_{jj}^{CatT}(Pt_{bjet}^{Had} bins)","PL");
  leg_mu->Draw();
  c1->cd(2);
  g_ele1->GetYaxis()->SetTitleOffset(1.73);
  g_ele1->Draw("ALP");
  g_ele2->Draw("LPsame");
  g_ele3->Draw("LPsame");
  g_ele4->Draw("LPsame");
  gPad->Modified();
  leg_ele->AddEntry(g_ele1,"ele+jets: M_{jj}","PL");
  leg_ele->AddEntry(g_ele2,"ele+jets: M_{jj}(Pt_{bjet}^{Had} bins)","PL");
  leg_ele->AddEntry(g_ele3,"ele+jets: M_{jj}^{CatL} + M_{jj}^{CatM}+ M_{jj}^{CatT}","PL");
  leg_ele->AddEntry(g_ele4,"ele+jets: M_{jj}^{CatL}(Pt_{bjet}^{Had} bins) + M_{jj}^{CatM}(Pt_{bjet}^{Had} bins)+ M_{jj}^{CatT}(Pt_{bjet}^{Had} bins)","PL");
  leg_ele->Draw();

  c1->cd(3);
  g_lep1->GetYaxis()->SetTitleOffset(1.73);
  g_lep1->Draw("ALP");
  g_lep2->Draw("LPsame");
  g_lep3->Draw("LPsame");
  g_lep4->Draw("LPsame");
  gPad->Modified();
  leg_lep->AddEntry(g_lep1,"lep+jets: M_{jj}","PL");
  leg_lep->AddEntry(g_lep2,"lep+jets: M_{jj}(Pt_{bjet}^{Had} bins)","PL");
  leg_lep->AddEntry(g_lep3,"lep+jets: M_{jj}^{CatL} + M_{jj}^{CatM}+ M_{jj}^{CatT}","PL");
  leg_lep->AddEntry(g_lep4,"lep+jets: M_{jj}^{CatL}(Pt_{bjet}^{Had} bins) + M_{jj}^{CatM}(Pt_{bjet}^{Had} bins)+ M_{jj}^{CatT}(Pt_{bjet}^{Had} bins)","PL");
  leg_lep->Draw();
  
  c1->cd(4);
  g_mu1->Draw("ALP");
  g_mu2->Draw("LPsame");
  g_mu3->Draw("LPsame");
  g_mu4->Draw("LPsame");
  g_ele1->Draw("LPsame");
  g_ele2->Draw("LPsame");
  g_ele3->Draw("LPsame");
  g_ele4->Draw("LPsame");
  g_lep1->Draw("LPsame");
  g_lep2->Draw("LPsame");
  g_lep3->Draw("LPsame");
  g_lep4->Draw("LPsame");
  gPad->Modified();
  leg_comb->AddEntry(g_mu1,"mu+jets: M_{jj}","PL");
  leg_comb->AddEntry(g_mu2,"mu+jets: M_{jj}(Pt_{bjet}^{Had} bins)","PL");
  leg_comb->AddEntry(g_mu3,"mu+jets: M_{jj}^{CatL} + M_{jj}^{CatM}+ M_{jj}^{CatT}","PL");
  leg_comb->AddEntry(g_mu4,"mu+jets: M_{jj}^{CatL}(Pt_{bjet}^{Had} bins) + M_{jj}^{CatM}(Pt_{bjet}^{Had} bins)+ M_{jj}^{CatT}(Pt_{bjet}^{Had} bins)","PL");
  leg_comb->AddEntry(g_ele1,"ele+jets: M_{jj}","PL");
  leg_comb->AddEntry(g_ele2,"ele+jets: M_{jj}(Pt_{bjet}^{Had} bins)","PL");
  leg_comb->AddEntry(g_ele3,"ele+jets: M_{jj}^{CatL} + M_{jj}^{CatM}+ M_{jj}^{CatT}","PL");
  leg_comb->AddEntry(g_ele4,"ele+jets: M_{jj}^{CatL}(Pt_{bjet}^{Had} bins) + M_{jj}^{CatM}(Pt_{bjet}^{Had} bins)+ M_{jj}^{CatT}(Pt_{bjet}^{Had} bins)","PL");
  leg_comb->AddEntry(g_lep1,"lep+jets: M_{jj}","PL");
  leg_comb->AddEntry(g_lep2,"lep+jets: M_{jj}(Pt_{bjet}^{Had} bins)","PL");
  leg_comb->AddEntry(g_lep3,"lep+jets: M_{jj}^{CatL} + M_{jj}^{CatM}+ M_{jj}^{CatT}","PL");
  leg_comb->AddEntry(g_lep4,"lep+jets: M_{jj}^{CatL}(Pt_{bjet}^{Had} bins) + M_{jj}^{CatM}(Pt_{bjet}^{Had} bins)+ M_{jj}^{CatT}(Pt_{bjet}^{Had} bins)","PL");
  leg_comb->Draw();
  
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
  c1->SaveAs("comb_overlay.pdf");
  //c1->Close();
}
