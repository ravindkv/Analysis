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

bool isMuChannel = false;
bool isEleChannel = true;

TString inFileDir="$PWD";
//TFile* fData      = TFile::Open(inFileDir+"/all_muData.root");
TFile* fData      = TFile::Open(inFileDir+"/all_EleData.root");
bool isSaveHisto = true;
//bkg
TFile* fTT      = TFile::Open(inFileDir+"/all_TTJetsP.root");
//signal
TFile *fWH80      = TFile::Open(inFileDir+"/all_Hplus80.root");
TFile *fWH90      = TFile::Open(inFileDir+"/all_Hplus90.root");
TFile *fWH100     = TFile::Open(inFileDir+"/all_Hplus100.root");
TFile *fWH120     = TFile::Open(inFileDir+"/all_Hplus120.root");
TFile *fWH140     = TFile::Open(inFileDir+"/all_Hplus140.root");
TFile *fWH150     = TFile::Open(inFileDir+"/all_Hplus150.root");
TFile *fWH155     = TFile::Open(inFileDir+"/all_Hplus155.root");
TFile *fWH160     = TFile::Open(inFileDir+"/all_Hplus160.root");


TProfile* getTProfile(TFile *profFile, TString profPath, TString dir, TString profName){
  TProfile* prof;
  if(!(profFile->Get(profPath+"/"+dir+"/"+profName))){
    prof = (TProfile*)(fTT->Get(profPath+"/"+dir+"/"+profName))->Clone(profName);
    prof->Add(prof, -1);
  }else prof = (TProfile*)(profFile->Get(profPath+"/"+dir+"/"+profName))->Clone(profName);
  prof->SetErrorOption("s");
  prof->GetXaxis()->SetRangeUser(25, 350);
  return prof;
}
TProfile* getTProfile2(TFile *profFile, TString profPath, TString profName){
  TProfile* prof;
  cout<<profPath+"/"+profName<<endl;
  if(!(profFile->Get(profPath+"/"+profName))){
    prof = (TProfile*)(fTT->Get(profPath+"/"+profName))->Clone(profName);
    prof->Add(prof, -1);
  }else prof = (TProfile*)(profFile->Get(profPath+"/"+profName))->Clone(profName);
  prof->SetErrorOption("s");
  return prof;
}

TProfile* decorateTProfile(TProfile *profile, TString xTitle, TString yTitle, TString myTitle, double yMin, double yMax, int color){
  profile->SetTitle(myTitle);
  profile->GetYaxis()->SetTitle(yTitle);
  profile->GetYaxis()->SetRangeUser(yMin, yMax);
  profile->GetXaxis()->SetRangeUser(25, 350);
  profile->GetXaxis()->SetTitle(xTitle);
  profile->GetXaxis()->SetTitleOffset(1.00);
  profile->GetYaxis()->SetTitleOffset(0.90);
  profile->GetYaxis()->SetTitleSize(0.08);   profile->GetXaxis()->SetTitleSize(0.07);
  profile->GetXaxis()->SetLabelSize(0.06);   profile->GetXaxis()->LabelsOption("u"); // extra
  profile->GetYaxis()->SetLabelSize(0.06);   profile->GetXaxis()->LabelsOption("u"); // extra
  profile->GetXaxis()->SetTickLength(0.05); 
  profile->GetYaxis()->SetTickLength(0.04); 
  profile->SetLineColor(color);
  profile->SetLineWidth(3);
  profile->SetMarkerSize(1);
  profile->SetMarkerStyle(20);
  profile->SetMarkerColor(color);
  return profile; 
}

void plotAllTProfile(){
  gStyle->SetOptStat(0);
  //gStyle->SetOptStat(11111111);
  gStyle->SetFrameLineWidth(3);
  auto c1 = new TCanvas();//"c1","Profile histogram example",200,10,700,500);
  c1->Divide(3,2);
  gPad->SetBottomMargin(0.20);
  gPad->SetLeftMargin(0.20);
  gPad->SetRightMargin(0.05);
 
  TLegend* leg = new TLegend(0.50,0.60,0.70,0.85,NULL,"brNDC");
  leg->SetBorderSize(0);
  TString bin_range = "";
  c1->cd(1); 
  TLegend* leg_tt_inc = new TLegend(0.50,0.60,0.70,0.85,NULL,"brNDC");
  leg_tt_inc->SetBorderSize(0);
  TProfile *prof_tt_inc = getTProfile(fTT, "base/Iso", "PtbJetInc", "mjj_kfit_pt_bjetH");
  decorateTProfile(prof_tt_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 40, 160, 1);
  prof_tt_inc->GetXaxis()->SetRangeUser(25, 350);
  prof_tt_inc->Draw();
  leg_tt_inc->AddEntry(prof_tt_inc, "ttbar","PL");
  leg_tt_inc->Draw();
 
  c1->cd(2);
  TLegend* leg_wh80_inc = new TLegend(0.50,0.60,0.70,0.85,NULL,"brNDC");
  leg_wh80_inc->SetBorderSize(0);
  TProfile *prof_wh80_inc = getTProfile(fWH80, "base/Iso", "PtbJetInc", "mjj_kfit_pt_bjetH");
  decorateTProfile(prof_wh80_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 40, 160, 2);
  prof_wh80_inc->Draw();
  leg_wh80_inc->AddEntry(prof_wh80_inc, "WH80","PL");
  leg_wh80_inc->Draw();

  c1->cd(3);
  TLegend* leg_wh100_inc = new TLegend(0.50,0.60,0.70,0.85,NULL,"brNDC");
  leg_wh100_inc->SetBorderSize(0);
  TProfile *prof_wh100_inc = getTProfile(fWH100, "base/Iso", "PtbJetInc", "mjj_kfit_pt_bjetH");
  decorateTProfile(prof_wh100_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 40, 160, 3);
  prof_wh100_inc->Draw();
  leg_wh100_inc->AddEntry(prof_wh100_inc, "WH100","PL");
  leg_wh100_inc->Draw();

  c1->cd(4);
  TLegend* leg_wh120_inc = new TLegend(0.50,0.60,0.70,0.85,NULL,"brNDC");
  leg_wh120_inc->SetBorderSize(0);
  TProfile *prof_wh120_inc = getTProfile(fWH120, "base/Iso", "PtbJetInc", "mjj_kfit_pt_bjetH");
  decorateTProfile(prof_wh120_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 40, 160, 4);
  prof_wh120_inc->Draw();
  leg_wh120_inc->AddEntry(prof_wh120_inc, "WH120","PL");
  leg_wh120_inc->Draw();
  
  c1->cd(5);
  TLegend* leg_wh150_inc = new TLegend(0.50,0.60,0.70,0.85,NULL,"brNDC");
  leg_wh150_inc->SetBorderSize(0);
  TProfile *prof_wh150_inc = getTProfile(fWH150, "base/Iso", "PtbJetInc", "mjj_kfit_pt_bjetH");
  decorateTProfile(prof_wh150_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 40, 160, 6);
  prof_wh150_inc->Draw();
  leg_wh150_inc->AddEntry(prof_wh150_inc, "WH150","PL");
  leg_wh150_inc->Draw();
  
  c1->cd(6);
  prof_tt_inc->Draw();
  prof_wh80_inc->Draw("Same");
  prof_wh100_inc->Draw("Same");
  prof_wh120_inc->Draw("Same");
  prof_wh150_inc->Draw("Same");
  leg->AddEntry(prof_tt_inc, "ttbar","PL");
  leg->AddEntry(prof_wh80_inc, "WH80","PL");
  leg->AddEntry(prof_wh100_inc, "WH100","PL");
  leg->AddEntry(prof_wh120_inc, "WH120","PL");
  leg->AddEntry(prof_wh150_inc, "WH150","PL");
  leg->Draw();
}
void overLayAllTProfile(){
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(3);
  auto c1 = new TCanvas();//"c1","Profile histogram example",200,10,700,500);
  gPad->SetBottomMargin(0.20);
  gPad->SetLeftMargin(0.20);
  gPad->SetRightMargin(0.05);
 
  TLegend* leg = new TLegend(0.50,0.60,0.70,0.85,NULL,"brNDC");
  leg->SetBorderSize(0);
  TString bin_range = "";
  TProfile *prof_tt_inc = getTProfile(fTT, "base/Iso", "PtbJetInc", "mjj_kfit_pt_bjetH");
  decorateTProfile(prof_tt_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 40, 160, 1);
 
  TProfile *prof_wh80_inc = getTProfile(fWH80, "base/Iso", "PtbJetInc", "mjj_kfit_pt_bjetH");
  decorateTProfile(prof_wh80_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 40, 160, 2);

  TProfile *prof_wh100_inc = getTProfile(fWH100, "base/Iso", "PtbJetInc", "mjj_kfit_pt_bjetH");
  decorateTProfile(prof_wh100_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 40, 160, 3);

  TProfile *prof_wh120_inc = getTProfile(fWH120, "base/Iso", "PtbJetInc", "mjj_kfit_pt_bjetH");
  decorateTProfile(prof_wh120_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 40, 160, 4);
  
  TProfile *prof_wh150_inc = getTProfile(fWH150, "base/Iso", "PtbJetInc", "mjj_kfit_pt_bjetH");
  decorateTProfile(prof_wh150_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 40, 160, 6);
  
  prof_tt_inc->Draw();
  prof_wh80_inc->Draw("Same");
  prof_wh100_inc->Draw("Same");
  prof_wh120_inc->Draw("Same");
  prof_wh150_inc->Draw("Same");
  leg->AddEntry(prof_tt_inc, "t#bar{t} + jets","PL");
  leg->AddEntry(prof_wh80_inc, "WH80","PL");
  leg->AddEntry(prof_wh100_inc, "WH100","PL");
  leg->AddEntry(prof_wh120_inc, "WH120","PL");
  leg->AddEntry(prof_wh150_inc, "WH150","PL");

  //pave text CMS box
  TPaveText *pt = new TPaveText(0.25,0.9354,0.90,0.9362, "brNDC"); // good_v1
  pt->SetBorderSize(1);
  pt->SetFillColor(19);
  pt->SetFillStyle(0);
  pt->SetTextSize(0.08);
  pt->SetLineColor(0);
  pt->SetTextFont(132);
  TText *text = pt->AddText("#sqrt{s}=13 TeV, 35.9 fb^{-1}; ");
  text->SetTextAlign(11);
  
  //pave text channel box
  TPaveText *ch = new TPaveText(1.00,0.9154898,0.7510067,0.9762187,"brNDC");
  ch->SetFillColor(19);
  ch->SetFillStyle(0);
  ch->SetLineColor(0);
  ch->SetTextSize(0.08);
  ch->SetBorderSize(1);
  if(isMuChannel) ch->AddText("#mu + jets");
  if(isEleChannel) ch->AddText("e + jets");
  pt->Draw();
  ch->Draw();
  leg->Draw();
  if(isMuChannel) c1->SaveAs("TProf_mu.pdf");
  if(isEleChannel) c1->SaveAs("TProf_ele.pdf");
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
void plotTProfile(){
  gStyle->SetFrameLineWidth(3);
  //pave text CMS box
  TPaveText *pt = paveText(0.09,0.9354,0.60,0.9362, 0, 19, 1, 0, 132);
  pt->SetTextSize(0.059);
  TText *text = pt->AddText("CMS Preliminary, #sqrt{s} = 13 TeV, 35.9 fb^{-1}");
  text->SetTextAlign(11);

  //pave text channel box
  TPaveText *ch = paveText(0.803,0.9154898,0.9010067,0.9762187, 0, 19, 1, 0, 132);
  ch->SetTextSize(0.08);
  if(isMuChannel) ch->AddText("#mu + jets");
  if(isEleChannel) ch->AddText("e + jets");
  gStyle->SetOptStat(0);
  //gStyle->SetOptStat(11111111);
  gStyle->SetFrameLineWidth(3);
  auto c1 = new TCanvas();//"c1","Profile histogram example",200,10,700,500);
  gPad->SetBottomMargin(0.20);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
 
  TLegend* leg = new TLegend(0.50,0.80,0.70,0.85,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.08);
  TString bin_range = "";
  
  TProfile *prof_tt_inc = getTProfile2(fData, "base", "RelIso_MET_TProf");
  if(isEleChannel) decorateTProfile(prof_tt_inc, "E_{T}^{miss} [GeV]", "I_{rel}^{e}", bin_range, -0.15, 0.15, 3);
  if(isMuChannel) decorateTProfile(prof_tt_inc, "E_{T}^{miss} [GeV]", "I_{rel}^{#mu}", bin_range, -0.50, 0.50, 3);
  prof_tt_inc->Draw();
  leg->AddEntry(prof_tt_inc, "data","PL");
  leg->Draw();
  leg->Draw();
  ch->Draw();
  pt->Draw();
  if(isSaveHisto){
    TString outFile("$PWD/");
    if(isMuChannel) outFile += "RelIso_MET_TProf_mu.pdf";
    if(isEleChannel) outFile += "RelIso_MET_TProf_ele.pdf";
    c1->SaveAs(outFile);
    //c1->Close();
  }
}
