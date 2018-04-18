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

TString inFileDir="$PWD";
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


TProfile* getTProfile(TFile *histFile, TString histPath, TString dir, TString histName){
  TProfile* hist;
  if(!(histFile->Get(histPath+"/"+dir+"/"+histName))){
    hist = (TProfile*)(fTT->Get(histPath+"/"+dir+"/"+histName))->Clone(histName);
    hist->Add(hist, -1);
  }else hist = (TProfile*)(histFile->Get(histPath+"/"+dir+"/"+histName))->Clone(histName);
  return hist;
}

TProfile* decorateTProfile(TProfile *graph, TString xTitle, TString yTitle, TString myTitle, double yMin, double yMax, int color){
  graph->SetTitle(myTitle);
  graph->GetYaxis()->SetTitle(yTitle);
  graph->GetYaxis()->SetRangeUser(yMin, yMax);
  graph->GetXaxis()->SetTitle(xTitle);
  graph->GetXaxis()->SetTitleOffset(1.15);
  graph->GetYaxis()->SetTitleSize(0.05);   graph->GetXaxis()->SetTitleSize(0.05);
  graph->GetXaxis()->SetLabelSize(0.05);   graph->GetXaxis()->LabelsOption("u"); // extra
  graph->GetYaxis()->SetLabelSize(0.05);   graph->GetXaxis()->LabelsOption("u"); // extra
  graph->GetXaxis()->SetTickLength(0.05); 
  graph->GetYaxis()->SetTickLength(0.04); 
  graph->GetYaxis()->SetTitleSize(0.06);   
  graph->GetYaxis()->SetTitleOffset(1.05);
  graph->SetLineColor(color);
  graph->SetLineWidth(4);
  graph->SetMarkerSize(20);
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(2);
  graph->SetMarkerColor(color);
  return graph; 
}

void plotTProfile(){
  gStyle->SetOptStat(11111111);
  gStyle->SetFrameLineWidth(4);
  auto c1 = new TCanvas("c1","Profile histogram example",200,10,700,500);
  c1->Divide(3,2);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
 
  TLegend* leg = new TLegend(0.50,0.60,0.70,0.85,NULL,"brNDC");
  TString bin_range = "";
  
  c1->cd(1); 
  TLegend* leg_tt_inc = new TLegend(0.50,0.60,0.70,0.85,NULL,"brNDC");
  TProfile *prof_tt_inc = getGraphForProc(fTT, "PtbJetInc");	
  decorateTProfile(prof_tt_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 0, 200, 1);
  prof_tt_inc->Draw("ALP");
  leg_tt_inc->AddEntry(prof_tt_inc, "ttbar","PL");
  leg_tt_inc->Draw();
 
  c1->cd(2);
  TLegend* leg_wh80_inc = new TLegend(0.50,0.60,0.70,0.85,NULL,"brNDC");
  TProfile *prof_wh80_inc = getGraphForProc(fWH80, "PtbJetInc");	
  decorateTProfile(prof_wh80_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 0, 200, 2);
  prof_wh80_inc->Draw("ALP");
  leg_wh80_inc->AddEntry(prof_wh80_inc, "WH80","PL");
  leg_wh80_inc->Draw();

  c1->cd(3);
  TLegend* leg_wh100_inc = new TLegend(0.50,0.60,0.70,0.85,NULL,"brNDC");
  TProfile *prof_wh100_inc = getGraphForProc(fWH100, "PtbJetInc");	
  decorateTProfile(prof_wh100_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 0, 200, 3);
  prof_wh100_inc->Draw("ALP");
  leg_wh100_inc->AddEntry(prof_wh100_inc, "WH100","PL");
  leg_wh100_inc->Draw();

  c1->cd(4);
  TLegend* leg_wh120_inc = new TLegend(0.50,0.60,0.70,0.85,NULL,"brNDC");
  TProfile *prof_wh120_inc = getGraphForProc(fWH120, "PtbJetInc");	
  decorateTProfile(prof_wh120_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 0, 200, 4);
  prof_wh120_inc->Draw("ALP");
  leg_wh120_inc->AddEntry(prof_wh120_inc, "WH120","PL");
  leg_wh120_inc->Draw();
  
  c1->cd(5);
  TLegend* leg_wh150_inc = new TLegend(0.50,0.60,0.70,0.85,NULL,"brNDC");
  TProfile *prof_wh150_inc = getGraphForProc(fWH150, "PtbJetInc");	
  decorateTProfile(prof_wh150_inc, "Pt_{bjet}^{Had}(GeV)", "Mean of M_{jj}^{Inc} (GeV)", bin_range, 0, 200, 6);
  prof_wh150_inc->Draw("ALP");
  leg_wh150_inc->AddEntry(prof_wh150_inc, "WH150","PL");
  leg_wh150_inc->Draw();
  
  c1->cd(6);
  prof_tt_inc->Draw("ALP");
  prof_wh80_inc->Draw("LPSame");
  prof_wh100_inc->Draw("LPSame");
  prof_wh120_inc->Draw("LPSame");
  prof_wh150_inc->Draw("LPSame");
  leg->AddEntry(prof_tt_inc, "ttbar","PL");
  leg->AddEntry(prof_wh80_inc, "WH80","PL");
  leg->AddEntry(prof_wh100_inc, "WH100","PL");
  leg->AddEntry(prof_wh120_inc, "WH120","PL");
  leg->AddEntry(prof_wh150_inc, "WH150","PL");
  leg->Draw();
}
