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
TFile* fitDiagOut = TFile::Open("fitDiagnostics.root");

double getStatUnc(TH1F* hCentral, double sError = 0.0){
  double  norm = hCentral->IntegralAndError(1, hCentral->GetNbinsX(), sError);
  //double statUnc = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00;
  double statUnc = sError;
  //2nd method to get stat unc
  double err2 = 0.0;
  for(int i = 1; i<hCentral->GetNbinsX()+1; i++){
    double binErr = hCentral->GetBinError(i);
    //cout<<"bin = "<<i<<endl;
    //cout<<"binCon = "<<hist->GetBinContent(i)<<endl;
    //cout<<"binErr = "<<hist->GetBinError(i)<<endl;
    err2 = err2 + binErr*binErr;
  }
  //statUnc = sqrt(err2);
  return statUnc;
}
void printYeild(TFile *histFile, TString chName, TString histName){
  TH1F* hPreFit = (TH1F*)(histFile->Get("shapes_prefit/hcs_13TeV/"+histName));
  double iPreFit = hPreFit->Integral();
  double ePreFit = 0.0;
  ePreFit = getStatUnc(hPreFit, ePreFit);

  TH1F* hPostFit_s = (TH1F*)(histFile->Get("shapes_fit_s/hcs_13TeV/"+histName));
  double iPostFit_s = hPostFit_s->Integral();
  double ePostFit_s = 0.0;
  ePostFit_s = getStatUnc(hPostFit_s, ePostFit_s);

  TH1F* hPostFit_b = (TH1F*)(histFile->Get("shapes_fit_b/hcs_13TeV/"+histName));
  double iPostFit_b = hPostFit_b->Integral();
  double ePostFit_b = 0.0;
  ePostFit_b = getStatUnc(hPostFit_b, ePostFit_b);
  string histName_(histName);
  printf("%s\t: %.2f +/- %.2f\t %.2f +/- %.2f\t %0.2f +/- %.2f \n", histName_.c_str(), iPreFit, ePreFit, iPostFit_s, ePostFit_s, iPostFit_b, ePostFit_b);
}

//--------------------------------------------//
void MyPrintPostFitYield(){
  printf("%s %s\t %s\t %s\n", "Process", "Prefit", "S+B", "B-only");
  printYeild(fitDiagOut, "", "ttbar");
  printYeild(fitDiagOut, "", "stop");
  printYeild(fitDiagOut, "", "qcd");
  printYeild(fitDiagOut, "", "wjet");
  printYeild(fitDiagOut, "", "zjet");
  printYeild(fitDiagOut, "", "vv");
} 


