#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <algorithm> 

using namespace std;

class MyHPlusDataCardMaker{
  public:
      
  TH1F* getHisto(TFile *inRootFile, TString histPath, TString histName, TFile* fTT);
  TH1F* readWriteHisto(TFile *inFile, TString histPath, TString inHistName, double sf, TFile *outFile, TFile *fTT, TString outHistName,  bool isWrite = false, double min_thres = 0, bool isNeffThreshold = false);
  double getBTagUnc(TH1F *hCentral, TH1F* hUp, TH1F* hDown);
  double getStatUnc(TH1F* hCentral, double sError = 0.0);
  
  private:
  double dont_use;  
};

//----------------------------------------//
//Variuos functions
//----------------------------------------//
TH1F*  MyHPlusDataCardMaker:: getHisto(TFile *inRootFile, TString histPath, TString histName, TFile *fTT){
  TH1F* hist;
  TString fullPath = histPath+histName;
  string exception_msg (inRootFile->GetName()+TString("/")+fullPath+", does not exist");
  try{
    if(!(inRootFile->Get(fullPath)))
       throw  exception_msg.c_str();
  }catch (const char *e){
    ///cout<<"WARNING:"<<e<<endl;
  }
  try{
    if(!(fTT->Get(fullPath)))
       throw  exception_msg.c_str();
  }catch (const char *e){
    cout<<"\033[01;31mERROR: \033[0m"<<e<< endl;
    exit(0);
  }
  if(!(inRootFile->Get(fullPath))){
    hist = (TH1F*)(fTT->Get(fullPath))->Clone(histName);
    hist->Reset();
  }else hist = (TH1F*)(inRootFile->Get(fullPath))->Clone(histName);
  return hist;
}

//Read histos from input file. Write to another file.
TH1F* MyHPlusDataCardMaker::readWriteHisto(TFile *inFile, TString histPath, TString inHistName, double sf, TFile *outFile, TFile *fTT, TString outHistName,  bool isWrite = false, double min_thres = 0, bool isNeffThreshold = false){
  TH1F* hist = (TH1F*) getHisto(inFile, histPath, inHistName, fTT)->Clone(outHistName);
  hist->Scale(sf);
  hist->GetXaxis()->SetRangeUser(0, 200);
  //---------------------------------------//
  // for bin-by-bin uncertainity
  //---------------------------------------//
  if(isNeffThreshold){
    double nentry_before = hist->GetEntries();
    double int_before = hist->Integral();
    //cout<<inHistName<<setw(15)<<"Before: "<<setw(10)<<outHistName<<setw(10)<<hist->GetEntries()<<setw(10)<<hist->Integral()<<setw(10)<<endl;
    double a1_before, a1_after, neff_before, neff_after;
    a1_before = 0; a1_after = 0; neff_before = 0; neff_after = 0;
    for(int ibin=1; ibin<41; ibin++){
      if(hist->GetBinContent(ibin) >0){
        float a1     = hist->GetBinContent(ibin);
        float a1_err = hist->GetBinError(ibin);
        double a2  = a1/a1_err;
        double neff = a2*a2;      //Effective number of events
        a1_before   = a1_before   + a1;
        neff_before = neff_before + neff;
        //cout<<ibin<<"\t"<<neff<<endl;
        if(neff<min_thres) {a1=0; a1_err=0; neff=0;}
        a1_after   = a1_after  + a1;
        neff_after = neff_after + neff;
        hist->SetBinContent(ibin, a1);
        hist->SetBinError(ibin, a1_err);
      }
    }
    //cout<<inHistName<<setw(15)<<"Before: "<<setw(10)<<outHistName<<setw(10)<<nentry_before<<setw(10)<<int_before<<setw(10)<<a1_before<<setw(10)<<neff_before<<setw(10)<<endl;
    //cout<<inHistName<<setw(15)<<"After: "<<setw(10)<<outHistName<<setw(10)<<hist->GetEntries()<<setw(10)<<hist->Integral()<<setw(10)<<a1_after<<setw(10)<<neff_after<<endl;
  }
  if(isWrite){
    outFile->cd();
    hist->Write(outHistName);
  }
  return hist;
}  

//get normalised uncertainity
double MyHPlusDataCardMaker::getBTagUnc(TH1F *hCentral, TH1F* hUp, TH1F* hDown){
  return 1 + max(fabs(hUp->Integral() - hCentral->Integral()), fabs(hCentral->Integral() - hDown->Integral()))/hCentral->Integral();
}

//get statistical uncertainity
double MyHPlusDataCardMaker::getStatUnc(TH1F* hCentral, double sError = 0.0){
  double  norm = hCentral->IntegralAndError(1, hCentral->GetNbinsX(), sError);
  double statUnc = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00; 
  return statUnc;
}

