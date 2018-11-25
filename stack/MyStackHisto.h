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
#include <sys/stat.h>
#include <stdlib.h>
#include <algorithm>
#include <iterator>

using namespace std;

//R. K. Verma
//Sat Jul 14 14:47:08 IST 2018

//INPUT FILES
TFile* fData = TFile::Open("all_muData.root");
//TFile* fData = TFile::Open("all_EleData.root");

TFile* fVV	= TFile::Open("all_VV.root");
TFile* fDY	= TFile::Open("all_DY.root");
TFile* fWJ	= TFile::Open("all_WJets.root");
TFile* fQCD	= TFile::Open("all_QCD.root");
TFile* fST	= TFile::Open("all_ST.root");
TFile* fTT	= TFile::Open("all_TTJetsP.root");
TFile *fSig     = TFile::Open("all_Hplus120.root");
TFile *fSig90     = TFile::Open("all_Hplus90.root");

class MyStackHisto{
  public : 
	//get histogram from root file. Return empty hist, if the hist does not exist.
	TH1F* getHisto(TFile *inRootFile, TString baseDir, TString isoDir, TString histDir, TString histName);
        //decorate histogram
        TH1F* decorateHisto(TH1F* hist, TString myTit, TString xTit, TString yTit);
	//function to stack histos
        void stackHisto(TFile *inRootFile, TString lable, TString baseDir, TString isoDir, TString histDir, TString histName, int color, double scale, THStack* MuptStack, TH1F* hMC, TLegend* leg);
	//qcd from data
	TH1F* getDataDrivenQCD(TString baseDir, TString isoDir, TString histDir, TString histName, double qcd_sf=1.0, double qcd_sf_err = 0.0);
	//function to add histograms
	TH1F*  addHistoForUnc(TString baseDir, TString isoDir, TString histDir, TString histName, bool isDataDrivenQCD = false);
	//Up/down error in the unc band
	double errBandUp(int iBin, TH1F *hCentral, TH1F *hJESPlus, TH1F *hJERPlus, TH1F *bTagPlus, TH1F *cTagPlus, TH1F *PileupPlus, TH1F* hQCD_dd);
	double errBandDown(int iBin, TH1F *hCentral, TH1F *hJESMinus, TH1F *hJERMinus, TH1F *bTagMinus, TH1F *cTagMinus, TH1F* PileupMinus, TH1F* hQCD_dd);
	//unc graph
	TGraphAsymmErrors *UNCGRAPH(TH1F *hCentral, TH1F *hJESPlus, TH1F *hJESMinus, TH1F *hJERPlus, TH1F *hJERMinus, TH1F *bTagPlus, TH1F *bTagMinus, TH1F *cTagPlus, TH1F *cTagMinus, TH1F *PileupPlus, TH1F* PileupMinus, TH1F* hQCD_dd, bool isFullGraph = false, bool isRatioGraph = false);
TPaveText *paveText(double minX, double minY, double maxX, double maxY, int lineColor, int fillColor, int size, int style, int font );

  private :
	int dont_use ;
};

//--------------------------------------------//
//define various functions
//--------------------------------------------//
TH1F*  MyStackHisto:: getHisto(TFile *inRootFile, TString baseDir, TString isoDir, TString histDir, TString histName){
  TH1F* hist;
  TString fullPath = baseDir+isoDir+histDir+histName;
  string exception_msg ("The histogram path, "+fullPath+", does not exist"); 
  try{
    if(!(inRootFile->Get(fullPath)))
       throw  exception_msg.c_str();
  }catch (const char *e){
    cout<<"WARNING:"<<e<<endl;
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
    hist->Add(hist, -1);
  }else hist = (TH1F*)(inRootFile->Get(fullPath))->Clone(histName);
  return hist;
}

TH1F* MyStackHisto:: decorateHisto(TH1F* hist, TString myTit, TString xTit, TString yTit){
  hist->SetTitle(myTit);
  hist->GetXaxis()->SetTitle(xTit);
  hist->GetYaxis()->SetTitle(yTit);
  hist->GetYaxis()->SetTitleOffset(1.00);
  hist->GetXaxis()->SetTitleOffset(1.00);
  hist->GetYaxis()->SetTitleSize(0.10);   
  hist->GetXaxis()->SetTitleSize(0.10);
  hist->GetXaxis()->SetLabelSize(0.10);   
  hist->GetYaxis()->SetLabelSize(0.10);   
  hist->GetXaxis()->SetTickLength(0.05);
  hist->GetXaxis()->SetNdivisions(10);
  hist->GetYaxis()->SetNdivisions(5);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTickLength(0.04);
  return hist;
}

void  MyStackHisto:: stackHisto(TFile *inRootFile, TString lable, TString baseDir, TString isoDir, TString histDir, TString histName, int color, double scale, THStack* MuptStack, TH1F* hMC, TLegend* leg){
  TH1F* hist = getHisto(inRootFile, baseDir, isoDir, histDir, histName);
  hist->Scale(scale);  
  hist->SetFillColor(color);
  leg->AddEntry(hist,lable,"F");
  MuptStack->Add(hist);
  hMC->Add(hist);
}

double MyStackHisto:: errBandUp(int iBin, TH1F *hCentral, TH1F *hJESPlus, TH1F *hJERPlus, TH1F *bTagPlus, TH1F *cTagPlus, TH1F* PileupPlus, TH1F *hQCD_dd){
  double errUp = sqrt(pow(fabs(hJESPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hJERPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(fabs(bTagPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(fabs(cTagPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(fabs(PileupPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(hCentral->GetBinError(iBin+1),2)+ pow(hQCD_dd->GetBinError(iBin+1),2));
  return errUp;
}

double MyStackHisto:: errBandDown(int iBin, TH1F *hCentral, TH1F *hJESMinus, TH1F *hJERMinus, TH1F *bTagMinus, TH1F *cTagMinus, TH1F* PileupMinus, TH1F *hQCD_dd){
  double errDown =sqrt(pow(fabs(hCentral->GetBinContent(iBin+1) - hJESMinus->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hCentral->GetBinContent(iBin+1) - hJERMinus->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hCentral->GetBinContent(iBin+1) - bTagMinus->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hCentral->GetBinContent(iBin+1) - cTagMinus->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hCentral->GetBinContent(iBin+1) - PileupMinus->GetBinContent(iBin+1)),2) + 
		  pow(hCentral->GetBinError(iBin+1),2)+pow(hQCD_dd->GetBinError(iBin+1),2));
  return errDown;
}

TGraphAsymmErrors * MyStackHisto:: UNCGRAPH(TH1F *hCentral, TH1F *hJESPlus, TH1F *hJESMinus, TH1F *hJERPlus, TH1F *hJERMinus, TH1F *bTagPlus, TH1F *bTagMinus, TH1F *cTagPlus, TH1F *cTagMinus, TH1F* PileupPlus, TH1F* PileupMinus, TH1F* hQCD_dd, bool isFullGraph = false, bool isRatioGraph = false){
  TGraphAsymmErrors *gr;
  int n1 = hCentral->GetNbinsX(); 
  double *Yval, *errorU, *errorD, *XerrorU, *XerrorD, *Xval ;
  Yval = new double[n1]; errorU = new double[n1]; errorD = new double[n1];
  XerrorU=new double[n1]; XerrorD=new double[n1]; Xval=new double[n1];
  //cout << "No. of bins= " << n1 << endl;
  for(int i=0; i<n1; i++){
    if(isFullGraph){
    Yval[i]   = hCentral->GetBinContent(i+1);
    errorU[i] = errBandUp(i, hCentral, hJESPlus, hJERPlus, bTagPlus, cTagPlus, PileupPlus, hQCD_dd); 
    errorD[i] = errBandDown(i, hCentral, hJESMinus, hJERMinus, bTagMinus, cTagMinus, PileupMinus, hQCD_dd); 
    //cout<<"bin = "<<i<<endl;
    cout<<i<<"\t"<<Yval[i]<<"\t"<<errorU[i]<<"\t"<<errorD[i]<<endl;
    }
    if(isRatioGraph){
    Yval[i]   = 1;
    errorU[i] = errBandUp(i, hCentral, hJESPlus, hJERPlus, bTagPlus, cTagPlus, PileupPlus, hQCD_dd); 
    errorD[i] = errBandDown(i, hCentral, hJESMinus, hJERMinus, bTagMinus, cTagMinus, PileupMinus, hQCD_dd); 
    errorU[i] = errorU[i]/hCentral->GetBinContent(i+1);
    errorD[i] = errorD[i]/hCentral->GetBinContent(i+1);
    //cout<<Yval[i]<<"\t"<<errorU[i]<<"\t"<<hCentral->GetBinContent(i+1)<<endl;
    }
    Xval[i]   = hCentral->GetBinCenter(i+1);
    XerrorU[i]= hCentral->GetBinWidth(i+1)/2;
    XerrorD[i]= hCentral->GetBinWidth(i+1)/2;
  }
  gr = new TGraphAsymmErrors(n1, Xval, Yval, XerrorD, XerrorU, errorD, errorU);
  return gr;
  delete [] Yval; delete [] errorU; delete [] errorD; delete [] XerrorU; delete [] XerrorD; delete [] Xval;
}

TH1F* MyStackHisto:: getDataDrivenQCD(TString baseDir, TString isoDir, TString histDir, TString histName, double qcd_sf=1.0, double qcd_sf_err = 0.0){
  TH1F* hVV =   getHisto(fVV,   baseDir, "NonIso/", histDir, histName);
  TH1F* hDY =   getHisto(fDY,   baseDir, "NonIso/", histDir, histName); 
  TH1F* hST =   getHisto(fST,   baseDir, "NonIso/", histDir, histName); 
  TH1F* hWJ =   getHisto(fWJ,   baseDir, "NonIso/", histDir, histName); 
  TH1F* hTT =   getHisto(fTT,   baseDir, "NonIso/", histDir, histName); 
  TH1F* hData = getHisto(fData, baseDir, "NonIso/", histDir, histName); 
  TH1F* hOtherMC = (TH1F*)hVV->Clone("hOtherMC"); 
  hOtherMC->Add(hDY); 
  hOtherMC->Add(hST); 
  hOtherMC->Add(hWJ); 
  hOtherMC->Add(hTT); 
  TH1F* hQCD = (TH1F*)hData->Clone(histName); 
  hQCD->Add(hOtherMC, -1);
  cout<<"-------------------------------------"<<endl;
  cout<<histDir<<"/"<<histName<<endl;
  double sError = 0.0;
  double  norm = hQCD->IntegralAndError(1, hQCD->GetNbinsX(), sError);
  cout<<"intB = "<<norm<<", intB_err = "<<sError<<endl;
  cout<<"qcdSF = "<<qcd_sf<<", qcdSF_err = "<<qcd_sf_err<<endl;
  double tot_bin_cont = 0.0; 
  double tot_bin_err = 0.0; 
  ///cout<<"bin"<<setw(10)<<"bin_cont"<<setw(15)<<"bin_err"<<setw(15)<<"new_bin_cont"<<setw(15)<<"new_bin_err"<<endl;
  for(int ibin=1; ibin<41; ibin++){
  //for(int ibin=1; ibin<hQCD->GetNbinsX(); ibin++){
      double bin_cont = hQCD->GetBinContent(ibin);
      double bin_err = hQCD->GetBinError(ibin);
      double new_bin_cont = qcd_sf*bin_cont;
      double new_bin_err = sqrt(pow(bin_cont*qcd_sf_err, 2) + pow(bin_err* qcd_sf, 2));
      ///cout<<ibin<<setw(10)<<bin_cont<<setw(15)<<bin_err<<setw(15)<<new_bin_cont<<setw(15)<<new_bin_err<<endl;
      tot_bin_cont = tot_bin_cont + new_bin_cont;
      tot_bin_err = tot_bin_err + new_bin_err*new_bin_err;

      hQCD->SetBinContent(ibin, new_bin_cont);
      hQCD->SetBinError(ibin, new_bin_err);
    }
  cout<<"tot_bin_cont= "<<tot_bin_cont<<", tot_bin_err = "<<sqrt(tot_bin_err)<<endl;
  ///hQCD->Scale(qcd_sf);
  return hQCD;
}

TH1F*  MyStackHisto:: addHistoForUnc(TString baseDir, TString isoDir, TString histDir, TString histName, bool isDataDrivenQCD = false){
  TH1F* hVV =   	getHisto(fVV,   baseDir, isoDir, histDir, histName);
  TH1F* hDY =   	getHisto(fDY,    baseDir, isoDir, histDir, histName);
  TH1F* hQCD_mc =   	getHisto(fQCD,  baseDir, isoDir, histDir, histName);
  TH1F* hWJ =   	getHisto(fWJ,   baseDir, isoDir, histDir, histName);
  TH1F* hST =   	getHisto(fST,   baseDir, isoDir, histDir, histName);
  TH1F* hTT =   	getHisto(fTT,   baseDir, isoDir, histDir, histName);
  TH1F* hAll = (TH1F*)hVV->Clone("hAllMC");
  hAll->Add(hDY);
  hAll->Add(hWJ);
  hAll->Add(hST);
  hAll->Add(hTT);
  if(!isDataDrivenQCD) hAll->Add(hQCD_mc);
  return hAll;
}

TPaveText * MyStackHisto:: paveText(double minX, double minY, double maxX, double maxY, int lineColor, int fillColor, int size, int style, int font ){
  TPaveText *pt = new TPaveText(minX, minY, maxX, maxY, "brNDC"); // good_v1
  pt->SetBorderSize(size);
  pt->SetFillColor(fillColor);
  pt->SetFillStyle(style);
  pt->SetLineColor(lineColor);
  pt->SetTextFont(font);
  return pt;
}
