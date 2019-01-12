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
TString baseDir = "baseLowMET";

class MyStackHisto{
  public : 
    //get histogram from root file. Return empty hist, if the hist does not exist.
    TH1F* getHisto(TFile *inRootFile, TString baseDir, TString isoDir, TString histDir, TString histName);
    //decorate histogram
    TH1F* decorateHisto(TH1F* hist, TString myTit, TString xTit, TString yTit);
    //function to stack histos
    void stackHisto(TFile *inRootFile, TString lable, TString baseDir, TString isoDir, TString histDir, TString histName, int color, double scale, THStack* MuptStack, TH1F* hMC, TLegend* leg);
    //qcd from data
    vector<double> getQcdSF(TFile* fData, TFile* fTT, TFile* fST, TFile* fWJ, TFile* fDY, TFile* fVV, TString histDir, TString histName);
    TH1F* getDDqcd(TString baseDir, TString isoDir, TString histDir, TString histName, double qcd_sf=1.0, double qcd_sf_err = 0.0);
    double getStatUnc(TH1F* hCentral, double sError = 0.0);
    //function to add histograms
    TH1F*  addHistoForUnc(TString baseDir, TString isoDir, TString histDir, TString histName, bool isDataDrivenQCD = false);
    //Up/down error in the unc band
    double errBandUp(int iBin, TH1F *hCentral, TH1F *hJESPlus, TH1F *hJERPlus, TH1F *bcTagPlus1, TH1F *bcTagPlus2, TH1F *bcTagPlus3, TH1F *PileupPlus, TH1F* hQCD_dd);
    double errBandDown(int iBin, TH1F *hCentral, TH1F *hJESMinus, TH1F *hJERMinus, TH1F *bcTagMinus1, TH1F *bcTagMinus2, TH1F *bcTagMinus3, TH1F* PileupMinus, TH1F* hQCD_dd);
    //unc graph
    TGraphAsymmErrors *UNCGRAPH(TH1F *hCentral, TH1F *hJESPlus, TH1F *hJESMinus, TH1F *hJERPlus, TH1F *hJERMinus, TH1F *bcTagPlus1, TH1F *bcTagMinus1, TH1F *bcTagPlus2, TH1F *bcTagMinus2,TH1F *bcTagPlus3, TH1F *bcTagMinus3, TH1F *PileupPlus, TH1F* PileupMinus, TH1F* hQCD_dd, bool isFullGraph = false, bool isRatioGraph = false);
TPaveText *paveText(double minX, double minY, double maxX, double maxY, int lineColor, int fillColor, int size, int style, int font );

  private :
    int dont_use ;
};

//--------------------------------------------//
//define various functions
//--------------------------------------------//
TH1F*  MyStackHisto:: getHisto(TFile *inRootFile, TString baseDir, TString isoDir, TString histDir, TString histName){
  TH1F* hist;
  TString inFile(inRootFile->GetName());
  TString fullPath = baseDir+isoDir+histDir+histName;
  string exception_msg ("The histogram path, "+inFile+"/"+fullPath+", does not exist"); 
  try{
    if(!(inRootFile->Get(fullPath)))
       throw  exception_msg.c_str();
  }catch (const char *e){
    //cout<<"WARNING:"<<e<<endl;
  }
  try{
    if(!(fTT->Get("base/Iso/"+histDir+histName))) //to initialise an empty hist
       throw  exception_msg.c_str();
  }catch (const char *e){
    cout<<"\033[01;31mERROR: \033[0m"<<e<< endl;
    exit(0);
  }
  if(!(inRootFile->Get(fullPath))){
    hist = (TH1F*)(fTT->Get("base/Iso/"+histDir+histName))->Clone(histName); //to initialise an empty hist
    hist->Reset();
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

double MyStackHisto:: errBandUp(int iBin, TH1F *hCentral, TH1F *hJESPlus, TH1F *hJERPlus, TH1F *bcTagPlus1, TH1F *bcTagPlus2, TH1F *bcTagPlus3, TH1F *PileupPlus, TH1F* hQCD_dd){
  double errUp = sqrt(pow(fabs(hJESPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hJERPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(fabs(bcTagPlus1->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(fabs(bcTagPlus2->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(fabs(bcTagPlus3->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(fabs(PileupPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(hCentral->GetBinError(iBin+1),2)+ pow(hQCD_dd->GetBinError(iBin+1),2));
  return errUp;
}

double MyStackHisto:: errBandDown(int iBin, TH1F *hCentral, TH1F *hJESMinus, TH1F *hJERMinus, TH1F *bcTagMinus1, TH1F *bcTagMinus2, TH1F *bcTagMinus3, TH1F* PileupMinus, TH1F* hQCD_dd){
  double errDown =sqrt(pow(fabs(hCentral->GetBinContent(iBin+1) - hJESMinus->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hCentral->GetBinContent(iBin+1) - hJERMinus->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hCentral->GetBinContent(iBin+1) - bcTagMinus1->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hCentral->GetBinContent(iBin+1) - bcTagMinus2->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hCentral->GetBinContent(iBin+1) - bcTagMinus3->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hCentral->GetBinContent(iBin+1) - PileupMinus->GetBinContent(iBin+1)),2) + 
		  pow(hCentral->GetBinError(iBin+1),2)+pow(hQCD_dd->GetBinError(iBin+1),2));
  return errDown;
}

TGraphAsymmErrors * MyStackHisto::UNCGRAPH(TH1F *hCentral, TH1F *hJESPlus, TH1F *hJESMinus, TH1F *hJERPlus, TH1F *hJERMinus, TH1F *bcTagPlus1, TH1F *bcTagMinus1, TH1F *bcTagPlus2, TH1F *bcTagMinus2,TH1F *bcTagPlus3, TH1F *bcTagMinus3, TH1F *PileupPlus, TH1F* PileupMinus, TH1F* hQCD_dd, bool isFullGraph = false, bool isRatioGraph = false){
  TGraphAsymmErrors *gr;
  int n1 = hCentral->GetNbinsX(); 
  double *Yval, *errorU, *errorD, *XerrorU, *XerrorD, *Xval ;
  Yval = new double[n1]; errorU = new double[n1]; errorD = new double[n1];
  XerrorU=new double[n1]; XerrorD=new double[n1]; Xval=new double[n1];
  //cout << "No. of bins= " << n1 << endl;
  for(int i=0; i<n1; i++){
    if(isFullGraph){
    Yval[i]   = hCentral->GetBinContent(i+1);
    errorU[i] = errBandUp(i, hCentral, hJESPlus, hJERPlus, bcTagPlus1, bcTagPlus2, bcTagPlus3, PileupPlus, hQCD_dd); 
    errorD[i] = errBandDown(i, hCentral, hJESMinus, hJERMinus, bcTagMinus1, bcTagMinus2, bcTagMinus3, PileupMinus, hQCD_dd); 
    }
    if(isRatioGraph){
    Yval[i]   = 1;
    errorU[i] = errBandUp(i, hCentral, hJESPlus, hJERPlus, bcTagPlus1, bcTagPlus2, bcTagPlus3, PileupPlus, hQCD_dd); 
    errorD[i] = errBandDown(i, hCentral, hJESMinus, hJERMinus, bcTagMinus1, bcTagMinus2, bcTagMinus3, PileupMinus, hQCD_dd); 
    //cout<<"bin = "<<i<<endl;
    //cout<<Yval[i]<<"\t"<<errorU[i]<<"\t"<<hCentral->GetBinContent(i+1)<<endl;
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

TH1F* MyStackHisto:: getDDqcd(TString baseDir, TString isoDir, TString histDir, TString histName, double qcd_sf=1.0, double qcd_sf_err = 0.0){
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
double MyStackHisto::getStatUnc(TH1F* hCentral, double sError = 0.0){
  double  norm = hCentral->IntegralAndError(1, hCentral->GetNbinsX(), sError);
  //double statUnc = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00;
  double statUnc = sError;
  return statUnc;
}

vector<double> MyStackHisto::getQcdSF(TFile* fData, TFile* fTT, TFile* fST, TFile* fWJ, TFile* fDY, TFile* fVV, TString histDir, TString histName){
  //RegionC = LowMET, Iso
  TH1F* hVV_RegC = getHisto(fVV,   "baseLowMET/", "NonIso/", histDir, histName);//Reg = Region
  TH1F* hDY_RegC = getHisto(fDY,   "baseLowMET/", "NonIso/", histDir, histName);
  TH1F* hWJ_RegC = getHisto(fWJ,   "baseLowMET/", "NonIso/", histDir, histName);
  TH1F* hST_RegC = getHisto(fST,   "baseLowMET/", "NonIso/", histDir, histName);
  TH1F* hTT_RegC = getHisto(fTT,   "baseLowMET/", "NonIso/", histDir, histName);
  TH1F* hMC_RegC = (TH1F*)hVV_RegC->Clone("hAllMC");
  hMC_RegC->Add(hDY_RegC);
  hMC_RegC->Add(hWJ_RegC);
  hMC_RegC->Add(hST_RegC);
  hMC_RegC->Add(hTT_RegC);
  TH1F* hData_RegC= (TH1F*) getHisto(fData, "baseLowMET/", "NonIso/", histDir, histName);

  //RegionD = LowMET, NonIso
  TH1F* hVV_RegD = getHisto(fVV,   "baseLowMET/", "Iso/", histDir, histName);
  TH1F* hDY_RegD = getHisto(fDY,   "baseLowMET/", "Iso/", histDir, histName);
  TH1F* hWJ_RegD = getHisto(fWJ,   "baseLowMET/", "Iso/", histDir, histName);
  TH1F* hST_RegD = getHisto(fST,   "baseLowMET/", "Iso/", histDir, histName);
  TH1F* hTT_RegD = getHisto(fTT,   "baseLowMET/", "Iso/", histDir, histName);
  TH1F* hMC_RegD = (TH1F*)hVV_RegD->Clone("hAllMC");
  hMC_RegD->Add(hDY_RegD);
  hMC_RegD->Add(hWJ_RegD);
  hMC_RegD->Add(hST_RegD);
  hMC_RegD->Add(hTT_RegD);
  TH1F* hData_RegD=  getHisto(fData, "baseLowMET/", "Iso/", histDir, histName);

  //(Data - MC) from RegionC
  double intMC_RegC   = hMC_RegC->Integral();
  double errMC_RegC   = getStatUnc(hMC_RegC, 0.0);
  double intData_RegC = hData_RegC->Integral();
  double errData_RegC = getStatUnc(hData_RegC, 0.0);
  //double intDiff_RegC = abs(intData_RegC - intMC_RegC);
  double intDiff_RegC = (intData_RegC - intMC_RegC);
  double errDiff_RegC = sqrt(errMC_RegC*errMC_RegC + errData_RegC*errData_RegC);

  //(Data - MC) from RegionD
  double intMC_RegD   = hMC_RegD->Integral();
  double errMC_RegD   = getStatUnc(hMC_RegD, 0.0);
  double intData_RegD = hData_RegD->Integral();
  double errData_RegD = getStatUnc(hData_RegD, 0.0);
  //double intDiff_RegD = abs(intData_RegD - intMC_RegD);
  double intDiff_RegD = (intData_RegD - intMC_RegD);
  double errDiff_RegD = sqrt(errMC_RegD*errMC_RegD + errData_RegD*errData_RegD);
  
  //Ratio of (Data-MC) from RegionD and RegionC
  double ratioDiffRegDC = intDiff_RegD/intDiff_RegC;
  double tmp_RegD = errDiff_RegD/intDiff_RegD;
  double tmp_RegC = errDiff_RegC/intDiff_RegC;
  double errDiffRegDC = ratioDiffRegDC*sqrt(tmp_RegD*tmp_RegD + tmp_RegC*tmp_RegC); 

  //QCD scale factor and error
  bool isRealSF = false;
  //if(intDiff_RegC >0 && intDiff_RegD >0 && intDiff_RegC > 0) isRealSF = true;
  //double sf = (isRealSF)?ratioDiffRegDC: 1.0;
  //double err = (isRealSF)?errDiffRegDC:1.0; 
  double sf =  ratioDiffRegDC;
  double err = errDiffRegDC; 
  cout<<"-------------------------------------"<<endl;
  cout<<"intMC_RegC   = "<<intMC_RegC<<endl;
  cout<<"intData_RegC = "<<intData_RegC<<endl;
  cout<<"intMC_RegD   = "<<intMC_RegD<<endl;
  cout<<"intData_RegD = "<<intData_RegD<<endl;
  cout<<"intDiff_RegC = "<<intDiff_RegC<<endl;
  cout<<"intDiff_RegD = "<<intDiff_RegD<<endl;
  cout<<"sf 	      = "<<sf<<endl;
  vector<double>sfAndErr;
  sfAndErr.push_back(sf);
  sfAndErr.push_back(err);
  return sfAndErr;
}
