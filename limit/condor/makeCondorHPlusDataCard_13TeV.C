#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <algorithm> 

/////////////////////////// USERS INPUT ///////////////////////////
///INPUT FILES
//data
double totLumi = LUMI;
TString inFileDir="inFileDir_";
//MAKE DATA CARD of mjj
TString histname="inShapeHistoDir_/inShapeHisto_";
bool isMuChannel = isMuonChannel_;
bool isEleChannel = isEleChannel_;

TFile* fData    = TFile::Open(inFileDir+"/all_Data.root");
//OUTPUT FILE
TString outShapeFile ="HplusShapes";
TFile *fout = new TFile(outShapeFile+"_CHANNEL_inShapeHistoDir__inShapeHisto__13TeV.root", "RECREATE");

//bkg
TFile* fVV      = TFile::Open(inFileDir+"/all_VV.root");
TFile* fDY      = TFile::Open(inFileDir+"/all_DY.root");
TFile* fWJ      = TFile::Open(inFileDir+"/all_WJets.root");
TFile* fQCD     = TFile::Open(inFileDir+"/all_QCD.root");
TFile* fST      = TFile::Open(inFileDir+"/all_ST.root");
TFile* fTT      = TFile::Open(inFileDir+"/all_TTJetsP.root");
TFile* fTT_up      	= TFile::Open(inFileDir+"/all_TTJetsP_up.root");
TFile* fTT_down      	= TFile::Open(inFileDir+"/all_TTJetsP_down.root");
TFile* fTT_mtop1715     = TFile::Open(inFileDir+"/all_TTJetsP_mtop1715.root");
TFile* fTT_mtop1735     = TFile::Open(inFileDir+"/all_TTJetsP_mtop1735.root");
TFile* fTT_hdampUP      = TFile::Open(inFileDir+"/all_TTJetsP_hdampUP.root");
TFile* fTT_hdampDOWN    = TFile::Open(inFileDir+"/all_TTJetsP_hdampDOWN.root");
//signal
TFile *fWH80      = TFile::Open(inFileDir+"/all_Hplus80.root");
TFile *fWH90      = TFile::Open(inFileDir+"/all_Hplus90.root");
TFile *fWH100     = TFile::Open(inFileDir+"/all_Hplus100.root");
TFile *fWH120     = TFile::Open(inFileDir+"/all_Hplus120.root");
TFile *fWH140     = TFile::Open(inFileDir+"/all_Hplus140.root");
TFile *fWH150     = TFile::Open(inFileDir+"/all_Hplus150.root");
TFile *fWH155     = TFile::Open(inFileDir+"/all_Hplus155.root");
TFile *fWH160     = TFile::Open(inFileDir+"/all_Hplus160.root");
//data driven qcd
TFile* fQCD_dd = TFile::Open(inFileDir+"/all_QCD_dd.root"); 


//////////////////////////////////////////////////////////////////

//----------------------------------------//
//Variuos functions
//----------------------------------------//
//Read histos from input file. Return empty hist if the hist does not exist. Write to another file.
TH1F* readWriteHisto(TFile *inFile, TString histPath, TString inHistName, double sf, TFile *outFile, TString outHistName, bool isWrite = false){
  TH1F* hist;
  if(!(inFile->Get(histPath+inHistName))){
    hist = (TH1F*)(fTT->Get(histPath+inHistName))->Clone(outHistName);
    hist->Add(hist, -1);
  }else hist = (TH1F*)(inFile->Get(histPath+inHistName))->Clone(outHistName);
  hist->Scale(sf);
  hist->GetXaxis()->SetRangeUser(0, 200);
  if(isWrite){
    fout->cd();
    hist->Write(outHistName);
  }
  return hist;
}  

//---------------------------------------//
// for bin-by-bin uncertainity
//---------------------------------------//
TString cat_name_ = "inShapeHistoDir__inShapeHisto_"; 
int bin_min = 4; 
int bin_max = 40;
double ttbar_unc = 0.00;
double wh_unc    = 0.00;
double stop_unc  = 0.00;
double wjet_unc  = 0.00;

void writeBinnedHisto(TFile *outFile, TH1F* hist, TString hist_label_, TString cat_name_, int bin_min, int bin_max, double min_unc){
  string hist_label(hist_label_);
  string cat_name(cat_name_);
  for(int ibin = bin_min; ibin <= bin_max; ibin++){
    if(hist->GetBinContent(ibin) > 0){
      if(hist->GetBinError(ibin)/hist->GetBinContent(ibin) > min_unc){
        outFile->cd();
        TH1F* hist_up = (TH1F*)hist->Clone(Form("%s_CMS_stat_%s_%s_bin_%dUp", hist_label.c_str(), hist_label.c_str(), cat_name.c_str(), ibin));
        hist_up->SetBinContent(ibin, hist->GetBinContent(ibin)+hist->GetBinError(ibin));
        hist_up->Write(Form("%s_CMS_stat_%s_%s_bin_%dUp", hist_label.c_str(), hist_label.c_str(), cat_name.c_str(), ibin));
 
        TH1F* hist_down = (TH1F*)hist->Clone(Form("%s_CMS_stat_%s_%s_bin_%dDown", hist_label.c_str(), hist_label.c_str(), cat_name.c_str(), ibin));
        hist_down->SetBinContent(ibin, hist->GetBinContent(ibin)-hist->GetBinError(ibin));
        hist_down->Write(Form("%s_CMS_stat_%s_%s_bin_%dDown", hist_label.c_str(), hist_label.c_str(), cat_name.c_str(), ibin));
      }
    }
  }
}

void writeBinnedCard(ofstream &outFile, TH1F* hist, TString hist_label_, TString cat_name_, int bin_min, int bin_max, double min_unc, string dash_space){
  string hist_label(hist_label_);
  string cat_name(cat_name_);
  for(int ibin = bin_min; ibin <= bin_max; ibin++){
    if(hist->GetBinContent(ibin) > 0){
      if(hist->GetBinError(ibin)/hist->GetBinContent(ibin) > min_unc){
        outFile<< Form("CMS_stat_%s_%s_bin_%d", hist_label.c_str(), cat_name.c_str(), ibin)<<dash_space<<endl;
      }
    }
  }
}
//get normalised uncertainity
double getBTagUnc(TH1F *hCentral, TH1F* hUp, TH1F* hDown){
  return 1 + max(fabs(hUp->Integral() - hCentral->Integral()), fabs(hCentral->Integral() - hDown->Integral()))/hCentral->Integral();
}

//get statistical uncertainity
double getStatUnc(TH1F* hCentral, double sError = 0.0){
  double  norm = hCentral->IntegralAndError(1, hCentral->GetNbinsX(), sError);
  double statUnc = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00; 
  return statUnc;
}

//----------------------------------------//
//Global variables
//----------------------------------------//
//ttbar
double sf_ttbar_lumi = 1; 
TH1F* ttbar_topPtweight = readWriteHisto(fTT, "base/", "SF_topPtWeights", 1, fout, "topPtWeights", false);
double sf_ttbar = sf_ttbar_lumi/ttbar_topPtweight->GetMean();
TH1F* ttbar = readWriteHisto(fTT, "base/Iso/", histname, sf_ttbar, fout, "ttbar", true);
TH1F* ttbar_JESUp = readWriteHisto(fTT, "JESPlus/Iso/", histname, sf_ttbar, fout, "ttbar_JESUp", true);
TH1F* ttbar_JERUp = readWriteHisto(fTT, "JERPlus/Iso/", histname, sf_ttbar, fout, "ttbar_JERUp", true);
TH1F* ttbar_JESDown = readWriteHisto(fTT, "JESMinus/Iso/", histname, sf_ttbar, fout, "ttbar_JESDown", true);
TH1F* ttbar_JERDown = readWriteHisto(fTT, "JERMinus/Iso/", histname, sf_ttbar, fout, "ttbar_JERDown", true);
TH1F* ttbar_topPtUp = readWriteHisto(fTT, "TopPtPlus/Iso/", histname, sf_ttbar, fout, "ttbar_topPtUp", true);
TH1F* ttbar_topPtDown = readWriteHisto(fTT, "TopPtMinus/Iso/", histname, sf_ttbar, fout, "ttbar_topPtDown", true);
TH1F* ttbar_bTagUp = readWriteHisto(fTT, "bTagPlus/Iso/", histname, sf_ttbar, fout, "ttbar_bTagUp", true);
TH1F* ttbar_bTagDown = readWriteHisto(fTT, "bTagMinus/Iso/", histname, sf_ttbar, fout, "ttbar_bTagDown", true);
TH1F* ttbar_cTagUp = readWriteHisto(fTT, "cTagPlus/Iso/", histname, sf_ttbar, fout, "ttbar_cTagUp", true);
TH1F* ttbar_cTagDown = readWriteHisto(fTT, "cTagMinus/Iso/", histname, sf_ttbar, fout, "ttbar_cTagDown", true);

//ttbar scaleUp
double sf_ttbar_scaleUp_lumi = 1; 
TH1F* ttbar_scaleUp_topPtweight = readWriteHisto(fTT_up, "base/", "SF_topPtWeights", 1, fout, "topPtWeights", false);
double sf_ttbar_scaleUp_ptRewt = ttbar_scaleUp_topPtweight->GetMean();
TH1F* ttbar_scaleUp_ = readWriteHisto(fTT_up, "base/Iso/", histname, 1, fout, "ttbar_scaleUp", false);
double sf_ttbar_scaleUp_norm = ttbar->Integral()/ttbar_scaleUp_->Integral();
double sf_ttbar_scaleUp = (sf_ttbar_scaleUp_lumi/sf_ttbar_scaleUp_ptRewt)*sf_ttbar_scaleUp_norm;
TH1F* ttbar_scaleUp = readWriteHisto(fTT_up, "base/Iso/", histname, sf_ttbar_scaleUp, fout, "ttbar_scaleUp", true);

//ttbar scaleDown
double sf_ttbar_scaleDown_lumi = 1; 
TH1F* ttbar_scaleDown_topPtweight = readWriteHisto(fTT_down, "base/", "SF_topPtWeights", 1, fout, "topPtWeights", false);
double sf_ttbar_scaleDown_ptRewt = ttbar_scaleDown_topPtweight->GetMean();
TH1F* ttbar_scaleDown_ = readWriteHisto(fTT_down, "base/Iso/", histname, 1, fout, "ttbar_scaleDown", false);
double sf_ttbar_scaleDown_norm = ttbar->Integral()/ttbar_scaleDown_->Integral();
double sf_ttbar_scaleDown = (sf_ttbar_scaleDown_lumi/sf_ttbar_scaleDown_ptRewt)*sf_ttbar_scaleDown_norm;
TH1F* ttbar_scaleDown = readWriteHisto(fTT_down, "base/Iso/", histname, sf_ttbar_scaleDown, fout, "ttbar_scaleDown", true);

//ttbar mtop1715
double sf_ttbar_mtop1715_lumi = 1; 
TH1F* ttbar_mtop1715_topPtweight = readWriteHisto(fTT_mtop1715, "base/", "SF_topPtWeights", 1, fout, "topPtWeights", false);
double sf_ttbar_mtop1715_ptRewt = ttbar_mtop1715_topPtweight->GetMean();
TH1F* ttbar_mtop1715_ = readWriteHisto(fTT_mtop1715, "base/Iso/", histname, 1, fout, "ttbar_massUp", false);
double sf_ttbar_mtop1715_norm = ttbar->Integral()/ttbar_mtop1715_->Integral();
double sf_ttbar_mtop1715 = (sf_ttbar_mtop1715_lumi/sf_ttbar_mtop1715_ptRewt)*sf_ttbar_mtop1715_norm;
TH1F* ttbar_mtop1715 = readWriteHisto(fTT_mtop1715, "base/Iso/", histname, sf_ttbar_mtop1715, fout, "ttbar_massUp", true);

//ttbar mtop1735
double sf_ttbar_mtop1735_lumi = 1; 
TH1F* ttbar_mtop1735_topPtweight = readWriteHisto(fTT_mtop1735, "base/", "SF_topPtWeights", 1, fout, "topPtWeights", false);
double sf_ttbar_mtop1735_ptRewt = ttbar_mtop1735_topPtweight->GetMean();
TH1F* ttbar_mtop1735_ = readWriteHisto(fTT_mtop1735, "base/Iso/", histname, 1, fout, "ttbar_massDown", false);
double sf_ttbar_mtop1735_norm = ttbar->Integral()/ttbar_mtop1735_->Integral();
double sf_ttbar_mtop1735 = (sf_ttbar_mtop1735_lumi/sf_ttbar_mtop1735_ptRewt)*sf_ttbar_mtop1735_norm;
TH1F* ttbar_mtop1735 = readWriteHisto(fTT_mtop1735, "base/Iso/", histname, sf_ttbar_mtop1735, fout, "ttbar_massDown", true);

//ttbar matchingUp
double sf_ttbar_matchingUp_lumi = 1; 
TH1F* ttbar_matchingUp_topPtweight = readWriteHisto(fTT_hdampUP, "base/", "SF_topPtWeights", 1, fout, "topPtWeights", false);
double sf_ttbar_matchingUp_ptRewt = ttbar_matchingUp_topPtweight->GetMean();
TH1F* ttbar_matchingUp_ = readWriteHisto(fTT_hdampUP, "base/Iso/", histname, 1, fout, "ttbar_matchingUp", false);
double sf_ttbar_matchingUp_norm = ttbar->Integral()/ttbar_matchingUp_->Integral();
double sf_ttbar_matchingUp = (sf_ttbar_matchingUp_lumi/sf_ttbar_matchingUp_ptRewt)*sf_ttbar_matchingUp_norm;
TH1F* ttbar_matchingUp = readWriteHisto(fTT_hdampUP, "base/Iso/", histname, sf_ttbar_matchingUp, fout, "ttbar_matchingUp", true);

//ttbar matchingDown
double sf_ttbar_matchingDown_lumi = 1; 
TH1F* ttbar_matchingDown_topPtweight = readWriteHisto(fTT_hdampDOWN, "base/", "SF_topPtWeights", 1, fout, "topPtWeights", false);
double sf_ttbar_matchingDown_ptRewt = ttbar_matchingDown_topPtweight->GetMean();
TH1F* ttbar_matchingDown_ = readWriteHisto(fTT_hdampDOWN, "base/Iso/", histname, 1, fout, "ttbar_matchingDown", false);
double sf_ttbar_matchingDown_norm = ttbar->Integral()/ttbar_matchingDown_->Integral();
double sf_ttbar_matchingDown = (sf_ttbar_matchingDown_lumi/sf_ttbar_matchingDown_ptRewt)*sf_ttbar_matchingDown_norm;
TH1F* ttbar_matchingDown = readWriteHisto(fTT_hdampDOWN, "base/Iso/", histname, sf_ttbar_matchingDown, fout, "ttbar_matchingDown", true);

//ttll
double sf_ttll = 0;
TH1F* ttll = readWriteHisto(fTT, "base/Iso/", histname, sf_ttll, fout, "ttll", true);
TH1F* ttll_JESUp = readWriteHisto(fTT, "JESPlus/Iso/", histname, sf_ttll, fout, "ttll_JESUp", true);
TH1F* ttll_JESDown = readWriteHisto(fTT, "JESMinus/Iso/", histname, sf_ttll, fout, "ttll_JESDown", true);

//w+jets
double sf_wjet = 1;
TH1F* wjet = readWriteHisto(fWJ, "base/Iso/", histname, sf_wjet, fout, "wjet", true);
TH1F* wjet_JESUp = readWriteHisto(fWJ, "JESPlus/Iso/", histname, sf_wjet, fout, "wjet_JESUp", true);
TH1F* wjet_JERUp = readWriteHisto(fWJ, "JERPlus/Iso/", histname, sf_wjet, fout, "wjet_JERUp", true);
TH1F* wjet_JESDown = readWriteHisto(fWJ, "JESMinus/Iso/", histname, sf_wjet, fout, "wjet_JESDown", true);
TH1F* wjet_JERDown = readWriteHisto(fWJ, "JERMinus/Iso/", histname, sf_wjet, fout, "wjet_JERDown", true);
TH1F* wjet_bTagUp = readWriteHisto(fWJ, "bTagPlus/Iso/", histname, sf_wjet, fout, "wjet_bTagUp", true);
TH1F* wjet_bTagDown = readWriteHisto(fWJ, "bTagMinus/Iso/", histname, sf_wjet, fout, "wjet_bTagDown", true); 
TH1F* wjet_cTagUp = readWriteHisto(fWJ, "cTagPlus/Iso/", histname, sf_wjet, fout, "wjet_cTagUp", true);
TH1F* wjet_cTagDown = readWriteHisto(fWJ, "cTagMinus/Iso/", histname, sf_wjet, fout, "wjet_cTagDown", true); 

//Z+Jets
double sf_zjet = 1;
TH1F* zjet = readWriteHisto(fDY, "base/Iso/", histname, sf_zjet, fout, "zjet", true);
TH1F* zjet_JESUp = readWriteHisto(fDY, "JESPlus/Iso/", histname, sf_zjet, fout, "zjet_JESUp", true);
TH1F* zjet_JERUp = readWriteHisto(fDY, "JERPlus/Iso/", histname, sf_zjet, fout, "zjet_JERUp", true);
TH1F* zjet_JESDown = readWriteHisto(fDY, "JESMinus/Iso/", histname, sf_zjet, fout, "zjet_JESDown", true);
TH1F* zjet_JERDown = readWriteHisto(fDY, "JERMinus/Iso/", histname, sf_zjet, fout, "zjet_JERDown", true);
TH1F* zjet_bTagUp = readWriteHisto(fDY, "bTagPlus/Iso/", histname, sf_zjet, fout, "zjet_bTagUp", true);
TH1F* zjet_bTagDown = readWriteHisto(fDY, "bTagMinus/Iso/", histname, sf_zjet, fout, "zjet_bTagDown", true);
TH1F* zjet_cTagUp = readWriteHisto(fDY, "cTagPlus/Iso/", histname, sf_zjet, fout, "zjet_cTagUp", true);
TH1F* zjet_cTagDown = readWriteHisto(fDY, "cTagMinus/Iso/", histname, sf_zjet, fout, "zjet_cTagDown", true);

//SingleTop
double sf_stop = 1;
TH1F* stop = readWriteHisto(fST, "base/Iso/", histname, sf_stop, fout, "stop", true);
TH1F* stop_JESUp = readWriteHisto(fST, "JESPlus/Iso/", histname, sf_stop, fout, "stop_JESUp", true);
TH1F* stop_JERUp = readWriteHisto(fST, "JERPlus/Iso/", histname, sf_stop, fout, "stop_JERUp", true);
TH1F* stop_JESDown = readWriteHisto(fST, "JESMinus/Iso/", histname, sf_stop, fout, "stop_JESDown", true);
TH1F* stop_JERDown = readWriteHisto(fST, "JERMinus/Iso/", histname, sf_stop, fout, "stop_JERDown", true);
TH1F* stop_bTagUp = readWriteHisto(fST, "bTagPlus/Iso/", histname, sf_stop, fout, "stop_bTagUp", true);
TH1F* stop_bTagDown = readWriteHisto(fST, "bTagMinus/Iso/", histname, sf_stop, fout, "stop_bTagDown", true);
TH1F* stop_cTagUp = readWriteHisto(fST, "cTagPlus/Iso/", histname, sf_stop, fout, "stop_cTagUp", true);
TH1F* stop_cTagDown = readWriteHisto(fST, "cTagMinus/Iso/", histname, sf_stop, fout, "stop_cTagDown", true);

//Dibosons
double sf_diboson = 1;
TH1F* diboson = readWriteHisto(fVV, "base/Iso/", histname, sf_diboson, fout, "diboson", true);
TH1F* diboson_JESUp = readWriteHisto(fVV, "JESPlus/Iso/", histname, sf_diboson, fout, "diboson_JESUp", true);
TH1F* diboson_JERUp = readWriteHisto(fVV, "JERPlus/Iso/", histname, sf_diboson, fout, "diboson_JERUp", true);
TH1F* diboson_JESDown = readWriteHisto(fVV, "JESMinus/Iso/", histname, sf_diboson, fout, "diboson_JESDown", true);
TH1F* diboson_JERDown = readWriteHisto(fVV, "JERMinus/Iso/", histname, sf_diboson, fout, "diboson_JERDown", true);
TH1F* diboson_bTagUp = readWriteHisto(fVV, "bTagPlus/Iso/", histname, sf_diboson, fout, "diboson_bTagUp", true);
TH1F* diboson_bTagDown = readWriteHisto(fVV, "bTagMinus/Iso/", histname, sf_diboson, fout, "diboson_bTagDown", true);
TH1F* diboson_cTagUp = readWriteHisto(fVV, "cTagPlus/Iso/", histname, sf_diboson, fout, "diboson_cTagUp", true);
TH1F* diboson_cTagDown = readWriteHisto(fVV, "cTagMinus/Iso/", histname, sf_diboson, fout, "diboson_cTagDown", true);

//QCD MC
double sf_qcd = 1;
TH1F* qcd = readWriteHisto(fQCD, "base/Iso/", histname, sf_qcd, fout, "qcd", true);
TH1F* qcd_JESUp = readWriteHisto(fQCD, "JESPlus/Iso/", histname, sf_qcd, fout, "qcd_JESUp", true);
TH1F* qcd_JERUp = readWriteHisto(fQCD, "JERPlus/Iso/", histname, sf_qcd, fout, "qcd_JERUp", true);
TH1F* qcd_JESDown = readWriteHisto(fQCD, "JESMinus/Iso/", histname, sf_qcd, fout, "qcd_JESDown", true);
TH1F* qcd_JERDown = readWriteHisto(fQCD, "JERMinus/Iso/", histname, sf_qcd, fout, "qcd_JERDown", true);
TH1F* qcd_bTagUp = readWriteHisto(fQCD, "bTagPlus/Iso/", histname, sf_qcd, fout, "qcd_bTagUp", true);
TH1F* qcd_bTagDown = readWriteHisto(fQCD, "bTagMinus/Iso/", histname, sf_qcd, fout, "qcd_bTagDown", true);
TH1F* qcd_cTagUp = readWriteHisto(fQCD, "cTagPlus/Iso/", histname, sf_qcd, fout, "qcd_cTagUp", true);
TH1F* qcd_cTagDown = readWriteHisto(fQCD, "cTagMinus/Iso/", histname, sf_qcd, fout, "qcd_cTagDown", true);
//QCD data driven
TH1F* qcd_dd = readWriteHisto(fQCD_dd, "base/Iso/", histname, sf_qcd, fout, "qcd", true);

//Data
double sf_data = 1; //should be 1, always
TH1F* data_obs = readWriteHisto(fData, "base/Iso/", histname, sf_data, fout, "data_obs", true);


//----------------------------------------//
//function to make data card for each mass
//----------------------------------------//
void makeOneDataCard(TFile *fWH, int mass=80, TString label="WH80", TString label2="HH80", TString histname = "mjj_kfit"){
  cout<<" mass point: "<<mass<<endl;
  //wh
  double sf_wh_lumi = 1; 
  TH1F* wh_topPtweight = readWriteHisto(fWH, "base/", "SF_topPtWeights", 1, fout, "topPtWeight"+label, false);
  double sf_wh = sf_wh_lumi/wh_topPtweight->GetMean();
  TH1F* wh = readWriteHisto(fWH, "base/Iso/", histname, sf_wh, fout, label, true);
  TH1F* wh_JESUp = readWriteHisto(fWH, "JESPlus/Iso/", histname, sf_wh, fout, label+"_JESUp", true);
  TH1F* wh_JERUp = readWriteHisto(fWH, "JERPlus/Iso/", histname, sf_wh, fout, label+"_JERUp", true);
  TH1F* wh_JESDown = readWriteHisto(fWH, "JESMinus/Iso/", histname, sf_wh, fout, label+"_JESDown", true);
  TH1F* wh_JERDown = readWriteHisto(fWH, "JERMinus/Iso/", histname, sf_wh, fout, label+"_JERDown", true);
  TH1F* wh_topPtUp = readWriteHisto(fWH,  "TopPtPlus/Iso/", histname, sf_wh, fout, label+"_topPtUp", true);
  TH1F* wh_topPtDown = readWriteHisto(fWH, "TopPtMinus/Iso/", histname, sf_wh, fout, label+"_topPtDown", true);
  TH1F* wh_bTagUp = readWriteHisto(fWH, "bTagPlus/Iso/", histname, sf_wh, fout, label+"_bTagUp", true);
  TH1F* wh_bTagDown = readWriteHisto(fWH, "bTagMinus/Iso/", histname, sf_wh, fout, label+"_bTagDown", true);
  TH1F* wh_cTagUp = readWriteHisto(fWH, "cTagPlus/Iso/", histname, sf_wh, fout, label+"_cTagUp", true);
  TH1F* wh_cTagDown = readWriteHisto(fWH, "cTagMinus/Iso/", histname, sf_wh, fout, label+"_cTagDown", true);
  writeBinnedHisto(fout, wh, label, cat_name_, 0, 40, 0.02);
  
  //hh
  double sf_hh = 0; 
  TH1F* hh = readWriteHisto(fWH, "base/Iso/", histname, sf_hh, fout, label2, true);
  TH1F* hh_JESUp = readWriteHisto(fWH, "JESPlus/Iso/", histname, sf_hh, fout, label2+"_JESUp", true);
  TH1F* hh_JERUp = readWriteHisto(fWH, "JERPlus/Iso/", histname, sf_hh, fout, label2+"_JERUp", true);
  TH1F* hh_JESDown = readWriteHisto(fWH, "JESMinus/Iso/", histname, sf_hh, fout, label2+"_JESDown", true);
  TH1F* hh_JERDown = readWriteHisto(fWH, "JERMinus/Iso/", histname, sf_hh, fout, label2+"_JERDown", true);
  TH1F* hh_topPtUp = readWriteHisto(fWH,  "TopPtPlus/Iso/", histname, sf_hh, fout, label2+"_topPtUp", true);
  TH1F* hh_topPtDown = readWriteHisto(fWH, "TopPtMinus/Iso/", histname, sf_hh, fout, label2+"_topPtDown", true);
  TH1F* hh_bTagUp = readWriteHisto(fWH, "bTagPlus/Iso/", histname, sf_hh, fout, label2+"_bTagUp", true);
  TH1F* hh_bTagDown = readWriteHisto(fWH, "bTagMinus/Iso/", histname, sf_hh, fout, label2+"_bTagDown", true);
  TH1F* hh_cTagUp = readWriteHisto(fWH, "cTagPlus/Iso/", histname, sf_hh, fout, label2+"_cTagUp", true);
  TH1F* hh_cTagDown = readWriteHisto(fWH, "cTagMinus/Iso/", histname, sf_hh, fout, label2+"_cTagDown", true);

  //open input template data card of 8 TeV
  ifstream in;
  char* c = new char[1000];
  in.open("template_datacard_csbar.log");
  //create output data card for 13 TeV
  string outDataCard = "datacard_csbar_13TeV_mH.txt";
  string histname_str("inShapeHistoDir__inShapeHisto_");
  if(isMuChannel) outDataCard = "datacard_csbar_mu_"+histname_str+"_13TeV_mH%d.txt"; 
  else outDataCard = "datacard_csbar_ele_"+histname_str+"_13TeV_mH%d.txt";
  ofstream out(Form(outDataCard.c_str(), mass));
  out.precision(8);

  time_t secs=time(0);
  tm *t=localtime(&secs);
  while (in.good()){
    in.getline(c,1000,'\n');
    if (in.good()){
      string line(c);
      if(line.find("Date")!=string::npos){
        string day = string(Form("%d",t->tm_mday));
        string month = string(Form("%d",t->tm_mon+1));
        string year = string(Form("%d",t->tm_year+1900));
        line.replace( line.find("XXX") , 3 , day+"/"+month+"/"+year);
        out << line << endl;
      }
      else if(line.find("Description")!=string::npos){
        line.replace( line.find("YYY") , 3 , string(Form("%d", mass)) );
        line.replace( line.find("ZZZ") , 3 , string(Form("%f", totLumi)) ); 
        out << line << endl;
      }
      else if(line.find("shapes")!=string::npos){
        line.replace( line.find("XXX") , 3 , string(outShapeFile+"_CHANNEL_inShapeHistoDir__inShapeHisto__13TeV"));
        out << line << endl;
      }
      else if(line.find("Observation")!=string::npos){
        line.replace( line.find("XXX") , 3 , string(Form("%.0f", data_obs->Integral())));
        out << line << endl;
      }
      else if(line.find("process")!=string::npos && line.find("WH")!=string::npos){
        line.replace( line.find("XXX") , 3 , string(Form("%d", mass)) );
        line.replace( line.find("YYY") , 3 , string(Form("%d", mass)) );
        out << line << endl;
      }
      else if(line.find("rate")!=string::npos){
        string rate = "rate               ";  
        string space = "     ";
        out << rate ;
        out << space << 0*wh->Integral()
            << space << wh->Integral()
            << space << ttbar->Integral()
            << space << ttll->Integral()
            << space << wjet->Integral()
            << space << zjet->Integral()
            << space << stop->Integral()
            << space << diboson->Integral()
            << space << qcd_dd->Integral()
            << endl;
      }
      else if(line.find("CMS_stat_wh")!=string::npos){ 
        line.replace( line.find("XXXX") , 4 , string(Form("%.3f", getStatUnc(wh,  0))));  
        out << line << endl;
      } 
      else if(line.find("CMS_stat_tt")!=string::npos){  
	line.replace( line.find("XXXX") , 4 , string(Form("%.3f", getStatUnc(ttbar,  0))));   
        out << line << endl;
      }  
      else if(line.find("CMS_stat_wjet")!=string::npos){  
        line.replace( line.find("XXXX") , 4 , string(Form("%.3f", getStatUnc(wjet,  0))));   
        out << line << endl;
      }  
      else if(line.find("CMS_stat_zjet")!=string::npos){ 
        line.replace( line.find("XXXX") , 4 , string(Form("%.3f", getStatUnc(zjet,  0))));  
        out << line << endl; 
      }
      else if(line.find("CMS_stat_stop")!=string::npos){ 
        line.replace( line.find("XXXX") , 4 , string(Form("%.3f", getStatUnc(stop,  0))));  
        out << line << endl; 
      } 
      else if(line.find("CMS_stat_vv")!=string::npos){  
        line.replace( line.find("XXXX") , 4 , string(Form("%.3f", getStatUnc(diboson,  0))));   
        out << line << endl;  
      }
      else if(line.find("CMS_stat_qcd")!=string::npos){  
        line.replace( line.find("XXXX") , 4 , string(Form("%.3f", getStatUnc(qcd_dd,  0))));   
        out << line << endl;  
      }
      else{ //default without changes
        out << line << endl;
      }
    }
  } 
  writeBinnedCard(out, wh, label, cat_name_, bin_min, bin_max, wh_unc,          "   shape   -       1.00    -         -         -        -        -        -      -");
  writeBinnedCard(out, ttbar, "ttbar", cat_name_, bin_min, bin_max, ttbar_unc,  "   shape   -       -       1.00      -         -        -        -        -      -");
  writeBinnedCard(out, stop,  "stop", cat_name_, bin_min, bin_max, stop_unc,    "   shape   -       -       -         -         -        -        1.00     -      -" );
  writeBinnedCard(out, wjet,  "wjet", cat_name_, bin_min, bin_max, wjet_unc,    "   shape   -       -       -         -         1.00     -        -        -      -");
  out.close();
  in.close();
}

//----------------------------------------//
//make data card for all masses
//----------------------------------------//
void makeCondorHPlusDataCard_13TeV(){
  writeBinnedHisto(fout, ttbar,"ttbar", cat_name_, bin_min, bin_max, ttbar_unc);
  writeBinnedHisto(fout, stop, "stop", cat_name_, bin_min, bin_max, stop_unc);
  writeBinnedHisto(fout, wjet, "wjet", cat_name_, bin_min, bin_max, wjet_unc);
  makeOneDataCard(fWH80,  80,  "WH80",  "HH80",  histname);
  makeOneDataCard(fWH90,  90,  "WH90",  "HH90",  histname);
  makeOneDataCard(fWH100, 100, "WH100", "HH100", histname);
  makeOneDataCard(fWH120, 120, "WH120", "HH120", histname);
  makeOneDataCard(fWH140, 140, "WH140", "HH140", histname);
  makeOneDataCard(fWH150, 150, "WH150", "HH150", histname);
  makeOneDataCard(fWH155, 155, "WH155", "HH155", histname);
  makeOneDataCard(fWH160, 160, "WH160", "HH160", histname);
  fout->Close();
}

