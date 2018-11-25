#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <algorithm> 
#include "MyHPlusDataCardMaker.h"

//----------------------------------------//
//make data card for each mass
//----------------------------------------//
MyHPlusDataCardMaker DC;
double totLumi = 35.9;
bool isNeffThreshold_ttbar = false;
bool isNeffThreshold_wjet = true;
bool isNeffThreshold_zjet = true;
bool isNeffThreshold_stop = true;
bool isNeffThreshold_vv = true;
bool isNeffThreshold_qcd_dd = true;
double minNeffThres = 2;

void MyHPlusDataCardMaker(TString inFileDir="stack_20180418_Mu_Sys_PreAppComent", 
        TString histSubDir_="KinFit", 
        TString histName="mjj_kfit", 
        TString channelName="mu", 
        int mass=80, 
        TString label="WH80", 
        TString hPlusFileName="all_Hplus80.root")
  {
  TString histSubDir = "Iso/"+histSubDir_+"/";
  bool isMuChannel = false; 
  if(channelName=="mu") isMuChannel = true;
  ///INPUT FILES
  TFile* fData    = TFile::Open(inFileDir+"/all_Data.root");
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
  TFile *fWH  = TFile::Open(inFileDir+"/"+hPlusFileName);
  //data driven qcd
  TFile* fQCD_dd = TFile::Open(inFileDir+"/all_QCD_dd.root"); 
  
  //OUTPUT FILE
  TFile *fout = new TFile(TString("HplusShapes_")+channelName+TString("_")+histSubDir_+TString("_")+histName+TString("_13TeV_")+label+TString(".root"), "RECREATE");

  //ttbar
  double sf_ttbar = 1; 
  TH1F* ttbar = DC.readWriteHisto(fTT, "base/"+histSubDir, histName, sf_ttbar, fout, fTT,  "ttbar", true, minNeffThres, isNeffThreshold_ttbar);
  TH1F* ttbar_JESUp = DC.readWriteHisto(fTT, "JESPlus/"+histSubDir, histName, sf_ttbar, fout, fTT, "ttbar_JESUp", true);
  TH1F* ttbar_JESDown = DC.readWriteHisto(fTT, "JESMinus/"+histSubDir, histName, sf_ttbar, fout, fTT, "ttbar_JESDown", true);
  TH1F* ttbar_PileupUp = DC.readWriteHisto(fTT, "PileupPlus/"+histSubDir, histName, sf_ttbar, fout, fTT, "ttbar_PileupUp", true);
  TH1F* ttbar_PileupDown = DC.readWriteHisto(fTT, "PileupMinus/"+histSubDir, histName, sf_ttbar, fout, fTT, "ttbar_PileupDown", true);
  TH1F* ttbar_JERUp = DC.readWriteHisto(fTT, "JERPlus/"+histSubDir, histName, sf_ttbar, fout, fTT, "ttbar_JERUp", true);
  TH1F* ttbar_JERDown = DC.readWriteHisto(fTT, "JERMinus/"+histSubDir, histName, sf_ttbar, fout, fTT, "ttbar_JERDown", true);
  TH1F* ttbar_topPtUp = DC.readWriteHisto(fTT, "TopPtPlus/"+histSubDir, histName, sf_ttbar, fout, fTT, "ttbar_topPtUp", true);
  TH1F* ttbar_topPtDown = DC.readWriteHisto(fTT, "TopPtMinus/"+histSubDir, histName, sf_ttbar, fout, fTT, "ttbar_topPtDown", true);
  TH1F* ttbar_bTagUp = DC.readWriteHisto(fTT, "bTagPlus/"+histSubDir, histName, sf_ttbar, fout, fTT, "ttbar_bTagUp", true);
  TH1F* ttbar_bTagDown = DC.readWriteHisto(fTT, "bTagMinus/"+histSubDir, histName, sf_ttbar, fout, fTT, "ttbar_bTagDown", true);
  TH1F* ttbar_cTagUp = DC.readWriteHisto(fTT, "cTagPlus/"+histSubDir, histName, sf_ttbar, fout, fTT, "ttbar_cTagUp", true);
  TH1F* ttbar_cTagDown = DC.readWriteHisto(fTT, "cTagMinus/"+histSubDir, histName, sf_ttbar, fout, fTT, "ttbar_cTagDown", true);


  //ttbar scaleUp
  double sf_ttbar_scaleUp = 1; 
  TH1F* ttbar_scaleUp_ = DC.readWriteHisto(fTT_up, "base/"+histSubDir, histName, 1, fout, fTT, "ttbar_scaleUp", false);
  double sf_ttbar_scaleUp_norm = ttbar->Integral()/ttbar_scaleUp_->Integral();
  sf_ttbar_scaleUp = (sf_ttbar_scaleUp)*sf_ttbar_scaleUp_norm;
  TH1F* ttbar_scaleUp = DC.readWriteHisto(fTT_up, "base/"+histSubDir, histName, sf_ttbar_scaleUp, fout, fTT, "ttbar_scaleUp", true);

  //ttbar scaleDown
  double sf_ttbar_scaleDown = 1; 
  TH1F* ttbar_scaleDown_ = DC.readWriteHisto(fTT_down, "base/"+histSubDir, histName, 1, fout, fTT, "ttbar_scaleDown", false);
  double sf_ttbar_scaleDown_norm = ttbar->Integral()/ttbar_scaleDown_->Integral();
  sf_ttbar_scaleDown = (sf_ttbar_scaleDown)*sf_ttbar_scaleDown_norm;
  TH1F* ttbar_scaleDown = DC.readWriteHisto(fTT_down, "base/"+histSubDir, histName, sf_ttbar_scaleDown, fout, fTT, "ttbar_scaleDown", true);

  //ttbar mtop1715
  double sf_ttbar_mtop1715 = 1; 
  TH1F* ttbar_mtop1715_ = DC.readWriteHisto(fTT_mtop1715, "base/"+histSubDir, histName, 1, fout, fTT, "ttbar_massUp", false);
  double sf_ttbar_mtop1715_norm = ttbar->Integral()/ttbar_mtop1715_->Integral();
  sf_ttbar_mtop1715 = (sf_ttbar_mtop1715)*sf_ttbar_mtop1715_norm;
  TH1F* ttbar_mtop1715 = DC.readWriteHisto(fTT_mtop1715, "base/"+histSubDir, histName, sf_ttbar_mtop1715, fout, fTT, "ttbar_massUp", true);

  //ttbar mtop1735
  double sf_ttbar_mtop1735 = 1; 
  TH1F* ttbar_mtop1735_ = DC.readWriteHisto(fTT_mtop1735, "base/"+histSubDir, histName, 1, fout, fTT, "ttbar_massDown", false);
  double sf_ttbar_mtop1735_norm = ttbar->Integral()/ttbar_mtop1735_->Integral();
  sf_ttbar_mtop1735 = (sf_ttbar_mtop1735)*sf_ttbar_mtop1735_norm;
  TH1F* ttbar_mtop1735 = DC.readWriteHisto(fTT_mtop1735, "base/"+histSubDir, histName, sf_ttbar_mtop1735, fout, fTT, "ttbar_massDown", true);

  //ttbar matchingUp
  double sf_ttbar_matchingUp = 1; 
  TH1F* ttbar_matchingUp_ = DC.readWriteHisto(fTT_hdampUP, "base/"+histSubDir, histName, 1, fout, fTT, "ttbar_matchingUp", false);
  double sf_ttbar_matchingUp_norm = ttbar->Integral()/ttbar_matchingUp_->Integral();
  sf_ttbar_matchingUp = (sf_ttbar_matchingUp)*sf_ttbar_matchingUp_norm;
  TH1F* ttbar_matchingUp = DC.readWriteHisto(fTT_hdampUP, "base/"+histSubDir, histName, sf_ttbar_matchingUp, fout, fTT, "ttbar_matchingUp", true);

  //ttbar matchingDown
  double sf_ttbar_matchingDown = 1; 
  TH1F* ttbar_matchingDown_ = DC.readWriteHisto(fTT_hdampDOWN, "base/"+histSubDir, histName, 1, fout, fTT, "ttbar_matchingDown", false);
  double sf_ttbar_matchingDown_norm = ttbar->Integral()/ttbar_matchingDown_->Integral();
  sf_ttbar_matchingDown = (sf_ttbar_matchingDown)*sf_ttbar_matchingDown_norm;
  TH1F* ttbar_matchingDown = DC.readWriteHisto(fTT_hdampDOWN, "base/"+histSubDir, histName, sf_ttbar_matchingDown, fout, fTT, "ttbar_matchingDown", true);

  //ttll
  double sf_ttll = 0;
  TH1F* ttll = DC.readWriteHisto(fTT, "base/"+histSubDir, histName, sf_ttll, fout, fTT, "ttll", true);
  TH1F* ttll_JESUp = DC.readWriteHisto(fTT, "JESPlus/"+histSubDir, histName, sf_ttll, fout, fTT, "ttll_JESUp", true);
  TH1F* ttll_JESDown = DC.readWriteHisto(fTT, "JESMinus/"+histSubDir, histName, sf_ttll, fout, fTT, "ttll_JESDown", true);

  //w+jets
  double sf_wjet = 1;
  TH1F* wjet = DC.readWriteHisto(fWJ, "base/"+histSubDir, histName, sf_wjet, fout, fTT, "wjet", true, minNeffThres, isNeffThreshold_wjet);
  TH1F* wjet_JESUp = DC.readWriteHisto(fWJ, "JESPlus/"+histSubDir, histName, sf_wjet, fout, fTT, "wjet_JESUp", true, minNeffThres, isNeffThreshold_wjet);
  TH1F* wjet_JESDown = DC.readWriteHisto(fWJ, "JESMinus/"+histSubDir, histName, sf_wjet, fout, fTT, "wjet_JESDown", true, minNeffThres, isNeffThreshold_wjet);
  TH1F* wjet_PileupUp = DC.readWriteHisto(fWJ, "PileupPlus/"+histSubDir, histName, sf_wjet, fout, fTT, "wjet_PileupUp", true, minNeffThres, isNeffThreshold_wjet);
  TH1F* wjet_PileupDown = DC.readWriteHisto(fWJ, "PileupMinus/"+histSubDir, histName, sf_wjet, fout, fTT, "wjet_PileupDown", true, minNeffThres, isNeffThreshold_wjet);
  TH1F* wjet_JERUp = DC.readWriteHisto(fWJ, "JERPlus/"+histSubDir, histName, sf_wjet, fout, fTT, "wjet_JERUp", true, minNeffThres, isNeffThreshold_wjet);
  TH1F* wjet_JERDown = DC.readWriteHisto(fWJ, "JERMinus/"+histSubDir, histName, sf_wjet, fout, fTT, "wjet_JERDown", true, minNeffThres, isNeffThreshold_wjet);
  TH1F* wjet_bTagUp = DC.readWriteHisto(fWJ, "bTagPlus/"+histSubDir, histName, sf_wjet, fout, fTT, "wjet_bTagUp", true, minNeffThres, isNeffThreshold_wjet);
  TH1F* wjet_bTagDown = DC.readWriteHisto(fWJ, "bTagMinus/"+histSubDir, histName, sf_wjet, fout, fTT, "wjet_bTagDown", true, minNeffThres, isNeffThreshold_wjet); 
  TH1F* wjet_cTagUp = DC.readWriteHisto(fWJ, "cTagPlus/"+histSubDir, histName, sf_wjet, fout, fTT, "wjet_cTagUp", true, minNeffThres, isNeffThreshold_wjet);
  TH1F* wjet_cTagDown = DC.readWriteHisto(fWJ, "cTagMinus/"+histSubDir, histName, sf_wjet, fout, fTT, "wjet_cTagDown", true, minNeffThres, isNeffThreshold_wjet); 

  //Z+Jets
  double sf_zjet = 1;
  TH1F* zjet = DC.readWriteHisto(fDY, "base/"+histSubDir, histName, sf_zjet, fout, fTT, "zjet", true, minNeffThres, isNeffThreshold_zjet);
  TH1F* zjet_JESUp = DC.readWriteHisto(fDY, "JESPlus/"+histSubDir, histName, sf_zjet, fout, fTT, "zjet_JESUp", true, minNeffThres, isNeffThreshold_zjet);
  TH1F* zjet_JESDown = DC.readWriteHisto(fDY, "JESMinus/"+histSubDir, histName, sf_zjet, fout, fTT, "zjet_JESDown", true, minNeffThres, isNeffThreshold_zjet);
  TH1F* zjet_PileupUp = DC.readWriteHisto(fDY, "PileupPlus/"+histSubDir, histName, sf_zjet, fout, fTT, "zjet_PileupUp", true, minNeffThres, isNeffThreshold_zjet);
  TH1F* zjet_PileupDown = DC.readWriteHisto(fDY, "PileupMinus/"+histSubDir, histName, sf_zjet, fout, fTT, "zjet_PileupDown", true, minNeffThres, isNeffThreshold_zjet);
  TH1F* zjet_JERUp = DC.readWriteHisto(fDY, "JERPlus/"+histSubDir, histName, sf_zjet, fout, fTT, "zjet_JERUp", true, minNeffThres, isNeffThreshold_zjet);
  TH1F* zjet_JERDown = DC.readWriteHisto(fDY, "JERMinus/"+histSubDir, histName, sf_zjet, fout, fTT, "zjet_JERDown", true, minNeffThres, isNeffThreshold_zjet);
  TH1F* zjet_bTagUp = DC.readWriteHisto(fDY, "bTagPlus/"+histSubDir, histName, sf_zjet, fout, fTT, "zjet_bTagUp", true, minNeffThres, isNeffThreshold_zjet);
  TH1F* zjet_bTagDown = DC.readWriteHisto(fDY, "bTagMinus/"+histSubDir, histName, sf_zjet, fout, fTT, "zjet_bTagDown", true, minNeffThres, isNeffThreshold_zjet);
  TH1F* zjet_cTagUp = DC.readWriteHisto(fDY, "cTagPlus/"+histSubDir, histName, sf_zjet, fout, fTT, "zjet_cTagUp", true, minNeffThres, isNeffThreshold_zjet);
  TH1F* zjet_cTagDown = DC.readWriteHisto(fDY, "cTagMinus/"+histSubDir, histName, sf_zjet, fout, fTT, "zjet_cTagDown", true, minNeffThres, isNeffThreshold_zjet);

  //SingleTop
  double sf_stop = 1;
  TH1F* stop = DC.readWriteHisto(fST, "base/"+histSubDir, histName, sf_stop, fout, fTT, "stop", true, minNeffThres, isNeffThreshold_stop);
  TH1F* stop_JESUp = DC.readWriteHisto(fST, "JESPlus/"+histSubDir, histName, sf_stop, fout, fTT, "stop_JESUp", true, minNeffThres, isNeffThreshold_stop);
  TH1F* stop_JESDown = DC.readWriteHisto(fST, "JESMinus/"+histSubDir, histName, sf_stop, fout, fTT, "stop_JESDown", true, minNeffThres, isNeffThreshold_stop);
  TH1F* stop_PileupUp = DC.readWriteHisto(fST, "PileupPlus/"+histSubDir, histName, sf_stop, fout, fTT, "stop_PileupUp", true, minNeffThres, isNeffThreshold_stop);
  TH1F* stop_PileupDown = DC.readWriteHisto(fST, "PileupMinus/"+histSubDir, histName, sf_stop, fout, fTT, "stop_PileupDown", true, minNeffThres, isNeffThreshold_stop);
  TH1F* stop_JERUp = DC.readWriteHisto(fST, "JERPlus/"+histSubDir, histName, sf_stop, fout, fTT, "stop_JERUp", true, minNeffThres, isNeffThreshold_stop);
  TH1F* stop_JERDown = DC.readWriteHisto(fST, "JERMinus/"+histSubDir, histName, sf_stop, fout, fTT, "stop_JERDown", true, minNeffThres, isNeffThreshold_stop);
  TH1F* stop_bTagUp = DC.readWriteHisto(fST, "bTagPlus/"+histSubDir, histName, sf_stop, fout, fTT, "stop_bTagUp", true, minNeffThres, isNeffThreshold_stop);
  TH1F* stop_bTagDown = DC.readWriteHisto(fST, "bTagMinus/"+histSubDir, histName, sf_stop, fout, fTT, "stop_bTagDown", true, minNeffThres, isNeffThreshold_stop);
  TH1F* stop_cTagUp = DC.readWriteHisto(fST, "cTagPlus/"+histSubDir, histName, sf_stop, fout, fTT, "stop_cTagUp", true, minNeffThres, isNeffThreshold_stop);
  TH1F* stop_cTagDown = DC.readWriteHisto(fST, "cTagMinus/"+histSubDir, histName, sf_stop, fout, fTT, "stop_cTagDown", true, minNeffThres, isNeffThreshold_stop);

  //Dibosons
  double sf_vv = 1;
  TH1F* vv = DC.readWriteHisto(fVV, "base/"+histSubDir, histName, sf_vv, fout, fTT, "vv", true, minNeffThres, isNeffThreshold_vv);
  TH1F* vv_JESUp = DC.readWriteHisto(fVV, "JESPlus/"+histSubDir, histName, sf_vv, fout, fTT, "vv_JESUp", true, minNeffThres, isNeffThreshold_vv);
  TH1F* vv_JESDown = DC.readWriteHisto(fVV, "JESMinus/"+histSubDir, histName, sf_vv, fout, fTT, "vv_JESDown", true, minNeffThres, isNeffThreshold_vv);
  TH1F* vv_PileupUp = DC.readWriteHisto(fVV, "PileupPlus/"+histSubDir, histName, sf_vv, fout, fTT, "vv_PileupUp", true, minNeffThres, isNeffThreshold_vv);
  TH1F* vv_PileupDown = DC.readWriteHisto(fVV, "PileupMinus/"+histSubDir, histName, sf_vv, fout, fTT, "vv_PileupDown", true, minNeffThres, isNeffThreshold_vv);
  TH1F* vv_JERUp = DC.readWriteHisto(fVV, "JERPlus/"+histSubDir, histName, sf_vv, fout, fTT, "vv_JERUp", true, minNeffThres, isNeffThreshold_vv);
  TH1F* vv_JERDown = DC.readWriteHisto(fVV, "JERMinus/"+histSubDir, histName, sf_vv, fout, fTT, "vv_JERDown", true, minNeffThres, isNeffThreshold_vv);
  TH1F* vv_bTagUp = DC.readWriteHisto(fVV, "bTagPlus/"+histSubDir, histName, sf_vv, fout, fTT, "vv_bTagUp", true, minNeffThres, isNeffThreshold_vv);
  TH1F* vv_bTagDown = DC.readWriteHisto(fVV, "bTagMinus/"+histSubDir, histName, sf_vv, fout, fTT, "vv_bTagDown", true, minNeffThres, isNeffThreshold_vv);
  TH1F* vv_cTagUp = DC.readWriteHisto(fVV, "cTagPlus/"+histSubDir, histName, sf_vv, fout, fTT, "vv_cTagUp", true, minNeffThres, isNeffThreshold_vv);
  TH1F* vv_cTagDown = DC.readWriteHisto(fVV, "cTagMinus/"+histSubDir, histName, sf_vv, fout, fTT, "vv_cTagDown", true, minNeffThres, isNeffThreshold_vv);

  //QCD MC
  double sf_qcd = 1;
  /*
  TH1F* qcd = DC.readWriteHisto(fQCD, "base/"+histSubDir, histName, sf_qcd, fout, fTT, "qcd", true);
  TH1F* qcd_JESUp = DC.readWriteHisto(fQCD, "JESPlus/"+histSubDir, histName, sf_qcd, fout, fTT, "qcd_JESUp", true);
  TH1F* qcd_JESDown = DC.readWriteHisto(fQCD, "JESMinus/"+histSubDir, histName, sf_qcd, fout, fTT, "qcd_JESDown", true);
  TH1F* qcd_PileupUp = DC.readWriteHisto(fQCD, "PileupPlus/"+histSubDir, histName, sf_qcd, fout, fTT, "qcd_PileupUp", true);
  TH1F* qcd_PileupDown = DC.readWriteHisto(fQCD, "PileupMinus/"+histSubDir, histName, sf_qcd, fout, fTT, "qcd_PileupDown", true);
  TH1F* qcd_JERUp = DC.readWriteHisto(fQCD, "JERPlus/"+histSubDir, histName, sf_qcd, fout, fTT, "qcd_JERUp", true);
  TH1F* qcd_JERDown = DC.readWriteHisto(fQCD, "JERMinus/"+histSubDir, histName, sf_qcd, fout, fTT, "qcd_JERDown", true);
  TH1F* qcd_bTagUp = DC.readWriteHisto(fQCD, "bTagPlus/"+histSubDir, histName, sf_qcd, fout, fTT, "qcd_bTagUp", true);
  TH1F* qcd_bTagDown = DC.readWriteHisto(fQCD, "bTagMinus/"+histSubDir, histName, sf_qcd, fout, fTT, "qcd_bTagDown", true);
  TH1F* qcd_cTagUp = DC.readWriteHisto(fQCD, "cTagPlus/"+histSubDir, histName, sf_qcd, fout, fTT, "qcd_cTagUp", true);
  TH1F* qcd_cTagDown = DC.readWriteHisto(fQCD, "cTagMinus/"+histSubDir, histName, sf_qcd, fout, fTT, "qcd_cTagDown", true);
  */
  //QCD data driven
  TH1F* qcd_dd = DC.readWriteHisto(fQCD_dd, "base/"+histSubDir, histName, sf_qcd, fout, fTT, "qcd", true, minNeffThres, isNeffThreshold_qcd_dd);

  //Data
  double sf_data = 1; //should be 1, always
  TH1F* data_obs = DC.readWriteHisto(fData, "base/"+histSubDir, histName, sf_data, fout, fTT, "data_obs", true);


  //wh
  double sf_wh = 1; 
  TH1F* wh = DC.readWriteHisto(fWH, "base/"+histSubDir, histName, sf_wh, fout, fTT, label, true);
  TH1F* wh_JESUp = DC.readWriteHisto(fWH, "JESPlus/"+histSubDir, histName, sf_wh, fout, fTT, label+"_JESUp", true);
  TH1F* wh_JESDown = DC.readWriteHisto(fWH, "JESMinus/"+histSubDir, histName, sf_wh, fout, fTT, label+"_JESDown", true);
  TH1F* wh_PileupUp = DC.readWriteHisto(fWH, "PileupPlus/"+histSubDir, histName, sf_wh, fout, fTT, label+"_PileupUp", true);
  TH1F* wh_PileupDown = DC.readWriteHisto(fWH, "PileupMinus/"+histSubDir, histName, sf_wh, fout, fTT, label+"_PileupDown", true);
  TH1F* wh_JERUp = DC.readWriteHisto(fWH, "JERPlus/"+histSubDir, histName, sf_wh, fout, fTT, label+"_JERUp", true);
  TH1F* wh_JERDown = DC.readWriteHisto(fWH, "JERMinus/"+histSubDir, histName, sf_wh, fout, fTT, label+"_JERDown", true);
  TH1F* wh_topPtUp = DC.readWriteHisto(fWH,  "TopPtPlus/"+histSubDir, histName, sf_wh, fout, fTT, label+"_topPtUp", true);
  TH1F* wh_topPtDown = DC.readWriteHisto(fWH, "TopPtMinus/"+histSubDir, histName, sf_wh, fout, fTT, label+"_topPtDown", true);
  TH1F* wh_bTagUp = DC.readWriteHisto(fWH, "bTagPlus/"+histSubDir, histName, sf_wh, fout, fTT, label+"_bTagUp", true);
  TH1F* wh_bTagDown = DC.readWriteHisto(fWH, "bTagMinus/"+histSubDir, histName, sf_wh, fout, fTT, label+"_bTagDown", true);
  TH1F* wh_cTagUp = DC.readWriteHisto(fWH, "cTagPlus/"+histSubDir, histName, sf_wh, fout, fTT, label+"_cTagUp", true);
  TH1F* wh_cTagDown = DC.readWriteHisto(fWH, "cTagMinus/"+histSubDir, histName, sf_wh, fout, fTT, label+"_cTagDown", true);
  
  //open input template data card of 8 TeV
  ifstream in;
  char* c = new char[1000];
  in.open("MyTemplateDataCard.txt");
  //create output data card for 13 TeV
  string outDataCard = "datacard_csbar_13TeV_WH.txt";
  string histName_str(histSubDir_+TString("_")+histName);
  if(isMuChannel) outDataCard = "datacard_csbar_mu_"+histName_str+"_13TeV_WH%d.txt"; 
  else outDataCard = "datacard_csbar_ele_"+histName_str+"_13TeV_WH%d.txt";
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
        line.replace( line.find("CCC") , 3 , string(Form("%s", string(channelName).c_str())) ); 
        out << line << endl;
      }
      else if(line.find("shapes")!=string::npos){
        line.replace( line.find("XXX") , 3 , string(TString("HplusShapes_")+channelName+TString("_")+histSubDir_+TString("_")+histName+TString("_13TeV_")+label));
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
            << space << 0*ttll->Integral()
            << space << wjet->Integral()
            << space << zjet->Integral()
            << space << stop->Integral()
            << space << vv->Integral()
            << space << qcd_dd->Integral()
            << endl;
      }
      else if(line.find("CMS_eff_b")!=string::npos){
        float bTagUnc_wh = (wh->Integral() > 0) ? DC.getBTagUnc(wh, wh_bTagUp, wh_bTagDown) : 1.00;
        line.replace( line.find("HHHH") , 4 , string(Form("%.3f", bTagUnc_wh)) );
        
        float bTagUnc_ttbar = (ttbar->Integral() > 0) ? DC.getBTagUnc(ttbar, ttbar_bTagUp, ttbar_bTagDown) : 1.00; 
        line.replace( line.find("TTTT") , 4 , string(Form("%.3f", bTagUnc_ttbar)) ); 
        
        float bTagUnc_wjet = (wjet->Integral() > 0) ? DC.getBTagUnc(wjet, wjet_bTagUp, wjet_bTagDown) : 1.00;
        line.replace( line.find("WWWW") , 4 , string(Form("%.3f", bTagUnc_wjet)) ); 
       
        float bTagUnc_zjet = (zjet->Integral() > 0) ? DC.getBTagUnc(zjet, zjet_bTagUp, zjet_bTagDown) : 1.00;
        line.replace( line.find("DDDD") , 4 , string(Form("%.3f", bTagUnc_zjet)) ); 

        float bTagUnc_stop = (stop->Integral() > 0) ? DC.getBTagUnc(stop, stop_bTagUp, stop_bTagDown) : 1.00; 
        line.replace( line.find("SSSS") , 4 , string(Form("%.3f", bTagUnc_stop)) ); 

        float bTagUnc_vv = (vv->Integral() > 0) ? DC.getBTagUnc(vv, vv_bTagUp, vv_bTagDown) : 1.00;
        line.replace( line.find("VVVV") , 4 , string(Form("%.3f", bTagUnc_vv)) ); 
        out << line << endl;
      }
      else if(line.find("CMS_eff_c")!=string::npos){
        float cTagUnc_wh = (wh->Integral() > 0) ? DC.getBTagUnc(wh, wh_cTagUp, wh_cTagDown) : 1.00;
        line.replace( line.find("HHHH") , 4 , string(Form("%.3f", cTagUnc_wh)) );
        
        float cTagUnc_ttbar = (ttbar->Integral() > 0) ? DC.getBTagUnc(ttbar, ttbar_cTagUp, ttbar_cTagDown) : 1.00; 
        line.replace( line.find("TTTT") , 4 , string(Form("%.3f", cTagUnc_ttbar)) ); 
        
        float cTagUnc_wjet = (wjet->Integral() > 0) ? DC.getBTagUnc(wjet, wjet_cTagUp, wjet_cTagDown) : 1.00;
        line.replace( line.find("WWWW") , 4 , string(Form("%.3f", cTagUnc_wjet)) ); 
       
        float cTagUnc_zjet = (zjet->Integral() > 0) ? DC.getBTagUnc(zjet, zjet_cTagUp, zjet_cTagDown) : 1.00;
        line.replace( line.find("DDDD") , 4 , string(Form("%.3f", cTagUnc_zjet)) ); 
        
        float cTagUnc_stop = (stop->Integral() > 0) ? DC.getBTagUnc(stop, stop_cTagUp, stop_cTagDown) : 1.00; 
        line.replace( line.find("SSSS") , 4 , string(Form("%.3f", cTagUnc_stop)) ); 

        float cTagUnc_vv = (vv->Integral() > 0) ? DC.getBTagUnc(vv, vv_cTagUp, vv_cTagDown) : 1.00;
        line.replace( line.find("VVVV") , 4 , string(Form("%.3f", cTagUnc_vv)) ); 
        out << line << endl;
      }
      else if(line.find("CMS_stat_wh")!=string::npos){ 
        line.replace( line.find("XXXX") , 4 , string(Form("%.3f", DC.getStatUnc(wh,  0))));  
        out << line << endl;
      } 
      else if(line.find("CMS_stat_tt")!=string::npos){  
	line.replace( line.find("XXXX") , 4 , string(Form("%.3f", DC.getStatUnc(ttbar,  0))));   
        out << line << endl;
      }  
      else if(line.find("CMS_stat_wjet")!=string::npos){  
        line.replace( line.find("XXXX") , 4 , string(Form("%.3f", DC.getStatUnc(wjet,  0))));   
        out << line << endl;
      }  
      else if(line.find("CMS_stat_zjet")!=string::npos){ 
        line.replace( line.find("XXXX") , 4 , string(Form("%.3f", DC.getStatUnc(zjet,  0))));  
        out << line << endl; 
      }
      else if(line.find("CMS_stat_stop")!=string::npos){ 
        line.replace( line.find("XXXX") , 4 , string(Form("%.3f", DC.getStatUnc(stop,  0))));  
        out << line << endl; 
      } 
      else if(line.find("CMS_stat_vv")!=string::npos){  
        line.replace( line.find("XXXX") , 4 , string(Form("%.3f", DC.getStatUnc(vv,  0))));   
        out << line << endl;  
      }
      else if(line.find("CMS_stat_qcd")!=string::npos){  
        line.replace( line.find("XXXX") , 4 , string(Form("%.3f", DC.getStatUnc(qcd_dd,  0))));   
        out << line << endl;  
      }
      else if(line.find("CMS_norm_qcd")!=string::npos){  
        if(isMuChannel) line.replace( line.find("QQQQ") , 4 , string(Form("%.3f", 1.72)));   
        else line.replace( line.find("QQQQ") , 4 , string(Form("%.3f", 1.40)));   
        out << line << endl;  
      }
      else{ //default without changes
        out << line << endl;
      }
    }
  } 
  out.close();
  in.close();
  fout->Close();
}
