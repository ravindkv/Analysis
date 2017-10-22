//#include "TROOT.h"    
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <math.h>
#include<string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <fstream>


// %%%%%%%%%%%% IMPORT THE FILE, PLOT THE HISTOGRAMS %%%%%%%%%%%%%%%

void getEvtProcessed(){

  bool isMC = true;
  bool isMuData = true;
  bool isEleData = false;

  TString mc_file[42] = {
  " DY1JetsToLL_Merged.root ", 
  " DY2JetsToLL_Merged.root ", 
  " DY3JetsToLL_Merged.root ",
  " DY4JetsToLL_Merged.root ", 
  " DYJetsToLL_Merged.root  ", 
  " HplusM100_Merged.root   ", 
  " HplusM120_Merged.root   ", 
  " HplusM140_Merged.root   ", 
  " HplusM150_Merged.root   ", 
  " HplusM155_Merged.root   ", 
  " HplusM160_Merged.root   ", 
  " HplusM80_Merged.root    ", 
  " HplusM90_Merged.root    ", 
  " QCD_Pt-15to20_Mu_Merged.root   ", 
  " QCD_Pt-20to30_Mu_Merged.root   ",
  " QCD_Pt-30to50_Mu_Merged.root   ",
  " QCD_Pt-50to80_Mu_Merged.root   ",
  " QCD_Pt-80to120_Mu_Merged.root  ",    
  " QCD_Pt-120to170_Mu_Merged.root ", 
  " QCD_Pt-170to300_Mu_Merged.root ", 
  " QCD_Pt-300to470_Mu_Merged.root ", 
  " ST_s_Merged.root       ",
  " ST_t__Merged.root      ", 
  " ST_tW_Merged.root      ", 
  " TTJetsM_Merged.root     ", 
  " TTJetsP_Merged.root     ", 
  " W1JetsToLNu_Merged.root", 
  " W2JetsToLNu_Merged.root", 
  " W3JetsToLNu_Merged.root", 
  " W4JetsToLNu_Merged.root", 
  " WJetsToLNu_Merged.root ", 
  " WW_Merged.root ", 
  " WZ_Merged.root ",
  " ZZ_Merged.root ", 
  " all_DY.root    ",
  " all_Hplus.root ",
  " all_QCD.root   ",
  " all_ST.root    ",
  " all_TTJetsM.root",
  " all_TTJetsP.root",
  " all_VV.root    ",
  " all_WJets.root ",
  };
  
  TString mu_data_file[9] = {
  " MuRunB2v2_Merged.root   ", 
  " MuRunCv1_Merged.root    ", 
  " MuRunDv1_Merged.root    ", 
  " MuRunEv1_Merged.root    ", 
  " MuRunFv1_Merged.root    ", 
  " MuRunGv1_Merged.root    ", 
  " MuRunH2v1_Merged.root   ", 
  " MuRunH3v1_Merged.root   ", 
  " all_muData.root"
  };
  
  TString ele_data_file[9] = {
  " EleRunBver2v2_Merged.root", 
  " EleRunCv1_Merged.root    ", 
  " EleRunDv1_Merged.root    ", 
  " EleRunEv1_Merged.root    ", 
  " EleRunFv1_Merged.root    ", 
  " EleRunGv1_Merged.root    ", 
  " EleRunHver2v1_Merged.root", 
  " EleRunHver3v1_Merged.root",
  " all_EleData.root"
  };
 
  if(isMC){ 
    cout<<endl;
    cout<<"=============================="<<endl;
    cout<<"        ALL MC SAMPLES        "<<endl;
    cout<<"=============================="<<endl;
    for(int i= 0; i<42; i++){
      TString inFile(mc_file[i]);
      TFile* ttbar= new TFile(inFile);
      TString path= "base/totalEvents";
      TH1F* hist= (TH1F*)(ttbar->Get(path));
      int entries= hist->GetBinContent(1);
      double mean= hist->GetMean();
      double events = entries*mean;
      cout<<events<<" : "<<"\t"<<inFile<<endl;
      //cout<<inFile<<endl;
      //cout<<" Total Events= "<<events<<endl;
      //hist->Draw();
    }
  }

  if(isMuData){ 
    cout<<endl;
    cout<<"=============================="<<endl;
    cout<<" ALL MUON DATA SAMPLES        "<<endl;
    cout<<"=============================="<<endl;
    for(int i= 0; i<9; i++){
      TString inFile(mu_data_file[i]);
      TFile* ttbar= new TFile(inFile);
      TString path= "base/totalEvents";
      TH1F* hist= (TH1F*)(ttbar->Get(path));
      int entries= hist->GetBinContent(1);
      double mean= hist->GetMean();
      double events = entries*mean;
      cout<<events<<" : "<<"\t"<<inFile<<endl;
      //cout<<inFile<<endl;
      //cout<<" Total Events= "<<events<<endl;
      //hist->Draw();
    }
  }

  if(isEleData){ 
    cout<<endl;
    cout<<"=============================="<<endl;
    cout<<" ALL ELECTRON DATA SAMPLES    "<<endl;
    cout<<"=============================="<<endl;
    for(int i= 0; i<9; i++){
      TString inFile(ele_data_file[i]);
      TFile* ttbar= new TFile(inFile);
      TString path= "base/totalEvents";
      TH1F* hist= (TH1F*)(ttbar->Get(path));
      int entries= hist->GetBinContent(1);
      double mean= hist->GetMean();
      double events = entries*mean;
      cout<<events<<" : "<<"\t"<<inFile<<endl;
      //cout<<" Total Events= "<<events<<endl;
      //hist->Draw();
    }
  }
}
