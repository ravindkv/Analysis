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
void plotNuisanceParam_13TeV( TString INPUT_FILE  = "diffNus.root", 
	   TString FIT_PARAM = "post_fit_errs",
	   TString CHANNEL_NAME  = "mu", 
	   TString HIST_NAME  = "mjj_kfit", 
	   TString MASS = "90")
  {
  TFile f(INPUT_FILE,"READ"); 
  if(f.IsZombie() || f.IsZombie()){
    cout<<"\033[01;31mfile not found: \033[0m"<<INPUT_FILE<< endl;
  }
  else{
    //gStyle->SetFrameLineWidth(3);
    gStyle->SetOptStat(0);
    TCanvas *c1 = (TCanvas*)f.Get(FIT_PARAM);
    //gPad->SetBottomMargin(0.35);
    
    TPaveText *pt = new TPaveText(0.35,0.90,0.50,0.95, "brNDC"); // good_v1
    pt->SetTextSize(0.02);
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->AddText(CHANNEL_NAME);
    pt->AddText(HIST_NAME);
    pt->AddText("M_{H^{+}}= "+MASS+" GeV");
    //pt->AddText(CHANNEL_NAME+"_"+HIST_NAME+"_"+MASS);
    pt->Draw();
    c1->Print(FIT_PARAM+"_"+CHANNEL_NAME+"_"+HIST_NAME+"_"+MASS+".pdf");
  }
}


