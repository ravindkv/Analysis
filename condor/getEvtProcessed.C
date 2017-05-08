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

    TString allFile[37] = {
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
    " MuRunB2v2_Merged.root   ", 
    " MuRunCv1_Merged.root    ", 
    " MuRunDv1_Merged.root    ", 
    " MuRunEv1_Merged.root    ", 
    " MuRunFv1_Merged.root    ", 
    " MuRunGv1_Merged.root    ", 
    " MuRunH2v1_Merged.root   ", 
    " MuRunH3v1_Merged.root   ", 
    " QCD_Pt-120to_Merged.root", 
    " QCD_Pt-15to2_Merged.root", 
    " QCD_Pt-170to_Merged.root", 
    " QCD_Pt-20to3_Merged.root", 
    " ST_s_Merged.root       " ,
    " ST_t__Merged.root       ", 
    " ST_tW_Merged.root      ", 
    " TTJets_Merged.root     ", 
    " W1JetsToLNu_Merged.root", 
    " W2JetsToLNu_Merged.root", 
    " W3JetsToLNu_Merged.root", 
    " W4JetsToLNu_Merged.root", 
    " WJetsToLNu_Merged.root ", 
    " WW_Merged.root ", 
    " WZ_Merged.root ",
    " ZZ_Merged.root " 
    }

    for(int i= 0; i<37; i++){
      TString inFile(allFile[i]);
      TFile* ttbar= new TFile(inFile);
      TString path= "base/totalEvents";
      TH1F* hist= (TH1F*)(ttbar->Get(path));
      int entries= hist->GetBinContent(1);
      double mean= hist->GetMean();
      double events = entries*mean;
      cout<<"-----------------------"<<endl;
      cout<<inFile<<endl;
      cout<<"Total events= "<<events<<endl;
      //hist->Draw();
    }
}
