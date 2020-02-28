#ifndef _reader_h_
#define _reader_h_

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TROOT.h"
#include <string>
#include <iostream>
#include <map>
#include <algorithm>

#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TObject.h"
#include "TString.h"
#include "TObjArray.h"
#include "TH2.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TVectorD.h"
#include "TFractionFitter.h"

#ifdef _STANDALONE
#include "MyEvent.h"
#else
#include "interface/MyEvent.h"
#endif

#endif

using namespace std;
class Reader 
{
 public:
  Reader() { }
 virtual ~Reader() { }
 ///~Reader() { }
  
  unsigned int AssignEventTreeFrom(TFile *f, const TString & path = "myMiniTreeProducer");
  unsigned int AssignEventTreeFromList(const char *file_list_name, const TString & path = "myMiniTreeProducer"); //for more than one files
  
  MyEvent *GetNewEvent(unsigned int ientry);
  MyEvent *GetNewEventFromList(unsigned int ientry); //for more than one files

  Long64_t LoadTree(Long64_t entry);
  vector<MyJet> getJets(MyEvent* ev, const string & algo, const int & jes, const int & jer, const bool & isData);
  vector<MyJet> getJetsNoCorr(MyEvent* ev, const string & algo);
  vector<MyElectron> getElectrons(MyEvent* ev, const string & algo);
  vector<MyMuon> getMuons(MyEvent* ev, const string & algo);
  MyMET getMET(MyEvent* ev, const string & algo);

  double getJERSF(const double & eta, const int & jer);
  double DeltaR(const MyLorentzVector & aV, const MyLorentzVector & bV);
 private :
  TTree *myTree;
  //TChain *chain;
  MyEvent *myEvent;
  
  ClassDef(Reader,1)
};
#endif
