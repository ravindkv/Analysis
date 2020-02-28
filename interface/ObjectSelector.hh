#ifndef _objectselector_h_
#define _objectselector_h_

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iomanip>
#include <iostream>
#include <fstream>

#include "TRandom2.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TTimeStamp.h"
#include <exception>

#ifdef _STANDALONE
#include "Reader.h"
#else
#include "interface/Reader.h"
#endif
#include "interface/UncertaintyComputer.hh"
#endif

class ObjectSelector : public UncertaintyComputer
{
public : 
  ObjectSelector() : UncertaintyComputer()
  {
  }
  virtual ~ObjectSelector(){}
  ///~ObjectSelector(){}

  // preselection of objects
  void preSelectElectrons(vector<int> * e_i, const vector<MyElectron> & vE , MyVertex & vertex, const bool & isPFlow=false);
  void preSelectMuons(TString url, vector<int> * m_i, const vector<MyMuon> & vM , MyVertex & vertex, const bool & isData=false, const double & random_u1=0.5, const double & random_u2 =0.4, const int & err_set=0, const int & err_member=0);
  void preSelectJets( string jetAlgo, vector<int> * j_i, const vector<MyJet> & vJ, const int & jes=0, const int & jer=0, const bool & isData=false);

  //Loose Lepton veto
  bool looseElectronVeto(unsigned long selectedElectron, const vector<MyElectron> & vE, MyVertex & vertex, const bool & isPFlow=false);
  bool looseMuonVeto( const int & selectedMuon, const vector<MyMuon> & vM, const bool & isPFlow=false);
  
  //Medium muon ID
  bool isMediumMuon(const MyMuon * m);
  bool isMediumMuonGH(const MyMuon * m);
  //egmGsfElectronIDs: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  bool cutBasedElectronID_Summer16_80X_V1_veto(const MyElectron *e); 
  bool cutBasedElectronID_Summer16_80X_V1_medium(const MyElectron *e, MyVertex & vertex); 

  // object cleaning
  void ElectronCleaning( const vector<MyElectron> & vE, const vector<MyMuon> & vM, vector<int> * e_old, vector<int> * e_new, vector<int> * mu, const double & DR );

  void JetCleaning(const vector<MyJet> & vJ, const vector<MyMuon> & vM, const vector<MyElectron> & vE, vector<int> * j_old, vector<int> * j_new, vector<int> * mu, vector<int> * el, const double & DR);
  double DeltaR(const MyLorentzVector & aV, const MyLorentzVector & bV); 
    
  
private :
  ClassDef(ObjectSelector, 1)
};
#endif
