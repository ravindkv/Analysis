#ifndef _uncertaintycomputer_h_
#define _uncertaintycomputer_h_

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
#include "interface/BTagSF.hh"
#include "interface/CTagSF.hh"
#include "BTagCalibrationStandalone.h"

#endif
class UncertaintyComputer{
public :
  UncertaintyComputer()
  {
    btsf = new BTagSF(12345);
    ctsf = new CTagSF(12345);
  }
   virtual ~UncertaintyComputer(){
   ///~UncertaintyComputer(){
     delete btsf;
     delete ctsf;
  }
  double muPtWithRochCorr(const MyMuon *mu, bool isData=false, double u1=0.5, double u2=0.4, int s=0, int m=0); 
  double metWithJESJER(const vector<MyJet> & vJ, vector<int> *j, MyMET MET, int jes=0, int jer=0, bool isData = false);
  void  openCSVfile(const std::string &filename); 
  double getJERSF(double eta, int jer);
  double DeltaR(MyLorentzVector aV, MyLorentzVector bV);
  
  // bTag SF, by event reweighting
  double getBTagPmcSys(TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, MyJet jet);
  double getBTagPdataSys(BTagCalibrationReader &reader, TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, MyJet jet, int scale);
  double getIncCTagPmcSys(TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, MyJet jet, bool isCTag);
  double getIncCTagPdataSys(BTagCalibrationReader &reader, TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, MyJet jet, bool isCTag, int scale);
  
private :
  BTagSF* btsf;
  CTagSF* ctsf;
  ClassDef(UncertaintyComputer, 1)
};
#endif
