#ifndef _btagsf_hh_
#define _btagsf_hh_

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TRandom3.h"
#include "TH2.h"
#include "TH2D.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "BTagCalibrationStandalone.h"
#endif
using namespace std;
class BtagSF{
  
public:
  BtagSF( int seed=0 )
  { 
    randm = new TRandom3(seed);};
  
  virtual ~BtagSF() {
    delete randm;
  };
  ///~BtagSF() {delete randm;};
  
  Bool_t isbtagged( BTagCalibrationReader &reader, TH2D *h2_BTaggingEff_Num, TH2D *h2_BTaggingEff_Denom, Float_t eta, Float_t pt, Float_t csv, Int_t jetflavor, Bool_t isdata, UInt_t btagsys); 
  Double_t getBTagSFb(BTagCalibrationReader &reader, Float_t eta, Float_t pt, Float_t csv, UInt_t btagsys);
  Double_t getBTagSFc(BTagCalibrationReader &reader, Float_t eta, Float_t pt, Float_t csv, UInt_t btagsys);
  Double_t getBTagSFl(BTagCalibrationReader &reader, Float_t eta, Float_t pt, Float_t csv, UInt_t btagsys);
Double_t getBTagEff(TH2D *h2_BTaggingEff_Num, TH2D *h2_BTaggingEff_Denom, Float_t pt, Float_t eta);

  enum { kNo, kDown, kUp };                     // systematic variations 


private:
 
  TRandom3* randm;
  ClassDef(BtagSF, 1)
};
#endif
