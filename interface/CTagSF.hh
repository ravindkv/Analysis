#ifndef _ctagsf_hh_
#define _ctagsf_hh_

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
class CTagSF{
  
public:
  CTagSF( int seed=0 )
  { 
    randm = new TRandom3(seed);};
  
  virtual ~CTagSF() {
    delete randm;
  };
  ///~CTagSF() {delete randm;};
  
  Bool_t isCtagged( BTagCalibrationReader &reader, TH2D *h2_CTagEff_Num, TH2D *h2_CTagEff_Denom, Float_t csv_OP, Float_t csv, Float_t eta, Float_t pt, Int_t jetflavor, Bool_t isdata, UInt_t CTagSys); 
  Double_t getCTagSFb(BTagCalibrationReader &reader, Float_t eta, Float_t pt, Float_t csv, UInt_t CTagSys);
  Double_t getCTagSFc(BTagCalibrationReader &reader, Float_t eta, Float_t pt, Float_t csv, UInt_t CTagSys);
  Double_t getCTagSFl(BTagCalibrationReader &reader, Float_t eta, Float_t pt, Float_t csv, UInt_t CTagSys);
Double_t getCTagEff(TH2D *h2_CTagEff_Num, TH2D *h2_CTagEff_Denom, Float_t pt, Float_t eta);

  enum { kNo, kDown, kUp};                     // systematic variations 

private:
 
  TRandom3* randm;
  ClassDef(CTagSF, 1)
};
#endif
