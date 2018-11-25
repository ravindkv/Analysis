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
  
  //For inclusive scale factors
  double getIncCTagPmc(TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, double eta, double pt, bool isCTag);
  double getIncCTagPdata(BTagCalibrationReader &reader, TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, double eta, double pt, double csv, bool isCTag, int jetFlavor, int cTagSys);
  //For exclusive scale factors
  double getExCTagPmc(TH2D *h2_qTagEff_NumL, TH2D *h2_qTagEff_NumM, TH2D *h2_qTagEff_NumT, TH2D *h2_qTagEff_Denom, double eta, double pt, double csv, bool isCTagLMT, bool isCTagMT, bool isCTagT);
  double getExCTagPdata(BTagCalibrationReader &readerL, BTagCalibrationReader &readerM, BTagCalibrationReader &readerT, TH2D *h2_qTagEff_NumL, TH2D *h2_qTagEff_NumM, TH2D *h2_qTagEff_NumT, TH2D *h2_qTagEff_Denom, double eta, double pt, double csv, int jetFlavor, int btagsys, bool isCTagLMT, bool isCTagMT, bool isCTagT);
  
  // For both
  double getCTagSF(BTagCalibrationReader &reader, double eta, double pt, double csv, double jetflavor, int CTagSys);
  double getCTagEff(TH2D *h2_CTagEff_Num, TH2D *h2_CTagEff_Denom, double eta, double pt);
  
  // For systematics
  enum { kNo, kDown, kUp};  

private:
 
  TRandom3* randm;
  ClassDef(CTagSF, 1)
};
#endif
