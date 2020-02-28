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
  double getIncCTagPmc(TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, const double & eta, const double & pt, const bool & isCTag);
  double getIncCTagPdata(BTagCalibrationReader &reader, TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, const double & eta, const double & pt, const double & csv, const bool & isCTag, const int & jetFlavor, const int & cTagSys);
  
  // For both
  double getCTagSF(BTagCalibrationReader &reader, const double & eta, const double & pt, const double & csv, const double & jetflavor, const int & CTagSys);
  double getCTagEff(TH2D *h2_CTagEff_Num, TH2D *h2_CTagEff_Denom, const double & eta, const double & pt);
  
private:
 
  TRandom3* randm;
  ClassDef(CTagSF, 1)
};
#endif
