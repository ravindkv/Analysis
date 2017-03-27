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
#include "interface/BtagSF.hh"
#include "interface/SVEffUnc.hh"

#endif

const double JEREtaMap[8] = {0., 0.5, 1.1, 1.7, 2.3, 2.8, 3.2, 5.0}; 
const double JERSF[7] = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};
const double JERSFUp[7] = {1.105, 1.127, 1.150, 1.254, 1.316, 1.458, 1.247}; 
const double JERSFDown[7] = {1.053, 1.071, 1.092, 1.162, 1.192, 1.332, 0.865}; 

class UncertaintyComputer{

public :
  UncertaintyComputer()
  {
    btsf = new BtagSF(12345);
    sveffunc = new SVEffUnc();
  }

   virtual ~UncertaintyComputer(){
   ///~UncertaintyComputer(){
     delete btsf;
     delete sveffunc;
  }
  
  double metWithJES(const vector<MyJet> & vJ, vector<int> *j, MyMET MET, int jes=0);
  double metWithJER(const vector<MyJet> & vJ, vector<int> *j, MyMET MET, int jer=0);
  double metWithJESJER(const vector<MyJet> & vJ, vector<int> *j, MyMET MET, int jes=0, int jer=0);
  double metWithUncl(const vector<MyJet> & vJ, vector<int> *j, const vector<MyMuon> &vMu, vector<int> *m, const vector<MyElectron> &vEle, vector<int> *el, MyMET MET, int unc=0);
  double getJERSF(double eta, int jer=0);
  double jetPtWithJESJER(MyJet jet, int jes=0, int jer=0); 
  bool getBtagWithSF(MyJet jet, bool isData, int scale, bool is2012);
  double EffUncOnSV(MyJet jet);
  
private :
  BtagSF* btsf;
  enum BVariation{kNo = 0, kDown = 1, kUp = 2};
  SVEffUnc* sveffunc;

  ClassDef(UncertaintyComputer, 1)
};
#endif
