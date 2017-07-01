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


//https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
//8TeV scale factors:
//const double JEREtaMap[8] = {0., 0.5, 1.1, 1.7, 2.3, 2.8, 3.2, 5.0}; 
//const double JERSF[7] = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};
//const double JERSFUp[7] = {1.105, 1.127, 1.150, 1.254, 1.316, 1.458, 1.247}; 
//const double JERSFDown[7] = {1.053, 1.071, 1.092, 1.162, 1.192, 1.332, 0.865}; 

//13TeV scale factors:
const double JEREtaMap[14] = {0., 0.5, 0.8, 1.1, 1.3, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3.0, 3.2, 5.0}; 
const double JERSF[13] = {1.109, 1.138, 1.114, 1.123, 1.084, 1.082, 1.140, 1.067, 1.177, 1.364, 1.857, 1.328, 1.16};
const double JERSFUp[13] = {1.109+0.008 , 1.138+0.013, 1.114+0.013, 1.123+0.024, 1.084+0.011, 1.082+0.035, 1.140+0.047, 1.067+0.053, 1.177+0.041, 1.364+0.039, 1.857+0.071, 1.328+0.022, 1.16+0.029};
const double JERSFDown[13] = {1.109-0.008 , 1.138-0.013, 1.114-0.013, 1.123-0.024, 1.084-0.011, 1.082-0.035, 1.140-0.047, 1.067-0.053, 1.177-0.041, 1.364-0.039, 1.857-0.071, 1.328-0.022, 1.16-0.029};

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
