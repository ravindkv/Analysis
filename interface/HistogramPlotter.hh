#ifndef _histogramplotter_h_ 
#define _histogramplotter_h_ 
 
#if !defined(__CINT__) || defined(__MAKECINT__) 
 
#include <iomanip> 
#include <iostream> 
#include <fstream> 
#include <map>
 
#include "TRandom2.h" 
#include "TMatrixD.h" 
#include "TF1.h" 
#include "TProfile.h" 
#include "TObjArray.h" 
#include "TMatrixD.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TTimeStamp.h" 
#include "TFile.h"
#include <exception> 
 
#endif 
 
class HistogramPlotter{ 

public :
  HistogramPlotter(){}
  virtual ~HistogramPlotter(){}
  ///~HistogramPlotter(){}
  void addHisto(const TString & name, const TString & dirname, const int & range, const double & min, const double & max);
  void addHisto2D(const TString & name, const TString & dirname, const int & range1, const double & min1, const double & max1, const int & range2, const double & min2, const double & max2);
  void fillHisto2d(const TString & name, const TString & dirname, const double & value1, const double & value2, const double & weight=1.0);
  void initHisto(TFile *file, const TString & dir, const TString & subdir, const TString & histName, const int & Nbin, const double & min, const double & max);
  void initHisto2D(TFile *file, const TString & dir, const TString & subdir, const TString & histName, const int & nBin1, const double & min1, const double & max1, const int & nBin2, const double & min2, const double & max2);
  void fillHisto2D(TFile *file, const TString & dir, const TString & subdir, const TString & histName, const int & Nbin1, const double & min1, const double & max1, const double & value1, const int & Nbin2, const double & min2, const double & max2, const double & value2, const double & weight);
  void fillHisto(TFile *file, const TString & dir, const TString & subdir, const TString & histName, const int & Nbin, const double & min, const double & max, const double & value, const double & weight);
  TH1* getHisto(const TString & name, const TString & dirname);
  TH2* getHisto2d(const TString & name, const TString & dirname);

  // TProfile
  void initTProfile(TFile *file, const TString & dir, const TString & subdir, const TString & histName, const int & nXBin, Float_t xBin_array[]);
  void addTProfile(const TString & name, const TString & dirname, const int & nXBin, Float_t xBin_array[]);
  void fillTProfile(TFile *file, const TString & dir, const TString & subdir, const TString & histName, const int & nXBin, Float_t xBin_array[], const double & value1, const double & value2, const double & weight);

private :
  std::map<TString, TH1*> histos1_;
  std::map<TString, TH2*> histos2_;
  std::map<TString, TProfile*> tprofile_;
  //static std::map<TString, histos1_> histDir1_;
  //static std::map<TString, histos2_> histDir2_;
  //std::map<TString, TDirectory*> hDir_;
  ClassDef(HistogramPlotter, 1)
};
#endif



