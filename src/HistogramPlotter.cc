#include "interface/HistogramPlotter.hh" 
#include <iostream> 
#include <iomanip> 
#include <math.h>
 
ClassImp(HistogramPlotter) 

void HistogramPlotter::addHisto(const TString & name, const TString & dirname, const int & range, const double & min, const double & max)
{
  //TString fullname = name+"_"+dirname; 
  TString fullname = dirname+"/"+name;
  std::string hname(fullname);
  histos1_[fullname] = new TH1F(name.Data(), hname.c_str(), range, min, max); 
  histos1_[fullname]->Sumw2();
}

void HistogramPlotter::initTProfile(TFile *file, const TString & dir, const TString & subdir, const TString & histName, const int & nXBin, Float_t xBin_array[]){
  std::string histPath;
  histPath = std::string(dir+"/"+subdir);
  TDirectory *d = file->GetDirectory(histPath.c_str());
  if(!d) file->mkdir(histPath.c_str());
  file->cd(histPath.c_str());
  addTProfile(histName, histPath, nXBin, xBin_array);
}

void HistogramPlotter::addTProfile(const TString & name, const TString & dirname, const int & nXBin, Float_t xBin_array[]){
  //TString fullname = name+"_"+dirname;
  TString fullname = dirname+"/"+name;
  std::string hname(fullname); 
  tprofile_[fullname] = new TProfile(name.Data(), hname.c_str(), nXBin, xBin_array);
  tprofile_[fullname]->Sumw2();
  tprofile_[fullname]->SetErrorOption("s");
}

void HistogramPlotter::fillTProfile(TFile *file, const TString & dir, const TString & subdir, const TString & histName, const int & nXBin, Float_t xBin_array[], const double & value1, const double & value2, const double & weight)
{
  TString fullname = dir+"/"+subdir+"/"+histName;
  if(!tprofile_[fullname]){
    initTProfile(file, dir, subdir, histName, nXBin, xBin_array); 
  }
  //then, fill TProfile
  tprofile_[fullname]->Fill(value1, value2, weight);
}

void HistogramPlotter::addHisto2D(const TString & name, const TString & dirname, const int & nBin1, const double & min1, const double & max1, const int & nBin2, const double & min2, const double & max2)
{
  //TString fullname = name+"_"+dirname;
  TString fullname = dirname+"/"+name;
  std::string hname(fullname); 
  histos2_[fullname] = new TH2D(name.Data(), hname.c_str(), nBin1, min1, max1, nBin2, min2, max2);
  histos2_[fullname]->Sumw2();
}

void HistogramPlotter::initHisto(TFile *file, const TString & dir, const TString & subdir, const TString & histName, const int & Nbin, const double & min, const double & max){
  std::string histPath;
  histPath = std::string(dir+"/"+subdir);
  TDirectory *d = file->GetDirectory(histPath.c_str());
  if(!d) file->mkdir(histPath.c_str());
  file->cd(histPath.c_str());
  addHisto(histName, histPath, Nbin, min, max);
  }

void HistogramPlotter::initHisto2D(TFile *file, const TString & dir, const TString & subdir, const TString & histName, const int & nBin1, const double & min1, const double & max1, const int & nBin2, const double & min2, const double & max2){
  std::string histPath;
  histPath = std::string(dir+"/"+subdir);
  TDirectory *d = file->GetDirectory(histPath.c_str());
  if(!d) file->mkdir(histPath.c_str());
  file->cd(histPath.c_str());
  addHisto2D(histName, histPath, nBin1, min1, max1, nBin2, min2, max2);
}

void HistogramPlotter::fillHisto(TFile *file, const TString & dir, const TString & subdir, const TString & histName, const int & Nbin, const double & min, const double & max, const double & value, const double & weight)
{
  TString fullname = dir+"/"+subdir+"/"+histName;
  if(!histos1_[fullname]){
    initHisto(file, dir, subdir, histName, Nbin, min,  max); 
  }
  //then, fill the histogram
  histos1_[fullname]->Fill(value, weight);
}

void HistogramPlotter::fillHisto2D(TFile *file, const TString & dir, const TString & subdir, const TString & histName, const int & nBin1, const double & min1, const double & max1, const double & value1, const int & nBin2, const double & min2, const double & max2, const double & value2, const double & weight)
{
  TString fullname = dir+"/"+subdir+"/"+histName;
  if(!histos2_[fullname]){
    initHisto2D(file, dir, subdir, histName, nBin1, min1, max1, nBin2, min2, max2); 
  }
  //then, fill the histogram
  histos2_[fullname]->Fill(value1, value2, weight);
}


TH2* HistogramPlotter::getHisto2d(const TString & name, const TString & dirname)
{
  TString fullname = dirname+"/"+name;
  TH2 * h = 0;
  if(histos2_.find(fullname) != histos2_.end())h = histos2_[fullname];
  return h;
}

TH1* HistogramPlotter::getHisto(const TString & name, const TString & dirname)
{
  //TString fullname = name+"_"+dirname;
  TString fullname = dirname+"/"+name;
  TH1 * h = 0;
  if(histos1_.find(fullname) != histos1_.end())h = histos1_[fullname];
  else if(histos2_.find(fullname) != histos2_.end())h = histos2_[fullname];
  return h;
}

