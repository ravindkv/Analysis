
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
#include "TH2.h"
#include "TTimeStamp.h"
#include "Math/VectorUtil.h"

#include "interface/Reader.h"
#include "interface/ObjectSelector.hh"
#include "interface/MomentumVec.h"
#include "interface/LumiReweighting.h"
#include "interface/UncertaintyComputer.hh"
#include "interface/HistogramPlotter.hh"
#include "interface/BTagCalibrationStandalone.h"

class hplusAnalyzer : public ObjectSelector, HistogramPlotter
{
public :
  hplusAnalyzer() : ObjectSelector(), HistogramPlotter()
  {
    DRMIN_JET = 0.4;
    DRMIN_ELE = 0.4;
    METCUT_   = 30.0;
    //---------------------------------------------------//
    //Pileup reweigting 
    //---------------------------------------------------//
    //PU info for Data: 
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData
    //PU info for MC:
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/Pileup_MC_Information
    LumiWeights_ = reweight::LumiReWeighting("stack/lumiRewgt/trueInTimePU_mcDY.root","stack/lumiRewgt/trueMinBiasPU_dataMu.root", "pileup", "pileup");
    PShiftDown_ = reweight::PoissonMeanShifter(-0.5);
    PShiftUp_ = reweight::PoissonMeanShifter(0.5);
    
    //---------------------------------------------------//
    //MC cross sections at 13 TeV 
    //---------------------------------------------------//
    //https://github.com/BristolTopGroup/AnalysisSoftware/blob/master/python/DataSetInfo_13TeV.py
    //https://github.com/BristolTopGroup/AnalysisSoftware/blob/master/python/DataSetInfo_13TeV_25ns.py
    //https://indico.cern.ch/event/617002/contributions/2490586/attachments/1419016/2173704/update_27022017.pdf
    //evtDBS= event at Data Base Server i.e in DAS (https://cmsweb.cern.ch/das/).
    xss["DY1JetsToLL"]       =  1016;          evtDBS["DY1JetsToLL"]       =  62285500;
  //xss["DY1JetsToLL"]       =  1016;          evtDBS["DY1JetsToLL"]       =  62627174;
    xss["DY2JetsToLL"]       =  331.3;         evtDBS["DY2JetsToLL"]       =  19970551;
    xss["DY3JetsToLL"]       =  96.6;          evtDBS["DY3JetsToLL"]       =  5856110;
    xss["DY4JetsToLL"]       =  51.4;          evtDBS["DY4JetsToLL"]       =  4197868;
    xss["DYJetsToLL"]        =  4895;          evtDBS["DYJetsToLL"]        =  48518500;
  //xss["DYJetsToLL"]        =  4895;          evtDBS["DYJetsToLL"]        =  49144274;
    xss["HplusM100"]         =  1;             evtDBS["HplusM100"]         =  996170; 
    xss["HplusM120"]         =  1;             evtDBS["HplusM120"]         =  994498; 
    xss["HplusM140"]         =  1;             evtDBS["HplusM140"]         =  987730; 
    xss["HplusM150"]         =  1;             evtDBS["HplusM150"]         =  990645;
    xss["HplusM155"]         =  1;             evtDBS["HplusM155"]         =  952984;
    xss["HplusM160"]         =  1;             evtDBS["HplusM160"]         =  992264;
    xss["HplusM80"]          =  1;             evtDBS["HplusM80"]          =  976710;
    xss["HplusM90"]          =  1;             evtDBS["HplusM90"]          =  988480;
    xss["QCD_Pt-15to20_Mu"]  =  3819570;       evtDBS["QCD_Pt-15to20_Mu"]  =  4141251;
    xss["QCD_Pt-20to30_Mu"]  =  2960198;       evtDBS["QCD_Pt-20to30_Mu"]  =  31475157;
    xss["QCD_Pt-30to50_Mu"]  =  1652471;       evtDBS["QCD_Pt-30to50_Mu"]  =  29954815;
    xss["QCD_Pt-50to80_Mu"]  =  437504;        evtDBS["QCD_Pt-50to80_Mu"]  =  19806915;
    xss["QCD_Pt-80to120_Mu"] =  106033;        evtDBS["QCD_Pt-80to120_Mu"] =  13786971;
    xss["QCD_Pt-120to170_Mu"]=  25190;         evtDBS["QCD_Pt-120to170_Mu"]=  8042721;
    xss["QCD_Pt-170to300_Mu"]=  8654;          evtDBS["QCD_Pt-170to300_Mu"]=  7947159;
    xss["QCD_Pt-300to470_Mu"]=  797;           evtDBS["QCD_Pt-300to470_Mu"]=  7937590;
    xss["ST_s"]              =  7.3;           evtDBS["ST_s"]              =  2989199;
    xss["ST_t"]              =  136.2;         evtDBS["ST_t"]              =  38811017;
    xss["ST_tW"]             =  35.6;          evtDBS["ST_tW"]             =  6933094;
    xss["TTJets"]            =  831.76;        evtDBS["TTJets"]            =  10139950;   
    xss["W1JetsToLNu"]       =  9493;          evtDBS["W1JetsToLNu"]       =  44981200;
  //xss["W1JetsToLNu"]       =  9493;          evtDBS["W1JetsToLNu"]       =  45367044;
    xss["W2JetsToLNu"]       =  3120;          evtDBS["W2JetsToLNu"]       =  29878415;
    xss["W3JetsToLNu"]       =  942.3;         evtDBS["W3JetsToLNu"]       =  19798117;
    xss["W4JetsToLNu"]       =  524.2;         evtDBS["W4JetsToLNu"]       =  9170576;
    xss["WJetsToLNu"]        =  50690;         evtDBS["WJetsToLNu"]        =  29705748;
  //xss["WJetsToLNu"]        =  50690;         evtDBS["WJetsToLNu"]        =  29705748;
    xss["WW"]                =  63.21;         evtDBS["WW"]                =  994012;
    xss["WZ"]                =  22.82;         evtDBS["WZ"]                =  1000000;
    xss["ZZ"]                =  10.32;         evtDBS["ZZ"]                =  990064; 
    
    //Lumis(inverse pb) of single muon DATA at 13TeV
  };
  ~hplusAnalyzer() {
    delete evR;
  };
  
  void CutFlowAnalysis(TString url,  string myKey="PFlow", string evtType="data");
  void CutFlowProcessor(TString url,  string myKey="PFlow", TString cutflowType="base", TFile *outFile_=0);
  //void CreateAnalHistos(TString flowType, TFile* outFile_);
  void processEvents();
  float reweightHEPNUPWJets(int hepNUP);
  float reweightHEPNUPDYJets(int hepNUP);
private :
  double DRMIN_JET, DRMIN_ELE, METCUT_;
  Reader *evR;
  
  reweight::LumiReWeighting LumiWeights_;
  reweight::PoissonMeanShifter PShiftUp_;   //pileup syst up
  reweight::PoissonMeanShifter PShiftDown_; //pileup syst down 
  std::map<string, double> xss;
  std::map<string, double> evtDBS;
  std::map<string, double> muSF;
  /*
  BTagCalibrationReader readCSVfile(const std::string &tagger, const std::string &filename);
  */
  ofstream outfile_;
  BTagCalibrationReader readCSVfile(const std::string &filename,const std::string &tagger, BTagEntry::OperatingPoint op, const std::string & measurementType, const std::string & sysType, const std::vector<std::string> & otherSysTypes, BTagEntry::JetFlavor jf);
  Double_t getMuonSF(TH2D *h2, double eta, double pt);
  Double_t getMuonTrigSF(TH2D *h2, double eta, double pt);
  Double_t getEleSF(TH2D *h2, double etaSC, double pt);

};

float hplusAnalyzer::reweightHEPNUPWJets(int hepNUP) {

  int nJets = hepNUP-5;
  if(nJets==0)      return 2.11;
  else if(nJets==1) return 0.23;
  else if(nJets==2) return 0.119;
  else if(nJets==3) return 0.0562;
  else if(nJets>=4) return 0.0671;
  else return 1 ;
}

float hplusAnalyzer::reweightHEPNUPDYJets(int hepNUP){

  int nJets = hepNUP-5;
  if(nJets==0)      return 0.120;
  else if(nJets==1) return 0.0164;
  else if(nJets==2) return 0.0168;
  else if(nJets==3) return 0.0167;
  else if(nJets>=4) return 0.0128;
  else return 1 ;
}

BTagCalibrationReader hplusAnalyzer::readCSVfile(const std::string &filename, const std::string &tagger,  
		BTagEntry::OperatingPoint op, 
		const std::string & measurementType,
		const std::string & sysType, 
		const std::vector<std::string> & otherSysTypes, 
		BTagEntry::JetFlavor jf
		)
{ 
  BTagCalibration calib(tagger, filename);
  BTagCalibrationReader reader(op, sysType, otherSysTypes);      
  reader.load(calib, jf, measurementType); 
  return reader;
}
//https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
Double_t hplusAnalyzer::getMuonSF(TH2D *h2, double eta, double pt){
  
  TAxis *xaxis = h2->GetXaxis();
  TAxis *yaxis = h2->GetYaxis();
  //since the Pt range of 2D histo is <120
  //for Pt >120, we use SF of Pt = 120
  if(pt<=120){
    Int_t binX = xaxis->FindBin(abs(eta));
    Int_t binY = yaxis->FindBin(pt);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }
  else{
    Int_t binX = xaxis->FindBin(abs(eta));
    Int_t binY = yaxis->FindBin(120);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }	  
}
Double_t hplusAnalyzer::getMuonTrigSF(TH2D *h2, double eta, double pt){
  
  TAxis *xaxis = h2->GetXaxis();
  TAxis *yaxis = h2->GetYaxis();
  //since the Pt range of 2D histo is <120
  //for Pt >120, we use SF of Pt = 120
  if(pt<=500){
    Int_t binX = xaxis->FindBin(abs(eta));
    Int_t binY = yaxis->FindBin(pt);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }
  else{
    Int_t binX = xaxis->FindBin(abs(eta));
    Int_t binY = yaxis->FindBin(500);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }	  
}

Double_t hplusAnalyzer::getEleSF(TH2D *h2, double etaSC, double pt){
  
  TAxis *xaxis = h2->GetXaxis();
  TAxis *yaxis = h2->GetYaxis();
  //since the Pt range of 2D histo is <500
  //for Pt >500, we use SF of Pt = 500
  if(pt<=500){
    Int_t binX = xaxis->FindBin(abs(etaSC));
    Int_t binY = yaxis->FindBin(pt);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }
  else{
    Int_t binX = xaxis->FindBin(abs(etaSC));
    Int_t binY = yaxis->FindBin(500);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }	  
}



