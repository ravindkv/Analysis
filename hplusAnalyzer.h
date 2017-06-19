
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
class hplusAnalyzer : public ObjectSelector, HistogramPlotter
{
public :
  hplusAnalyzer() : ObjectSelector(), HistogramPlotter()
  {
    DRMIN_JET = 0.5;
    DRMIN_ELE = 0.5;
    METCUT_   = 30.0;
    //---------------------------------------------------//
    //Pileup reweigting 
    //---------------------------------------------------//
    //PU info for Data: 
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData
    //PU info for MC:
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/Pileup_MC_Information
    LumiWeights_ = reweight::LumiReWeighting("stack/trueInTimePU_mcDY.root","stack/trueMinBiasPU_dataMu.root", "pileup", "pileup");
    PShiftDown_ = reweight::PoissonMeanShifter(-0.5);
    PShiftUp_ = reweight::PoissonMeanShifter(0.5);
    
    //---------------------------------------------------//
    //MC cross sections at 13 TeV 
    //---------------------------------------------------//
    //https://github.com/BristolTopGroup/AnalysisSoftware/blob/master/python/DataSetInfo_13TeV.py
    //https://github.com/BristolTopGroup/AnalysisSoftware/blob/master/python/DataSetInfo_13TeV_25ns.py
    //https://indico.cern.ch/event/617002/contributions/2490586/attachments/1419016/2173704/update_27022017.pdf
    //evtDBS= event at Data Base Server i.e in DAS (https://cmsweb.cern.ch/das/).
    xss["DY1JetsToLL"]       =  1016;          evtDBS["DY1JetsToLL"]       =  62472600;
    //xss["DY1JetsToLL"]       =  1016;          evtDBS["DY1JetsToLL"]       =  62627174;
    xss["DY2JetsToLL"]       =  331.3;         evtDBS["DY2JetsToLL"]       =  19970551;
    xss["DY3JetsToLL"]       =  96.6;          evtDBS["DY3JetsToLL"]       =  5856110;
    xss["DY4JetsToLL"]       =  51.4;          evtDBS["DY4JetsToLL"]       =  4197868;
    xss["DYJetsToLL"]        =  4895;          evtDBS["DYJetsToLL"]        =  48906200;
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
    xss["TTJets"]            =  831.76;        evtDBS["TTJets"]            =  9676510;   
    //xss["TTJets"]            =  831.76;        evtDBS["TTJets"]            =  10139950;   
    xss["W1JetsToLNu"]       =  9493;          evtDBS["W1JetsToLNu"]       =  45230942;
    //xss["W1JetsToLNu"]       =  9493;          evtDBS["W1JetsToLNu"]       =  45367044;
    xss["W2JetsToLNu"]       =  3120;          evtDBS["W2JetsToLNu"]       =  29878415;
    xss["W3JetsToLNu"]       =  942.3;         evtDBS["W3JetsToLNu"]       =  19798117;
    xss["W4JetsToLNu"]       =  524.2;         evtDBS["W4JetsToLNu"]       =  9170576;
    xss["WJetsToLNu"]        =  50690;         evtDBS["WJetsToLNu"]        =  29531700;
    //xss["WJetsToLNu"]        =  50690;         evtDBS["WJetsToLNu"]        =  29705748;
    xss["WW"]                =  63.21;         evtDBS["WW"]                =  994012;
    xss["WZ"]                =  22.82;         evtDBS["WZ"]                =  1000000;
    xss["ZZ"]                =  10.32;         evtDBS["ZZ"]                =  990064; 
    
    //Lumis(inverse pb) of single muon DATA at 13TeV
    //https://docs.google.com/spreadsheets/d/1lQyfcY0gnG_IgFrtnbBES_HV1M7ARQM9qCw01vsnxSk/edit?usp=sharing
    double lumiB = 5403; 
    double lumiC = 2395;
    double lumiD = 4255; 
    double lumiE = 4053;
    double lumiF = 3105;
    double lumiG = 7544;
    double lumiH = 8529+216;
    double lumiTotal = lumiB+ lumiC+ lumiD+ lumiE+ lumiF+ lumiG+ lumiH;
    
    //muon Trigger/ID/ISo SFs, in bins of eta (from muon POG)
    //SFs for different lumi period are weighted by lumi fraction.
    //Trigger SF for HLT_IsoMu24_eta2p1
    double sfEta1 = (lumiE*0.956+lumiB*0.9798+lumiC*0.9841+lumiD*0.98151)/lumiTotal; // 0<|eta|<0.9 
    double sfEta2 = (lumiE*0.9528+lumiB*0.9618+lumiC*0.9688+lumiD*0.96156)/lumiTotal; // 0.9<|eta|<1.2
    double sfEta3 = (lumiE*0.9809+lumiB*0.9814+lumiC*1.0021+lumiD*0.99721)/lumiTotal; // 1.2<|eta|<2.1
    //multiply mu ID/Iso SFs
    sfEta1 = sfEta1*0.9939*1.0004;
    sfEta2 = sfEta2*0.9902*1.0031;
    sfEta3 = sfEta3*0.9970*1.0050;
    muSF["sfEta1"] = 1; //sfEta1;
    muSF["sfEta2"] = 1; //sfEta2;
    muSF["sfEta3"] = 1; //sfEta3;
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
  ofstream outfile_;
};
