#include "interface/CTagSF.hh"

ClassImp(CTagSF)

double CTagSF::getIncCTagPmc(TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, double eta, double pt, bool isCTag){
  double pMC = 1.0; 
  if(isCTag) pMC = getCTagEff(h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt);
  else pMC = 1 - getCTagEff(h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt);
  return pMC;
}

double CTagSF::getIncCTagPdata(BTagCalibrationReader &reader, TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, double eta, double pt, double csv, bool isCTag, int jetFlavor, int cTagSys){
  double pData = 1.0; 
  double sf = 1.0;
  double eff = 1.0;
  if(isCTag){//tagged 
    sf = getCTagSF(reader, eta, pt, csv, jetFlavor, cTagSys);
    eff = getCTagEff(h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt);
    pData = sf*eff;
  }
  else{// untagged
    sf = getCTagSF(reader, eta, pt, csv, jetFlavor, cTagSys);
    eff = getCTagEff(h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt);
    pData = 1.0 - sf*eff;
  }
  return pData;
}

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods#1a)%20Event%20reweighting%20using%20scal
double CTagSF::getExCTagPmc(TH2D *h2_qTagEff_NumL, TH2D *h2_qTagEff_NumM, TH2D *h2_qTagEff_NumT, TH2D *h2_qTagEff_Denom, double eta, double pt, double csv, bool isCTagLMT, bool isCTagMT, bool isCTagT){
  double pMC = 1.0;
  if(isCTagT){
    pMC = getCTagEff(h2_qTagEff_NumT, h2_qTagEff_Denom, eta, pt);
  }
  else if(isCTagMT){
    double effM = getCTagEff(h2_qTagEff_NumM, h2_qTagEff_Denom, eta, pt);
    double effT = getCTagEff(h2_qTagEff_NumT, h2_qTagEff_Denom, eta, pt);
    pMC = effM - effT; 
  }
  else if(isCTagLMT){
    double effL = getCTagEff(h2_qTagEff_NumL, h2_qTagEff_Denom, eta, pt);
    double effM = getCTagEff(h2_qTagEff_NumM, h2_qTagEff_Denom, eta, pt);
    double effT = getCTagEff(h2_qTagEff_NumT, h2_qTagEff_Denom, eta, pt);
    pMC = (effL - effM)* (effM - effT);
  }
  return pMC;
}

double CTagSF::getExCTagPdata(BTagCalibrationReader &readerL, BTagCalibrationReader &readerM, BTagCalibrationReader &readerT, TH2D *h2_qTagEff_NumL, TH2D *h2_qTagEff_NumM, TH2D *h2_qTagEff_NumT, TH2D *h2_qTagEff_Denom, double eta, double pt, double csv, int jetFlavor, int btagsys, bool isCTagLMT, bool isCTagMT, bool isCTagT){
  double pData = 1.0;
  double sfL = 1.0;
  double sfM = 1.0;
  double sfT = 1.0;
  double effL = 1.0;
  double effM = 1.0;
  double effT = 1.0;
  if(isCTagT){//tagged 
    sfT  = getCTagSF(readerT, eta, pt, csv, jetFlavor, btagsys);
    effT = getCTagEff(h2_qTagEff_NumT, h2_qTagEff_Denom, eta, pt);
    pData = sfT*effT;
  }
  else if(isCTagMT){
    effM = getCTagEff(h2_qTagEff_NumM, h2_qTagEff_Denom, eta, pt);
    effT = getCTagEff(h2_qTagEff_NumT, h2_qTagEff_Denom, eta, pt);
    sfM  = getCTagSF(readerM, eta, pt, csv, jetFlavor, btagsys); //csv value of the jet will be same for MT ?
    sfT  = getCTagSF(readerT, eta, pt, csv, jetFlavor, btagsys); 
    pData  = (sfM*effM - sfT*effT);
  }
  else if(isCTagLMT){
    effL = getCTagEff(h2_qTagEff_NumL, h2_qTagEff_Denom, eta, pt);
    effM = getCTagEff(h2_qTagEff_NumM, h2_qTagEff_Denom, eta, pt);
    effT = getCTagEff(h2_qTagEff_NumT, h2_qTagEff_Denom, eta, pt);
    sfL  = getCTagSF(readerL, eta, pt, csv, jetFlavor, btagsys); 
    sfM  = getCTagSF(readerM, eta, pt, csv, jetFlavor, btagsys); 
    sfT  = getCTagSF(readerT, eta, pt, csv, jetFlavor, btagsys); 
    pData = (sfL*effL - sfM*effM)*(sfM*effM - sfT*effT);
  }
  return pData;
}

//https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration#Additional_scripts
//https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation80XReReco/CSVv2_Moriond17_B_H.csv
double CTagSF::getCTagSF(BTagCalibrationReader &reader, double eta, double pt, double csv, double jetFlavor, int cTagSys){
  double sf     =1.0;
  double sfUp   =1.0;
  double sfDown =1.0;
  if(abs(jetFlavor) ==5){
    sf     = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, eta, pt,csv);
    sfUp   = reader.eval_auto_bounds("up",      BTagEntry::FLAV_B, eta, pt, csv);
    sfDown = reader.eval_auto_bounds("down",    BTagEntry::FLAV_B, eta, pt, csv);
  }
  else if(abs(jetFlavor)==4){
    sf     = reader.eval_auto_bounds("central", BTagEntry::FLAV_C, eta, pt,csv);
    sfUp   = reader.eval_auto_bounds("up",      BTagEntry::FLAV_C, eta, pt, csv);
    sfDown = reader.eval_auto_bounds("down",    BTagEntry::FLAV_C, eta, pt, csv);
  }
  else{
    sf     = reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta, pt,csv);
    sfUp   = reader.eval_auto_bounds("up",      BTagEntry::FLAV_UDSG, eta, pt, csv);
    sfDown = reader.eval_auto_bounds("down",    BTagEntry::FLAV_UDSG, eta, pt, csv);
  }
  double scalefactor = 1.0;
  if(cTagSys == kNo)   scalefactor = sf; 
  if(cTagSys == kUp)   scalefactor = sfUp;
  if(cTagSys == kDown) scalefactor = sfDown;
  return scalefactor;
}

//https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
double CTagSF::getCTagEff(TH2D *h2_CTagEff_Num, TH2D *h2_CTagEff_Denom, double eta, double pt){
  double eff = 0.0;
  double bin_num = h2_CTagEff_Num->FindBin(pt, double(eta));
  double bin_denom = h2_CTagEff_Denom->FindBin(pt, double(eta));
  double num = h2_CTagEff_Num->GetBinContent(bin_num); 
  double denom = h2_CTagEff_Denom->GetBinContent(bin_denom); 
  eff = num/denom;
  return eff;  
}
