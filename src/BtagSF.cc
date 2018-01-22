#include "interface/BtagSF.hh"

ClassImp(BtagSF)

Bool_t BtagSF::isbtagged(BTagCalibrationReader &reader, TH2D *h2_BTagEff_Num, TH2D *h2_BTagEff_Denom, Float_t eta, Float_t pt, Float_t csv, Int_t jetflavor, Bool_t isdata, UInt_t btagsys) 
{ 
  //----------------------------------------------
  //for DATA     
  //----------------------------------------------
  //cout<<"isbtagged called "<<endl;
  //double csvL = 0.5426;
  double csvM = 0.8484;
  //double csvT = 0.9535;
  double csvOP = csvM;
  Bool_t btagged = kFALSE; 
  if(isdata){ 
    if(csv>csvOP) btagged = kTRUE; 
    else          btagged = kFALSE; 
    return btagged; 
  } 
  //----------------------------------------------
  //for b-quark jet     
  //----------------------------------------------
  Double_t SFb = 0.0; 
  Double_t eff_b = 0.0; 
  Double_t promoteProb_btag=0; // ~probability to promote to tagged 
  Double_t demoteProb_btag=0; // ~probability to demote from tagged 
  if(fabs(jetflavor) == 5) {                // real b-jet 
    SFb = getBTagSFb(reader, eta, pt, csv, btagsys);
    eff_b = getBTagEff(h2_BTagEff_Num, h2_BTagEff_Denom, pt, eta);
    if(SFb < 1) demoteProb_btag = fabs(1.0 - SFb); 
    else promoteProb_btag = fabs(SFb - 1.0)/((1/eff_b) - 1.0); 
    if(csv > csvOP){
      btagged = kTRUE; // if tagged 
      if(demoteProb_btag > 0 && randm->Uniform() < demoteProb_btag) btagged = kFALSE;  // demote it to untagged  
      else                  					    btagged = kTRUE; // leave it tagged 	
    }	      
    else{ 
      btagged = kFALSE; 
      if(promoteProb_btag > 0 && randm->Uniform() < promoteProb_btag) btagged = kTRUE;  // promote it to tagged 
      else                                                            btagged = kFALSE; // leave it untagged 
    } 
    return btagged; 
  } 
  
  //----------------------------------------------
  //for c-quark jet     
  //----------------------------------------------
  Double_t SFl = 0, eff_l = 0;
  Double_t promoteProb_mistag=0; // ~probability to promote to tagged 
  Double_t demoteProb_mistag=0; // ~probability to demote from tagged 
  if(fabs(jetflavor) == 4) {
    SFl = getBTagSFc(reader, eta, pt, csv, btagsys);
    eff_l = getBTagEff(h2_BTagEff_Num, h2_BTagEff_Denom, pt, eta);
    if(SFl > 1) promoteProb_mistag = fabs(SFl - 1.0)/((1/eff_l) - 1.0); 
    else demoteProb_mistag = SFl; 
    if(csv > csvOP) {         // if tagged 
      btagged = kTRUE; 
      if(demoteProb_mistag > 0 && randm->Uniform() > demoteProb_mistag) btagged = kFALSE; // demote it to untagged 
      else                                                              btagged = kTRUE;  // leave it tagged 
    }
    else {                    // not tagged 
      btagged = kFALSE; 
      if(promoteProb_mistag > 0 && randm->Uniform() < promoteProb_mistag) btagged = kTRUE;  // promote it to tagged 
      else                                                                btagged = kFALSE; // leave it untagged 
    } 
  } 
  //----------------------------------------------
  //for other(light) quarks and gluon jet     
  //----------------------------------------------
  else{
    SFl = getBTagSFl(reader, eta, pt, csv, btagsys);
    
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //
    //One doubt: what, if SFl = 0 ??
    //
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    eff_l = getBTagEff(h2_BTagEff_Num, h2_BTagEff_Denom, pt, eta);
    ///if(SFl!=0) cout<<"SFl = "<<SFl<<endl;
    if(SFl > 1) promoteProb_mistag = fabs(SFl - 1.0)/((1/eff_l) - 1.0); 
    else demoteProb_mistag = SFl; 
    if(csv > csvOP) {         // if tagged 
      btagged = kTRUE; 
      if(demoteProb_mistag > 0 && randm->Uniform() > demoteProb_mistag) btagged = kFALSE; // demote it to untagged 
      else                                                              btagged = kTRUE;  // leave it tagged 
    }
    else {                    // not tagged 
      btagged = kFALSE; 
      if(promoteProb_mistag > 0 && randm->Uniform() < promoteProb_mistag) btagged = kTRUE;  // promote it to tagged 
      else                                                                btagged = kFALSE; // leave it untagged 
    } 
  }
  return btagged; 
} 

//https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration#Additional_scripts
//https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation80XReReco/CSVv2_Moriond17_B_H.csv
//scale factors for b-quark
Double_t BtagSF::getBTagSFb(BTagCalibrationReader &reader, Float_t eta, Float_t pt, Float_t csv, UInt_t btagsys){
  double SFb     = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, eta, pt,csv);
  double SFbUp   = reader.eval_auto_bounds("up", BTagEntry::FLAV_B, eta, pt, csv);
  double SFbDown = reader.eval_auto_bounds("down", BTagEntry::FLAV_B, eta, pt, csv);
  Double_t scalefactor = 1.0;
  if(btagsys == kNo)   scalefactor = SFb; 
  if(btagsys == kUp)   scalefactor = SFbUp;
  if(btagsys == kDown) scalefactor = SFbDown;
  return scalefactor;
}
//scale factors for c-quark
Double_t BtagSF::getBTagSFc(BTagCalibrationReader &reader, Float_t eta, Float_t pt, Float_t csv, UInt_t btagsys){
  double SFc     = reader.eval_auto_bounds("central", BTagEntry::FLAV_C, eta, pt,csv);
  double SFcUp   = reader.eval_auto_bounds("up", BTagEntry::FLAV_C, eta, pt, csv);
  double SFcDown = reader.eval_auto_bounds("down", BTagEntry::FLAV_C, eta, pt, csv);
  Double_t scalefactor = 1.0;
  if(btagsys == kNo)   scalefactor = SFc; 
  if(btagsys == kUp)   scalefactor = SFcUp;
  if(btagsys == kDown) scalefactor = SFcDown;
  return scalefactor;
}

//scale factors for u, d, s-quark, gluon
Double_t BtagSF::getBTagSFl(BTagCalibrationReader &reader, Float_t eta, Float_t pt, Float_t csv, UInt_t btagsys){
  double SFl     = reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta, pt,csv);
  double SFlUp   = reader.eval_auto_bounds("up", BTagEntry::FLAV_UDSG, eta, pt, csv);
  double SFlDown = reader.eval_auto_bounds("down", BTagEntry::FLAV_UDSG, eta, pt, csv);
  Double_t scalefactor = 1.0;
  if(btagsys == kNo)   scalefactor = SFl; 
  if(btagsys == kUp)   scalefactor = SFlUp;
  if(btagsys == kDown) scalefactor = SFlDown;
  return scalefactor;
}

//https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
Double_t BtagSF::getBTagEff(TH2D *h2_BTagEff_Num, TH2D *h2_BTagEff_Denom, Float_t pt, Float_t eta){
  double eff = 0.0;
  double bin_num = h2_BTagEff_Num->FindBin(pt, double(eta));
  double bin_denom = h2_BTagEff_Denom->FindBin(pt, double(eta));
  double num = h2_BTagEff_Num->GetBinContent(bin_num); 
  double denom = h2_BTagEff_Denom->GetBinContent(bin_denom); 
  eff = num/denom;
  return eff;  
}
