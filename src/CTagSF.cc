#include "interface/CTagSF.hh"

ClassImp(CTagSF)

Bool_t CTagSF::isCtagged(BTagCalibrationReader &reader, TH2D *h2_CTagEff_Num, TH2D *h2_CTagEff_Denom, Float_t csv_OP, Float_t csv, Float_t eta, Float_t pt , Int_t jetflavor, Bool_t isdata, UInt_t cTagSys) 
{ 
  
  /*
  cout<<"-------------------------"<<endl;
  cout<<"csv_OP = "<<csv_OP<<endl;
  cout<<"csv    = "<<csv<<endl;
  cout<<"eta    = "<<eta<<endl;
  cout<<"pt     = "<<pt<<endl;
  cout<<"jetflavor ="<<jetflavor<<endl;
  cout<<"cTagSys= "<<cTagSys<<endl;
  */
  
  //----------------------------------------------
  //for DATA     
  //----------------------------------------------
  Bool_t cTagged = kFALSE; 
  if(isdata){ 
    if(csv > csv_OP) cTagged = kTRUE; 
    else          cTagged = kFALSE; 
    return cTagged; 
  } 
  //----------------------------------------------
  //for b-quark jet     
  //----------------------------------------------
  Double_t SFb = 0.0; 
  Double_t eff_b = 0.0; 
  Double_t promoteProb_btag=0; // ~probability to promote to tagged 
  Double_t demoteProb_btag=0; // ~probability to demote from tagged 
  if(fabs(jetflavor) == 5) {                // real b-jet 
    SFb = getCTagSFb(reader, eta, pt, csv, cTagSys);
    eff_b = getCTagEff(h2_CTagEff_Num, h2_CTagEff_Denom, pt, eta);
    if(SFb < 1) demoteProb_btag = fabs(1.0 - SFb); 
    else promoteProb_btag = fabs(SFb - 1.0)/((1/eff_b) - 1.0); 
    if(csv > csv_OP){
      cTagged = kTRUE; // if tagged 
      if(demoteProb_btag > 0 && randm->Uniform() < demoteProb_btag) cTagged = kFALSE;  // demote it to untagged  
      else                  					    cTagged = kTRUE; // leave it tagged 	
    }	      
    else{ 
      cTagged = kFALSE; 
      if(promoteProb_btag > 0 && randm->Uniform() < promoteProb_btag) cTagged = kTRUE;  // promote it to tagged 
      else                                                            cTagged = kFALSE; // leave it untagged 
    } 
    return cTagged; 
  } 
  
  //----------------------------------------------
  //for c-quark jet     
  //----------------------------------------------
  Double_t SFl = 0, eff_l = 0;
  Double_t promoteProb_mistag=0; // ~probability to promote to tagged 
  Double_t demoteProb_mistag=0; // ~probability to demote from tagged 
  if(fabs(jetflavor) == 4) {
    SFl = getCTagSFc(reader, eta, pt, csv, cTagSys);
    eff_l = getCTagEff(h2_CTagEff_Num, h2_CTagEff_Denom, pt, eta);
    if(SFl > 1) promoteProb_mistag = fabs(SFl - 1.0)/((1/eff_l) - 1.0); 
    else demoteProb_mistag = SFl; 
    if(csv > csv_OP) {         // if tagged 
      cTagged = kTRUE; 
      if(demoteProb_mistag > 0 && randm->Uniform() > demoteProb_mistag) cTagged = kFALSE; // demote it to untagged 
      else                                                              cTagged = kTRUE;  // leave it tagged 
    }
    else {                    // not tagged 
      cTagged = kFALSE; 
      if(promoteProb_mistag > 0 && randm->Uniform() < promoteProb_mistag) cTagged = kTRUE;  // promote it to tagged 
      else                                                                cTagged = kFALSE; // leave it untagged 
    } 
  } 
  //----------------------------------------------
  //for other(light) quarks and gluon jet     
  //----------------------------------------------
  else{
    SFl = getCTagSFl(reader, eta, pt, csv, cTagSys);
    
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //
    //One doubt: what, if SFl = 0 ??
    //
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    eff_l = getCTagEff(h2_CTagEff_Num, h2_CTagEff_Denom, pt, eta);
    ///if(SFl!=0) cout<<"SFl = "<<SFl<<endl;
    if(SFl > 1) promoteProb_mistag = fabs(SFl - 1.0)/((1/eff_l) - 1.0); 
    else demoteProb_mistag = SFl; 
    if(csv > csv_OP) {         // if tagged 
      cTagged = kTRUE; 
      if(demoteProb_mistag > 0 && randm->Uniform() > demoteProb_mistag) cTagged = kFALSE; // demote it to untagged 
      else                                                              cTagged = kTRUE;  // leave it tagged 
    }
    else {                    // not tagged 
      cTagged = kFALSE; 
      if(promoteProb_mistag > 0 && randm->Uniform() < promoteProb_mistag) cTagged = kTRUE;  // promote it to tagged 
      else                                                                cTagged = kFALSE; // leave it untagged 
    } 
  }
  return cTagged; 
} 

//https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration#Additional_scripts
//https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation80XReReco/CSVv2_Moriond17_B_H.csv
//scale factors for b-quark
Double_t CTagSF::getCTagSFb(BTagCalibrationReader &reader, Float_t eta, Float_t pt, Float_t csv, UInt_t cTagSys){
  double SFb     = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, eta, pt,csv);
  double SFbUp   = reader.eval_auto_bounds("up", BTagEntry::FLAV_B, eta, pt, csv);
  double SFbDown = reader.eval_auto_bounds("down", BTagEntry::FLAV_B, eta, pt, csv);
  Double_t scalefactor = 1.0;
  if(cTagSys == kNo)   scalefactor = SFb; 
  if(cTagSys == kUp)   scalefactor = SFbUp;
  if(cTagSys == kDown) scalefactor = SFbDown;
  return scalefactor;
}
//scale factors for c-quark
Double_t CTagSF::getCTagSFc(BTagCalibrationReader &reader, Float_t eta, Float_t pt, Float_t csv, UInt_t cTagSys){
  double SFc     = reader.eval_auto_bounds("central", BTagEntry::FLAV_C, eta, pt,csv);
  double SFcUp   = reader.eval_auto_bounds("up", BTagEntry::FLAV_C, eta, pt, csv);
  double SFcDown = reader.eval_auto_bounds("down", BTagEntry::FLAV_C, eta, pt, csv);
  Double_t scalefactor = 1.0;
  if(cTagSys == kNo)   scalefactor = SFc; 
  if(cTagSys == kUp)   scalefactor = SFcUp;
  if(cTagSys == kDown) scalefactor = SFcDown;
  return scalefactor;
}

//scale factors for u, d, s-quark, gluon
Double_t CTagSF::getCTagSFl(BTagCalibrationReader &reader, Float_t eta, Float_t pt, Float_t csv, UInt_t cTagSys){
  double SFl     = reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta, pt,csv);
  double SFlUp   = reader.eval_auto_bounds("up", BTagEntry::FLAV_UDSG, eta, pt, csv);
  double SFlDown = reader.eval_auto_bounds("down", BTagEntry::FLAV_UDSG, eta, pt, csv);
  Double_t scalefactor = 1.0;
  if(cTagSys == kNo)   scalefactor = SFl; 
  if(cTagSys == kUp)   scalefactor = SFlUp;
  if(cTagSys == kDown) scalefactor = SFlDown;
  return scalefactor;
}

//https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
Double_t CTagSF::getCTagEff(TH2D *h2_CTagEff_Num, TH2D *h2_CTagEff_Denom, Float_t pt, Float_t eta){
  double eff = 0.0;
  double bin_num = h2_CTagEff_Num->FindBin(pt, double(eta));
  double bin_denom = h2_CTagEff_Denom->FindBin(pt, double(eta));
  double num = h2_CTagEff_Num->GetBinContent(bin_num); 
  double denom = h2_CTagEff_Denom->GetBinContent(bin_denom); 
  eff = num/denom;
  //cout<<"eff = "<<eff<<endl;
  return eff;  
}
