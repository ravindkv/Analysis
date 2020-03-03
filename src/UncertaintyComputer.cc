#include "interface/UncertaintyComputer.hh"
#include <iostream>
#include <iomanip>
#include "TRandom3.h"
#include <stack/Roch/RoccoR.cc>

ClassImp(UncertaintyComputer)

using namespace std;
RoccoR  rc("stack/Roch/rcdata.2016.v3");
double UncertaintyComputer::muPtWithRochCorr(const MyMuon *mu, bool isData, double u1, double u2, int s, int m){
  double charge = mu->charge;
  double pt 	= mu->p4.pt();
  double eta 	= mu->p4.eta();
  double phi 	= mu->p4.phi();
  int nl 	= mu->nTrackerLayers;
  double dataSF = rc.kScaleDT(charge, pt, eta, phi, s, m); 
  double mcSF 	= rc.kScaleAndSmearMC(charge, pt, eta, phi, nl, u1, u2, s, m); 
  double SF = 1.0; 
  if(isData)SF = dataSF;
  else SF = mcSF;
  ///cout<<pt<<"\t"<<SF<<"\t"<<SF*pt<<"\t"<<charge<<"\t"<<eta<<"\t"<<phi<<"\t"<<nl<<"\t"<<u1<<"\t"<<u2<<"\t"<<s<<"\t"<<m<<endl;
  return SF*pt;
}

double UncertaintyComputer::getJERSF(double eta, int jer){
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
  //New SF with stat + sys unc
  //https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Spring16_25nsV10a_MC/Spring16_25nsV10a_MC_SF_AK4PFchs.txt
  std::vector<double> JEREtaMap;
  std::vector<double> JERSF;
  std::vector<double> JERSFUp;
  std::vector<double> JERSFDown;
  //eta values
  JEREtaMap.push_back(0.0); 
  JEREtaMap.push_back(0.5); 
  JEREtaMap.push_back(0.8); 
  JEREtaMap.push_back(1.1); 
  JEREtaMap.push_back(1.3); 
  JEREtaMap.push_back(1.7); 
  JEREtaMap.push_back(1.9); 
  JEREtaMap.push_back(2.1); 
  JEREtaMap.push_back(2.3); 
  JEREtaMap.push_back(2.5); 
  JEREtaMap.push_back(2.8); 
  JEREtaMap.push_back(3.0); 
  JEREtaMap.push_back(3.2); 
  JEREtaMap.push_back(4.7);
  //SF values
  JERSF.push_back(1.109); JERSFDown.push_back(1.044); JERSFUp.push_back(1.174);
  JERSF.push_back(1.138); JERSFDown.push_back(1.072); JERSFUp.push_back(1.204);
  JERSF.push_back(1.114); JERSFDown.push_back(1.050); JERSFUp.push_back(1.178);
  JERSF.push_back(1.123); JERSFDown.push_back(1.022); JERSFUp.push_back(1.224);
  JERSF.push_back(1.084); JERSFDown.push_back(0.985); JERSFUp.push_back(1.183);
  JERSF.push_back(1.082); JERSFDown.push_back(0.973); JERSFUp.push_back(1.191);
  JERSF.push_back(1.140); JERSFDown.push_back(1.020); JERSFUp.push_back(1.260);
  JERSF.push_back(1.067); JERSFDown.push_back(0.953); JERSFUp.push_back(1.181);
  JERSF.push_back(1.177); JERSFDown.push_back(0.967); JERSFUp.push_back(1.387);
  JERSF.push_back(1.364); JERSFDown.push_back(1.203); JERSFUp.push_back(1.525);
  JERSF.push_back(1.857); JERSFDown.push_back(1.654); JERSFUp.push_back(2.060);
  JERSF.push_back(1.328); JERSFDown.push_back(1.203); JERSFUp.push_back(1.453);
  JERSF.push_back(1.160); JERSFDown.push_back(1.013); JERSFUp.push_back(1.307);
  double SF = 1.0;
  for(size_t i = 0; i < 13; i++){
    if(TMath::Abs(eta) >= JEREtaMap[i] && TMath::Abs(eta) < JEREtaMap[i+1]){
      if(jer == 0)SF = JERSF[i];
      else if (jer == 1) SF = JERSFUp[i];
      else if(jer == -1) SF = JERSFDown[i];  
    }
  }
  return SF;
}
//https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections
double UncertaintyComputer::metWithJESJER(const vector<MyJet> & vJ, vector<int> *j, MyMET MET, int jes, int jer, bool isData, const int & seed) 
{ 
  double metX = MET.p4.px(); 
  double metY = MET.p4.py(); 
  if(!isData){
    std::mt19937 m_random_generator;
    m_random_generator = std::mt19937(seed);
    for(size_t i = 0; i < j->size(); i++){ 
      int j_ind = j->at(i); 
      double gen_pt = vJ[j_ind].Genp4.pt(); 
      double jet_pt = vJ[j_ind].p4.pt();
      double sigmaJER = vJ[j_ind].resolution;
      //apply JER uncert, scaling
      double delR = DeltaR(vJ[j_ind].Genp4, vJ[j_ind].p4);
      double rCone = 0.4;
      //https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
      MyLorentzVector rawJet = vJ[j_ind].p4; 
      metX += rawJet.px(); 
      metY += rawJet.py(); 
      double SF = getJERSF(vJ[j_ind].p4.eta(), jer); 
      double jerF = 1.0;
      std::normal_distribution<> d(0, sigmaJER);
      double N0sigma = d(m_random_generator);
      //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
      if(delR<rCone/2 && abs(jet_pt -gen_pt)<3*sigmaJER*jet_pt ){
        jerF = max(0.0, 1.0 + (SF - 1)*(jet_pt - gen_pt)/jet_pt); 
      }
      else{
        jerF   = 1+ N0sigma* std::sqrt(std::max(SF*SF -1, 0.0));
      }
      metX -= jerF*rawJet.px(); 
      metY -= jerF*rawJet.py();
      //apply JES scaling
      double jesF = vJ[j_ind].JECUncertainty*double(jes);
      metX -= jesF*rawJet.px(); 
      metY -= jesF*rawJet.py();
    } 
  }  
  return sqrt(metX*metX + metY*metY); 
}

//bottom mistagging, by event re-weighting 
double UncertaintyComputer::getBTagPmcSys(TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, MyJet jet){
  double csv =jet.bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
  double pMC = 1.0; 
  pMC = btsf->getBTagPmc(h2_qTagEff_Num, h2_qTagEff_Denom, jet.p4.eta(), jet.p4.pt(), csv);
  return pMC;
}
double UncertaintyComputer::getBTagPdataSys(BTagCalibrationReader &reader, TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, MyJet jet, int scale){
  double pData = 1.0;
  double csv =jet.bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
  double eta = jet.p4.eta();
  double pt = jet.p4.pt();
  int flavor = abs(jet.partonFlavour);
  if(scale == 0) pData = btsf->getBTagPdata(reader, h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt, csv, flavor ,0);
  else if(scale == 1) pData = btsf->getBTagPdata(reader, h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt, csv, flavor ,1);
  else if(scale == -1) pData = btsf->getBTagPdata(reader, h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt, csv, flavor ,-1);
  return pData;
}

//charm mistagging, for inclusive categories, by event re-weighting 
double UncertaintyComputer::getIncCTagPmcSys(TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, MyJet jet, bool isCTag){
  double pMC = 1.0; 
  pMC = ctsf->getIncCTagPmc(h2_qTagEff_Num, h2_qTagEff_Denom, jet.p4.eta(), jet.p4.pt(), isCTag);
  return pMC;
}
double UncertaintyComputer::getIncCTagPdataSys(BTagCalibrationReader &reader, TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, MyJet jet, bool isCTag, int scale){
  double pData = 1.0;
  //double csv =jet.bDiscriminator["pfCombinedCvsLJetTags"]; 
  double csv =abs(jet.bDiscriminator["pfCombinedCvsBJetTags"]); //SF is calculated as a func of pT
  double eta = jet.p4.eta();
  double pt = jet.p4.pt();
  int flavor = abs(jet.partonFlavour);
  if(scale == 0) pData = ctsf->getIncCTagPdata(reader, h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt, csv, isCTag, flavor ,0);
  else if(scale == 1) pData = ctsf->getIncCTagPdata(reader, h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt, csv, isCTag, flavor ,1);
  else if(scale == -1) pData = ctsf->getIncCTagPdata(reader, h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt, csv, isCTag, flavor ,-1);
  return pData;
}

double UncertaintyComputer::DeltaR(MyLorentzVector aV, MyLorentzVector bV){
  double deta = TMath::Abs(aV.eta() - bV.eta());
  double dphi = TMath::Abs(aV.phi() - bV.phi());
  if(dphi > M_PI) dphi = 2*M_PI - dphi;
  double delR = sqrt(deta*deta + dphi*dphi);
  return delR;
}
