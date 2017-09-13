#include "interface/ObjectSelector.hh"
#include <iostream>
#include <iomanip>

ClassImp(ObjectSelector)

using namespace std;
void ObjectSelector::preSelectElectrons(vector<int> * e_i, const vector<MyElectron> & vE , MyVertex & vertex, bool isPFlow){
 
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria
  for(unsigned int i=0;i<vE.size();i++){
    const MyElectron * e   = &vE[i];
    double sigmaIetaIeta   = e->sigmaIetaIeta;
    double dEtaInSeed      = e->dEtaInSeed; 
    double dPhiIn          = e->dPhiIn;     
    double hadOverEm       = e->hadOverEm;  
    double iEminusiP       = abs(e->iEminusiP);  
    double nInnerHits      = e->nInnerHits; 
    double eRelIso  	   = e->relCombPFIsoEA; 

    double eEta     	   = TMath::Abs(e->p4.eta());
    double ePt     	   = TMath::Abs(e->p4.pt());
    double d0      	   = fabs(e->D0);
    double zvertex   	   = vertex.XYZ.z();
    double zelectron 	   = e->vertex.z();
    double dz 		   = fabs(zvertex - zelectron);

    std::map<std::string, float>idWPs = e->eidWPs;
    float eid = (isPFlow) ? idWPs["eidLoose"] : idWPs["cutBasedElectronID-Spring15-25ns-V1-standalone-veto"];
    //float eid = idWPs["eidLoose"];
    bool passId = (int(eid) & 0x1);

    bool isPassConVeto = e->isPassConVeto;
    //barrel cuts ( |eta supercluster| <= 1.479) 
    if(abs(eEta) <= 1.479 
       && ePt  			>30 
       && sigmaIetaIeta 	< 0.00998 
       && abs(dEtaInSeed) 	< 0.00308
       && abs(dPhiIn) 		< 0.0816 
       && hadOverEm 		< 0.0414 
       && abs(iEminusiP) 	< 0.0129 
       && eRelIso 		< 0.0588
       && nInnerHits       	<= 1
       && isPassConVeto    
       && d0  			< 0.05 
       && dz 			< 0.1 
       && passId    
      )e_i->push_back(i);

    //endcap cuts ( |eta supercluster| > 1.479) |
    if(abs(eEta) > 1.479 
       && ePt  			>30 
       && sigmaIetaIeta 	< 0.0292 
       && abs(dEtaInSeed) 	< 0.00605 
       && abs(dPhiIn) 		< 0.0394
       && hadOverEm 		< 0.0641
       && abs(iEminusiP) 	< 0.0129
       && eRelIso 		< 0.0571
       && nInnerHits       	<= 1
       && isPassConVeto    
       && d0  			< 0.1 
       && dz 			< 0.3 
       && passId    
      )e_i->push_back(i);
    //To print all eidWPs
    /*
    map <std::string, float> :: iterator itr;
    for(itr = idWPs.begin(); itr != idWPs.end(); ++itr){
      cout<<'\t'<< itr->first<<'\t'<<itr->second<<'\n';
    }
    //All possible names: 	
        cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose	0
	cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium	0
	cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight	0
	cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto	0
	cutBasedElectronID-Spring15-25ns-V1-standalone-loose	0
	cutBasedElectronID-Spring15-25ns-V1-standalone-medium	0
	cutBasedElectronID-Spring15-25ns-V1-standalone-tight	0
	cutBasedElectronID-Spring15-25ns-V1-standalone-veto	0
	cutBasedElectronID-Spring15-50ns-V2-standalone-loose	0
	cutBasedElectronID-Spring15-50ns-V2-standalone-medium	0
	cutBasedElectronID-Spring15-50ns-V2-standalone-tight	0
	cutBasedElectronID-Spring15-50ns-V2-standalone-veto	0
	eidLoose	4
	eidRobustHighEnergy	4
	eidRobustLoose	4
	eidRobustTight	4
	eidTight	4
	heepElectronID-HEEPV60	0
	mvaEleID-Spring15-25ns-Trig-V1-wp80	0
	mvaEleID-Spring15-25ns-Trig-V1-wp90	0
	mvaEleID-Spring15-25ns-nonTrig-V1-wp80	0
	mvaEleID-Spring15-25ns-nonTrig-V1-wp90	0
	mvaEleID-Spring15-25ns-nonTrig-V1-wpLoose	0
	mvaEleID-Spring15-50ns-Trig-V1-wp80	0
	mvaEleID-Spring15-50ns-Trig-V1-wp90	0
     */
  }
}

void ObjectSelector::preSelectMuons(vector<int> * m_i, const vector<MyMuon> & vM , MyVertex & vertex, bool isPFlow){
  
  for( int i=0;i< (int) vM.size();i++){
    
    const MyMuon * m = &vM[i];
    double mEta     = TMath::Abs(m->p4.eta());
    double mPt      = TMath::Abs(m->p4.pt());
    double mD0      = fabs(m->D0);
    //double mRelIso  = (isPFlow) ? m->RelIso : m->UserPFRelIso;
    double mRelIso  = m->pfRelIso;

    //bool isGlobalMuon = (m->type & (1<<1)); // 1 = 01, 1<<1 = 10 = 2
    //bool isPFMuon = (m->type & (1<<5)); //1 = 000001, 1<<5 = 100000 = 32
    bool isGlobalMuon = m->isGlobalMuon; 
    bool isPFMuon = m->isPFMuon; 
    //bool passId = (m->GlobalMuonPromptTight > 0);
    bool passId = (isGlobalMuon && isPFMuon && m->nMuonHits >=1 
		   && m->nPixelHits >= 1 && m->nMatchedStations >=2 
		   && m->nTrackerLayers >= 6 && m->normChi2 < 10); 

    double zvertex   = vertex.XYZ.z();
    double zmuon     = m->vertex.z();
    double dz =  fabs(zvertex-zmuon);
    
    if(passId && mPt > M_PT_MIN_ && mEta < M_ETA_MAX_ && 
		    mD0 < M_D0_MAX_&& dz < ZMAX_ && mRelIso < M_RELISO_MAX_ ){ 
      m_i->push_back(i);
    }
  }
  
}

void ObjectSelector::preSelectMuonsNoIso(vector<int> * m_i, const vector<MyMuon> & vM , MyVertex & vertex, bool isPFlow){
  
  for( int i=0;i< (int) vM.size();i++){
    
    const MyMuon * m = &vM[i];
    double mEta     = TMath::Abs(m->p4.eta());
    double mPt      = TMath::Abs(m->p4.pt());
    double mD0      = fabs(m->D0);

    bool isGlobalMuon = m->isGlobalMuon; 
    bool isPFMuon = m->isPFMuon; 
    bool passId = (isGlobalMuon && isPFMuon && m->nMuonHits >=1 
		   && m->nPixelHits >= 1 && m->nMatchedStations >=2 
		   && m->nTrackerLayers >= 6 && m->normChi2 < 10); 

    double zvertex   = vertex.XYZ.z();
    double zmuon     = m->vertex.z();
    double dz =  fabs(zvertex-zmuon);
    
    if(passId && mPt > M_PT_MIN_ && mEta < M_ETA_MAX_ && mD0 < M_D0_MAX_&& dz < ZMAX_ ){ 
      m_i->push_back(i);
    }
  }
  
}


void ObjectSelector::preSelectJets( string jetAlgo, vector<int> * j_i, const vector<MyJet> & vJ, int jes, int jer, double sigmaJER){
 
  for(unsigned int i=0;i<vJ.size();i++){
    const MyJet *jet = &vJ[i];
    double jetEta     = TMath::Abs(jet->p4.eta());
    //double jetPt      = TMath::Abs(jet->p4.pt());
    double jetPt      = jetPtWithJESJER(vJ[i], jes, jer, sigmaJER); 
    double neutralHadEnFrac = jet->neutralHadronEnergyFraction;
    double neutralEmEnFrac = jet->neutralEmEnergyFraction;
    double chargedHadEnFrac = jet->chargedHadronEnergyFraction;

    ///double pujetid    = int(jet->puIDMVALoose);
    ///if(jetPt > JET_PT_MIN_ && jetEta < JET_ETA_MAX_ && pujetid==1.0  )
    if(jetPt > JET_PT_MIN_ 
      && jetEta < JET_ETA_MAX_
      && neutralHadEnFrac < 0.9
      && neutralEmEnFrac  < 0.9
      && chargedHadEnFrac > 0
      ){
      j_i->push_back(i);
    }
  }
}

//https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
bool ObjectSelector::isMediumMuon(const MyMuon * m, bool isPFlow){
  
  bool isMedium(false);
  bool goodGlob = m->isGlobalMuon && 
	  m->normChi2 && 
	  m->chi2LocalPosition < 12 && 
	  m->trkKink < 20; 
  bool isLooseMuon = m->isPFMuon && 
          (m->isGlobalMuon || m->isTrackerMuon);
  isMedium =  isLooseMuon &&  
	    m->validFraction > 0.8 && 
	    m->segmentCompatibility >(goodGlob ? 0.303 : 0.451); 
  return isMedium; 
}


bool ObjectSelector::looseMuonVeto( int selectedMuon, const vector<MyMuon> & vM, bool isPFlow){

  bool looseVeto(false);
  
  for(int i=0;i< (int)vM.size();i++){
    
    if( i==selectedMuon ){continue;}
    
    const MyMuon * m = &vM[i];
    
    bool isGlobalMuon = m->isGlobalMuon; 
    double mEta     = TMath::Abs(m->p4.eta());
    double mPt      = TMath::Abs(m->p4.pt());
    double mRelIso  = m->pfRelIso;
    
    if(! isGlobalMuon) continue;
    if( mEta<LOOSE_M_ETA_MAX_  && mPt> LOOSE_M_PT_MIN_ && mRelIso < LOOSE_M_RELISO_MAX_ ){ looseVeto = true; }
  }
  return looseVeto;
    
}

bool ObjectSelector::loose2ndMuonVeto( int firstMuon, int secondMuon, const vector<MyMuon> & vM, bool isPFlow){

  bool looseVeto(false);
  
  for(int i=0;i< (int)vM.size();i++){
    
    if( i==firstMuon ){continue;}
    if( i==secondMuon   ){continue;}
    
    const MyMuon * m = &vM[i];
    double mEta     = TMath::Abs(m->p4.eta());
    double mPt      = TMath::Abs(m->p4.pt());
    double mRelIso  = m->pfRelIso;
    
    bool isGlobalMuon = m->isGlobalMuon; 
    if(! isGlobalMuon ) continue;
    if( mEta<LOOSE_M_ETA_MAX_  && mPt> LOOSE_M_PT_MIN_ && mRelIso < LOOSE_M_RELISO_MAX_ ){ looseVeto = true; }
  }
  return looseVeto;
    
}

bool ObjectSelector::looseElectronVeto(unsigned long selectedElectron, const vector<MyElectron> & vE, MyVertex & vertex, bool isPFlow){

  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria
  bool looseVeto(false);
  for(unsigned long i=0;i<vE.size();i++){
    const MyElectron * e   = &vE[i];
    double sigmaIetaIeta   = e->sigmaIetaIeta;
    double dEtaInSeed      = e->dEtaInSeed; 
    double dPhiIn          = e->dPhiIn;     
    double hadOverEm       = e->hadOverEm;  
    double iEminusiP       = abs(e->iEminusiP);  
    double nInnerHits      = e->nInnerHits; 
    double eRelIso  	   = e->relCombPFIsoEA; 

    double eEta     	   = TMath::Abs(e->p4.eta());
    double ePt     	   = TMath::Abs(e->p4.pt());
    double d0      	   = fabs(e->D0);
    double zvertex   	   = vertex.XYZ.z();
    double zelectron 	   = e->vertex.z();
    double dz 		   = fabs(zvertex - zelectron);

    std::map<std::string, float>idWPs = e->eidWPs;
    //float eid = (isPFlow) ? idWPs["eidLoose"] : idWPs["cutBasedElectronID-Spring15-25ns-V1-standalone-veto"];
    float eid = idWPs["eidLoose"];
    bool passId = (int(eid) & 0x1);

    bool isPassConVeto = e->isPassConVeto;
    //barrel cuts ( |eta supercluster| <= 1.479) 
    if( i==selectedElectron) continue; 
    if(abs(eEta) <= 1.479 
       && ePt  			> 15 
       && sigmaIetaIeta 	< 0.011 
       && abs(dEtaInSeed) 	< 0.00477
       && abs(dPhiIn) 		< 0.222
       && hadOverEm 		< 0.298
       && abs(iEminusiP) 	< 0.241 
       && eRelIso 		< 0.0994
       && nInnerHits       	<= 1
       && isPassConVeto    
       && d0  			< 0.05 
       && dz 			< 0.1 
       && passId    
      )looseVeto = true;

    //endcap cuts ( |eta supercluster| > 1.479) |
    if(abs(eEta) > 1.479 
       && ePt  			> 15 
       && sigmaIetaIeta 	< 0.0314 
       && abs(dEtaInSeed) 	< 0.00868
       && abs(dPhiIn) 		< 0.213
       && hadOverEm 		< 0.101
       && abs(iEminusiP) 	< 0.14
       && eRelIso 		< 0.107
       && nInnerHits       	<=1
       && isPassConVeto    
       && d0  			< 0.1 
       && dz 			< 0.3 
       && passId    
      )looseVeto = true;
  }
  return looseVeto;
}

void ObjectSelector::ElectronCleaning( const vector<MyElectron> & vE, const vector<MyMuon> & vM, vector<int> * e_old, vector<int> * e_new, vector<int> * mu, double DR ){
  
  for(size_t i = 0; i < e_old->size(); i++){
    int iele = (*e_old)[i];
    double delR2Mu = 5.0;
    for(size_t j = 0; j < mu->size(); j++){
      int imu = (*mu)[j];
      double delR = DeltaR(vE[iele].p4, vM[imu].p4);
      if(delR < delR2Mu)delR2Mu = delR;
    }

    if(delR2Mu > DR) e_new->push_back(iele);
  }
}

void ObjectSelector::JetCleaning(const vector<MyJet> & vJ, const vector<MyMuon> & vM, const vector<MyElectron> & vE,vector<int> * j_old, vector<int> * j_new, vector<int> * mu, vector<int> * el, double DR){

  for(size_t i = 0; i < j_old->size(); i++){
    int ijet = (*j_old)[i];

    double delR2Mu = 5.0, delR2Ele = 5.0;
    
    for(size_t j = 0; j < mu->size(); j++){
      int imu = (*mu)[j];
      double delR = DeltaR(vJ[ijet].p4, vM[imu].p4);
      if(delR < delR2Mu)delR2Mu = delR;
    }

    for(size_t k = 0; k < el->size(); k++){
      int iele = (*el)[k];
      double delR = DeltaR(vJ[ijet].p4, vE[iele].p4);
      if(delR < delR2Ele)delR2Ele = delR;
    }
    if(delR2Mu > DR && delR2Ele > DR )
    {
        j_new->push_back(ijet);
    }
    }
}

double ObjectSelector::DeltaR(MyLorentzVector aV, MyLorentzVector bV){
  
  double deta = TMath::Abs(aV.eta() - bV.eta());
  double dphi = TMath::Abs(aV.phi() - bV.phi());
  if(dphi > M_PI) dphi = 2*M_PI - dphi;

  double delR = sqrt(deta*deta + dphi*dphi);

  return delR;
}

