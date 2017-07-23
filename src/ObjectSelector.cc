#include "interface/ObjectSelector.hh"
#include <iostream>
#include <iomanip>

ClassImp(ObjectSelector)

using namespace std;
void ObjectSelector::preSelectElectrons(vector<int> * e_i, const vector<MyElectron> & vE , MyVertex & vertex, bool isPFlow){
  
  for(unsigned int i=0;i<vE.size();i++){

    const MyElectron * e = &vE[i];
    ///double scEta    = e->electronSCEta;
    double eEta     = TMath::Abs(e->p4.eta());
    double eEt      = TMath::Abs(e->p4.Et());
    double ePt      = TMath::Abs(e->p4.pt());
    ///double eRelIso  = e->RelIso;
    double d0       = fabs(e->D0);
    std::map<std::string, float>idWPs = e->eidWPs;
    float eid = (isPFlow) ? idWPs["eidTightMC"] : idWPs["simpleEleId70cIso"];
    bool passId = (int(eid) & 0x1);
    bool isNotFromConversion = ((int(eid) >> 2) & 0x1);
    ///bool isEcalDriven = (e->isEcalDriven > 0 || e->isPFlow > 0);
    
    double zvertex   = vertex.XYZ.z();
    double zelectron = e->vertex.z();
    double dz = fabs(zvertex - zelectron);
    
    ///if( defaultSelection_  && ( fabs(scEta)> 1.4442 && fabs(scEta)<1.5660)  ){ continue; }
    if( defaultSelection_ ){ continue; }
    if( defaultSelection_  && eEt < E_ET_MIN_ ){ continue; }
    if( !defaultSelection_ && ePt < E_PT_MIN_ ){ continue; }

    if(dz > ZMAX_)continue;

    ///if(passId && isNotFromConversion && isEcalDriven  && d0 < E_D0_MAX_ && eEta < E_ETA_MAX_ && eRelIso < E_RELISO_MAX_){ e_i->push_back(i);}
    if(passId && isNotFromConversion && d0 < E_D0_MAX_ && eEta < E_ETA_MAX_){ e_i->push_back(i);}

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

    bool isGlobalMuon = (m->type & (1<<1));
    bool isPFMuon = (m->type & (1<<5));
    //bool passId = (m->GlobalMuonPromptTight > 0);
    bool passId = (isGlobalMuon && isPFMuon && m->nMuonHits > 0 
		   && m->nPixelHits > 0 && m->nMatchedStations > 1 
		   && m->nTrackerLayers > 5 && m->normChi2 < 10); 

    double zvertex   = vertex.XYZ.z();
    double zmuon     = m->vertex.z();
    double dz =  fabs(zvertex-zmuon);
    
    if(dz > ZMAX_)continue;
    if(passId && mD0 < M_D0_MAX_&&  mPt > M_PT_MIN_ && mEta < M_ETA_MAX_ && mRelIso < M_RELISO_MAX_ ){ m_i->push_back(i);}
    
  }
  
}

void ObjectSelector::preSelectMuonsNoIso(vector<int> * m_i, const vector<MyMuon> & vM , MyVertex & vertex, bool isPFlow){
  
  for( int i=0;i< (int) vM.size();i++){
    
    const MyMuon * m = &vM[i];
    
    double mEta     = TMath::Abs(m->p4.eta());
    double mPt      = TMath::Abs(m->p4.pt());
    double mD0      = fabs(m->D0);
    
    bool isGlobalMuon = (m->type & (1<<1));
    bool isPFMuon = (m->type & (1<<5));
    
    bool passId = (isGlobalMuon && isPFMuon && m->nMuonHits > 0
                   && m->nPixelHits > 0 && m->nMatchedStations > 1
                   && m->nTrackerLayers > 5 && m->normChi2 < 10);
    
    double zvertex   = vertex.XYZ.z();
    double zmuon     = m->vertex.z();
    double dz =  fabs(zvertex-zmuon);
    
    if(dz > ZMAX_)continue;
    
    if(passId && mD0 < M_D0_MAX_&&  mPt > M_PT_MIN_ && mEta < M_ETA_MAX_ ){ m_i->push_back(i);}
    
  }
  
}

void ObjectSelector::preSelectJets( string jetAlgo, vector<int> * j_i, const vector<MyJet> & vJ, int jes, int jer){
  
  for(unsigned int i=0;i<vJ.size();i++){
    
    const MyJet *jet = &vJ[i];
    double jetEta     = TMath::Abs(jet->p4.eta());
    //double jetPt      = TMath::Abs(jet->p4.pt());
    double jetPt      = jetPtWithJESJER(vJ[i], jes, jer); 
    ///double pujetid    = int(jet->puIDMVALoose);

    ///if(jetPt > JET_PT_MIN_ && jetEta < JET_ETA_MAX_ && pujetid==1.0  )
    if(jetPt > JET_PT_MIN_ && jetEta < JET_ETA_MAX_){
      //    if(jetPt > JET_PT_MIN_ && jetEta < JET_ETA_MAX_ ){ 
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
    
    double mEta     = TMath::Abs(m->p4.eta());
    double mPt      = TMath::Abs(m->p4.pt());
    //double mRelIso  = (isPFlow) ? m->RelIso : m->UserPFRelIso;
    double mRelIso  = m->pfRelIso;
    //see if this muon is glogal
    static const unsigned int GlobalMuon     =  1<<1;
    bool isGlobal = (m->type & GlobalMuon);
    
    if(! isGlobal ) continue;
    if( mEta<LOOSE_M_ETA_MAX_  && mPt> LOOSE_M_PT_MIN_ && mRelIso < LOOSE_M_RELISO_MAX_ ){ looseVeto = true; }
    
  }
  
  return looseVeto;
    
}

bool ObjectSelector::looseElectronVeto(int selectedElectron, const vector<MyElectron> & vE, bool isPFlow){
  bool looseVeto(false);
  for( int i=0;i< (int) vE.size();i++ ){
    if(i==selectedElectron){ continue; }
    
    const MyElectron * e = &vE[i];
    
    double eEta      = TMath::Abs(e->p4.eta());
    double eEt       = TMath::Abs(e->p4.Et());
    ///double eRelIso   = e->RelIso;
    std::map<std::string, float>idWPs = e->eidWPs;
    float eid = (isPFlow) ? idWPs["eidLooseMC"] : idWPs["simpleEleId90cIso"];
    bool id = (int(eid) & 0x1);
    double minDR2mu  = 0.4;
        
    ///if(id && eEt > LOOSE_E_ET_MIN_ && eEta < LOOSE_E_ETA_MAX_ && eRelIso < LOOSE_E_RELISO_MAX_ && minDR2mu >0.1){ looseVeto = true;}
    if(id && eEt > LOOSE_E_ET_MIN_ && eEta < LOOSE_E_ETA_MAX_ && minDR2mu >0.1){ looseVeto = true;}
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

