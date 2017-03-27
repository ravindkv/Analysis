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

using namespace std;

class KineFit : public ObjectSelector
{
public :
  KineFit() : ObjectSelector()
  {
    DRMIN_JET = 0.5;
    DRMIN_ELE = 0.5;
    METCUT_   = 50.0;

    TAUPRONGS_ = 0;
    
    std::vector<float> MCPUDist, DataPUDist;
    float dataDist_2011A[25] = {
      0.00592951, 0.0255428, 0.0591468, 0.097016, 0.126287, 0.138848, 0.134117, 0.11691, 0.0937398, 0.0700927, 0.0493627, 0.0329741, 0.0209976, 0.0127917, 0.00747402, 0.00419649, 0.00226774, 0.00118102, 0.000593481, 0.000288109, 0.000135272, 6.14976e-05, 2.71017e-05, 1.15906e-05, 4.81577e-06};
    
    float mc_Dist[25] = {0.112605, 0.0643133, 0.0695342, 0.0698206, 0.06943, 0.0689683, 0.0676831, 0.065551, 0.0628637, 0.0585737, 0.053444, 0.0473381, 0.0410107, 0.0346135, 0.0281938, 0.0224807, 0.0175619, 0.0132723, 0.00987218, 0.00720604, 0.00513216, 0.00355089, 0.00243336, 0.00161019, 0.00107677}; // wjet

    for(int i = 0; i < 25; i++){
      DataPUDist.push_back(dataDist_2011A[i]);
      MCPUDist.push_back(mc_Dist[i]);
    }
    
    LumiWeights_ = reweight::LumiReWeighting(MCPUDist, DataPUDist);
    
  };
  ~KineFit() { };
  
  void SelectMuonJet(TString url,  string myKey="PFlow", string tauIsoLabel = "Loose", bool isData = true, TString histoFiles_="MuonJet_data", TString text_ = "output.txt");
  
  void DefineHistos();
  void processEvents();

private :
  int    TAUPRONGS_;
  double DRMIN_JET, DRMIN_ELE, METCUT_;
  Reader *evR;
  map<string, TH1D*> histos_;
  map<string, TH2D*> histos2_;
  
  reweight::LumiReWeighting LumiWeights_;

};

/// TTbar  ///////
void KineFit::SelectMuonJet(TString url,  string myKey, string tauIsoLabel, bool isData, TString histoFiles_, TString text_){
  
  ofstream myfile1;
  myfile1.open (text_+"_debug.txt");

  string eAlgo("PFlow"), mAlgo("PFlow"), jAlgo("PFlow"), tAlgo("PFlow");
  if(myKey.find("PFlow") != string::npos)jAlgo = "PFlow", tAlgo = "PFlow";
  else if(myKey.find("HPS") != string::npos)jAlgo = "JetsAK5PF", tAlgo = "TausHpsPFTau";

  std::map<string, double> xss;
  xss["WJETS"] = 31314.0;
  xss["TTBAR"] = 165.4;
  xss["ZJETS"] = 3048.0;
  xss["QCD"]   = 84679.3;

  double Lumi = 1000.0;
  
  
  evR = new Reader();
  
  TFile *f = TFile::Open(url);
  if(f==0) return ;
  if(f->IsZombie()) { f->Close(); return; }

  int nEntries = evR->AssignEventTreeFrom(f);
  cout << nEntries << " events are available" << endl;
  if( nEntries == 0) {return; }

  //get initial number of events

  TH1F* inputcf = (TH1F*)(f->Get("allEventsFilter/totalEvents"))->Clone("inputcf");
  double initialEvents = inputcf->GetBinContent(1);

  myfile1<<"Input file : "<<url<<endl;
  myfile1<<"Available input sample : "<<initialEvents<<endl;

  //create the output file ////////////////////////////////

  TString Filename_ = histoFiles_ +"_"+myKey+"_"+tauIsoLabel+".root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel( 9 );
  //////////////////////////////////////////////////

  double sampleWeight(1.0);

  //define histograms
  DefineHistos();
  histos_["w_mt"] = new TH1D("w_mt", "mt (mu, met)", 25, 0., 200.);
  histos_["w_mt_gen"] = new TH1D("w_mt_gen","mt (genmu, genmet)",25, 0., 200.);
  histos_["w_mt_kinfit"] = new TH1D("w_mt_kinfit","mt (kinfitmu, kinfitmet)",25, 0., 200.);
  MyEvent *ev;
  int nTriggEvent = 0, nSelEvents = 0, matchjetcount= 0, threepairjet = 0;
  double nVerticesFailCount = 0.0;
  for(int i=0; i<nEntries; ++i){
    Long64_t ientry = evR->LoadTree(i);
    if (ientry < 0) break;

    ev = evR->GetNewEventFromList(i);
    if(ev==0) continue;

    // apply PU re-weighting

    double evtWeight = 1.0;
    if(!isData){
      //get sample information
      if(i < 1){
        string sampleName = ev->sampleInfo.sampleName;
        sampleWeight = xss[sampleName] * Lumi / initialEvents;
        myfile1<<"Scale factor for lumi "<<Lumi<<" pb is "<< sampleWeight<<endl;
      }
      evtWeight *= sampleWeight; // upto this only sigma*lumi weight is applied
      
      //vector<double>puweights = ev->sampleInfo.puWeights;  // this line and below this line applies the pile up corrections
      
      //if(puweights.size() > 0) evtWeight *= puweights[0];
    }
        
    //trigger
    bool passTrig = false;
    vector<string> trig = ev->hlt;
    for(size_t it = 0; it < trig.size(); it++){
      if(trig[it].find("Mu") != string::npos) passTrig = true;
    }
    if(!passTrig){
      //      cout << "not satisfying trigger" << endl;
      continue;
    }
    nTriggEvent++;
        
    //get objects //
    vector<MyVertex> Vertices = ev->PrimaryVtxs;
    
    if(Vertices.size() <= 0){
      nVerticesFailCount+=evtWeight;
      cout<<" no vertexes , exit"<<endl;
      continue;
    }
    
    vector<MyMuon> pfMuons = evR->getMuons(ev, mAlgo);
    vector<MyElectron> pfElectrons = evR->getElectrons(ev, eAlgo);
    vector<MyJet> pfJets = evR->getJets(ev, jAlgo);

    // preselect objects //
    vector<int> m_init; m_init.clear();
    preSelectMuons(&m_init, pfMuons, Vertices[0]);
    vector<int> e_init; e_init.clear();
    preSelectElectrons(&e_init, pfElectrons, Vertices[0]);
    vector<int> j_init; j_init.clear();
    preSelectJets(jAlgo, &j_init, pfJets);
    
    // clean objects //
    vector<int> e_final; e_final.clear();
    ElectronCleaning( pfElectrons, pfMuons, &e_init, &e_final, &m_init, DRMIN_ELE);
    vector<int> j_final; j_final.clear();
    JetCleaning(pfJets, pfMuons, pfElectrons, &j_init, &j_final, &m_init, &e_final, DRMIN_JET, false);
    
    // Muon properties plot before selection 

    int nMuon = m_init.size();
    histos_["Muon_mult"]->Fill(nMuon, evtWeight);

    vector<MyLorentzVector> HadP; HadP.clear();
    vector<MyLorentzVector> HadQ; HadQ.clear();
    vector<MyLorentzVector> HadB; HadB.clear();
    vector<MyLorentzVector> Lepton; Lepton.clear();
    vector<MyLorentzVector> Met; Met.clear();

    vector<MyKineFitParticle> allKineFitParticles = ev->KineFitParticles;
    for(size_t imk=0; imk < allKineFitParticles.size(); imk++){
      if(allKineFitParticles[imk].partName.find("PartonsHadP") != string::npos ){
        HadP.push_back(allKineFitParticles[imk].p4);
      }
      if(allKineFitParticles[imk].partName.find("PartonsHadQ") != string::npos ){
        HadQ.push_back(allKineFitParticles[imk].p4);
      }
      if(allKineFitParticles[imk].partName.find("Leptons") != string::npos ){
        Lepton.push_back(allKineFitParticles[imk].p4);
      }
      if(allKineFitParticles[imk].partName.find("Neutrinos") != string::npos ){
        Met.push_back(allKineFitParticles[imk].p4);
      }
      if(allKineFitParticles[imk].partName.find("PartonsHadB") !=string::npos ){
	HadB.push_back(allKineFitParticles[imk].p4);
      }
      //myfile1 << "Particle Type:    "<< allKineFitParticles[imk].partName << "  Pt:  "<< allKineFitParticles[imk].p4.pt()<<endl; 
      
    }
    
    /*    
    // Applying condition to get the parton level information
    vector<MyLorentzVector> quarks; quarks.clear();
    vector<MyLorentzVector> muons; muons.clear();
    vector<MyLorentzVector> bquark; bquark.clear();
    vector<MyLorentzVector> cquark; cquark.clear();
    vector<MyLorentzVector> squark; squark.clear();
    if(!(ev->isData)){
      vector<MyMCParticle>allMCParticles = ev->mcParticles;
      MyMET mcMET = ev->mcMET;
      histos_["mcMET"]->Fill(mcMET.p4.pt(),evtWeight);

      for(size_t imc=0; imc < allMCParticles.size(); ++imc){
        if(allMCParticles[imc].status != 3) continue;
        if(abs(allMCParticles[imc].pid) == 13 && allMCParticles[imc].mother.size() > 0 && abs(allMCParticles[imc].mother[0])==24 ){
          myfile1 << "Particle type for muons:  " << allMCParticles[imc].pid << "       Mother id for muons:   " << allMCParticles[imc].mother[0] << endl;
          muons.push_back(allMCParticles[imc].p4Gen);
        }
	if(abs(allMCParticles[imc].pid) == 5 && allMCParticles[imc].mother.size() > 0 && abs(allMCParticles[imc].mother[0])==6 ){
	  myfile1 << "Particle type for b-quark:  " << allMCParticles[imc].pid << "       Mother id for b-quark:   " << allMCParticles[imc].mother[0] << endl;
	  bquark.push_back(allMCParticles[imc].p4Gen);
	}
	if(abs(allMCParticles[imc].pid) == 4 && allMCParticles[imc].mother.size() > 0 && abs(allMCParticles[imc].mother[0])==37 ){
          myfile1 << "Particle type for c-quark:  " << allMCParticles[imc].pid << "       Mother id for c-quark:   " << allMCParticles[imc].mother[0] << endl;
          cquark.push_back(allMCParticles[imc].p4Gen);
        }
	if(abs(allMCParticles[imc].pid) == 3 && allMCParticles[imc].mother.size() > 0 && abs(allMCParticles[imc].mother[0])==37 ){
          myfile1 << "Particle type for s-quark:  " << allMCParticles[imc].pid << "       Mother id for s-quark:   " << allMCParticles[imc].mother[0] << endl;
          squark.push_back(allMCParticles[imc].p4Gen);
        }


      }
      // here calculating mt using generated muon and met
      double   genleptonPt(0), gendeltaPhi(0);
      double genmetPt = mcMET.p4.pt();
      if(genmetPt > 30){
	
	//	int m_j = m_init[0];
	genleptonPt = TMath::Abs(muons[0].pt());
	gendeltaPhi = ROOT::Math::VectorUtil::DeltaPhi(muons[0], mcMET.p4);
	double genmt = sqrt (  2*genleptonPt*genmetPt*(1 - cos(gendeltaPhi) ) ) ;
	histos_["w_mt_gen"]->Fill(genmt, evtWeight);
      }
      
      for(size_t imc=0; imc < allMCParticles.size(); ++imc){
	if(allMCParticles[imc].status != 3) continue;
	if(abs(allMCParticles[imc].pid) <= 4 && allMCParticles[imc].mother.size() > 0 && abs(allMCParticles[imc].mother[0])==37 ){
	  myfile1 << "Particle type:  " << allMCParticles[imc].pid << "       Mother id:   " << allMCParticles[imc].mother[0] << endl;
	  quarks.push_back(allMCParticles[imc].p4Gen);
	}
      }
    }
    
    */

    //Apply Lepton selection//////////////////////////////
    
    int nLepton = m_init.size();  // this condition proof that only muon + jet events
    if(nLepton != 1)continue;
    if( looseMuonVeto( m_init[0],pfMuons) ) continue;
    if( looseElectronVeto(-1,pfElectrons) ) continue;
    int count_muon = m_init.size();
    int m_i = m_init[0];
    double delta = DeltaR(pfMuons[m_i].p4 , Lepton[0]);
    if(delta > 0.2) continue;
    histos_["Muon_mult_final"]->Fill(count_muon, evtWeight);
    
    // Fill histogram after trigger and one offline isolated muon and applied 2nd lepton veto
    
    histos_["pt_muon"]->Fill(pfMuons[m_i].p4.pt(),evtWeight);
    histos_["eta_muon"]->Fill(pfMuons[m_i].p4.eta(),evtWeight);
    histos_["phi_muon"]->Fill(pfMuons[m_i].p4.phi(),evtWeight);
    /*
    // generator level pt information
    histos_["pt_bjet_gen"]->Fill(bquark[0].pt(), evtWeight);
    histos_["pt_cjet_gen"]->Fill(cquark[0].pt(), evtWeight);
    histos_["pt_sjet_gen"]->Fill(squark[0].pt(), evtWeight);
   */
    // vertex just after one lepton selection
    
    double pri_vtxs = Vertices.size();
    histos_["pri_vtxs"]->Fill(pri_vtxs, evtWeight);
    
    // Jet properties before applying jet criteria
    int prenJet = j_final.size();
    histos_["pre_jet_mult"]->Fill(prenJet, evtWeight);
    

    // Apply Jet Selection
    int nJet = j_final.size();
    if(nJet < 4)continue;  // this condition implies events should contains more than 4 jets
    histos_["Jet_mult"]->Fill(nJet, evtWeight);
        
    int threejet = 0;
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
//       histos_["pt_alljet"]->Fill(pfJets[ind_jet].p4.pt(), evtWeight);   
//       histos_["eta_alljet"]->Fill(pfJets[ind_jet].p4.eta(), evtWeight);
//       histos_["phi_alljet"]->Fill(pfJets[ind_jet].p4.phi(), evtWeight); 
      
      if(pfJets[ind_jet].p4.pt() > 30 )threejet++;
    }
    if(threejet < 3) continue;
    
    // Met distribution   
    
    vector<MyMET> mets = ev->mets;
    if(mets.size() <= 0)continue;
    MyMET met;
    for(size_t imet = 0; imet < mets.size(); imet++){
      if(mets[imet].metName.find("PFlow") != string::npos)met = mets[imet];
    }
    histos_["pre_Met"]->Fill(met.p4.pt(),evtWeight);

    //count number of b-Jets with combinedSecondaryVertex Medium working point after selecting the events contains atleast 4 jets
    
    int Post_nBTaggedJets = 0;
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      histos_["btag_value_pre"]->Fill(pfJets[ind_jet].bDiscriminator["combinedSecondaryVertexBJetTags"], evtWeight);
      if(pfJets[ind_jet].bDiscriminator["combinedSecondaryVertexBJetTags"] > 0.679)Post_nBTaggedJets++;
    }

    double   leptonPt(0), deltaPhi(0);
    double metPt = met.p4.pt();
    if(metPt < 30) continue;  // Missing transverse energy cut 30 GeV

    histos_["Final_Met"]->Fill(metPt,evtWeight);
    int m_j = m_init[0];
    leptonPt = TMath::Abs(pfMuons[m_j].p4.pt());
    deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(pfMuons[m_j].p4, met.p4);
    double mt = sqrt (  2*leptonPt*metPt*(1 - cos(deltaPhi) ) ) ;
    histos_["w_mt"]->Fill(mt, evtWeight);
    
    // Here putting the criteria that event should contain at least one b-tag jet  
    int count = 0; 
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      if(pfJets[ind_jet].bDiscriminator["combinedSecondaryVertexBJetTags"] > 0.679) count++;
      histos_["btag_value_final"]->Fill(pfJets[ind_jet].bDiscriminator["combinedSecondaryVertexBJetTags"], evtWeight);
      
    }
    if(count <= 0) continue;
    histos_["btag_jet_mult"]->Fill(count,evtWeight);

    nSelEvents++;  // this is the counter for counting final number of events  

    // matching  KinFit Jets to PAT jets
    //    double smalldelta1 = 0.1; 
    vector<int> IndexForBtag; IndexForBtag.clear();
    vector<double> BtagForKinFit; BtagForKinFit.clear();
    bool first_flag = false;
    bool sec_flag = false;
    bool third_flag = false;
    for(size_t ijet1 = 0; ijet1 < j_final.size(); ijet1++){
      double tmp_delta1 = DeltaR(pfJets[j_final[ijet1]].p4, HadP[0]);
      if(tmp_delta1 < 0.2){
	first_flag = true;
	double btag_HadP = pfJets[j_final[ijet1]].bDiscriminator["combinedSecondaryVertexBJetTags"];
	//	myfile1 <<"indexForHadP:  "<< ijet1 <<"   btag_HadP::   " << btag_HadP << endl;
	IndexForBtag.push_back(ijet1);
	BtagForKinFit.push_back(btag_HadP);
      }
    }
    if(first_flag){
      for(size_t ijet2 = 0; ijet2 < j_final.size(); ijet2++){
	double tmp_delta2 = DeltaR(pfJets[j_final[ijet2]].p4, HadQ[0]);
	if(tmp_delta2 < 0.2){
	  sec_flag = true;
	  double btag_HadQ = pfJets[j_final[ijet2]].bDiscriminator["combinedSecondaryVertexBJetTags"];
	  //	  myfile1 <<"indexForHadQ:  "<< ijet2 <<"   btag_HadQ::   " << btag_HadQ << endl;
	  IndexForBtag.push_back(ijet2);
	  BtagForKinFit.push_back(btag_HadQ);
	}
      }
    }
    if(first_flag && sec_flag){
      for(size_t ijet3 = 0; ijet3 < j_final.size(); ijet3++){
        double tmp_delta3 = DeltaR(pfJets[j_final[ijet3]].p4, HadB[0]);
        if(tmp_delta3 < 0.2){
	  third_flag = true;
          double btag_HadB = pfJets[j_final[ijet3]].bDiscriminator["combinedSecondaryVertexBJetTags"];
	  //          myfile1 <<"indexForHadB:  "<< ijet3 <<"   btag_HadB::   " << btag_HadB << endl;
	  IndexForBtag.push_back(ijet3);
	  BtagForKinFit.push_back(btag_HadB);
        }
      }
    }
    
    myfile1 << endl;
    myfile1 << endl;

    if(first_flag && sec_flag && third_flag && IndexForBtag.size()==3 && BtagForKinFit.size()==3){
      myfile1 << "found three correct jets " << endl;
      threepairjet++; 
      myfile1 << "BtagForKinFit index are:   "<<BtagForKinFit[0]<<"  "<<BtagForKinFit[1]<<"  "<<BtagForKinFit[2]<< endl;
      //      myfile1 << "IndexForBtag index are:   "<<IndexForBtag[0]  <<IndexForBtag[1]  <<IndexForBtag[2]<< endl;
      if (BtagForKinFit[0] >= BtagForKinFit[1] && BtagForKinFit[0] >= BtagForKinFit[2]){
	MyLorentzVector mjj_kinfit_btag = HadQ[0] + HadB[0] ;
	histos_["mjj_kinfit_btag"]->Fill(mjj_kinfit_btag.M(), evtWeight);
	myfile1 << "largest number is::   " << BtagForKinFit[0] << endl;
      }
      else if (BtagForKinFit[1] >= BtagForKinFit[0] && BtagForKinFit[1] >= BtagForKinFit[2]){
	MyLorentzVector mjj_kinfit_btag = HadP[0] + HadB[0] ;
        histos_["mjj_kinfit_btag"]->Fill(mjj_kinfit_btag.M(), evtWeight);
	myfile1 << "largest number is::   " << BtagForKinFit[1] << endl;
      }
      else if (BtagForKinFit[2] >= BtagForKinFit[0] && BtagForKinFit[2] >= BtagForKinFit[1]){
	MyLorentzVector mjj_kinfit_btag = HadP[0] + HadQ[0] ;
	histos_["mjj_kinfit_btag"]->Fill(mjj_kinfit_btag.M(), evtWeight);
	myfile1 << "largest number is::   " << BtagForKinFit[2] << endl;
      }
    }



    //    myfile1 << "total number of selected events " << nSelEvents <<endl;
    if((nSelEvents%50)==0){
      cout << "Events number:    " << nSelEvents << endl;
    }
    
    //loop over jets 
    
    vector<int> j_final_nob; j_final_nob.clear();
    double nvtxs = Vertices.size();
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      //b-tag jet selection
      if(pfJets[ind_jet].bDiscriminator["combinedSecondaryVertexBJetTags"] > 0.679){
	histos_["pt_bjet"]->Fill(pfJets[ind_jet].p4.pt(), evtWeight);
	histos_["eta_bjet"]->Fill(pfJets[ind_jet].p4.eta(), evtWeight);
	continue;
      }
      j_final_nob.push_back(ind_jet);
      histos_["pt_alljet"]->Fill(pfJets[ind_jet].p4.pt(), evtWeight);
      histos_["eta_alljet"]->Fill(pfJets[ind_jet].p4.eta(), evtWeight);
      histos_["phi_alljet"]->Fill(pfJets[ind_jet].p4.phi(), evtWeight);
      histos_["nvtx_alljet"]->Fill(nvtxs, evtWeight);
    }
    

    /*
    if(quarks.size() >= 2){
      MyLorentzVector mjj_quark = quarks[0] + quarks[1];
      histos_["mjj_quark"]->Fill(mjj_quark.M(), evtWeight);
    }
    */
    //    double delta = DeltaR(pfMuons[m_i].p4 , Lepton[0]);
    //    myfile1 << "delta::     " << delta << endl;

    double   kinfitleptonPt(0), kinfitdeltaPhi(0);
    double kinfitmetPt = Met[0].pt();
    if( kinfitmetPt > 30){

      kinfitleptonPt = TMath::Abs(Lepton[0].pt());
      kinfitdeltaPhi = ROOT::Math::VectorUtil::DeltaPhi(Lepton[0], Met[0]);
      double kinfitmt = sqrt (  2*kinfitleptonPt*kinfitmetPt*(1 - cos(kinfitdeltaPhi) ) ) ;
      histos_["w_mt_kinfit"]->Fill(kinfitmt, evtWeight);
    }

    if( HadP.size() > 0 && HadQ.size() > 0 && HadB.size() > 0){
      MyLorentzVector mjj_kinfit = HadP[0] + HadQ[0] ;
      histos_["mjj_kinfit"]->Fill(mjj_kinfit.M(), evtWeight);
      MyLorentzVector trijet_kinfit = HadP[0] + HadQ[0] + HadB[0];
      histos_["trijet_kinfit"]->Fill(trijet_kinfit.M(), evtWeight);
    }
    
    // matching with jet in cone 0.3
    if(HadP.size() > 0 && HadQ.size() > 0){
      for(size_t index1 = 0; index1 < j_final.size(); index1++){
	  double q1dr = DeltaR(pfJets[j_final[index1]].p4, HadP[0]); 
	  if(q1dr < 0.3){ 
	    for(size_t index2 = 0; index2 < j_final.size() ; index2++){
	      if(index2 != index1){
		double q2dr = DeltaR(pfJets[j_final[index2]].p4, HadQ[0]); 
		if(q2dr < 0.3){
		  matchjetcount++;
		  MyLorentzVector mjj = HadP[0] + HadQ[0] ;
		  histos_["mjj"]->Fill(mjj.M(), evtWeight);
		} 
	      } // if(index2
	    } //index2 for loop
	  }  
      } // index1 for loop                                                                                                                                 
    } // first if loop      
    
    if(j_final_nob.size() >= 2){
      int first_index = j_final_nob[0];
      int sec_index = j_final_nob[1];  
      histos_["pt_cjet"]->Fill(pfJets[first_index].p4.pt(), evtWeight);
      histos_["pt_sjet"]->Fill(pfJets[sec_index].p4.pt(), evtWeight);
      MyLorentzVector diJet = pfJets[first_index].p4 + pfJets[sec_index].p4;
      histos_["diJet_Mass"]->Fill(diJet.M(), evtWeight); 
      //  jet parton matching
      /*
	if(quarks.size() >= 2){      
	for(size_t index1 = 0; index1 < j_final_nob.size(); index1++){
	for(size_t qq1 = 0; qq1 < quarks.size(); qq1++){
	double q1dr = DeltaR(pfJets[j_final_nob[index1]].p4, quarks[qq1]);
	if(q1dr < 0.4){
	for(size_t index2 = index1+1; index2 < j_final_nob.size(); index2++){
	for(size_t qq2 = qq1+1; qq2 < quarks.size(); qq2++){
	double q2dr = DeltaR(pfJets[j_final_nob[index2]].p4, quarks[qq2]);
	if(q2dr < 0.4){
	MyLorentzVector mjj = pfJets[j_final_nob[index1]].p4 + pfJets[j_final_nob[index2]].p4;
	histos_["mjj"]->Fill(mjj.M(), evtWeight);
	}
	} // qq2 for loop
	} //index2 for loop
	}
	} // qq1 for loop
	} // index1 for loop
	} // first if loop
      */
    }
  }

  myfile1 << "Number of times HadP and HadQ matches with jets" << matchjetcount << endl;
  myfile1 << "total number of selected events " << nSelEvents <<endl; 
  myfile1 << "No of times three pair jet matched:    " << threepairjet << endl;
  outFile_->Write();
  outFile_->Close();
  myfile1.close();
  
}

// new function 

void KineFit::DefineHistos()
{

  histos_["Muon_mult"] = new TH1D("Muon_mult","Muon_multiplicity",10,0,10);
  // Event Yield

  histos_["btag_value_pre"] = new TH1D("btag_value_pre","btag discriminator for CSV",40,0,1.5);
  histos_["btag_value_final"] = new TH1D("btag_value_final","btag discriminator for CSV",40,0,1.5);
  histos_["pre_btag_jet_mult"] = new TH1D("pre_btag_jet_mult","pre_btag_jet_mult",10,0,10);
  histos_["pre_jet_mult"] = new TH1D("pre_jet_mult","Jet multiplicity before selection ",10,0,10);
  histos_["pre_Met"] = new TH1D("pre_Met","pre_Met",40,0,200);
  histos_["Muon_mult_final"] = new TH1D("Muon_mult_final","Muon_multiplicity",10,0,10);
  // now all the histograms are after applying various cut
  
  histos_["Final_Met"] = new TH1D("Final_Met","Final MET after selection",50,0,200);
  histos_["mcMET"] = new TH1D("mcMET","generator level met",50,0,200);
  histos_["btag_jet_mult"] = new TH1D("btag_jet_mult","b-tagging multiplicity",10,0,10);
  histos_["pt_muon"] = new TH1D("pt_muon", "pt muon",200,0,200);
  histos_["eta_muon"] = new TH1D("eta_muon","eta muon",50,-2.5,2.5);
  histos_["phi_muon"] = new TH1D("phi_muon","phi muon",20,-3.5,3.5);
  histos_["pt_alljet"] = new TH1D("pt_alljet", "pt alljet", 200,0,200);
  histos_["eta_alljet"] = new TH1D("eta_alljet", "eta alljet", 50, -2.5, 2.5);

  histos_["pt_bjet_gen"] = new TH1D("pt_bjet_gen", "pt_bjet_gen",40,0,200);
  histos_["pt_cjet_gen"] = new TH1D("pt_cjet_gen", "pt_cjet_gen",40,0,200); 
  histos_["pt_sjet_gen"] = new TH1D("pt_sjet_gen", "pt_sjet_gen",40,0,200);
  histos_["pt_bjet"] = new TH1D("pt_bjet", "pt_bjet",40,0,200);
  histos_["pt_cjet"] = new TH1D("pt_cjet","pt_cjet",40,0,200);
  histos_["pt_sjet"] = new TH1D("pt_sjet","pt_sjet",40,0,200);
  histos_["eta_bjet"] = new TH1D("eta_bjet","eta_bjet",50,-2.5,2.5);
  histos_["phi_alljet"] = new TH1D("phi_alljet", "phi alljet", 14, -3.5, 3.5);

  histos_["Jet_mult"] = new TH1D("Jet_mult","Jet multiplicity after selection",10,0,10);
  
  histos_["diJet_Mass"] = new TH1D("diJet_Mass","DiJet Mass",50,0,200);
  histos_["mjj"] = new TH1D("mjj","Dijet_Mass_with_quark_matching",50,0,200);
  histos_["mjj_quark"]= new TH1D("mjj_quark","Dijet_Mass_with_quark_only",25,0,200);
  
  histos_["mjj_kinfit_btag"] = new TH1D("mjj_kinfit_btag","mjj_kinfit_btag",50,0,200);
  histos_["mjj_kinfit"]= new TH1D("mjj_kinfit","mjj_kinfit",50,0,200);
  histos_["trijet_kinfit"]= new TH1D("trijet_kinfit","trijet_kinfit",50,0,200);
  histos_["pri_vtxs"] = new TH1D("pri_vtxs", "pri_vtxs", 30, 0, 30.);
  histos_["nvtx_alljet"] = new TH1D("nvtx_alljet", "nvtx_alljet", 30, 0, 30.);

  
}

void KineFit::processEvents(){
  //  SelectMuonJet("/tmp/gkole/Muon_all/Muon_all.root", "PFlow", "Loose",  true, "MuonJet_data", "MuonJet_data");
  //  SelectMuonJet("/tmp/gkole/wjet-su11/wjet_1_9.root", "PFlow", "Loose",  false, "MuonJet_wjet", "MuonJet_wjet");
  //  SelectMuonJet("/tmp/gkole/wjet-su11/wjet_all.root","PFlow", "Loose",  false, "MuonJet_wjet", "MuonJet_wjet");
  //  SelectMuonJet("/tmp/gkole/zjet-su11/zjet_all.root", "PFlow", "Loose", false, "MuonJet_zjet", "MuonJet_zjet");
  //  SelectMuonJet("/tmp/gkole/ttbar-su11/ttbar_all.root", "PFlow", "Loose", false, "MuonJet_ttbar", "MuonJet_ttbar");
  //  SelectMuonJet("/tmp/gkole/qcd-su11/qcd_all.root", "PFlow","Loose", false, "MuonJet_qcd", "MuonJet_qcd");
  //  SelectMuonJet("/tmp/gkole/signal/sample/chargedHiggs_all_new.root", "PFlow","Loose", false, "MuonJet_hp", "MuonJet_hp");  
  //SelectMuonJet("/tmp/gkole/signal/sample/mc_charged_higgs_all.root", "PFlow","Loose", false, "MuonJet_hp_test", "MuonJet_hp_test");
  //  SelectMuonJet("/tmp/gkole/mc_charged_higgs.root", "PFlow","Loose", false, "MuonJet_hp_test", "MuonJet_hp_test");
  //  SelectMuonJet("/tmp/gkole/mc_charged_higgs_csbar_wmu.root", "PFlow","Loose", false, "MuonJet_hp_test", "MuonJet_hp_test");
  //  SelectMuonJet("/tmp/gkole/8TeV/ttjet/mc_ttjet_all.root", "PFlow","Loose", false, "TTbar_jet", "TTbar_jet");
  //  SelectMuonJet("/tmp/gkole/8TeV/ttjet_cons234/mc_ttjet_cons234_all.root", "PFlow","Loose", false, "TTbar_jet", "TTbar_jet");
  SelectMuonJet("/tmp/gkole/8TeV/52x/ttjet_all_nobtag.root", "PFlow","Loose", false, "TTbar_jet", "TTbar_jet");
}
