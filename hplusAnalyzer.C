
#include "hplusAnalyzer.h"

using namespace std;
void hplusAnalyzer::CutFlowAnalysis(TString url, string myKey, string evtType){
  
  TString outFile("13TeV/outputDir/");
  TString Filename_ = outFile+evtType+"_Anal.root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel( 9 );
  
  TString debug_Filename_ = Filename_+"_debug.txt";
  string debug_file(debug_Filename_);
  outfile_.open(debug_file.c_str());

  CutFlowProcessor(url, myKey, "base", outFile_);
  outfile_.close();
  outFile_->Write(); 
  outFile_->Close(); 
}

//---------------------------------------------------//
//Process the cuts, event by event
//---------------------------------------------------//  
void hplusAnalyzer::CutFlowProcessor(TString url,  string myKey, TString cutflowType, TFile *outFile_){
  
  int input_count = 0;
  bool isPFlow = (myKey.find("PFlow") != string::npos) ? true : false;
  string eAlgo("Electrons"), mAlgo("Muons"), jAlgo("Jets"), metAlgo("METs");
  
  //Uncertainty variations, JES, JER, MET unclustered, bTag
  ///int jes = 0, jer = 0, metuc = 0, bscale = 0;
  int jes = 1, jer = 1; 
  if(cutflowType.Contains("JESPlus"))jes = 1;
  else if (cutflowType.Contains("JESMinus"))jes = -1;
  else if (cutflowType.Contains("JERPlus"))jer = 1;
  else if (cutflowType.Contains("JERMinus"))jer = -1;
  //else if (cutflowType.Contains("METUCPlus"))metuc = 1;
  //else if (cutflowType.Contains("METUCMinus"))metuc = -1;
  //else if (cutflowType.Contains("bTagPlus"))bscale = 1;
  //else if (cutflowType.Contains("bTagMinus"))bscale = -1; 

  evR = new Reader();
  TFile *f = TFile::Open(url);
  if(f==0) return ;
  if(f->IsZombie()) { f->Close(); return; }
 
  //---------------------------------------------------//
  //get initial number of events, from ntuples
  //store initial informations, in a txt file
  //---------------------------------------------------//
  double lumiTotal = 35500;
  int nEntries = evR->AssignEventTreeFrom(f);
  if( nEntries == 0) {return; }
  TH1F* inputcf = (TH1F*)(f->Get("allEventsFilter/totalEvents"))->Clone("inputcf");
  TH1F* intimepu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/intimepu"))->Clone("intimepu");
  TH1F* outoftimepu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/outoftimepu"))->Clone("outoftimepu");
  TH1F* totalpu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/totalpu"))->Clone("totalpu");
  TH1F* trueintimepu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/trueintimepu"))->Clone("trueintimepu");
  TH1F* trueoutoftimepu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/trueoutoftimepu"))->Clone("trueoutoftimepu");
  TH1F* truetotalpu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/truetotalpu"))->Clone("truetotalpu");
  double initialEvents = inputcf->GetBinContent(1);
  cout<<"input file: "<<url<<endl;
  outfile_<<"input file: "<<url<<endl;
  outfile_<<"totalEvents: "<<initialEvents<<endl;
  CreateAnalHistos(cutflowType, outFile_);
  for(int nbin=1; nbin<6001; nbin++){
    fillHistoPU("intimepu", cutflowType, nbin, intimepu->GetBinContent(nbin));  
    fillHistoPU("outoftimepu", cutflowType, nbin, outoftimepu->GetBinContent(nbin));  
    fillHistoPU("totalpu", cutflowType, nbin, totalpu->GetBinContent(nbin));  
    fillHistoPU("trueintimepu", cutflowType, nbin, trueintimepu->GetBinContent(nbin));  
    fillHistoPU("trueoutoftimepu", cutflowType, nbin, trueoutoftimepu->GetBinContent(nbin));  
    fillHistoPU("truetotalpu", cutflowType, nbin, truetotalpu->GetBinContent(nbin));  
  }
  fillHisto("totalEvents", cutflowType, initialEvents, 1);
  MyEvent *ev;
  int nTriggEvent = 0, nSelEvents = 0, matchjetcount= 0, threepairjet = 0;
  double nVerticesFailCount = 0.0;
  ///double corr_jet_pair_cs = 0.0, corr_jet_pair_ud = 0.0;
  double matchedJet_q = 0.0, matchedJet_b = 0.0, matched_quark_eta_pt = 0.0, not_matchedJet_q = 0.0;
  double TotalTopPtWeights = 0, TotalLplusJEvents = 0; 
  double TotalSVEff = 0, TotalEvtsWithSV = 0.;
  
  //---------------------------------------------------//
  //loop over each event, of the ntuple
  //---------------------------------------------------//
  for(int i=0; i<nEntries; ++i){
    Long64_t ientry = evR->LoadTree(i);
    if (ientry < 0) break;
    ev = evR->GetNewEvent(i);
    if(ev==0) continue;
    if(i%1000==0) cout<<"\033[01;32mEvent number = \033[00m"<< i << endl;
 
    //---------------------------------------------------//
    //apply lumi, k factor and pileup weight
    //---------------------------------------------------//
    double evtWeight = 1.0;
    if(!ev->isData){
      string sampleName = ev->sampleInfo.sampleName;
      //k-factor weight (along with lumi weight) 
      if(sampleName.find("WJetsToLNu") != string::npos || sampleName.find("W1JetsToLNu") != string::npos || sampleName.find("W2JetsToLNu") != string::npos || sampleName.find("W3JetsToLNu") != string::npos || sampleName.find("W4JetsToLNu") != string::npos){
	    int hepNUP = ev->sampleInfo.hepNUP;
        double weightK = reweightHEPNUPWJets(hepNUP) * (lumiTotal/1000.0);
        evtWeight *= weightK;  
        fillHisto("hepNUP", cutflowType, hepNUP, evtWeight);
      }
      else if(sampleName.find("DYJetsToLL") != string::npos || sampleName.find("DY1JetsToLL") != string::npos || sampleName.find("DY2JetsToLL") != string::npos || sampleName.find("DY3JetsToLL") != string::npos || sampleName.find("DY4JetsToLL") != string::npos){
        int hepNUP = ev->sampleInfo.hepNUP;
        std::vector<int> hepIDUP = ev->sampleInfo.hepIDUP;
        std::vector<int> hepISTUP = ev->sampleInfo.hepISTUP;
        int countZ = 0;
        for(size_t p=0; p<hepIDUP.size(); p++){
          if(hepIDUP[p]==23 && hepISTUP[p]==2)
            countZ = countZ + 1;
        }
        if(countZ==0) hepNUP = hepNUP+1;
        double weightK = reweightHEPNUPDYJets(hepNUP) * (lumiTotal/1000.0);
        evtWeight *= weightK;  
        fillHisto("hepNUP", cutflowType, hepNUP, evtWeight);
      }
      //lumi weight
      else {
      double sampleWeight(1.0);
      sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
      evtWeight *= sampleWeight; 
      if(i < 1){
        outfile_<<"Sample weight for "<<sampleName<<" is= "<< sampleWeight<<endl;
        }
      }
      //pileup weight 
      vector<double>pu = ev->sampleInfo.truepileup;
      if(pu.size() > 0) {
	  float npu = pu[0];
	  double weightPU = LumiWeights_.weight(npu);
      evtWeight *= weightPU;  
      }
    } 

    //---------------------------------------------------//
    // apply top re-weighting weight
    //---------------------------------------------------//
    double topPtWeights_offline = 1.0;
    if(!ev->isData){
      string sampleName = ev->sampleInfo.sampleName;
      //(we have to update it ==========)
      if(sampleName.find("TTJets") != string::npos || sampleName.find("WH") != string::npos || sampleName.find("HplusM120") != string::npos){
        vector<double>topptweights = ev->sampleInfo.topPtWeights;
        if(topptweights.size() > 0){
          topPtWeights_offline = topptweights[0];
          if(cutflowType.Contains("TopPtPlus"))
            topPtWeights_offline = topPtWeights_offline*topPtWeights_offline;
          else if(cutflowType.Contains("TopPtMinus"))
            topPtWeights_offline = 1.0;
        }
      }
    }
    TotalTopPtWeights += topPtWeights_offline;
    //cout<<"topPtWeights_offline = "<<topPtWeights_offline<<endl;
    TotalLplusJEvents++;
    evtWeight *= topPtWeights_offline; //Multiply to the total weights
    
    //---------------------------------------------------//
    //apply muon triggers
    //---------------------------------------------------//
    double nCutPass = 0.0;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);
    bool passTrig = false;
    vector<string> trig = ev->hlt;
    for(size_t it = 0; it < trig.size(); it++){
      if(trig[it].find("Mu") != string::npos) passTrig = true;
    }
    if(!passTrig){
            //cout << "not satisfying trigger" << endl;
      continue;
    }
    nTriggEvent++;
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);
   
    //---------------------------------------------------//
    //get all objets e.g. leptons, jets, vertices etc.
    //---------------------------------------------------//
    vector<MyVertex> Vertices = ev->PrimaryVtxs;
    if(Vertices.size() <= 0){
      nVerticesFailCount+=evtWeight;
      cout<<" no vertexes , exit"<<endl;
      continue;
    }
    vector<MyMuon> pfMuons = evR->getMuons(ev, mAlgo);
    vector<MyMuon> pfMuons_noiso = evR->getMuons(ev, mAlgo);
    vector<MyElectron> pfElectrons = evR->getElectrons(ev, eAlgo);
    vector<MyJet> pfJets = evR->getJets(ev, jAlgo);
    MyMET met = evR->getMET(ev, metAlgo);

    // preselect objects 
    vector<int> m_init; m_init.clear();
    preSelectMuons(&m_init, pfMuons, Vertices[0], isPFlow);
    vector<int> m_init_noiso; m_init_noiso.clear();
    preSelectMuonsNoIso(&m_init_noiso, pfMuons_noiso, Vertices[0], isPFlow);
    vector<int> e_init; e_init.clear();
    preSelectElectrons(&e_init, pfElectrons, Vertices[0], isPFlow);
    vector<int> j_init; j_init.clear();
    preSelectJets(jAlgo, &j_init, pfJets, jes, jer);
    
    // clean objects //
    vector<int> e_final; e_final.clear();
    ElectronCleaning( pfElectrons, pfMuons, &e_init, &e_final, &m_init, DRMIN_ELE);
    vector<int> j_final; j_final.clear();
    vector<int> t_final; t_final.clear();
    JetCleaning(pfJets, pfMuons, pfElectrons,  &j_init, &j_final, &m_init, &e_final, DRMIN_JET);
    
    //Get MC partons
    vector<MyLorentzVector> bquarks; bquarks.clear();
    vector<MyLorentzVector> lquarks; lquarks.clear();
    MyLorentzVector mcTop, mcAntiTop;
    mcTop.SetPxPyPzE(0., 0., 0., 0.); mcAntiTop.SetPxPyPzE(0., 0., 0., 0.);
    if(!ev->isData){
      vector<MyMCParticle>allMCParticles = ev->mcParticles;
      for(size_t imc=0; imc < allMCParticles.size(); ++imc){
	    if(abs(allMCParticles[imc].pid) == 5 && allMCParticles[imc].mother.size() > 0 && (abs(allMCParticles[imc].mother[0])==6) )
	      bquarks.push_back(allMCParticles[imc].p4Gen);
	    else if(abs(allMCParticles[imc].pid) <= 4 && allMCParticles[imc].mother.size() > 0 && (abs(allMCParticles[imc].mother[0])==24 || abs(allMCParticles[imc].mother[0])==37) )
	      lquarks.push_back(allMCParticles[imc].p4Gen); 
	    else if(allMCParticles[imc].pid == 6 && allMCParticles[imc].status == 3)
	      mcTop = allMCParticles[imc].p4Gen;
	    else if(allMCParticles[imc].pid == -6 && allMCParticles[imc].status == 3)
	      mcAntiTop = allMCParticles[imc].p4Gen;
      }
    }
    
    //---------------------------------------------------//
    //apply selection cuts on leptons
    //---------------------------------------------------//
    if(m_init_noiso.size() > 0){
      int m_i_noiso = m_init_noiso[0];
      double mRelIso_no_iso = pfMuons_noiso[m_i_noiso].pfRelIso;
      fillHisto("pre_RelIso_mu",cutflowType, mRelIso_no_iso, evtWeight);
    }
    double pri_vtxs = Vertices.size();
    int nLepton = m_init.size();  // this condition proof that only muon + jet events
    if(nLepton != 1)continue;
    if( looseMuonVeto( m_init[0],pfMuons, isPFlow) ) continue;
    if( looseElectronVeto(-1,pfElectrons, isPFlow) ) continue;
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);
    
    //apply muon SF to eventWeights 
    int m_i = m_init[0];
/* 
    double musfWeight = 1.0;
    if(fabs(pfMuons[m_i].p4.eta()) < 0.9)musfWeight = muSF["sfEta1"];
    else if(fabs(pfMuons[m_i].p4.eta()) > 0.9 && fabs(pfMuons[m_i].p4.eta()) < 1.2)musfWeight = muSF["sfEta2"];
    else musfWeight = muSF["sfEta3"];
    evtWeight *= musfWeight;
    //cout<<"evtWeight musf = "<<evtWeight<<endl;
*/   
    int count_muon = m_init.size();
    ///int muCharge = pfMuons[m_i].charge;

    ///double nCutPass = 0.0;// double nCutPass_plus = 0.0; double nCutPass_minus = 0.0;
    double mRelIso = pfMuons[m_i].pfRelIso;
    fillHisto("pt_mu", cutflowType, pfMuons[m_i].p4.pt(), evtWeight);
    fillHisto("final_RelIso_mu",cutflowType, mRelIso, evtWeight);
    fillHisto("final_multi_mu",cutflowType, count_muon, evtWeight);

    // Fill histogram after trigger and one offline isolated muon and applied 2nd lepton veto
    int nJet = j_final.size();
    fillHisto("eta_mu", cutflowType, pfMuons[m_i].p4.eta(), evtWeight);
    fillHisto("phi_mu", cutflowType, pfMuons[m_i].p4.phi(), evtWeight);
    // vertex just after one lepton selection
    //double pri_vtxs = Vertices.size();
    fillHisto("nvtx", cutflowType, pri_vtxs, evtWeight);
    fillHisto("nvtx_6Kbins", cutflowType, pri_vtxs, evtWeight);
    fillHisto("rhoAll0", cutflowType, Vertices[0].rhoAll, evtWeight);
    for( std::size_t n=0; n<Vertices.size(); n++){
      fillHisto("rhoAll", cutflowType, Vertices[n].rhoAll, evtWeight);
      fillHisto("chi2", cutflowType, Vertices[n].chi2, evtWeight);
      fillHisto("ndof", cutflowType, Vertices[n].ndof, evtWeight);
      }
    //---------------------------------------------------//
    // Apply Jet Selection
    //---------------------------------------------------//

    fillHisto("multi_jet", cutflowType, nJet, evtWeight);
    if(nJet < 4)continue;  // this condition implies event should contain at least 4 jets
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      fillHisto("pt_jet", cutflowType, pfJets[ind_jet].p4.pt(), evtWeight);
      fillHisto("eta_jet", cutflowType, pfJets[ind_jet].p4.eta(), evtWeight);
      fillHisto("phi_jet", cutflowType, pfJets[ind_jet].p4.phi(), evtWeight);
    }
    fillHisto("final_multi_jet", cutflowType, nJet, evtWeight);
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);
    
    //---------------------------------------------------//
    //apply MET selection   
    //---------------------------------------------------//
    double   leptonPt(0), deltaPhi(0);
    double metPt = 0; 
    metPt= met.p4.pt();
    ///if(!metuc)metPt = metWithJESJER(pfJets, &j_final, met, jes, jer);
    ///else metPt = metWithUncl(pfJets, &j_final, pfMuons, &m_init, pfElectrons, &e_final, met, metuc);
    fillHisto("pt_met", cutflowType, metPt, evtWeight);
    fillHisto("phi_met", cutflowType, met.p4.phi(), evtWeight);
    if(metPt < 20) continue;  // Missing transverse energy cut 30 GeV(CMS) for ATLAS 20 GeV 
    fillHisto("final_pt_met", cutflowType, metPt, evtWeight);
    fillHisto("final_phi_met", cutflowType, met.p4.phi(), evtWeight);
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);

    int m_j = m_init[0];
    leptonPt = TMath::Abs(pfMuons[m_j].p4.pt());
    deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(pfMuons[m_j].p4, met.p4);
    double mt = sqrt (  2*leptonPt*metPt*(1 - cos(deltaPhi) ) ) ;
    fillHisto("wmt", cutflowType, mt, evtWeight);
    /*
    if((mt + metPt ) < 30 ) continue;  // extra condition to reduce multijet bkgs
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);
    */
    
    //---------------------------------------------------//
    //apply B-tagging, C-tagging
    //---------------------------------------------------//
    int count_CSVL = 0; int count_CSVM = 0; //int count_CSVT = 0; 
    vector<int> j_final_nob; j_final_nob.clear();
    double pfCISV = 0.0; //pfCombinedInclusiveSecondaryVertexV2BJetTags
    double pfCMVA = 0.0; //pfCombinedMVAV2BJetTags
    double pfCCvsL = 0.0;//pfCombinedCvsLJetTags
    double pfCCvsB = 0.0; //pfCombinedCvsBJetTags
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      pfCISV = pfJets[ind_jet].bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
      pfCMVA = pfJets[ind_jet].bDiscriminator["pfCombinedMVAV2BJetTags"];
      pfCCvsL= pfJets[ind_jet].bDiscriminator["pfCombinedCvsLJetTags"];
      pfCCvsB = pfJets[ind_jet].bDiscriminator["pfCombinedCvsBJetTags"];
      fillHisto("pfCISV", cutflowType, pfCISV , evtWeight); 
      fillHisto("pfCMVA", cutflowType, pfCMVA , evtWeight); 
      fillHisto("pfCCvsL", cutflowType, pfCCvsL, evtWeight); 
      fillHisto("pfCCvsB", cutflowType, pfCCvsB , evtWeight); 
      if(pfCISV > 0.5426){
        count_CSVL++;
        fillHisto("bDiscr_Loose", cutflowType+"/BTag", pfCISV, evtWeight); 
      }
      else j_final_nob.push_back(ind_jet);  
    }
    fillHisto("CSVL_count", cutflowType, count_CSVL, evtWeight);
    
    double pfCCvsL_0 = 0.0;
    double pfCCvsL_1 = 0.0;
    double pfCCvsB_0 = 0.0;
    double pfCCvsB_1 = 0.0;
    //Get the invariant mass of csbar pair
    if(j_final_nob.size() >= 2){
      int first_index = j_final_nob[0];
      int sec_index = j_final_nob[1];
      pfCCvsL_0 = pfJets[first_index].bDiscriminator["pfCombinedCvsLJetTags"];
      pfCCvsL_1 = pfJets[sec_index].bDiscriminator["pfCombinedCvsLJetTags"];
      pfCCvsB_0 = pfJets[first_index].bDiscriminator["pfCombinedCvsBJetTags"];
      pfCCvsB_1 = pfJets[sec_index].bDiscriminator["pfCombinedCvsBJetTags"];
      fillHisto("pfCCvsL_0", cutflowType, pfCCvsL_0, evtWeight); 
      fillHisto("pfCCvsL_1", cutflowType, pfCCvsL_1, evtWeight); 
      fillHisto("pfCCvsB_0", cutflowType, pfCCvsB_0 , evtWeight); 
      fillHisto("pfCCvsB_1", cutflowType, pfCCvsB_1 , evtWeight); 
      MyLorentzVector diJet = pfJets[first_index].p4 + pfJets[sec_index].p4;
      fillHisto("mjj", cutflowType, diJet.M(), evtWeight);
    }
    /*
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
            if(pfJets[ind_jet].bDiscriminator["combinedSecondaryVertexBJetTags"] > 0.679)count_CSVM++;
      bool isBtag = getBtagWithSF(pfJets[ind_jet], isData, bscale, true); 
      if(isBtag)count_CSVM++; 
    }// 13TeV, it is 0.8484, the tight one is 0.9535 

    //    if(count_CSVL <= 2) continue; // Atfirst demanding for ATLAS 3L and then 0M
    */
    if(count_CSVL <= 1) continue; // Demanding for 2M b-tagged jets
    fillHisto("CSVM_count", cutflowType, count_CSVM, evtWeight);
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight); 
    nSelEvents++;  // this is the counter for counting final number of events  

    //---------------------------------------------------//
    // add set of plots after BTag:
    //---------------------------------------------------//
    fillHisto("pt_mu", cutflowType+"/BTag", pfMuons[m_i].p4.pt(), evtWeight);
    fillHisto("eta_mu", cutflowType+"/BTag", pfMuons[m_i].p4.eta(), evtWeight);
    fillHisto("phi_mu", cutflowType+"/BTag", pfMuons[m_i].p4.phi(), evtWeight);
    fillHisto("final_RelIso_mu",cutflowType+"/BTag", mRelIso, evtWeight);
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      //double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      double jetPt = pfJets[ind_jet].p4.pt(); 
      fillHisto("pt_jet", cutflowType+"/BTag", jetPt, evtWeight);
      fillHisto("eta_jet", cutflowType+"/BTag", pfJets[ind_jet].p4.eta(), evtWeight);
      fillHisto("phi_jet", cutflowType+"/BTag", pfJets[ind_jet].p4.phi(), evtWeight);
    }
    fillHisto("final_multi_jet", cutflowType+"/BTag", nJet, evtWeight);
    fillHisto("final_pt_met", cutflowType+"/BTag", metPt, evtWeight);
    fillHisto("nvtx", cutflowType+"/BTag", pri_vtxs, evtWeight);
    fillHisto("nvtx_6Kbins", cutflowType+"/BTag", pri_vtxs, evtWeight);
    fillHisto("rhoAll0", cutflowType+"/BTag", Vertices[0].rhoAll, evtWeight);
    for(std::size_t n=0; n<Vertices.size(); n++){
      fillHisto("rhoAll", cutflowType+"/BTag", Vertices[n].rhoAll, evtWeight);
      fillHisto("chi2", cutflowType+"/BTag", Vertices[n].chi2, evtWeight);
      fillHisto("ndof", cutflowType+"/BTag", Vertices[n].ndof, evtWeight);
      }
    fillHisto("wmt", cutflowType+"/BTag", mt, evtWeight);

    input_count++;
    if(input_count%10==0)
    cout << "input count: "<< input_count << endl;
    //if(i > 5000) break;
    }//event loop
      
  outfile_ << "Number of times HadP and HadQ matches with jets" << matchjetcount << endl;
  outfile_ << "total number of selected events " << nSelEvents <<endl; 
  outfile_ << "No of times three pair jet matched:    " << threepairjet << endl;
  outfile_ << "No of correct jet pair assign to W:   " << matchedJet_q << endl; 
  outfile_ << "No of not correct jet pair assign to W:   " << not_matchedJet_q << endl;
  outfile_ << "No of correct jet pair not assign to W but |eta| < 2.5 and pt > 25:   " << matched_quark_eta_pt << endl;
  outfile_ << "No of correct jet matched to b (from top) :   " << matchedJet_b << endl;
  outfile_ << "Average weights for top-reweighting : "<<TotalTopPtWeights<<"/"<<TotalLplusJEvents<<" = "<<TotalTopPtWeights/TotalLplusJEvents<<endl; 
  outfile_ << "Efficiency Uncertainty on the SV finding : "<<TotalSVEff<<"/"<<TotalEvtsWithSV<<" = "<<TotalSVEff/TotalEvtsWithSV<<endl;
  fillHisto("AvTopPtWeight", cutflowType, 1.0, TotalTopPtWeights/TotalLplusJEvents);
  fillHisto("SVEffUncert", cutflowType, 2.0, TotalSVEff/TotalEvtsWithSV);
  fillHisto("SVEffUncert", cutflowType, 1.0, TotalEvtsWithSV);
  fillHisto("SVEffUncert", cutflowType, 0.0, TotalSVEff);
  f->Close(); 
  delete f;
}

void hplusAnalyzer::processEvents(){ 
  
  //Data, MC sample from lxplus and T2
  //CutFlowAnalysis("TTJets_ntuple_MuChannel.root", "PF", "TTJets_MuMC_check"); 
  
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/rverma/ntuple_MuMC_nokfit_20170528/MuMC_20170528/DY1JetsToLL_MuMC_20170528/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DY1JetsToLL_MuMC_20170528/170528_112855/0000/DY1JetsToLL_MuMC_20170528_Ntuple_1.root", "PF","DY1");
  
  //condor submission
  CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", "PF", "outputFile");
} 

float hplusAnalyzer::reweightHEPNUPWJets(int hepNUP) {

  int nJets = hepNUP-5;
  if(nJets==0)      return 2.07;
  else if(nJets==1) return 0.226;
  else if(nJets==2) return 0.119;
  else if(nJets==3) return 0.0562;
  else if(nJets>=4) return 0.0671;
  else return 1 ;
}

float hplusAnalyzer::reweightHEPNUPDYJets(int hepNUP){

  int nJets = hepNUP-5;
  if(nJets==0)      return 0.117;
  else if(nJets==1) return 0.0164;
  else if(nJets==2) return 0.0167;
  else if(nJets==3) return 0.0167;
  else if(nJets>=4) return 0.0128;
  else return 1 ;
}

