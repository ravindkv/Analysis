
///////////////////////
// Electron Channel New
///////////////////////

#include "hplusAnalyzer.h"
#include <map>

using namespace std;
void hplusAnalyzer::CutFlowAnalysis(TString url, string myKey, string evtType){
  
  TString outFile("13TeV/outputDir/");
  TString Filename_ = outFile+evtType+"_Anal.root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel( 9 );
  
  TString debug_Filename_ = Filename_+"_debug.txt";
  string debug_file(debug_Filename_);
  
  //check if the input file is MC or Data  
  Reader *evR_;  
  evR_ = new Reader();
  TFile *f_ = TFile::Open(url);
  int nEntries = evR_->AssignEventTreeFrom(f_);
  MyEvent *ev_;
  ev_ = evR_->GetNewEvent(1);

  CutFlowProcessor(url, myKey, "base", outFile_);
  //CutFlowProcessor(url, myKey, "baseLowMET", outFile_);
  /*
  //to estimate unc in the data-driven qcd 
  CutFlowProcessor(url, myKey, "baseIso20HighMET", outFile_);
  CutFlowProcessor(url, myKey, "baseIso20LowMET", outFile_);
  //---------------------------------------------------//
  //for systematics (all sys in one go)
  //---------------------------------------------------//  
  if(!ev_->isData){ 
    CutFlowProcessor(url, myKey, "JESPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "JESMinus", 	outFile_);
    CutFlowProcessor(url, myKey, "JERPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "JERMinus", 	outFile_);
    //CutFlowProcessor(url, myKey, "METUCPlus", 	outFile_);
    //CutFlowProcessor(url, myKey, "METUCMinus", 	outFile_);
    CutFlowProcessor(url, myKey, "bTagPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "bTagMinus", 	outFile_);
    CutFlowProcessor(url, myKey, "TopPtPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "TopPtMinus", 	outFile_);
  }
  */
  outFile_->Write(); 
  outFile_->Close();
  f_->Close();
  delete f_;
}

//---------------------------------------------------//
//Process the cuts, event by event
//---------------------------------------------------//  
void hplusAnalyzer::CutFlowProcessor(TString url,  string myKey, TString cutflowType, TFile *outFile_){
  int input_count = 0;
  bool isPFlow = (myKey.find("PFlow") != string::npos) ? true : false;
  string eAlgo("Electrons"), mAlgo("Muons"), jAlgo("Jets"), metAlgo("METs");
  
  //Uncertainty variations, JES, JER, MET unclustered, bTag
  int jes = 0, jer = 0, metuc = 0, bscale = 0, minMET =20, minMT =0;
  //to estimate unc in the data-driven qcd 
  bool isLowMET = false, isIso20 = false;

  if(cutflowType.Contains("JESPlus"))jes = 1;
  else if (cutflowType.Contains("JESMinus"))jes = -1;
  else if (cutflowType.Contains("JERPlus"))jer = 1;
  else if (cutflowType.Contains("JERMinus"))jer = -1;
  else if (cutflowType.Contains("METUCPlus"))metuc = 1;
  else if (cutflowType.Contains("METUCMinus"))metuc = -1;
  else if (cutflowType.Contains("bTagPlus"))bscale = 1;
  else if (cutflowType.Contains("bTagMinus"))bscale = -1; 
  //to estimate unc in the data-driven qcd 
  else if (cutflowType.Contains("baseIso")){ 
    if (cutflowType.Contains("Iso20HighMET"))isIso20 = true; 
    if (cutflowType.Contains("Iso20LowMET")){isIso20 = true; isLowMET= true;}
  }
  else if (cutflowType.Contains("LowMET"))isLowMET = true;
  
  evR = new Reader();
  TFile *f = TFile::Open(url);
  if(f==0) return ;
  if(f->IsZombie()) { f->Close(); return; }
  
  //---------------------------------------------------//
  //get initial number of events, from ntuples
  //store initial informations, in a txt file
  //---------------------------------------------------//
  double lumiTotal = 34698+758;
  int nEntries = evR->AssignEventTreeFrom(f);
  if(nEntries == 0) {return; }
  TH1F* inputcf = (TH1F*)(f->Get("allEventsFilter/totalEvents"));
  double initialEvents = inputcf->GetBinContent(1);
  cout<<"\033[01;32m input file: \033[00m"<<url<<"\n"<<endl;
  fillHisto(outFile_, cutflowType, "", "totalEvents", 10, 0, 10000000000, initialEvents, 1 );
  MyEvent *ev;
  int nTriggEvent = 0, nSelEvents = 0, matchjetcount= 0, threepairjet = 0;
  double nVerticesFailCount = 0.0;
  double matchedJet_q = 0.0, matchedJet_b = 0.0, matched_quark_eta_pt = 0.0, not_matchedJet_q = 0.0;
  double TotalTopPtWeights = 0, TotalLplusJEvents = 0; 
  double TotalSVEff = 0, TotalEvtsWithSV = 0.;
  
  //---------------------------------------------------//
  //BTag SF: read CSV file for SF, 2D histos for eff 
  //---------------------------------------------------//      
  //https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
  const std::string & filename 		= "stack/CSVv2_Moriond17_B_H.csv";
  const std::string & tagger 		= "CSVv2";
  const std::string & measurementType 	= "comb";
  const std::string & sysType 		= "central"; 
  if(bscale==1) const std::string &sysType 		= "up"; 
  if(bscale==-1)const std::string &sysType 		= "down"; 
  const std::vector<std::string> & otherSysTypes = {"up", "down"};
  //b-quark
  BTagCalibrationReader readCSVbL= readCSVfile(filename, tagger, BTagEntry::OP_MEDIUM,
    	      measurementType, sysType, otherSysTypes, BTagEntry::FLAV_B);
  //c-quark
  BTagCalibrationReader readCSVcL= readCSVfile(filename, tagger, BTagEntry::OP_MEDIUM,
    	      measurementType, sysType, otherSysTypes, BTagEntry::FLAV_C);
  //other(light) quarks and gluon
  BTagCalibrationReader readCSVlL= readCSVfile(filename, tagger, BTagEntry::OP_MEDIUM,
    	      "incl", sysType, otherSysTypes, BTagEntry::FLAV_UDSG);
  
  //getBTagEffHistos(f);
  TString histPath("myMiniTreeProducer/Jets/");
  TH2D* h2_BTaggingEff_Denom_b 		= (TH2D*)(f->Get(histPath+"h2_BTaggingEff_Denom_b"));
  TH2D* h2_BTaggingEff_Denom_c 		= (TH2D*)(f->Get(histPath+"h2_BTaggingEff_Denom_c"));
  TH2D* h2_BTaggingEff_Denom_udsg 	= (TH2D*)(f->Get(histPath+"h2_BTaggingEff_Denom_udsg")); 
  TH2D* h2_BTaggingEff_Num_bL 		= (TH2D*)(f->Get(histPath+"h2_BTaggingEff_Num_bL"));
  TH2D* h2_BTaggingEff_Num_cL 		= (TH2D*)(f->Get(histPath+"h2_BTaggingEff_Num_cL"));
  TH2D* h2_BTaggingEff_Num_udsgL 	= (TH2D*)(f->Get(histPath+"h2_BTaggingEff_Num_udsgL")); 
  
  //---------------------------------------------------//
  //loop over each event, of the ntuple
  //---------------------------------------------------//
  double kfCount = 0;
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
        if(i < 1){
        }
        evtWeight *= weightK;  
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
        if(i < 1){
        }
      }
      //lumi weight
      else {
      double sampleWeight(1.0);
      sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
      evtWeight *= sampleWeight; 
      fillHisto(outFile_, cutflowType, "", "MuPt30JetPt25MET20MT20", 10, 0, 1000, sampleWeight, 1 );
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
    // apply top re-weighting
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
    //---------------------------------------------------//
    double topPtWeights_offline = 1.0;
    if(!ev->isData){
      string sampleName = ev->sampleInfo.sampleName;
      if(sampleName.find("Hplus") != string::npos ||
		      sampleName.find("TTJetsM") != string::npos || 
		      sampleName.find("TTJetsP") != string::npos)
      {
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
    fillHisto(outFile_, cutflowType, "", "SF_topPtWeights", 1000, 0, 3, topPtWeights_offline, 1 );
    //cout<<"topPtWeights_offline = "<<topPtWeights_offline<<endl;
    TotalLplusJEvents++;
    evtWeight *= topPtWeights_offline; //Multiply to the total weights
    
    //---------------------------------------------------//
    //apply muon triggers
    //---------------------------------------------------//
    
    bool passTrig = false;
    vector<string> trig = ev->hlt;
    for(size_t it = 0; it < trig.size(); it++){
      if(trig[it].find("HLT_Ele27_WPTight_Gsf") != string::npos) passTrig = true;
    }
    if(!passTrig){
    //cout << "not satisfying trigger" << endl;
      continue;
    }
    nTriggEvent++;
    double nCutPass = 1.0;
    double nCutPass_NonIso = 1.0;
    fillHisto(outFile_, cutflowType+"/Iso", "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    fillHisto(outFile_, cutflowType+"/NonIso", "", "cutflow", 20, 0.5, 20.5, nCutPass_NonIso, evtWeight );
   
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
    vector<MyElectron> pfElectrons = evR->getElectrons(ev, eAlgo);
    vector<MyJet> pfJets = evR->getJets(ev, jAlgo);
    MyMET met = evR->getMET(ev, metAlgo);

    // preselect objects 
    vector<int> m_init; m_init.clear();
    double u1 	= gRandom->Rndm();//used for rochester corrections
    double u2 	= gRandom->Rndm();
    preSelectMuons(&m_init, pfMuons, Vertices[0], ev->isData, u1, u2, 0, 0);
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
    //get KinFit objects
    //---------------------------------------------------//
    vector<MyLorentzVector> kfJets; kfJets.clear();
    vector<MyLorentzVector> kfJetsLepB; kfJetsLepB.clear();
    vector<MyLorentzVector> kfLepton; kfLepton.clear();
    vector<MyLorentzVector> kfMet; kfMet.clear();
  
    double chi2OfKinFit=999.;
    double statusOfKinFit=-99;
    double probOfKinFit=-99;
    vector<MyKineFitParticle> allKineFitParticles = ev->KineFitParticles;
    for(size_t imk=0; imk < allKineFitParticles.size(); imk++){
      string labelName = "";
      if(cutflowType.Contains("JESPlus"))labelName="JESUp";
      else if(cutflowType.Contains("JESMinus"))labelName="JESDown";
      else if(cutflowType.Contains("JERPlus"))labelName="JERUp";
      else if(cutflowType.Contains("JERMinus"))labelName="JERDown";
      if(labelName=="JESUp" || labelName=="JESDown" || labelName=="JERUp" ||labelName=="JERDown"){
        if(allKineFitParticles[imk].labelName.find(labelName) != string::npos ){
          if(allKineFitParticles[imk].partName.find("PartonsHadP") != string::npos ||
             allKineFitParticles[imk].partName.find("PartonsHadQ") != string::npos ||
             allKineFitParticles[imk].partName.find("PartonsHadB") != string::npos
             ){
            kfJets.push_back(allKineFitParticles[imk].p4);
          }
          else if(allKineFitParticles[imk].partName.find("Leptons") != string::npos ){
            kfLepton.push_back(allKineFitParticles[imk].p4);
          }
          else if(allKineFitParticles[imk].partName.find("Neutrinos") != string::npos ){
            kfMet.push_back(allKineFitParticles[imk].p4);
          }
          else if(allKineFitParticles[imk].partName.find("PartonsLepB") != string::npos ){
            kfJetsLepB.push_back(allKineFitParticles[imk].p4);
          }
          ///if(imk<1){
          chi2OfKinFit = allKineFitParticles[imk].chi2OfFit;
          statusOfKinFit = allKineFitParticles[imk].statusOfFit;
          probOfKinFit = allKineFitParticles[imk].probOfFit;
	  ///}
        }
      }
      else{
        if(allKineFitParticles[imk].labelName.find("JESUp") == string::npos && allKineFitParticles[imk].labelName.find("JESDown") == string::npos && allKineFitParticles[imk].labelName.find("JERUp") == string::npos && allKineFitParticles[imk].labelName.find("JERDown") == string::npos){
          if(allKineFitParticles[imk].partName.find("PartonsHadP") != string::npos ||
             allKineFitParticles[imk].partName.find("PartonsHadQ") != string::npos ||
             allKineFitParticles[imk].partName.find("PartonsHadB") != string::npos
             ){
            kfJets.push_back(allKineFitParticles[imk].p4);
          }
          else if(allKineFitParticles[imk].partName.find("Leptons") != string::npos ){
            kfLepton.push_back(allKineFitParticles[imk].p4);
          }
          else if(allKineFitParticles[imk].partName.find("Neutrinos") != string::npos ){
            kfMet.push_back(allKineFitParticles[imk].p4);
          }
          else if(allKineFitParticles[imk].partName.find("PartonsLepB") != string::npos ){
            kfJetsLepB.push_back(allKineFitParticles[imk].p4);
          }
          if(imk<1){
            chi2OfKinFit = allKineFitParticles[imk].chi2OfFit;
            statusOfKinFit = allKineFitParticles[imk].statusOfFit;
            probOfKinFit = allKineFitParticles[imk].probOfFit;
          }
        }
      }
    }
    
    //---------------------------------------------------//
    //apply selection cuts on leptons
    //---------------------------------------------------//
    if(e_final.size() > 0){
      int e_i = e_final[0];
      double eRelIso = pfElectrons[e_i].relCombPFIsoEA/pfElectrons[e_i].p4.pt();
      fillHisto(outFile_, cutflowType+"/Iso", "", "RelIso", 100, 0, 1, eRelIso, evtWeight);
      fillHisto(outFile_, cutflowType+"/NonIso", "", "RelIso", 100, 0, 1, eRelIso, evtWeight);
    }
    int nEle = e_final.size();
    double pri_vtxs = Vertices[0].totVtx;
    //select only one electron
    if(nEle != 1)continue;
    int e_i = e_final[0];
    //loose electron veto
    if(looseElectronVeto(e_final[0], pfElectrons, Vertices[0], isPFlow)) continue;
    nCutPass = 2; 
    nCutPass_NonIso++;
    fillHisto(outFile_, cutflowType+"/Iso", "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    fillHisto(outFile_, cutflowType+"/NonIso", "", "cutflow", 20, 0.5, 20.5, nCutPass_NonIso, evtWeight );
     
    //events should not have any muon
    if(looseMuonVeto( -1, pfMuons, isPFlow) ) continue;
    nCutPass = 3; 
    nCutPass_NonIso++;
    fillHisto(outFile_, cutflowType+"/Iso", "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    fillHisto(outFile_, cutflowType+"/NonIso", "", "cutflow", 20, 0.5, 20.5, nCutPass_NonIso, evtWeight );
    int count_ele = e_final.size();
    ///int muCharge = pfMuons[m_i].charge;
    
    //---------------------------------------------------//
    //apply Electron SF to eventWeights 
    //---------------------------------------------------//
    double eleSF =1.0;
    //Reco, ID, trigger	
    double ele_recoSF 		= getEleSF(h2_ele_recoSF, pfElectrons[e_i].eleSCEta, pfElectrons[e_i].p4.pt());
    double ele_medium_idSF  	= getEleSF(h2_ele_medium_idSF, pfElectrons[e_i].eleSCEta, pfElectrons[e_i].p4.pt());
    double ele_trigSF 		= getEleTrigSF(h2_ele_trigSF, pfElectrons[e_i].p4.pt(), pfElectrons[e_i].eleSCEta);
    if(!ev->isData){
      eleSF = ele_recoSF*ele_medium_idSF*ele_trigSF;
    }
    evtWeight *= eleSF;
    nCutPass = 4;
    nCutPass_NonIso++;
    fillHisto(outFile_, cutflowType+"/Iso", "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    fillHisto(outFile_, cutflowType+"/NonIso", "", "cutflow", 20, 0.5, 20.5, nCutPass_NonIso, evtWeight );
    //---------------------------------------------------//
    // Iso(<0.12) and Non-iso(>0.12) region 
    //---------------------------------------------------//
    bool noisofound = false;
    bool isofound = false;
    double tmp_iso = pfElectrons[e_i].relCombPFIsoEA;
    fillHisto(outFile_, cutflowType, "", "RelIso_ele", 100, 0, 10, tmp_iso, evtWeight );
    string cutflowType_(cutflowType);
    if(isIso20){
      if(tmp_iso <= 0.20) cutflowType_ = cutflowType+"/Iso";
      if(tmp_iso > 0.20 && tmp_iso <= 0.4) cutflowType_ = cutflowType+"/NonIso";
    }
    else{
      if(tmp_iso <= 0.08) cutflowType_ = cutflowType+"/Iso";
      if(tmp_iso > 0.08 && tmp_iso <= 0.3) cutflowType_ = cutflowType+"/NonIso";
    }
    if(tmp_iso > 0.3) continue;
    nCutPass = 5;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    
    double eRelIso = pfElectrons[e_i].relCombPFIsoEA;
    //double muonPt = pfMuons[m_i].p4.pt();
    double elePt = pfElectrons[e_i].p4.pt();
    
    //---------------------------------------------------//
    // Apply Jet Selection
    //---------------------------------------------------//
    int count_jets = j_final.size();
    if(count_jets < 4)continue;  // events should have 4 or more jets
    nCutPass = 6;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    
    //---------------------------------------------------//
    //apply MET selection   
    //---------------------------------------------------//
    double metPt = 0; 
    metPt = metWithJESJER(pfJets, &j_final, met, jes, jer);
    //if(metuc==0)metPt = metWithJESJER(pfJets, &j_final, met, jes, jer);
    //else metPt = metWithUncl(pfJets, &j_final, pfMuons, &m_init, pfElectrons, &e_final, met, metuc);
    
    if(isLowMET){
      if(metPt > minMET) continue;  
    }
    else if(metPt < minMET) continue;  // Missing transverse energy cut 30 GeV(CMS) for ATLAS 20 GeV 
    nCutPass = 7;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    
    //Transverse mass b/w lepton and MET
    double deltaPhi(0);
    //leptonPt = TMath::Abs(pfMuons[m_i].p4.pt());
    deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(pfElectrons[e_i].p4, met.p4);
    double mt = sqrt (  2*elePt*metPt*(1 - cos(deltaPhi) ) ) ;
    if(mt < minMT) continue; // mt cuts
    nCutPass = 8;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    
    //---------------------------------------------------//
    //apply B-tagging, C-tagging
    //---------------------------------------------------//
    vector<int> j_final_nob; j_final_nob.clear();
    vector<int> j_final_b; j_final_b.clear();
    vector<double> bdiscr; bdiscr.clear();
    double pfCISV = 0.0; //pfCombinedInclusiveSecondaryVertexV2BJetTags
    double pfCMVA = 0.0; //pfCombinedMVAV2BJetTags
    double pfCCvsL = 0.0;//pfCombinedCvsLJetTags
    double pfCCvsB = 0.0; //pfCombinedCvsBJetTags
    double pfCCvsL_0 = 0.0;
    double pfCCvsL_1 = 0.0;
    double pfCCvsB_0 = 0.0;
    double pfCCvsB_1 = 0.0;

    //LOOSE BTAG
    int count_CSVL_SF = 0; 
    int count_CSVL = 0;
    ///bool isBtagL = true; bool isBtagM = true; isBtagT = true; 
    bool isBtag = false;
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      pfCISV = pfJets[ind_jet].bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
      pfCMVA = pfJets[ind_jet].bDiscriminator["pfCombinedMVAV2BJetTags"];
      pfCCvsL= pfJets[ind_jet].bDiscriminator["pfCombinedCvsLJetTags"];
      pfCCvsB = pfJets[ind_jet].bDiscriminator["pfCombinedCvsBJetTags"];
      fillHisto(outFile_, cutflowType_, "BTag", "pfCISV", 50, 0, 2, pfCISV, evtWeight );
      fillHisto(outFile_, cutflowType_, "BTag", "pfCMVA", 100, -2, 2, pfCMVA, evtWeight );
      fillHisto(outFile_, cutflowType_, "BTag", "pfCCvsL", 100, -2, 2, pfCCvsL, evtWeight );
      fillHisto(outFile_, cutflowType_, "BTag", "pfCCvsB", 100, -2, 2, pfCCvsB, evtWeight );
      //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
      //b-quark
      if(abs(pfJets[ind_jet].partonFlavour) ==5)
        isBtag = getBtagWithSF(readCSVbL, h2_BTaggingEff_Num_bL, h2_BTaggingEff_Denom_b, pfJets[ind_jet], ev->isData, bscale); 
      //c-quark
      else if(abs(pfJets[ind_jet].partonFlavour) ==4) 
        isBtag = getBtagWithSF(readCSVcL, h2_BTaggingEff_Num_cL, h2_BTaggingEff_Denom_c, pfJets[ind_jet], ev->isData, bscale); 
      //other quarks and gluon
      else isBtag = getBtagWithSF(readCSVlL, h2_BTaggingEff_Num_udsgL, h2_BTaggingEff_Denom_udsg, pfJets[ind_jet], ev->isData, bscale); 
      
      if(isBtag){
        count_CSVL_SF++; 
        double jetPt = jetPtWithJESJER(pfJets[ijet], jes, jer);
        fillHisto(outFile_, cutflowType_, "BTag", "pt_bjet", 50, 0, 500, jetPt, evtWeight );
        fillHisto(outFile_, cutflowType_, "BTag", "eta_bjet", 50, -5, 5, pfJets[ijet].p4.eta(), evtWeight );
        j_final_b.push_back(ind_jet);
        bdiscr.push_back(pfCISV);
      }
      else j_final_nob.push_back(ind_jet); 
      //if(pfCISV >0.5426) count_CSVL ++;
      if(pfCISV >0.8484) count_CSVL ++;
      //if(pfCISV >0.9535) count_CSVL ++;
    }
    if(count_CSVL >=2 ){
    nCutPass = 9;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    }
    fillHisto(outFile_, cutflowType_, "BTag", "CSVL_count", 50, 0, 10, count_CSVL_SF, evtWeight );
    if(count_CSVL_SF <= 1) continue; // Demanding for 2L b-tagged jets
    nCutPass = 10;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    
    //---------------------------------------------------//
    //invariant mass of c sbar
    //---------------------------------------------------//
    //sort j_final_b w.r.t b-discriminator value(ascending order)
    std::map<double, int> bdiscr_sorted_bjets;
    for(unsigned long k=0; k<j_final_b.size(); k++){
      bdiscr_sorted_bjets.insert(pair <double, int> (bdiscr[k],j_final_b[k])); 
    }
    map <double, int> :: iterator bdiscr_itr;
    int index_of_2nd_bjet;
    int index_of_1st_bjet;
    vector<int> index_of_other_bjets;
    int total_bjets = j_final_b.size();
    for(bdiscr_itr = bdiscr_sorted_bjets.begin(); bdiscr_itr != bdiscr_sorted_bjets.end(); ++bdiscr_itr){
       total_bjets --;
       if(total_bjets==1) index_of_2nd_bjet = bdiscr_itr->second;  
       else if(total_bjets==0) index_of_1st_bjet = bdiscr_itr->second;  
       else index_of_other_bjets.push_back(bdiscr_itr->second);
    }
    //mjj will involve 2 non-bjet, highest pt jets
    if(j_final_b.size()==2){
      if(j_final_nob.size() >= 2){
        int index_of_1st_mjj = j_final_nob[0];
        int index_of_2nd_mjj = j_final_nob[1];
        MyLorentzVector diJet = pfJets[index_of_1st_mjj].p4 + pfJets[index_of_2nd_mjj].p4;
        fillHisto(outFile_, cutflowType_, "BTag", "mjj", 200, 0, 1000, diJet.M(), evtWeight );
      }
    }
    //Arrange other bjets and non-bjets in pt order in a list
    //mjj will involve 2 highest pt jets in this list
    else{ 
      std::map<double, int> pt_sorted_jets;
      for(unsigned long k=0; k<j_final_nob.size(); k++){
        double pt_of_nobjet = pfJets[j_final_nob[k]].p4.pt(); 
        pt_sorted_jets.insert(pair <double, int> (pt_of_nobjet, j_final_nob[k])); 
      }
      for(unsigned long k=0; k<index_of_other_bjets.size(); k++){
        double pt_of_other = pfJets[index_of_other_bjets[k]].p4.pt();
        pt_sorted_jets.insert(pair <double, int> (pt_of_other, index_of_other_bjets[k])); 
      }
      //select two highest pt jets
      int index_of_1st_mjj = 0;
      int index_of_2nd_mjj = 0;
      int total_jets_for_mjj = pt_sorted_jets.size();
      map <double, int> :: iterator itr_pt;
      for(itr_pt = pt_sorted_jets.begin(); itr_pt != pt_sorted_jets.end(); ++itr_pt){
         total_jets_for_mjj --;
         if(total_jets_for_mjj==1) index_of_2nd_mjj = itr_pt->second;  
         if(total_jets_for_mjj==0) index_of_1st_mjj = itr_pt->second;  
      }
      MyLorentzVector diJet = pfJets[index_of_1st_mjj].p4 + pfJets[index_of_2nd_mjj].p4;
      fillHisto(outFile_, cutflowType_, "BTag", "mjj", 200, 0, 1000, diJet.M(), evtWeight );
    }	

    //---------------------------------------------------//
    // add set of plots after BTag:
    //---------------------------------------------------//
    //fillHisto("pt_mu", cutflowType_+"/BTag", pfMuons[m_i].p4.pt(), evtWeight);
    fillHisto(outFile_, cutflowType_, "BTag","pt_ele", 50, 0, 500, elePt, evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","eta_ele", 50, -5, 5, pfElectrons[e_i].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","phi_ele", 50, -5, 5, pfElectrons[e_i].p4.phi(), evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","final_RelIso_ele", 100, 0, 1, eRelIso, evtWeight );
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      fillHisto(outFile_, cutflowType_, "BTag","pt_jet", 50, 0, 500, jetPt, evtWeight );
      fillHisto(outFile_, cutflowType_, "BTag","eta_jet", 50, -5, 5, pfJets[ind_jet].p4.eta(), evtWeight );
      fillHisto(outFile_, cutflowType_, "BTag","phi_jet", 50, -5, 5, pfJets[ind_jet].p4.phi(), evtWeight );
    }
    fillHisto(outFile_, cutflowType_, "BTag","final_multi_jet", 15, 0, 15, count_jets, evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","final_pt_met", 50, 0, 500, metPt, evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","wmt", 50, 0, 500, mt, evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","nvtx", 100, 0, 100, pri_vtxs, evtWeight );
    for(std::size_t n=0; n<Vertices.size(); n++){
      fillHisto(outFile_, cutflowType_, "BTag","rhoAll", 100, 0, 100, Vertices[n].rhoAll, evtWeight );
    }
    input_count++;
    if(input_count%10==0)
    cout << "input count iso: "<< input_count << endl;
    //if(i > 2000) break;
  

    //kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk//
    // 		add set of plots after KinFit: 		    //
    //kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk//
    //make sure that the fit converges 
    if(statusOfKinFit !=0) continue ; 
    nCutPass =11;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
  
    //---------------------------------------------------//
    //get KF muon by matching with PF muons
    //---------------------------------------------------//
    bool foundkfEle = false;
    double dR =DeltaR(pfElectrons[e_i].p4 , kfLepton[0]);
    if(kfLepton.size()>0){
      if(dR< 0.2)foundkfEle = true;
    }
    if(!foundkfEle) continue; 
    
    //---------------------------------------------------//
    //Apply same Pt selections on PF and KF jets
    //---------------------------------------------------//
    if(kfJets[0].pt() < 25) 	continue;
    if(kfJets[1].pt() < 25) 	continue;
    if(kfJets[2].pt() < 25) 	continue;
    if(kfJetsLepB[0].pt() < 25) continue;
    nCutPass = 12;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    
    //---------------------------------------------------//
    //select maximum b-tag discriminator jet in KF
    //---------------------------------------------------//
    unsigned long maxBtagJet = -1;
    double maxBDiscr = -999.;
    int count_kfJets = 0;
    fillHisto(outFile_, cutflowType_, "KinFit","multi_kfjet", 100, 0, 20, kfJets.size(), evtWeight );
    for(unsigned long ik = 0; ik < kfJets.size(); ik++){
      fillHisto(outFile_, cutflowType_, "KinFit","pt_kf_jets", 500, 0, 500, kfJets[ik].pt(), evtWeight );
      for(size_t ij = 0; ij < j_final.size(); ij++){
        int ind_ij = j_final[ij];
        if(DeltaR(kfJets[ik], pfJets[ind_ij].p4) < 0.2){
          double discr = pfJets[ind_ij].bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
          if(discr > maxBDiscr){
            maxBDiscr = discr;
            maxBtagJet = ik;
          }
        }
      }
    }
    //---------------------------------------------------//
    //make sure that each event has 2 light jets
    //---------------------------------------------------//
    vector<MyLorentzVector> kfLightJets; kfLightJets.clear();
    unsigned long zero = 0; 
    if(kfJets.size() >=3 && maxBtagJet >= zero){
      for(unsigned long ik = 0; ik < kfJets.size(); ik++){
        if(ik != maxBtagJet)kfLightJets.push_back(kfJets[ik]);
	else fillHisto(outFile_, cutflowType_, "KinFit","pt_bjet", 50, 0, 500, kfJets[ik].pt(), evtWeight );
      }
    }
    fillHisto(outFile_, cutflowType_, "KinFit","pt_bjetL", 50, 0, 500, kfJetsLepB[0].pt(), evtWeight );
    fillHisto(outFile_, cutflowType_, "KinFit", "chi2OfKinFit", 100, 0, 100, chi2OfKinFit, evtWeight );
    fillHisto(outFile_, cutflowType_, "KinFit", "probOfKinFit", 100, 0, 1, probOfKinFit, evtWeight );
    if(kfLightJets.size() < 2) continue; 
    bool match_j1 = false, match_j2 = false;
    int indexForCTag0 = 0, indexForCTag1 = 0;
    for(size_t ij = 0; ij < j_final.size(); ij++){
      int ind_ij = j_final[ij];
      if(DeltaR(kfLightJets[0], pfJets[ind_ij].p4) < 0.2){
        match_j1=true;
        indexForCTag0 = ind_ij;
      }
      if(DeltaR(kfLightJets[1], pfJets[ind_ij].p4) < 0.2){
        match_j2=true;
        indexForCTag1 = ind_ij;
      }
    }
    //---------------------------------------------------//
    //fill histos after c, sbar are found
    //---------------------------------------------------//
    if(!match_j1) continue;
    if(!match_j2) continue;
    kfCount++;
    nCutPass = 13;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    MyLorentzVector diJet = kfLightJets[0]+kfLightJets[1];
    fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit", 200, 0, 1000, diJet.mass(), evtWeight );
    fillHisto(outFile_, cutflowType_, "KinFit","pt_ele", 50, 0, 500, elePt, evtWeight );
    fillHisto(outFile_, cutflowType_, "KinFit","eta_ele", 50, -5, 5, pfElectrons[e_i].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "KinFit","phi_ele", 50, -5, 5, pfElectrons[e_i].p4.phi(), evtWeight );
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      fillHisto(outFile_, cutflowType_, "KinFit","pt_jet", 50, 0, 500, jetPt, evtWeight );
      fillHisto(outFile_, cutflowType_, "KinFit","eta_jet", 50, -5, 5, pfJets[ind_jet].p4.eta(), evtWeight );
      fillHisto(outFile_, cutflowType_, "KinFit","phi_jet", 50, -5, 5, pfJets[ind_jet].p4.phi(), evtWeight );
    }
    for( std::size_t n=0; n<Vertices.size(); n++){
      fillHisto(outFile_, cutflowType_, "KinFit","rhoAll", 100, 0, 100, Vertices[n].rhoAll, evtWeight );
    }
    fillHisto(outFile_, cutflowType_, "KinFit","final_multi_jet", 15, 0, 15, count_jets, evtWeight );
    fillHisto(outFile_, cutflowType_, "KinFit","final_pt_met", 50, 0, 500, metPt, evtWeight );
    fillHisto(outFile_, cutflowType_, "KinFit","wmt", 50, 0, 500, mt, evtWeight );
    fillHisto(outFile_, cutflowType_, "KinFit","nvtx", 100, 0, 100, pri_vtxs, evtWeight );
    //c jet 
    fillHisto(outFile_, cutflowType_, "KinFit","pt_kfjet0", 50, 0, 500, kfLightJets[0].pt(), evtWeight );
    fillHisto(outFile_, cutflowType_, "KinFit","eta_kfjet0", 50, -5, 5, kfLightJets[0].eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "KinFit","phi_kfjet0", 50, -5, 5, kfLightJets[0].phi(), evtWeight );
    //s bar jet
    fillHisto(outFile_, cutflowType_, "KinFit","pt_kfjet1", 50, 0, 500, kfLightJets[1].pt(), evtWeight );
    fillHisto(outFile_, cutflowType_, "KinFit","eta_kfjet1", 50, -5, 5, kfLightJets[1].eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "KinFit","phi_kfjet1", 50, -5, 5, kfLightJets[1].phi(), evtWeight );
    
    //---------------------------------------------------//
    // category depending on C-tagging
    //---------------------------------------------------//
    double pfCCvsL0 = pfJets[indexForCTag0].bDiscriminator["pfCombinedCvsLJetTags"];
    double pfCCvsL1 = pfJets[indexForCTag1].bDiscriminator["pfCombinedCvsLJetTags"];
    double pfCCvsB0 = pfJets[indexForCTag0].bDiscriminator["pfCombinedCvsBJetTags"]; 
    double pfCCvsB1 = pfJets[indexForCTag1].bDiscriminator["pfCombinedCvsBJetTags"];
    //CharmL-> (CvsL: -0.48 and CvsB: -0.17)
    MyLorentzVector diJet_tag = kfLightJets[0]+kfLightJets[1];
    /*
    if((pfCCvsL0 > dCvsL || pfCCvsL1 >dCvsL) && (pfCCvsB0 >dCvsB || pfCCvsB1 >dCvsB)){
    */
    bool isCTagL = false; //loose
    double dCCvsL = -0.5426;
    bool isCTagM1 = false; //medium
    double dCCvsM1 = -0.4426;
    bool isCTagM2 = false; //medium
    double dCCvsM2 = -0.2426;
    bool isCTagM3 = false; //medium
    double dCCvsM3 = 0.0;
    bool isCTagM4 = false; //medium
    double dCCvsM4 = 0.2426;
    bool isCTagM5 = false; //medium
    double dCCvsM5 = 0.4426;
    bool isCTagT = false; //tight
    if((pfCCvsL0 > -0.5426 &&  pfCCvsB0 > -0.5426) ||(pfCCvsL1 > -0.5426 && pfCCvsB1 > -0.5426))isCTagL = true;
    if((pfCCvsL0 > dCCvsM1 && pfCCvsB0 > dCCvsM1) ||(pfCCvsL1 > dCCvsM1 && pfCCvsB1 > dCCvsM1))isCTagM1 = true;
    if((pfCCvsL0 > dCCvsM2 && pfCCvsB0 > dCCvsM2) ||(pfCCvsL1 > dCCvsM2 && pfCCvsB1 > dCCvsM2))isCTagM2 = true;
    if((pfCCvsL0 > dCCvsM3 && pfCCvsB0 > dCCvsM3) ||(pfCCvsL1 > dCCvsM3 && pfCCvsB1 > dCCvsM3))isCTagM3 = true;
    if((pfCCvsL0 > dCCvsM4 && pfCCvsB0 > dCCvsM4) ||(pfCCvsL1 > dCCvsM4 && pfCCvsB1 > dCCvsM4))isCTagM4 = true;
    if((pfCCvsL0 > dCCvsM5 && pfCCvsB0 > dCCvsM5) ||(pfCCvsL1 > dCCvsM5 && pfCCvsB1 > dCCvsM5))isCTagM5 = true;
    if((pfCCvsL0 > 0.5 && pfCCvsB0 > 0.5) ||(pfCCvsL1 > 0.5 && pfCCvsB1 > 0.5))isCTagT = true;
    //loose
    if(isCTagL)fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_CTagL", 200, 0, 1000, diJet_tag.mass(), evtWeight );
    else fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_noCTagL", 200, 0, 1000, diJet_tag.mass(), evtWeight );
    //medium
    if(isCTagM1)fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_CTagM1", 200, 0, 1000, diJet_tag.mass(), evtWeight );
    else fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_noCTagM1", 200, 0, 1000, diJet_tag.mass(), evtWeight );
    //medium
    if(isCTagM2)fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_CTagM2", 200, 0, 1000, diJet_tag.mass(), evtWeight );
    else fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_noCTagM2", 200, 0, 1000, diJet_tag.mass(), evtWeight );
    //medium
    if(isCTagM3)fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_CTagM3", 200, 0, 1000, diJet_tag.mass(), evtWeight );
    else fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_noCTagM3", 200, 0, 1000, diJet_tag.mass(), evtWeight );
    //medium
    if(isCTagM4)fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_CTagM4", 200, 0, 1000, diJet_tag.mass(), evtWeight );
    else fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_noCTagM4", 200, 0, 1000, diJet_tag.mass(), evtWeight );
    //medium
    if(isCTagM5)fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_CTagM5", 200, 0, 1000, diJet_tag.mass(), evtWeight );
    else fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_noCTagM5", 200, 0, 1000, diJet_tag.mass(), evtWeight );
    //tight
    if(isCTagT)fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_CTagT", 200, 0, 1000, diJet_tag.mass(), evtWeight );
    else fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_noCTagT", 200, 0, 1000, diJet_tag.mass(), evtWeight );
    
    //---------------------------------------------------//
    // cuts on chi2 and prob of KinFit
    //---------------------------------------------------//
    if(probOfKinFit > 0.1){
      if(chi2OfKinFit <100)
      fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_Chi100", 200, 0, 1000, diJet.mass(), evtWeight );
      if(chi2OfKinFit <10)
      fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_Chi10", 200, 0, 1000, diJet.mass(), evtWeight );
      if(chi2OfKinFit <5)
      fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_Chi5", 200, 0, 1000, diJet.mass(), evtWeight );
      if(chi2OfKinFit <2)
      fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_Chi2", 200, 0, 1000, diJet.mass(), evtWeight );
      if(chi2OfKinFit <1)
      fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_Chi1", 200, 0, 1000, diJet.mass(), evtWeight );
      if(chi2OfKinFit <0.5)
      fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_Chi0p5", 200, 0, 1000, diJet.mass(), evtWeight );
      if(chi2OfKinFit <0.2)
      fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_Chi0p2", 200, 0, 1000, diJet.mass(), evtWeight );
      if(chi2OfKinFit <0.1)
      fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_Chi0p1", 200, 0, 1000, diJet.mass(), evtWeight );
    }
  }//event loop
  cout<<"kfCount = "<<kfCount<<endl;
  f->Close(); 
  delete f;
}

void hplusAnalyzer::processEvents(){ 
  
  //Data, MC sample from lxplus and T2
  ///CutFlowAnalysis("TTJetsP_MuMC_20171104_Ntuple_1.root", "PF", ""); 
  //CutFlowAnalysis("outFile_.root", "PF", ""); 
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/", "PF", "");
  
  CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/rverma/ntuple_EleMC_kfitM_20171215/EleMC_20171215/TTJetsP_EleMC_20171215/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/TTJetsP_EleMC_20171215/171215_130507/0000/TTJetsP_EleMC_20171215_Ntuple_162.root", "PF", "");

  //====================================
  //condor submission
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", "PF", "outputFile");
  //====================================
} 
