
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
  /*
  CutFlowProcessor(url, myKey, "baseLowMET", outFile_);
  //to estimate unc in the data-driven qcd 
  CutFlowProcessor(url, myKey, "baseIso20HighMET", outFile_);
  CutFlowProcessor(url, myKey, "baseIso20LowMET", outFile_);
  //---------------------------------------------------//
  //for systematics (all sys in one go)
  //---------------------------------------------------//  
  if(!ev_->isData){
    CutFlowProcessor(url, myKey, "PileupPlus",     outFile_); 
    CutFlowProcessor(url, myKey, "PileupMinus",    outFile_);
    CutFlowProcessor(url, myKey, "JESPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "JESMinus", 	outFile_);
    CutFlowProcessor(url, myKey, "JERPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "JERMinus", 	outFile_);
    //CutFlowProcessor(url, myKey, "METUCPlus", 	outFile_);
    //CutFlowProcessor(url, myKey, "METUCMinus", 	outFile_);
    CutFlowProcessor(url, myKey, "bTagPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "bTagMinus", 	outFile_);
    CutFlowProcessor(url, myKey, "cTagPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "cTagMinus", 	outFile_);
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
  int jes = 0, jer = 0, metuc = 0, bScale = 0, cScale =0, minMET =20, minMT =0;
  //to estimate unc in the data-driven qcd 
  bool isLowMET = false, isIso20 = false;

  if(cutflowType.Contains("JESPlus"))jes = 1;
  else if (cutflowType.Contains("JESMinus"))jes = -1;
  else if (cutflowType.Contains("JERPlus"))jer = 1;
  else if (cutflowType.Contains("JERMinus"))jer = -1;
  else if (cutflowType.Contains("METUCPlus"))metuc = 1;
  else if (cutflowType.Contains("METUCMinus"))metuc = -1;
  else if (cutflowType.Contains("bTagPlus"))bScale = 1;
  else if (cutflowType.Contains("bTagMinus"))bScale = -1; 
  else if (cutflowType.Contains("cTagPlus"))cScale = 1;
  else if (cutflowType.Contains("cTagMinus"))cScale = -1; 
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
  double lumiTotal = 35860;
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
  
  //---------------------------------------------------//
  //BTag SF: read CSV file for SF, 2D histos for eff 
  //---------------------------------------------------//      
  //https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
  const std::string & bTagCSVfile 	= "stack/CSVv2_Moriond17_B_H.csv";
  const std::string & bTagName 		= "CSVv2";
  const std::string & bTagSys 		= "central"; 
  if(bScale==1) const std::string &bTagSys 		= "up"; 
  if(bScale==-1)const std::string &bTagSys 		= "down"; 
  const std::vector<std::string> & otherSysTypes = {"up", "down"};
  //b-quark
  BTagCalibrationReader readBTagCSV_bM= readCSV(bTagCSVfile, bTagName, BTagEntry::OP_MEDIUM,
    	      "comb", bTagSys, otherSysTypes, BTagEntry::FLAV_B);
  //c-quark
  BTagCalibrationReader readBTagCSV_cM= readCSV(bTagCSVfile, bTagName, BTagEntry::OP_MEDIUM,
    	      "comb", bTagSys, otherSysTypes, BTagEntry::FLAV_C);
  //other(light) quarks and gluon
  BTagCalibrationReader readBTagCSV_lM= readCSV(bTagCSVfile, bTagName, BTagEntry::OP_MEDIUM,
    	      "incl", bTagSys, otherSysTypes, BTagEntry::FLAV_UDSG);
  
  //getBTagEffHistos(f);
  TString histPath("myMiniTreeProducer/Jets/");
  TH2D* h2_BTagEff_Denom_b 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Denom_b"));
  TH2D* h2_BTagEff_Denom_c 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Denom_c"));
  TH2D* h2_BTagEff_Denom_udsg 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Denom_udsg")); 
  TH2D* h2_BTagEff_Num_bM 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Num_bM"));
  TH2D* h2_BTagEff_Num_cM 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Num_cM"));
  TH2D* h2_BTagEff_Num_udsgM 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Num_udsgM")); 
  
  //---------------------------------------------------//
  //CTag SF: read CSV file for SF, 2D histos for eff 
  //---------------------------------------------------//      
  //https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
  const std::string & cTagCSVfile 		= "stack/ctagger_Moriond17_B_H.csv";
  const std::string & cTagName 			= "cTag";
  const std::string & cTagSys 			= "central"; 
  if(cScale==1) const std::string &cTagSys 		= "up"; 
  if(cScale==-1)const std::string &cTagSys 		= "down"; 
  
  //LOOSE WP
  //b-quark
  BTagCalibrationReader readCTagCSV_bL= readCSV(cTagCSVfile, cTagName, BTagEntry::OP_LOOSE,
    	      "TnP", cTagSys, otherSysTypes, BTagEntry::FLAV_B);
  //c-quark
  BTagCalibrationReader readCTagCSV_cL= readCSV(cTagCSVfile, cTagName, BTagEntry::OP_LOOSE,
    	      "comb", cTagSys, otherSysTypes, BTagEntry::FLAV_C);
  //other(light) quarks and gluon
  BTagCalibrationReader readCTagCSV_lL= readCSV(cTagCSVfile, cTagName, BTagEntry::OP_LOOSE,
    	      "incl", cTagSys, otherSysTypes, BTagEntry::FLAV_UDSG);
  
  //MEDIUM WP
  //b-quark
  BTagCalibrationReader readCTagCSV_bM= readCSV(cTagCSVfile, cTagName, BTagEntry::OP_MEDIUM,
    	      "TnP", cTagSys, otherSysTypes, BTagEntry::FLAV_B);
  //c-quark
  BTagCalibrationReader readCTagCSV_cM= readCSV(cTagCSVfile, cTagName, BTagEntry::OP_MEDIUM,
    	      "comb", cTagSys, otherSysTypes, BTagEntry::FLAV_C);
  //other(light) quarks and gluon
  BTagCalibrationReader readCTagCSV_lM= readCSV(cTagCSVfile, cTagName, BTagEntry::OP_MEDIUM,
    	      "incl", cTagSys, otherSysTypes, BTagEntry::FLAV_UDSG);
  
  //TIGHT WP
  //b-quark
  BTagCalibrationReader readCTagCSV_bT= readCSV(cTagCSVfile, cTagName, BTagEntry::OP_TIGHT,
    	      "TnP", cTagSys, otherSysTypes, BTagEntry::FLAV_B);
  //c-quark
  BTagCalibrationReader readCTagCSV_cT= readCSV(cTagCSVfile, cTagName, BTagEntry::OP_TIGHT,
    	      "comb", cTagSys, otherSysTypes, BTagEntry::FLAV_C);
  //other(light) quarks and gluon
  BTagCalibrationReader readCTagCSV_lT= readCSV(cTagCSVfile, cTagName, BTagEntry::OP_TIGHT,
    	      "incl", cTagSys, otherSysTypes, BTagEntry::FLAV_UDSG);
  
  //getCTagEffHistos(f);
  TH2D* h2_CTagEff_Denom_b 		= (TH2D*)(f->Get(histPath+"h2_CTagEff_Denom_b"));
  TH2D* h2_CTagEff_Denom_c 		= (TH2D*)(f->Get(histPath+"h2_CTagEff_Denom_c"));
  TH2D* h2_CTagEff_Denom_udsg 		= (TH2D*)(f->Get(histPath+"h2_CTagEff_Denom_udsg")); 
  //LOOSE
  TH2D* h2_CTagEff_Num_bL 		= (TH2D*)(f->Get(histPath+"h2_CTagEff_Num_bL"));
  TH2D* h2_CTagEff_Num_cL 		= (TH2D*)(f->Get(histPath+"h2_CTagEff_Num_cL"));
  TH2D* h2_CTagEff_Num_udsgL 		= (TH2D*)(f->Get(histPath+"h2_CTagEff_Num_udsgL")); 
  //MEDIUM
  TH2D* h2_CTagEff_Num_bM 		= (TH2D*)(f->Get(histPath+"h2_CTagEff_Num_bM"));
  TH2D* h2_CTagEff_Num_cM 		= (TH2D*)(f->Get(histPath+"h2_CTagEff_Num_cM"));
  TH2D* h2_CTagEff_Num_udsgM 		= (TH2D*)(f->Get(histPath+"h2_CTagEff_Num_udsgM")); 
  //TIGHT
  TH2D* h2_CTagEff_Num_bT 		= (TH2D*)(f->Get(histPath+"h2_CTagEff_Num_bT"));
  TH2D* h2_CTagEff_Num_cT		= (TH2D*)(f->Get(histPath+"h2_CTagEff_Num_cT"));
  TH2D* h2_CTagEff_Num_udsgT 		= (TH2D*)(f->Get(histPath+"h2_CTagEff_Num_udsgT")); 
   
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
      fillHisto(outFile_, cutflowType, "", "LumiScaleFactor", 10, 0, 1000, sampleWeight, 1 );
      }
      //pileup weight
      double weightPU = 1.0;
      vector<double>pu = ev->sampleInfo.truepileup;
      if(pu.size() > 0) {
        float npu = pu[0];
        weightPU = LumiWeights_.weight(npu);
        if(cutflowType.Contains("PileupPlus"))
          weightPU = LumiWeights_Up.weight(npu);
        else if(cutflowType.Contains("PileupMinus"))
          weightPU = LumiWeights_Down.weight(npu);
      evtWeight *= weightPU;  
      }
    } 
    
    //---------------------------------------------------//
    // apply top re-weighting
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
    //---------------------------------------------------//
    double topPtWt = 1.0;
    if(!ev->isData){
      string sampleName = ev->sampleInfo.sampleName;
      if(sampleName.find("Hplus") != string::npos ||
		      sampleName.find("TTJetsM") != string::npos || 
		      sampleName.find("TTJetsP") != string::npos)
      {
        vector<double>topptweights = ev->sampleInfo.topPtWeights;
        if(topptweights.size() > 0){
            //topPtWt = topptweights[0];
          if(cutflowType.Contains("TopPtPlus"))
            topPtWt = topPtWt*topPtWt;
          else if(cutflowType.Contains("TopPtMinus"))
            topPtWt = 1.0;
        }
      }
    }
    fillHisto(outFile_, cutflowType, "", "SF_topPtWeights", 1000, 0, 3, topPtWt, 1 );
    evtWeight *= topPtWt; //Multiply to the total weights
    
    //---------------------------------------------------//
    //apply electron triggers
    //---------------------------------------------------//
    
    bool passTrig = false;
    vector<string> trig = ev->hlt;
    for(size_t it = 0; it < trig.size(); it++){
      if(trig[it].find("HLT_Ele27_WPTight_Gsf") != string::npos) passTrig = true;
    }
    if(!passTrig) continue;
    nTriggEvent++;
    double nCutPass = 1.0;
    fillHisto(outFile_, cutflowType+"/Iso", "", "cutflow", 10, 0.5, 10.5, nCutPass, evtWeight );
    fillHisto(outFile_, cutflowType+"/NonIso", "", "cutflow", 10, 0.5, 10.5, nCutPass, evtWeight );
   
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
    
    int nEle = e_final.size();
    double pri_vtxs = Vertices[0].totVtx;
    //select only one electron
    if(nEle != 1)continue;
    int e_i = e_final[0];
    //loose electron veto
    if(looseElectronVeto(e_final[0], pfElectrons, Vertices[0], isPFlow)) continue;
    double tmp_iso = pfElectrons[e_i].relCombPFIsoEA;
     
    //events should not have any muon
    if(looseMuonVeto( -1, pfMuons, isPFlow) ) continue;

    //---------------------------------------------------//
    // Iso(<0.08) and Non-iso(>0.08) region 
    //---------------------------------------------------//
    bool noisofound = false;
    bool isofound = false;
    string cutflowType_(cutflowType);
    if(isIso20){
      if(tmp_iso <= 0.11) cutflowType_ = cutflowType+"/Iso";
      if(tmp_iso > 0.11 && tmp_iso <= 0.30) cutflowType_ = cutflowType+"/NonIso";
    }
    else{
      if(tmp_iso <= 0.08) cutflowType_ = cutflowType+"/Iso";
      if(tmp_iso > 0.08 && tmp_iso <= 0.30) cutflowType_ = cutflowType+"/NonIso";
    }
    double eRelIso = pfElectrons[e_i].relCombPFIsoEA;
    double elePt = pfElectrons[e_i].p4.pt();
    
    //---------------------------------------------------//
    //apply Electron SF to eventWeights 
    //---------------------------------------------------//
    //Reco, ID, trigger	
    double ele_recoSF 		= getEleSF(h2_ele_recoSF, pfElectrons[e_i].eleSCEta, pfElectrons[e_i].p4.pt());
    double ele_medium_idSF  	= getEleSF(h2_ele_medium_idSF, pfElectrons[e_i].eleSCEta, pfElectrons[e_i].p4.pt());
    double ele_trigSF 		= getEleTrigSF(h2_ele_trigSF, pfElectrons[e_i].p4.pt(), pfElectrons[e_i].eleSCEta);
    
    //combined SF
    double eleSF =1.0;
    if(!ev->isData) eleSF = ele_recoSF*ele_medium_idSF*ele_trigSF;
    evtWeight *= eleSF;
    double metPt = 0; 
    metPt = metWithJESJER(pfJets, &j_final, met, jes, jer);
    Float_t xBinIso_array[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410, 420, 430, 440, 450, 460, 470, 480, 490, 500};
    fillHisto(outFile_, cutflowType, "", "RelIso_1Ele", 50, 0, 3, tmp_iso, evtWeight);
    fillHisto(outFile_, cutflowType, "", "pt_met_1Ele", 100, 0, 1000, metPt, evtWeight );
    fillHisto2D(outFile_, cutflowType, "", "RelIso_MET_1Ele", 100, 0, 1000, metPt, 50, 0, 3, tmp_iso,evtWeight );
    fillTProfile(outFile_, cutflowType, "", "RelIso_MET_TProf_1Ele", 50, xBinIso_array, metPt, tmp_iso, evtWeight); 
    
    if(tmp_iso > 0.30) continue;
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 10, 0.5, 10.5, nCutPass, evtWeight );

    //---------------------------------------------------//
    // Apply Jet Selection
    //---------------------------------------------------//
    int count_jets = j_final.size();
    if(count_jets < 4)continue;  // events should have 4 or more jets
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 10, 0.5, 10.5, nCutPass, evtWeight );
    
    //---------------------------------------------------//
    //apply MET selection   
    //---------------------------------------------------//
    //if(metuc==0)metPt = metWithJESJER(pfJets, &j_final, met, jes, jer);
    //else metPt = metWithUncl(pfJets, &j_final, pfMuons, &m_init, pfElectrons, &e_final, met, metuc);
    
    if(isLowMET){
      if(metPt > minMET) continue;  
    }
    else if(metPt < minMET) continue;  // Missing transverse energy cut 30 GeV(CMS) for ATLAS 20 GeV 
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 10, 0.5, 10.5, nCutPass, evtWeight );
    //Transverse mass b/w lepton and MET
    double deltaPhi(0);
    //leptonPt = TMath::Abs(pfMuons[m_i].p4.pt());
    deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(pfElectrons[e_i].p4, met.p4);
    double mt = sqrt (  2*elePt*metPt*(1 - cos(deltaPhi) ) ) ;

    //---------------------------------------------------//
    //apply B-tagging, C-tagging
    //---------------------------------------------------//
    vector<int> j_final_nob; j_final_nob.clear();
    vector<int> j_final_b; j_final_b.clear();
    vector<double> bdiscr; bdiscr.clear();
    double pfCISV = 0.0; //pfCombinedInclusiveSecondaryVertexV2BJetTags
    double pfCMVA = 0.0; //pfCombinedMVAV2BJetTags
    int count_CSVM_SF = 0; 
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      pfCISV = pfJets[ind_jet].bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
      pfCMVA = pfJets[ind_jet].bDiscriminator["pfCombinedMVAV2BJetTags"];
      fillHisto(outFile_, cutflowType_, "BTag", "pfCISV", 100, -2, 2, pfCISV, evtWeight );
      fillHisto(outFile_, cutflowType_, "BTag", "pfCMVA", 100, -2, 2, pfCMVA, evtWeight );
      //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
      if(pfCISV > 0.8484){
        count_CSVM_SF++; 
        double jetPt = jetPtWithJESJER(pfJets[ijet], jes, jer);
        fillHisto(outFile_, cutflowType_, "BTag", "pt_bjet", 100, 0, 1000, jetPt, evtWeight );
        fillHisto(outFile_, cutflowType_, "BTag", "eta_bjet", 50, -5, 5, pfJets[ijet].p4.eta(), evtWeight );
        j_final_b.push_back(ind_jet);
        bdiscr.push_back(pfCISV);
      }
      else j_final_nob.push_back(ind_jet);
    }
    if(count_CSVM_SF <= 1) continue; // Demanding for 2M b-tagged jets

    //Apply b-tag SF
    //https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1a)%20Event%20reweighting%20using%20scal
    double pmc_btag = 1.0;
    double pdata_btag = 1.0;
    double bTagWt = 1.0; 
    if(!ev->isData){
      for(size_t ijet = 0; ijet < j_final.size(); ijet++){
        int ind_jet = j_final[ijet];
        double pMC_ = 1.0;
        double pData_ = 1.0;
        //b-quark
        if(abs(pfJets[ind_jet].partonFlavour) ==5){
          pMC_ = getBTagPmcSys(h2_BTagEff_Num_bM, h2_BTagEff_Denom_b, pfJets[ind_jet]); 
          pData_ = getBTagPdataSys(readBTagCSV_bM, h2_BTagEff_Num_bM, h2_BTagEff_Denom_b, pfJets[ind_jet],bScale);
        }
        //c-quark
        else if(abs(pfJets[ind_jet].partonFlavour) ==4){ 
          pMC_ = getBTagPmcSys(h2_BTagEff_Num_cM, h2_BTagEff_Denom_c, pfJets[ind_jet]); 
          pData_ = getBTagPdataSys(readBTagCSV_cM, h2_BTagEff_Num_cM, h2_BTagEff_Denom_c, pfJets[ind_jet],bScale);
        }
        //other quarks and gluon
        else{ 
          pMC_ = getBTagPmcSys(h2_BTagEff_Num_udsgM, h2_BTagEff_Denom_udsg, pfJets[ind_jet]); 
          pData_ = getBTagPdataSys(readBTagCSV_lM, h2_BTagEff_Num_udsgM, h2_BTagEff_Denom_udsg, pfJets[ind_jet], bScale); 
        }
        pmc_btag = pmc_btag*pMC_;
        pdata_btag = pdata_btag*pData_;
      }
    bTagWt = pdata_btag/pmc_btag;
    }
    evtWeight *= bTagWt;
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 10, 0.5, 10.5, nCutPass, evtWeight);
    fillHisto(outFile_, cutflowType, "", "bTagWeight", 100, 0, 2, bTagWt, 1);

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
        fillHisto(outFile_, cutflowType_, "BTag", "mjj", 100, 0, 500, diJet.M(), evtWeight );
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
      fillHisto(outFile_, cutflowType_, "BTag", "mjj", 100, 0, 500, diJet.M(), evtWeight );
    }	
    //---------------------------------------------------//
    // add set of plots after BTag:
    //---------------------------------------------------//
    //fillHisto("pt_mu", cutflowType_+"/BTag", pfMuons[m_i].p4.pt(), evtWeight);
    fillHisto(outFile_, cutflowType_, "BTag","pt_ele", 100, 0, 1000, elePt, evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","eta_ele", 50, -5, 5, pfElectrons[e_i].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","phi_ele", 50, -5, 5, pfElectrons[e_i].p4.phi(), evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","final_RelIso_ele", 100, 0, 1, eRelIso, evtWeight );
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      fillHisto(outFile_, cutflowType_, "BTag","pt_jet", 100, 0, 1000, jetPt, evtWeight );
      fillHisto(outFile_, cutflowType_, "BTag","eta_jet", 50, -5, 5, pfJets[ind_jet].p4.eta(), evtWeight );
      fillHisto(outFile_, cutflowType_, "BTag","phi_jet", 50, -5, 5, pfJets[ind_jet].p4.phi(), evtWeight );
    }
    fillHisto(outFile_, cutflowType_, "BTag","final_multi_jet", 15, 0.5, 15.5, count_jets, evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag", "CSVL_count", 10, 0.5, 10.5, count_CSVM_SF, evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","final_pt_met", 100, 0, 1000, metPt, evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","wmt", 100, 0, 1000, mt, evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","nvtx", 100, 0, 100, pri_vtxs, evtWeight );
    for(std::size_t n=0; n<Vertices.size(); n++){
      fillHisto(outFile_, cutflowType_, "BTag","rhoAll", 100, 0, 100, Vertices[n].rhoAll, evtWeight );
    }
    fillHisto(outFile_, cutflowType_, "BTag", "chi2OfKinFit", 100, 0, 100, chi2OfKinFit, evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag", "probOfKinFit", 100, 0, 3, probOfKinFit, evtWeight );
    input_count++;
    if(input_count%10==0)
    cout << "input count iso: "<< input_count << endl;
    //if(i > 2000) break;
  
    //---------------------------------------------------//
    //get KF lepton by matching with PF leptons
    //---------------------------------------------------//
    bool status_ = false; 
    if(statusOfKinFit ==0) status_=true ;
    bool foundkfEle = false;
    if(status_){
      double dR =DeltaR(pfElectrons[e_i].p4 , kfLepton[0]);
      if(kfLepton.size()>0){
        if(dR< 0.2)foundkfEle = true;
      }
    }
    bool foundkfEle_ = false; 
    bool kfJetsSel_=false; 
    if(status_ && foundkfEle && kfLepton[0].pt() >25){
      foundkfEle_=true;
      if(kfJets[0].pt() > 25 && 
              kfJets[1].pt() > 25 && 
              kfJets[2].pt() > 25 && 
              kfJetsLepB[0].pt() > 25 )kfJetsSel_ =true;
    }
    //---------------------------------------------------//
    //select maximum b-tag discriminator jet in KF
    //---------------------------------------------------//
    vector<MyLorentzVector> kfLightJets; kfLightJets.clear();
    double pt_bjetHad = 0;
    if(status_ && foundkfEle_ && kfJetsSel_){
      unsigned long maxBtagJet = -1;
      double maxBDiscr = -999.;
      int count_kfJets = 0;
      for(unsigned long ik = 0; ik < kfJets.size(); ik++){
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
      unsigned long zero = 0;
      if(kfJets.size() >=3 && maxBtagJet >= zero){
        for(unsigned long ik = 0; ik < kfJets.size(); ik++){
          if(ik != maxBtagJet)kfLightJets.push_back(kfJets[ik]);
      else pt_bjetHad = kfJets[ik].pt();
        }
      }
    }
    bool kfLightJetSel_=false; 
    if(kfLightJets.size() >= 2) kfLightJetSel_ = true;
    bool match_j1 = false, match_j2 = false;
    int indexForCTag0 = 0, indexForCTag1 = 0;
    if(status_ && foundkfEle_ && kfJetsSel_ && kfLightJetSel_){
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
    }
    //---------------------------------------------------//
    // All KinFit selections
    //make sure that the fit converges 
    bool matchJets_ =false;
    if(match_j1 && match_j2) matchJets_ = true;
    bool allKinFitSel_ = false; 
    if(status_ && foundkfEle_ && kfJetsSel_ && matchJets_) allKinFitSel_=true;
    //---------------------------------------------------//
    if(allKinFitSel_){
      kfCount++;
      nCutPass++;
      fillHisto(outFile_, cutflowType_, "", "cutflow", 10, 0.5, 10.5, nCutPass, evtWeight );
      MyLorentzVector diJet = kfLightJets[0]+kfLightJets[1];
      fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit", 100, 0, 500, diJet.mass(), evtWeight );
      fillHisto(outFile_, cutflowType_, "KinFit","pt_ele", 100, 0, 1000, elePt, evtWeight );
      fillHisto(outFile_, cutflowType_, "KinFit","eta_ele", 50, -5, 5, pfElectrons[e_i].p4.eta(), evtWeight );
      fillHisto(outFile_, cutflowType_, "KinFit","phi_ele", 50, -5, 5, pfElectrons[e_i].p4.phi(), evtWeight );
      fillHisto(outFile_, cutflowType_, "KinFit","final_RelIso_ele", 100, 0, 1, eRelIso, evtWeight );
      for(size_t ijet = 0; ijet < j_final.size(); ijet++){
        int ind_jet = j_final[ijet];
        double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
        fillHisto(outFile_, cutflowType_, "KinFit","pt_jet", 100, 0, 1000, jetPt, evtWeight );
        fillHisto(outFile_, cutflowType_, "KinFit","eta_jet", 50, -5, 5, pfJets[ind_jet].p4.eta(), evtWeight );
        fillHisto(outFile_, cutflowType_, "KinFit","phi_jet", 50, -5, 5, pfJets[ind_jet].p4.phi(), evtWeight );
      }
      for( std::size_t n=0; n<Vertices.size(); n++){
        fillHisto(outFile_, cutflowType_, "KinFit","rhoAll", 100, 0, 100, Vertices[n].rhoAll, evtWeight );
      }
      fillHisto(outFile_, cutflowType_, "KinFit","final_multi_jet", 15, 0.5, 15.5, count_jets, evtWeight );
      fillHisto(outFile_, cutflowType_, "KinFit", "CSVL_count", 10, 0.5, 10.5, count_CSVM_SF, evtWeight );
      fillHisto(outFile_, cutflowType_, "KinFit","final_pt_met", 100, 0, 1000, metPt, evtWeight );
      fillHisto(outFile_, cutflowType_, "KinFit","wmt", 100, 0, 1000, mt, evtWeight );
      fillHisto(outFile_, cutflowType_, "KinFit","nvtx", 100, 0, 100, pri_vtxs, evtWeight );

      //---------------------------------------------------//
      // b-jet pT categorization from inclusive events
      //---------------------------------------------------//
      fillHisto(outFile_, cutflowType_, "KinFit","pt_bjetH", 100, 0, 1000, pt_bjetHad, evtWeight );
      fillHisto(outFile_, cutflowType_, "KinFit","pt_bjetL", 100, 0, 1000, kfJetsLepB[0].pt(), evtWeight );
      Float_t xBin_array[] = {25,  42,  57,  74,  99,  500};
      fillTProfile(outFile_, cutflowType_, "PtbJetInc", "mjj_kfit_pt_bjetH", 5, xBin_array, pt_bjetHad, 
              diJet.mass(), evtWeight); 
      if(pt_bjetHad >= 25 && pt_bjetHad < 42)
        fillHisto(outFile_, cutflowType_, "PtbJetInc", "mjj_kfit_25To42", 100, 0, 500, diJet.mass(), evtWeight );
      if(pt_bjetHad >= 42 && pt_bjetHad < 57)
        fillHisto(outFile_, cutflowType_, "PtbJetInc", "mjj_kfit_42To57", 100, 0, 500, diJet.mass(), evtWeight );
      if(pt_bjetHad >= 57 && pt_bjetHad < 74)
        fillHisto(outFile_, cutflowType_, "PtbJetInc", "mjj_kfit_57To74", 100, 0, 500, diJet.mass(), evtWeight );
      if(pt_bjetHad >= 74 && pt_bjetHad < 99)
        fillHisto(outFile_, cutflowType_, "PtbJetInc", "mjj_kfit_74To99", 100, 0, 500, diJet.mass(), evtWeight );
      if(pt_bjetHad >= 99 && pt_bjetHad < 500)
        fillHisto(outFile_, cutflowType_, "PtbJetInc", "mjj_kfit_99To500", 100, 0, 500, diJet.mass(), evtWeight );
      
      //---------------------------------------------------//
      //apply CTagging
      //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
      //---------------------------------------------------//
      int count_cJetsIncL   = 0;
      int count_cJetsIncM   = 0;
      int count_cJetsIncT   = 0;
      for(size_t ijet = 0; ijet < kfLightJets.size(); ijet++){
        if(ijet ==0) ijet = indexForCTag0;
        if(ijet ==1) ijet = indexForCTag1;
        double pfCCvsL = pfJets[ijet].bDiscriminator["pfCombinedCvsLJetTags"];
        double pfCCvsB = pfJets[ijet].bDiscriminator["pfCombinedCvsBJetTags"]; 
        fillHisto(outFile_, cutflowType_, "KinFit", "pfCCvsL", 100, -2, 2, pfCCvsL, evtWeight );
        fillHisto(outFile_, cutflowType_, "KinFit", "pfCCvsB", 100, -2, 2, pfCCvsB, evtWeight );
        if(pfCCvsL > -0.48  && pfCCvsB > -0.17)  count_cJetsIncL++;
        if(pfCCvsL > -0.1   && pfCCvsB > 0.08)   count_cJetsIncM++;
        if(pfCCvsL > 0.69   && pfCCvsB > -0.45)  count_cJetsIncT++;
      }
      //---------------------------------------------------//
      //Charm mistag scale factors : Inclusive 
      //---------------------------------------------------//
      //https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1a)%20Event%20reweighting%20using%20scal
      double cTagWt_IncL  = 1.0;
      double cTagWt_IncM  = 1.0;
      double cTagWt_IncT  = 1.0;
      if(!ev->isData){
        bool isCTagIncL = false;
        if(count_cJetsIncL > 0){
          isCTagIncL = true; 
          double pMC_IncL     = 1.0;
          double pData_IncL   = 1.0;
          for(size_t ijet = 0; ijet < kfLightJets.size(); ijet++){
            if(ijet ==0) ijet = indexForCTag0;
            if(ijet ==1) ijet = indexForCTag1;
            double pmc   = 1.0;
            double pdata = 1.0;
            //b-quark
            if(abs(pfJets[ijet].partonFlavour) ==5){
              pmc = getIncCTagPmcSys(h2_CTagEff_Num_bL, h2_CTagEff_Denom_b, pfJets[ijet], isCTagIncL);   
              pdata = getIncCTagPdataSys(readCTagCSV_bL, h2_CTagEff_Num_bL, h2_CTagEff_Denom_b, 
                      pfJets[ijet], isCTagIncL, cScale);
            }
            //c-quark
            else if(abs(pfJets[ijet].partonFlavour) ==4){ 
              pmc = getIncCTagPmcSys(h2_CTagEff_Num_cL, h2_CTagEff_Denom_c, pfJets[ijet], isCTagIncL);   
              pdata = getIncCTagPdataSys(readCTagCSV_cL, h2_CTagEff_Num_cL, h2_CTagEff_Denom_c, 
                      pfJets[ijet], isCTagIncL, cScale);
            }
            //other quarks and gluon
            else{ 
              pmc = getIncCTagPmcSys(h2_CTagEff_Num_udsgL, h2_CTagEff_Denom_udsg, pfJets[ijet], isCTagIncL);   
              pdata = getIncCTagPdataSys(readCTagCSV_lL, h2_CTagEff_Num_udsgL, h2_CTagEff_Denom_udsg, 
                      pfJets[ijet], isCTagIncL, cScale);
            }
            pMC_IncL = pMC_IncL*pmc;
            pData_IncL = pData_IncL*pdata;
          }
          cTagWt_IncL  = pData_IncL/pMC_IncL; 
        }//Inc loose scale factor
        bool isCTagIncM = false;
        if(count_cJetsIncM > 0){
          isCTagIncM = true;
          double pMC_IncM     = 1.0;
          double pData_IncM   = 1.0;
          for(size_t ijet = 0; ijet < kfLightJets.size(); ijet++){
            if(ijet ==0) ijet = indexForCTag0;
            if(ijet ==1) ijet = indexForCTag1;
            double pmc   = 1.0;
            double pdata = 1.0;
            //b-quark
            if(abs(pfJets[ijet].partonFlavour) ==5){
              pmc = getIncCTagPmcSys(h2_CTagEff_Num_bM, h2_CTagEff_Denom_b, pfJets[ijet], isCTagIncM);   
              pdata = getIncCTagPdataSys(readCTagCSV_bM, h2_CTagEff_Num_bM, h2_CTagEff_Denom_b, 
                      pfJets[ijet], isCTagIncM, cScale);
            }
            //c-quark
            else if(abs(pfJets[ijet].partonFlavour) ==4){ 
              pmc = getIncCTagPmcSys(h2_CTagEff_Num_cM, h2_CTagEff_Denom_c, pfJets[ijet], isCTagIncM);   
              pdata = getIncCTagPdataSys(readCTagCSV_cM, h2_CTagEff_Num_cM, h2_CTagEff_Denom_c, 
                      pfJets[ijet], isCTagIncM, cScale);
            }
            //other quarks and gluon
            else{ 
              pmc = getIncCTagPmcSys(h2_CTagEff_Num_udsgM, h2_CTagEff_Denom_udsg, pfJets[ijet], isCTagIncM);   
              pdata = getIncCTagPdataSys(readCTagCSV_lM, h2_CTagEff_Num_udsgM, h2_CTagEff_Denom_udsg, 
                      pfJets[ijet], isCTagIncM, cScale);
            }
            pMC_IncM = pMC_IncM*pmc;
            pData_IncM = pData_IncM*pdata;
          }
          cTagWt_IncM  = pData_IncM/pMC_IncM; 
        }//Inc medium scale factor
        bool isCTagT = false;
        if(count_cJetsIncT > 0){
          isCTagT = true;
          double pMC_IncT     = 1.0;
          double pData_IncT   = 1.0;
          for(size_t ijet = 0; ijet < kfLightJets.size(); ijet++){
            if(ijet ==0) ijet = indexForCTag0;
            if(ijet ==1) ijet = indexForCTag1;
            double pmc   = 1.0;
            double pdata = 1.0;
            //b-quark
            if(abs(pfJets[ijet].partonFlavour) ==5){
              pmc = getIncCTagPmcSys(h2_CTagEff_Num_bT, h2_CTagEff_Denom_b, pfJets[ijet], isCTagT);   
              pdata = getIncCTagPdataSys(readCTagCSV_bT, h2_CTagEff_Num_bT, h2_CTagEff_Denom_b, 
                      pfJets[ijet], isCTagT, cScale);
            }
            //c-quark
            else if(abs(pfJets[ijet].partonFlavour) ==4){ 
              pmc = getIncCTagPmcSys(h2_CTagEff_Num_cT, h2_CTagEff_Denom_c, pfJets[ijet], isCTagT);   
              pdata = getIncCTagPdataSys(readCTagCSV_cT, h2_CTagEff_Num_cT, h2_CTagEff_Denom_c, 
                      pfJets[ijet], isCTagT, cScale);
            }
            //other quarks and gluon
            else{ 
              pmc = getIncCTagPmcSys(h2_CTagEff_Num_udsgT, h2_CTagEff_Denom_udsg, pfJets[ijet], isCTagT);   
              pdata = getIncCTagPdataSys(readCTagCSV_lT, h2_CTagEff_Num_udsgT, h2_CTagEff_Denom_udsg, 
                      pfJets[ijet], isCTagT, cScale);
            }
            pMC_IncT = pMC_IncT*pmc;
            pData_IncT = pData_IncT*pdata;
          }
          cTagWt_IncT  = pData_IncT/pMC_IncT; 
        }
      }
      fillHisto(outFile_, cutflowType, "", "cTagWeight_incL", 100, 0, 2, cTagWt_IncL, 1);
      fillHisto(outFile_, cutflowType, "", "cTagWeight_incM", 100, 0, 2, cTagWt_IncM, 1);
      fillHisto(outFile_, cutflowType, "", "cTagWeight_incT", 100, 0, 2, cTagWt_IncT, 1);
      double evtWtIncL = evtWeight*cTagWt_IncL;
      double evtWtIncM = evtWeight*cTagWt_IncM;
      double evtWtIncT = evtWeight*cTagWt_IncT;
      if(count_cJetsIncL > 0) 
          fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_CTagIncL", 100, 0, 500, diJet.mass(), evtWtIncL);
      if(count_cJetsIncM > 0) 
          fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_CTagIncM", 100, 0, 500, diJet.mass(), evtWtIncM);
      if(count_cJetsIncT > 0) 
          fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_CTagIncT", 100, 0, 500, diJet.mass(), evtWtIncT);
      if(count_cJetsIncL > 0){
      nCutPass++;
      fillHisto(outFile_, cutflowType_, "", "cutflow", 10, 0.5, 10.5, nCutPass, evtWtIncL);
      }
      //---------------------------------------------------//
      //ctag - categorisation 
      //Charm mistag scale factors : Exclusive 
      //https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1a)%20Event%20reweighting%20using%20scal
      //---------------------------------------------------//
      bool isCTagL = false;
      bool isCTagM = false;
      bool isCTagT = false;
      double cTagWt_ExL  = 1.0;
      double cTagWt_ExM  = 1.0;
      double cTagWt_ExT  = 1.0;
      if(!ev->isData){
        if(count_cJetsIncL > 0){
          isCTagL = true; 
          isCTagM = false; 
          isCTagT = false; 
          double pMC_ExL     = 1.0;
          double pData_ExL   = 1.0;
          for(size_t ijet = 0; ijet < kfLightJets.size(); ijet++){
            if(ijet ==0) ijet = indexForCTag0;
            if(ijet ==1) ijet = indexForCTag1;
            double pmc   = 1.0;
            double pdata = 1.0;
            //b-quark
            if(abs(pfJets[ijet].partonFlavour) ==5){
              pmc = getExCTagPmcSys(h2_CTagEff_Num_bL, h2_CTagEff_Num_bM, h2_CTagEff_Num_bT, 
                      h2_CTagEff_Denom_b, pfJets[ijet], isCTagL, isCTagM, isCTagT);   
              pdata = getExCTagPdataSys(readCTagCSV_bL, readCTagCSV_bM, readCTagCSV_bT, 
                      h2_CTagEff_Num_bL, h2_CTagEff_Num_bM, h2_CTagEff_Num_bT, h2_CTagEff_Denom_b, 
                      pfJets[ijet], isCTagL, isCTagM, isCTagT, cScale);
            }
            //c-quark
            else if(abs(pfJets[ijet].partonFlavour) ==4){ 
              pmc = getExCTagPmcSys(h2_CTagEff_Num_cL, h2_CTagEff_Num_cM, h2_CTagEff_Num_cT, 
                      h2_CTagEff_Denom_c, pfJets[ijet], isCTagL, isCTagM, isCTagT);   
              pdata = getExCTagPdataSys(readCTagCSV_cL, readCTagCSV_cM, readCTagCSV_cT, 
                      h2_CTagEff_Num_cL, h2_CTagEff_Num_cM, h2_CTagEff_Num_cT, h2_CTagEff_Denom_c, 
                      pfJets[ijet], isCTagL, isCTagM, isCTagT, cScale);
            }
            //other quarks and gluon
            else{ 
              pmc = getExCTagPmcSys(h2_CTagEff_Num_udsgL, h2_CTagEff_Num_udsgM, h2_CTagEff_Num_udsgT, 
                      h2_CTagEff_Denom_udsg, pfJets[ijet], isCTagL, isCTagM, isCTagT);   
              pdata = getExCTagPdataSys(readCTagCSV_lL, readCTagCSV_lM, readCTagCSV_lT, 
                      h2_CTagEff_Num_udsgL, h2_CTagEff_Num_udsgM, h2_CTagEff_Num_udsgT, h2_CTagEff_Denom_udsg, 
                      pfJets[ijet], isCTagL, isCTagM, isCTagT, cScale);
            }
            pMC_ExL = pMC_ExL*pmc;
            pData_ExL = pData_ExL*pdata;
          }
          cTagWt_ExL  = pData_ExL/pMC_ExL; 
          cout<<"-------------------------"<<endl;
          cout<<"pData_ExL = "<<pData_ExL<<endl;
          cout<<"pMC_ExL = "<<pMC_ExL<<endl;
          cout<<"cTagWt_ExL = "<<cTagWt_ExL<<endl;
        }//Ex loose scale factor
        if(count_cJetsIncM > 0){
          isCTagL = false; 
          isCTagM = true; 
          isCTagT = false; 
          double pMC_ExM     = 1.0;
          double pData_ExM   = 1.0;
          for(size_t ijet = 0; ijet < kfLightJets.size(); ijet++){
            if(ijet ==0) ijet = indexForCTag0;
            if(ijet ==1) ijet = indexForCTag1;
            //b-quark
            double pmc   = 1.0;
            double pdata = 1.0;
            if(abs(pfJets[ijet].partonFlavour) ==5){
              pmc = getExCTagPmcSys(h2_CTagEff_Num_bL, h2_CTagEff_Num_bM, h2_CTagEff_Num_bT, 
                      h2_CTagEff_Denom_b, pfJets[ijet], isCTagL, isCTagM, isCTagT);   
              pdata = getExCTagPdataSys(readCTagCSV_bL, readCTagCSV_bM, readCTagCSV_bT, 
                      h2_CTagEff_Num_bL, h2_CTagEff_Num_bM, h2_CTagEff_Num_bT, h2_CTagEff_Denom_b, 
                      pfJets[ijet], isCTagL, isCTagM, isCTagT, cScale);
            }
            //c-quark
            else if(abs(pfJets[ijet].partonFlavour) ==4){ 
              pmc = getExCTagPmcSys(h2_CTagEff_Num_cL, h2_CTagEff_Num_cM, h2_CTagEff_Num_cT, 
                      h2_CTagEff_Denom_c, pfJets[ijet], isCTagL, isCTagM, isCTagT);   
              pdata = getExCTagPdataSys(readCTagCSV_cL, readCTagCSV_cM, readCTagCSV_cT, 
                      h2_CTagEff_Num_cL, h2_CTagEff_Num_cM, h2_CTagEff_Num_cT, h2_CTagEff_Denom_c, 
                      pfJets[ijet], isCTagL, isCTagM, isCTagT, cScale);
            }
            //other quarks and gluon
            else{ 
              pmc = getExCTagPmcSys(h2_CTagEff_Num_udsgL, h2_CTagEff_Num_udsgM, h2_CTagEff_Num_udsgT, 
                      h2_CTagEff_Denom_udsg, pfJets[ijet], isCTagL, isCTagM, isCTagT);   
              pdata = getExCTagPdataSys(readCTagCSV_lL, readCTagCSV_lM, readCTagCSV_lT, 
                      h2_CTagEff_Num_udsgL, h2_CTagEff_Num_udsgM, h2_CTagEff_Num_udsgT, h2_CTagEff_Denom_udsg, 
                      pfJets[ijet], isCTagL, isCTagM, isCTagT, cScale);
            }
            pMC_ExM = pMC_ExM*pmc;
            pData_ExM = pData_ExM*pdata;
          }
          cTagWt_ExM  = pData_ExM/pMC_ExM; 
        }//Ex medium scale factor
        if(count_cJetsIncT > 0){
          isCTagL = false; 
          isCTagM = false; 
          isCTagT = true; 
          double pMC_ExT     = 1.0;
          double pData_ExT   = 1.0;
          for(size_t ijet = 0; ijet < kfLightJets.size(); ijet++){
            if(ijet ==0) ijet = indexForCTag0;
            if(ijet ==1) ijet = indexForCTag1;
            double pmc   = 1.0;
            double pdata = 1.0;
            //b-quark
            if(abs(pfJets[ijet].partonFlavour) ==5){
              pmc = getExCTagPmcSys(h2_CTagEff_Num_bL, h2_CTagEff_Num_bM, h2_CTagEff_Num_bT, 
                      h2_CTagEff_Denom_b, pfJets[ijet], isCTagL, isCTagM, isCTagT);   
              pdata = getExCTagPdataSys(readCTagCSV_bL, readCTagCSV_bM, readCTagCSV_bT, 
                      h2_CTagEff_Num_bL, h2_CTagEff_Num_bM, h2_CTagEff_Num_bT, h2_CTagEff_Denom_b, 
                      pfJets[ijet], isCTagL, isCTagM, isCTagT, cScale);
            }
            //c-quark
            else if(abs(pfJets[ijet].partonFlavour) ==4){ 
              pmc = getExCTagPmcSys(h2_CTagEff_Num_cL, h2_CTagEff_Num_cM, h2_CTagEff_Num_cT, 
                      h2_CTagEff_Denom_c, pfJets[ijet], isCTagL, isCTagM, isCTagT);   
              pdata = getExCTagPdataSys(readCTagCSV_cL, readCTagCSV_cM, readCTagCSV_cT, 
                      h2_CTagEff_Num_cL, h2_CTagEff_Num_cM, h2_CTagEff_Num_cT, h2_CTagEff_Denom_c, 
                      pfJets[ijet], isCTagL, isCTagM, isCTagT, cScale);
            }
            //other quarks and gluon
            else{ 
              pmc = getExCTagPmcSys(h2_CTagEff_Num_udsgL, h2_CTagEff_Num_udsgM, h2_CTagEff_Num_udsgT, 
                      h2_CTagEff_Denom_udsg, pfJets[ijet], isCTagL, isCTagM, isCTagT);   
              pdata = getExCTagPdataSys(readCTagCSV_lL, readCTagCSV_lM, readCTagCSV_lT, 
                      h2_CTagEff_Num_udsgL, h2_CTagEff_Num_udsgM, h2_CTagEff_Num_udsgT, h2_CTagEff_Denom_udsg, 
                      pfJets[ijet], isCTagL, isCTagM, isCTagT, cScale);
            }
            pMC_ExT = pMC_ExT*pmc;
            pData_ExT = pData_ExT*pdata;
          }
        cTagWt_ExT  = pData_ExT/pMC_ExT; 
        }
      }
      fillHisto(outFile_, cutflowType, "", "cTagWeight_ExL", 100, 0, 2, cTagWt_ExL, 1);
      fillHisto(outFile_, cutflowType, "", "cTagWeight_ExM", 100, 0, 2, cTagWt_ExM, 1);
      fillHisto(outFile_, cutflowType, "", "cTagWeight_ExT", 100, 0, 2, cTagWt_ExT, 1);
      double evtWtExL = evtWeight*cTagWt_ExL;
      double evtWtExM = evtWeight*cTagWt_ExM;
      double evtWtExT = evtWeight*cTagWt_ExT;

      //---------------------------------------------------//
      // b-jet pT categorization from exclusive events
      //---------------------------------------------------//
      if(count_cJetsIncT> 0){ 
        fillTProfile(outFile_, cutflowType_, "PtbJetCTagExT", "mjj_kfit_pt_bjetH", 5, xBin_array, pt_bjetHad, 
                diJet.mass(), evtWtExT); 
        fillHisto(outFile_, cutflowType_, "PtbJetCTagExT","pt_bjetH_CTagExT", 100, 0, 1000, pt_bjetHad, evtWtExT );
        fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_CTagExT", 100, 0, 500, diJet.mass(), evtWtExT );
        if(pt_bjetHad >= 25 && pt_bjetHad < 42)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExT", "mjj_kfit_25To42", 100, 0, 500, diJet.mass(), evtWtExT );
        if(pt_bjetHad >= 42 && pt_bjetHad < 57)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExT", "mjj_kfit_42To57", 100, 0, 500, diJet.mass(), evtWtExT );
        if(pt_bjetHad >= 57 && pt_bjetHad < 74)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExT", "mjj_kfit_57To74", 100, 0, 500, diJet.mass(), evtWtExT );
        if(pt_bjetHad >= 74 && pt_bjetHad < 99)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExT", "mjj_kfit_74To99", 100, 0, 500, diJet.mass(), evtWtExT );
        if(pt_bjetHad >= 99 && pt_bjetHad < 500)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExT", "mjj_kfit_99To500", 100, 0, 500, diJet.mass(), evtWtExT );
        }
      else if(count_cJetsIncM> 0){ 
        fillTProfile(outFile_, cutflowType_, "PtbJetCTagExM", "mjj_kfit_pt_bjetH", 5, xBin_array, pt_bjetHad, 
                diJet.mass(), evtWtExM); 
        fillHisto(outFile_, cutflowType_, "PtbJetCTagExM","pt_bjetH_CTagExM", 100, 0, 1000, pt_bjetHad, evtWtExM );
        fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_CTagExM", 100, 0, 500, diJet.mass(), evtWtExM );
        if(pt_bjetHad >= 25 && pt_bjetHad < 42)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExM", "mjj_kfit_25To42", 100, 0, 500, diJet.mass(), evtWtExM );
        if(pt_bjetHad >= 42 && pt_bjetHad < 57)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExM", "mjj_kfit_42To57", 100, 0, 500, diJet.mass(), evtWtExM );
        if(pt_bjetHad >= 57 && pt_bjetHad < 74)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExM", "mjj_kfit_57To74", 100, 0, 500, diJet.mass(), evtWtExM );
        if(pt_bjetHad >= 74 && pt_bjetHad < 99)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExM", "mjj_kfit_74To99", 100, 0, 500, diJet.mass(), evtWtExM );
        if(pt_bjetHad >= 99 && pt_bjetHad < 500)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExM", "mjj_kfit_99To500", 100, 0, 500, diJet.mass(), evtWtExM );
      }
      else if(count_cJetsIncL > 0){ 
        fillTProfile(outFile_, cutflowType_, "PtbJetCTagExL", "mjj_kfit_pt_bjetH", 5, xBin_array, pt_bjetHad, 
                diJet.mass(), evtWtExL); 
        fillHisto(outFile_, cutflowType_, "PtbJetCTagExL","pt_bjetH_CTagExL", 100, 0, 1000, pt_bjetHad, evtWtExL );
        fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_CTagExL", 100, 0, 500, diJet.mass(), evtWtExL );
        if(pt_bjetHad >= 25 && pt_bjetHad < 42)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExL", "mjj_kfit_25To42", 100, 0, 500, diJet.mass(), evtWtExL );
        if(pt_bjetHad >= 42 && pt_bjetHad < 57)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExL", "mjj_kfit_42To57", 100, 0, 500, diJet.mass(), evtWtExL );
        if(pt_bjetHad >= 57 && pt_bjetHad < 74)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExL", "mjj_kfit_57To74", 100, 0, 500, diJet.mass(), evtWtExL );
        if(pt_bjetHad >= 74 && pt_bjetHad < 99)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExL", "mjj_kfit_74To99", 100, 0, 500, diJet.mass(), evtWtExL );
        if(pt_bjetHad >= 99 && pt_bjetHad < 500)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExL", "mjj_kfit_99To500", 100, 0, 500, diJet.mass(), evtWtExL );
      }
      else{
        fillTProfile(outFile_, cutflowType_, "PtbJetCTagExO", "mjj_kfit_pt_bjetH", 5, xBin_array, pt_bjetHad, 
                diJet.mass(), evtWeight); 
        fillHisto(outFile_, cutflowType_, "PtbJetCTagExO","pt_bjetH_CTagExO", 100, 0, 1000, pt_bjetHad, evtWeight );
        fillHisto(outFile_, cutflowType_, "KinFit", "mjj_kfit_CTagExO", 100, 0, 500, diJet.mass(), evtWeight );
        if(pt_bjetHad >= 25 && pt_bjetHad < 42)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExO", "mjj_kfit_25To42", 100, 0, 500, diJet.mass(), evtWeight );
        if(pt_bjetHad >= 42 && pt_bjetHad < 57)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExO", "mjj_kfit_42To57", 100, 0, 500, diJet.mass(), evtWeight );
        if(pt_bjetHad >= 57 && pt_bjetHad < 74)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExO", "mjj_kfit_57To74", 100, 0, 500, diJet.mass(), evtWeight );
        if(pt_bjetHad >= 74 && pt_bjetHad < 99)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExO", "mjj_kfit_74To99", 100, 0, 500, diJet.mass(), evtWeight );
        if(pt_bjetHad >= 99 && pt_bjetHad < 500)
          fillHisto(outFile_, cutflowType_, "PtbJetCTagExO", "mjj_kfit_99To500", 100, 0, 500, diJet.mass(), evtWeight );
      }
    }//allKinFitSel
  else{
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
        fillHisto(outFile_, cutflowType_, "nonKinFit", "mjj_non_kfit", 100, 0, 500, diJet.M(), evtWeight );
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
      fillHisto(outFile_, cutflowType_, "nonKinFit", "mjj_non_kfit", 100, 0, 500, diJet.M(), evtWeight );
      }	
      fillHisto(outFile_, cutflowType_, "nonKinFit","pt_ele", 100, 0, 1000, elePt, evtWeight );
      fillHisto(outFile_, cutflowType_, "nonKinFit","eta_ele", 50, -5, 5, pfElectrons[e_i].p4.eta(), evtWeight );
      fillHisto(outFile_, cutflowType_, "nonKinFit","phi_ele", 50, -5, 5, pfElectrons[e_i].p4.phi(), evtWeight );
      fillHisto(outFile_, cutflowType_, "nonKinFit","final_RelIso_ele", 100, 0, 1, eRelIso, evtWeight );
      for(size_t ijet = 0; ijet < j_final.size(); ijet++){
        int ind_jet = j_final[ijet];
        double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
        fillHisto(outFile_, cutflowType_, "nonKinFit","pt_jet", 100, 0, 1000, jetPt, evtWeight );
        fillHisto(outFile_, cutflowType_, "nonKinFit","eta_jet", 50, -5, 5, pfJets[ind_jet].p4.eta(), evtWeight );
        fillHisto(outFile_, cutflowType_, "nonKinFit","phi_jet", 50, -5, 5, pfJets[ind_jet].p4.phi(), evtWeight );
      }
      fillHisto(outFile_, cutflowType_, "nonKinFit","final_multi_jet", 15, 0.5, 15.5, count_jets, evtWeight );
      fillHisto(outFile_, cutflowType_, "nonKinFit", "CSVL_count", 10, 0.5, 10.5, count_CSVM_SF, evtWeight );
      fillHisto(outFile_, cutflowType_, "nonKinFit","final_pt_met", 100, 0, 1000, metPt, evtWeight );
      fillHisto(outFile_, cutflowType_, "nonKinFit","wmt", 100, 0, 1000, mt, evtWeight );
      fillHisto(outFile_, cutflowType_, "nonKinFit","nvtx", 100, 0, 100, pri_vtxs, evtWeight );
      for(std::size_t n=0; n<Vertices.size(); n++){
        fillHisto(outFile_, cutflowType_, "nonKinFit","rhoAll", 100, 0, 100, Vertices[n].rhoAll, evtWeight );
      }
    }//nonKinFitSel
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

  CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/rverma/ntuple_EleMC_kfitM_20180419/EleMC_20180419/TTJetsP_EleMC_20180419/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/TTJetsP_EleMC_20180419/180419_065438/0000/TTJetsP_EleMC_20180419_Ntuple_1.root", "PF", "");

  //CutFlowAnalysis("/home/rverma/t3store2/ele_unProcNtuples/EleRunDv1_EleData_20180418_Ntuple_1.root", "PF", ""); 
  //====================================
  //condor submission
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", "PF", "outputFile");
  //====================================
} 
